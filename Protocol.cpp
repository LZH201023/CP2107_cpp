#include <KZG.hpp>
#include <NTT.hpp>
#include <mcl/bls12_381.hpp>
#include <iostream>

//univariate zero test PIOP
bool UZPIOP(bool tag = 1);
//univariate sum check PIOP
bool USCPIOP(bool tag = 1);

static mcl::Fr evalPoly(const mcl::Fr* poly, int df, const mcl::Fr& i);
static mcl::Fr evalZx(int order, const mcl::Fr& x);
static void poly_div_Z(const mcl::Fr* f, int df, const int H_size, mcl::Fr* Q, mcl::Fr* r);
static mcl::Fr* gen_H(int pow);

bool UZPIOP(bool tag) //1 for honest P, 0 for malicious P
{
    clock_t timer;
    clock_t P_time = 0, V_time = 0;
    //initialize H
    int pow = 10;
    int l = 1;
    l = l << pow;
    mcl::Fr* H = gen_H(pow);

    //initialize f in ZeroTest, using randomly generated polynomial
    int t = 2047;
    if (t < 3) t += 10;
    mcl::Fr* q = new mcl::Fr[t + 1];
    for (int i = 0; i <= t; i++)
        q[i].setByCSPRNG();
    mcl::Fr* Z = new mcl::Fr[l + 1];
    for (int k = 1; k < l; k++)
        Z[k] = 0;
    Z[0] = -1; Z[l] = 1;
    mcl::Fr* f = poly_mult(q, t, Z, l);
    std::cout << "Univariate ZeroTest PIOP\n\n";
    int df = t + l;
    if (!tag)
    {
        mcl::Fr r;
        r.setByCSPRNG();
        if (r == 0) r += 1;
        f[rand() % (df + 1)] += r;
    }

    //setup KZG
    timer = clock();
    PK pk = Setup(128, t + l);
    timer = clock() - timer;
    std::cout << "Setup time: " << timer * 1000.0 / CLOCKS_PER_SEC << "ms\n\n";

    //P commits to f
    timer = clock();
    mcl::G1 PCS_f = Commit(pk, f, df);
    P_time += clock() - timer;
    std::cout << "Prover's commitment:\nf: " << PCS_f.getStr() << "\n";

    //P commits to q (f = Z * q)
    mcl::Fr* r = new mcl::Fr[df + 1];
    timer = clock();
    poly_div_Z(f, df, l, q, r);
    mcl::G1 PCS_q = Commit(pk, q, t);
    P_time += clock() - timer;
    std::cout << "q: " << PCS_q.getStr() << "\n\n";

    //V sends challenge x, randomness ra
    timer = clock();
    mcl::Fr x, ra; x.setByCSPRNG(); ra.setByCSPRNG();
    V_time += clock() - timer;
    std::cout << "Verifier: x = " << x.getStr() << ", r = " << ra.getStr() << "\n\n";

    //P sends evaluations and witness
    timer = clock();
    mcl::Fr** polys = new mcl::Fr*[2];
    polys[0] = f; polys[1] = q;
    int degs[2];
    degs[0] = df; degs[1] = t;
    Witness_B2 w = CreateWitness_B2(pk, polys, 2, degs, x, ra);
    P_time += clock() - timer;
    std::cout << "Prover:\nf(x) = " << w.evals[0] << "\nq(x) = " << w.evals[1]
        << "\nWitness: " << w.witness.getStr() << "\n\n";

    //V verifies evaluations
    timer = clock();
    mcl::Fr Zx;
    bool check = true;
    mcl::G1 cs[2];
    cs[0] = PCS_f; cs[1] = PCS_q;
    check = check && VerifyEval_B2(pk, cs, 2, w, ra);
    Zx = evalZx(l, x);
    check = check && (w.evals[0] == w.evals[1] * Zx);
    V_time += clock() - timer;
    if (check)
        std::cout << "Verifier accepts\n";
    else std::cout << "Verifier rejects\n";

    std::cout << "Protocol ends\n\n";

    std::cout << "Prover time: " << P_time * 1000.0 / CLOCKS_PER_SEC << "ms\n"
        << "Verifier time: " << V_time * 1000.0 / CLOCKS_PER_SEC << "ms\n";

    delete[] H;
    delete[] q; delete[] Z; delete[] f;
    delete[] polys;
    return check;
}

bool USCPIOP(bool tag) //1 for honest P, 0 for malicious P
{
    clock_t timer;
    clock_t P_time = 0, V_time = 0;
    //initialize H
    int pow = 10;
    int l = 1;
    l = l << pow;
    mcl::Fr* H = gen_H(pow);

    //initialize f in SumCheck, using randomly generated polynomial
    int t = 2047;
    if (t < l) t += l;
    mcl::Fr* f = new mcl::Fr[t + 1];
    for (int i = 0; i <= t; i++)
        f[i].setByCSPRNG();
    mcl::Fr res = 0;
    for (int k = 0; k < l; k++)
        res += evalPoly(f, t, H[k]);
    f[0] = f[0] - res / l + (1 - tag) * rand();
    std::cout << "Univariate SumCheck PIOP\n\n";

    //setup KZG
    timer = clock();
    PK pk = Setup(128, t);
    timer = clock() - timer;
    std::cout << "Setup time: " << timer * 1000.0 / CLOCKS_PER_SEC << "ms\n\n";

    //P commits to f
    timer = clock();
    mcl::G1 PCS_f = Commit(pk, f, t);
    P_time += clock() - timer;
    std::cout << "Prover's commitment:\nf: " << PCS_f.getStr() << "\n";
 
    timer = clock();
    mcl::Fr* Z = new mcl::Fr[l + 1];
    Z[0] = -1; Z[l] = 1;
    mcl::Fr* q = new mcl::Fr[t + 1 - l];
    mcl::Fr* p = new mcl::Fr[t + 1];
    poly_div_Z(f, t, l, q, p);
    
    //P commits to q (f = Z * q + x * p)
    mcl::G1 PCS_q = Commit(pk, q, t - l);
    P_time += clock() - timer;
    std::cout << "q: " << PCS_q.getStr() << "\n";

    //P commits to p (f = Z * q + x * p)
    timer = clock();
    mcl::G1 PCS_p = Commit(pk, p + 1, l - 2);
    P_time += clock() - timer;
    std::cout << "p: " << PCS_p.getStr() << "\n\n";

    //V sends challenge x, randomness ra
    timer = clock();
    mcl::Fr x, ra; x.setByCSPRNG(); ra.setByCSPRNG();
    V_time += clock() - timer;
    std::cout << "Verifier: x = " << x.getStr() << ", r = " << ra.getStr() << "\n\n";

    //P sends evaluations and witness
    timer = clock();
    mcl::Fr** polys = new mcl::Fr*[3];
    polys[0] = f; polys[1] = q; polys[2] = p + 1;
    int degs[3];
    degs[0] = t; degs[1] = t - l; degs[2] = l - 2;
    Witness_B2 w = CreateWitness_B2(pk, polys, 3, degs, x, ra);
    P_time += clock() - timer;
    std::cout << "Prover:\nf(x) = " << w.evals[0] << "\nq(x) = " << w.evals[1]
        << "\nWitness: " << w.witness.getStr() << "\n\n";

    //V verifies evaluations
    timer = clock();
    mcl::Fr Zx;
    bool check = true;
    mcl::G1 cs[3];
    cs[0] = PCS_f; cs[1] = PCS_q; cs[2] = PCS_p;
    check = check && VerifyEval_B2(pk, cs, 3, w, ra);
    Zx = evalZx(l, x);
    check = check && (w.evals[0] == w.evals[1] * Zx + x * w.evals[2]);
    V_time += clock() - timer;
    if (check)
        std::cout << "Verifier accepts\n";
    else std::cout << "Verifier rejects\n";

    std::cout << "Protocol ends\n\n";

    std::cout << "Prover time: " << P_time * 1000.0 / CLOCKS_PER_SEC << "ms\n"
        << "Verifier time: " << V_time * 1000.0 / CLOCKS_PER_SEC << "ms\n";

    delete[] H;
    delete[] f; delete[] Z; delete[] q; delete[] p;
    delete[] polys;
    return check;
}

static mcl::Fr evalZx(int order, const mcl::Fr& x)
{
    mcl::Fr p;
    mcl::Fr::pow(p, x, order);
    return p - 1;
}

static void poly_div_Z(const mcl::Fr* f, int df, int m, mcl::Fr* Q, mcl::Fr* r)
{
    int n = df + 1;
    std::copy(f, f + n, r);
    for (int i = n - m - 1; i >= 0; --i) {
        Q[i] = r[i + m];
        r[i] = r[i] + Q[i];
    }
}

static mcl::Fr evalPoly(const mcl::Fr* poly, int df, const mcl::Fr& i)
{
    mcl::Fr x = 1;
    mcl::Fr sum = 0;
    int n = df + 1;
    for (int j = 0; j < n; j++)
    {
        sum += x * poly[j];
        x *= i;
    }
    return sum;
}

//pow >= 0
static mcl::Fr* gen_H(int pow)
{
    if (pow == 0) return new mcl::Fr{1};
    mcl::Fr root = -1;
    for (int i = 1; i < pow; i++)
        mcl::Fr::squareRoot(root, root);
    int id = 1;
    id = id << pow;
    mcl::Fr* res = new mcl::Fr[id];
    res[0] = 1;
    for (int i = 1; i < id; i++)
        res[i] = res[i - 1] * root;
    return res;
}