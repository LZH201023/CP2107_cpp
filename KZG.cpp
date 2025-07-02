#include <KZG.hpp>
#include <NTT.hpp>

static mcl::Fr* calc_B_poly(const mcl::Fr*, int, int);
static void poly_div(const mcl::Fr* f, int df, const mcl::Fr* Z, int dZ, mcl::Fr* Q, mcl::Fr* r);
static mcl::Fr evalPoly(const mcl::Fr* poly, int df, const mcl::Fr& i);

PK Setup(int k, int t)
{
    if (k > 128)
    {
        std::cout << "Unable to provide required security level with BN_SNARK1!";
        exit(-1);
    }

    mcl::G1 g; mcl::Fr alpha, aux;
    mcl::Fp r; r.setByCSPRNG();
    mapToG1(g, r);
    alpha.setByCSPRNG();
    aux.setByCSPRNG();
    mcl::G2 h; mcl::hashAndMapToG2(h, aux.getStr());
    if (alpha == 0) alpha = alpha + 1;
    mcl::G1* v1 = new mcl::G1[t + 1];
    mcl::G2* v2 = new mcl::G2[t + 1];
    mcl::Fr p = 1;
    for (int i = 0; i <= t; i++)
    {
        v1[i] = g * p;
        v2[i] = h * p;
        p = p * alpha;
    }

    return PK(v1, v2, t);
}

mcl::G1 Commit(const PK& pk, const mcl::Fr* poly, int df)
{
    if (df > pk.max_deg)
    {
        std::cout << "Invalid input polynomial!";
        exit(-1);
    }
    mcl::G1 c;
    mcl::G1::mulVec(c, pk.srs, poly, df + 1);
    return c;
}

Witness CreateWitness(const PK& pk, const mcl::Fr* poly, int df, const mcl::Fr& i)
{
    mcl::Fr phi_i = evalPoly(poly, df, i);
    if (df == 0) return Witness(i, phi_i, mcl::getG1basePoint());
    mcl::Fr* ps1 = new mcl::Fr[df];
    ps1[df - 1] = poly[df];
    for (int j = df - 2; j >= 0; j--)
        ps1[j] = poly[j + 1] + i * ps1[j + 1];
    mcl::G1 witness;
    mcl::G1::mulVec(witness, pk.srs, ps1, df);
    return Witness(i, phi_i, witness);
}

bool VerifyEval(const PK& pk, const mcl::G1& c, const Witness& w)
{
    mcl::Fp12 e1;
    mcl::pairing(e1, w.witness, pk.hs[1] - (pk.hs[0] * w.i));
    mcl::Fp12 e2;
    mcl::pairing(e2, c - (pk.srs[0] * w.phi_i), pk.hs[0]);
    return e1 == e2; 
}

Witness_B1 CreateWitness_B1(const PK& pk, const mcl::Fr* poly, int df, const mcl::Fr* B, int B_size)
{
    if (B_size > df)
    {
        mcl::Fr* r = new mcl::Fr[df + 1];
        std::copy(poly, poly + df + 1, r);
        return Witness_B1(B, r, mcl::getG1basePoint(), B_size);
    }
    mcl::Fr* B_poly = calc_B_poly(B, 0, B_size);
    mcl::Fr* Q = new mcl::Fr[df - B_size + 1];
    mcl::Fr* r = new mcl::Fr[df + 1];
    poly_div(poly, df, B_poly, B_size, Q, r);
    mcl::G1 witness; 
    mcl::G1::mulVec(witness, pk.srs, Q, df - B_size + 1);
    delete[] B_poly; delete[] Q;
    return Witness_B1(B, r, witness, B_size);
}

bool VerifyEval_B1(const PK& pk, const mcl::G1& c, const Witness_B1& w)
{
    mcl::Fp12 e1, e2, e3;
    mcl::pairing(e1, c, pk.hs[0]);
    mcl::G2 r_alpha;
    mcl::G2::mulVec(r_alpha, pk.hs, w.r, std::min(w.B_size, pk.max_deg + 1));
    mcl::pairing(e2, pk.srs[0], r_alpha);
    mcl::Fr* Z = calc_B_poly(w.B, 0, w.B_size);
    mcl::G2 h_pi;
    mcl::G2::mulVec(h_pi, pk.hs, Z, std::min(w.B_size, pk.max_deg) + 1);
    mcl::pairing(e3, w.witness, h_pi);
    delete[] Z;
    return e1 == e2 * e3;
}

Witness_B2 CreateWitness_B2(const PK& pk, mcl::Fr** polys, int num, const int* degs, const mcl::Fr& i, mcl::Fr& r)
{
    mcl::Fr* evals = new mcl::Fr[num];
    for (int j = 0; j < num; j++)
        evals[j] = evalPoly(polys[j], degs[j], i);
    int max = 0;
    for (int k = 0; k < num; k++)
        if (degs[k] > max) max = degs[k];
    max++;
    mcl::Fr* p = new mcl::Fr[max];
    for (int k = 0; k < max; k++)
        p[k] = 0;
    mcl::Fr temp = 1;
    for (int k = 0; k < num; k++)
    {
        for (int j = 0; j < degs[k] + 1; j++)
            p[j] += polys[k][j] * temp;
        temp *= r;
    }
    mcl::Fr* ps1 = new mcl::Fr[max - 1];
    ps1[max - 2] = p[max - 1];
    for (int j = max - 3; j >= 0; j--)
        ps1[j] = p[j + 1] + i * ps1[j + 1];
    mcl::G1 witness;
    mcl::G1::mulVec(witness, pk.srs, ps1, max - 1);
    delete[] p; delete[] ps1;
    return Witness_B2(i, evals, witness, num);
}

bool VerifyEval_B2(const PK& pk, const mcl::G1* cs, int num, const Witness_B2& w, mcl::Fr& r)
{
    mcl::Fp12 e1;
    mcl::pairing(e1, w.witness, pk.hs[1] - (pk.hs[0] * w.i));
    mcl::Fp12 e2;
    mcl::G1 c = mcl::getG1basePoint();
    mcl::Fr temp = 1;
    mcl::Fr y = evalPoly(w.evals, w.num - 1, r);
    //optimize?
    for (int i = 0; i < num; i++)
    {
        c += cs[i] * temp;
        temp *= r;
    }
    mcl::pairing(e2, c - (pk.srs[0] * y), pk.hs[0]);
    return e1 == e2;
}

static mcl::Fr* calc_B_poly(const mcl::Fr* B, int begin, int end)
{
    if (end - begin == 1)
    {
        mcl::Fr* v = new mcl::Fr[2];
        v[0] = -1 * B[begin]; v[1] = 1;
        return v;
    }
    else
    {
        int mid = begin + (end - begin) / 2;
        mcl::Fr* v1, *v2, *vt;
        v1 = calc_B_poly(B, begin, mid);
        v2 = calc_B_poly(B, mid, end);
        vt = poly_mult(v1, mid - begin, v2, end - mid);
        delete[] v1; delete[] v2;
        return vt;
    }
}

static void poly_div(const mcl::Fr* f, int df, const mcl::Fr* Z, int dZ, mcl::Fr* Q, mcl::Fr* r)
{
    int n = df + 1, m = dZ + 1;
    std::copy(f, f + n, r);
    for (int i = n - m; i >= 0; --i) {
        mcl::Fr t = r[i + m - 1] / Z[m - 1];
        Q[i] = t;
        for (int j = i; j < i + m; ++j)
            r[j] = r[j] - t * Z[j - i];
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