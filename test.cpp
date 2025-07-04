#include <mcl/bls12_381.hpp>
#include <NTT.hpp>
#include <KZG.hpp>
#include <ctime>

extern bool UZPIOP(bool tag = 1);
extern bool USCPIOP(bool tag = 1);

void NTT_test();
void NTT_time();
void KZG_test();
void protocol_test();

int man()
{
    mcl::initPairing(mcl::BN_SNARK1);
    NTT_time();
}

static mcl::Fr* naive_poly_mult(mcl::Fr* f, int df, mcl::Fr* g, int dg)
{
    int fs = df + 1;
    int gs = dg + 1;
    int N = fs + gs - 1;
    mcl::Fr* res = new mcl::Fr[N];
    for (int i = 0; i < N; i++)
        res[i] = 0;
    for (int i = 0; i < N; i++)
    {
        int j = i - gs + 1 > 0 ? i - gs + 1 : 0;
        for (;j <= i && j < fs; j++)
            res[i] = (res[i] + f[j] * g[i - j]);
    }
    return res;
}

void NTT_test()
{
    for (int i = 0; i < 5; i++)
    {
        int df = rand() % 2048 + 1;
        int dg = rand() % 2048 + 1;
        mcl::Fr* f = new mcl::Fr[df + 1];
        mcl::Fr* g = new mcl::Fr[dg + 1];
        for (int j = 0; j <= df; j++)
            f[j].setByCSPRNG();
        for (int j = 0; j <= dg; j++)
            g[j].setByCSPRNG();
        mcl::Fr* result = naive_poly_mult(f, df, g, dg);
        mcl::Fr* got = poly_mult(f, df, g, dg);
        delete[] f; delete[] g;
        bool check = true;
        int N = df + dg + 1;
        for (int i = 0; i < N; i++)
            check = check && (result[i] == got[i]);
        delete[] result; delete[] got;
        std::cout << "NTT test " << i + 1 << ": " << (check ? "Success" : "Failure") << "\n";
    }
}

void NTT_time()
{
    clock_t timer;
    clock_t total = 0;
    for (int i = 0; i < 6; i++)
    {
        int df = 65535;
        int dg = 65536;
        mcl::Fr* f = new mcl::Fr[df + 1];
        mcl::Fr* g = new mcl::Fr[dg + 1];
        for (int j = 0; j <= df; j++)
            f[j].setByCSPRNG();
        for (int j = 0; j <= dg; j++)
            g[j].setByCSPRNG();
        timer = clock();
        mcl::Fr* got = poly_mult(f, df, g, dg);
        timer = clock() - timer;
        delete[] f; delete[] g; delete[] got;
        if (i > 0) total += timer;
    }
    std::cout << "NTT average: " << total * 1000.0 / 5 / CLOCKS_PER_SEC << "ms\n";
}
//avg. 471.6ms

void KZG_test()
{
    mcl::initPairing(mcl::BN_SNARK1);
    for (int m = 0; m < 6; m++)
    {
        printf("KZG (%s) test%d: ", (m < 3 ? "completeness" : "soundness"), m + 1);
        bool tag = m < 3;
        int t = rand() % 100;
        if (t < 3) t += 10;
        mcl::Fr* f = new mcl::Fr[t + 1];
        for (int i = 0; i <= t; i++)
            f[i].setByCSPRNG();
        if (f[t] == 0) f[t] = f[t] + 1;
        PK pk = Setup(128, t);
        mcl::G1 c = Commit(pk, f, t);
        mcl::Fr i; i.setByCSPRNG();
        Witness w = CreateWitness(pk, f, t, i);
        if (!tag)
        {
            int ran = rand() % 2;
            mcl::Fr r; r.setByCSPRNG();
            if (ran == 0)
                w.phi_i = w.phi_i + r;
            else
                w.witness = w.witness * r;
        }
        printf(tag == VerifyEval(pk, c, w) ? "Success\n" : "Failure\n");
        delete[] f;
    }
}


void protocol_test()
{
    std::cout << "Univariate ZeroTest PIOP:\n";
    for (int i = 0; i < 10; i++)
        std::cout << "Test " << i + 1 << (i < 5 ? " (completeness)" : " (soundness)") << ": " << ((i < 5) == UZPIOP(i < 5) ? "Success\n" : "Failure\n");
    std::cout << "\nUnivariate SumCheck PIOP:\n";
    for (int i = 0; i < 10; i++)
        std::cout << "Test " << i + 1 << (i < 5 ? " (completeness)" : " (soundness)") << ": " << ((i < 5) == USCPIOP(i < 5) ? "Success\n" : "Failure\n");
}