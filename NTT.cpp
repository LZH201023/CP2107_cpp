#include <NTT.hpp>

mcl::Fr* poly_mult(const mcl::Fr* f, int df, const mcl::Fr* h, int dh)
{
    int n = df + dh + 1;
    int N = 1;
    while (N < n) N = N << 1;

    mcl::Fr* fN = new mcl::Fr[N];
    mcl::Fr* hN = new mcl::Fr[N];
    std::copy(f, f + df + 1, fN);
    for (int i = df + 1; i < N; i++) fN[i] = 0;
    std::copy(h, h + dh + 1, hN);
    for (int i = dh + 1; i < N; i++) hN[i] = 0;

    fN = NTT<mcl::Fr>(fN, N); hN = NTT<mcl::Fr>(hN, N);
    for (int i = 0; i < N; i++) fN[i] = fN[i] * hN[i];
    mcl::Fr* res = iNTT<mcl::Fr>(fN, N);
    delete[] fN; delete[] hN;
    return res;
}