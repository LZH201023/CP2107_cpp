#include <NTT.hpp>
#include <ctime>

mcl::Fr* iter_NTT(const mcl::Fr* f, int N)
{
    if (N == 1)
    {
        g = inv_g = 1;
        return new mcl::Fr[1] { f[0] };
    }
    mcl::Fr* arr = new mcl::Fr[N];
    int* arr1 = new int[N];
    int T = N / 2;
    int H = N;
    int temp = 1;
    int logN = log2(N);

    arr1[0] = 0;

    for (int i = 0; i < logN; i++)
    {
        for (int j = 0; j < N; j += H)
            arr1[j + T] = arr1[j] + temp;
        T /= 2;
        H /= 2;
        temp *= 2;
    }

    for (int i = 0; i < N; i++)
        arr[i] = f[arr1[i]];
    T = 2;
    mcl::Fr* temp_arr = new mcl::Fr[N];
    mcl::Fr TthRoot = -1;
    for (int i = 0; i < logN; i++)
    {
        for (int j = 0; j < N; j += T)
        {
            mcl::Fr t = 1;
            for (int k = 0; k < H; k++)
            {
                temp_arr[k] = arr[j + k] + arr[j + k + H] * t;
                temp_arr[k + H] = arr[j + k] - arr[j + k + H] * t;
                t *= TthRoot;
            }
            for (int k = 0; k < T; k++)
                arr[j + k] = temp_arr[k];
        }
        T *= 2;
        H *= 2;
        mcl::Fr::squareRoot(TthRoot, TthRoot);
    }
    g = TthRoot * TthRoot;
    mcl::Fr::inv(inv_g, g);
    delete[] arr1; delete[] temp_arr;
    delete[] f;
    return arr;
}

mcl::Fr* inv_NTT(const mcl::Fr* f, int N)
{
    mcl::Fr* arr = new mcl::Fr[N];
    std::copy(f, f + N, arr);
    mcl::Fr* res = new mcl::Fr[N];
    int logN = log2(N);
    int T = N; int H = N / 2;
    mcl::Fr w = inv_g;

    for (int i = 0; i < logN; i++) {
        for (int j = 0; j < N; j += T)
        {
            mcl::Fr a = 1;
            for (int k = 0; k < H; k++)
            {
                res[j + k] = (arr[j + k] + arr[j + k + H]) / 2;
                res[j + k + H] = (arr[j + k] - arr[j + k + H]) / 2 * a;
                a *= w;
            }
        }
        w *= w;
        T /= 2;
        H /= 2;
        std::swap(arr, res);
    }

    int* temp = new int[N];
    temp[0] = 0; int t = 1;
    for (H = N / 2; H > 0; H /= 2) {
        for (int j = 0; j < t; j++)
            temp[j + t] = temp[j] + H;
        t *= 2;
    }
    for (int i = 0; i < N; i++)
        res[i] = arr[temp[i]];
    delete[] temp; delete[] arr;
    return res;
}

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

    fN = iter_NTT(fN, N); hN = iter_NTT(hN, N);
    for (int i = 0; i < N; i++) fN[i] = fN[i] * hN[i];
    mcl::Fr* res = inv_NTT(fN, N);
    delete[] fN; delete[] hN;
    return res;
}