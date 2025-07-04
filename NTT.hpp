#pragma once
#include <mcl/bls12_381.hpp>

static mcl::Fr g, inv_g;
static int order = 2;
static bool started = false;

template<typename S>
S* NTT(const S* f, int N)
{
    started = true;
    if (N == 1)
    {
        g = 1;
        return new S[1]{ f[0] };
    }
    S* arr = new S[N];
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
    S* temp_arr = new S[N];
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
    order = N;
    mcl::Fr::inv(inv_g, g);
    delete[] arr1; delete[] temp_arr;
    return arr;
}

template<typename S>
S* iNTT(const S* f, int N, mcl::Fr ge = 0)
{
    if (!started)
    {
        g = inv_g = -1;
        started = true;
    }
    if (ge == 0)
    {
        if (order != N)
        {
            if (order > N)
            {
                while (order > N)
                {
                    order /= 2;
                    g *= g;
                }
                mcl::Fr::inv(inv_g, g);
            }
            else
            {
                while (order < N)
                {
                    if (order == 1)
                    {
                        order *= 2;
                        g = -1;
                        continue;
                    }
                    order *= 2;
                    mcl::Fr::squareRoot(g, g);
                }
                mcl::Fr::inv(inv_g, g);
            }
        }
    }
    else
    {
        order = N;
        g = ge;
        mcl::Fr::inv(inv_g, ge);
    }
    S* arr = new S[N];
    std::copy(f, f + N, arr);
    S* res = new S[N];
    int logN = log2(N);
    int T = N; int H = N / 2;
    mcl::Fr w = inv_g;
    mcl::Fr half = 2; half = 1 / half;
    for (int i = 0; i < logN; i++) {
        for (int j = 0; j < N; j += T)
        {
            mcl::Fr a = 1;
            for (int k = 0; k < H; k++)
            {
                res[j + k] = (arr[j + k] + arr[j + k + H]) * half;
                res[j + k + H] = (arr[j + k] - arr[j + k + H]) * half * a;
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

mcl::Fr* poly_mult(const mcl::Fr*, int, const mcl::Fr*, int);