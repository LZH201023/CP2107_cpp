#pragma once

#include <mcl/bls12_381.hpp>
#include <vector>

struct Witness
{
    mcl::Fr i;
    mcl::Fr phi_i;
    mcl::G1 witness;

    Witness(mcl::Fr i, mcl::Fr phi_i, mcl::G1 witness)
    {
        this->i = i;
        this->phi_i = phi_i;
        this->witness = witness;
    }

};

struct PK
{
    mcl::G1* srs;
    mcl::G2* hs;
    int max_deg;
    
    PK(mcl::G1* srs, mcl::G2* hs, int t)
    {
        this->srs = srs;
        this->hs = hs;
        this->max_deg = t;
    }

    ~PK()
    {
        delete[] srs;
        delete[] hs;
    }
};

struct Witness_B1
{
    const mcl::Fr* B;
    mcl::Fr* r;
    mcl::G1 witness;
    int B_size;

    Witness_B1(const mcl::Fr* B, mcl::Fr* r, mcl::G1 witness, int B_size)
    {
        this->B = B;
        this->r = r;
        this->witness = witness;
        this->B_size = B_size;
    }

    ~Witness_B1()
    {
        delete[] B;
        delete[] r;
    }
};

struct Witness_B2
{
    mcl::Fr i;
    mcl::Fr* evals;
    mcl::G1 witness;
    int num;

    Witness_B2(mcl::Fr i, mcl::Fr* evals, mcl::G1 witness, int num)
    {
        this->i = i;
        this->evals = evals;
        this->witness = witness;
        this->num = num;
    }

    ~Witness_B2()
    {
        delete[] this->evals;
    }

};

PK Setup(int k, int t);
mcl::G1 Commit(const PK& pk, const mcl::Fr* poly, int df);
Witness CreateWitness(const PK& pk, const mcl::Fr* poly, int df, const mcl::Fr& i);
bool VerifyEval(const PK& pk, const mcl::G1& c, const Witness& w);

Witness_B1 CreateWitness_B1(const PK& pk, const mcl::Fr* poly, int df, const mcl::Fr* B, int B_size);
bool VerifyEval_B1(const PK& pk, const mcl::G1& c, const Witness_B1& w);

Witness_B2 CreateWitness_B2(const PK& pk, mcl::Fr** polys, int num, const int* degs, const mcl::Fr& i, mcl::Fr& r);
bool VerifyEval_B2(const PK& pk, const mcl::G1* cs, int num, const Witness_B2& w, mcl::Fr& r);