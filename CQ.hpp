#pragma once
#include <mcl/bls12_381.hpp>

struct SRS
{
	int N;
	mcl::G1* srs1;
	mcl::G2* srs2;
	mcl::G2 Zx2;
	mcl::G2 Tx2;
	mcl::G1* qxs1;
	mcl::G1* Lxs1;
	mcl::G1* Lps1;

	SRS(int N, mcl::G1* srs1, mcl::G2* srs2, mcl::G2 Zx2, mcl::G2 Tx2, mcl::G1* qxs1, mcl::G1* Lxs1, mcl::G1* Lps1)
	{
		this->N = N;
		this->srs1 = srs1;
		this->srs2 = srs2;
		this->Zx2 = Zx2;
		this->Tx2 = Tx2;
		this->qxs1 = qxs1;
		this->Lxs1 = Lxs1;
		this->Lps1 = Lps1;
	}

	~SRS()
	{
		delete[] srs1;
		delete[] srs2;
		delete[] qxs1;
		delete[] Lxs1;
		delete[] Lps1;
	}
};