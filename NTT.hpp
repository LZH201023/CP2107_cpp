#pragma once
#include <mcl/bls12_381.hpp>

static mcl::Fr g, inv_g;

mcl::Fr* iter_NTT(const mcl::Fr*, int);
mcl::Fr* inv_NTT(const mcl::Fr*, int);
mcl::Fr* poly_mult(const mcl::Fr*, int, const mcl::Fr*, int);