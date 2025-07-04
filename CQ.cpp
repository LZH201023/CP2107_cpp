#include <mcl/bls12_381.hpp>
#include <unordered_map>
#include <ctime>
#include <CQ.hpp>
#include <NTT.hpp>
#include <KZG.hpp>

SRS gen(int N, mcl::Fr* t);
void IsInTable(mcl::G1& cm, mcl::Fr* t, int N, SRS& srs, int n, mcl::Fr* f);

static mcl::Fr evalPoly(const mcl::Fr* poly, int df, const mcl::Fr& i);
static mcl::G1* fast_KZG(mcl::Fr* T, int N, const mcl::Fr& ge, mcl::G1* srs);
static void poly_div_Z(const mcl::Fr* f, int df, int m, mcl::Fr* Q);

SRS gen(int N, mcl::Fr* t)
{
	clock_t timer = clock();
	mcl::Fr x;
	do
		x.setByCSPRNG();
	while (x == 0 || x == 1);
	mcl::Fr a = 1;
	mcl::G1* srs1 = new mcl::G1[N];
	mcl::G2* srs2 = new mcl::G2[N + 1];
	mcl::G1 g; mcl::Fr aux;
	mcl::Fp r; r.setByCSPRNG();
	mapToG1(g, r);
	aux.setByCSPRNG();
	mcl::G2 h; mcl::hashAndMapToG2(h, aux.getStr());
	for (int i = 0; i < N; i++)
	{
		srs1[i] = g * a;
		srs2[i] = h * a;
		a *= x;
	}
	srs2[N] = h * a;

	mcl::G2 Zx2 = srs2[N] - h;

	mcl::Fr ge = -1;
	for (int i = 1; i < log2(N); i++)
		mcl::Fr::squareRoot(ge, ge);

	mcl::Fr* T = iNTT<mcl::Fr>(t, N, ge);
	mcl::G2 Tx2 = h * evalPoly(T, N - 1, x);

	mcl::G1* qxs1 = fast_KZG(T, N, ge, srs1);

	mcl::G1* Lxs1 = NTT<mcl::G1>(srs1, N);
	int hN = N / 2; mcl::G1 temp;
	mcl::Fr N_inv = N; mcl::Fr::inv(N_inv, N_inv);
	for (int i = 1; i < hN; i++)
	{
		temp = Lxs1[i] * N_inv; 
		Lxs1[i] = Lxs1[N - i] * N_inv;
		Lxs1[N - i] = temp;
	}
	Lxs1[0] *= N_inv; Lxs1[hN] *= N_inv;
	mcl::G1* Lps1 = new mcl::G1[N];
	mcl::Fr gi = 1;
	for (int i = 0; i < N; i++)
	{
		Lps1[i] = Lxs1[i] * gi - srs1[N - 1] * N_inv;
		gi /= ge;
	}
	timer = clock() - timer;
	std::cout << "Setup time: " << timer * 1000.0 / CLOCKS_PER_SEC << "ms\n";
	return SRS(N, srs1, srs2, Zx2, Tx2, qxs1, Lxs1, Lps1);
}

void IsInTable(mcl::G1& cm, mcl::Fr* t, int N, SRS& srs, int n, mcl::Fr* f)
{
	clock_t timer = clock();
	clock_t P_time = 0, V_time = 0;
	//precompute table
	std::unordered_map<mcl::Fr, int> table;
	for (int i = 0; i < N; i++)
		table[t[i]] = i;
	timer = clock() - timer;

	std::cout << "Protocol begins\n\n";

	//Prover computes m(X) and sends commitment m
	timer = clock();
	std::unordered_map<int, int> m_table;
	for (int i = 0; i < n; i++)
		m_table[table[f[i]]]++;
	mcl::G1 m = mcl::getG1basePoint();
	mcl::Fr coeff;
	for (const auto& pair : m_table)
	{
		coeff = pair.second;
		m += srs.Lxs1[pair.first] * coeff;
	}
	P_time += clock() - timer;
	std::cout << "Prover's commitment m: " << m.getStr() << "\n";

	//Verifer sends challenge beta
	timer = clock();
	mcl::Fr beta; beta.setByCSPRNG();
	V_time += clock() - timer;
	std::cout << "Verifier: beta = " << beta << "\n";

	//Prover computes and commits to A(X)
	timer = clock();
	std::unordered_map<int, mcl::Fr> A;
	for (const auto& pair : m_table)
	{
		coeff = pair.second;
		A[pair.first] = coeff / (t[pair.first] + beta);
	}
	mcl::G1 a = mcl::getG1basePoint();
	for (const auto& pair : A)
	{
		coeff = pair.second;
		a += srs.Lxs1[pair.first] * coeff;
	}
	P_time += clock() - timer;
	std::cout << "Prover's commitment a: " << a.getStr() << "\n";

	//Prover computes and commits to Q_A(X)
	timer = clock();
	mcl::G1 qa = mcl::getG1basePoint();
	for (const auto& pair : A)
	{
		coeff = pair.second;
		qa += srs.qxs1[pair.first] * coeff;
	}
	P_time += clock() - timer;
	std::cout << "Prover's commitment qa: " << qa.getStr() << "\n";

	//Prover computes and commits to B0(X)
	timer = clock();
	mcl::Fr* B = new mcl::Fr[n];
	for (int i = 0; i < n; i++)
		B[i] = 1 / (f[i] + beta);
	mcl::G1 b0;
	mcl::Fr* Bx = iNTT<mcl::Fr>(B, n);
	mcl::G1::mulVec(b0, srs.srs1, Bx + 1, n - 1);
	P_time += clock() - timer;
	std::cout << "Prover's commitment b0: " << b0.getStr() << "\n";

	//Prover computes and commits to Q_B(X)
	timer = clock();
	mcl::Fr* fx = iNTT<mcl::Fr>(f, n);
	fx[0] += beta;
	mcl::Fr* pdt = poly_mult(Bx, n - 1, fx, n - 1);
	fx[0] -= beta;
	mcl::Fr* Q_B = new mcl::Fr[n - 1];
	poly_div_Z(pdt, 2 * n - 2, n, Q_B);
	mcl::G1 qb;
	mcl::G1::mulVec(qb, srs.srs1, Q_B, n - 1);
	P_time += clock() - timer;
	std::cout << "Prover's commitment qb: " << qb.getStr() << "\n";
	delete[] B;
	delete[] pdt;
	
	//Prover computes and commits to P(X)
	timer = clock();
	mcl::G1 p;
	mcl::G1::mulVec(p, srs.srs1 + N - n + 1, Bx + 1, n - 1);
	P_time += clock() - timer;
	std::cout << "Prover's commitment p: " << p.getStr() << "\n";

	//Verfier checks A and B0
	timer = clock();
	mcl::Fp12 e1, e2, e3;
	mcl::pairing(e1, a, srs.Tx2);
	mcl::pairing(e2, qa, srs.Zx2);
	mcl::pairing(e3, m - a * beta, srs.srs2[0]);
	bool check = (e1 == e2 * e3);
	V_time += clock() - timer;
	if (!check)
	{
		std::cout << "Verifier rejects\nProtocol ends\n\n";
		std::cout << "Prover time: " << P_time * 1000.0 / CLOCKS_PER_SEC << "ms\n"
			<< "Verifier time: " << V_time * 1000.0 / CLOCKS_PER_SEC << "ms\n";
		delete[] Bx; delete[] fx; delete[] Q_B;
		return;
	}

	//Verifer sends challenge gamma
	timer = clock();
	mcl::Fr gamma;
	gamma.setByCSPRNG();
	V_time += clock() - timer;
	std::cout << "Verifier: gamma = " << gamma << "\n";

	//Prover computes B0(gamma), f(gamma), A(0)
	timer = clock();
	mcl::Fr b0g, fg, a_0 = 0;
	b0g = evalPoly(Bx + 1, n - 2, gamma);
	fg = evalPoly(fx, n - 1, gamma);
	mcl::Fr Ni = N; mcl::Fr::inv(Ni, Ni);
	for (const auto& pair : A)
		a_0 += pair.second * Ni;
	P_time += clock() - timer;
	std::cout << "Prover: B0(gamma) = " << b0g
		<< "\nf(gamma) = " << fg
		<< "\nA(0) = " << a_0 << "\n";

	//Verifier's computations and challenege eta
	timer = clock();
	mcl::Fr ni = n; mcl::Fr::inv(ni, ni);
	mcl::Fr b_0 = N * a_0 * ni;
	mcl::Fr::pow(coeff, gamma, n);
	mcl::Fr Qbg = ((b0g * gamma + b_0) * (fg + beta) - 1) / (coeff - 1);
	mcl::Fr eta; eta.setByCSPRNG();
	mcl::Fr v_v = b0g + eta * fg + eta * eta * Qbg;
	V_time += clock() - timer;
	std::cout << "Verifier: eta = " << eta << "\n";

	//Prover computes batched KZH proof pi
	timer = clock();
	mcl::Fr::pow(coeff, gamma, n);
	mcl::Fr v_p = b0g + eta * fg + eta * eta * ((b0g * gamma + b_0) * (fg + beta) - 1) / (coeff - 1);
	coeff = eta * eta;

	for (int i = 0; i < n - 1; i++)
	{
		fx[i] *= eta;
		fx[i] += Bx[i + 1];
		fx[i] += coeff * Q_B[i];
	}
	fx[n - 1] *= eta;
	fx[0] -= v_p;

	mcl::Fr rem = fx[n - 1];
	for (int i = n - 2; i >= 0; i--) {
		mcl::Fr temp = fx[i];
		fx[i] = rem;
		rem = temp + rem * gamma;
	}

	mcl::G1 pi;
	mcl::G1::mulVec(pi, srs.srs1, fx, n - 1);
	P_time += clock() - timer;
	std::cout << "Prover's commitment pi: " << pi.getStr() << "\n";

	//Verifier's batched check
	timer = clock();
	mcl::G1 c = b0 + cm * eta + qb * eta * eta;
	mcl::pairing(e1, c - srs.srs1[0] * v_v + pi * gamma, srs.srs2[0]);
	mcl::pairing(e2, pi, srs.srs2[1]);
	check = (e1 == e2);
	V_time += clock() - timer;
	if (!check)
	{
		std::cout << "Verifier rejects\nProtocol ends\n\n";
		std::cout << "Prover time: " << P_time * 1000.0 / CLOCKS_PER_SEC << "ms\n"
			<< "Verifier time: " << V_time * 1000.0 / CLOCKS_PER_SEC << "ms\n";
		delete[] Bx; delete[] fx; delete[] Q_B;
		return;
	}

	//Prover computes commitment a0
	timer = clock();
	mcl::G1 a0 = mcl::getG1basePoint();
	for (const auto& pair : A)
		a0 += srs.Lps1[pair.first] * pair.second;
	P_time += clock() - timer;
	std::cout << "Prover's commitment a0: " << a0.getStr() << "\n";

	//Verifier checks a0
	timer = clock();
	mcl::pairing(e1, a - srs.srs1[0] * a_0, srs.srs2[0]);
	mcl::pairing(e2, a0, srs.srs2[1]);
	check = (e1 == e2);
	V_time += clock() - timer;
	if (!check)
	{
		std::cout << "Verifier rejects\nProtocol ends\n\n";
		std::cout << "Prover time: " << P_time * 1000.0 / CLOCKS_PER_SEC << "ms\n"
			<< "Verifier time: " << V_time * 1000.0 / CLOCKS_PER_SEC << "ms\n";
		delete[] Bx; delete[] fx; delete[] Q_B;
		return;
	}

	std::cout << "Verifier accepts\nProtocol ends\n\n";
	std::cout << "Prover time: " << P_time * 1000.0 / CLOCKS_PER_SEC << "ms\n"
		<< "Verifier time: " << V_time * 1000.0 / CLOCKS_PER_SEC << "ms\n";
	delete[] Bx; delete[] fx; delete[] Q_B;
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

int main()
{
	mcl::initPairing(mcl::BN_SNARK1);
	mcl::Fr* t = new mcl::Fr[128];
	for (int i = 0; i < 128; i++)
		t[i].setByCSPRNG();
	mcl::Fr* f = new mcl::Fr[32];
	for (int i = 0; i < 32; i++)
		f[i] = t[i];
	mcl::Fr* fx = iNTT(f, 32);
	SRS srs = gen(128, t);

	//Commitment to f
	mcl::G1 cm;
	mcl::G1::mulVec(cm, srs.srs1, fx, 32);

	IsInTable(cm, t, 128, srs, 32, f);
	delete[] t; delete[] f; delete[] fx;
}

static mcl::G1* fast_KZG(mcl::Fr* T, int N, const mcl::Fr& ge, mcl::G1* srs)
{

	int n = 2 * N;
	mcl::G1* s = new mcl::G1[n];
	for (int i = 0; i < N; i++)
		s[i] = srs[N - 1 - i];
	for (int i = N; i < n; i++)
		s[i] = mcl::getG1basePoint();
	s = NTT<mcl::G1>(s, n);
	mcl::Fr* v = new mcl::Fr[n];
	for (int i = 0; i <= N; i++)
		v[i] = 0;
	for (int i = N + 1; i < n; i++)
		v[i] = T[i - N];
	v = NTT<mcl::Fr>(v, n);
	for (int i = 0; i < n; i++)
		s[i] *= v[i];
	mcl::G1* h = iNTT<mcl::G1>(s, n);
	mcl::G1* proofs = NTT<mcl::G1>(h, N);
	mcl::Fr x = 1;
	mcl::Fr Ni = N; mcl::Fr::inv(Ni, Ni);
	for (int i = 0; i < N; i++)
	{
		proofs[i] *= x * Ni;
		x *= ge;
	}
	return proofs;
}

static void poly_div_Z(const mcl::Fr* f, int df, int m, mcl::Fr* Q)
{
	int n = df + 1;
	mcl::Fr* r = new mcl::Fr[n];
	std::copy(f, f + n, r);
	for (int i = n - m - 1; i >= 0; --i) {
		Q[i] = r[i + m];
		r[i] = r[i] + Q[i];
	}
}