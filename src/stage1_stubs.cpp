#include "myutils.h"
#include "flat_base.h"
#include "flat_dist_types.h"
#include "flat_helpers.h"
#include "chaq.h"
#include "alpha.h"

// Stage-1 host stubs / adapters for symbols that live in reseek-only units.
// alpha.cpp is not in muscle.vcxproj (alpha2.cpp owns SetAlpha), so Mu tables live here.

unsigned char g_LetterToCharMu[256] =
	{
	'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
	'L', 'K', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
	'U', 'V', 'W', 'X', 'Y', 'Z',
	'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'
	};

unsigned char g_CharToLetterMu[256];

static struct Init_g_CharToLetterMu
	{
	Init_g_CharToLetterMu()
		{
		memset(g_CharToLetterMu, 0xff, 256);
		for (uint i = 0; i < 36; ++i)
			g_CharToLetterMu[(unsigned char) g_LetterToCharMu[i]] = (unsigned char) i;
		}
	} s_Init_g_CharToLetterMu;

#if TRACK_ACTIVE
atomic<int64_t> g_flat_creates[FE_N];
atomic<int64_t> g_flat_destroys[FE_N];
atomic<int64_t> g_flat_bytes[FE_N];
#endif

#if TRACK_SRC
list<void *> g_flat_obj_list;
mutex g_flat_obj_list_lock;
#endif

void log_flat_stats(const string &msg)
	{
#if TRACK_ACTIVE
	Log("log_flat_stats(%s)\n", msg.c_str());
#else
	(void) msg;
#endif
	}

const ic_t sid2ic[65536] = {};

uint8_t s_sec32_to_sec4[32] = {
	0, 0, 1, 1, 2, 2, 1, 3, 1, 2, 2, 2, 2, 1, 2, 3,
	3, 2, 2, 2, 3, 3, 3, 3, 3, 1, 2, 3, 2, 2, 2, 2
	};

uint8_t g_nucode_to_kappacode[256] = {
	 0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,
	 8,  9, 10, 11, 12, 13, 14, 15,  8,  9, 10, 11, 12, 13, 14, 15,
	16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
	 8,  9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31,
	 8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
	16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
	16, 17, 18, 19, 20, 21, 22, 23,  8,  9, 10, 11, 12, 13, 14, 15,
	16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
	24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23,
	16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
	24, 25, 26, 27, 28, 29, 30, 31, 24, 25, 26, 27, 28, 29, 30, 31,
	24, 25, 26, 27, 28, 29, 30, 31, 24, 25, 26, 27, 28, 29, 30, 31,
	24, 25, 26, 27, 28, 29, 30, 31,  8,  9, 10, 11, 12, 13, 14, 15,
	16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
	16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
	16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
	};

void chaq::sec32_codeseq_to_sec4(
	const uint8_t *codeseq_sec32, uint L,
	uint8_t *codeseq_sec4)
	{
	for (uint i = 0; i < L; ++i)
		{
		uint8_t sec32code = codeseq_sec32[i];
		assert(sec32code < 32);
		codeseq_sec4[i] = s_sec32_to_sec4[sec32code];
		}
	}

void chaq::codeseq_nu_to_kappa_inplace(uint8_t *codeseq, uint L)
	{
	for (uint i = 0; i < L; ++i)
		{
		uint8_t nu_code = codeseq[i];
		uint8_t kappa_code = g_nucode_to_kappacode[nu_code];
		assert(kappa_code < 32);
		codeseq[i] = kappa_code;
		}
	}

void chaq::codeseq_nu_to_kappa(
	const uint8_t *codeseq_nu, uint L,
	uint8_t *codeseq_kappa, size_t codeseq_kappa_bytes)
	{
	asserta(L <= codeseq_kappa_bytes);
	for (uint i = 0; i < L; ++i)
		{
		uint8_t nu_code = codeseq_nu[i];
		uint8_t kappa_code = g_nucode_to_kappacode[nu_code];
		assert(kappa_code < 32);
		codeseq_kappa[i] = kappa_code;
		}
	}

void chaq::fill_codeseq_nu(
	const char *charseq_aa20,
	const uint8_t *codeseq_pm2,
	const uint8_t *codeseq_sec32,
	const uint L,
	uint8_t *codeseq_nu,
	size_t codeseq_nu_bytes)
	{
	asserta(L <= codeseq_nu_bytes);
	for (uint32_t pos = 0; pos < L; ++pos)
		{
		char c = charseq_aa20[pos];
		uint8_t code_aa20 = g_CharToLetterAmino[(unsigned char) c];
		if (code_aa20 >= 20) code_aa20 = 0;

		const uint8_t code_pm2 = codeseq_pm2[pos];
		const uint8_t code_sec32 = codeseq_sec32[pos];

		assert(code_aa20 < 20);
		assert(code_pm2 < 2);
		assert(code_sec32 < 32);

		const uint8_t code_aa4 = chaq::m_aacode2aa4code[code_aa20];
		const uint8_t code_nu = uint8_t(code_aa4 + 4*code_pm2 + 4*2*code_sec32);
		codeseq_nu[pos] = code_nu;
		}
	}

void chaq::fill_codeseq_nu_from_chain(
	const flat_chain_t *chain,
	sid_t *distmx_buffer,
	chaq_vecs2 *cv_buffer,
	uint8_t *codeseq_nu,
	uint buffer_L)
	{
	const uint L = chain->get_length();
	asserta(L <= buffer_L);

	sid_t *distmx = distmx_buffer;
	chaq_vecs2 *cv = cv_buffer;
	const char *charseq_aa20 = chain->m_aa->m_data;

	chaq::fill_distmx(chain, distmx);
	chaq::fill_chaq_vecs2(distmx, L, *cv);
	chaq::fill_codeseq_nu(
		charseq_aa20, cv->pm2_codeseq, cv->sec32_codeseq,
		L, codeseq_nu, buffer_L);
	}

uint32_t get_flat_pssm_feature_block_offsets(
	const uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	uint32_t * __restrict feature_block_offsets)
	{
	uint32_t offset = 0;
	for (uint32_t fi = 0; fi < nfeat; ++fi)
		{
		feature_block_offsets[fi] = offset;
		offset += alpha_sizes[fi];
		}
	return offset;
	}

void trunc_label(const string &Label, string &TruncatedLabel)
	{
	TruncatedLabel = Label;
	size_t n = TruncatedLabel.find(' ');
	if (n != string::npos)
		TruncatedLabel.resize(n);
	n = TruncatedLabel.find('|');
	if (n != string::npos)
		TruncatedLabel.resize(n);
	n = TruncatedLabel.find('/');
	if (n != string::npos)
		TruncatedLabel.resize(n);
	}

void trunc_label(string &Label)
	{
	string TruncatedLabel;
	trunc_label(Label, TruncatedLabel);
	Label = TruncatedLabel;
	}

void fill_pattern_offsets(const string &Str, uint8_t *offsets)
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(Str); ++i)
		{
		char c = Str[i];
		asserta(c == '0' || c == '1');
		if (c == '1')
			offsets[n++] = (uint8_t) i;
		}
	}

uint get_nr_pattern_ones(const string &Str)
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(Str); ++i)
		{
		char c = Str[i];
		asserta(c == '0' || c == '1');
		if (c == '1')
			++n;
		}
	return n;
	}

int kappa_max_pos_logodds()
	{
	return 0;
	}

// g_alpha_collect_lines is defined in alpha_collect.cpp
