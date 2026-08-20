#pragma once

#include "flat_base.h"
#include "flat_chain.h"
#include "fan.h"
#include "alpha.h"
#include "scratch_mem.h"

class sec_kmeans;
class flat_params;

struct chaq_vecs
	{
	p_uint16_t nens = 0;
	p_uint16_t rens = 0;
	p_uint16_t pens = 0;
	p_uint16_t mens = 0;
	p_sid_t nensids = 0;
	p_sid_t rensids = 0;
	p_sid_t pensids = 0;
	p_sid_t mensids = 0;
	p_uint8_t sec32_codeseq = 0;
	};

// To keep a vector, copy pointer and set it 
// to zero in chaq_vecs2 object
struct chaq_vecs2
	{
	p_uint16_t nens = 0;
	p_uint16_t rens = 0;
	p_uint16_t pens = 0;
	p_uint16_t mens = 0;
	p_sid_t nensids = 0;
	p_sid_t rensids = 0;
	p_sid_t pensids = 0;
	p_sid_t mensids = 0;
	p_uint8_t pm2_codeseq = 0;
	p_uint8_t sec32_codeseq = 0;
	uint maxL = 0;
	};

/***
Chain quantizer / quantifier
All member functions are static.
Function arguments are pointers to data, not to flat_base objects.
Consuming code responsible for memory ownership and lifefime.
***/
class chaq
	{
private:
	chaq() = delete;

public:
	static uint8_t m_aacode2aa4code[20];

public:
	static void fill_distmx(
		cp_ic_t xyz,
		uint L,
		sid_t *distmx);

	static void fill_distmx(
		const flat_chain_t *chain,
		sid_t *distmx);

	static void fill_nen_vecs(
		cp_sid_t distmx,
		uint L,
		p_uint16_t nnvec,
		p_uint16_t nnsidvec);

	static void fill_nen_ren_vecs(
		cp_uint16_t pens,
		cp_uint16_t mens,
		cp_sid_t pensids,
		cp_sid_t mensids,
		uint L,
		p_uint16_t nens,
		p_uint16_t rens,
		p_sid_t nensids,
		p_sid_t rensids);

	static void fill_pen_men_vecs(
		cp_sid_t distmx,
		uint L,
		p_uint16_t penvec,
		p_uint16_t pensidvec,
		p_uint16_t menvec,
		p_uint16_t mensidvec);

	static void fill_fen_vecs(
		cp_sid_t distmx,
		uint L,
		p_uint16_t fenvec,
		p_uint16_t fensidvec);

	static uint8_t get_ss3(cp_sid_t distmx, uint L, uint pos);
	static uint8_t get_ss4(cp_sid_t distmx, uint L, uint pos);
	static void get_ss4_str(cp_sid_t distmx, uint L, string &ss);

	static void get_ss3_codeseq(cp_sid_t distmx, uint L, p_uint8_t codeseq);
	static void get_ss4_codeseq(cp_sid_t distmx, uint L, p_uint8_t codeseq);

	static void get_aa3_codeseq(const char *aacharseq, uint L, p_uint8_t codeseq);
	static void get_aa4_codeseq(const char *aacharseq, uint L, p_uint8_t codeseq);
	static void get_aa20_codeseq(const flat_chain_t *chain, p_uint8_t codeseq);

	static void get_pm_codeseq(cp_sid_t pensids, cp_sid_t mensids, uint L, p_uint8_t codeseq);

	static void get_packing_values(cp_sid_t distmx, uint L,
		uint maxsid, bool include_plus, bool include_minus,
		p_uint16_t values);

	static void get_turnd_values(cp_sid_t distmx, uint L,
		uint16_t undef_value, p_uint16_t values);

	static void get_packing_codeseq(cp_sid_t distmx, uint L, 
		uint maxsid, bool include_plus, bool include_minus, p_uint8_t codeseq);

	static void slow_get_values(
		const flat_chain_t *chain,
		FAN fan,
		uint alpha_size,
		p_uint16_t values);

	static void slow_get_codeseq(
		const flat_params &params,
		const flat_chain_t *chain,
		FAN fan,
		uint alpha_size,
		p_uint8_t codeseq);

	static void slow_get_codeseq_binned(
		const flat_params &params,
		const flat_chain_t *chain,
		FAN fan,
		uint alpha_size,
		p_uint8_t codeseq);

	static void slow_get_codeseq_discrete(
		const flat_params &params,
		const flat_chain_t *chain,
		FAN fan,
		uint alpha_size,
		p_uint8_t codeseq);

	static void slow_get_charseq_binned(
		const flat_params &params,
		const flat_chain_t *chain,
		FAN fan,
		uint8_t alpha_size,
		cp_uint16_t thresholds,
		uint16_t undef_value,
		char *charseq);

	static void slow_get_charseq_discrete(
		const flat_params &params,
		const flat_chain_t *chain,
		FAN fan,
		uint8_t alpha_size,
		uint8_t undef_code,
		char *charseq);

	static void get_sec_codeseq(
		uint alpha_size,
		cp_sid_t distmx,
		uint L,
		p_uint8_t codeseq);

	static const uint8_t *get_char2letter(uint alpha_size)
		{
		return (alpha_size == 20 ? g_CharToLetterAmino : g_CharToLetterMu);
		}

	static const unsigned char *get_letter2char(uint alpha_size)
		{
		return (alpha_size == 20 ? g_LetterToCharAmino : g_LetterToCharMu);
		}

	static void codeseq2charseq(
		const uint8_t *codeseq,
		uint L,
		uint alpha_size,
		char *charseq);

	static void charseq2codeseq(
		const char *charseq,
		uint L,
		uint alpha_size,
		uint8_t *codeseq);

	static void slow_get_angle_values(
		const flat_chain_t *chain,
		uint n,
		uint alpha_size,
		p_uint16_t values);

	static cp_uint16_t get_thresholds(
		const flat_params &params, FAN fan, uint alpha_size);
	static uint16_t get_undef_value(
		const flat_params &params, FAN fan, uint alpha_size);
	static uint8_t get_undef_code(
		FAN fan, uint alpha_size);

	static void set_aagroups(const string &aagroups);

	static size_t get_fill_chaq_vecs_bytes_per_pos();
	static size_t get_fast_get_codeseq_scratch_bytes_per_pos();

	// To keep a vector, copy pointer and set it 
	// to zero in chaq_vecs2 object
	static void alloc_chaq_vecs2(chaq_vecs2 &cv, uint L);
	static void free_chaq_vecs2(chaq_vecs2 &cv);

	static void fill_mega_prof(
		const flat_params &params,
		const flat_chain_t *chain,
		const sid_t *distmx,
		uint8_t *mega_prof,
		chaq_vecs2 *scratch_cv,
		uint8_t *scratch_buffer,
		uint scratch_buffer_bytes);

	static uint8_t *make_mega_prof(
		const flat_params &params,
		const flat_chain_t *chain,
		const sid_t *distmx,
		chaq_vecs2 *scratch_cv,
		uint8_t *scratch_buffer,
		uint scratch_buffer_bytes);

	static void fill_chaq_vecs2(
		cp_sid_t distmx,
		uint L,
		chaq_vecs2 &cv);

	static void fill_chaq_vecs_scratch_mem(
		cp_sid_t distmx,
		uint L,
		chaq_vecs &cv,
		scratch_mem &mem);

	static void fast_get_codeseq(
		const flat_params &params,
		const flat_chain_t *chain,
		const sid_t *distmx,
		const chaq_vecs2 *cv,
		FAN fan,
		uint alpha_size,
		p_uint8_t codeseq,
		uint8_t *scratch_buffer,
		uint scratch_buffer_bytes);

	static void fast_get_values(
		const flat_params &params,
		const sid_t *distmx,
		const chaq_vecs2 *cv,
		const flat_chain_t *chain,
		FAN fan,
		uint alpha_size,
		p_uint16_t values);

	static void sec32_codeseq_to_sec4(
		const uint8_t *codeseq_sec32, uint L,
		uint8_t *codeseq_sec4);

	static void codeseq_nu_to_kappa(
		const uint8_t *codeseq_nu, uint L,
		uint8_t *codeseq_kappa, size_t codeseq_kappa_bytes);

	static void codeseq_nu_to_kappa_inplace(
		uint8_t *codeseq, uint L);

	static void fill_codeseq_nu(
		const char *charseq_aa20,
		const uint8_t *codeseq_pm2,
		const uint8_t *codeseq_sec32,
		const uint L,
		uint8_t *codeseq_nu,
		size_t codeseq_nu_bytes);

	static void fill_codeseq_nu_from_chain(
		const flat_chain_t *chain,
		sid_t *distmx_buffer,
		chaq_vecs2 *cv_buffer,
		uint8_t *codeseq_nu,
		uint buffer_L);

	};
