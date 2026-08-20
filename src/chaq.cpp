#include "myutils.h"
#include "flat_base.h"
#include "flat_params.h"
#include "flat_distmx.h"
#include "sec_kmeans.h"
#include "scratch_mem.h"
#include "quantize.h"
#include "abcxyz.h"
#include "alpha.h"
#include "chaq.h"

/***
Groups sweep by $src/reseek_tune2/bash/aa4s_mix.bash
     0      1       2 3
ACHPST-DEKNQR-FILMVWY-G
***/
uint8_t chaq::m_aacode2aa4code[20] =
	{
	0,	//'A'
	0,  //'C'
	1,  //'D'
	1,  //'E'
	2,  //'F'
	3,  //'G'
	0,  //'H'
	2,  //'I'
	1,  //'K'
	2,  //'L'
	2,  //'M'
	1,  //'N'
	0,  //'P'
	1,  //'Q'
	1,  //'R'
	0,  //'S'
	0,  //'T'
	2,  //'V'
	2,  //'W'
	2,  //'Y'
	};

static const uint16_t s_packing_maxsid = dist2sid(15.0f);

static inline uint16_t radians_to_uint16(float theta)
	{
	const float MY_PI = 3.1415926535f; // M_PI not in MSVC <cmath>?
    // clamp just in case of small FP drift
    if (theta < 0.0f)
        theta = 0.0f;
    else if (theta > float(MY_PI))
        theta = float(MY_PI);

    // normalize to [0,1]
    float u = theta * (1.0f / float(MY_PI));

    // map to [0,65535] with 0 centered at ~32767
    float val = (2.0f * u) * 65535.0f * 0.5f; // same as u * 65535

    return uint16_t(val + 0.5f);
	}

static uint8_t get_packing_code(uint n)
	{
	return n;
	}

uint8_t chaq::get_undef_code(FAN fan, uint alpha_size)
	{
	switch (fan)
		{
	case FAN_sec:
	case FAN_nensec:
	case FAN_rensec:
	case FAN_pensec:
	case FAN_mensec:
	// lowest frequency code
		return alpha_size - 1;

	//case FAN_pack:
	//case FAN_mpack:
	//case FAN_ppack:
	//case FAN_angle:
	//// roughly median value, TODO?
	//	return alpha_size/2;

	case FAN_aa:
	case FAN_pm:
	case FAN_nu:
	case FAN_kappa:
		return 0;
		}

	// roughly median value, TODO?
	if (is_quantized(fan))
		return alpha_size/2;
	Die("get_undef_code(%s)", FAN2str(fan));
	return 0;
	}

static uint8_t get_aa4code(char c)
	{
	uint8_t aacode = g_CharToLetterAmino[c];
	if (aacode >= 20) return 0;
	return chaq::m_aacode2aa4code[aacode];
	}

static uint8_t get_aa3code(char c)
	{
	c = toupper(c);
	if (c == 'G')
		return 0;
	if (strchr("CFILMVWY", c) != 0)
		return 1;
	return 2;
	}

void chaq::fill_distmx(
	cp_ic_t xyz,
	uint L,
	sid_t *distmx)
	{
	fill_flat_distmx(xyz, L, distmx);
	}

void chaq::fill_distmx(
	const flat_chain_t *chain,
	sid_t *distmx)
	{
	fill_flat_distmx(chain->m_xyz->m_data,
		chain->get_length(), distmx);
	}

// 0=helix 1=strand 2=turn 3=loop
// Method from sec_str() in TMalign.cpp Zhang & Skolnick 2005
uint8_t chaq::get_ss4(const sid_t *distmx, uint L, uint pos)
	{
	if (pos < 2 || pos + 2 >= L)
		return 3;

	float dis13 = sid2dist(distmx[banded_i_lt_j_to_k(pos-2, pos)]);
	float dis14 = sid2dist(distmx[banded_i_lt_j_to_k(pos-2, pos+1)]);
	float dis15 = sid2dist(distmx[banded_i_lt_j_to_k(pos-2, pos+2)]);
	float dis24 = sid2dist(distmx[banded_i_lt_j_to_k(pos-1, pos+1)]);
	float dis25 = sid2dist(distmx[banded_i_lt_j_to_k(pos-1, pos+2)]);
	float dis35 = sid2dist(distmx[banded_i_lt_j_to_k(pos, pos+2)]);

	const float DH = 2.1f;
	if (fabs(dis15 - 6.37f) < DH && fabs(dis14 - 5.18f) < DH &&
		fabs(dis25 - 5.18f) < DH && fabs(dis13 - 5.45f) < DH &&
		fabs(dis24 - 5.45f) < DH && fabs(dis35 - 5.45f) < DH)
		return 0;	// helix

	const float DS = 1.42f;
	if (fabs(dis15 - 13.0f) < DS && fabs(dis14 - 10.4f) < DS &&
		fabs(dis25 - 10.4f) < DS && fabs(dis13 - 6.1f) < DS &&
		fabs(dis24 - 6.1f) < DS && fabs(dis35 - 6.1f) < DS)
		return 1;	// strand

	if (dis15 < 8.2)
		return 2;	// turn

	return 3; // loop
	}

// 0=helix 1=strand 2=other
uint8_t chaq::get_ss3(const sid_t * __restrict distmx, uint L, uint pos)
	{
	uint8_t ss4 = get_ss4(distmx, L, pos);
	return ss4 <= 2 ? ss4 : 2;
	}

void chaq::get_ss4_str(const sid_t * __restrict distmx, uint L, string &ss)
	{
	ss.clear();
	ss.reserve(L);
	for (uint pos = 0; pos < L; ++pos)
		{
		uint8_t letter = get_ss4(distmx, L, pos);
		assert(letter < 4);
		ss += "hst~"[letter];
		}
	}

void chaq::get_ss4_codeseq(cp_sid_t distmx, uint L, p_uint8_t intseq)
	{
	for (uint pos = 0; pos < L; ++pos)
		intseq[pos] = get_ss4(distmx, L, pos);
	}

void chaq::fill_nen_vecs(
	cp_sid_t distmx,
	uint L,
	p_uint16_t nen,
	p_uint16_t nensid)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	for (uint i = 0; i < L; ++i)
		{
		nen[i] = UINT16_MAX;
		nensid[i] = UINT16_MAX;
		}

	for (uint i = 0; i < L; ++i)
		{
		uint dmax = L - 1 - i;
		if (dmax > M)
			dmax = M;

		for (uint d = m; d <= dmax; ++d)
			{
			uint j = i + d;
			sid_t sid = distmx[banded_ij_to_k(i, j)];

			if (sid < nensid[i])
				{
				nensid[i] = sid;
				nen[i] = uint16_t(j);
				}

			if (sid < nensid[j])
				{
				nensid[j] = sid;
				nen[j] = uint16_t(i);
				}
			}
		}
	}

void chaq::fill_fen_vecs(
	cp_sid_t distmx,
	uint L,
	p_uint16_t fen,
	p_uint16_t fensid)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;
	for (uint i = 0; i < L; ++i)
		{
		fen[i] = UINT16_MAX;
		fensid[i] = 0;
		}

	for (uint i = 0; i < L; ++i)
		{
		uint dmax = L - 1 - i;
		if (dmax > M)
			dmax = M;

		for (uint d = m; d <= dmax; ++d)
			{
			uint j = i + d;
			sid_t sid = distmx[banded_ij_to_k(i, j)];

			if (sid > fensid[i])
				{
				fensid[i] = sid;
				fen[i] = uint16_t(j);
				}

			if (sid > fensid[j])
				{
				fensid[j] = sid;
				fen[j] = uint16_t(i);
				}
			}
		}
	}

void chaq::fill_nen_ren_vecs(
	cp_uint16_t pens,
	cp_uint16_t mens,
	cp_sid_t pensids,
	cp_sid_t mensids,
	uint L,
	p_uint16_t nens,
	p_uint16_t rens,
	p_sid_t nensids,
	p_sid_t rensids)
	{
	for (uint i = 0; i < L; ++i)
		{
		sid_t pensid = pensids[i];
		sid_t mensid = mensids[i];
		if (pensid <= mensid)
			{
			nens[i] = pens[i];
			rens[i] = mens[i];
			nensids[i] = pensids[i];
			rensids[i] = mensids[i];
			}
		else
			{
			nens[i] = mens[i];
			rens[i] = pens[i];
			nensids[i] = mensids[i];
			rensids[i] = pensids[i];
			}
		}
	}

void chaq::fill_pen_men_vecs(
	cp_sid_t distmx,
	uint L,
	p_uint16_t pen,
	p_uint16_t pensid,
	p_uint16_t men,
	p_uint16_t mensid)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	for (uint i = 0; i < L; ++i)
		{
		pen[i]  = UINT16_MAX;
		pensid[i] = UINT16_MAX;
		men[i]  = UINT16_MAX;
		mensid[i] = UINT16_MAX;
		}

	for (uint i = 0; i < L; ++i)
		{
		uint dmax = min(L - 1 - i, M);
		for (uint d = m; d <= dmax; ++d)
			{
			uint j = i + d;
			sid_t sid = distmx[banded_ij_to_k(i, j)];

			// forward for i: j > i
			if (sid < pensid[i])
				{
				pensid[i] = sid;
				pen[i] = uint16_t(j);
				}

			// reverse for j: i < j
			if (sid < mensid[j])
				{
				mensid[j] = sid;
				men[j] = uint16_t(i);
				}
			}
		}
	}

void chaq::get_aa3_codeseq(const char *aacharseq, uint L, p_uint8_t codeseq)
	{
	for (uint i = 0; i < L; ++i)
		codeseq[i] = get_aa3code(aacharseq[i]);
	}

void chaq::get_aa4_codeseq(const char *aacharseq, uint L, p_uint8_t codeseq)
	{
	for (uint i = 0; i < L; ++i)
		codeseq[i] = get_aa4code(aacharseq[i]);
	}

void chaq::get_aa20_codeseq(const flat_chain_t *chain, p_uint8_t codeseq)
	{
	uint8_t undef_code = get_undef_code(FAN_aa, 20);
	const char *charseq = chain->m_aa->m_data;
	const uint L = chain->get_length();
	for (uint i = 0; i < L; ++i)
		{
		char c = charseq[i];
		uint8_t code = g_CharToLetterAmino[c];
		if (code == 0xff)
			code = 0;
		assert(code < 20);
		codeseq[i] = code;
		}
	}

void chaq::get_pm_codeseq(cp_sid_t pensids, cp_sid_t mensids, uint L, p_uint8_t codeseq)
	{
	for (uint i = 0; i < L; ++i)
		codeseq[i] = (pensids[i] <= mensids[i] ? 0 : 1);
	}

void chaq::get_packing_codeseq(cp_sid_t distmx, uint L,
	uint maxsid, bool include_plus, bool include_minus, p_uint8_t codeseq)
	{
	Die("TODO");
	}

void chaq::get_turnd_values(cp_sid_t distmx, uint L,
	uint16_t undef_value, p_uint16_t values)
	{
	const uint w = flat_params::m_turnd_w;

	// Avoid unsigned underflow of L-w when L < w (and OOB writes when L < w).
	if (L <= 2*w)
		{
		for (uint pos = 0; pos < L; ++pos)
			values[pos] = undef_value;
		return;
		}

	for (uint pos = 0; pos < w; ++pos)
		values[pos] = undef_value;

	for (uint pos = w; pos < L-w; ++pos)
		{
		sid_t sid = distmx[banded_ij_to_k(pos-w, pos+w)];
		values[pos] = sid;
		}

	for (uint pos = L-w; pos < L; ++pos)
		values[pos] = undef_value;
	}

void chaq::get_packing_values(cp_sid_t distmx, uint L,
	uint maxsid, bool include_plus, bool include_minus,
	p_uint16_t values)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	for (uint i = 0; i < L; ++i)
		{
		uint32_t n = 0;
		if (include_minus)
			{
			int jmin = int(i) - int(M);
			if (jmin < 0) jmin = 0;
			for (uint j = jmin; j < i; ++j)
				{
				sid_t sid = distmx[banded_ij_to_k(i, j)];
				if (sid <= maxsid)
					++n;
				}
			}

		if (include_plus)
			{
			uint jmax = i + M;
			if (jmax >= L) jmax = L - 1;
			for (uint j = i+1; j < jmax; ++j)
				{
				sid_t sid = distmx[banded_ij_to_k(i, j)];
				if (sid <= maxsid)
					++n;
				}
			}

		values[i] = n;
		}
	}

void chaq::get_sec_codeseq(
	uint alpha_size,
	cp_sid_t distmx,
	uint L,
	p_uint8_t codeseq)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	const sec_kmeans *SK = sec_kmeans::get_SK(alpha_size, M);
	assert(SK->m_K == alpha_size);
	SK->get_codeseq(distmx, L, codeseq);
	}


void chaq::slow_get_codeseq_discrete(
	const flat_params &params,
	const flat_chain_t *chain,
	FAN fan,
	uint alpha_size,
	p_uint8_t codeseq)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	const uint L = chain->get_length();
	if (fan == FAN_aa)
		{
		uint8_t undef_code = get_undef_code(fan, alpha_size);
		switch (alpha_size)
			{
		case 3:	get_aa3_codeseq(chain->m_aa->m_data, L, codeseq); return;
		case 4:	get_aa4_codeseq(chain->m_aa->m_data, L, codeseq); return;
		case 20:
			const char *charseq = chain->m_aa->m_data;
			for (uint i = 0; i < L; ++i)
				{
				char c = charseq[i];
				uint8_t code = g_CharToLetterAmino[c];
				if (code == 0xff)
					code = undef_code;
				asserta(code < alpha_size);
				codeseq[i] = code;
				}
			return;
			}
		Die("chaq::slow_get_codeseq_discrete(aa, alpha_size=%u)", alpha_size);
		}

	uint16_t *distmx = myalloc(sid_t, L*M);
	uint16_t *pens = myalloc(uint16_t, L);
	uint16_t *mens = myalloc(uint16_t, L);
	uint16_t *nens = myalloc(uint16_t, L);
	uint16_t *rens = myalloc(uint16_t, L);
	sid_t *pensids = myalloc(sid_t, L);
	sid_t *mensids = myalloc(sid_t, L);
	sid_t *nensids = myalloc(sid_t, L);
	sid_t *rensids = myalloc(sid_t, L);
	uint8_t *sec_codeseq = myalloc(uint8_t, L);

	chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);

	chaq::fill_pen_men_vecs(
		distmx, L, pens, pensids, mens, mensids);

	chaq::fill_nen_ren_vecs(
		pens, mens, pensids, mensids, L,
		nens, rens, nensids, rensids);

	bool need_sec_codeseq = false;
	switch (fan)
		{
	case FAN_sec:
	case FAN_nensec:
	case FAN_rensec:
	case FAN_pensec:
	case FAN_mensec:
	case FAN_nu:
	case FAN_kappa:
		need_sec_codeseq = true;
		}

	if (need_sec_codeseq)
		{
		chaq::get_sec_codeseq(alpha_size, distmx, L, sec_codeseq);
#if DEBUG
		{
		for (uint i = 0; i < L; ++i)
			asserta(sec_codeseq[i] < alpha_size);
		}
#endif
		}

	// Only "discrete" features, not *dist etc.
	const size_t bytes = L*sizeof(uint16_t);
	switch (fan)
		{
	case FAN_nu:
		{
		asserta(alpha_size == 256);
		const char *aacharseq = chain->m_aa->m_data;
		uint8_t *pm2_codeseq = myalloc(uint8_t, L);
		chaq::get_pm_codeseq(pensids, mensids, L, pm2_codeseq);
		for (uint32_t pos = 0; pos < L; ++pos)
			{
			uint8_t aa = aacharseq[pos];
			uint8_t code_aa20 = g_CharToLetterAmino[aa];
			if (code_aa20 >= 20) code_aa20 = 0;
			const uint8_t code_pm2 = pm2_codeseq[pos];
			const uint8_t code_sec32 = sec_codeseq[pos];

			assert(code_aa20 < 20);
			assert(code_pm2 < 2);
			assert(code_sec32 < 32);

			const uint8_t code_aa4 = chaq::m_aacode2aa4code[code_aa20];
			const uint8_t code_nu = uint8_t(code_aa4 + 4*code_pm2 + 4*2*code_sec32);
			assert(code_nu < 256);

			codeseq[pos] = code_nu;
			}
		myfree(pm2_codeseq);
		return;
		}

	case FAN_kappa:
		{
		extern uint8_t g_nucode_to_kappacode[256];
		asserta(alpha_size == 32);
		const char *aacharseq = chain->m_aa->m_data;
		uint8_t *pm2_codeseq = myalloc(uint8_t, L);
		chaq::get_pm_codeseq(pensids, mensids, L, pm2_codeseq);
		for (uint32_t pos = 0; pos < L; ++pos)
			{
			uint8_t aa = aacharseq[pos];
			uint8_t code_aa20 = g_CharToLetterAmino[aa];
			if (code_aa20 >= 20) code_aa20 = 0;
			const uint8_t code_pm2 = pm2_codeseq[pos];
			const uint8_t code_sec32 = sec_codeseq[pos];

			assert(code_aa20 < 20);
			assert(code_pm2 < 2);
			assert(code_sec32 < 32);

			const uint8_t code_aa4 = chaq::m_aacode2aa4code[code_aa20];
			const uint8_t code_nu = uint8_t(code_aa4 + 4*code_pm2 + 4*2*code_sec32);
			assert(code_nu < 256);
			const uint8_t code_kappa = g_nucode_to_kappacode[code_nu];
			assert(code_kappa < 32);
			codeseq[pos] = code_kappa;
			}
		myfree(pm2_codeseq);
		return;
		}

	case FAN_sec:
		if (alpha_size == 32)
			memcpy(codeseq, sec_codeseq, L);
		else if (alpha_size == 4)
			sec32_codeseq_to_sec4(sec_codeseq, L, codeseq);
		else
			Die("chaq::slow_get_codeseq(sec) alpha_size=%u", alpha_size);
		break;

	case FAN_nensec:
		{
		assert(need_sec_codeseq);
		uint8_t undef_code = get_undef_code(fan, alpha_size);
		for (uint i = 0; i < L; ++i)
			{
			uint16_t xen = nens[i];
			asserta(xen < L || xen == UINT16_MAX);
			uint8_t code = (xen == UINT16_MAX) ?
				undef_code : sec_codeseq[xen];
			asserta(undef_code < alpha_size);
			codeseq[i] = code;
			}
		break;
		}

	case FAN_rensec:
		{
		assert(need_sec_codeseq);
		uint8_t undef_code = get_undef_code(fan, alpha_size);
		for (uint i = 0; i < L; ++i)
			{
			uint16_t xen = rens[i];
			asserta(xen < L || xen == UINT16_MAX);
			uint8_t code = (xen == UINT16_MAX) ?
				undef_code : sec_codeseq[xen];
			asserta(undef_code < alpha_size);
			codeseq[i] = code;
			}
		break;
		}

	case FAN_pensec:
		{
		uint8_t undef_code = get_undef_code(fan, alpha_size);
		assert(need_sec_codeseq);
		for (uint i = 0; i < L; ++i)
			{
			uint16_t xen = pens[i];
			asserta(xen < L || xen == UINT16_MAX);
			uint8_t code = (xen == UINT16_MAX) ?
				undef_code : sec_codeseq[xen];
			asserta(undef_code < alpha_size);
			codeseq[i] = code;
			}
		break;
		}

	case FAN_mensec:
		{
		uint8_t undef_code = get_undef_code(fan, alpha_size);
		assert(need_sec_codeseq);
		for (uint i = 0; i < L; ++i)
			{
			uint16_t xen = mens[i];
			asserta(xen < L || xen == UINT16_MAX);
			uint8_t code = (xen == UINT16_MAX) ?
				undef_code : sec_codeseq[xen];
			asserta(undef_code < alpha_size);
			codeseq[i] = code;
			}
		break;
		}

	case FAN_pm:
		asserta(alpha_size == 2);
		chaq::get_pm_codeseq(pensids, mensids, L, codeseq);
		break;

	case FAN_pack:
	case FAN_ppack:
	case FAN_mpack:
		chaq::slow_get_codeseq_binned(
			params, chain, fan, alpha_size, codeseq);
		break;

	default:	Die("slow_get_codeseq_discrete(%s)", FAN2str(fan));
		}
#if DEBUG
	{
	for (uint i = 0; i < L; ++i)
		asserta(codeseq[i] < alpha_size);
	}
#endif

	myfree(distmx);
	myfree(pens);
	myfree(mens);
	myfree(nens);
	myfree(rens);
	myfree(pensids);
	myfree(mensids);
	myfree(nensids);
	myfree(rensids);
	myfree(sec_codeseq);
	}

void chaq::slow_get_values(
	const flat_chain_t *chain,
	FAN fan,
	uint alpha_size,
	p_uint16_t values)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	if (fan == FAN_angle)
		{
		const uint n = flat_params::m_angle_n;
		slow_get_angle_values(chain, n, alpha_size, values);
		return;
		}

	const uint L = chain->get_length();

	uint16_t *distmx = myalloc(sid_t, L*M);
	uint16_t *pens = myalloc(uint16_t, L);
	uint16_t *mens = myalloc(uint16_t, L);
	uint16_t *nens = myalloc(uint16_t, L);
	uint16_t *rens = myalloc(uint16_t, L);
	uint16_t *fens = myalloc(uint16_t, L);
	sid_t *pensids = myalloc(sid_t, L);
	sid_t *mensids = myalloc(sid_t, L);
	sid_t *nensids = myalloc(sid_t, L);
	sid_t *rensids = myalloc(sid_t, L);
	sid_t *fensids = myalloc(sid_t, L);

	chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);

	chaq::fill_pen_men_vecs(
		distmx, L, pens, pensids, mens, mensids);

	chaq::fill_nen_ren_vecs(
		pens, mens, pensids, mensids, L,
		nens, rens, nensids, rensids);

	chaq::fill_fen_vecs(
		distmx, L, fens, fensids);

	// Only "float" features, not aa, ss3 etc.
	const size_t bytes = L*sizeof(uint16_t);
	switch (fan)
		{
	case FAN_nendist:	memcpy(values, nensids, bytes); break;
	case FAN_rendist:	memcpy(values, rensids, bytes); break;
	case FAN_pendist:	memcpy(values, pensids, bytes); break;
	case FAN_mendist:	memcpy(values, mensids, bytes); break;
	case FAN_fendist:	memcpy(values, fensids, bytes); break;

	case FAN_pmdiff:
		{
		for (uint i = 0; i < L; ++i)
			{
			const sid_t ic_40A = 20400;//TODO param?
			sid_t pensid = pensids[i];
			sid_t mensid = mensids[i];
			ic_t penic = sid2ic[pensid];
			ic_t menic = sid2ic[mensid];
			sid_t pmdiff = UINT16_MAX;
			if (pensid == UINT16_MAX || mensid == UINT16_MAX)
				pmdiff = ic_40A; // pensid==mensid
			else
				{
				pensid = min(pensid, ic_40A);
				mensid = min(mensid, ic_40A);
				pmdiff = ic_40A + pensid - mensid;
				assert(pmdiff <= 2*ic_40A);
				}
			values[i] = pmdiff;
			}
		break;
		}

	case FAN_pmdd:
		{
		// # python /mnt/c/src/py/angstroms_to_ic_and_sid.py 20
		// 20.0 Angstroms = 10200 ic
		// Squared distance sid 2500

		for (uint i = 0; i < L; ++i)
			{
			const sid_t sid_20A = 2500;//TODO param?
			sid_t pensid = pensids[i];
			sid_t mensid = mensids[i];
			sid_t pmdd = UINT16_MAX;
			if (pensid == UINT16_MAX || mensid == UINT16_MAX)
				pmdd = sid_20A; // pensid==mensid
			else
				{
				pensid = min(pensid, sid_20A);
				mensid = min(mensid, sid_20A);
				pmdd = sid_20A + pensid - mensid;
				assert(pmdd <= 2*sid_20A);
				}
			values[i] = pmdd;
			}
		break;
		}

	case FAN_pm:
		{
		// Special-case hack, convert 8- to 16-bit.
		uint8_t *codeseq = myalloc(uint8_t, L);
		chaq::get_pm_codeseq(pensids, mensids, L, codeseq);
		for (uint i = 0; i < L; ++i)
			values[i] = codeseq[i];
		myfree(codeseq);
		break;
		}

	case FAN_pack:
		{
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
		static const uint16_t maxsid = dist2sid(15.0f);//TODO param for 15.0f
		chaq::get_packing_values(distmx, L, maxsid, true, true, values);
		break;
		}

	case FAN_ppack:
		{
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
		static const uint16_t maxsid = dist2sid(15.0f);//TODO param for 15.0f
		chaq::get_packing_values(distmx, L, maxsid, true, false, values);
		break;
		}

	case FAN_mpack:
		{
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
		static const uint16_t maxsid = dist2sid(15.0f);//TODO param for 15.0f
		chaq::get_packing_values(distmx, L, maxsid, false, true, values);
		break;
		}

	case FAN_turnd:
		{
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
		static const uint16_t w = 5; //TODO param for w=5
		static const uint16_t undef_value =  1540; // measured median
		chaq::get_turnd_values(distmx, L, undef_value, values);
		break;
		}

	default:	Die("slow_get_values(%s)", FAN2str(fan));
		}

	myfree(distmx);
	myfree(pens);
	myfree(mens);
	myfree(nens);
	myfree(rens);
	myfree(pensids);
	myfree(mensids);
	myfree(nensids);
	myfree(rensids);
	}

void chaq::codeseq2charseq(
	const uint8_t *codeseq,
	uint L,
	uint alpha_size,
	char *charseq)
	{
	const unsigned char *code2char = get_letter2char(alpha_size);
	for (uint i = 0; i < L; ++i)
		{
		uint8_t code = codeseq[i];
		assert(code < alpha_size);
		charseq[i] = code2char[code];
		}
	}

void chaq::charseq2codeseq(
	const char *charseq,
	uint L,
	uint alpha_size,
	uint8_t *codeseq)
	{
	const uint8_t *char2code = get_char2letter(alpha_size);
	for (uint i = 0; i < L; ++i)
		{
		char c = charseq[i];
		uint8_t code = char2code[c];
		assert(code < alpha_size);
		codeseq[i] = code;
		}
	}

void chaq::slow_get_charseq_discrete(
	const flat_params &params,
	const flat_chain_t *chain,
	FAN fan,
	uint8_t alpha_size,
	uint8_t undef_code,
	char *charseq)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	const uint L = chain->get_length();
	assert(L > 0);
	uint8_t *codeseq = myalloc(uint8_t, L);
	slow_get_codeseq_discrete(params, chain, fan, alpha_size, codeseq);
	codeseq2charseq(codeseq, L, alpha_size, charseq);
	myfree(codeseq);
	}

void chaq::slow_get_codeseq_binned(
	const flat_params &params,
	const flat_chain_t *chain,
	FAN fan,
	uint alpha_size,
	p_uint8_t codeseq)
	{
	const uint L = chain->get_length();

	cp_uint16_t thresholds = get_thresholds(params, fan, alpha_size);
	const uint16_t undef_value = get_undef_value(params, fan, alpha_size);

	uint16_t *values = myalloc(uint16_t, L);
	chaq::slow_get_values(chain, fan, alpha_size, values);
	for (uint i = 0; i < L; ++i)
		{
		uint16_t value = values[i];
		if (value == UINT16_MAX)
			value = undef_value;
		uint8_t code = get_bin(value, alpha_size, thresholds);
		assert(code < alpha_size);
		codeseq[i] = code;
		}
	myfree(values);
	}

void chaq::slow_get_codeseq(
	const flat_params &params,
	const flat_chain_t *chain,
	FAN fan,
	uint alpha_size,
	p_uint8_t codeseq)
	{
	const uint L = chain->get_length();

	if (is_quantized(fan))
		chaq::slow_get_codeseq_binned(params, chain, fan, alpha_size, codeseq);
	else
		chaq::slow_get_codeseq_discrete(params, chain, fan, alpha_size, codeseq);
	}

void chaq::slow_get_charseq_binned(
	const flat_params &params,
	const flat_chain_t *chain,
	FAN fan,
	uint8_t alpha_size,
	cp_uint16_t thresholds,
	uint16_t undef_value,
	char *charseq)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	const uint L = chain->get_length();
	uint16_t *values = myalloc(uint16_t, L);
	chaq::slow_get_values(chain, fan, alpha_size, values);
	for (uint i = 0; i < L; ++i)
		{
		uint16_t value = values[i];
		if (value == UINT16_MAX)
			value = undef_value;
		uint8_t code = get_bin(value, alpha_size, thresholds);
		assert(code < alpha_size);
		charseq[i] = g_LetterToCharMu[code];
		}
	myfree(values);
	}

void chaq::slow_get_angle_values(
	const flat_chain_t *chain,
	uint n,
	uint alpha_size,
	p_uint16_t values)
	{
	const uint16_t undef_value = radians_to_uint16(0);

	const uint L = chain->get_length();

	// Avoid unsigned underflow of L-n-1 when L < n+1 (and OOB writes when L < n).
	if (L <= 2*n)
		{
		for (uint pos = 0; pos < L; ++pos)
			values[pos] = undef_value;
		return;
		}

	for (uint pos = 0; pos < n; ++pos)
		values[pos] = undef_value;

	const ic_t *xyz = chain->m_xyz->m_data;
	for (uint pos = n; pos < L - n - 1; ++pos)
		{
		uint posa = pos - n;
		uint posc = pos + n;

		ic_t xa = xyz[3*posa];
		ic_t ya = xyz[3*posa+1];
		ic_t za = xyz[3*posa+2];

		ic_t xb = xyz[3*pos];
		ic_t yb = xyz[3*pos+1];
		ic_t zb = xyz[3*pos+2];

		ic_t xc = xyz[3*posc];
		ic_t yc = xyz[3*posc+1];
		ic_t zc = xyz[3*posc+2];

		float theta = GetTheta3D_3pts<float>(
			xa, ya, za,
			xb, yb, zb,
			xc, yc, zc);
		values[pos] = radians_to_uint16(theta);
		}

	for (uint pos = L - n - 1; pos < L; ++pos)
		values[pos] = undef_value;
	}

void chaq::set_aagroups(const string &aagroups)
	{
	vector<string> flds;
	Split(aagroups, flds, '-');
	asserta(flds.size() == 4);

	for (uint i = 0; i < 4; ++i)
		m_aacode2aa4code[i] = 0xff;

	for (uint aa4code = 0; aa4code < 4; ++aa4code)
		{
		const string &group = flds[aa4code];
		for (auto c : group)
			{
			uint8_t aacode = g_CharToLetterAmino[c];
			asserta(aacode < 20);
			m_aacode2aa4code[aacode] = aa4code;
			}
		}

	for (uint i = 0; i < 4; ++i)
		asserta(m_aacode2aa4code[i] < 4);
	}

size_t chaq::get_fast_get_codeseq_scratch_bytes_per_pos()
	{
	return sizeof(uint16_t);	// values
	}

size_t chaq::get_fill_chaq_vecs_bytes_per_pos()
	{
	return 4*sizeof(uint16_t)	// Xenss
		+ 4*sizeof(sid_t)		// Xensids
		+ sizeof(uint8_t);		// sec32_codeseq
	}

void chaq::alloc_chaq_vecs2(chaq_vecs2 &cv, uint maxL)
	{
	cv.nens = myalloc(uint16_t, maxL);
	cv.rens = myalloc(uint16_t, maxL);
	cv.pens = myalloc(uint16_t, maxL);
	cv.mens = myalloc(uint16_t, maxL);
	cv.nensids = myalloc(sid_t, maxL);
	cv.rensids = myalloc(sid_t, maxL);
	cv.pensids = myalloc(sid_t, maxL);
	cv.mensids = myalloc(sid_t, maxL);
	cv.sec32_codeseq = myalloc(uint8_t, maxL);
	cv.pm2_codeseq = myalloc(uint8_t, maxL);
	cv.maxL = maxL;
	}

void chaq::free_chaq_vecs2(chaq_vecs2 &cv)
	{
	myfree(cv.nens);
	myfree(cv.rens);
	myfree(cv.pens);
	myfree(cv.mens);
	myfree(cv.nensids);
	myfree(cv.rensids);
	myfree(cv.pensids);
	myfree(cv.mensids);
	myfree(cv.sec32_codeseq);
	myfree(cv.pm2_codeseq);
	cv.maxL = 0;
	}

void chaq::fill_chaq_vecs_scratch_mem(
	cp_sid_t distmx,
	uint L,
	chaq_vecs &cv,
	scratch_mem &scratch)
	{
	cv.nens = scratch.get<uint16_t>(L);
	cv.rens = scratch.get<uint16_t>(L);
	cv.pens = scratch.get<uint16_t>(L);
	cv.mens = scratch.get<uint16_t>(L);
	cv.nensids = scratch.get<sid_t>(L);
	cv.rensids = scratch.get<sid_t>(L);
	cv.pensids = scratch.get<sid_t>(L);
	cv.mensids = scratch.get<sid_t>(L);
	cv.sec32_codeseq = scratch.get<uint8_t>(L);

	fill_pen_men_vecs(distmx, L,
		cv.pens, cv.pensids,
		cv.mens, cv.mensids);

	fill_nen_ren_vecs(
		cv.pens, cv.mens,
		cv.pensids, cv.mensids, L,
		cv.nens, cv.rens,
		cv.nensids, cv.rensids);

	get_sec_codeseq(32, distmx, L, cv.sec32_codeseq);
	}

void chaq::fill_chaq_vecs2(
	cp_sid_t distmx,
	uint L,
	chaq_vecs2 &cv)
	{
	asserta(L <= cv.maxL); // TODO=maxL
	fill_pen_men_vecs(distmx, L,
		cv.pens, cv.pensids,
		cv.mens, cv.mensids);

	fill_nen_ren_vecs(
		cv.pens, cv.mens,
		cv.pensids, cv.mensids, L,
		cv.nens, cv.rens,
		cv.nensids, cv.rensids);

	get_sec_codeseq(32, distmx, L, cv.sec32_codeseq);
	get_pm_codeseq(cv.pensids, cv.mensids, L, cv.pm2_codeseq);
	}

uint8_t *chaq::make_mega_prof(
	const flat_params &params,
	const flat_chain_t *chain,
	const sid_t *distmx,
	chaq_vecs2 *scratch_cv,
	uint8_t *scratch_buffer,
	uint scratch_buffer_bytes)
	{
	const uint L = chain->get_length();
	asserta(L > 0);
	const uint nfeat = params.m_nfeat;
	asserta(nfeat > 0);
	uint8_t *mega_prof = myalloc(uint8_t, nfeat*L);
	fill_mega_prof(params, chain, distmx, mega_prof,
	  scratch_cv, scratch_buffer, scratch_buffer_bytes);
	return mega_prof;
	}

void chaq::fill_mega_prof(
	const flat_params &params,
	const flat_chain_t *chain,
	const sid_t *distmx,
	uint8_t *mega_prof,
	chaq_vecs2 *scratch_cv,
	uint8_t *scratch_buffer,
	uint scratch_buffer_bytes)
	{
	const uint L = chain->get_length();
	asserta(L > 0);

	const uint nfeat = params.m_nfeat;
	asserta(nfeat > 0);
	asserta(mega_prof != 0);
	asserta(scratch_cv != 0);

	chaq::fill_chaq_vecs2(distmx, L, *scratch_cv);

#if DEBUG
	memset(mega_prof, 0xff, nfeat*L);
#endif

	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const FAN fan = params.m_fans[fi];
		const uint alpha_size = params.m_alpha_sizes[fi];

		uint8_t *codeseq = mega_prof + fi*L;
		uint8_t undef_code = chaq::get_undef_code(fan, alpha_size);
		(void) undef_code;
		chaq::fast_get_codeseq(
			params, chain, distmx, scratch_cv, fan, alpha_size,
			codeseq, scratch_buffer, scratch_buffer_bytes);

#if DEBUG
		for (uint pos = 0; pos < L; ++pos)
			assert(codeseq[pos] < alpha_size);
#endif
		}
	}
