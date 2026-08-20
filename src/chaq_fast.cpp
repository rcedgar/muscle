#include "myutils.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "chaq.h"
#include "quantize.h"

void chaq::fast_get_values(
	const flat_params &params,
	const sid_t *distmx,
	const chaq_vecs2 *cv,
	const flat_chain_t *chain,
	FAN fan,
	uint alpha_size,
	p_uint16_t values)
	{
	const uint L = chain->get_length();
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;
	const size_t bytes16 = L*sizeof(uint16_t);
	switch (fan)
		{
	case FAN_angle:
		{
		const uint n = flat_params::m_angle_n;
		slow_get_angle_values(chain, n, alpha_size, values);
		return;
		}

	case FAN_nendist:	memcpy(values, cv->nensids, bytes16); return;
	case FAN_rendist:	memcpy(values, cv->rensids, bytes16); return;
	case FAN_pendist:	memcpy(values, cv->pensids, bytes16); return;
	case FAN_mendist:	memcpy(values, cv->mensids, bytes16); return;

	case FAN_pmdiff:
		{
		for (uint i = 0; i < L; ++i)
			{
			const sid_t ic_40A = 20400;//TODO param?
			sid_t pensid = cv->pensids[i];
			sid_t mensid = cv->mensids[i];
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
		return;
		}

	case FAN_pmdd:
		{
		// # python /mnt/c/src/py/angstroms_to_ic_and_sid.py 20
		// 20.0 Angstroms = 10200 ic
		// Squared distance sid 2500

		for (uint i = 0; i < L; ++i)
			{
			const sid_t sid_20A = 2500;//TODO param?
			sid_t pensid = cv->pensids[i];
			sid_t mensid = cv->mensids[i];
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
		return;
		}

	case FAN_pack:
		{
		static const uint16_t maxsid = dist2sid(15.0f);//TODO param for 15.0f
		chaq::get_packing_values(distmx, L, maxsid, true, true, values);
		return;
		}

	case FAN_ppack:
		{
		static const uint16_t maxsid = dist2sid(15.0f);//TODO param for 15.0f
		chaq::get_packing_values(distmx, L, maxsid, true, false, values);
		return;
		}

	case FAN_mpack:
		{
		static const uint16_t maxsid = dist2sid(15.0f);//TODO param for 15.0f
		chaq::get_packing_values(distmx, L, maxsid, false, true, values);
		return;
		}

	case FAN_turnd:
		{
		static const uint16_t w = 5; //TODO param for w=5
		static const uint16_t undef_value =  1540; // measured median
		chaq::get_turnd_values(distmx, L, undef_value, values);
		return;
		}

	default:
		Die("chaq::fast_get_values(fan=%d=%s)", fan, FAN2str(fan));
		}
	}

void chaq::fast_get_codeseq(
	const flat_params &params,
	const flat_chain_t *chain,
	const sid_t *distmx,
	const chaq_vecs2 *cv,
	FAN fan,
	uint alpha_size,
	p_uint8_t codeseq,
	uint8_t *scratch_buffer,
	uint scratch_buffer_bytes)
	{
	const uint L = chain->get_length();
	switch (fan)
		{
	case FAN_aa:
		{
		switch (alpha_size)
			{
		case 3:		get_aa3_codeseq(chain->m_aa->m_data, L, codeseq);	return;
		case 4:		get_aa4_codeseq(chain->m_aa->m_data, L, codeseq);	return;
		case 20:	get_aa20_codeseq(chain, codeseq);					return;
			}
		Die("chaq::fast_get_codeseq(aa, alpha_size=%u)", alpha_size);	return;
		}

	case FAN_sec:
		{
		if (alpha_size == 32)
			memcpy(codeseq, cv->sec32_codeseq, L);
		else if (alpha_size == 4)
			sec32_codeseq_to_sec4(cv->sec32_codeseq, L, codeseq);
		else
			Die("chaq::fast_get_codeseq(sec) alpha_size=%u", alpha_size);
		return;
		}

	case FAN_nensec:
		{
		asserta(alpha_size == 32);
		uint8_t undef_code = get_undef_code(fan, alpha_size);
		for (uint pos = 0; pos < L; ++pos)
			{
			uint16_t xen = cv->nens[pos];
			asserta(xen < L || xen == UINT16_MAX);
			uint8_t code = (xen == UINT16_MAX) ?
				undef_code : cv->sec32_codeseq[xen];
			asserta(undef_code < alpha_size);
			codeseq[pos] = code;
			}
		return;
		}

	case FAN_rensec:
		{
		asserta(alpha_size == 32);
		uint8_t undef_code = get_undef_code(fan, alpha_size);
		for (uint pos = 0; pos < L; ++pos)
			{
			uint16_t xen = cv->rens[pos];
			asserta(xen < L || xen == UINT16_MAX);
			uint8_t code = (xen == UINT16_MAX) ?
				undef_code : cv->sec32_codeseq[xen];
			asserta(undef_code < alpha_size);
			codeseq[pos] = code;
			}
		return;
		}

	case FAN_pensec:
		{
		asserta(alpha_size == 32);
		uint8_t undef_code = get_undef_code(fan, alpha_size);
		for (uint pos = 0; pos < L; ++pos)
			{
			uint16_t xen = cv->pens[pos];
			asserta(xen < L || xen == UINT16_MAX);
			uint8_t code = (xen == UINT16_MAX) ?
				undef_code : cv->sec32_codeseq[xen];
			asserta(undef_code < alpha_size);
			codeseq[pos] = code;
			}
		return;
		}

	case FAN_mensec:
		{
		asserta(alpha_size == 32);
		uint8_t undef_code = get_undef_code(fan, alpha_size);
		for (uint pos = 0; pos < L; ++pos)
			{
			uint16_t xen = cv->mens[pos];
			asserta(xen < L || xen == UINT16_MAX);
			uint8_t code = (xen == UINT16_MAX) ?
				undef_code : cv->sec32_codeseq[xen];
			asserta(undef_code < alpha_size);
			codeseq[pos] = code;
			}
		return;
		}

	case FAN_pm:
		asserta(alpha_size == 2);
		chaq::get_pm_codeseq(cv->pensids, cv->mensids, L, codeseq);
		return;

	///////////////////////////////////////////////////////////////
	// Quantized features
	///////////////////////////////////////////////////////////////
	case FAN_nendist:
	case FAN_rendist:
	case FAN_pendist:
	case FAN_mendist:
	case FAN_fendist:
	case FAN_pmdd:
	case FAN_pmdiff:
	case FAN_pack:
	case FAN_ppack:
	case FAN_mpack:
	case FAN_angle:
	case FAN_turnd:
		{
		cp_uint16_t thresholds = get_thresholds(params, fan, alpha_size);
		const uint16_t undef_value = get_undef_value(params, fan, alpha_size);
		asserta(scratch_buffer_bytes >= sizeof(uint16_t)*L);
		uint16_t *values = (uint16_t *) scratch_buffer;
		chaq::fast_get_values(params, distmx, cv, chain, fan, alpha_size, values);
		for (uint i = 0; i < L; ++i)
			{
			uint16_t value = values[i];
			if (value == UINT16_MAX)
				value = undef_value;
			uint8_t code = get_bin(value, alpha_size, thresholds);
			assert(code < alpha_size);
			codeseq[i] = code;
			}
		return;
		}
	///////////////////////////////////////////////////////////////

	case FAN_kappa:
		{
		asserta(alpha_size == 32);//@@TODO
		slow_get_codeseq(params, chain, FAN_kappa, 32, codeseq);//@@TODO
		return;
		}

	default:
		Die("chaq::fast_get_codeseq(fan=%d %s)", fan, FAN2str(fan));
		}
	}

#if 0
void cmd_test_chaq_fast()
	{
	asserta(optset_alphadir);
	asserta(optset_feature);
	asserta(optset_alpha_size);
	asserta(optset_varstr);
	vector<string> param_names;
	vector<float> param_values;

	FAN fan = str2FAN(opt(feature));
	uint alpha_size = opt(alpha_size);

	parse_varstr(opt(varstr), param_names, param_values);
	const uint32_t M = flat_params::m_distmx_bandwidth;

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	const string &alphadir = opt(alphadir);
	flat_params params;
	params.init_from_alphadir(alphadir, alpha_names);
	Paralign::set_final_nu();

	vector<flat_chain_t *> chains;
	read_flat_chains(g_Arg1, chains);
	const uint nchain = uint(chains.size());
	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		ProgressStep(chainidx, nchain, "Working");
		const flat_chain_t *chain = chains[chainidx];
		const uint L = chain->get_length();

		uint8_t *slow_codeseq = myalloc(uint8_t, L);
		uint8_t *fast_codeseq = myalloc(uint8_t, L);

		size_t mem_bytes1 = chaq::get_fill_chaq_vecs_bytes_per_pos();
		scratch_mem mem1(mem_bytes1*L);
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain, distmx);

		chaq_vecs2 cv;
		chaq::fill_chaq_vecs2(distmx, L, cv);

		//size_t scratch_bytes2 = chaq::get_fast_get_codeseq_scratch_bytes_per_pos();
		uint scratch_buffer_bytes = 2*L;
		uint8_t *scratch_buffer = myalloc(uint8_t, scratch_buffer_bytes);

		chaq::slow_get_codeseq(
			params, chain, fan, alpha_size, slow_codeseq);
		chaq::fast_get_codeseq(
			params, chain, distmx, &cv, fan, alpha_size, fast_codeseq,
			scratch_buffer, scratch_buffer_bytes);

		for (uint pos = 0; pos < L; ++pos)
			{
			uint8_t slow_code = slow_codeseq[pos];
			asserta(slow_code < alpha_size);

			uint8_t fast_code = fast_codeseq[pos];
			asserta(fast_code < alpha_size);

			asserta(fast_code == slow_code);
			}

		myfree(slow_codeseq);
		myfree(fast_codeseq);
		}
	}
#endif
