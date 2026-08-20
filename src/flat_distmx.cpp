#include "myutils.h"
#include "flat_chain.h"
#include "flat_distmx.h"
#include "getticks.h"

void test_indexing(uint32_t M)
	{
	const int L = 250;
	const uint K = 250*M;
	vector<bool> touched_plus(K);
	for (int i = 0; i < L; ++i)
		{
		for (int j = i+1; j < L; ++j)
			{
			if (abs(i-j) > int(M))
				continue;
			uint32_t k = banded_ij_to_k(i, j);
			asserta(k < K);

			asserta(!touched_plus[k]);
			touched_plus[k] = true;

			uint32_t i2, j2;
			banded_k_to_ij(k, i2, j2);
			if (i2 != i)
				Die("i2=%u i=%u j=%u", i2, i, j);
			if (j2 != j)
				Die("j2=%u i=%u j=%u", j2, i, j);
			}
		}

	vector<bool> touched_minus(K);
	for (int i = 0; i < L; ++i)
		{
		for (int j = i+1; j < L; ++j)
			{
			if (abs(i-j) > int(M))
				continue;
			uint32_t k = banded_ij_to_k(i, j);
			asserta(k < K);

			asserta(!touched_minus[k]);
			touched_minus[k] = true;
			uint32_t i2, j2;
			banded_k_to_ij(k, i2, j2);
			if (i2 != i)
				Die("i2=%u i=%u j=%u", i2, i, j);
			if (j2 != j)
				Die("j2=%u i=%u j=%u", j2, i, j);
			}
		}
	ProgressLog("test_indexing OK\n");
	}

#if TRACE
static const uint trace_i = 0;
static const uint trace_j = 1;
#endif

static inline void fill_pen(sid_t *__restrict sdmx,
	uint32_t L, uint32_t M, uint32_t m, uint16_t *pen)
	{
	for (uint32_t i = 0; i < L; ++i)
		{
		const uint32_t jend = min(i+M, L-1);
		uint16_t pen_i = UINT16_MAX;
		uint32_t min_sd = UINT32_MAX;
		uint32_t k = i*M + m - 1;
		for (uint32_t j = i + m; j <= jend; ++j)
			{
			assert(k == banded_ij_to_k(i, j));
			uint32_t sd = sdmx[k++];
			if (sd < min_sd)
				{
				min_sd = sd;
				pen_i = j;
				}
			}
		pen[i] = pen_i;
		}
	}

static inline void fill_men(sid_t *__restrict sdmx,
	uint32_t L, uint16_t *men)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;
	for (int i = 0; i < int(L); ++i)
		{
		const uint32_t jstart = min(i+M, L-1);
		uint16_t men_i = UINT16_MAX;
		uint32_t min_sd = UINT32_MAX;
		for (int j = max(0,i-int(M)); j <= i-int(m); ++j)
			{
			uint32_t k = banded_ij_to_k(i, j);
			uint32_t sd = sdmx[k];
			if (sd < min_sd)
				{
				min_sd = sd;
				men_i = j;
				}
			}
		men[i] = men_i;
		}
	}

static uint compare_fill(const flat_chain_t &chain, const sid_t *sdmx)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;
	const uint L = chain.get_length();
	uint diffs = 0;
	for (int i = 0; i < int(L); ++i)
		{
		ic_t icx_i, icy_i, icz_i;
		chain.get_ic_xyz(i, icx_i, icy_i, icz_i);
		for (int j = i+1; j < int(L); ++j)
			{
			if (i==j || abs(i-j) > int(M))
				continue;

			ic_t icx_j, icy_j, icz_j;
			chain.get_ic_xyz(j, icx_j, icy_j, icz_j);

			uint32_t sid = icxyzpair2sid(
				icx_i, icy_i, icz_i,
				icx_j, icy_j, icz_j);

			uint k = banded_ij_to_k(i, j);
			uint32_t sid2 = sdmx[k];
			if (sid2 != sid)
				++diffs;

			if (opt_verbose)
				{
				Log("\n");
				Log("i=%u j=%u\n", i, j);
				Log(" xyz(%u) = %u,%u,%u", i, icx_i, icy_i, icz_i);
				Log(" = %.1f, %.1f, %.1f\n", ic2coord(icx_i), ic2coord(icy_i), ic2coord(icz_i));
				Log(" xyz(%u) = %u,%u,%u", j, icx_j, icy_j, icz_j);
				Log(" = %.1f, %.1f, %.1f\n", ic2coord(icx_j), ic2coord(icy_j), ic2coord(icz_j));
				Log(" sid = %u = %.1f A", sid, sid2dist(sid));
				Log("\n");
				}
			}
		}
	return diffs;
	}

static uint compare_pen(const flat_chain_t &chain, uint M, uint m,
	const uint16_t *pen)
	{
	const uint L = chain.get_length();
	uint diffs = 0;
	for (int i = 0; i < int(L); ++i)
		{
		const int jend = min(i+int(M), int(L)-1);
		float MinDist = FLT_MAX;
		uint16_t pen_i = UINT16_MAX;
		for (int j = i + int(m); j <= jend; ++j)
			{
			if (i==j || abs(i-j) < int(m) || abs(i-j) > int(M))
				continue;
			float d = chain.slow_float_dist(i, j);
			if (d < MinDist)
				{
				MinDist = d;
				pen_i = j;
				}
			}
		if (pen_i != pen[i])
			++diffs;
		}
	return diffs;
	}

static uint compare_men(const flat_chain_t &chain, const uint16_t *men)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;
	const uint L = chain.get_length();
	uint diffs = 0;
	for (int i = 0; i < int(L); ++i)
		{
		const int jend = min(i+int(M), int(L)-1);
		float MinDist = FLT_MAX;
		uint16_t men_i = UINT16_MAX;
		for (int j = 0; j < i; ++j)
			{
			if (i==j || abs(i-j) < int(m) || abs(i-j) > int(M))
				continue;
			float d = chain.slow_float_dist(i, j);
			if (d < MinDist)
				{
				MinDist = d;
				men_i = j;
				}
			}
		if (men_i != men[i])
			++diffs;
		}
	return diffs;
	}

static void test_distmx(const vector<flat_chain_t *> &chains, uint M)
	{
	const uint ChainCount = SIZE(chains);
	vector<uint16_t> ICs;
	uint total_diffs = 0;
	TICKS total_ticks = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "working diffs %u", total_diffs);
		const flat_chain_t &chain = *chains[ChainIdx];
		const uint L = chain.get_length();
		const uint K = L*M;
		sid_t *sdmx = myalloc(sid_t, K);
		TICKS t1 = GetClockTicks();
		fill_flat_distmx(chain.m_xyz->m_data, L, sdmx);
		TICKS t2 = GetClockTicks();
		total_ticks += t2 - t1;
		uint diffs = compare_fill(chain, sdmx);
		total_diffs += diffs;
		}
	ProgressLog("%.3g ticks, %u diffs sd\n", double(total_ticks), total_diffs);
	}

static void test_pen(const vector<flat_chain_t *> &chains, uint M, uint m)
	{
	const uint ChainCount = SIZE(chains);
	vector<uint16_t> ICs;
	uint total_diffs = 0;
	TICKS total_ticks = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "working diffs %u", total_diffs);
		const flat_chain_t &chain = *chains[ChainIdx];
		const uint L = chain.get_length();
		const uint K = L*M;
		const uint16_t *xyz = chain.m_xyz->m_data;
		sid_t *sdmx = myalloc(sid_t, K);
		fill_flat_distmx(xyz, L, sdmx);

		uint16_t *pen = myalloc(uint16_t, L);

		TICKS t1 = GetClockTicks();
		fill_pen(sdmx, L, m, M, pen);
		TICKS t2 = GetClockTicks();

		total_ticks += t2 - t1;
		uint diffs = compare_pen(chain, m, M, pen);
		total_diffs += diffs;
		}
	ProgressLog("%.3g ticks, %u diffs pen\n", double(total_ticks), total_diffs);
	}

static void test_men(const vector<flat_chain_t *> &chains, uint M, uint m)
	{
	const uint ChainCount = SIZE(chains);
	vector<uint16_t> ICs;
	uint total_diffs = 0;
	TICKS total_ticks = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "working diffs %u", total_diffs);
		const flat_chain_t &chain = *chains[ChainIdx];
		const uint L = chain.get_length();
		const uint K = L*M;
		const uint16_t *xyz = chain.m_xyz->m_data;
		sid_t *sdmx = myalloc(sid_t, K);
		fill_flat_distmx(xyz, L, sdmx);

		uint16_t *men = myalloc(uint16_t, L);

		TICKS t1 = GetClockTicks();
		fill_men(sdmx, L, men);
		TICKS t2 = GetClockTicks();

		total_ticks += t2 - t1;
		uint diffs = compare_men(chain, men);
		total_diffs += diffs;
		}
	ProgressLog("%.3g ticks, %u diffs men\n", double(total_ticks), total_diffs);
	}

#if 0
void cmd_test_flat_distmx()
	{
	const uint M = flat_params::m_distmx_bandwidth;
	test_indexing(M);
	vector<flat_chain_t *> chains;
	read_flat_chains(g_Arg1, chains);
	test_distmx(chains, M);
	test_pen(chains, M, 16);
	test_men(chains, M, 16);
	}
#endif