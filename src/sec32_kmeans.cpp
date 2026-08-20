#include "myutils.h"
#include "flat_chain.h"
#include "flat_params.h"
#include "sec_kmeans.h"
#include "chaq.h"
#include "alpha.h"

struct sec32_best_t
	{
	uint32_t dist = UINT32_MAX;
	const flat_chain_t *chain = 0;
	uint pos_0based = UINT_MAX;
	};

static string match_window_aa(const flat_chain_t *chain, uint pos, int w)
	{
	asserta(chain != 0);
	asserta(chain->m_aa != 0);
	const uint lo = pos - uint(w);
	const uint hi = pos + uint(w);
	asserta(hi < chain->get_length());
	return string(chain->m_aa->m_data + lo, size_t(2*w + 1));
	}

void cmd_match_sec32()
	{
	asserta(optset_output);
	vector<flat_chain_t *> chains;
	read_flat_chains(g_Arg1, chains);

	sec_kmeans SK;
	SK.from_sec_n(32);
	asserta(SK.m_K == 32);

	const uint M = flat_params::m_distmx_bandwidth;
	const int w = SK.m_w;
	sec32_best_t best[32];
	sid_t tmpv[32];
	asserta(SK.m_D <= 32);

	for (uint chain_idx = 0; chain_idx < SIZE(chains); ++chain_idx)
		{
		const flat_chain_t *chain = chains[chain_idx];
		const uint L = chain->get_length();
		if (int(L) < 2*w + 1)
			continue;
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
		for (int pos = w; pos <= int(L) - w - 1; ++pos)
			{
			SK.get_v(distmx, pos, L, tmpv);
			for (uint k = 0; k < SK.m_K; ++k)
				{
				uint32_t d = SK.get_dist_early_quit(
					tmpv, SK.m_means + k*SK.m_D, best[k].dist);
				if (d < best[k].dist)
					{
					best[k].dist = d;
					best[k].chain = chain;
					best[k].pos_0based = uint(pos);
					}
				}
			}
		myfree(distmx);
		}

	FILE *fOut = CreateStdioFile(opt(output));
	fprintf(fOut, "sec32\tchain\tpos\tseq\n");
	for (uint k = 0; k < SK.m_K; ++k)
		{
		if (best[k].dist == UINT32_MAX)
			Die("cmd_match_sec32: no valid position for sec32 letter %c",
				char(g_LetterToCharMu[k]));
		asserta(best[k].chain != 0);
		const string seq = match_window_aa(best[k].chain, best[k].pos_0based, w);
		asserta(seq.size() == size_t(2*w + 1));
		fprintf(fOut, "%c\t%s\t%u\t%s\n",
			char(g_LetterToCharMu[k]),
			best[k].chain->m_label.c_str(),
			best[k].pos_0based + 1,
			seq.c_str());
		}
	CloseStdioFile(fOut);
	}
