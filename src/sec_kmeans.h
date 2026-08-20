#pragma once

#include "flat_distmx.h"
#include "chaq.h"
#include "sort.h"
#include "tabbedlines.h"
#include "alpha.h"

// Cluster subset of local distance
//	matrix by k-meeans clustering
class sec_kmeans
	{
public:
	// Parameters
	/////////////
	uint m_K = 0;				// number of clusters for K-means
	uint m_D = 0;				// dimension of feature vector, length of m_i/jvalues
	int m_w = 0;				// band width for sec (e.g. 3), max index in m_i/jvalues
	int* m_off1s = 0;			// +/- offsets from position
	int* m_off2s = 0;			// +/- offsets from position
	sid_t *m_means = 0;			// flat matrix of current means size m_K x m_D

	// Training data
	////////////////
	const vector<flat_chain_t *> *m_chains = 0;
	uint m_N = 0;				// number of residues, size of m_vs
	sid_t *m_vs = 0;			// flat matrix of feature vectors size m_N x m_D
	uint* m_cluster_idxs = 0;	// current cluster assignments
	uint* m_cluster_sizes = 0;	// cluster sizes
	uint* m_size_order = 0;
	uint m_zero_count = 0;
	uint m_nrchanges = 0;

public:
	void clear_params()
		{
		m_K = 0;
		m_D = 0;
		m_w = 0;
		myfree(m_off1s);
		myfree(m_off2s);
		m_off1s = 0;
		m_off2s = 0;		// +/- offsets from position
		}

	void clear_training_data()
		{
		myfree(m_vs);
		myfree(m_cluster_idxs);
		myfree(m_cluster_sizes);
		myfree(m_size_order);

		m_N = 0;
		m_vs = 0;
		m_cluster_idxs = 0;
		m_cluster_sizes = 0;
		m_size_order = 0;
		m_zero_count = 0;
		m_nrchanges = 0;
		}

	void clear()
		{
		clear_params();
		clear_training_data();
		}

	void to_tsv(const string &fn) const
		{
		vector<string> lines;
		to_lines(lines);
		tabbedlines::to_tsv(fn, lines);
		}

	void from_tsv(const string &fn)
		{
		vector<string> lines;
		ReadLinesFromFile(fn, lines);
		from_lines(lines);
		}

	void to_lines(vector<string> &lines) const
		{
		tabbedlines tl;

		tl.put_int("sec", m_K);
		tl.put_int("dim", m_D);
		tl.put_signed_int_flat_vec("offs1", m_off1s, m_D);
		tl.put_signed_int_flat_vec("offs2", m_off2s, m_D);
		sid_t *sorted_means = get_sorted_means();
		tl.put_int16_flat_vec("mean", sorted_means, m_K*m_D);
		myfree(sorted_means);

		lines = tl.m_lines;
		}

	uint calc_w() const
		{
		uint max_off = 0;
		for (uint i = 0; i < m_D; ++i)
			{
			max_off = max(max_off, uint(abs(m_off1s[i])));
			max_off = max(max_off, uint(abs(m_off2s[i])));
			}
		return max_off;
		}

	void alloc_DK()
		{
		asserta(m_D > 0);
		asserta(m_K > 0);
		}

	void from_sec_n(uint alpha_size);

	void from_lines(const vector<string> &lines)
		{
		clear();
		tabbedlines tl(lines);

		m_K = tl.get_int("sec");
		m_D = tl.get_int("dim");
		m_off1s = tl.get_signed_int_flat_vec("offs1", m_D);
		m_off2s = tl.get_signed_int_flat_vec("offs2", m_D);
		m_means = tl.get_int16_flat_vec("mean", m_K*m_D);
		tl.get_eof();

		m_w = calc_w();

		alloc_DK();
		}

	void log_params() const
		{
		Log("K %u, N %u, D %u, w %u\n", m_K, m_N, m_D, m_w);
		Log("off1s[%u] =", m_D);
		for (uint i = 0; i < m_D; ++i)
			{
			if (i > 0)
				Log(",");
			Log(" %2d", m_off1s[i]);
			}
		Log("\n");

		Log("off2s[%u] =", m_D);
		for (uint i = 0; i < m_D; ++i)
			{
			if (i > 0)
				Log(",");
			Log(" %2d", m_off2s[i]);
			}
		Log("\n");
		}

	void log_v(const sid_t *v) const
		{
		for (uint i = 0; i < m_D; ++i)
			Log(" %5u(%5.2f)", v[i], sid2dist(v[i]));
		Log("\n");
		}

	void log_means() const
		{
		Log("\nmeans:\n");
		for (uint i = 0; i < m_K; ++i)
			{
			uint cluster_idx = (m_size_order == 0 ? i : m_size_order[i]);
			if (m_cluster_sizes == 0)
				Log("%3u [-] ", cluster_idx);
			else
				{
				double pct = GetPct(m_cluster_sizes[cluster_idx], m_N);
				Log("%3u [%6.1f%%] ", cluster_idx, pct);
				}
			log_v(m_means + cluster_idx*m_D);
			}
		}

	void log_head_vs(uint n=10) const
		{
		if (m_vs == 0)
			{
			Log("m_vs=nullptr\n");
			return;
			}
		Log("\nhead_vs(%u):\n", n);
		for (uint residue_idx = 0; residue_idx < min(n, m_N); ++residue_idx)
			{
			Log("[%3u] ", residue_idx);
			log_v(m_vs + residue_idx*m_D);
			}
		}

	void log_random_vs(uint n=100) const
		{
		Log("\nrandom_vs(%u):\n", n);
		for (uint i = 0; i < n; ++i)
			{
			uint residue_idx = randu32()%m_N;
			Log("[%7u] ", residue_idx);
			log_v(m_vs + residue_idx*m_D);
			}
		}

	void logme() const
		{
		log_params();
		log_head_vs();
		log_means();
		}

	// Adjacent residues in the backbone should have distance ~3.81 A
	// 5162 / 184690 bad backbones (3%)
	bool check_backbone(uint chain_idx, const sid_t* distmx, int pos, int L)
		{
		const sid_t backbone_sid = dist2sid(3.81f);
		const sid_t min_backbone_sid = backbone_sid - 10;
		const sid_t max_backbone_sid = backbone_sid + 10;
		assert(pos >= m_w && pos + m_w < L);
		for (int i = -int(m_w); i < int(m_w); ++i)
			{
			int ifirst_pos = int(pos)+i;
			assert(ifirst_pos >= 0 && ifirst_pos + 1 < L);
			uint32_t first_pos = uint32_t(ifirst_pos);
			uint k = banded_ij_to_k(first_pos, first_pos+1);
			sid_t sid = distmx[k];
			if (sid < min_backbone_sid || sid > max_backbone_sid)
				return false;
			}
		return true;
		}

	void get_v(const sid_t* distmx, int pos, int L, sid_t* v) const
		{
		assert(pos >= m_w && pos + m_w < L);
		for (uint m = 0; m < m_D; ++m)
			{
			int off1 = m_off1s[m];
			int off2 = m_off2s[m];
			uint k = banded_ij_to_k(pos+off1, pos+off2);
			sid_t sid = distmx[k];
			v[m] = sid;
			}
		}

	// Euclidean squared distance (no need to sqrt)
	uint32_t get_dist(const sid_t* v1, const sid_t* v2) const
		{
		uint32_t sum2 = 0;
		for (uint m = 0; m < m_D; ++m)
			{
			int32_t diff = int32_t(v1[m]) - int32_t(v2[m]);
			sum2 += uint32_t(diff*diff);
			}
		return sum2;
		}

	// When looking for best match can give up early
	uint32_t get_dist_early_quit(const sid_t* v1, const sid_t* v2,
		uint32_t smallest_so_far) const
		{
		uint32_t sum2 = 0;
		for (uint m = 0; m < m_D; ++m)
			{
			int32_t diff = int32_t(v1[m]) - int32_t(v2[m]);
			sum2 += diff*diff;
			if (sum2 >= smallest_so_far)
				return sum2;
			}
		return sum2;
		}

	uint8_t assign_cluster(const sid_t* v) const
		{
		uint best_cluster = UINT32_MAX;
		uint32_t min_dist = UINT32_MAX;
		for (uint cluster_idx = 0; cluster_idx < m_K; ++cluster_idx)
			{
			uint32_t d = get_dist_early_quit(v, m_means + cluster_idx*m_D, min_dist);
			if (d < min_dist)
				{
				min_dist = d;
				best_cluster = cluster_idx;
				}
			}
		assert(best_cluster != UINT32_MAX);
		assert(best_cluster < UINT8_MAX);
		return uint8_t(best_cluster);
		}

	void assign_random_means()
		{
		for (uint cluster_idx = 0; cluster_idx < m_K; ++cluster_idx)
			{
			uint residue_idx = randu32()%m_N;
			memcpy(m_means + cluster_idx*m_D, m_vs + residue_idx*m_D, m_D*sizeof(sid_t));
			}
		}

	sid_t *get_sorted_means() const
		{
		sid_t *sorted_means = myalloc(sid_t, m_D*m_K);
		for (uint i = 0; i < m_K; ++i)
			{
			uint j = m_size_order[i];
			memcpy(sorted_means + i*m_D, m_means + j*m_D, m_D*sizeof(sid_t));
			}
		return sorted_means;
		}

	uint assign_clusters()
		{
		myfree(m_cluster_sizes);
		myfree(m_size_order);
		m_cluster_sizes = myalloc(uint, m_K);
		m_size_order = myalloc(uint, m_K);
		zero_array(m_cluster_sizes, m_K);
		uint nrchanges = 0;
		for (uint residue_idx = 0; residue_idx < m_N; ++residue_idx)
			{
			uint old_cluster_idx = m_cluster_idxs[residue_idx];
			uint new_cluster_idx = assign_cluster(m_vs + residue_idx*m_D);
			if (new_cluster_idx != old_cluster_idx)
				{
				++nrchanges;
				m_cluster_idxs[residue_idx] = new_cluster_idx;
				}
			++m_cluster_sizes[new_cluster_idx];
			}
		QuickSortOrderDesc(m_cluster_sizes, m_K, m_size_order);
		return nrchanges;
		}

	sid_t get_random_value(uint d) const
		{
		assert(d < m_D);
		uint residue_idx = randu32()%m_N;
		return m_vs[residue_idx*m_D + d];
		}

	uint calc_means()
		{
		uint zero_count = 0;
		uint n = m_K*m_D;
		uint64_t *sums = myalloc(uint64_t, n);
		uint32_t *residue_counts = myalloc(uint32_t, m_K);
		zero_array(sums, n);
		zero_array(residue_counts, m_K);
#if DEBUG
		uint32_t *check_counts = myalloc(uint32_t, n);
		zero_array(check_counts, n);
#endif

		for (uint residue_idx = 0; residue_idx < m_N; ++residue_idx)
			{
			uint cluster_idx = m_cluster_idxs[residue_idx];
			++residue_counts[cluster_idx];
			for (uint d = 0; d < m_D; ++d)
				{
				sums[cluster_idx*m_D + d] += m_vs[residue_idx*m_D + d];
#if DEBUG
				check_counts[cluster_idx*m_D + d] += 1;
#endif
				}
			}

		uint sum_residue_count = 0;
		for (uint cluster_idx = 0; cluster_idx < m_K; ++cluster_idx)
			{
			uint residue_count = residue_counts[cluster_idx];
			sum_residue_count += residue_count;
			if (residue_count == 0)
				{
				++zero_count;
				uint random_residue_idx = randu32()%m_N;
				for (uint d = 0; d < m_D; ++d)
					m_means[cluster_idx*m_D + d] = m_vs[random_residue_idx*m_D + d];
				continue;
				}

			for (uint d = 0; d < m_D; ++d)
				{
#if DEBUG
				uint32_t check_count = check_counts[cluster_idx*m_D + d];
				assert(check_count == residue_count);
#endif
				uint64_t mean64 = sums[cluster_idx*m_D + d]/residue_count;
				sid_t mean = sid_t(mean64);
				asserta(uint64_t(mean) == mean64);
				m_means[cluster_idx*m_D + d] = mean;
				}
			}
		assert(sum_residue_count == m_N);
		myfree(sums);
		myfree(residue_counts);
#if DEBUG
		myfree(check_counts);
#endif
		return zero_count;
		}

	void init(uint K, uint M,
		const vector<int> &off1s,
		const vector<int> &off2s)
		{
		m_K = K;

		m_D = SIZE(off1s);
		asserta(SIZE(off2s) == m_D);

		m_off1s = myalloc(int, m_D);
		m_off2s = myalloc(int, m_D);

		memcpy(m_off1s, off1s.data(), m_D*sizeof(int));
		memcpy(m_off2s, off2s.data(), m_D*sizeof(int));

		m_w = 0;
		for (uint i = 0; i < m_D; ++i)
			{
			m_w = max(m_w, abs(off1s[i]));
			m_w = max(m_w, abs(off2s[i]));
			}
		myfree(m_means);
		m_means = myalloc(sid_t, m_K*m_D);
		}

	void set_vs(const vector<flat_chain_t *> &chains)
		{
		const uint M = flat_params::m_distmx_bandwidth;
		m_N = 0;
		m_chains = &chains;
		asserta(m_D > 0);
		const uint nrchains = SIZE(chains);

		uint total_length = 0;
		for (uint chainidx = 0; chainidx < nrchains; ++chainidx)
			{
			uint L = chains[chainidx]->get_length();
			total_length += L;
			}

		// Will skip some residues, total_length is > m_N
		myfree(m_vs);
		m_vs = myalloc(sid_t, m_D*total_length);
#if DEBUG
		memset(m_vs, 0xff, m_D*total_length*sizeof(sid_t));
#endif

		uint residue_idx = 0;
		uint bad_backbones = 0;
		for (uint chain_idx = 0; chain_idx < nrchains; ++chain_idx)
			{
			const flat_chain_t* chain = chains[chain_idx];
			const int L = (int) chain->get_length();

			sid_t *distmx = myalloc(sid_t, L*M);
			chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
			for (int pos = m_w; pos < L - m_w; ++pos)
				{
				bool ok = check_backbone(chain_idx, distmx, pos, L);
				if (!ok)
					{
					++bad_backbones;
					continue;
					}
				get_v(distmx, pos, L, m_vs + m_D*residue_idx++);
				}
			}
		m_N = residue_idx;
		myfree(m_cluster_idxs);
		m_cluster_idxs = myalloc(uint, m_N);
		ProgressLog("%u / %u bad backbones\n", bad_backbones, m_N);
		}

	void get_codeseq(
		const sid_t *distmx,
		uint L,
		uint8_t *codeseq) const;

	void run_iter()
		{
		m_zero_count = calc_means();
		m_nrchanges = assign_clusters();
		}

	void train(uint niter)
		{
		assign_random_means();
		assign_clusters();
		for (uint iter = 0; iter < niter; ++iter)
			{
			run_iter();
			ProgressLog("iter %u, changes %u",
				iter, m_nrchanges);
			if (m_zero_count > 0)
				ProgressLog(", zero %u", m_zero_count);
			ProgressLog("\n");
			if (m_nrchanges == 0)
				{
				ProgressLog("Converged\n");
				break;
				}
			}

		}

	void ss4stats()
		{
		const uint M = flat_params::m_distmx_bandwidth;
		assert(m_chains);
		vector<vector<uint> > countmx(m_K);
		for (uint i = 0; i < m_K; ++i)
			countmx[i].resize(4);

		uint nrchains = SIZE(*m_chains);
		for (uint chainidx = 0; chainidx < nrchains; ++chainidx)
			{
			const flat_chain_t* chain =(*m_chains)[chainidx];
			const uint L = chain->get_length();
			sid_t *distmx = myalloc(sid_t, L*M);
			chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);

			uint8_t *codeseq = myalloc(uint8_t, L);
			get_codeseq(distmx, L, codeseq);

			uint8_t *ss4codeseq = myalloc(uint8_t, L);
			chaq::get_ss4_codeseq(distmx, L, ss4codeseq);

			for (uint pos = 2; pos < L - 2; ++pos)
				{
				uint8_t letter = codeseq[pos];
				uint8_t ss4letter = ss4codeseq[pos];
				countmx[letter][ss4letter] += 1;
				}
			myfree(ss4codeseq);
			}

		ProgressLog("X    Helix   Strand     Turn     Loop\n");
		for (uint j = 0; j < m_K; ++j)
			{
			uint i = m_size_order[j];
			ProgressLog("%c", g_LetterToCharMu[j]);
			for (uint j = 0; j < 4; ++j)
				ProgressLog("  %7u", countmx[i][j]);
			ProgressLog("  [%2u]", j);
			ProgressLog("  %6.1f%%", GetPct(m_cluster_sizes[i], m_N));
			ProgressLog("\n");
			}
		}
public:
	static sec_kmeans *m_SK2;
	static sec_kmeans *m_SK3;
	static sec_kmeans *m_SK4;
	static sec_kmeans *m_SK8;
	static sec_kmeans *m_SK16;
	static sec_kmeans *m_SK32;
	static void get_sec2_lines(vector<string> &lines);
	static void get_sec3_lines(vector<string> &lines);
	static void get_sec4_lines(vector<string> &lines);
	static void get_sec8_lines(vector<string> &lines);
	static void get_sec16_lines(vector<string> &lines);
	static void get_sec32_lines(vector<string> &lines);
	static void get_sec_lines(uint alpha_size, vector<string> &lines);
	static void get_conf_lines(vector<string> &lines);
	static sec_kmeans *get_SK(uint alpha_size, uint M);
	};
