#include <vector>
#include <cstdint>
#include <cassert>
#include <limits>
#include <algorithm>

using std::vector;
using std::uint16_t;
using std::uint32_t;
using std::uint64_t;

struct QuantizeResult
	{
	// thresholds[i] is the largest value in bin i, for i=0..K-2.
	// Bin 0: x <= thresholds[0]
	// Bin 1: thresholds[0] < x <= thresholds[1]
	// ...
	// Bin K-1: x > thresholds[K-2]
	vector<uint16_t> thresholds;

	// Actual population in each bin, size K.
	vector<uint64_t> bin_counts;

	// Sum over bins of (bin_count - N/K)^2, multiplied by K^2.
	// This scaling avoids floating point in the DP.
	uint64_t scaled_sse = 0;
	uint64_t sum_count = 0;
	uint16_t median_value = UINT16_MAX;
	};

static inline QuantizeResult quantize_histogram_equal_mass_dp(
	const vector<uint16_t>& counts,
	unsigned K)
	{
	assert(counts.size() == 65536);
	assert(K >= 1);

	QuantizeResult qr;

	// Compress histogram to nonzero masses.
	vector<uint16_t> vals;
	vector<uint32_t> w;
	vals.reserve(65536);
	w.reserve(65536);

	uint64_t N = 0;
	for (uint32_t v = 0; v < 65536; ++v)
		{
		uint32_t c = counts[v];
		if (c)
			{
			vals.push_back((uint16_t)v);
			w.push_back(c);
			N += c;
			}
		}
	qr.sum_count = N;

	uint32_t sum_count = 0;
	for (uint32_t v = 0; v < 65536; ++v)
		{
		sum_count += counts[v];
		if (sum_count >= N/2)
			{
			qr.median_value = v;
			break;
			}
		}

	const unsigned M = (unsigned)w.size();

	qr.thresholds.resize(K > 0 ? K - 1 : 0, 0);
	qr.bin_counts.resize(K, 0);

	if (K == 1)
		{
		qr.bin_counts[0] = N;
		qr.scaled_sse = 0; // single bin, trivial
		return qr;
		}

	if (N == 0)
		{
		// Degenerate: no data
		std::fill(qr.thresholds.begin(), qr.thresholds.end(), uint16_t(0));
		std::fill(qr.bin_counts.begin(), qr.bin_counts.end(), uint64_t(0));
		qr.scaled_sse = 0;
		return qr;
		}

	if (M == 0)
		{
		// Same as N==0, just for completeness
		std::fill(qr.thresholds.begin(), qr.thresholds.end(), uint16_t(0));
		std::fill(qr.bin_counts.begin(), qr.bin_counts.end(), uint64_t(0));
		qr.scaled_sse = 0;
		return qr;
		}

	// Prefix sums of masses: pref[i] = sum of first i masses, pref[0]=0.
	vector<uint64_t> pref(M + 1, 0);
	for (unsigned i = 0; i < M; ++i)
		pref[i + 1] = pref[i] + w[i];

	// Cost of grouping masses [a, b) into one bin:
	//   ( sum(w[a:b]) - N/K )^2
	// To avoid fractions, minimize:
	//   (K*sum - N)^2
	auto scaled_bin_cost = [&](unsigned a, unsigned b) -> uint64_t
		{
		uint64_t s = pref[b] - pref[a];
		uint64_t ks = uint64_t(K) * s;
		uint64_t d = (ks >= N ? ks - N : N - ks);
		return d * d;
		};

	const uint64_t INF = std::numeric_limits<uint64_t>::max() / 4;

	// dp[k][i] = minimum scaled SSE for partitioning first i masses into k bins.
	// i ranges 0..M, k ranges 0..K.
	vector<vector<uint64_t>> dp(K + 1, vector<uint64_t>(M + 1, INF));
	vector<vector<unsigned>> prev(K + 1, vector<unsigned>(M + 1, 0));

	dp[0][0] = 0;

	// Allow empty bins.
	// That is important when K > number of useful splits.
	// Recurrence:
	//   dp[k][i] = min over j in [0..i] of dp[k-1][j] + cost(j,i)
	// If j==i, this creates an empty final bin with cost (0 - N/K)^2.
	//
	// Note: this is O(K*M^2). For uint16_t histograms M <= 65536, but in real
	// usage M is often much smaller. For moderate K this is usually fine.
	for (unsigned k = 1; k <= K; ++k)
		{
		for (unsigned i = 0; i <= M; ++i)
			{
			uint64_t best = INF;
			unsigned best_j = 0;

			for (unsigned j = 0; j <= i; ++j)
				{
				if (dp[k - 1][j] == INF)
					continue;

				uint64_t c = scaled_bin_cost(j, i);
				uint64_t cand = dp[k - 1][j] + c;
				if (cand < best)
					{
					best = cand;
					best_j = j;
					}
				}

			dp[k][i] = best;
			prev[k][i] = best_j;
			}
		}

	qr.scaled_sse = dp[K][M];

	// Recover bin boundaries in compressed index space.
	// cuts[k] = start index of bin k, with cuts[0]=0 and cuts[K]=M.
	vector<unsigned> cuts(K + 1, 0);
	cuts[K] = M;
	{
	unsigned i = M;
	for (unsigned k = K; k >= 1; --k)
		{
		unsigned j = prev[k][i];
		cuts[k - 1] = j;
		i = j;
		if (k == 1)
			break;
		}
	}

	// Bin counts
	for (unsigned k = 0; k < K; ++k)
		qr.bin_counts[k] = pref[cuts[k + 1]] - pref[cuts[k]];

	// Convert cuts to uint16 thresholds.
	//
	// Bin k contains compressed indices [cuts[k], cuts[k+1]).
	// Threshold between bin k and k+1 is the largest original value in bin k.
	//
	// If bin k is empty, use the previous threshold if possible; otherwise 0.
	uint16_t last_threshold = 0;
	for (unsigned k = 0; k + 1 < K; ++k)
		{
		unsigned a = cuts[k];
		unsigned b = cuts[k + 1];

		uint16_t t;
		if (a < b)
			{
			t = vals[b - 1];
			}
		else
			{
			// Empty bin: repeated threshold
			t = last_threshold;
			}

		qr.thresholds[k] = t;
		last_threshold = t;
		}

	return qr;
	}

static inline uint8_t get_bin(uint16_t x, uint8_t alpha_size,
	const uint16_t *thresholds)
	{
	uint8_t lo = 0, hi = alpha_size - 1;
	while (lo < hi)
		{
		uint8_t mid = (lo + hi) >> 1;
		if (x <= thresholds[mid])
			hi = mid;
		else
			lo = mid + 1;
		}
	assert(lo < alpha_size);
	return lo;
	}