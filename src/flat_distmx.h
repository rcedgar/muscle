#pragma once

#include "flat_dist_types.h"
#include "flat_params.h"

static inline uint32_t banded_i_lt_j_to_k(uint32_t i, uint32_t j)
    {
	const uint M = flat_params::m_distmx_bandwidth;
	assert(i < j);
	assert(abs(int(i)-int(j)) <= int(M));
    uint32_t offset = j - i;
    return M*i + offset - 1;
    }

static inline uint32_t banded_ij_to_k(uint32_t i, uint32_t j)
    {
	const uint M = flat_params::m_distmx_bandwidth;
	assert(i != j);
	assert(abs(int(i)-int(j)) <= int(M));
    if (j < i) std::swap(i, j);
    uint32_t offset = j - i;
    return M*i + offset - 1;
    }

static inline void banded_k_to_ij(uint32_t k, uint32_t& i, uint32_t& j) 
    {
	const uint M = flat_params::m_distmx_bandwidth;
    i = k / M;
    uint32_t offset = (k % M) + 1;
    j = i + offset;
    }

static inline void fill_flat_distmx(
	cp_ic_t xyz,
	uint32_t L,
	p_sid_t sdmx)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	uint i3 = 0;
	for (uint32_t i = 0; i < L; ++i)
		{
		int32_t icx_i = xyz[i3++];
		int32_t icy_i = xyz[i3++];
		int32_t icz_i = xyz[i3++];
		uint32_t k = i*M;
		const uint32_t jend = min(i+M, L-1);
		for (uint32_t j = i + 1; j <= jend; ++j)
			{
			int32_t icx_j = xyz[3*j];
			int32_t icy_j = xyz[3*j+1];
			int32_t icz_j = xyz[3*j+2];
			sid_t sd = icxyzpair2sid(
				icx_i, icy_i, icz_i,
				icx_j, icy_j, icz_j);
			assert(k == banded_ij_to_k(i, j));
			sdmx[k++] = sd;
			}
		}
	}
