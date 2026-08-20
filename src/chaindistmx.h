#pragma once

#include "flat_base.h"

class chaindistmx_t : public flat_mx<sid_t, FE_chaindistmx>
	{
protected:
	chaindistmx_t(uint32_t n)
        : flat_mx<sid_t, FE_chaindistmx>(n, 0)
		{
		}

	chaindistmx_t(uint32_t rows, uint32_t cols)
        : flat_mx<sid_t, FE_chaindistmx>(rows, cols)
		{
		}

public:
#if TRACK_SRC
	static chaindistmx_t *newflat_src(uint32_t n,
		const char *srcfile, int srcline)
		{
		chaindistmx_t *p = new chaindistmx_t(n, 0);
		p->m_srcfile = srcfile;
		p->m_srcline = srcline;
		return p;
		}

	static chaindistmx_t *newflat_src(uint32_t L, uint32_t M,
		const char *srcfile, int srcline)
		{
		chaindistmx_t *p = new chaindistmx_t(L, M);
		p->m_srcfile = srcfile;
		p->m_srcline = srcline;
		return p;
		}
#else
	static chaindistmx_t *newflat(uint32_t n)
		{
		chaindistmx_t *p = new chaindistmx_t(n, 0);
		return p;
		}

	static chaindistmx_t *newflat(uint32_t L, uint32_t M)
		{
		chaindistmx_t *p = new chaindistmx_t(L, M);
		return p;
		}
#endif
	};
