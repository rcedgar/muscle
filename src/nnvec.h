#pragma once

#include "flat_base.h"

class nnvec_t : public flat_vec<uint16_t, FE_nnvec>
	{
protected:
    nnvec_t(uint32_t L) : flat_vec<uint16_t, FE_nnvec>(L) {}

public:
#if TRACK_SRC
	static nnvec_t *newflat_src(uint32_t L,
		const char *srcfile, int srcline)
		{
		nnvec_t *p = new nnvec_t(L);
		p->m_srcfile = srcfile;
		p->m_srcline = srcline;
		return p;
		}
#else
	static nnvec_t *newflat(uint32_t L)
		{
		nnvec_t *p = new nnvec_t(L);
		return p;
		}
#endif
	};
