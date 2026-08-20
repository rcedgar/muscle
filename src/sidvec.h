#pragma once

#include "flat_base.h"

class sidvec_t : public flat_vec<sid_t, FE_sidvec>
	{
protected:
    sidvec_t(uint32_t L) : flat_vec<sid_t, FE_sidvec>(L) {}

public:
#if TRACK_SRC
	static sidvec_t *newflat_src(uint32_t L,
		const char *srcfile, int srcline)
		{
		sidvec_t *p = new sidvec_t(L);
		p->m_srcfile = srcfile;
		p->m_srcline = srcline;
		return p;
		}
#else
	static sidvec_t *newflat(uint32_t L)
		{
		sidvec_t *p = new sidvec_t(L);
		return p;
		}
#endif
	};
