#pragma once

#include "flat_base.h"

class chainaa_t : public flat_vec<char, FE_chainaa>
	{
protected:
    chainaa_t(uint32_t L) : flat_vec<char, FE_chainaa>(L) {}

public:
#if TRACK_SRC
	static chainaa_t *newflat_src(uint32_t L,
		const char *srcfile, int srcline)
		{
		chainaa_t *p = new chainaa_t(L);
		p->m_srcfile = srcfile;
		p->m_srcline = srcline;
		return p;
		}
#else
	static chainaa_t *newflat(uint32_t L)
		{
		chainaa_t *p = new chainaa_t(L);
		return p;
		}
#endif
	};
