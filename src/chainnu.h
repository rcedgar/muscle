#pragma once

#include "flat_base.h"

class chainnu_t : public flat_vec<uint8_t, FE_chainnu>
	{
protected:
	chainnu_t(uint32_t L) : flat_vec<uint8_t, FE_chainnu>(L) {}

public:
#if TRACK_SRC
	static chainnu_t *newflat_src(uint32_t L,
		const char *srcfile, int srcline)
		{
		chainnu_t *p = new chainnu_t(L);
		p->m_srcfile = srcfile;
		p->m_srcline = srcline;
		return p;
		}
#else
	static chainnu_t *newflat(uint32_t L)
		{
		chainnu_t *p = new chainnu_t(L);
		return p;
		}
#endif
	};
