#pragma once

#include "flat_base.h"

class chainxyz_t : public flat_mx<ic_t, FE_chainxyz>
	{
protected:
	chainxyz_t(uint32_t rows, uint32_t cols)
        : flat_mx<ic_t, FE_chainxyz>(rows, cols)
		{
		}

public:
#if TRACK_SRC
	static chainxyz_t *newflat_src(uint32_t L,
		const char *srcfile, int srcline)
		{
		chainxyz_t *p = new chainxyz_t(L, 3);
		p->m_srcfile = srcfile;
		p->m_srcline = srcline;
		return p;
		}
#else
	static chainxyz_t *newflat(uint32_t L)
		{
		chainxyz_t *p = new chainxyz_t(L, 3);
		return p;
		}
#endif
	};
