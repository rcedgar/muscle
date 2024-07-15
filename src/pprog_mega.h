#pragma once

#include "pprog.h"
#include "mega.h"

class PProg_mega : public PProg
	{
// Members which know about HMM are virtual to
//  allow subclass override in PProg_mega.
public:
	virtual void CalcFwdFlat_PProg(uint GSI1, uint L1, 
	  uint GSI2, uint L2, float *Flat);

	virtual void CalcBwdFlat_PProg(uint GSI1, uint L1, 
	  uint GSI2, uint L2, float *Flat);
	};
