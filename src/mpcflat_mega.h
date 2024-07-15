#pragma once

// Created by Igor Tolstoy 2024-06-19.
// Updated by Igor Tolstoy & Robert Edgar

#include "mpcflat.h"
#include "mega.h"

class MPCFlat_mega : public MPCFlat
    {
public:
	virtual void CalcFwdFlat_MPCFlat(uint GSIX, uint LX,
	  uint GSIY, uint LY, float *Flat);
	virtual void CalcBwdFlat_MPCFlat(uint GSIX, uint LX,
	  uint GSIY, uint LY, float *Flat);
    };
