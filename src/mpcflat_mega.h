#pragma once

// Created by Igor Tolstoy 2024-06-19.
// Updated by Igor Tolstoy & Robert Edgar

#include "mpcflat.h"
#include "mega.h"

class MPCFlat_mega : public MPCFlat
    {
public:
    Mega *m_MM = 0;

public:
    virtual void CalcPosterior(uint PairIndex);
    };
