#pragma once

// Created by Igor Tolstoy 2024-06-19.
// Updated by Igor Tolstoy & Robert Edgar

#include "mpcflat.h"
#include "mega.h"

class MPCFlat_mega : public MPCFlat
    {
public:
    Mega *m_MM = 0;
	const vector<vector<vector<byte> > *> *m_ProfilePtrVec = 0;

private:
	void Run(MultiSequence *InputSeqs) { Die("MPCFlat_mega::Run()"); }

public:
	void Run(MultiSequence *InputSeqs,
	  const vector<vector<vector<byte> > *> &m_ProfilePtrVec);

protected:
    virtual void CalcPosterior(uint PairIndex);

protected:
	const vector<vector<byte> > &GetProfile(uint SMI) const;
    };
