#ifndef MSADist_h
#define MSADist_h

#include <math.h>

double GetScoreDist(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2);

class MSADist
	{
public:
	MSADist(DISTANCE Distance)
		{
		m_Distance = Distance;
		}

	double ComputeDist(const MSA &msa, unsigned uSeqIndex1, unsigned uSeqIndex2)
		{
		if (m_Distance == DISTANCE_ScoreDist)
			return GetScoreDist(msa, uSeqIndex1, uSeqIndex2);

		double dPctId = msa.GetPctIdentityPair(uSeqIndex1, uSeqIndex2);
		switch(m_Distance)
			{
		case DISTANCE_PctIdKimura:
			return KimuraDist(dPctId);
		case DISTANCE_PctIdLog:
			if (dPctId < 0.05)
				dPctId = 0.05;
			return -log(dPctId);
			}
		Quit("MSADist::ComputeDist, invalid DISTANCE_%u", m_Distance);
		return 0;
		}

private:
	DISTANCE m_Distance;
	};

#endif	// MSADist_h
