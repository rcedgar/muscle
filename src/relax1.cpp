#include "myutils.h"
#include "probcons.h"

/////////////////////////////////////////////////////////////////
// Relax1()
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////
void Relax1(SparseMatrix* matZX, SparseMatrix* matZY, vector<float>& posterior)
	{
	asserta(matZX != 0);
	asserta(matZY != 0);

	int lengthZ = matZX->GetSeq1Length();
	int lengthY = matZY->GetSeq2Length();

// for every z[k]
	for (int k = 1; k <= lengthZ; k++)
		{
		vector<PIF>::iterator ZXptr = matZX->GetRowPtr(k);
		vector<PIF>::iterator ZXend = ZXptr + matZX->GetRowSize(k);

	// iterate through all z[k]-x[i]
		while (ZXptr != ZXend)
			{
			vector<PIF>::iterator ZYptr = matZY->GetRowPtr(k);
			vector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(k);
			const float ZXval = ZXptr->second;
			vector<float>::iterator base = posterior.begin() +
			  ZXptr->first * (lengthY + 1);

		// iterate through all z[k]-y[j]
			while (ZYptr != ZYend)
				{
				base[ZYptr->first] += ZXval * ZYptr->second;
				ZYptr++;
				}
			ZXptr++;
			}
		}
	}
