#include "myutils.h"
#include "probcons.h"

/////////////////////////////////////////////////////////////////
// Relax()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

void Relax (SparseMatrix *matXZ, SparseMatrix *matZY, vector<float> &posterior){

  assert (matXZ);
  assert (matZY);

  int lengthX = matXZ->GetSeq1Length();
  int lengthY = matZY->GetSeq2Length();
  assert (matXZ->GetSeq2Length() == matZY->GetSeq1Length());

  // for every x[i]
  for (int i = 1; i <= lengthX; i++){
    vector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
    vector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

    vector<float>::iterator base = posterior.begin() + i * (lengthY + 1);

    // iterate through all x[i]-z[k]
    while (XZptr != XZend){
      vector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
      vector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(XZptr->first);
      const float XZval = XZptr->second;

      // iterate through all z[k]-y[j]
      while (ZYptr != ZYend){
        base[ZYptr->first] += XZval * ZYptr->second;
        ZYptr++;
      }
      XZptr++;
    }
  }
}
