#ifndef tracebit_h
#define tracebit_h

#include "xdpmem.h"

static const byte TRACEBITS_DM = 0x01;
static const byte TRACEBITS_IM = 0x02;
static const byte TRACEBITS_MD = 0x04;
static const byte TRACEBITS_MI = 0x08;
static const byte TRACEBITS_SM = 0x10;
static const byte TRACEBITS_UNINIT = ~0x1f;

void TraceBackBitSW(XDPMem &Mem,
  uint LA, uint LB, uint Besti, uint Bestj,
  uint &Leni, uint &Lenj, string &Path);
void LogTBSW(const char *Msg, XDPMem &Mem, uint LA, uint LB);

#endif // tracebit_h
