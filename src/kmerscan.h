#pragma once

typedef uint fn_OnKmer(uint32 Code, uint32 Pos, void *UserData);
void SyncmerScanNt(const byte *Seq, uint Lo, uint Len, uint k, uint d,
  fn_OnKmer OnKmer, void *UserData);
void SyncmerScanAa(const byte *Seq, uint Lo, uint Len, uint k, uint d,
  fn_OnKmer OnKmer, void *UserData);

uint32 GetKmerMaskNt(uint k);
uint32 GetKmerMaskAa(uint k);

const char *CodeToStrNt(uint32 Code, uint k);
const char *CodeToStrAa(uint32 Code, uint k);

uint32 SeqToCodeNt(const byte *Seq, uint k);
uint32 SeqToCodeAa(const byte *Seq, uint k);
uint32 SeqToCode(bool Nucleo, const byte *Seq, uint k);

static const uint BITS_PER_LETTER_SEB8 = 3;
static const uint ALPHA_SIZE_SEB8 = 8;
