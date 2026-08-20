#pragma once

#include "flat_chain.h"
#include "linereader2.h"
#include "pdbfilescanner.h"
#include "bcadata.h"
#include "flat_dist_types.h"
#include <list>

struct chaq_vecs2;

class flat_chain_reader
	{
public:
	enum STATE
		{
		STATE_Closed,
		STATE_PendingFile,
		STATE_ReadingCALFile,
		STATE_ReadingCANFile,
		STATE_ReadingBCAFile,
		STATE_ReadingPDBFile,
		STATE_ReadingCIFFile,
		};

public:
	mutex m_CRPerThreadLock;
	mutex m_CRGlobalLock;
	STATE m_State = STATE_Closed;
	FILE *m_fCal = 0;
	LineReader2 m_LR;
	string m_Line;
	vector<string> m_Lines;
	vector<flat_chain_t *> m_Chains_PDB;
	vector<flat_chain_t *> m_Chains_CIF;
	vector<flat_chain_t *> *m_ptrChains = 0;
	uint m_ChainIdx_PDB = 0;
	uint m_ChainIdx_CIF = 0;
	uint m_ChainIdx_Vec = 0;
	string m_Label_PDB;
	BCAData m_BCA;
	uint64 m_ChainIdx_BCA = 0;
	string m_CurrentFN;
	bool m_Trace = false;
	bool m_SaveLines = false;
// When true, GetNext() fills chain nu codes (BCB read or chaq). Default
// false so MSA / structure-load paths do not pull chaq into the read loop.
	bool m_ComputeNu = false;

	sid_t *m_distmx = 0;
	uint8_t *m_codeseq_nu_scratch = 0;
	chaq_vecs2 *m_cv = 0;
	bool m_NuScratchInited = false;
	uint64 m_LastBCAChainIdx = UINT64_MAX;

// FS object shared with other threads
	PDBFileScanner *m_ptrFS = 0;

public:
	static uint m_CRGlobalChainCount;
	static uint m_CRGlobalFormatErrors;

public:
	void Open(const string &FileName);
	void Open(PDBFileScanner &FS);
	void Open(vector<flat_chain_t *> &Chains);
	flat_chain_t *GetNext();

private:
	void Close();
	flat_chain_t *GetNextLo1();
	flat_chain_t *GetFirst(const string &FN);

	flat_chain_t *GetFirst_BCA(const string &FN);
	flat_chain_t *GetNext_BCA();

	flat_chain_t *GetFirst_CAL(const string &FN);
	flat_chain_t *GetNext_CAL();

	flat_chain_t *GetFirst_CAN(const string &FN);
	flat_chain_t *GetNext_CAN();

	flat_chain_t *GetFirst_PDB(const string &FN);
	flat_chain_t *GetNext_PDB();

	flat_chain_t *GetFirst_CIF(const string &FN);
	flat_chain_t *GetNext_CIF();

	flat_chain_t *GetNext_Vec();

	void ChainsFromLines_PDB(const vector<string> &Lines,
		vector<flat_chain_t *> &Chains, const string &Label) const;
	void ChainsFromLines_CIF(const vector<string> &Lines,
		vector<flat_chain_t *> &Chains, const string &FallbackLabel);
	bool IsATOMLine_PDB(const string &Line) const;
	bool IsChainEndLine_PDB(const string &Line) const;
	uint GetCIFFieldIdx(const map<string, uint> &FieldToIdx, const string &Name);
	void IncFormatErrors();
	void InitNuScratch();
	void FreeNuScratch();
	void CacheNuOnChain(flat_chain_t *chain);
	static uint8_t ParseNuHexField(const string &hex, const string &FN,
		const string &line);
	};

void GetFallbackLabelFromFN(const string &FN, string &Label);
