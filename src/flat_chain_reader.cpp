#include "myutils.h"
#include "flat_chain.h"
#include "flat_chain_reader.h"
#include "struct_desc.h"
#include "chaq.h"
#include "flat_params.h"

uint flat_chain_reader::m_CRGlobalChainCount;
uint flat_chain_reader::m_CRGlobalFormatErrors;

void flat_chain_reader::InitNuScratch()
	{
	if (m_NuScratchInited)
		return;
	asserta(m_ComputeNu);
	const uint M = flat_params::m_distmx_bandwidth;
	m_distmx = myalloc(sid_t, flat_params::m_maxL*M);
	m_codeseq_nu_scratch = myalloc(uint8_t, flat_params::m_maxL);
	asserta(m_cv == 0);
	m_cv = new chaq_vecs2;
	chaq::alloc_chaq_vecs2(*m_cv, flat_params::m_maxL);
	m_NuScratchInited = true;
	}

void flat_chain_reader::FreeNuScratch()
	{
	if (!m_NuScratchInited)
		return;
	myfree(m_distmx);
	myfree(m_codeseq_nu_scratch);
	if (m_cv != 0)
		{
		chaq::free_chaq_vecs2(*m_cv);
		delete m_cv;
		m_cv = 0;
		}
	m_distmx = 0;
	m_codeseq_nu_scratch = 0;
	m_NuScratchInited = false;
	}

uint8_t flat_chain_reader::ParseNuHexField(
	const string &hex, const string &FN, const string &line)
	{
	if (hex.size() != 2)
		Die("%s: Expected 2-digit hex, got '%s' in '%s'",
		  FN.c_str(), hex.c_str(), line.c_str());
	char *endptr = 0;
	long nu = strtol(hex.c_str(), &endptr, 16);
	if (endptr == hex.c_str())
		Die("%s: Expected hex digits, got '%s' in '%s'",
		  FN.c_str(), hex.c_str(), line.c_str());
	if (nu < 0 || nu >= 256)
		Die("%s: Expected hex digits in [0,255], got '%s'=%ld in '%s'",
		  FN.c_str(), hex.c_str(), nu, line.c_str());
	uint8_t nu_code = uint8_t(nu);
	assert(long(nu_code) == nu);
	return nu_code;
	}

void flat_chain_reader::CacheNuOnChain(flat_chain_t *chain)
	{
	if (chain->has_nu())
		return;
	const uint L = chain->get_length();
	if (L == 0)
		return;
	asserta(L <= flat_params::m_maxL);
	InitNuScratch();
	if (m_State == STATE_ReadingBCAFile && m_BCA.m_HasNuSequences)
		{
		asserta(m_LastBCAChainIdx != UINT64_MAX);
		uint nL = m_BCA.read_codeseq_nu(m_codeseq_nu_scratch,
			uint(m_LastBCAChainIdx), flat_params::m_maxL);
		asserta(nL == L);
		chain->set_nu_codes(m_codeseq_nu_scratch, L);
		return;
		}
	asserta(m_cv != 0);
	chaq::fill_codeseq_nu_from_chain(
		chain, m_distmx, m_cv,
		m_codeseq_nu_scratch, flat_params::m_maxL);
	chain->set_nu_codes(m_codeseq_nu_scratch, L);
	}

void flat_chain_reader::Close()
	{
	m_CRGlobalLock.lock();
	if (m_Trace) Log("flat_chain_reader::Close()\n");
	if (m_State != STATE_Closed)
		{
		m_State = STATE_Closed;
		FreeNuScratch();
		if (m_ptrFS != 0)
			delete m_ptrFS;
		m_ptrFS = 0;
		}
	m_CRGlobalLock.unlock();
	}

void flat_chain_reader::Open(const string &FileName)
	{
	asserta(m_State == STATE_Closed);
	asserta(m_ptrFS == 0);
	PDBFileScanner *FS = new PDBFileScanner;
	FS->Open(FileName);
	Open(*FS);
	}

void flat_chain_reader::Open(PDBFileScanner &FS)
	{
	asserta(m_State == STATE_Closed);
	m_ptrFS = &FS;
	m_Trace = opt(trace_chainreader2);
	m_ptrFS->m_Trace = opt(trace_chainreader2);
	if (m_Trace) Log("flat_chain_reader::Open()\n");
	m_State = STATE_PendingFile;
	m_CRGlobalChainCount = 0;
	}

void flat_chain_reader::Open(vector<flat_chain_t *> &Chains)
	{
	asserta(m_State == STATE_Closed);
	m_ptrChains = &Chains;
	m_ChainIdx_Vec = 0;
	}

// Files first, then directories to reduce queue
flat_chain_t* flat_chain_reader::GetFirst(const string &FN)
	{
	m_CurrentFN = FN;

	string Ext;
	GetExtFromPathName(FN, Ext);
	ToLower(Ext);
	if (m_Trace) Log("  FN=%s Ext=%s\n", FN.c_str(), Ext.c_str());
	if (Ext == "cal")
		{
		m_State = STATE_ReadingCALFile;
		flat_chain_t* Chain = GetFirst_CAL(FN);
		if (Chain)
			return Chain;
		}
	else if (Ext == "can")
		{
		m_State = STATE_ReadingCANFile;
		flat_chain_t* Chain = GetFirst_CAN(FN);
		if (Chain)
			return Chain;
		}
	else if (Ext == "bca" || Ext == "bcb")
		{
		m_State = STATE_ReadingBCAFile;
		flat_chain_t* Chain = GetFirst_BCA(FN);
		if (Chain)
			return Chain;
		}
	else if (Ext == "pdb" || Ext == "pdb.gz" || Ext == "ent" || Ext == "ent.gz")
		{
		m_State = STATE_ReadingPDBFile;
		flat_chain_t* Chain = GetFirst_PDB(FN);
		if (Chain)
			return Chain;
		}
	else if (Ext == "cif" || Ext == "cif.gz" || Ext == "mmcif" || Ext == "mmcif.gz")
		{
		m_State = STATE_ReadingCIFFile;
		flat_chain_t* Chain = GetFirst_CIF(FN);
		if (Chain)
			return Chain;
		}
	else
		Die("flat_chain_reader::GetNext(%s), unknown extension", FN.c_str());
	return 0;
	}

flat_chain_t* flat_chain_reader::GetNext()
	{
	for (uint SanityCounter = 0; ; ++SanityCounter)
		{
		if (SanityCounter > 100)
			Warning("Excessive looping in flat_chain_reader::GetNext()");

		m_CRPerThreadLock.lock();
		flat_chain_t* Chain = GetNextLo1();
		m_CRPerThreadLock.unlock();

		if (!Chain)
			return 0;

		uint L = Chain->get_length();
		if (L == 0)
			{
			delete Chain;
			continue;
			}
		if (L > flat_params::m_maxL)
			Chain->truncate(flat_params::m_maxL);
		if (m_ComputeNu)
			CacheNuOnChain(Chain);
		return Chain;
		}
	}

flat_chain_t* flat_chain_reader::GetNextLo1()
	{
	for (uint SanityCounter = 0; SanityCounter < 100; ++SanityCounter)
		{
		if (m_Trace) Log("GetNextLo1() state=%d\n", m_State);
		switch (m_State)
			{
		case STATE_Closed:
			return 0;

		case STATE_PendingFile:
			{
			string FN;
			bool Ok = m_ptrFS->GetNext(FN);
			if (!Ok)
				return 0;
			flat_chain_t* Chain = GetFirst(FN);
			if (Chain)
				return Chain;
			continue;
			}

		case STATE_ReadingCALFile:
			{
			flat_chain_t* Chain = GetNext_CAL();
			if (Chain)
				return Chain;
			if (m_Trace) Log("GetNext_CAL()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingCANFile:
			{
			flat_chain_t* Chain = GetNext_CAN();
			if (Chain)
				return Chain;
			if (m_Trace) Log("GetNext_CAN()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingBCAFile:
			{
			flat_chain_t* Chain = GetNext_BCA();
			if (Chain)
				return Chain;
			if (m_Trace) Log("GetNext_BCA()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingPDBFile:
			{
			flat_chain_t* Chain = GetNext_PDB();
			if (Chain)
				return Chain;
			if (m_Trace) Log("GetNext_PDB()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingCIFFile:
			{
			flat_chain_t* Chain = GetNext_CIF();
			if (Chain)
				return Chain;
			if (m_Trace) Log("GetNext_CIF()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		default:
			asserta(false);
			}
		}
	Die("Excessive looping in flat_chain_reader::GetNextLo1()");
	return 0;
	}

flat_chain_t* flat_chain_reader::GetFirst_BCA(const string &FN)
	{
	m_BCA.Open(FN);
	m_ChainIdx_BCA = 0;
	return GetNext_BCA();
	}

flat_chain_t* flat_chain_reader::GetNext_BCA()
	{
	uint64 ChainCount = m_BCA.GetChainCount();
	if (m_ChainIdx_BCA >= ChainCount)
		{
		m_BCA.Close();
		m_LastBCAChainIdx = UINT64_MAX;
		return 0;
		}
	uint64 idx = m_ChainIdx_BCA++;
	flat_chain_t* chain = m_BCA.read_flat_chain(idx);
	m_LastBCAChainIdx = idx;
	return chain;
	}

flat_chain_t* flat_chain_reader::GetFirst_CAL(const string &FN)
	{
	m_LR.Open(FN);
	bool Ok = m_LR.ReadLine(m_Line);
	if (!Ok)
		Die("Failed to read first line of CAL file '%s'",
		  FN.c_str());
	return GetNext_CAL();
	}

flat_chain_t* flat_chain_reader::GetNext_CAL()
	{
	if (m_LR.m_EOF)
		{
		m_LR.Close();
		return 0;
		}
	if (m_Line.empty() || m_Line[0] != '>')
		Die("%s: Expected '>' in CAL file",
		  m_CurrentFN.c_str());

	const string Label = m_Line.substr(1);
	if (m_Trace) Log("flat_chain_reader::GetNext_CAL() Label=%s\n", Label.c_str());
	m_Lines.clear();
	while (m_LR.ReadLine(m_Line))
		{
		if (m_Line.c_str()[0] == '>')
			break;
		m_Lines.push_back(m_Line);
		}

/***
>102l
M       43.619  -1.924  8.869
N       40.445  -0.876  10.670
I       38.254  2.240   11.220
F       40.340  3.621   14.036
***/
	const uint N = SIZE(m_Lines);
	vector<string> Fields;
	vector<char> aas;
	vector<float> Xs, Ys, Zs;
	aas.reserve(RESERVE_CHAIN_LENGTH);
	Xs.reserve(RESERVE_CHAIN_LENGTH);
	Ys.reserve(RESERVE_CHAIN_LENGTH);
	Zs.reserve(RESERVE_CHAIN_LENGTH);
	for (uint LineNr = 0; LineNr < N; ++LineNr)
		{
		const string &Line = m_Lines[LineNr];
		if (Line.empty())
			continue;
		Split(Line, Fields, '\t');
		if (Fields.size() != 4 || Fields[0].size() != 1)
			Die("%s: Invalid CAL record '%s'",
			  m_CurrentFN.c_str(), Line.c_str());

		char aa = Fields[0][0];
		float X = StrToFloatf(Fields[1]);
		float Y = StrToFloatf(Fields[2]);
		float Z = StrToFloatf(Fields[3]);

		aas.push_back(aa);
		Xs.push_back(X);
		Ys.push_back(Y);
		Zs.push_back(Z);
		}
	auto chain = flat_chain_t::newflat(Label, aas, Xs, Ys, Zs);
	return chain;
	}

flat_chain_t* flat_chain_reader::GetFirst_CAN(const string &FN)
	{
	m_LR.Open(FN);
	bool Ok = m_LR.ReadLine(m_Line);
	if (!Ok)
		Die("Failed to read first line of CAN file '%s'",
		  FN.c_str());
	return GetNext_CAN();
	}

flat_chain_t* flat_chain_reader::GetNext_CAN()
	{
	if (m_LR.m_EOF)
		{
		m_LR.Close();
		return 0;
		}
	if (m_Line.empty() || m_Line[0] != '>')
		Die("%s: Expected '>' in CAN file",
		  m_CurrentFN.c_str());

	const string Label = m_Line.substr(1);
	if (m_Trace) Log("flat_chain_reader::GetNext_CAN() Label=%s\n", Label.c_str());
	m_Lines.clear();
	while (m_LR.ReadLine(m_Line))
		{
		if (m_Line.c_str()[0] == '>')
			break;
		m_Lines.push_back(m_Line);
		}

	const uint N = SIZE(m_Lines);
	vector<string> Fields;
	vector<char> aas;
	vector<float> Xs, Ys, Zs;
	vector<uint8_t> nu_codes;
	aas.reserve(RESERVE_CHAIN_LENGTH);
	Xs.reserve(RESERVE_CHAIN_LENGTH);
	Ys.reserve(RESERVE_CHAIN_LENGTH);
	Zs.reserve(RESERVE_CHAIN_LENGTH);
	nu_codes.reserve(RESERVE_CHAIN_LENGTH);
	for (uint LineNr = 0; LineNr < N; ++LineNr)
		{
		const string &Line = m_Lines[LineNr];
		if (Line.empty())
			continue;
		Split(Line, Fields, '\t');
		if (Fields.size() != 5 || Fields[0].size() != 1)
			Die("%s: Invalid CAN record '%s'",
			  m_CurrentFN.c_str(), Line.c_str());

		char aa = Fields[0][0];
		float X = StrToFloatf(Fields[1]);
		float Y = StrToFloatf(Fields[2]);
		float Z = StrToFloatf(Fields[3]);
		uint8_t nu_code = ParseNuHexField(Fields[4], m_CurrentFN, Line);

		aas.push_back(aa);
		Xs.push_back(X);
		Ys.push_back(Y);
		Zs.push_back(Z);
		nu_codes.push_back(nu_code);
		}
	auto chain = flat_chain_t::newflat(Label, aas, Xs, Ys, Zs);
	if (!nu_codes.empty())
		chain->set_nu_codes(nu_codes.data(), SIZE(nu_codes));
	return chain;
	}

flat_chain_t* flat_chain_reader::GetFirst_PDB(const string &FN)
	{
	ReadLinesFromFile(FN, m_Lines);
	string Label;
	GetFallbackLabelFromFN(FN, Label);
	ChainsFromLines_PDB(m_Lines, m_Chains_PDB, Label);
	m_ChainIdx_PDB = 0;
	return GetNext_PDB();
	}

flat_chain_t* flat_chain_reader::GetFirst_CIF(const string &FN)
	{
	ReadLinesFromFile(FN, m_Lines);

	string FallbackLabel;
	GetFallbackLabelFromFN(FN, FallbackLabel);
	ChainsFromLines_CIF(m_Lines, m_Chains_CIF, FallbackLabel);
	m_ChainIdx_CIF = 0;
	return GetNext_CIF();
	}

flat_chain_t* flat_chain_reader::GetNext_PDB()
	{
	const uint N = SIZE(m_Chains_PDB);
	if (m_ChainIdx_PDB == N)
		return 0;
	asserta(m_ChainIdx_PDB < N);
	flat_chain_t* Chain = m_Chains_PDB[m_ChainIdx_PDB++];
	if (m_Trace) Log("flat_chain_reader::GetNext_PDB() %u/%u Label=%s\n", m_ChainIdx_PDB, N, Chain->m_label.c_str());
	return Chain;
	}

flat_chain_t* flat_chain_reader::GetNext_Vec()
	{
	asserta(m_ptrChains != 0);
	const uint N = SIZE(*m_ptrChains);
	if (m_ChainIdx_Vec == N)
		return 0;
	asserta(m_ChainIdx_CIF < N);
	flat_chain_t* Chain = (*m_ptrChains)[m_ChainIdx_Vec++];
	if (m_Trace) Log("flat_chain_reader::GetNext_Vec() %u/%u Label=%s\n", m_ChainIdx_Vec, N, Chain->m_label.c_str());
	return Chain;
	}

flat_chain_t* flat_chain_reader::GetNext_CIF()
	{
	const uint N = SIZE(m_Chains_CIF);
	if (m_ChainIdx_CIF == N)
		return 0;
	asserta(m_ChainIdx_CIF < N);
	flat_chain_t* Chain = m_Chains_CIF[m_ChainIdx_CIF++];
	if (m_Trace) Log("flat_chain_reader::GetNext_CIF() %u/%u Label=%s\n", m_ChainIdx_CIF, N, Chain->m_label.c_str());
	return Chain;
	}

bool flat_chain_reader::IsATOMLine_PDB(const string &Line) const
	{
	if (SIZE(Line) < 27)
		return false;
	if (strncmp(Line.c_str(), "ATOM  ", 6) == 0)
		return true;
	return false;
	}

bool flat_chain_reader::IsChainEndLine_PDB(const string &Line) const
	{
	// ENDMDL ends the first model (later models ignored).
	// Do not treat TER as chain end: mid-chain TER is common before
	// HETATM chromophores (e.g. GFP CRO), with the same chain ID
	// continuing afterward. Chain splits are by chain ID only.
	if (StartsWith(Line, "ENDMDL"))
		return true;
	return false;
	}

void flat_chain_reader::ChainsFromLines_PDB(const vector<string> &Lines,
  vector<flat_chain_t *> &Chains, const string &Label) const
	{
	Chains.clear();

	string Entry = Label;
	string Title;
	map<string, string> MolByChain;
	map<string, vector<string> > RefsByChain;
	ExtractPdbMeta(Lines, Entry, Title, MolByChain, RefsByChain);

	const uint N = SIZE(Lines);
	vector<string> ChainLines;
	char CurrChainChar = 0;
	bool AnyAtoms = false;
	bool EndOfChainFound = false;

	auto FinishChain = [&]()
		{
		if (!AnyAtoms || ChainLines.empty())
			return;
		flat_chain_t* Chain = flat_chain_t::newflat(0);
		bool Ok = Chain->from_pdb_lines(Entry, ChainLines, m_SaveLines);
		if (Ok)
			{
			string ChainStr;
			ChainStr.push_back(CurrChainChar);
			if (ChainStr == " " || ChainStr.empty())
				ChainStr = "_";
			string db_ref;
			map<string, vector<string> >::const_iterator rit =
				RefsByChain.find(ChainStr);
			if (rit != RefsByChain.end())
				db_ref = PreferUnp(rit->second);
			string molecule = PickMoleculeByChain(MolByChain, ChainStr);
			AppendStructDescToLabel(Chain->m_label, Entry, db_ref,
			  molecule, Title);
			Chains.push_back(Chain);
			}
		else
			delete Chain;
		ChainLines.clear();
		EndOfChainFound = false;
		AnyAtoms = false;
		};

	for (uint i = 0; i < N; ++i)
		{
		const string &Line = Lines[i];
		if (IsChainEndLine_PDB(Line))
			EndOfChainFound = true;
		if (IsATOMLine_PDB(Line))
			{
			if (Line.size() < 54)
				continue;
			char ChainChar = Line[21];
			if (ChainChar != CurrChainChar)
				{
				FinishChain();
				CurrChainChar = ChainChar;
				}
			if (!EndOfChainFound)
				ChainLines.push_back(Line);
			AnyAtoms = true;
			}
		}

	FinishChain();
	}
