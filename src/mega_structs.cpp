#include "muscle.h"
#include "mega.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "flat_chain.h"
#include "flat_chain_reader.h"
#include "pdbfilescanner.h"
#include "chaq.h"
#include "alpha.h"

bool Mega::IsStructsInput(const string &FileName)
	{
	if (FileName.empty())
		return false;
	if (IsDirectory(FileName))
		return true;

	string Ext;
	GetExtFromPathName(FileName, Ext);
	ToLower(Ext);
	if (Ext == "pdb" || Ext == "ent" || Ext == "cif" || Ext == "mmcif" ||
	  Ext == "cal" || Ext == "can" || Ext == "bca" || Ext == "bcb" ||
	  Ext == "files")
		return true;

	// Bare chain lists / scanners also accept paths without these
	// extensions when -structs is set (handled by caller).
	return false;
	}

static void CopyLogOddsMx(const float *FlatMx, uint AlphaSize,
  vector<vector<float> > &Mx)
	{
	Mx.clear();
	Mx.resize(AlphaSize);
	for (uint i = 0; i < AlphaSize; ++i)
		{
		Mx[i].resize(AlphaSize);
		for (uint j = 0; j < AlphaSize; ++j)
			Mx[i][j] = FlatMx[i*AlphaSize + j];
		}
	}

static void AppendChainProfile(const flat_params &params,
  const flat_chain_t *chain, const uint8_t *mega_prof)
	{
	const uint L = chain->get_length();
	const uint nfeat = params.m_nfeat;
	asserta(L > 0);
	asserta(nfeat > 0);

	const string &Label = chain->m_label;
	string TruncLabel = Label;
	trunc_label(TruncLabel);
	if (Mega::m_LabelToIdx.find(TruncLabel) != Mega::m_LabelToIdx.end())
		Die("Duplicate label in STRUCTS >%s", TruncLabel.c_str());

	const uint ProfileIdx = SIZE(Mega::m_Profiles);
	Mega::m_LabelToIdx[TruncLabel] = ProfileIdx;
	Mega::m_Labels.push_back(TruncLabel);

	vector<vector<byte> > Profile;
	Profile.resize(L);
	string Seq;
	Seq.reserve(L);

	for (uint Pos = 0; Pos < L; ++Pos)
		{
		vector<byte> &Col = Profile[Pos];
		Col.resize(nfeat);
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			byte Letter = mega_prof[fi*L + Pos];
			asserta(Letter < params.m_alpha_sizes[fi]);
			Col[fi] = Letter;
			}
		char aa = chain->get_aa(Pos);
		Seq.push_back(aa);
		}

	Mega::m_SeqToIdx[Seq] = ProfileIdx;
	Mega::m_Profiles.push_back(Profile);
	Mega::m_Seqs.push_back(Seq);
	}

void Mega::FromStructs(const string &FileName)
	{
	if (FileName == "")
		Die("Missing STRUCTS filename");
	RejectLegacyMega(FileName);
	if (m_Loaded)
		Die("Mega already loaded");

	m_Loaded = true;
	m_FileName = FileName;

	flat_params params;
	params.init_from_varstr("=sf");

	m_FeatureCount = params.m_nfeat;
	asserta(m_FeatureCount > 0);
	m_FeatureNames = params.m_alpha_names;
	m_AlphaSizes.resize(m_FeatureCount);
	m_Weights.resize(m_FeatureCount);
	m_LogOddsMxVec.resize(m_FeatureCount);
	m_LogProbsVec.resize(m_FeatureCount);

	for (uint fi = 0; fi < m_FeatureCount; ++fi)
		{
		const uint AS = params.m_alpha_sizes[fi];
		m_AlphaSizes[fi] = AS;
		m_Weights[fi] = params.m_weights[fi];
		CopyLogOddsMx(params.m_unweighted_logoddsvec[fi], AS,
		  m_LogOddsMxVec[fi]);
		m_LogProbsVec[fi].assign(AS, 0.0f);
		}

	m_GapOpen = params.m_open;
	m_GapExt = params.m_ext;

	const uint M = params.m_distmx_bandwidth;
	asserta(M > 0);

	PDBFileScanner FS;
	FS.Open(FileName);

	flat_chain_reader CR;
	CR.m_ComputeNu = false;
	CR.Open(FS);

	chaq_vecs2 cv;
	sid_t *distmx = 0;
	uint distmx_L = 0;
	uint8_t *scratch = 0;
	uint scratch_bytes = 0;
	uint8_t *mega_prof = 0;
	uint mega_prof_bytes = 0;

	uint ChainCount = 0;
	for (;;)
		{
		flat_chain_t *chain = CR.GetNext();
		if (chain == 0)
			break;

		const uint L = chain->get_length();
		if (L == 0)
			{
			delete chain;
			continue;
			}
		if (L > flat_params::m_maxL)
			Die("Chain '%s' length %u exceeds maxL %u",
			  chain->m_label.c_str(), L, flat_params::m_maxL);

		if (L > distmx_L)
			{
			myfree(distmx);
			distmx = myalloc(sid_t, L*M);
			distmx_L = L;

			if (cv.maxL > 0)
				chaq::free_chaq_vecs2(cv);
			chaq::alloc_chaq_vecs2(cv, L);

			const uint need_scratch = uint(
			  chaq::get_fast_get_codeseq_scratch_bytes_per_pos()*L);
			if (need_scratch > scratch_bytes)
				{
				myfree(scratch);
				scratch = myalloc(uint8_t, need_scratch);
				scratch_bytes = need_scratch;
				}
			}

		const uint need_prof = params.m_nfeat*L;
		if (need_prof > mega_prof_bytes)
			{
			myfree(mega_prof);
			mega_prof = myalloc(uint8_t, need_prof);
			mega_prof_bytes = need_prof;
			}

		chaq::fill_distmx(chain, distmx);
		chaq::fill_mega_prof(params, chain, distmx, mega_prof,
		  &cv, scratch, scratch_bytes);
		AppendChainProfile(params, chain, mega_prof);
		++ChainCount;
		delete chain;
		}

	if (cv.maxL > 0)
		chaq::free_chaq_vecs2(cv);
	myfree(distmx);
	myfree(scratch);
	myfree(mega_prof);

	if (ChainCount == 0)
		Die("No chains found in STRUCTS '%s'", FileName.c_str());
	ProgressLog("STRUCTS %s: %u chains, %u features (sf)\n",
	  FileName.c_str(), ChainCount, m_FeatureCount);
	}
