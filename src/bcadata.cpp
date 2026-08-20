#include "myutils.h"
#include "flat_chain.h"
#include "flat_params.h"
#include "bcadata.h"
#include "chaq.h"

void BCAData::Close()
	{
	if (m_Reading && !m_Writing)
		CloseReader();
	else if (m_Writing && !m_Reading)
		CloseWriter();
	else
		Die("BCAData::Close(), not open");
	}

void BCAData::Create(const string &FN, bool WithNu)
	{
	if (FN == "")
		Die("Empty BCA filename");
	asserta(!m_Writing && !m_Reading);
	m_HasNuSequences = WithNu;
	m_FN = FN;
	m_f = CreateStdioFile(FN);
	const uint Magic = (WithNu ? BCB_MAGIC : BCA_MAGIC);
	WriteStdioFile(m_f, &Magic, sizeof(Magic));

// Placeholder #1 overwritten with number of chains
// Placeholder #2 overwritten with address of labels in Close()
// Placeholder #3 overwritten with size of labels data
	uint64_t Placeholder = 0;
	WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	if (WithNu)
		{
// Placeholder #4 overwritten with start of contiguous nu section (BCB only)
		WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
		}
	m_Writing = true;

// nu bytes are streamed to a temp file during writing, then appended as
// one contiguous section in CloseWriter (keeps peak memory bounded).
	if (m_HasNuSequences)
		{
		m_NuTmpFN = FN + ".nutmp";
		m_nu_tmp_f = CreateStdioFile(m_NuTmpFN);
		}

	const uint N = 1024*1024;
	m_Labels.reserve(N);
	m_Offsets.reserve(N);
	m_SeqLengths.reserve(N);
	}

uint64 BCAData::get_offset_aaseq(uint idx) const
	{
	asserta(idx < m_Offsets.size());
	uint64 offset = m_Offsets[idx];
	return offset;
	}

// Contiguous nu layout: nu records are packed back-to-back in chain-index
// order in a section starting at m_NuSeqPos64. Offset of chain idx is the
// section base plus the sum of lengths of all preceding chains.
uint64 BCAData::get_offset_nuseq(uint idx) const
	{
	asserta(m_HasNuSequences);
	asserta(m_NuSeqPos64 != UINT64_MAX);
	asserta(idx < m_NuPrefix.size());
	return m_NuSeqPos64 + m_NuPrefix[idx];
	}

uint64 BCAData::get_offset_ICs(uint idx) const
	{
	asserta(idx < m_SeqLengths.size());
	uint64 offset_aaseq = get_offset_aaseq(idx);
	uint L = m_SeqLengths[idx];
	uint64 offset_ICs = offset_aaseq + L;
	return offset_ICs;
	}

void BCAData::write_flat_chain(const flat_chain_t *chain, chaq_vecs2 *cv)
	{
	asserta(m_Writing && !m_Reading);
	uint L = chain->get_length();
	if (L == 0) return;
	uint64_t Offset = GetStdioFilePos64(m_f);
	size_t n = m_Offsets.size();
	asserta(m_SeqLengths.size() == n);
// nu is no longer interleaved, so every chain occupies 7*L bytes
// (AA + ICs) in the main stream regardless of m_HasNuSequences.
	if (n > 0)
		{
		uint Ln_1 = m_SeqLengths[n-1];
		asserta(Offset == m_Offsets[n-1] + 7*Ln_1);
		}
	const char *seq = chain->m_aa->m_data;
	uint Idx = SIZE(m_Labels);
	asserta(SIZE(m_SeqLengths) == Idx);

	m_Labels.push_back(chain->m_label);
	m_SeqLengths.push_back(L);
	m_Offsets.push_back(Offset);
	vector<uint16_t> ICs;
	chain->get_ICs(ICs);
	asserta(SIZE(ICs) == 3*L);
	assert(GetStdioFilePos64(m_f) == get_offset_aaseq(Idx));
	WriteStdioFile64(m_f, seq, L);
	assert(GetStdioFilePos64(m_f) == get_offset_ICs(Idx));
	WriteStdioFile64(m_f, ICs.data(), 6*L);
	if (m_HasNuSequences)
		append_codeseq_nu(chain, cv);
	}

void BCAData::append_codeseq_nu(
	const flat_chain_t *chain,
	chaq_vecs2 *cv)
	{
	if (m_distmx == 0)
		{
		m_distmx = myalloc(sid_t,
			flat_params::m_distmx_bandwidth*flat_params::m_maxL);
		m_codeseq_nu = myalloc(uint8_t, flat_params::m_maxL);
		chaq::alloc_chaq_vecs2(*cv, flat_params::m_maxL);
		}
	const uint L = chain->get_length();
	asserta(L <= flat_params::m_maxL); // TODO=maxL
	asserta(m_nu_tmp_f != 0);
	if (chain->has_nu())
		{
		WriteStdioFile64(m_nu_tmp_f, chain->get_nu_data(), L);
		return;
		}
	chaq::fill_codeseq_nu_from_chain(
		chain, m_distmx, cv, m_codeseq_nu, flat_params::m_maxL);
	WriteStdioFile64(m_nu_tmp_f, m_codeseq_nu, L);
	}

//void BCAData::WriteChain(const PDBChain &Chain)
//	{
//	asserta(m_Writing && !m_Reading);
//	asserta(!m_HasNuSequences);
//	uint64_t Offset = GetStdioFilePos64(m_f);
//	size_t n = m_Offsets.size();
//	asserta(m_SeqLengths.size() == n);
//	uint L = Chain.GetSeqLength();
//	if (n > 0)
//		{
//		uint Ln_1 = m_SeqLengths[n-1];
//		asserta(Offset == m_Offsets[n-1] + 7*Ln_1);
//		}
//	const string &Seq = Chain.m_Seq;
//	uint Idx = SIZE(m_Labels);
//	asserta(SIZE(m_SeqLengths) == Idx);
//
//	m_Labels.push_back(Chain.m_Label);
//	m_SeqLengths.push_back(L);
//	m_Offsets.push_back(Offset);
//	vector<uint16_t> ICs;
//	Chain.GetICs(ICs);
//	asserta(SIZE(ICs) == 3*L);
//	WriteStdioFile64(m_f, Seq.c_str(), L);
//	WriteStdioFile64(m_f, ICs.data(), 6*L);
//	}

void BCAData::Open(const string &FN)
	{
	if (FN == "")
		Die("Empty BCA filename");
	m_FN = FN;
	asserta(!m_Writing && !m_Reading);
	asserta(m_f == 0);

	m_f = OpenStdioFile(FN);

	//m_scratch_buffer_bytes = 2*m_maxL;
	//m_scratch_buffer = myalloc(uint8_t, m_scratch_buffer_bytes);
	//chaq::alloc_chaq_vecs2(m_cv, m_maxL);

	uint32_t Magic;
	ReadStdioFile(m_f, &Magic, sizeof(Magic));
	if (Magic == BCA_MAGIC)
		m_HasNuSequences = false;
	else if (Magic == BCB_MAGIC)
		m_HasNuSequences = true;
	else
		Die("Bad magic %08lx, invalid .bcx file '%s'",
		  Magic, FN.c_str());

// Placeholder #1 overwritten with number of chains
// Placeholder #2 overwritten with address of labels in Close()
// Placeholder #3 overwritten with size of labels data
// Placeholder #4 overwritten with start of contiguous nu section
	uint64_t ChainCount64;
	ReadStdioFile(m_f, &ChainCount64, sizeof(uint64_t));
	ReadStdioFile(m_f, &m_SeqLengthsPos64, sizeof(uint64_t));
	ReadStdioFile(m_f, &m_LabelDataSize64, sizeof(uint64_t));
	if (m_HasNuSequences)
		ReadStdioFile(m_f, &m_NuSeqPos64, sizeof(uint64_t));
	else
		m_NuSeqPos64 = UINT64_MAX;
	uint64 Offset = GetStdioFilePos64(m_f);

	uint ChainCount = uint(ChainCount64);
	asserta(ChainCount == ChainCount64);
	uint64 SeqLengthsBytes = sizeof(uint32_t)*ChainCount;

	m_SeqLengths.resize(ChainCount64);

	SetStdioFilePos64(m_f, m_SeqLengthsPos64);
	ReadStdioFile64NoPos(m_f, m_SeqLengths.data(), SeqLengthsBytes);

// AA+ICs occupy 7*L per chain in the main stream; nu (if present) lives
// in a separate contiguous section, indexed by m_NuPrefix.
	m_NuPrefix.resize(ChainCount64);
	uint64 NuPrefix = 0;
	for (uint64 i = 0; i < ChainCount64; ++i)
		{
		uint L = m_SeqLengths[i];
		m_Offsets.push_back(Offset);
		Offset += 7*uint64(L);
		m_NuPrefix[i] = NuPrefix;
		NuPrefix += L;
		}

	uint LabelDataSize = uint(m_LabelDataSize64);
	asserta(LabelDataSize == m_LabelDataSize64);
	char *LabelData = myalloc(char, LabelDataSize);
	ReadStdioFile(m_f, LabelData, LabelDataSize);
	m_Labels.clear();
	uint n = 0;
	m_Labels.push_back(LabelData);
	for (uint i = 0; i + 1 < LabelDataSize; ++i)
		if (LabelData[i] == 0)
			m_Labels.push_back(LabelData + i + 1);
	uint LabelCount = SIZE(m_Labels);
	if (LabelCount != ChainCount64)
		Die("Bad BCA file, %u chains %u labels",
		  ChainCount, LabelCount);

	m_Reading = true;
	}

void BCAData::Clear()
	{
	m_Labels.clear();
	m_Offsets.clear();
	m_SeqLengths.clear();
	m_NuPrefix.clear();
	m_FN.clear();
	if (m_f != 0)
		CloseStdioFile(m_f);
	m_f = 0;
	if (m_nu_tmp_f != 0)
		{
		CloseStdioFile(m_nu_tmp_f);
		m_nu_tmp_f = 0;
		}
	m_NuTmpFN.clear();
	m_Writing = false;
	m_Reading = false;
	m_SeqLengthsPos64 = UINT64_MAX;
	m_LabelDataSize64 = UINT64_MAX;
	m_NuSeqPos64 = UINT64_MAX;
	}

void BCAData::CloseReader()
	{
	asserta(m_Reading && !m_Writing);
	Clear();
	}

// Stream a write-then-read temp file's full contents to the end of m_f.
void BCAData::AppendTempFileToMain(FILE *tmp_f)
	{
	asserta(tmp_f != 0);
	fflush(tmp_f);
	SetStdioFilePos64(tmp_f, 0);
	const size_t BufBytes = 4*1024*1024;
	uint8_t *buf = myalloc(uint8_t, BufBytes);
	for (;;)
		{
		size_t n = fread(buf, 1, BufBytes, tmp_f);
		if (n == 0)
			break;
		WriteStdioFile64(m_f, buf, n);
		}
	myfree(buf);
	}

// Write the seq-length table and label data at the current file position,
// then rewind and patch the header placeholders. Assumes the main stream
// (and nu section / m_NuSeqPos64) are already in place.
void BCAData::WriteSeqLengthsLabelsHeader()
	{
	const uint ChainCount = GetChainCount();

	m_SeqLengthsPos64 = GetStdioFilePos64(m_f);
	WriteStdioFile64(m_f, m_SeqLengths.data(), sizeof(uint32_t)*ChainCount);

	m_LabelDataSize64 = 0;
	for (uint i = 0; i < ChainCount; ++i)
		{
		const string &Label = m_Labels[i];
		const uint n = SIZE(Label) + 1;
		WriteStdioFile(m_f, Label.c_str(), n);
		m_LabelDataSize64 += n;
		}

// Re-wind to overwrite placeholders in file header
	SetStdioFilePos64(m_f, sizeof(BCA_MAGIC));
	uint64 ChainCount64 = ChainCount;

// #1 number of chains
// #2 address of labels
// #3 size of labels data
// #4 start of contiguous nu section
	WriteStdioFile(m_f, &ChainCount64, sizeof(ChainCount64));
	WriteStdioFile(m_f, &m_SeqLengthsPos64, sizeof(m_SeqLengthsPos64));
	WriteStdioFile(m_f, &m_LabelDataSize64, sizeof(m_LabelDataSize64));
	if (m_HasNuSequences)
		WriteStdioFile(m_f, &m_NuSeqPos64, sizeof(m_NuSeqPos64));
	}

void BCAData::CloseWriter()
	{
	asserta(m_Writing && !m_Reading);

// The AA+ICs section is now complete; the contiguous nu section (if any)
// is appended here by streaming the temp file in.
	if (m_HasNuSequences)
		{
		asserta(m_nu_tmp_f != 0);
		m_NuSeqPos64 = GetStdioFilePos64(m_f);
		AppendTempFileToMain(m_nu_tmp_f);
		CloseStdioFile(m_nu_tmp_f);
		m_nu_tmp_f = 0;
		DeleteStdioFile(m_NuTmpFN);
		}
	else
		m_NuSeqPos64 = UINT64_MAX;

	WriteSeqLengthsLabelsHeader();
	Clear();
	}

// ---- Parallel sharded writing -------------------------------------------

void BCAData::CreateSharded(const string &FN, bool WithNu, uint nshard)
	{
	if (FN == "")
		Die("Empty BCA filename");
	asserta(!m_Writing && !m_Reading);
	asserta(nshard >= 1);
	m_HasNuSequences = WithNu;
	m_FN = FN;
	m_f = CreateStdioFile(FN);
	const uint Magic = (WithNu ? BCB_MAGIC : BCA_MAGIC);
	WriteStdioFile(m_f, &Magic, sizeof(Magic));

// Placeholders #1..#3 (and #4 for BCB), patched in WriteSeqLengthsLabelsHeader.
	uint64_t Placeholder = 0;
	WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	if (WithNu)
		WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	m_Writing = true;

	const uint Reserve = 1024*1024;
	m_Labels.reserve(Reserve);
	m_Offsets.reserve(Reserve);
	m_SeqLengths.reserve(Reserve);

	m_Shards.resize(nshard);
	for (uint s = 0; s < nshard; ++s)
		{
		Shard *sh = new Shard;
		sh->aa_fn = FN + ".aatmp." + std::to_string(s);
		sh->aa_f = CreateStdioFile(sh->aa_fn);
		if (WithNu)
			{
			sh->nu_fn = FN + ".nutmp." + std::to_string(s);
			sh->nu_f = CreateStdioFile(sh->nu_fn);
			}
		m_Shards[s] = sh;
		}
	}

void BCAData::write_flat_chain_shard(uint s, const flat_chain_t *chain)
	{
	asserta(m_Writing && !m_Reading);
	asserta(s < SIZE(m_Shards));
	Shard *sh = m_Shards[s];
	const uint L = chain->get_length();
	if (L == 0)
		return;
	const char *seq = chain->m_aa->m_data;

	sh->labels.push_back(chain->m_label);
	sh->seqlengths.push_back(L);

	vector<uint16_t> ICs;
	chain->get_ICs(ICs);
	asserta(SIZE(ICs) == 3*L);
	WriteStdioFile64(sh->aa_f, seq, L);
	WriteStdioFile64(sh->aa_f, ICs.data(), 6*L);

	if (m_HasNuSequences)
		{
		asserta(sh->nu_f != 0);
		asserta(L <= flat_params::m_maxL);
		if (chain->has_nu())
			WriteStdioFile64(sh->nu_f, chain->get_nu_data(), L);
		else
			{
			if (sh->distmx == 0)
				{
				sh->distmx = myalloc(sid_t,
					flat_params::m_distmx_bandwidth*flat_params::m_maxL);
				sh->codeseq_nu = myalloc(uint8_t, flat_params::m_maxL);
				chaq::alloc_chaq_vecs2(sh->cv, flat_params::m_maxL);
				sh->cv_inited = true;
				}
			chaq::fill_codeseq_nu_from_chain(
				chain, sh->distmx, &sh->cv,
				sh->codeseq_nu, flat_params::m_maxL);
			WriteStdioFile64(sh->nu_f, sh->codeseq_nu, L);
			}
		}
	}

void BCAData::CloseSharded()
	{
	asserta(m_Writing && !m_Reading);
	asserta(!m_Shards.empty());

// Main stream = AA+ICs of all shards concatenated in shard order, with the
// global label/length tables built in the same order.
	for (uint s = 0; s < SIZE(m_Shards); ++s)
		{
		Shard *sh = m_Shards[s];
		AppendTempFileToMain(sh->aa_f);
		const uint ns = SIZE(sh->seqlengths);
		asserta(SIZE(sh->labels) == ns);
		for (uint i = 0; i < ns; ++i)
			{
			m_SeqLengths.push_back(sh->seqlengths[i]);
			m_Labels.push_back(sh->labels[i]);
			}
		}

// Contiguous nu section, same shard order.
	if (m_HasNuSequences)
		{
		m_NuSeqPos64 = GetStdioFilePos64(m_f);
		for (uint s = 0; s < SIZE(m_Shards); ++s)
			AppendTempFileToMain(m_Shards[s]->nu_f);
		}
	else
		m_NuSeqPos64 = UINT64_MAX;

	WriteSeqLengthsLabelsHeader();

// Release shard temp files and scratch.
	for (uint s = 0; s < SIZE(m_Shards); ++s)
		{
		Shard *sh = m_Shards[s];
		if (sh->aa_f != 0)
			{
			CloseStdioFile(sh->aa_f);
			DeleteStdioFile(sh->aa_fn);
			}
		if (sh->nu_f != 0)
			{
			CloseStdioFile(sh->nu_f);
			DeleteStdioFile(sh->nu_fn);
			}
		myfree(sh->distmx);
		myfree(sh->codeseq_nu);
		if (sh->cv_inited)
			chaq::free_chaq_vecs2(sh->cv);
		delete sh;
		}
	m_Shards.clear();
	Clear();
	}

//uint64 BCAData::GetICsOffset(uint64 ChainIdx) const
//	{
//	asserta(ChainIdx < SIZE(m_Offsets));
//	asserta(ChainIdx < SIZE(m_SeqLengths));
//	uint64 Offset = m_Offsets[ChainIdx];
//	uint L = m_SeqLengths[ChainIdx];
//	return Offset + L;
//	}

uint64 BCAData::GetSeqOffset(uint64 ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_Offsets));
	return m_Offsets[ChainIdx];
	}

uint BCAData::GetSeqLength(uint64 ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_SeqLengths));
	return m_SeqLengths[ChainIdx];
	}

// convert flat x0,y0,z0, x1,y1,z1 ...
// to flat x0,x1 ... y0,y1 ... z0,y1
static inline void aos_to_soa_u16(
	const uint16_t* __restrict aos,
	uint16_t* __restrict soa,
	uint32_t L
){
	for (uint32_t i = 0; i < L; ++i) {
		soa[i] = aos[3*i + 0];
		soa[L + i] = aos[3*i + 1];
		soa[2*L + i] = aos[3*i + 2];
	}
}

uint BCAData::read_codeseq_nu(
	uint8_t *codeseq_nu, uint idx, uint buffer_length) const
	{
	asserta(m_Reading && !m_Writing);
	uint L = GetSeqLength(idx);
	if (L > buffer_length)
		L = buffer_length;
	uint64 offset = get_offset_nuseq(idx);
	uint64 nL = ReadStdioFile64_NoFail(m_f, offset, codeseq_nu, L);
	if (nL != L)
		{
		Log("FN=%s\n", m_FN.c_str());
		Log("ChainIdx=%u\n", idx);
		Log("Chains=%u\n", SIZE(m_SeqLengths));
		Log("L=%u\n", L);
		Log("SeqOffset=%llu\n", (unsigned long long) offset);
		Log("nL=%llu\n", (unsigned long long) nL);
		Die("BCAData::read_codeseq_nu()");
		}
	return L;
	}

flat_chain_t* BCAData::read_flat_chain(uint64 ChainIdx) const
	{
	asserta(m_Reading && !m_Writing);
	uint L = flat_chain_cap_L(GetSeqLength(ChainIdx));
	auto chain = flat_chain_t::newflat(L);
	uint64 SeqOffset = GetSeqOffset(ChainIdx);
	uint64 nL = ReadStdioFile64_NoFail(m_f, SeqOffset, chain->m_aa->m_data, L);
	if (nL != L)
		{
		Log("FN=%s\n", m_FN.c_str());
		Log("ChainIdx=%u\n", ChainIdx);
		Log("Chains=%u\n", SIZE(m_SeqLengths));
		Log("L=%u\n", L);
		Log("SeqOffset=%llu\n", (unsigned long long) SeqOffset);
		Log("nL=%llu\n", (unsigned long long) nL);
		Die("BCAData::ReadChain(#2)");
		}

	uint64 BytesToRead = 6*L;
	uint64 nIC = ReadStdioFile64_NoFail(m_f, SeqOffset + L,
		chain->m_xyz->m_data, BytesToRead);
	if (nIC != BytesToRead)
		{
		Log("FN=%s\n", m_FN.c_str());
		Log("ChainIdx=%u\n", ChainIdx);
		Log("Chains=%u\n", SIZE(m_SeqLengths));
		Log("L=%u\n", L);
		Log("SeqOffset=%llu\n", (unsigned long long) SeqOffset);
		Log("nIC=%llu\n", (unsigned long long) nIC);
		Die("BCAData::ReadChain(#2)");
		}

	asserta(ChainIdx < SIZE(m_Labels));
	chain->m_label = m_Labels[ChainIdx];
	return chain;
	}

//void BCAData::ReadChain(uint64 ChainIdx, PDBChain &Chain) const
//	{
//	asserta(m_Reading && !m_Writing);
//	asserta(!m_HasNuSequences);
//	Chain.Clear();
//	uint L = GetSeqLength(ChainIdx);
//	uint64 SeqOffset = GetSeqOffset(ChainIdx);
//	char *Seq = myalloc(char, L+1);
//	uint64 nL = ReadStdioFile64_NoFail(m_f, SeqOffset, Seq, L);
//	if (nL != L)
//		{
//		Log("FN=%s\n", m_FN.c_str());
//		Log("ChainIdx=%u\n", ChainIdx);
//		Log("Chains=%u\n", SIZE(m_SeqLengths));
//		Log("L=%u\n", L);
//		Log("SeqOffset=%llu\n", (unsigned long long) SeqOffset);
//		Log("nL=%llu\n", (unsigned long long) nL);
//		Die("BCAData::ReadChain(#2)");
//		}
//
//	Seq[L] = 0;
//	Chain.m_Seq = string(Seq);
//	myfree(Seq);
//
//	uint16_t *ICs = myalloc(uint16_t, 3*L);
//	uint64 BytesToRead = 6*L;
//	uint64 nIC = ReadStdioFile64_NoFail(m_f, SeqOffset + L, ICs, BytesToRead);
//	if (nIC != BytesToRead)
//		{
//		Log("FN=%s\n", m_FN.c_str());
//		Log("ChainIdx=%u\n", ChainIdx);
//		Log("Chains=%u\n", SIZE(m_SeqLengths));
//		Log("L=%u\n", L);
//		Log("SeqOffset=%llu\n", (unsigned long long) SeqOffset);
//		Log("nIC=%llu\n", (unsigned long long) nIC);
//		Die("BCAData::ReadChain(#2)");
//		}
//
//	Chain.CoordsFromICs(ICs, L);
//	myfree(ICs);
//	asserta(ChainIdx < SIZE(m_Labels));
//	Chain.m_Label = m_Labels[ChainIdx];
//	}

void BCAData::make_kappa_codeseqs(
	uint8_t ***ptr_kappa_codeseqs,
	uint **ptr_lengths) const
	{
	const uint nchain = GetChainCount();
	uint8_t **kappa_codeseqs = myalloc(uint8_t *, nchain);
	uint *lengths = myalloc(uint, nchain);

	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		const uint L = GetSeqLength(chainidx);
		lengths[chainidx] = L;
		uint8_t *kappa_codeseq = myalloc(uint8_t, L);
		read_codeseq_nu(kappa_codeseq, chainidx, L);
		chaq::codeseq_nu_to_kappa_inplace(kappa_codeseq, L);
		kappa_codeseqs[chainidx] = kappa_codeseq;
		}
	*ptr_kappa_codeseqs = kappa_codeseqs;
	*ptr_lengths = lengths;
	}

void BCAData::make_nu_and_kappa_codeseqs(
	uint8_t ***ptr_nu_codeseqs,
	uint8_t ***ptr_kappa_codeseqs,
	uint **ptr_lengths) const
	{
	const uint nchain = GetChainCount();
	uint8_t **nu_codeseqs = myalloc(uint8_t *, nchain);
	uint8_t **kappa_codeseqs = myalloc(uint8_t *, nchain);
	uint *lengths = myalloc(uint, nchain);

	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		const uint L = GetSeqLength(chainidx);
		lengths[chainidx] = L;
		uint8_t *nu_codeseq = myalloc(uint8_t, L);
		uint8_t *kappa_codeseq = myalloc(uint8_t, L);
		read_codeseq_nu(nu_codeseq, chainidx, L);
		chaq::codeseq_nu_to_kappa(nu_codeseq, L, kappa_codeseq, L);
		nu_codeseqs[chainidx] = kappa_codeseq;
		kappa_codeseqs[chainidx] = kappa_codeseq;
		}
	*ptr_nu_codeseqs = nu_codeseqs;
	*ptr_kappa_codeseqs = kappa_codeseqs;
	*ptr_lengths = lengths;
	}

void cmd_bca_stats()
	{
	BCAData BCA;
	BCA.Open(g_Arg1);
	uint ChainCount = BCA.GetChainCount();
	ProgressLog("%10u  Chains with_nu=%c (%s)\n",
		ChainCount,
		tof(BCA.m_HasNuSequences),
		FloatToStr(ChainCount));
	uint64 SumL = 0;
	for (uint i = 0; i < ChainCount; ++i)
		SumL += BCA.m_SeqLengths[i];
	ProgressLog("%10u  Residues (%s)\n",
		SumL,
		FloatToStr(double(SumL)));
	ProgressLog("%10u  Mean length\n", SumL/ChainCount);
	ProgressLog("%10.0f  Label data bytes (%s)\n",
		(double) BCA.m_LabelDataSize64,
		FloatToStr((double) BCA.m_LabelDataSize64));
	}
