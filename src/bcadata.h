#pragma once

class PDBChain;
#include "flat_chain.h"
#include "chaq.h"

// Binary C-alpha
class BCAData
	{
public:
	vector<string> m_Labels;
	vector<uint64_t> m_Offsets; // start of IC vector in file
	vector<uint32_t> m_SeqLengths;
	bool m_HasNuSequences = false;
	string m_FN;
	FILE *m_f = 0;
	bool m_Writing = false;
	bool m_Reading = false;
	uint64 m_SeqLengthsPos64 = UINT64_MAX;
	uint64 m_LabelDataSize64 = UINT64_MAX;
// Start of the contiguous nu-sequence section (BCB only; UINT64_MAX if none).
	uint64 m_NuSeqPos64 = UINT64_MAX;
// Prefix sums of m_SeqLengths: byte offset of chain idx's nu within the
// contiguous nu section (m_NuPrefix[idx]); built at Open.
	vector<uint64_t> m_NuPrefix;
	sid_t *m_distmx = 0;
	//chaq_vecs2 m_cv;
	//uint8_t *m_scratch_buffer = 0;
	//uint m_scratch_buffer_bytes = 0;
	uint8_t *m_codeseq_nu = 0;
// Temp file accumulating nu bytes during write; appended contiguously at close.
	FILE *m_nu_tmp_f = 0;
	string m_NuTmpFN;

// Sharded (parallel) writing. Each shard is owned by a single thread and
// streams its AA+ICs and nu bytes to private temp files plus local
// label/length tables. CloseSharded concatenates shards in shard order;
// because the reader recomputes all offsets from m_SeqLengths, the only
// requirement is that the four per-shard streams agree on chain order.
	struct Shard
		{
		FILE *aa_f = 0;			// AA + ICs bytes (7*L per chain)
		FILE *nu_f = 0;			// nu bytes (L per chain, BCB only)
		string aa_fn;
		string nu_fn;
		vector<string> labels;
		vector<uint32_t> seqlengths;
		sid_t *distmx = 0;		// per-shard scratch for nu computation
		uint8_t *codeseq_nu = 0;
		chaq_vecs2 cv;
		bool cv_inited = false;
		};
	vector<Shard *> m_Shards;

public:
	void Clear();
	void Create(const string &FN, bool WithNu = false);
	void Open(const string &FN);
	//void WriteChain(const PDBChain &Chain);

	const string &GetLabel(uint idx) const
		{
		asserta(idx < m_Labels.size());
		return m_Labels[idx];
		}

	uint GetSeqLength(uint idx) const
		{
		asserta(idx < m_SeqLengths.size());
		return m_SeqLengths[idx];
		}

	//void ReadChain(uint64 ChainIdx, PDBChain &Chain) const;
	flat_chain_t* read_flat_chain(uint64 ChainIdx) const;
	void write_flat_chain(const flat_chain_t *chain,
		chaq_vecs2 *cv);

// Parallel sharded writing (single-writer-per-shard, merged at close).
	void CreateSharded(const string &FN, bool WithNu, uint nshard);
	void write_flat_chain_shard(uint shard, const flat_chain_t *chain);
	void CloseSharded();

	void Close();
	uint GetChainCount() const { return SIZE(m_Labels); }
	uint64 GetSeqOffset(uint64 ChainIdx) const;
	uint GetSeqLength(uint64 ChainIdx) const;
	void append_codeseq_nu(const flat_chain_t *chain,
		chaq_vecs2 *cv);
	uint read_codeseq_nu(
		uint8_t *codeseq_nu, uint idx, uint buffer_length) const;
	uint64 get_offset_ICs(uint idx) const;
	uint64 get_offset_aaseq(uint idx) const;
	uint64 get_offset_nuseq(uint idx) const;
	void make_kappa_codeseqs(
		uint8_t ***ptr_kappa_codeseqs,
		uint **ptr_lengths) const;
	void make_nu_and_kappa_codeseqs(
		uint8_t ***ptr_nu_codeseqs,
		uint8_t ***ptr_kappa_codeseqs,
		uint **ptr_lengths) const;

private:
	void CloseWriter();
	void CloseReader();
// Shared by CloseWriter and CloseSharded.
	void AppendTempFileToMain(FILE *tmp_f);
	void WriteSeqLengthsLabelsHeader();
	};

const uint32_t BCA_MAGIC = 0xBCABCA;
// BCB layout v2: contiguous nu section + 4th header field (m_NuSeqPos64).
// Value bumped from 0xBCBBCB so pre-v2 .bcb files are rejected, not misread.
const uint32_t BCB_MAGIC = 0xBCBBC2;
