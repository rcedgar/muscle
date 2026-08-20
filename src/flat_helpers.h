#pragma once
#include "alpha.h"
#include "flat_dist_types.h"
#include "flat_distmx.h"
#include "fan.h"

class flat_aligner;

using colscorefn = float(uint i, uint j);
float sw_flat(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB,
	uint LA, uint LB, colscorefn sf,
	float Open, float Ext, uint &Loi, uint &Loj,
	char *path_buffer, uint &ncol);

float sw_flat_pssm(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profA, uint LA,
	const float *__restrict pssm, uint LB,
	const uint32_t * __restrict feature_block_offsets,
	uint nfeat,
	float Open, float Ext, uint &Loi, uint &Loj,
	char *path_buffer, uint &ncol);

float sw_flat_pssm_scoreonly(
	float *__restrict scratch_rows,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profQ, uint LQ,
	const float *__restrict pssmT, uint LT,
	const uint32_t * __restrict feature_block_offsets,
	uint nfeat, float Open, float Ext);

void log_profile(
	const string &label,
	const uint8_t *prof,
	uint nfeat,
	uint L);

void profiles2faprof(
	const string &fn,
	const vector<string> &feature_names,
	const vector<uint> &alpha_sizes,
	const vector<string> &labels,
	const vector<vector<uint8_t> > &profiles);

void read_profiles_faprof(
	const string &fn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles);

void read_logoddsvec(
	const vector<string> &fns,
	vector<vector<float> > &logoddsvec);

void read_logoddsvec_pattern(
	const string &fnpattern,
	const vector<string> &feature_names,
	const vector<uint> &alpha_sizes,
	vector<vector<float> > &logoddsvec);

void read_profiles_and_logoddsvec(
	const string &specfn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles,
	vector<vector<float> > &logoddsvec);

void check_profiles(
	vector<vector<uint8_t> > &profiles,
	vector<uint> &alpha_sizes);

uint32_t get_flat_pssm_feature_block_offsets(
	const uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	uint32_t * __restrict feature_block_offsets);

void fill_flat_pssm(
	const uint8_t * __restrict profQ,
	uint32_t LQ,
	uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	const uint32_t * __restrict feature_block_offsets,
	const float *const * __restrict weighted_logoddsvec,
	float * __restrict pssm);

void fill_flat_pssm_reversed(
	const uint8_t * __restrict profQ,
	uint32_t LQ,
	uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	const uint32_t * __restrict feature_block_offsets,
	const float *const * __restrict weighted_logoddsmxvec,
	float * __restrict pssm);

void fill_smx_using_flat_pssm(
	const uint8_t * __restrict profA,
	uint32_t LA,
	uint32_t LB,
	uint32_t nfeat,
	const uint32_t * __restrict feature_block_offsets,
	const float * __restrict pssm,
	float * __restrict smx);

void write_flat_aln(
	FILE *f,
	const string &labelQ, const uint8_t *profQ, uint LQ,
	const string &labelT, const uint8_t *profT, uint LT,
	uint LoQ, uint LoT, const string &path,
	const vector<string> &feature_names,
	const vector<uint> &alpha_sizes,
	const vector<string> &symbolsvec,
	float score,
	const string &style = "");

void WriteLocalAln(FILE *f, const string &LabelA, const uint8_t *A,
	const string &LabelB, const uint8_t *B,
	uint Loi, uint Loj, const char *Path);

void flat_logodds_symbols(
	const float *logodds,
	uint alpha_size,
	string &symbols);

void read_fasta_label2idx(
	const string &fafn,
	unordered_map<string, uint> &label2idx);

void read_feature_fasta(
	const string &fafn,
	uint alpha_size,
	const unordered_map<string, uint> &label2idx,
	vector<vector<uint8_t> > &codeseqs);

void make_fn_pattern(
	const string &fnpattern,
	const string &feature_name,
	string &fn);

uint read_logodds(
	const string &fn,
	vector<float> &logoddsmx);

void write_logoddsmx(FILE *f,
	const vector<vector<double> > &logoddsmx,
	bool asintegers);

void write_flat_logoddsmx(FILE *f,
	const vector<float> &logoddsmx,
	uint alpha_size,
	bool asintegers);

uint read_logodds_and_freqmx(
	const string &fn,
	vector<double> &logoddsmx,
	vector<double> &freqmx);

double get_expected_score(
	const vector<double> &freqs,
	vector<vector<double> > &scoremx);

double get_expected_score_flat(
	const vector<double> &freqs,
	const vector<double> &scoremx);

double get_relative_entropy_flat(
	const vector<double> &freqmx,
	const vector<double> &scoremx,
	uint alpha_size);

void get_logoddsmx_from_flat_freqmx(
	const vector<double> &freqmx,
	uint alpha_size,
	vector<double> &logoddsmx);

void trunc_label(const string &Label,
	string &TruncatedLabel);

void trunc_label(string &Label);

void ChainizeLabel(string &Label, const string &ChainStr);
char ExtractChainIdFromLabel(const string &Label);

uint32_t get_alpha_size_from_feature_name(const string &name);

static inline const uint8_t *get_letter2char(uint alpha_size)
	{
	return (alpha_size == 20 ? g_LetterToCharAmino : g_LetterToCharMu);
	}

static inline const uint8_t *get_char2letter(uint alpha_size)
	{
	return (alpha_size == 20 ? g_CharToLetterAmino : g_CharToLetterMu);
	}

static inline const uint8_t *get_char2letter(const string &feature_name)
	{
	uint alpha_size = get_alpha_size_from_feature_name(feature_name);
	return (alpha_size == 20 ? g_CharToLetterAmino : g_CharToLetterMu);
	}

float flat_get_dali(
	const string &labelQ, const string &labelT,
	const string &path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const sid_t *distmxQ,
	const sid_t *distmxT);

float flat_get_dalix(
	const char *path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float *colscores);

float flat_get_dalix3(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float *colscores);

float flat_get_dali3(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT);

float flat_getlddt_muscle_some_floats(
	const uint32_t *posQs,
	const uint32_t LQ,
	const uint32_t *posTs,
	const uint32_t LT,
	const uint ncol,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint32_t *nr_considered_vec,
	uint32_t *nr_preserved_vec);

float flat_getlddt_muscle_some_floats3(
	const string &labelQ, const string &labelT,
	const string &path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const sid_t *distmxQ,
	const sid_t *distmxT);

float flat_getlddt_old(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT);

float flat_getlddt_muscle_some_floats4(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT);

float flat_get_entropy2(
	const flat_aligner &fa,
	const uint8_t *profQ,
	const uint8_t *profT,
	uint nfeat, uint fi);

void flat_reverse_profile(
	const uint8_t *prof,
	uint32_t L,
	uint32_t nfeat,
	uint8_t *revprof);

void flat_reverse_distmx(
	cp_sid_t distmx, uint32_t L, p_sid_t reversed_distmx);

void trunc_label(const string &label, string &tlabel);

float flat_getlddt_muscle_some_floats2(
	const string &labelQ, const string &labelT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint LQ, uint LT,
	const vector<uint32_t> &posQs,
	const vector<uint32_t> &posTs);

float flat_getlddt_old_some_floats(
	const string &labelQ, const string &labelT,
	const uint32_t *posQs,
	const uint32_t LQ,
	const uint32_t *posTs,
	const uint32_t LT,
	const uint ncol,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint32_t *nr_considered_vec,
	uint32_t *nr_preserved_vec);

uint16_t *read_quantize(
	const string &fn,
	uint alpha_size,
	uint16_t &median);

void get_alpha_names_from_peaker_spec_file_lines(
	vector<string> &lines,
	vector<string> &alpha_names);

void GetFeatures(
	const string &varstr,
	vector<string> &feature_names,
	vector<float> &weights);

void path2posvecs(
	const string &labelQ, const string &labelT,
	const string &path,
	uint loQ, uint LQ,
	uint loT, uint LT,
	vector<uint> &posQs,
	vector<uint> &posTs);

float flat_get_entropy(
	const string &labelQ, const string &labelT,
	const string &path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const uint8_t *profQ,
	const uint8_t *profT,
	uint nfeat, uint fi);

void parse_varstr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

void flat_classify_params(
	const vector<string> &names,
	const vector<float> &values,
	vector<string> &alphan_ames,
	vector<float> &weights,
	vector<string> &scalar_names,
	vector<float> &scalar_values);

FAN parse_alpha_name(
	const string &alpha_name,
	uint &alpha_size);

void flat_make_varstr(
	const flat_params &params,
	string &varstr);

void flat_make_peaker_spec_const(
	const flat_params &params,
	vector<string> &lines);

void flat_make_peaker_spec_range(
	const flat_params &params,
	vector<string> &lines,
	double rate = 1.1);

void validate_flat_peaker_spec(
	const vector<string> &SpecLines,
	vector<string> &alpha_names);

uint path2posvecs3(
	const char *path, uint ncol,
	uint loQ, uint LQ,
	uint loT, uint LT,
	uint *posQs,
	uint *posTs,
	uint bufnpos);

float flat_get_dali4(
	const uint *posQs, uint LQ, 
	const uint *posTs, uint LT, uint nmatch,
	const sid_t *distmxQ, const sid_t *distmxT);

/* Bytes for float[nmatch*nmatch] pair-term cache (includes float alignment padding). */
size_t dali_greedy_terms_bytes(uint nmatch);

/* Bytes for float[nmatch] column scores + uint8_t[nmatch] active flags
   (includes float alignment padding). */
size_t dali_greedy_work_bytes(uint nmatch);

/* Adaptive greedy DALI column pruning. No allocation.
   Fills retained_cols[0..nretained) with original column indices in ascending order.
   Caller must supply terms_scratch and work_scratch sized via the helpers above;
   retained_cols must have capacity >= nmatch (except nmatch==0). */
float dali_greedy(
	const uint *posQs, uint LQ,
	const uint *posTs, uint LT, uint nmatch,
	const sid_t *distmxQ, const sid_t *distmxT,
	void *terms_scratch, size_t terms_bytes,
	void *work_scratch, size_t work_bytes,
	uint *retained_cols, uint retained_cols_capacity,
	uint &nretained);

void codeseq_to_hexfasta(FILE *f, const string &label,
	const uint8_t *codeseq, uint L);

void codeseq_to_fasta(FILE *f, const string &label,
	const uint8_t *codeseq, uint L, uint alpha_size);

uint logodds_and_freqmx_from_lines(
	const vector<string> &lines,
	vector<double> &logoddsmx,
	vector<double> &freqmx);

void set_sec4_groups(const string &sec4_groups);
void decide_query_or_db_kmer_neighborhood(uint QSeqCount, uint DBSeqCount);
class kappa_dex;
class kappa_mermx;
void setup_kappa_qkmer_index(kappa_dex &QKmerIndex, const kappa_mermx &ScoreMx);
static uint16_t CoordToIC(float X) { return uint16_t((X + 1000)*10 + 0.5); }
static float ICToCoord(uint16_t IC) { return float(IC/10.0f) - 1000; }

static float GetFloatFromString(const string &s, uint Pos, uint n)
	{
	string t = s.substr(Pos, n);
	StripWhiteSpace(t);
	float Value = (float) StrToFloat(t);
	return Value;
	}

static void GetXYZFromATOMLine(const string &InputLine,
  float &x, float &y, float &z)
	{
	x = GetFloatFromString(InputLine, 30, 8);
	y = GetFloatFromString(InputLine, 38, 8);
	z = GetFloatFromString(InputLine, 46, 8);
	}

static void SetXYZInATOMLine(const string &InputLine,
  float x, float y, float z, string &OutputLine)
	{
	string sx;
	string sy;
	string sz;
	Ps(sx, "%8.3f", x);
	Ps(sy, "%8.3f", y);
	Ps(sz, "%8.3f", z);
	asserta(SIZE(sx) == 8);
	asserta(SIZE(sy) == 8);
	asserta(SIZE(sz) == 8);

	OutputLine = InputLine;
	for (uint i = 0; i < 8; ++i)
		{
		OutputLine[30+i] = sx[i];
		OutputLine[38+i] = sy[i];
		OutputLine[46+i] = sz[i];
		}
	}

static bool IsATOMLine(const string &Line)
	{
	if (SIZE(Line) < 27)
		return false;
// Insertion code is old PDB hack to preserve numbering from
// a reference sequence like NAST or ICTV numbers -- should
// just ignore residue numbering.
	//if (Line[26] != ' ') // insertion code
	//	return false;
	if (strncmp(Line.c_str(), "ATOM  ", 6) == 0)
		return true;
	return false;
	}
