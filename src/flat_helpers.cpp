#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "tabbedlines.h"
#include "flat_helpers.h"
#include "features.h"
#include "triangle.h"
#include "kappa_dex.h"
#include "kappa_mermx.h"
#include "kappa_filter_params.h"

void flat_reverse_profile(
	const uint8_t *prof,
	uint32_t L,
	uint32_t nfeat,
	uint8_t *revprof)
	{
	for (uint32_t fi = 0; fi < nfeat; ++fi)
		{
		const uint8_t *row = prof + fi*L;
		uint8_t *revrow = revprof + fi*L;
		for (uint32_t pos = 0; pos < L; ++pos)
			{
			assert(fi*L + pos < L*nfeat);
			assert(fi*L + (L-1-pos) < L*nfeat);
			uint8_t code = row[L - 1 - pos];
			revrow[pos] = code;
			}
		}
	}

void flat_reverse_distmx(
	cp_sid_t distmx,
	uint32_t L,
	uint32_t M,
	p_sid_t reversed_distmx)
	{
	for (uint32_t i = 0; i < L; ++i)
		{
		const uint32_t jend = min(i + M, L - 1);
		for (uint32_t j = i + 1; j <= jend; ++j)
			{
			uint32_t k_src = banded_ij_to_k(i, j);
			uint32_t ir = L - 1 - j;
			uint32_t jr = L - 1 - i;
			uint32_t k_dst = banded_ij_to_k(ir, jr);
			reversed_distmx[k_dst] = distmx[k_src];
			}
		}
	}

void read_fasta_label2idx(
	const string &fafn,
	unordered_map<string, uint> &label2idx)
	{
	label2idx.clear();
	SeqDB db_fa;
	db_fa.FromFasta(fafn, false);
	const uint nseqs = db_fa.GetSeqCount();
	for (uint i = 0; i < nseqs; ++i)
		{
		const string &label = db_fa.GetLabel(i);
		label2idx[label] = i;
		}
	}

void read_feature_fasta(
	const string &fafn,
	uint alpha_size,
	const unordered_map<string, uint> &label2idx,
	vector<vector<uint8_t> > &codeseqs)
	{
	codeseqs.clear();

	const uint8_t *char2letter = get_char2letter(alpha_size);

	SeqDB db_fa;
	db_fa.FromFasta(fafn, false);
	const uint db_seq_count = db_fa.GetSeqCount();
	for (uint i = 0; i < db_seq_count; ++i)
		trunc_label(db_fa.m_Labels[i]);
	db_fa.SetLabelToIndex();
	const uint nseqs = SIZE(label2idx);
	codeseqs.resize(nseqs);
	for (unordered_map<string, uint>::const_iterator iter =
		label2idx.begin(); iter != label2idx.end(); ++iter)
		{
		string label = iter->first;
		trunc_label(label);
		uint idx = iter->second;
		uint seqidx = db_fa.GetSeqIndex(label);
		uint L = db_fa.GetSeqLength(seqidx);
		asserta(idx < nseqs);
		asserta(codeseqs[idx].size() == 0);
		codeseqs[idx].resize(L);
		const byte *byteseq = db_fa.GetByteSeq(seqidx);
		for (uint pos = 0; pos < L; ++pos)
			{
			uint8_t code = char2letter[byteseq[pos]];
			if (code < alpha_size)
				codeseqs[idx][pos] = code;
			else if (code == 0xff)
				codeseqs[idx][pos] = 0;
			else
				asserta(false);
			}
		}
	}

void read_profiles_from_fastas(
	const vector<string> &fafns,
	const vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles)
	{
	profiles.clear();
	labels.clear();

	const uint nfeat = SIZE(fafns);
	asserta(SIZE(alpha_sizes) == nfeat);
	asserta(nfeat > 0);

	unordered_map<string, uint> label2idx;
	read_fasta_label2idx(fafns[0], label2idx);
	const uint nprof = SIZE(label2idx);

	vector<vector<vector<uint8_t> > > codeseqsvec(nfeat);

	read_feature_fasta(fafns[0], alpha_sizes[0], label2idx, codeseqsvec[0]);
	for (uint fi = 1; fi < nfeat; ++fi)
		read_feature_fasta(fafns[fi], alpha_sizes[fi], label2idx, codeseqsvec[fi]);

	profiles.resize(nprof);
	for (auto iter : label2idx)
		{
		void TruncLabel(string &lab);
		string label = iter.first;
		uint idx = iter.second;
		labels.push_back(label);
		vector<uint8_t> &profile = profiles[idx];
		uint L = SIZE(codeseqsvec[0][idx]);
		profile.resize(nfeat*L, 0xff);
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			const vector<uint8_t> &codeseq = codeseqsvec[fi][idx];
			for (uint i = 0; i < L; ++i)
				asserta(codeseq[i] != 0xff);
			for (uint pos = 0; pos < L; ++pos)
				profile[fi*L + pos] = codeseq[pos];
			}
		for (uint i = 0; i < nfeat*L; ++i)
			asserta(profile[i] != 0xff);
		}
	}

uint lines2logoddsmx(
	const vector<string> &lines,
	vector<float> &logoddsmx)
	{
	logoddsmx.clear();
	tabbedlines tl(lines);
	uint alpha_size = tl.get_int("logodds");
	asserta(alpha_size != 0);
	logoddsmx.resize(alpha_size*alpha_size);
	tl.get_float_flat_square_mx(alpha_size, logoddsmx.data());
	return alpha_size;
	}

uint read_logodds(
	const string &fn,
	vector<float> &logoddsmx)
	{
	vector<string> lines;
	ReadLinesFromFile(fn, lines);
	return lines2logoddsmx(lines, logoddsmx);
	}

uint logodds_and_freqmx_from_lines(
	const vector<string> &lines,
	vector<double> &logoddsmx,
	vector<double> &freqmx)
	{
	logoddsmx.clear();
	tabbedlines tl(lines);
	uint alpha_size = tl.get_int("logodds");
	asserta(alpha_size != 0);
	logoddsmx.resize(alpha_size*alpha_size);
	freqmx.resize(alpha_size*alpha_size);
	tl.get_double_flat_square_mx(alpha_size, logoddsmx.data());

	uint alpha_size2 = tl.get_int("freqs");
	asserta(alpha_size2 == alpha_size);
	tl.get_double_flat_square_mx(alpha_size, freqmx.data());
	return alpha_size;
	}

uint read_logodds_and_freqmx(
	const string &fn,
	vector<double> &logoddsmx,
	vector<double> &freqmx)
	{
	vector<string> lines;
	ReadLinesFromFile(fn, lines);
	uint alpha_size = logodds_and_freqmx_from_lines(
		lines, logoddsmx, freqmx);
	return alpha_size;
	}

void read_logoddsvec(
	const vector<string> &fns,
	vector<vector<float> > &logoddsvec)
	{
	logoddsvec.clear();
	for (auto fn : fns)
		{
		vector<float> logoddsmx;
		read_logodds(fn, logoddsmx);
		logoddsvec.push_back(logoddsmx);
		}
	}

uint32_t get_alpha_size_from_feature_name(const string &name)
	{
	if (name == "aa" || name == "AA")
		return 20;
	uint n = 0;
	for (auto c : name)
		{
		if (isdigit(c))
			n = 10*n + c - '0';
		else
			n = 0;
		}
	if (n == 0)
		Die("get_alpha_size_from_feature_name(%s)", name.c_str());
	return n;
	}

// @=name endswith AS e.g. AA20
void make_fn_pattern(
	const string &fnpattern,
	const string &feature_name,
	string &fn)
	{
	fn.clear();
	for (auto c : fnpattern)
		{
		if (c == '@')
			fn += feature_name;
		else
			fn += c;
		}
	}

void read_logoddsvec_pattern(
	const string &fnpattern,
	const vector<string> &feature_names,
	const vector<uint> &alpha_sizes,
	vector<vector<float> > &logoddsvec)
	{
	const uint nfeat = SIZE(feature_names);
	asserta(SIZE(alpha_sizes) == nfeat);
	vector<string> fns(nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		make_fn_pattern(
			fnpattern,
			feature_names[fi],
			fns[fi]);
	read_logoddsvec(fns, logoddsvec);
	}

void log_flat_square_mx(const float *mx, uint n)
	{
	for (uint i = 0; i < n; ++i)
		{
		Log("%2u  |", i);
		for (uint j = 0; j < n; ++j)
			Log(" %7.3g", mx[n*i + j]);
		Log("\n");
		}
	}

void check_profile(
	vector<uint8_t> &profile,
	vector<uint> &alpha_sizes)
	{
	const uint n = SIZE(profile);
	const uint nfeat = SIZE(alpha_sizes);
	asserta(n%nfeat == 0);
	const uint L = n/nfeat;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = alpha_sizes[fi];
		for (uint i = 0; i < L; ++i)
			asserta(profile[fi*L + i] < AS);
		}
	}

void check_profiles(
	vector<vector<uint8_t> > &profiles,
	vector<uint> &alpha_sizes)
	{
	const uint nprof = SIZE(profiles);
	for (uint i = 0; i < nprof; ++i)
		check_profile(profiles[i], alpha_sizes);
	}

void profiles2faprof(
	const string &fn,
	const vector<string> &feature_names,
	const vector<uint> &alpha_sizes,
	const vector<string> &labels,
	const vector<vector<uint8_t> > &profiles)
	{
	FILE *fap = CreateStdioFile(fn);
		const uint nprof = SIZE(labels);
	asserta(SIZE(profiles) == nprof);
	const uint nfeat = SIZE(feature_names);
	asserta(SIZE(alpha_sizes) == nfeat);
	for (uint seqidx = 0; seqidx < nprof; ++seqidx)
		{
		const string &label = labels[seqidx];
		const vector<uint8_t> &profile = profiles[seqidx];
		asserta(SIZE(profile)%nfeat == 0);
		const uint L = SIZE(profile)/nfeat;
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			uint alpha_size = alpha_sizes[fi];
			const uint8_t *letter2char = get_letter2char(alpha_size);

			string seq;
			for (uint pos = 0; pos < L; ++pos)
				{
				uint8_t code = profile[fi*L + pos];
				seq += letter2char[code];
				}
			string label_feat;
			Psa(label_feat, "%s:%s",
				label.c_str(),
				feature_names[fi].c_str());
			SeqToFasta(fap, label_feat, seq, L);
			}
		}
	CloseStdioFile(fap);
	}

void read_profiles_faprof(
	const string &fn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles)
	{
	if (fn == "")
		return;

	feature_names.clear();
	alpha_sizes.clear();
	labels.clear();
	profiles.clear();

	SeqDB DB;
	DB.FromFasta(fn);
	const uint ndbseq = DB.GetSeqCount();
	asserta(ndbseq > 0);
	vector<string> flds;
	string acc;
	uint nfeat = 0;
	for (uint i = 0; i < ndbseq; ++i)
		{
		const string &label = DB.GetLabel(i);
		Split(label, flds, ':');
		asserta(SIZE(flds) == 2);
		if (i == 0)
			{
			acc = flds[0];
			continue;
			}
		else
			{
			if (flds[0] != acc)
				{
				nfeat = i;
				break;
				}
			}
		}
	asserta(nfeat > 0);
	for (uint i = 0; i < nfeat; ++i)
		{
		const string &label = DB.GetLabel(i);
		Split(label, flds, ':');
		asserta(flds.size() == 2);
		string alpha_annot = flds[1];
		Split(alpha_annot, flds, '*');
		asserta(flds.size() == 2);
		feature_names.push_back(flds[0]);
		alpha_sizes.push_back(StrToUint(flds[1]));
		}
	ProgressLog("%u features", nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		ProgressLog(" %s*%u",
			feature_names[fi].c_str(),
			alpha_sizes[fi]);
	ProgressLog("\n");
	asserta(ndbseq%nfeat == 0);
	uint nprof = ndbseq/nfeat;
	profiles.resize(nprof);
	for (uint profidx = 0; profidx < nprof; ++profidx)
		{
		uint L = DB.GetSeqLength(nfeat*profidx);
		uint profile_length = nfeat*L;
		vector<uint8_t> &profile = profiles[profidx];
		profile.resize(profile_length);
		uint k = 0;
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			if (fi == 0)
				{
				const string &label = DB.GetLabel(nfeat*profidx);
				Split(label, flds, ':');
				asserta(flds.size() == 2);
				const string &acc = flds[0];
				labels.push_back(acc);
				}
			uint alpha_size = alpha_sizes[fi];
			const uint8_t *char2letter = get_char2letter(alpha_size);
			uint L_fi = DB.GetSeqLength(nfeat*profidx + fi);
			asserta(L_fi == L);
			const string &seq = DB.GetSeq(nfeat*profidx + fi);
			asserta(SIZE(seq) == L);
			for (uint k = 0; k < L; ++k)
				profile[fi*L + k] = char2letter[seq[k]];
			}
		}
	}

void read_profiles_and_logoddsvec(
	const string &specfn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles,
	vector<vector<float> > &logoddsvec)
	{
	vector<string> lines;
	ReadLinesFromFile(specfn, lines);
	uint nfeat = SIZE(lines);

	vector<string> fafns;
	vector<string> logoddsfns;
	vector<float> weights;
	float sumw = 0;
	for (auto line : lines)
		{
		vector<string> flds;
		Split(line, flds, '\t');
		asserta(SIZE(flds) == 5);
		feature_names.push_back(flds[0]);
		fafns.push_back(flds[1]);
		alpha_sizes.push_back(StrToUint(flds[2]));
		logoddsfns.push_back(flds[3]);
		float w = (float) StrToFloat(flds[4]);
		weights.push_back(w);
		sumw += w;
		}
	asserta(SIZE(weights) == nfeat);

	asserta(sumw > 0);
	float sumw2 = 0;
	for (uint i = 0; i < nfeat; ++i)
		{
		float w = weights[i]/sumw;
		sumw2 += w;
		weights[i] = w;
		}
	asserta(sumw2 > 0.99 && sumw2 < 1.01);

	read_logoddsvec(logoddsfns, logoddsvec);
	asserta(SIZE(logoddsvec) == nfeat);

	read_profiles_from_fastas(fafns, alpha_sizes, labels, profiles);
	}

void log_profile(
	const string &label,
	const uint8_t *prof,
	uint nfeat,
	uint L)
	{
	Log("log_profile(%s) L=%u nfeat=%u\n", label.c_str(), L, nfeat);
	Log("  pos  ");
	for (uint fi = 0; fi < nfeat; ++fi)
		Log(" %2u", fi);
	Log("\n");
	for (uint pos = 0; pos < L; ++pos)
		{
		Log("[%4u] ", pos);
		for (uint fi = 0; fi < nfeat; ++fi)
			Log(" %2x", prof[fi*L + pos]);
		Log("\n");
		}
	}

#if 0
void cmd_flat_profiles()
	{
	const string &specfn = g_Arg1;
	vector<string> labels;
	vector<vector<uint8_t> > profiles;
	vector<string> feature_names;
	vector<uint> alpha_sizes;
	vector<vector<float> > logoddsvec;
	read_profiles_and_logoddsvec(
		specfn,
		feature_names,
		alpha_sizes,
		labels,
		profiles,
		logoddsvec);

	const uint nfeat = SIZE(feature_names);
	const uint nprof = SIZE(labels);

	asserta(SIZE(alpha_sizes) == nfeat);
	asserta(SIZE(profiles) == nprof);

	check_profiles(profiles, alpha_sizes);

	for (uint fi = 0; fi < nfeat; ++fi)
		{
		Log("\n%s\n", feature_names[fi].c_str());
		log_flat_square_mx(logoddsvec[fi].data(), alpha_sizes[fi]);
		}

	// Convert profiles back to FASTA for correctness checking
	if (optset_output2)
		{
		const uint nprof = SIZE(labels);
		asserta(SIZE(profiles) == nprof);
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			uint alpha_size = alpha_sizes[fi];
			const uint8_t *letter2char = (alpha_size == 20 ? g_LetterToCharAmino : g_LetterToCharMu);

			string fn = opt(output2) + feature_names[fi];
			Progress("FASTA %s\n", fn.c_str());
			FILE *ffa = CreateStdioFile(fn);
			for (uint seqidx = 0; seqidx < nprof; ++seqidx)
				{
				const string &label = labels[seqidx];
				const vector<uint8_t> &profile = profiles[seqidx];
				asserta(SIZE(profile)%nfeat == 0);
				const uint L = SIZE(profile)/nfeat;
				string seq;
				for (uint pos = 0; pos < L; ++pos)
					{
					uint8_t code = profile[fi*L + pos];
					seq += letter2char[code];
					}
				SeqToFasta(ffa, label, seq);
				}
			CloseStdioFile(ffa);
			}
		}

	// Convert profiles to faprof
	if (optset_faprof)
		profiles2faprof(opt(faprof), feature_names, alpha_sizes, labels, profiles);
	}
#endif

void write_flat_aln(
	FILE *f,
	const string &labelQ, const uint8_t *profQ, uint LQ,
	const string &labelT, const uint8_t *profT, uint LT,
	uint LoQ, uint LoT, const string &path,
	const vector<string> &feature_names,
	const vector<uint> &alpha_sizes,
	const vector<string> &symbolsvec,
	float score,
	const string &style)
	{
	if (f == 0)
		return;
	fprintf(f, "\n");
	uint nfeat = SIZE(feature_names);
	asserta(alpha_sizes.size() == nfeat);

	size_t ncol = path.size();
	vector<string> feature_rowsQ(nfeat);
	vector<string> feature_rowsT(nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		string feature_rowQ = feature_rowsQ[fi];
		string feature_rowT = feature_rowsT[fi];
		string annot_row;

		feature_rowQ.reserve(ncol);
		feature_rowT.reserve(ncol);
		annot_row.reserve(ncol);
		uint alpha_size = alpha_sizes[fi];
		uint posQ = 0;
		uint posT = 0;
		const uint8_t *letter2char = get_letter2char(alpha_size);
		for (uint col = 0; col < ncol; ++col)
			{
			uint8_t codeQ = profQ[fi*LQ + posQ];
			uint8_t codeT = profT[fi*LT + posT];
			char c = path[col];
			if (c == 'M')
				annot_row += (codeQ == codeT) ? '|' :
					symbolsvec[fi][codeQ*alpha_size + codeT];
			else
				annot_row += ' ';

			if (c == 'M' || c == 'D')
				{
				feature_rowQ += letter2char[codeQ];
				++posQ;
				}
			else
				feature_rowQ += '-';

			if (c == 'M' || c == 'I')
				{
				feature_rowT += letter2char[codeT];
				++posT;
				}
			else
				feature_rowT += '-';
			}
		fprintf(f, "\n");
		fprintf(f, "%s", feature_rowQ.c_str());
		fprintf(f, "  %8.8s*%2u", feature_names[fi].c_str(), alpha_sizes[fi]);
		fprintf(f, "  %s\n", labelQ.c_str());

		fprintf(f, "%s\n", annot_row.c_str());

		fprintf(f, "%s", feature_rowT.c_str());
		fprintf(f, "  %8.8s*%2u", feature_names[fi].c_str(), alpha_sizes[fi]);
		fprintf(f, "  %s\n", labelT.c_str());
		}
	fprintf(f, "score %.1f\n", score);
	}

#if 0
void cmd_test_faprof()
	{
	vector<string> feature_names;
	vector<uint> alpha_sizes;
	vector<string> labels;
	vector<vector<uint8_t> > profiles;
	read_profiles_faprof(g_Arg1,
		feature_names,
		alpha_sizes,
		labels,
		profiles);
	profiles2faprof(opt(output), feature_names, alpha_sizes, labels, profiles);
	}
#endif

void flat_logodds_symbols(
	const float *logodds,
	uint alpha_size,
	string &symbols)
	{
	symbols.clear();
	float min_score = FLT_MAX;
	float max_score = FLT_MAX;
	for (uint i = 0; i < alpha_size*alpha_size; ++i)
		{
		float score = logodds[i];
		min_score = (i == 0 ? score : min(min_score, score));
		max_score = (i == 0 ? score : max(max_score, score));
		}

	// __. +*^
	// 0123456
	static const char s[7] = { 'V', '_', '.', ' ', '+', '*', '^' };
	for (uint i = 0; i < alpha_size*alpha_size; ++i)
		{
		float score = logodds[i];
		uint k = uint(7*(score - min_score)/(max_score - min_score + max_score/7));
		symbols += s[k];
		}
	}

void cmd_flat_logodds_info()
	{
	vector<float> logodds;
	uint AS = read_logodds(g_Arg1, logodds);
	asserta(logodds.size() == AS*AS);

	float sum_diag = 0;
	float min_diag = 999;
	float max_diag = -999;
	for (uint i = 0; i < AS; ++i)
		{
		float score = logodds[AS*i + i];
		sum_diag += score;
		max_diag = max(score, max_diag);
		min_diag = min(score, min_diag);
		}
	float sum = 0;
	float min_offdiag_score = 999;
	float max_offdiag_score = -999;
	for (uint i = 0; i < AS; ++i)
		{
		for (uint j = 0; j < AS; ++j)
			{
			float score = logodds[i*AS + j];
			if (i != j)
				{
				min_offdiag_score = min(score, min_offdiag_score);
				max_offdiag_score = max(score, max_offdiag_score);
				}
			sum += score;
			}
		}

	float mean = sum/(AS*AS);
	float mean_off_diag = (sum - sum_diag)/(AS*(AS-1));

	string symbols;
	flat_logodds_symbols(logodds.data(), AS, symbols);

	ProgressLog("%s\n", g_Arg1.c_str());
	ProgressLog("%8.3f  min off-diag\n", min_offdiag_score);
	ProgressLog("%8.3f  min diag\n", min_diag);
	ProgressLog("%8.3f  max off-diag\n", max_offdiag_score);
	ProgressLog("%8.3f  max diag\n", max_diag);
	ProgressLog("%8.3f  mean\n", mean);
	ProgressLog("%8.3f  mean_off_diag\n", mean_off_diag);

	for (uint i = 0; i < AS; ++i)
		{
		for (uint j = 0; j < AS; ++j)
			Log("%c", symbols[i*AS + j]);
		Log("\n");
		}
	}

static bool LabelAlreadyHasChain(const string &Label, 
  const string &ChainStr)
	{
	uint chn = SIZE(ChainStr);
	if (chn != 1)
		return false;
	uint labn = SIZE(Label);
	if (labn < 6)
		return false;
	if (tolower(Label[labn-1]) != tolower(ChainStr[chn-1]))
		return false;
	char c = Label[labn-2];
	if (c == '_' || c == ':' || c == '.')
		return true;
	return false;
	}

void ChainizeLabel(string &Label, const string &_ChainStr)
	{
	if (opt(nochainchar))
		return;
	string ChainStr = _ChainStr;
	if (ChainStr == "" || ChainStr == " ")
		ChainStr = '_';
	if (LabelAlreadyHasChain(Label, _ChainStr))
		return;
	Label += (optset_chainsep ? string(opt(chainsep)) : "_");
	Label += ChainStr;
	}

char ExtractChainIdFromLabel(const string &Label)
	{
	string lab = Label;
	trunc_label(lab);
	if (opt(nochainchar))
		return 'A';
	const string sep = optset_chainsep ? opt(chainsep) : "_";
	const uint sepLen = SIZE(sep);
	const uint labLen = SIZE(lab);
	if (sepLen > 0 && labLen >= sepLen + 1)
		{
		if (lab.substr(labLen - sepLen - 1, sepLen) == sep)
			return lab[labLen - 1];
		}
	if (labLen >= 2)
		{
		char sep2 = lab[labLen - 2];
		if (sep2 == '_' || sep2 == ':' || sep2 == '.')
			return lab[labLen - 1];
		}
	return 'A';
	}

void GetFallbackLabelFromFN(const string &FN, string &Label)
	{
	GetStemName(FN, Label);

/***
Special case for anomalous SCOP40 domain names

	d1dnu.1.pdb
	01234567890

d1dnu.1 d3n55.1 d1o7d.2 d1r8o.1 d1dy9.1 d1ko6.1 d1f8v.1 d1o7d.3 d1pyo.1
d1qtn.1 d1sc3.1 d1f2t.1 d1xew.1 d1q7l.1 d1wht.1 d1w2w.1 d1k3b.1 d1or0.1
d1gk9.1 d1k2x.1 d1apy.1 d2dg5.1 d1pya.1 d1mtp.1
***/

	string Ext;
	GetExtFromPathName(FN, Ext);
	ToLower(Ext);

// Special-case for downloaded PDB files e.g. pdb1iv1.ent
	if (Ext == "pdb" || Ext == "ent" || Ext == "pdb.gz" || Ext == "ent.gz")
		{
		if (Label.size() == 7 && Label[0] == 'p' && Label[1] == 'd' && Label[2] == 'b')
			{
			Label = Label.substr(3, string::npos);
			ToUpper(Label);
			}
		}
	}

void GetPathCounts(const char *Path, uint ncol, uint &M, uint &D, uint &I)
	{
	M = 0;
	D = 0;
	I = 0;
	for (uint i = 0; i < ncol; ++i)
		{
		char c = Path[i];
		if (c == 'M')
			++M;
		else if (c == 'D')
			++D;
		else if (c == 'I')
			++I;
		}
	}


void GetPathCounts(const string &Path, uint &M, uint &D, uint &I)
	{
	GetPathCounts(Path.c_str(), uint(Path.size()), M, D, I);
	}

void decide_query_or_db_kmer_neighborhood(uint QSeqCount, uint DBSeqCount)
	{
	extern bool g_QueryNeighborhood;
	g_QueryNeighborhood = true;
	if (opt(idxq))
		g_QueryNeighborhood = true;
	else if (opt(idxt))
		g_QueryNeighborhood = false;
	ProgressLog("Index %s k-mer neighborhoods\n", 
		g_QueryNeighborhood ? "query" : "target");
	return;

	//static const uint MAX_QUERY_CHAINS_FOR_QUERY_NEIGHBORHOOD = 100;
	//if (opt(idxq))
	//	g_QueryNeighborhood = true;
	//else if (opt(idxt))
	//	g_QueryNeighborhood = false;
	//else
	//	{
	//	if (QSeqCount <= MAX_QUERY_CHAINS_FOR_QUERY_NEIGHBORHOOD)
	//		g_QueryNeighborhood = true;
	//	else
	//		g_QueryNeighborhood = false;
	//	}
	//Log("g_QueryNeighborhood=%c\n", tof(g_QueryNeighborhood));
	}

void setup_kappa_qkmer_index(kappa_dex &QKmerIndex, const kappa_mermx &ScoreMx)
	{
	QKmerIndex.m_AddNeighborhood = g_QueryNeighborhood;
	if (g_QueryNeighborhood)
		QKmerIndex.m_ptrScoreMx = &ScoreMx;
	}

uint path2posvecs3(
	const char *path, uint ncol,
	uint loQ, uint LQ,
	uint loT, uint LT,
	uint *posQs,
	uint *posTs,
	uint bufnpos)
	{
	if (ncol == 0)
		return 0;
	if (loQ >= LQ || loT >= LT)
		return 0;
	uint posQ = loQ;
	uint posT = loT;
	uint nmatch = 0;
	for (uint col = 0; col < ncol; ++col)
		{
		char c = path[col];
		if (c == 'M')
			{
			if (nmatch >= bufnpos || posQ >= LQ || posT >= LT)
				return 0;
			posQs[nmatch] = posQ;
			posTs[nmatch] = posT;
			++nmatch;
			}
		else if (c != 'D' && c != 'I')
			return 0;
		if (c == 'M' || c == 'D')
			{
			++posQ;
			if (posQ > LQ)
				return 0;
			}
		if (c == 'M' || c == 'I')
			{
			++posT;
			if (posT > LT)
				return 0;
			}
		}
	if (nmatch == 0 || nmatch >= bufnpos)
		return 0;
	return nmatch;
	}

void path2posvecs(
	const string &labelQ, const string &labelT,
	const string &path,
	uint loQ, uint LQ,
	uint loT, uint LT,
	vector<uint> &posQs,
	vector<uint> &posTs)
	{
	posQs.clear();
	posTs.clear();
	const uint colcount = uint(path.size());
	posQs.reserve(colcount);
	posTs.reserve(colcount);
	uint posQ = loQ;
	uint posT = loT;
	for (uint col = 0; col < colcount; ++col)
		{
		char c = path[col];
		if (c == 'M')
			{
			assert(posQ < LQ);
			assert(posT < LT);
			posQs.push_back(posQ);
			posTs.push_back(posT);
			}
		if (c == 'M' || c == 'D')
			posQ++;
		if (c == 'M' || c == 'I')
			posT++;
		}
	}
