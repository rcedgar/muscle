#include "myutils.h"
#include "pdbfilescanner.h"
#include "flat_chain.h"
#include "flat_chain_reader.h"
#include "flat_helpers.h"
#include "flat_params.h"

atomic<uint> g_flat_n_truncated_chains;

uint flat_chain_cap_L(uint L)
	{
	if (L > flat_params::m_maxL)
		{
		g_flat_n_truncated_chains.fetch_add(1, memory_order_relaxed);
		return flat_params::m_maxL;
		}
	return L;
	}

void log_flat_n_truncated_chains()
	{
	uint n = g_flat_n_truncated_chains.load();
	if (n > 0)
		ProgressLogPrefix("%u chains truncated to max length %u\n",
			n, flat_params::m_maxL);
	}

void ChainizeLabel(string &Label, const string &_ChainStr);
void GetThreeFromOne(char aa, string &AAA);

flat_chain_t::~flat_chain_t()
	{
	delete m_aa;
	delete m_xyz;
	delete m_nu;
	m_aa = 0;
	m_xyz = 0;
	m_nu = 0;
	}

static bool GetFieldsFromATOMLine(const string &Line,
  float &X, float &Y, float &Z, char &aa)
	{
	aa = 'X';
	X = -999;
	Y = -999;
	Z = -999;
	string AtomName = Line.substr(12, 4);
	StripWhiteSpace(AtomName);
	if (AtomName != "CA")
		return false;
	char AltLoc = Line[16];
	if (AltLoc != ' ' && AltLoc != 'A' && AltLoc != '1')
		return false;

	string AAA = Line.substr(17, 3);
	aa = GetOneFromThree(AAA);

	string sX, sY, sZ;
	sX = Line.substr(30, 8);
	sY = Line.substr(38, 8);
	sZ = Line.substr(46, 8);

	StripWhiteSpace(sX);
	StripWhiteSpace(sY);
	StripWhiteSpace(sZ);

	X = StrToFloatf(sX);
	Y = StrToFloatf(sY);
	Z = StrToFloatf(sZ);

	return true;
	}

void flat_chain_t::truncate(uint L)
	{
	if (L >= m_L) return;
	if (m_aa) m_aa->truncate(L);
	if (m_xyz) m_xyz->truncate(L);
	if (m_nu) m_nu->truncate(L);
	m_L = L;
	}

void flat_chain_t::set_xyz(const vector<float> &Xs,
	const vector<float> &Ys, const vector<float> &Zs)
	{
	const uint32_t L = SIZE(Xs);
	//m_xyz->falloc(L);
	for (uint32_t i = 0; i < L; ++i)
		{
		uint16_t ic_x = CoordToIC(Xs[i]);
		uint16_t ic_y = CoordToIC(Ys[i]);
		uint16_t ic_z = CoordToIC(Zs[i]);
		m_xyz->set(i, 0, ic_x);
		m_xyz->set(i, 1, ic_y);
		m_xyz->set(i, 2, ic_z);
		}
	}

void flat_chain_t::set_aa(const vector<char> &aas)
	{
	const uint32_t L = SIZE(aas);
	//m_aa->falloc(L);
	memcpy(m_aa->m_data, aas.data(), L);
	}

void flat_chain_t::set_nu_codes(const uint8_t *codes, uint L)
	{
	asserta(L == m_L);
	asserta(codes != 0);
	if (m_nu == 0)
		{
#if TRACK_SRC
		m_nu = chainnu_t::newflat_src(L, m_srcfile, m_srcline);
#else
		m_nu = chainnu_t::newflat(L);
#endif
		}
	memcpy(m_nu->m_data, codes, L);
	}

bool flat_chain_t::from_pdb_lines(const string &label,
	const vector<string> &lines, bool save_lines)
	{
	asserta(m_L == 0);
	if (save_lines)
		m_lines = lines;
	m_label = label;
	const uint N = SIZE(lines);
	uint ResidueCount = 0;
	int CurrentResidueNumber = INT_MAX;
	string ChainStr;
	vector<char> aas;
	vector<float> Xs, Ys, Zs;
	aas.reserve(RESERVE_CHAIN_LENGTH);
	Xs.reserve(RESERVE_CHAIN_LENGTH);
	Ys.reserve(RESERVE_CHAIN_LENGTH);
	Zs.reserve(RESERVE_CHAIN_LENGTH);
	for (uint lineNr = 0; lineNr < N; ++lineNr)
		{
		const string &line = lines[lineNr];
	// First model only; ignore TER (may appear mid-chain before HETATM).
		if (StartsWith(line, "ENDMDL"))
			break;
		const size_t L = line.size();

		char lineChainChar = line[21];
		string lineChainStr;
		lineChainStr.push_back(lineChainChar);
		if (ChainStr == "")
			ChainStr = lineChainStr;
		else if (ChainStr != lineChainStr)
			Die("flat_chain_t::from_pdb_lines() two chains %s, %s",
			  ChainStr.c_str(), lineChainStr.c_str());

		char aa;
		float X, Y, Z;
		bool IsCA = GetFieldsFromATOMLine(line, X, Y, Z, aa);
		if (!IsCA)
			continue;

		aas.push_back(aa);
		Xs.push_back(X);
		Ys.push_back(Y);
		Zs.push_back(Z);
		}
	const uint L = uint(aas.size());
	falloc(L);
	set_xyz(Xs, Ys, Zs);
	set_aa(aas);

	ChainizeLabel(m_label, ChainStr);
	bool Ok = (SIZE(Xs) > 0);
	return Ok;
	}

void read_flat_chains(const string &fn, vector<flat_chain_t *> &chains)
	{
	Progress("Read chains %s ...", fn.c_str());
	PDBFileScanner FS;
	FS.Open(fn);

	flat_chain_reader CR;
	CR.Open(FS);
	for (;;)
		{
		flat_chain_t *chain = CR.GetNext();
		if (!chain)
			break;
		chains.push_back(chain);
		}
	Progress("done\n");
	}

void read_flat_chains_idx(
	const string &fn,
	vector<flat_chain_t *> &chains,
	unordered_map<string, uint> &label2idx)
	{
	read_flat_chains(fn, chains);
	size_t nchain = chains.size();
	for (size_t i = 0; i < nchain; ++i)
		{
		const string &label = chains[i]->m_label;
		label2idx[label] = uint(i);
		}
	}

void read_flat_chains_idx_trunclabel(
	const string &fn,
	vector<flat_chain_t *> &chains,
	unordered_map<string, uint> &label2idx)
	{
	void trunc_label(string &Label);
	read_flat_chains(fn, chains);
	size_t nchain = chains.size();
	for (size_t i = 0; i < nchain; ++i)
		{
		string label = chains[i]->m_label;
		trunc_label(label);
		label2idx[label] = uint(i);
		}
	}

void flat_chain_t::to_fasta(const string &fn) const
	{
	FILE *f = CreateStdioFile(fn);
	to_fasta(f);
	CloseStdioFile(f);
	}

void flat_chain_t::to_cal(const string &fn) const
	{
	FILE *f = CreateStdioFile(fn);
	to_cal(f);
	CloseStdioFile(f);
	}

void flat_chain_t::to_fasta(FILE *f) const
	{
	if (f == 0)
		return;
	string seq(m_aa->m_data, m_aa->m_size);
	SeqToFasta(f, m_label, seq);
	}

void flat_chain_t::to_cal(FILE *f) const
	{
	if (f == 0)
		return;
	const uint L = get_length();
	string buf;
	buf.reserve(m_label.size() + 2 + size_t(L)*40);
	buf.push_back('>');
	buf += m_label;
	buf.push_back('\n');
	char line[64];
	for (uint i = 0; i < L; ++i)
		{
		char aa = get_aa(i);
		float x, y, z;
		get_coords(i, x, y, z);
		int n = snprintf(line, sizeof(line), "%c\t%.1f\t%.1f\t%.1f\n",
			aa, x, y, z);
		asserta(n > 0 && n < (int) sizeof(line));
		buf.append(line, size_t(n));
		}
	WriteStdioFile(f, buf.data(), uint32(buf.size()));
	}

void flat_chain_t::to_pdb(const string &fn, char chainId) const
	{
	if (fn == "")
		return;
	FILE *f = CreateStdioFile(fn);
	to_pdb(f, chainId);
	CloseStdioFile(f);
	}

void flat_chain_t::to_pdb(FILE *f, char chainId) const
	{
	if (f == 0)
		return;
	const uint L = get_length();
	for (uint i = 0; i < L; ++i)
		{
		char aa = get_aa(i);
		string sAAA;
		GetThreeFromOne(aa, sAAA);
		const char *AAA = sAAA.c_str();
		float x, y, z;
		get_coords(i, x, y, z);

		fprintf(f, "ATOM  ");
		fprintf(f, "%5u", i+1);
		fprintf(f, " ");
		fprintf(f, " CA ");
		fprintf(f, " ");
		fprintf(f, "%3.3s", AAA);
		fprintf(f, " ");
		fprintf(f, "%c", chainId);
		fprintf(f, "%4u", i+1);
		fprintf(f, " ");
		fprintf(f, "   ");
		fprintf(f, "%8.3f", x);
		fprintf(f, "%8.3f", y);
		fprintf(f, "%8.3f", z);
		fprintf(f, "%6.2f", 1.0);
		fprintf(f, "%6.2f", 0.0);
		fprintf(f, "          ");
		fprintf(f, " C");
		fprintf(f, "  ");
		fprintf(f, "\n");
		}
	}

#if 0
void cmd_test()
	{
	vector<flat_chain_t *>chains;
	read_flat_chains(g_Arg1, chains);
	FILE *f = CreateStdioFile(opt(output));
	uint n = SIZE(chains);
	ProgressLog("%u chains\n", n);
	for (uint i = 0; i < n; ++i)
		chains[i]->to_cal(f);
	log_flat_stats();
	CloseStdioFile(f);
	}
#endif
