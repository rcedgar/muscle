#include "myutils.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "sort.h"

int flat_params::m_kappa_min_kmerpairscore = 65;
int flat_params::m_kappa_min_diagscore = 0;
uint flat_params::m_kappa_min_chainlength =
	flat_params::DEFAULT_MIN_CHAINLENGTH;
string flat_params::m_kappa_pattern = "1010011";

uint flat_params::get_min_chainlength()
	{
	if (optset_minchainlength)
		return opt(minchainlength);
	return DEFAULT_MIN_CHAINLENGTH;
	}

void flat_params::sync_min_chainlength()
	{
	m_kappa_min_chainlength = get_min_chainlength();
	}

uint flat_params::m_kappa_kmer_nrones = KAPPA_NRONES;
uint flat_params::m_kappa_kmer_width = 7;
uint flat_params::m_kappa_dict_size = myipow(KAPPA_AS, KAPPA_NRONES);
uint8_t *flat_params::m_kappa_kmer_onesoffsets;
uint flat_params::m_rsb_size = 1500;
bool flat_params::m_kappa_hsp_rsb_prune = false;
int flat_params::m_kappa_max_pos_logodds = 0;

/***
C:\src\reseek_bench3\bash\onehitdiag_twohitdiag_test.bash

SEPQ0.1=0.298 SEPQ1=0.407 SEPQ10=0.462 Sum3=1.667 sf Kappa filter 13 secs	9.80Gb	01:16	default	
SEPQ0.1=0.298 SEPQ1=0.407 SEPQ10=0.462 Sum3=1.667 sf Kappa filter 19 secs	9.68Gb	01:19	onehitdiag
SEPQ0.1=0.287 SEPQ1=0.383 SEPQ10=0.424 Sum3=1.573 sf Kappa filter  7 secs	9.70Gb	01:02	twohitdiag
 
PDB vs. 100
Kappa filter 39 secs	4.91Gb	2:09	default
Kappa filter 71 secs	4.12Gb	3:19	onehitdiag
Kappa filter 22 secs	4.52Gb	1:48	twohitdiag

AFDB50 vs. 100
Kappa filter 252 secs	10.2Gb	4:48	default
Kappa filter 472 secs	11.3Gb	08:14	onehitdiag
Kappa filter 121 secs	10.2Gb	02:46	twohitdiag
Kappa filter 118 secs	 8.3Gb	02:30	twohitdiag kappa_hsp_rsb_prune 
Kappa filter 113 secs	 8.3Gb	02:30	twohitdiag kappa_hsp_rsb_prune kappa_mindiagscore 200
***/
bool flat_params::m_kappa_onehitdiag = false;
bool flat_params::m_kappa_twohitdiag = false;

bool flat_params::m_hsp_align = true;
bool flat_params::m_hsp_align_check = false;
uint flat_params::m_hsp_align_min_length =
	flat_params::DEFAULT_HSP_ALIGN_MIN_LENGTH;
float flat_params::m_hsp_x2 = 20.0f;

bool flat_params::want_hsp_seeds()
	{
	return m_hsp_align || m_hsp_align_check;
	}

/////////////////////
// Chain quantization
// MUST RE-TRAIN THRESHOLDS AND LOGODDS
///////////////////////////////////////
uint32_t const flat_params::m_nn_min_offset = 12;
uint32_t const flat_params::m_distmx_bandwidth = 256;
uint32_t const flat_params::m_turnd_w = 5;
uint32_t const flat_params::m_angle_n = 4;
uint32_t const flat_params::m_maxL = 4000;
////////////////////////////////////

///////////////////////
// LDDT -- special case
float const flat_params::m_LDDT_R0 = 15;
static const float thresholds[] = { 0.5, 1, 2, 4 };
const float *flat_params::m_LDDT_thresholds = thresholds;
const uint flat_params::m_LDDT_nr_thresholds
	= sizeof(thresholds)/sizeof(thresholds[0]);
///////////////////////////////////////////////

uint flat_params::m_max_nu_filter_accepts = 0;

void flat_params::init_kappa()
	{
	if (optset_rsb_size)
		flat_params::m_rsb_size = opt(rsb_size);
	if (optset_kappa_pattern)
		flat_params::m_kappa_pattern = opt(kappa_pattern);

	uint get_nr_pattern_ones(const string &Str);
	uint k = get_nr_pattern_ones(flat_params::m_kappa_pattern);
	uint K = uint(flat_params::m_kappa_pattern.size());
	flat_params::m_kappa_kmer_onesoffsets = myalloc(uint8_t, k);
	flat_params::m_kappa_kmer_nrones = k; 
	flat_params::m_kappa_kmer_width = K;
	flat_params::m_kappa_dict_size = myipow(KAPPA_AS, k);
	void fill_pattern_offsets(const string &Str, uint8_t *offsets);
	fill_pattern_offsets(flat_params::m_kappa_pattern,
		flat_params::m_kappa_kmer_onesoffsets);

	if (optset_kappa_minkmerscore)
		flat_params::m_kappa_min_kmerpairscore = opt(kappa_minkmerscore);
	if (optset_kappa_mindiagscore)
		flat_params::m_kappa_min_diagscore = opt(kappa_mindiagscore);
	if (optset_kappa_hsp_rsb_prune)
		flat_params::m_kappa_hsp_rsb_prune = true;

	if (opt(onehitdiag) && opt(twohitdiag))
		Die("-onehitdiag and -twohitdiag are mutually exclusive");
	flat_params::m_kappa_onehitdiag = opt(onehitdiag);
	flat_params::m_kappa_twohitdiag = opt(twohitdiag);

	int kappa_max_pos_logodds();
	flat_params::m_kappa_max_pos_logodds = kappa_max_pos_logodds();
	flat_params::sync_min_chainlength();

	if (opt(hsp_align_check))
		flat_params::m_hsp_align_check = true;
	if (optset_mkfl)
		flat_params::m_hsp_align_min_length = opt(mkfl);
	if (optset_hsp_x2)
		flat_params::m_hsp_x2 = float(opt(hsp_x2));
	}

void flat_params::init_from_cmdline()
	{
	if (!optset_stats)
		Die("Must set -stats LEVEL (family, superfamily or fold)");
	uint nmode = 0;
	string smode;
	if (opt(fast))
		{
		flat_params::m_kappa_min_kmerpairscore = 65;
		flat_params::m_rsb_size = 1500;
		smode = "fast";
		}
	else if (opt(sensitive))
		{
		flat_params::m_kappa_min_kmerpairscore = 55;
		flat_params::m_rsb_size = 3000;
		smode = "sensitive";
		}
	else if (opt(verysensitive))
		{
		flat_params::m_kappa_min_kmerpairscore = 50;
		flat_params::m_rsb_size = 5000;
		smode = "verysensitive";
		}
	else
		Die("Must set -fast, -sensitive or -verysensitive");

	sync_min_chainlength();

	if (opt(onehitdiag) && opt(twohitdiag))
		Die("-onehitdiag and -twohitdiag are mutually exclusive");
	flat_params::m_kappa_onehitdiag = opt(onehitdiag);
	flat_params::m_kappa_twohitdiag = opt(twohitdiag);

	if (optset_pvalue)
		{
		m_max_pvalue = opt(pvalue);
		if (m_max_pvalue > 1 || m_max_pvalue <= 0)
			Die("Invalid -pvalue, must be >0 and <= 1");
		}

	const string stats = opt(stats);
	if (stats == "family" || stats == "fam")
		init_from_varstr("=fam");
	else if (stats == "superfamily" || stats == "sf")
		init_from_varstr("=sf");
	else if (stats == "fold")
		init_from_varstr("=fold");
	else
		Die("Invalid -stats '%s', must be family, fam, superfamily, sf or fold", stats.c_str());
	ProgressLog("%s %s\n", smode.c_str(), stats.c_str());
	init_kappa();
	}

void flat_params::init_from_varstr(const string &varstr)
	{
	vector<string> param_names;
	vector<float> param_values;
	parse_varstr(varstr, param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_classify_params(
		param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	const string &alphadir = opt(alphadir);
	set_scalars(scalar_names, scalar_values);
	init_from_alphadir(alphadir, alpha_names);
	uint n = SIZE(alpha_names);
	asserta(SIZE(weights) == n);
	unordered_map<string, float> NameToWeight;
	for (uint i = 0; i < n; ++i)
		{
		const string &name = alpha_names[i];
		if (NameToWeight.find(name) != NameToWeight.end())
			Die("Dupe name in spec '%s'", name.c_str());
		NameToWeight[name] = weights[i];
		}
	apply_weights(NameToWeight);
	}

// non-alpha
void flat_params::set_scalars(
	const vector<string> &names,
	const vector<float> &values)
	{
	assert(names.size() == values.size());

	// test statistic
	m_self_w = FLT_MAX;
	m_rev_w = FLT_MAX;
	m_lddt_w = FLT_MAX;
	m_lddtx_w = FLT_MAX;
	m_dali_w = FLT_MAX;
	m_dalix_w = FLT_MAX;
	m_nurev_w = FLT_MAX;
	m_pvm = PVM_invalid;

	// filters
	m_mega_filter_min_fwd = FLT_MAX;
	m_nu_filter_self_w = FLT_MAX;
	m_nu_filter_rev_w = FLT_MAX;
	m_nu_filter_min_fwd_score = FLT_MAX;
	m_nu_filter_min_combined_score = FLT_MAX;

	for (size_t i = 0; i < names.size(); ++i)
		{
		const string &name = names[i];
		float value = values[i];

		if (name == "gap2")
			{
			m_open = value;
			m_ext = value/10;
			}
		else if (name == "pv")
			{
			if (value == 1)
				m_pvm = PVM_fam;
			else if (value == 2)
				m_pvm = PVM_sf;
			else if (value == 3)
				m_pvm = PVM_fold;
			else
				Die("invalid pvm=%.3g in varstr", value);
			}
#define x(param_name, m_name)	else if (name == #param_name) m_name = value;
#include "tunable_flat_params.h"
		else
			Die("flat_params::setparams() %s=%.3g",
				name.c_str(), value);
		}
	asserta(m_open != FLT_MAX);
	asserta(m_ext != FLT_MAX);
	asserta(m_self_w != FLT_MAX);
	asserta(m_rev_w != FLT_MAX);
	asserta(m_lddt_w != FLT_MAX);
	asserta(m_lddtx_w != FLT_MAX);
	asserta(m_dali_w != FLT_MAX);
	asserta(m_dalix_w != FLT_MAX);
	asserta(m_nurev_w != FLT_MAX);
	asserta(m_mega_filter_min_fwd != FLT_MAX);
	asserta(m_nu_filter_self_w != FLT_MAX);
	asserta(m_nu_filter_rev_w != FLT_MAX);
	asserta(m_nu_filter_min_fwd_score != FLT_MAX);
	asserta(m_nu_filter_min_combined_score != FLT_MAX);
	asserta(m_pvm != PVM_invalid);
	}

bool flat_params::need_distmx()
	{
	return
		flat_params::m_dalix_w > 0 ||
		flat_params::m_dali_w > 0 ||
		flat_params::m_lddt_w > 0 ||
		flat_params::m_lddtx_w > 0;
	}

bool flat_params::need_self()
	{
	return flat_params::m_self_w > 0;
	}

bool flat_params::need_nu_self()
	{
	return flat_params::m_nu_filter_self_w > 0;
	}

bool flat_params::need_reverse()
	{
	return flat_params::m_rev_w > 0;
	}

bool flat_params::need_alignx()
	{
	return
		need_self() ||
		need_reverse();
	}

void flat_params::logme()
	{
	Log("\n");
#define w(x)	Log("%10.3g  %s\n", m_##x, #x)
	w(open);
	w(ext);
	w(self_w);
	w(rev_w);
	w(lddt_w);
	w(lddtx_w);
	w(dali_w);
	w(dalix_w);
	w(nu_filter_self_w);
	w(nu_filter_rev_w);
	w(nu_filter_min_fwd_score);
	w(nu_filter_min_combined_score);
#undef w

#define w(x)	Log("%10d  %s\n", m_##x, #x)
	w(kappa_kmer_nrones);
	w(kappa_kmer_width);
	w(kappa_dict_size);
	w(kappa_min_kmerpairscore);
	w(kappa_min_diagscore);
	w(kappa_min_chainlength);
	w(rsb_size);
	Log("%10d  kappa_hsp_rsb_prune\n", m_kappa_hsp_rsb_prune);
	Log("%10d  kappa_max_pos_logodds\n", m_kappa_max_pos_logodds);
	Log("%10d  kappa_onehitdiag\n", m_kappa_onehitdiag);
	Log("%10d  kappa_twohitdiag\n", m_kappa_twohitdiag);
	Log("%10d  hsp_align\n", m_hsp_align);
	Log("%10d  hsp_align_check\n", m_hsp_align_check);
	Log("%10u  hsp_align_min_length\n", m_hsp_align_min_length);
	Log("%10.3g  hsp_x2\n", m_hsp_x2);
	{
	const char *diag_mode = "unique_fine";
	if (m_kappa_onehitdiag)
		diag_mode = "onehit_insert";
	else if (m_kappa_twohitdiag)
		diag_mode = "twohit_dupes";
	Log("           kappa_diag_mode  %s\n", diag_mode);
	}
#undef w

#define w(x)	Log("%10u  %s\n", m_##x, #x)
	w(nn_min_offset);
	w(distmx_bandwidth);
	w(turnd_w);
	w(angle_n);
#undef x

	Log("LDDT: R0=%.3g thresholds", m_LDDT_R0);
	for (uint i = 0; i < m_LDDT_nr_thresholds; ++i)
		Log(" %.1f", m_LDDT_thresholds[i]);
	Log("\n");
	vector<uint> order(m_nfeat);

	QuickSortOrderDesc(m_weights, m_nfeat, order.data());
	Log("\n");
	Log("%u alphas, sum_sizes=%u\n", m_nfeat, m_sum_alpha_sizes);
	float sumw = 0;
	for (uint k = 0; k < m_nfeat; ++k)
		{
		uint i = order[k];
		float w = m_weights[i];
		sumw += w;
		Log("%10.10s  %7.3f  ", m_alpha_names[i].c_str(), w);
		uint H = uint(w*80);
		for (uint h = 0; h < H; ++h)
			Log("■");
		if (H == 0) Log("o");
		Log("\n");
		}
	Log("%10.10s  %7.3f\n", "TOTAL", sumw);
	}
