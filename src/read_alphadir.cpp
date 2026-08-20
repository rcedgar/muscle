#include "myutils.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "tabbedlines.h"
#include "collect.h"
#include "chaq.h"
#include "fan.h"

uint16_t *read_quantize(const string &fn, uint alpha_size, uint16_t &median)
	{
	tabbedlines tl(fn);
	string s = tl.get_str("fan");
	FAN fan = str2FAN(s);
	uint alpha_size2 = tl.get_int("alpha_size");
	asserta(alpha_size2 == alpha_size);
	median = tl.get_int("median");
	uint16_t *thresholds = tl.get_int16_flat_vec("thresholds", alpha_size-1);
	return thresholds;
	}

uint16_t *quantize_from_lines(const vector<string> &lines,
	FAN &fan, uint &alpha_size, uint16_t &median)
	{
	uint16_t *thresholds = 0;
	alpha_size = UINT_MAX;
	median = UINT16_MAX;
	fan = FAN_COUNT;

	const size_t n = lines.size();
	vector<string> flds;
	for (size_t i = 0; i < n; ++i)
		{
		const string &line = lines[i];
		if (StartsWith(line, "#"))
			continue;
		Split(line, flds, '\t');
		const string &f0 = flds[0];
		if (f0 == "fan")
			{
			asserta(flds.size() == 2);
			fan = str2FAN(flds[1]);
			}
		else if (f0 == "alpha_size")
			{
			asserta(flds.size() == 2);
			alpha_size = StrToUint(flds[1]);
			}
		else if (f0 == "median")
			{
			asserta(flds.size() == 2);
			uint umedian = StrToUint(flds[1]);
			asserta(umedian < UINT16_MAX);
			median = uint16_t(umedian);
			}
		else if (f0 == "thresholds")
			{
			asserta(flds.size() == alpha_size + 1);
			asserta(StrToUint(flds[1]) == alpha_size - 1);
			thresholds = myalloc(uint16_t, alpha_size - 1);
			for (uint j = 0; j < alpha_size - 1; ++j)
				thresholds[j] = StrToUint(flds[j+2]);
			}
		else
			Die("quantize_from_lines(f0=%s)", f0.c_str());
		}

	asserta(thresholds != 0);
	asserta(fan != FAN_COUNT);
	asserta(alpha_size != UINT_MAX);
	asserta(median != UINT16_MAX);

	return thresholds;
	}

void flat_params::init_from_collect(
	const collect &C,
	const vector<string> &alpha_names)
	{
	asserta(!alpha_names.empty());

	set_alpha_names(alpha_names);

	string compound;
	for (uint fi = 0; fi < m_nfeat; ++fi)
		{
		const string &alpha_name = alpha_names[fi];
		uint alpha_size = m_alpha_sizes[fi];
		FAN fan = m_fans[fi];

		if (fi > 0)
			compound += "+";
		compound += alpha_name;

		string logoddsfn;
		Ps(logoddsfn, "%s.logodds", alpha_name.c_str());

		vector<float> logodds;
		const vector<string> &logodds_lines = C.get_lines(logoddsfn);
		uint alpha_size2 = lines2logoddsmx(logodds_lines, logodds);
		asserta(alpha_size2 == alpha_size);

		const uint n = alpha_size*alpha_size;
		m_unweighted_logoddsvec[fi] = myalloc(float, n);
		m_weighted_logoddsvec[fi] = myalloc(float, n);
		for (uint k = 0; k < n; ++k)
			{
			const float score = logodds[k];
			asserta(score >= MIN_SANE_SCORE && score <= MAX_SANE_SCORE);
			m_unweighted_logoddsvec[fi][k] = score;
			m_weighted_logoddsvec[fi][k] = BAD_SCORE;
			}

		uint16_t median = UINT16_MAX;
		uint16_t *thresholds = 0;
		if (is_quantized(fan))
			{
			string quantizefn;
			Ps(quantizefn, "%s.quantize", alpha_name.c_str());
			const vector<string> &quantize_lines = C.get_lines(quantizefn);
			//thresholds = read_quantize(quantizefn, alpha_size, median);
			FAN fan2 = FAN_COUNT;
			uint alpha_size2 = UINT_MAX;
			thresholds = quantize_from_lines(quantize_lines, fan2, alpha_size2, median);
			asserta(fan2 == fan);
			asserta(alpha_size2 == alpha_size);
			}

		m_fans[fi] = fan;
		m_alpha_sizes[fi] = alpha_size;
		m_medians[fi] = median;
		m_thresholds[fi] = thresholds;
		}

	apply_unit_weights();
	set_feature_block_offsets();
	set_symbolsvec();

	Log("%s: %s\n", C.m_name.c_str(), compound.c_str());
	}

void flat_params::init_from_fnprefixes(
	const vector<string> &alpha_names,
	const vector<string> &fnprefixes)
	{
	asserta(!alpha_names.empty());
	asserta(fnprefixes.size() == alpha_names.size());

	set_alpha_names(alpha_names);

	string compound;
	for (uint fi = 0; fi < m_nfeat; ++fi)
		{
		const string &alpha_name = alpha_names[fi];
		const string &fnprefix = fnprefixes[fi];
		FAN fan = m_fans[fi];
		uint alpha_size = m_alpha_sizes[fi];

		if (fi > 0)
			compound += "+";
		compound += alpha_name;

		string logoddsfn;
		Ps(logoddsfn, "%s.logodds", fnprefix.c_str());

		vector<float> logodds;
		uint alpha_size2 = read_logodds(logoddsfn, logodds);
		asserta(alpha_size2 == alpha_size);

		const uint n = alpha_size*alpha_size;
		m_unweighted_logoddsvec[fi] = myalloc(float, n);
		m_weighted_logoddsvec[fi] = myalloc(float, n);
		for (uint k = 0; k < n; ++k)
			{
			const float score = logodds[k];
			asserta(score >= MIN_SANE_SCORE && score <= MAX_SANE_SCORE);
			m_unweighted_logoddsvec[fi][k] = score;
			m_weighted_logoddsvec[fi][k] = BAD_SCORE;
			}

		uint16_t median = UINT16_MAX;
		uint16_t *thresholds = 0;
		if (is_quantized(fan))
			{
			string quantizefn;
			Ps(quantizefn, "%s.quantize", fnprefix.c_str());
			thresholds = read_quantize(quantizefn, alpha_size, median);
			}

		m_fans[fi] = fan;
		m_alpha_sizes[fi] = alpha_size;
		m_medians[fi] = median;
		m_thresholds[fi] = thresholds;
		}

	apply_unit_weights();
	set_feature_block_offsets();
	set_symbolsvec();

	ProgressLog("Loaded %s\n", compound.c_str());
	}

void flat_params::init_from_alphadir(
	const string &arg_alphadir,
	const vector<string> &alpha_names)
	{
	asserta(!alpha_names.empty());

	bool has_kappa = false;
	for (size_t i = 0; i < alpha_names.size(); ++i)
		if (alpha_names[i] == "kappa32")
			{
			has_kappa = true;
			break;
			}
	extern const vector<string> g_alpha_collect_lines;
	if (arg_alphadir == "")
		{
		collect C;
		C.from_lines(g_alpha_collect_lines);
		C.m_name = "[default_alphadir]";
		flat_params::init_from_collect(C, alpha_names);
		return;
		}

	asserta(!has_kappa);
	if (StartsWith(arg_alphadir, "@"))
		{
		const string fn = arg_alphadir.substr(1);
		collect C;
		C.from_file(fn);
		flat_params::init_from_collect(C, alpha_names);
		return;
		}

	string alphadir = arg_alphadir;
	Dirize(alphadir);

	vector<string> fnprefixes;
	for (uint fi = 0; fi < uint(alpha_names.size()); ++fi)
		{
		string fnprefix;
		Ps(fnprefix, "%s%s", alphadir, alpha_names[fi].c_str());
		fnprefixes.push_back(fnprefix);
		}
	init_from_fnprefixes(alpha_names, fnprefixes);
	}

uint16_t chaq::get_undef_value(const flat_params &params, FAN fan, uint alpha_size)
	{
	for (size_t i = 0; i < params.m_nfeat; ++i)
		{
		if (params.m_fans[i] == fan &&
			params.m_alpha_sizes[i] == alpha_size)
			{
			uint16_t median = params.m_medians[i];
			if (median == UINT16_MAX)
				Die("chaq::get_undef_value(%s, %u) median=UINT16_MAX",
					FAN2str(fan), alpha_size);
			return median;
			}
		}
	Die("chaq::get_undef_value(%s, %u) not found",
		FAN2str(fan), alpha_size);
	return UINT16_MAX;
	}

cp_uint16_t chaq::get_thresholds(const flat_params &params, FAN fan, uint alpha_size)
	{
	const size_t n = params.m_fans.size();
	for (size_t i = 0; i < n; ++i)
		{
		if (params.m_fans[i] == fan && params.m_alpha_sizes[i] == alpha_size)
			return params.m_thresholds[i];
		}
	Die("chaq::get_thresholds(%s, %u)", FAN2str(fan), alpha_size);
	return 0;
	}

FAN parse_alpha_name(const string &alpha_name, uint &alpha_size)
	{
	size_t n = alpha_name.size();
	asserta(n > 2);
	string ssize;
	for (int j = int(n)-1; j > 0; --j)
		{
		char c = alpha_name[j];
		if (isdigit(c))
			ssize = c + ssize;
		else
			{
			alpha_size = StrToUint(ssize);
			string fan = alpha_name.substr(0, j+1);
			return str2FAN(fan);
			}
		}
	Die("parse_alpha_name(%s)", alpha_name.c_str());
	return FAN_COUNT;
	}

static void load_alphadir(flat_params &params, const string &arg_alphadir)
	{
	string alphadir = arg_alphadir;
	Dirize(alphadir);

	vector<string> fns;
	vector<bool> subdirs;
	mylistdir(alphadir, fns, subdirs);

	vector<string> alpha_names;
	for (size_t i = 0; i < fns.size(); ++i)
		{
		if (subdirs[i]) continue;
		const string &fn = fns[i];
		if (EndsWith(fn, ".logodds"))
			{
			size_t n = fn.size() - strlen(".logodds");
			string alpha_name = fn.substr(0, n);
			alpha_names.push_back(alpha_name);
			}
		}
	params.init_from_alphadir(alphadir, alpha_names);
	}

void cmd_read_alphadir()
	{
	flat_params params;
	load_alphadir(params, g_Arg1);
	}