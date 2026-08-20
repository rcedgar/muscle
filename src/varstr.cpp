#include "myutils.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "peaker.h"
#include "fan.h"
#include <unordered_set>

/***
[scalar]     0.89  gap2
[scalar]     0.95  selfw
[scalar]     0.64  revw
[scalar]    0.058  nurevw
[scalar] 0.0001755  dali
[scalar]    0.045  lddt
[scalar]    16.19  minfwd
[scalar]      0.5  nfselfw
[scalar]   0.2899  nfrevw
[scalar]    113.8  nfminfwd
[scalar]    33.82  nfmincmb

   ppack32    0.002  o
    pmdd32    0.008  o
 mendist32    0.008  o
       pm2    0.014  ■
 pendist32    0.015  ■
  nensec32    0.016  ■
 rendist32    0.022  ■
   turnd32    0.035  ■■
     sec32    0.059  ■■■■
   angle32    0.087  ■■■■■■
    pack32    0.088  ■■■■■■■
 nendist32    0.118  ■■■■■■■■■
      aa20    0.528  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
     Total    1.000
***/

/***
secs          sum3  C:\src\reseek_tune2\flat_bench2_nufilter_sweep2\report2.txt
16	1.7760	1.7800	minmufwd140.minmucmb70.minmgfwd20	1.7%
124	1.8110	1.8110	minmufwd80.minmucmb30.minmgfwd15	0.0%
***/

static const string default_varstr =
// "gap2=8.9E-01;selfw=9.5E-01;revw=6.4E-01;nurevw=5.8E-02;dali=0.00017545;lddt=4.5E-02;minfwd=16.191;nfselfw=5.0E-01;nfrevw=0.28992;nfminfwd=113.85;nfmincmb=33.82;aa20=5.2848E-01;angle32=8.6751E-02;mendist32=7.9881E-03;nendist32=1.1838E-01;nensec32=1.5953E-02;pack32=8.8342E-02;pendist32=1.4577E-02;pm2=1.3725E-02;pmdd32=7.9770E-03;ppack32=1.6439E-03;rendist32=2.1936E-02;sec32=5.9344E-02;turnd32=3.4900E-02;";
"pv=2;gap2=8.9E-01;selfw=9.5E-01;revw=6.4E-01;nurevw=5.8E-02;dali=0.00017545;lddt=4.5E-02;minfwd=16.191;nfselfw=5.0E-01;nfrevw=0.28992;nfminfwd=113.85;nfmincmb=40;aa20=5.2848E-01;angle32=8.6751E-02;mendist32=7.9881E-03;nendist32=1.1838E-01;nensec32=1.5953E-02;pack32=8.8342E-02;pendist32=1.4577E-02;pm2=1.3725E-02;pmdd32=7.9770E-03;ppack32=1.6439E-03;rendist32=2.1936E-02;sec32=5.9344E-02;turnd32=3.4900E-02;lddtx=0;dalix=0;";

static const string fast_varstr_sf =
// "gap2=8.9E-01;selfw=9.5E-01;revw=6.4E-01;nurevw=5.8E-02;dali=0.00017545;lddt=4.5E-02;minfwd=16.191;nfselfw=5.0E-01;nfrevw=0.28992;nfminfwd=113.85;nfmincmb=33.82;aa20=5.2848E-01;angle32=8.6751E-02;mendist32=7.9881E-03;nendist32=1.1838E-01;nensec32=1.5953E-02;pack32=8.8342E-02;pendist32=1.4577E-02;pm2=1.3725E-02;pmdd32=7.9770E-03;ppack32=1.6439E-03;rendist32=2.1936E-02;sec32=5.9344E-02;turnd32=3.4900E-02;";
"pv=2;minfwd=20;nfminfwd=140;nfmincmb=70;gap2=8.9E-01;selfw=9.5E-01;revw=6.4E-01;nurevw=5.8E-02;dali=0.00017545;lddt=4.5E-02;nfselfw=5.0E-01;nfrevw=0.28992;aa20=5.2848E-01;angle32=8.6751E-02;mendist32=7.9881E-03;nendist32=1.1838E-01;nensec32=1.5953E-02;pack32=8.8342E-02;pendist32=1.4577E-02;pm2=1.3725E-02;pmdd32=7.9770E-03;ppack32=1.6439E-03;rendist32=2.1936E-02;sec32=5.9344E-02;turnd32=3.4900E-02;lddtx=0;dalix=0;";
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
static const string best_varstr_fold =
"pv=3;minfwd=16.191;nfselfw=5.0E-01;nfrevw=0.28992;nfminfwd=113.85;nfmincmb=40;aa20=0.0001;pm2=0.0001;revw=9.96E-01;selfw=9.16E-01;gap2=7.07E-01;lddt=1.10E-01;dali=6.52E-03;nensec32=2.71E-01;aa4=1.84E-01;sec32=1.00E-01;nendist32=8.82E-02;mendist32=1.13E-01;ppack32=1.01E-01;pmdd32=1.32E-02;mpack32=6.72E-02;angle32=1.65E-02;pendist32=4.49E-02;lddtx=0;dalix=0;nurevw=0;";

static const string fast_varstr_fold =
"pv=3;minfwd=20;nfselfw=5.0E-01;nfrevw=0.28992;nfminfwd=140;nfmincmb=70;aa20=0.0001;pm2=0.0001;revw=9.96E-01;selfw=9.16E-01;gap2=7.07E-01;lddt=1.10E-01;dali=6.52E-03;nensec32=2.71E-01;aa4=1.84E-01;sec32=1.00E-01;nendist32=8.82E-02;mendist32=1.13E-01;ppack32=1.01E-01;pmdd32=1.32E-02;mpack32=6.72E-02;angle32=1.65E-02;pendist32=4.49E-02;lddtx=0;dalix=0;nurevw=0;";
//^^^^^^^^                                ^^^^^^^^^^^^^^^^^^^^^^^^

// TODO nurevw=0 added by hand for fam, seems not optimized in TS
static const string best_varstr_fam =
//"minfwd=16.191;nfselfw=5.0E-01;nfrevw=0.28992;nfminfwd=113.85;nfmincmb=40;selfw=6.49E-01;lddt=1.05E+00;revw=6.36E-01;gap2=4.48E-01;dali=3.47E-04;aa20=6.07E-01;nendist32=1.04E-01;mendist32=1.03E-01;sec32=4.53E-02;pendist32=3.09E-02;pmdiff32=3.19E-02;nensec32=3.57E-02;turnd32=2.05E-02;mensec32=2.15E-02;";
"pv=1;minfwd=16.191;nfselfw=5.0E-01;nurevw=0;nfrevw=0.28992;nfminfwd=113.85;nfmincmb=40;pm2=0.0001;selfw=6.49E-01;lddt=1.05E+00;revw=6.36E-01;gap2=4.48E-01;dali=3.47E-04;aa20=6.07E-01;nendist32=1.04E-01;mendist32=1.03E-01;sec32=4.53E-02;pendist32=3.09E-02;pmdiff32=3.19E-02;nensec32=3.57E-02;turnd32=2.05E-02;mensec32=2.15E-02;lddtx=0;dalix=0";

static const string fast_varstr_fam =
//"minfwd=16.191;nfselfw=5.0E-01;nfrevw=0.28992;nfminfwd=113.85;nfmincmb=40;selfw=6.49E-01;lddt=1.05E+00;revw=6.36E-01;gap2=4.48E-01;dali=3.47E-04;aa20=6.07E-01;nendist32=1.04E-01;mendist32=1.03E-01;sec32=4.53E-02;pendist32=3.09E-02;pmdiff32=3.19E-02;nensec32=3.57E-02;turnd32=2.05E-02;mensec32=2.15E-02;";
"pv=1;minfwd=20;nfselfw=5.0E-01;nurevw=0;nfrevw=0.28992;nfminfwd=140;nfmincmb=70;pm2=0.0001;selfw=6.49E-01;lddt=1.05E+00;revw=6.36E-01;gap2=4.48E-01;dali=3.47E-04;aa20=6.07E-01;nendist32=1.04E-01;mendist32=1.03E-01;sec32=4.53E-02;pendist32=3.09E-02;pmdiff32=3.19E-02;nensec32=3.57E-02;turnd32=2.05E-02;mensec32=2.15E-02;lddtx=0;dalix=0;";
//^^^^^^^^                                         ^^^^^^^^^^^^^^^^^^^^^^^^

void parse_varstr(
	const string &arg_VarStr,
	vector<string> &Names,
	vector<float> &Values)
	{
	Names.clear();
	Values.clear();

	string VarStr;
	if (arg_VarStr == "" || arg_VarStr == "=sf")
		VarStr = default_varstr;
	else if (arg_VarStr == "=fastsf")
		VarStr = fast_varstr_sf;
	else if (arg_VarStr == "=fold")
		VarStr = best_varstr_fold;
	else if (arg_VarStr == "=fastfold")
		VarStr = fast_varstr_fold;
	else if (arg_VarStr == "=fam")
		VarStr = best_varstr_fam;
	else if (arg_VarStr == "=fastfam")
		VarStr = fast_varstr_fam;
	else
		{
		if (StartsWith(arg_VarStr, "@"))
			{
			const string fn = arg_VarStr.substr(1);
			vector<string> lines;
			ReadLinesFromFile(fn, lines);
			for (auto line : lines)
				{
				if (StartsWith(line, "#"))
					continue;
				StripWhiteSpace(line);
				asserta(EndsWith(line, ";"));
				VarStr += line;
				}
			}
		else
			VarStr = arg_VarStr;
		}

	StripAllWhiteSpace(VarStr);

	vector<string> Fields;
	Split(VarStr, Fields, ';');

	const uint n = SIZE(Fields);
	for (uint i = 0; i < n; ++i)
		{
		const string &NameEqValue = Fields[i];
		vector<string> Fields2;
		Split(NameEqValue, Fields2, '=');
		if (SIZE(Fields2) != 2)
			Die("SubsetBench::parse_varstr(%s) not name=value '%s'",
				VarStr.c_str(), Fields[i].c_str());
		const string &Name = Fields2[0];
		const string &ValueStr = Fields2[1];
		float Value = StrToFloatf(ValueStr);
		Names.push_back(Name);
		Values.push_back(Value);
		}
	}

void flat_classify_params(
	const vector<string> &names,
	const vector<float> &values,
	vector<string> &alpha_names,
	vector<float> &alpaha_weights,
	vector<string> &scalar_names,
	vector<float> &scalar_values)
	{
	for (uint i = 0; i < SIZE(names); ++i)
		{
		const string &name = names[i];
		float Value = values[i];
		bool is_scalar = false;

		if (name == "gap2" || name == "pv") { is_scalar = true; }
#define x(param_name, m_name)	else if (name == #param_name) {is_scalar = true; }
#include "tunable_flat_params.h"
		if (is_scalar)
			{
			scalar_names.push_back(name);
			scalar_values.push_back(Value);
			}
		else
			{
			uint alpha_size = 0;
			FAN fan = parse_alpha_name(name, alpha_size);
			alpha_names.push_back(name);
			alpaha_weights.push_back(Value);
			}
		}
	}

void flat_make_varstr(const flat_params &params, string &varstr)
	{
	varstr.clear();

	if (params.m_pvm == PVM_fam)
		Psa(varstr, "pv=1;\n");
	else if (params.m_pvm == PVM_sf)
		Psa(varstr, "pv=2;\n");
	else if (params.m_pvm == PVM_fold)
		Psa(varstr, "pv=3;\n");
	else
		Die("flat_make_varstr: invalid pvm");

	if (feq(params.m_open, params.m_ext*10))
		Psa(varstr, "gap2=%.4g;\n", params.m_open);
	else
		{
		Psa(varstr, "open=%.4g;\n", params.m_open);
		Psa(varstr, "ext=%.4g;\n", params.m_ext);
		}

#define x(param_name, member_name)	\
	if (string(#param_name) != "open" && string(#param_name) != "ext") \
		Psa(varstr, "%s=%.4g;\n", #param_name, params.member_name);
#include "tunable_flat_params.h"

	for (uint fi = 0; fi < params.m_nfeat; ++fi)
		{
		Psa(varstr, "%s=%.4g;\n",
			params.m_alpha_names[fi].c_str(),
			params.m_weights[fi]);
		}
	}

static float pvm_to_float(PVALUE_MODE pvm)
	{
	if (pvm == PVM_fam) return 1;
	if (pvm == PVM_sf) return 2;
	if (pvm == PVM_fold) return 3;
	Die("pvm_to_float: invalid pvm");
	return 0;
	}

void flat_make_peaker_spec_const(
	const flat_params &params, vector<string> &lines)
	{
	lines.clear();

	string line;
	Ps(line, "var=pv;constant=%.4g;", pvm_to_float(params.m_pvm));
	lines.push_back(line);

	if (feq(params.m_open, params.m_ext*10))
		{
		Ps(line, "var=gap2;constant=%.4g;", params.m_open);
		lines.push_back(line);
		}
	else
		{
		Ps(line, "var=open;constant=%.4g;", params.m_open);
		lines.push_back(line);

		Ps(line, "var=ext;constant=%.4g;", params.m_ext);
		lines.push_back(line);
		}

#define x(param_name, member_name)	\
	if (string(#param_name) != "open" && string(#param_name) != "ext") { \
		Ps(line, "var=%s;constant=%.4g;", #param_name, params.member_name); \
		lines.push_back(line); }
#include "tunable_flat_params.h"

	for (uint fi = 0; fi < params.m_nfeat; ++fi)
		{
		Ps(line, "var=%s;constant=%.4g;isalpha=yes;weight=yes;",
			params.m_alpha_names[fi].c_str(),
			params.m_weights[fi]);
		lines.push_back(line);
		}
	}

void flat_make_peaker_spec_range(
	const flat_params &params,
	vector<string> &lines,
	double rate)
	{
	lines.clear();

	string line;
	// pv is discrete; keep constant in range templates
	Ps(line, "var=pv;constant=%.4g;", pvm_to_float(params.m_pvm));
	lines.push_back(line);

	if (feq(params.m_open, params.m_ext*10))
		{
		Ps(line, "var=gap2;min=%.4g;max=%.4g;",
			params.m_open/rate, params.m_open*rate);
		lines.push_back(line);
		}
	else
		{
		Ps(line, "var=open;min=%.4g;max=%.4g;",
			params.m_open/rate, params.m_open*rate);
		lines.push_back(line);

		Ps(line, "var=ext;min=%.4g;max=%.4g;",
			params.m_ext/rate, params.m_ext*rate);
		lines.push_back(line);
		}

#define x(param_name, member_name)	\
	if (string(#param_name) != "open" && string(#param_name) != "ext") { \
		if (params.member_name == 0) \
			Ps(line, "var=%s;constant=0;", #param_name); \
		else \
			Ps(line, "var=%s;min=%.4g;max=%.4g;", #param_name, \
		params.member_name/rate, params.member_name*rate); \
		lines.push_back(line); }
#include "tunable_flat_params.h"

	for (uint fi = 0; fi < params.m_nfeat; ++fi)
		{
		float w = params.m_weights[fi];
		if (w == 0)
			Ps(line, "var=%s;constant=0;isalpha=yes;weight=yes;",
				params.m_alpha_names[fi].c_str());
		else
			Ps(line, "var=%s;min=%.4g;max=%.4g;isalpha=yes;weight=yes;",
				params.m_alpha_names[fi].c_str(),
				w/rate, w*rate);
		lines.push_back(line);
		}
	}

static bool is_known_scalar_name(const string &name)
	{
	if (name == "gap2" || name == "pv")
		return true;
#define x(param_name, member_name) if (name == #param_name) return true;
#include "tunable_flat_params.h"
	return false;
	}

void validate_flat_peaker_spec(
	const vector<string> &SpecLines,
	vector<string> &alpha_names)
	{
	alpha_names.clear();

	unordered_set<string> seen_names;
	unordered_set<string> scalar_names;
	bool has_gap2 = false;
	bool has_open = false;
	bool has_ext = false;

	for (uint i = 0; i < SIZE(SpecLines); ++i)
		{
		const string &Line = SpecLines[i];
		if (Line.empty() || StartsWith(Line, "#"))
			continue;
		if (!StartsWith(Line, "var="))
			continue;

		string name;
		Peaker::SpecGetStr(Line, "var", name, "");
		if (name.empty())
			Die("peaker spec line missing var= name: %s", Line.c_str());
		if (seen_names.find(name) != seen_names.end())
			Die("Duplicate peaker var '%s'", name.c_str());
		seen_names.insert(name);

		const bool isalpha = Peaker::SpecGetBool(Line, "isalpha", false);
		const bool isweight = Peaker::SpecGetBool(Line, "weight", false);
		string constant;
		Peaker::SpecGetStr(Line, "constant", constant, "");
		string minv;
		string maxv;
		Peaker::SpecGetStr(Line, "min", minv, "");
		Peaker::SpecGetStr(Line, "max", maxv, "");
		const bool has_const = (constant != "");
		const bool has_range = (minv != "" && maxv != "");
		if (has_const == has_range)
			Die("Peaker var '%s' must have either constant= or min=+max= (not both/neither)",
				name.c_str());
		if ((minv != "") != (maxv != ""))
			Die("Peaker var '%s' needs both min= and max=", name.c_str());

		if (isalpha)
			{
			if (is_known_scalar_name(name))
				Die("Peaker var '%s' is a scalar but has isalpha=yes", name.c_str());
			uint alpha_size = 0;
			FAN fan = parse_alpha_name(name, alpha_size);
			if (fan == FAN_COUNT || alpha_size == 0)
				Die("Invalid alphabet name in peaker spec '%s'", name.c_str());
			alpha_names.push_back(name);
			}
		else
			{
			if (isweight)
				Die("weight=yes is only legal with isalpha=yes (var=%s)", name.c_str());
			if (!is_known_scalar_name(name))
				Die("Unknown peaker scalar '%s' (alphabet names require isalpha=yes)",
					name.c_str());
			scalar_names.insert(name);
			if (name == "gap2") has_gap2 = true;
			if (name == "open") has_open = true;
			if (name == "ext") has_ext = true;
			}
		}

	if (SIZE(alpha_names) == 0)
		Die("Peaker spec must include at least one alphabet (isalpha=yes)");

	if (has_gap2)
		{
		if (has_open || has_ext)
			Die("Peaker spec: use gap2 or open+ext, not both");
		}
	else
		{
		if (!has_open || !has_ext)
			Die("Peaker spec: missing gap2, or missing open/ext");
		}

	if (scalar_names.find("pv") == scalar_names.end())
		Die("Peaker spec missing required scalar pv");

#define x(param_name, member_name)	\
	if (string(#param_name) != "open" && string(#param_name) != "ext") { \
		if (scalar_names.find(#param_name) == scalar_names.end()) \
			Die("Peaker spec missing required scalar %s", #param_name); }
#include "tunable_flat_params.h"

	for (unordered_set<string>::const_iterator it = scalar_names.begin();
		it != scalar_names.end(); ++it)
		{
		if (!is_known_scalar_name(*it))
			Die("Unexpected scalar '%s'", it->c_str());
		}
	}

