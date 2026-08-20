#ifndef x
#error "x not defined"
#endif

/***
Peaker-tunable scalars for flat_hjmega / flat_params::set_scalars.
Every peaker spec must list each of these as constant= or var= (min/max),
except open+ext may be replaced by a single gap2 (ext=open/10).

Also required in peaker specs but not in this X-macro:
  pv   (1=fam, 2=sf, 3=fold) — special-cased in set_scalars.

Alphabet feature weights are NOT listed here; they appear in peaker
specs as var=<alpha_name>;isalpha=yes;[weight=yes;]...

Excluded from peaker (not optimized here):
  - Kappa filter statics (m_kappa_*, pattern, min kmer score, diag modes)
  - Nu SW open/ext (Paralign::m_Open / m_Ext, hardcoded outside flat_params)
  - Log-odds matrices and quantize thresholds (fixed; only weights tune)
  - Geometry statics (m_LDDT_*, m_distmx_bandwidth, m_maxL, ...)
  - m_max_pvalue, m_max_nu_filter_accepts, m_nu_only
***/

x(open,		m_open)
x(ext,		m_ext)
x(selfw,	m_self_w)
x(revw,		m_rev_w)
x(nurevw,	m_nurev_w)
x(dali,		m_dali_w)
x(dalix,	m_dalix_w)
x(lddt,		m_lddt_w)
x(lddtx,	m_lddtx_w)
x(minfwd,	m_mega_filter_min_fwd)
x(nfselfw,	m_nu_filter_self_w)
x(nfrevw,	m_nu_filter_rev_w)
x(nfminfwd,	m_nu_filter_min_fwd_score)
x(nfmincmb,	m_nu_filter_min_combined_score)

#undef x
