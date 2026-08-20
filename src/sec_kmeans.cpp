#include "myutils.h"
#include "flat_chain.h"
#include "sec_kmeans.h"
#include "chaq.h"
#include "quarts.h"
#include "alpha.h"

sec_kmeans *sec_kmeans::m_SK2 = 0;
sec_kmeans *sec_kmeans::m_SK3 = 0;
sec_kmeans *sec_kmeans::m_SK4 = 0;
sec_kmeans *sec_kmeans::m_SK8 = 0;
sec_kmeans *sec_kmeans::m_SK16 = 0;
sec_kmeans *sec_kmeans::m_SK32 = 0;

// C:\src\2025-10_reseek_tune\bash\make_sec32.bash
void sec_kmeans::get_sec32_lines(vector<string> &lines)
	{
	lines.clear();
	lines.push_back("sec	32");
	lines.push_back("dim	6");
	lines.push_back("offs1	6	-3	-3	0	1	2	-3");
	lines.push_back("offs2	6	1	0	2	3	0	3");
	lines.push_back("mean	192	234	162	186	188	186	604	257	179	187	191	187	662	1076	619	274	276	274	2347	1077	624	268	240	268	2153	222	177	227	225	227	456	449	433	211	212	211	843	1024	605	264	232	264	1918	799	473	275	269	275	1842	1155	660	283	292	283	2543	291	187	234	248	234	895	454	216	206	214	206	698	267	425	203	206	203	637	424	210	250	251	250	1069	908	528	275	280	275	2089	538	432	239	236	239	1093	656	456	259	246	259	1361	662	392	271	270	271	1616	819	550	195	206	195	528	750	489	209	215	209	819	204	186	221	246	221	246	986	596	228	230	228	1330	816	524	231	227	231	1112	837	529	267	242	267	1567	510	267	269	266	269	1319	984	598	197	205	197	720	1038	611	249	239	249	1665	522	455	230	225	230	546	1007	603	205	212	205	985	248	182	269	245	269	703	483	261	193	210	193	310	275	411	254	252	254	244	614	547	204	249	204	247");
	}

// C:\src\2025-10_reseek_tune\bash\make_sec8.bash
void sec_kmeans::get_sec8_lines(vector<string> &lines)
	{
	lines.clear();
	lines.push_back("sec	8");
	lines.push_back("dim	6");
	lines.push_back("offs1	6	-3	-3	0	1	2	-3");
	lines.push_back("offs2	6	1	0	2	3	0	3");
	lines.push_back("mean	48	244	177	194	195	194	613	400	321	222	226	222	825	943	554	268	253	268	1933	1097	630	276	274	276	2365	829	512	254	245	254	1482	871	558	204	211	204	766	367	328	221	237	221	300	535	344	253	250	253	1201");
	}

// # -sec_kmeans ../data/scop40c.bca -log sec_kmeans2.log -randseed 1  -alpha_size 2 -output sec2.tsv
void sec_kmeans::get_sec2_lines(vector<string> &lines)
	{
	lines.clear();
	lines.push_back("sec	2");
	lines.push_back("dim	6");
	lines.push_back("offs1	6	-3	-3	0	1	2	-3");
	lines.push_back("offs2	6	1	0	2	3	0	3");
	lines.push_back("mean	12	401	289	209	213	209	688	951	561	267	258	267	1932");
	}


// C:\src\2025-10_reseek_tune\2026-03-25_logodds_and_bins\sec_3.kmeans
void sec_kmeans::get_sec3_lines(vector<string> &lines)
	{
	lines.clear();
	lines.push_back("sec	3");
	lines.push_back("dim	5");
	lines.push_back("offs1	5	-2	-2	-1	-1	-3");
	lines.push_back("offs2	5	1	2	1	2	3");
	lines.push_back("mean	15	233	305	197	233	598	599	1024	276	599	2105	475	742	250	474	1162");
	}

// C:\src\2025-10_reseek_tune\2026-03-25_logodds_and_bins\sec_4.kmeans
void sec_kmeans::get_sec4_lines(vector<string> &lines)
	{
	lines.clear();
	lines.push_back("sec	4");
	lines.push_back("dim	6");
	lines.push_back("offs1	6	-3	-3	0	1	-2	-3");
	lines.push_back("offs2	6	1	0	2	3	0	3");
	lines.push_back("mean	24	274	213	203	207	194	591	698	471	218	222	243	860	1053	609	273	265	276	2229	789	482	260	251	255	1544");
	}

// C:\src\2025-10_reseek_tune\2026-03-25_logodds_and_bins\sec_16.kmeans
void sec_kmeans::get_sec16_lines(vector<string> &lines)
	{
	lines.clear();
	lines.push_back("sec	16");
	lines.push_back("dim	5");
	lines.push_back("offs1	5	-3	-3	0	1	-3");
	lines.push_back("offs2	5	1	0	2	3	3");
	lines.push_back("mean	80	246	168	192	193	613	1029	596	271	258	2177	1126	645	279	284	2463	326	200	229	234	832	233	224	228	241	289	533	421	216	216	827	792	470	274	270	1822	1018	601	262	236	1869	863	561	198	207	619	646	403	268	261	1504	936	575	211	215	984	958	583	241	235	1426	435	232	253	256	1132	619	450	245	238	1181	301	429	212	212	635	542	433	208	227	325");
	}

// Conf from reseek v2 converted to flat
void sec_kmeans::get_conf_lines(vector<string> &lines)
	{
	lines.clear();
	lines.push_back("sec	16");
	lines.push_back("dim	9");
	lines.push_back("offs1	9	-2	-2	-2	-1	-1	0	-3	0	-3");
	lines.push_back("offs2	9	0	1	2	1	2	2	3	3	0");
	lines.push_back("mean	144	187	170	248	187	171	188	620	177	169	269	578	977	271	569	264	1989	531	575	295	658	1129	294	655	292	2416	640	645	225	408	677	252	496	268	1484	560	488	210	428	710	256	540	268	929	530	219	271	579	1011	274	605	287	1505	578	369	265	549	765	271	427	201	801	204	499	194	192	278	187	183	189	759	197	429	202	358	545	235	449	221	403	311	204	258	447	350	203	203	188	648	206	517	275	558	894	257	510	221	1356	282	574	187	204	346	195	268	227	900	456	196	194	180	240	202	356	234	287	496	246	255	474	499	211	257	228	1099	461	513	215	246	455	214	433	262	732	514	455	265	358	250	201	194	200	276	331	482");
	}

void sec_kmeans::get_sec_lines(uint alpha_size, vector<string> &lines)
	{
	switch (alpha_size)
		{
	case 2:
		get_sec2_lines(lines); return;
	case 3:
		get_sec3_lines(lines); return;
	case 4:
		get_sec4_lines(lines); return;
	case 8:
		get_sec8_lines(lines); return;
	case 16:
		get_sec16_lines(lines); return;
	case 32:
		get_sec32_lines(lines); return;
		}
	Die("get_sec_lines(%u)", alpha_size);
	}

sec_kmeans *sec_kmeans::get_SK(uint alpha_size, uint M)
	{
	sec_kmeans **ptrSK = 0;
	sec_kmeans *SK = 0;
	switch (alpha_size)
		{
	case 2:		ptrSK = &m_SK2; break;
	case 3:		ptrSK = &m_SK3; break;
	case 4:		ptrSK = &m_SK4; break;
	case 8:		ptrSK = &m_SK8; break;
	case 16:	ptrSK = &m_SK16; break;
	case 32:	ptrSK = &m_SK32; break;
	default:	Die("sec_kmeans::get_SK(%u)", alpha_size);
		}
	if (*ptrSK == 0)
		{
		SK = new sec_kmeans;
		vector<string> lines;
		get_sec_lines(alpha_size, lines);
		SK->from_lines(lines);
		assert(SK->m_K == alpha_size);
		*ptrSK = SK;
		}
	SK = *ptrSK;
	return SK;
	}

void sec_kmeans::from_sec_n(uint alpha_size)
	{
	vector<string> lines;
	switch (alpha_size)
		{
	case 3:		get_sec3_lines(lines); break;
	case 4:		get_sec4_lines(lines); break;
	case 8:		get_sec8_lines(lines); break;
	case 16:	get_sec16_lines(lines); break;
	case 32:	get_sec32_lines(lines); break;
	default: Die("from_sec_n(%u)", alpha_size);
		}
	from_lines(lines);
	}

static bool check_backbone(const flat_chain_t* chain)
	{
	const uint L = chain->get_length();
	for (uint i = 1; i < L; ++i)
		{
		float d = chain->slow_float_dist(i-1, i);
		if (d < 3.7 || d > 3.9)
			return false;
		}
	return true;
	}

void cmd_sec_variance()
	{
	asserta(optset_output);
	FILE *fOut = CreateStdioFile(opt(output));
	vector<flat_chain_t *>chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);
	const uint M = 32;
	const int w = 16;
	vector<vector<float> > dists(w+1);
	for (uint chain_idx = 0; chain_idx < nrchains; ++chain_idx)
		{
		const flat_chain_t* chain = chains[chain_idx];
		if (!check_backbone(chain))
			continue;
		const uint L = chain->get_length();
		for (uint pos = 0; pos < L; ++pos)
		for (uint dij = 1; dij < w; ++dij)
			{
			uint pos2 = pos + dij;
			if (pos2 >= L)
				break;
			float d = chain->slow_float_dist(pos, pos2);
			dists[dij].push_back(d);
			}
		}

	fprintf(fOut, "dij");
	fprintf(fOut, "\tN");
	fprintf(fOut, "\tMin");
	fprintf(fOut, "\tLoQ");
	fprintf(fOut, "\tMed");
	fprintf(fOut, "\tHiQ");
	fprintf(fOut, "\tMax");
	fprintf(fOut, "\tAvg");
	fprintf(fOut, "\tStdDev");
	fprintf(fOut, "\n");
	for (uint dij = 1; dij < w; ++dij)
		{
		const vector<float> &ds = dists[dij];
		QuartsFloat QF;
		GetQuartsFloat(ds, QF);
		ProgressLog("dij=%2u  ", dij);
		QF.ProgressLogMe();

		fprintf(fOut, "%u\t", dij);
		QF.ToTsv(fOut);
		}
	CloseStdioFile(fOut);
	}

void codeseq2fasta(FILE *f, const string &label, const uint8_t *codeseq, uint L)
	{
	if (f == 0)
		return;
	fprintf(f, ">%s\n", label.c_str());
	for (uint i = 0; i < L; ++i)
		{
		if (i > 0 && i%80 == 0)
			fputc('\n', f);
		char c = g_LetterToCharMu[codeseq[i]];
		fputc(c, f);
		}
	fputc('\n', f);
	}

void cmd_sec_fasta()
	{
	const uint M = flat_params::m_distmx_bandwidth;

	vector<flat_chain_t *>chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);

	sec_kmeans SK;
	SK.from_tsv(opt(input));

	FILE *ffa = CreateStdioFile(opt(output));

	for (uint chainidx = 0; chainidx < nrchains; ++chainidx)
		{
		const flat_chain_t* chain = chains[chainidx];
		const uint L = chain->get_length();
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
		uint8_t *codeseq = myalloc(uint8_t, L);
		SK.get_codeseq(distmx, L, codeseq);
		codeseq2fasta(ffa, chain->m_label, codeseq, L);
		}
	CloseStdioFile(ffa);
	}

static void validate_offs(
	const vector<int> &off1s,
	const vector<int> &off2s)
	{
	size_t n = off1s.size();
	asserta(off2s.size() == n);
	for (uint i = 0; i < n; ++i)
		{
		int off1 = off1s[i];
		int off2 = off2s[i];
		int dij = max(off1, off2) - min(off1, off2);
		asserta(dij > 1);
		for (uint j = 0; j < i; ++j)
			{
			if (off1s[j] == off1 && off2s[j] == off2)
				Die("dupe %d,%d", off1, off2);
			}
		}
	}

static void get_rand_off12(uint M, int &off1, int &off2)
	{
	for (uint iter = 0; iter < 100; ++iter)
		{
		uint u = randu32()%(M+1);
		bool plus = (randu32()%2 == 0);
		off1 = (plus ? u : -int(u));

		u = randu32()%(M+1);
		plus = (randu32()%2 == 0);
		off2 = (plus ? u : -int(u));
		int dij = max(off1, off2) - min(off1, off2);
		if (dij > 1)
			return;
		}
	asserta(false);
	}

static bool ok_to_append(
	const vector<int> &off1s,
	const vector<int> &off2s,
	int off1,
	int off2)
	{
	size_t n = off1s.size();
	asserta(off2s.size() == n);
	for (size_t i = 0; i < n; ++i)
		{
		if (off1s[i] == off1 && off2s[i] == off2)
			return false;
		}
	return true;
	}

// r_D_M where D=dimension M=max offset
static void get_random_offs(const string &spec,
	vector<int> &off1s, vector<int> &off2s)
	{
	vector<string> flds;
	Split(spec, flds, '_');
	asserta(flds.size() == 3);
	assert(flds[0] == "r");
	uint D = StrToUint(flds[1]);
	uint M = StrToUint(flds[2]);
	asserta(D > 1 && D < 32);
	asserta(M > 2 && M < 32);
	off1s.clear();
	off2s.clear();

	for (uint d = 0; d < D; ++d)
		{
		bool ok = false;
		for (uint iter = 0; iter < 100; ++iter)
			{
			int off1, off2;
			get_rand_off12(M, off1, off2);
			if (ok_to_append(off1s, off2s, off1, off2))
				{
				off1s.push_back(min(off1,off2));
				off2s.push_back(max(off1,off2));
				ok = true;
				break;
				}
			}
		asserta(ok);
		}
	asserta(off1s.size() == D);
	asserta(off2s.size() == D);
	}

void sec_kmeans::get_codeseq(
	const sid_t *distmx,
	uint L, 
	uint8_t *codeseq) const
	{
	if (int(L) < 2*m_w + 1)
		{
		memset(codeseq, m_K-1, L);
		return;
		}
	sid_t tmpv[32];
	assert(m_D <= 32);

#if DEBUG
	memset(codeseq, UINT8_MAX, L);
#endif

	get_v(distmx, m_w, L, tmpv);
	int8_t letter_lo = assign_cluster(tmpv);
	for (int pos = 0; pos <= m_w; ++pos)
		{
#if DEBUG
		assert(codeseq[pos] == UINT8_MAX);
#endif
		codeseq[pos] = letter_lo;
		}

	int pos_hi = L - m_w - 1;
	for (int pos = m_w + 1; pos < pos_hi; ++pos)
		{
		get_v(distmx, pos, L, tmpv);
#if DEBUG
		assert(codeseq[pos] == UINT8_MAX);
#endif
		codeseq[pos] = assign_cluster(tmpv);
		}

	get_v(distmx, pos_hi, L, tmpv);
	int8_t letter_hi = assign_cluster(tmpv);
	for (int pos = pos_hi; pos < int(L); ++pos)
		{
#if DEBUG
		assert(codeseq[pos] == UINT8_MAX);
#endif
		codeseq[pos] = letter_hi;
		}

#if DEBUG
	for (uint pos = 0; pos < L; ++pos)
		assert(codeseq[pos] < m_K);
#endif
	}

void cmd_sec_kmeans()
	{
	asserta(optset_alpha_size);
	vector<flat_chain_t *>chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);

/***
     _________________________________________
	|     |     |     | [7] |     | o-o |  o  |
	|-3,+3|-2,+3|-1,+3| 0,+3|+1,+3|+2,+3|+3,+3|    o = zero
	|_____|_____|_____|_____|_____|_____|_____|    o-o = 3.81 A (adjacent)
	|     | [2] | [4] | [5] | o-o |  o  | o-o |
	|-3,+2|-2,+2|-1,+2| 0,+2|+1,+2|+2,+2|+3,+2|
	|_____|_____|_____|_____|_____|_____|_____|
	|     | [1] | [3] | o-o |  o  | o-o |     |
	|-3,+1|-2,+1|-1,+1| 0,+1|+1,+1|+2,+1|+3,+1|
	|_____|_____|_____|_____|_____|_____|_____|
	| [8] | [0] | o-o |  o  | o-o |     |     |
	|-3,0 |-2,0 |-1,0 | 0,0 |+1,0 |+2,0 |+3,0 |
	|_____|_____|_____|_____|_____|_____|_____|
	|     | o-o |  o  | o-o |     |     |     |
	|-3,-1|-2,-1|-1,-1| 0,-1|+1,-1|+2,-1|+3,-1|
	|_____|_____|_____|_____|_____|_____|_____|
	| o-o |  o  | o-o |     |     |     |     |
	|-3,-2|-2,-2|-1,-2| 0,-2|+1,-2|+2,-2|+3,-2|
	|_____|_____|_____|_____|_____|_____|_____|
	|  o  | o-o |     |     |     |     | [6] |
	|-3,-3|-2,-3|-1,-3| 0,-3|+1,-3|+2,-3|+3,-3|
	|_____|_____|_____|_____|_____|_____|_____|

***/

#if 0	/////// OLD DEFAULTS
	//                                   0   1   2   3   4   5   6   7   8
	const vector<int> default_off1s = { -2, -2, -2, -1, -1,  0, -3,  0, -3 };
	const vector<int> default_off2s = {  0,  1,  2,  1,  2,  2,  3,  3,  0 };
	//                                   3   3   5   3   3   3   7   4   4  dij
#endif

///////////////////////////////////////////////////////////////
// $src/2025-10_reseek_tune/2026-03-23_sec_variants_best_sweep
// nbr16_16      1.303   -3,1,-3,0,0,2,1,3,2,0,-3,3 <<<<<< best
///////////////////////////////////////////////////////////////
    //                                   0   1   2   3   4   5
	const vector<int> default_off1s = { -3, -3,  0,  1,  2, -3 };
	const vector<int> default_off2s = {  1,  0,  2,  3,  0,  3 };
	//                                   4   3   2   2   2   6  dij

	vector<int> off1s;
	vector<int> off2s;
	uint D = UINT_MAX;
	if (optset_spec)
		{
		const string &spec = opt(spec);
		if (spec[0] == 'r')
			get_random_offs(spec, off1s, off2s);
		else
			{
			vector<string> flds;
			Split(spec, flds, ',');
			size_t n = flds.size();
			asserta(n%2 == 0);
			n /= 2;
			for (uint i = 0; i < n; ++i)
				{
				int off1 = StrToInt(flds[2*i]);
				int off2 = StrToInt(flds[2*i+1]);
				off1s.push_back(off1);
				off2s.push_back(off2);
				}
			}
		}
	else
		{
		off1s = default_off1s;
		off2s = default_off2s;
		}
	validate_offs(off1s, off2s);

	Log("off1s "); for (uint i = 0; i < off1s.size(); ++i) Log(" %3d", off1s[i]); Log("\n");
	Log("off2s "); for (uint i = 0; i < off2s.size(); ++i) Log(" %3d", off2s[i]); Log("\n");
	Log("-spec ");
	for (uint i = 0; i < off1s.size(); ++i)
		{
		if (i > 0) Log(",");
		Log("%d,%d", off1s[i], off2s[i]);
		}
	Log("\n");
	
	const uint K = opt(alpha_size);
	const uint32_t M = 32; // dist mx band width

	uint niter = 1000;
	if (optset_iters) niter = opt(iters);

	sec_kmeans SK;
	SK.init(K, M, off1s, off2s);

	sid_t *vs = myalloc(sid_t, SK.m_N*SK.m_D);
	SK.m_cluster_idxs = myalloc(uint, SK.m_N);
	SK.set_vs(chains);
	SK.train(niter);
	SK.logme();
	SK.ss4stats();
	SK.to_tsv(opt(output));
	if (optset_fasta)
		{
		FILE *ffa = CreateStdioFile(opt(fasta));
		for (uint chainidx = 0; chainidx < nrchains; ++chainidx)
			{
			const flat_chain_t* chain = chains[chainidx];
			const uint L = chain->get_length();
			sid_t *distmx = myalloc(sid_t, L*M);
			chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
			uint8_t *codeseq = myalloc(uint8_t, L);
			SK.get_codeseq(distmx, L, codeseq);
			codeseq2fasta(ffa, chain->m_label, codeseq, L);
			}
		CloseStdioFile(ffa);
		}

#if DEBUG
	sec_kmeans SK2;
	SK2.from_tsv(opt(output));
	SK2.logme();
#endif
	log_flat_stats();
	}

void cmd_flat_secn()
	{
	asserta(optset_alpha_size);
	asserta(optset_fasta);
	vector<flat_chain_t *>chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);
	FILE *ffa = CreateStdioFile(opt(fasta));
	
	const uint alpha_size = opt(alpha_size);
	const uint K = alpha_size;
	const uint M = 32;

	sec_kmeans SK;
	SK.from_sec_n(alpha_size);

	for (uint chainidx = 0; chainidx < nrchains; ++chainidx)
		{
		const flat_chain_t* chain = chains[chainidx];
		const uint L = chain->get_length();
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
		uint8_t *codeseq = myalloc(uint8_t, L);
		SK.get_codeseq(distmx, L, codeseq);
		codeseq2fasta(ffa, chain->m_label, codeseq, L);
		}
	CloseStdioFile(ffa);
	}
