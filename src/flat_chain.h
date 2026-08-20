#pragma once

#include "flat_base.h"
#include "chainaa.h"
#include "chainnu.h"
#include "chainxyz.h"

static const size_t RESERVE_CHAIN_LENGTH = 400;
//static const uint M = 64;

class flat_chain_t : public flat_base<uint32_t, FE_flat_chain>
	{
public:
	uint32_t m_L = 0;
	string m_label = "_null_";
	chainxyz_t *m_xyz = 0;
	chainaa_t *m_aa = 0;
	chainnu_t *m_nu = 0;
	vector<string> m_lines;

protected:
	flat_chain_t(uint32_t L) :
		flat_base<uint32_t, FE_flat_chain>(L)
		{
		if (L == 0)
			{
			m_L = 0;
			m_label = "_null_";
			m_xyz = 0;
			m_aa = 0;
			}
		else
			{
			m_L = L;
#if TRACK_SRC
			m_xyz = chainxyz_t::newflat_src(L, m_srcfile, m_srcline);
			m_aa = chainaa_t::newflat_src(L, m_srcfile, m_srcline);
#else
			m_xyz = chainxyz_t::newflat(L);
			m_aa = chainaa_t::newflat(L);
#endif
			}
		}

public:
	~flat_chain_t();

	void falloc(uint L)
		{
		asserta(m_L == 0);
		m_L = L;
#if TRACK_SRC
		m_aa = chainaa_t::newflat_src(m_L, m_srcfile, m_srcline);
		m_xyz = chainxyz_t::newflat_src(m_L, m_srcfile, m_srcline);
#else
		m_aa = chainaa_t::newflat(m_L);
		m_xyz = chainxyz_t::newflat(m_L);
#endif
		}

	void truncate(uint L);
	void set_xyz(const vector<float> &Xs,
		const vector<float> &Ys, const vector<float> &Zs);
	void set_aa(const vector<char> &aas);
	bool from_pdb_lines(const string &label,
		const vector<string> &lines, bool save_lines);
	void to_fasta(const string &fn) const;
	void to_fasta(FILE *f) const;
	void to_cal(const string &fn) const;
	void to_cal(FILE *f) const;
	void to_pdb(const string &fn, char chainId) const;
	void to_pdb(FILE *f, char chainId) const;

	uint32_t get_length() const
		{
		return m_L;
		}

	void get_ICs(vector<uint16_t> &ICs) const
		{
		uint L = get_length();
		ICs.clear();
		ICs.reserve(3*L);
		const ic_t *data = m_xyz->m_data;
		for (uint pos = 0; pos < L; ++pos)
			{
			uint k = 3*pos;

			ic_t ic_x = data[k];
			ic_t ic_y = data[k+1];
			ic_t ic_z = data[k+2];

			ICs.push_back(ic_x);
			ICs.push_back(ic_y);
			ICs.push_back(ic_z);
			}
		}

	void get_ic_xyz(uint i, ic_t &ic_x, ic_t &ic_y, ic_t &ic_z) const
		{
		const ic_t *data = m_xyz->m_data;
		uint k = 3*i;
		ic_x = data[k];
		ic_y = data[k+1];
		ic_z = data[k+2];
		}

	void get_coords(uint i, float &x, float &y, float &z) const
		{
		const uint16_t *data = m_xyz->m_data;
		uint k = 3*i;
		uint ic_x = data[k];
		uint ic_y = data[k+1];
		uint ic_z = data[k+2];
		x = ic2coord(ic_x);
		y = ic2coord(ic_y);
		z = ic2coord(ic_z);
		}

	uint16_t slow_sd(uint i, uint j) const
		{
		assert(i < m_aa->m_size);
		assert(j < m_aa->m_size);
		const uint16_t *data = m_xyz->m_data;

		int32_t dx = int32_t(data[3*i]) - int32_t(data[3*j]);
		int32_t dy = int32_t(data[3*i+1]) - int32_t(data[3*j+1]);
		int32_t dz = int32_t(data[3*i+2]) - int32_t(data[3*j+2]);
		return dx*dx + dy*dy + dz*dz;
		}

	float slow_float_dist(uint i, uint j) const
		{
		assert(i < m_aa->m_size);
		assert(j < m_aa->m_size);
		const uint16_t *data = m_xyz->m_data;

		int32_t dx = int32_t(data[3*i]) - int32_t(data[3*j]);
		int32_t dy = int32_t(data[3*i+1]) - int32_t(data[3*j+1]);
		int32_t dz = int32_t(data[3*i+2]) - int32_t(data[3*j+2]);
		
		float d2 = float(dx*dx + dy*dy + dz*dz);
		float d = sqrtf(d2)/10;
		return d;
		}

	char get_aa(uint i) const
		{
		assert(i < m_aa->m_size);
		return m_aa->m_data[i];
		}

	bool has_nu() const
		{
		return m_nu != 0;
		}

	const uint8_t *get_nu_data() const
		{
		asserta(m_nu != 0);
		return m_nu->m_data;
		}

	void set_nu_codes(const uint8_t *codes, uint L);

public:
	static uint16_t coord2ic(float X) { return uint16_t((X + 1000)*10 + 0.5); }
	static float ic2coord(uint16_t IC) { return float(IC/10.0f) - 1000; }

#if TRACK_SRC
	static flat_chain_t* newflat_src(uint32_t L,
		const char *srcfile, int srcline)
		{
		flat_chain_t *chain = new flat_chain_t(L);
		chain->m_srcfile = srcfile;
		chain->m_srcline = srcline;

		chain->m_aa->m_srcfile = srcfile;
		chain->m_aa->m_srcline = srcline;

		chain->m_xyz->m_srcfile = srcfile;
		chain->m_xyz->m_srcline = srcline;
		return chain;
		}

	static flat_chain_t* newflat_src(
		const string &label,
		const vector<char> &aas,
		const vector<float> &Xs,
		const vector<float> &Ys,
		const vector<float> &Zs,
		const char *srcfile, int srcline)
		{
		uint32_t L = SIZE(aas);
		assert(SIZE(Xs) == L);
		assert(SIZE(Ys) == L);
		assert(SIZE(Zs) == L);

		flat_chain_t *chain = newflat(L);
		chain->m_label = label;
		chain->set_aa(aas);
		chain->set_xyz(Xs, Ys, Zs);
		chain->m_srcfile = srcfile;
		chain->m_srcline = srcline;
		return chain;
		}
#else
	static flat_chain_t* newflat(uint32_t L)
		{
		flat_chain_t *chain = new flat_chain_t(L);
		return chain;
		}

	static flat_chain_t* newflat(
		const string &label,
		const vector<char> &aas,
		const vector<float> &Xs,
		const vector<float> &Ys,
		const vector<float> &Zs)
		{
		uint32_t L = SIZE(aas);
		assert(SIZE(Xs) == L);
		assert(SIZE(Ys) == L);
		assert(SIZE(Zs) == L);

		flat_chain_t *chain = newflat(L);
		chain->m_label = label;
		chain->set_aa(aas);
		chain->set_xyz(Xs, Ys, Zs);
		return chain;
		}
#endif
	};

extern atomic<uint> g_flat_n_truncated_chains;

uint flat_chain_cap_L(uint L);
void log_flat_n_truncated_chains();

void read_flat_chains(const string &fn, vector<flat_chain_t *> &chains);

void read_flat_chains_idx(
	const string &fn,
	vector<flat_chain_t *> &chains,
	unordered_map<string, uint> &label2idx);

void read_flat_chains_idx_trunclabel(
	const string &fn,
	vector<flat_chain_t *> &chains,
	unordered_map<string, uint> &label2idx);
