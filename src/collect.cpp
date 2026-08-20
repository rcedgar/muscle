#include "myutils.h"
#include "collect.h"
#include "flat_helpers.h"

void collect::from_lines(const vector<string> &lines)
	{
	const size_t N = lines.size();
	size_t i = 0;
	vector<string> flds;
	while (i < N)
		{
		const string &line = lines[i];
		Split(lines[i], flds, '\t');
		asserta(flds.size() == 3);
		asserta(flds[0] == "@");
		const string &name = flds[1];
		uint n = StrToUint(flds[2]);
		asserta(i + n <= N);
		vector<string> name_lines;
		for (uint j = 0; j < n; ++j)
			name_lines.push_back(lines[i+j+1]);
		m_name2lines[name] = name_lines;
		i += n + 1;
		}
	}

void collect::from_file(const string &fn)
	{
	m_name = fn;
	vector<string> lines;
	ReadLinesFromFile(fn, lines);
	from_lines(lines);
	}

const vector<string> &collect::get_lines(const string &name) const
	{
	unordered_map<string, vector<string> >::const_iterator
		iter = m_name2lines.find(name);
	if (iter == m_name2lines.end())
		Die("collect::get_lines(%s)", name.c_str());
	return iter->second;
	}

//void collect::append_kappa()
//	{
//	asserta(optset_sec4_groups);
//	asserta(optset_kappa_logodds);
//	set_sec4_groups(opt(sec4_groups));
//	vector<string> lines;
//	ReadLinesFromFile(opt(kappa_logodds), lines);
//	m_name2lines["kappa32.logodds"] = lines;
//	}
