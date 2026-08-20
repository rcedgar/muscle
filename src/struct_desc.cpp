#include "myutils.h"
#include "struct_desc.h"
#include <cctype>
#include <set>

static const char *g_UnpDbs[] =
	{
	"UNP", "UNIPROT", "SWS", "SWISSPROT", "SWISS-PROT", 0
	};

bool CifValuePresent(const string &value)
	{
	string s = value;
	StripWhiteSpace(s);
	return s != "" && s != "." && s != "?";
	}

void CleanStructText(string &s)
	{
	string out;
	bool in_ws = true;
	for (size_t i = 0; i < s.size(); ++i)
		{
		unsigned char c = (unsigned char) s[i];
		if (isspace(c))
			{
			if (!in_ws && !out.empty())
				{
				out.push_back(' ');
				in_ws = true;
				}
			}
		else
			{
			out.push_back((char) c);
			in_ws = false;
			}
		}
	if (!out.empty() && out[out.size() - 1] == ' ')
		out.resize(out.size() - 1);
	if (out == "." || out == "?")
		out.clear();
	s = out;
	}

string CleanStructTextCopy(const string &s)
	{
	string t = s;
	CleanStructText(t);
	return t;
	}

static bool IsUnpDb(const string &db_name)
	{
	string u = db_name;
	ToUpper(u);
	for (unsigned i = 0; g_UnpDbs[i]; ++i)
		{
		if (u == g_UnpDbs[i])
			return true;
		}
	return false;
	}

string FormatDbRef(const string &db_name, const string &accession)
	{
	string db = db_name;
	string acc = accession;
	StripWhiteSpace(db);
	StripWhiteSpace(acc);
	if (!CifValuePresent(db) && !CifValuePresent(acc))
		return "";
	if (!CifValuePresent(acc))
		return "";
	if (IsUnpDb(db) || !CifValuePresent(db))
		return "UNP:" + acc;
	return db + ":" + acc;
	}

string PreferUnp(const vector<string> &refs)
	{
	if (refs.empty())
		return "";
	for (size_t i = 0; i < refs.size(); ++i)
		{
		if (StartsWith(refs[i], "UNP:"))
			return refs[i];
		}
	return refs[0];
	}

void AppendStructDescToLabel(string &Label, const string &Entry,
	const string &DbRef, const string &Molecule, const string &Title)
	{
	if (opt(trunclabels))
		return;

	const string seq_id = Label;
	string db_ref = CleanStructTextCopy(DbRef);
	string molecule = CleanStructTextCopy(Molecule);
	string title = CleanStructTextCopy(Title);

	if (!db_ref.empty())
		{
		Label += " ";
		Label += db_ref;
		}
	if (!molecule.empty())
		{
		Label += " ";
		Label += molecule;
		}
	if (!title.empty())
		{
		string title_lc = title;
		string mol_lc = molecule;
		ToLower(title_lc);
		ToLower(mol_lc);
		if (title_lc != mol_lc)
			{
			string sid_lc = seq_id;
			string entry_lc = Entry;
			ToLower(sid_lc);
			ToLower(entry_lc);
			if (sid_lc.find(title_lc) == string::npos && title_lc != entry_lc)
				{
				Label += " ";
				Label += title;
				}
			}
		}
	}

string PickRefByEntity(const map<string, vector<string> > &by_entity,
	const string &entity_id)
	{
	if (!entity_id.empty())
		{
		map<string, vector<string> >::const_iterator it =
			by_entity.find(entity_id);
		if (it != by_entity.end())
			return PreferUnp(it->second);
		}
	map<string, vector<string> >::const_iterator it_empty =
		by_entity.find("");
	if (it_empty != by_entity.end())
		return PreferUnp(it_empty->second);
	if (by_entity.size() == 1)
		return PreferUnp(by_entity.begin()->second);
	return "";
	}

string PickMoleculeByEntity(const map<string, string> &mol_by_entity,
	const string &entity_id)
	{
	if (!entity_id.empty())
		{
		map<string, string>::const_iterator it = mol_by_entity.find(entity_id);
		if (it != mol_by_entity.end())
			return it->second;
		}
	if (mol_by_entity.size() == 1)
		return mol_by_entity.begin()->second;
	return "";
	}

string PickMoleculeByChain(const map<string, string> &mol_by_chain,
	const string &chain)
	{
	map<string, string>::const_iterator it = mol_by_chain.find(chain);
	if (it != mol_by_chain.end())
		return it->second;
	if (mol_by_chain.empty())
		return "";
	set<string> uniq;
	for (map<string, string>::const_iterator i = mol_by_chain.begin();
	  i != mol_by_chain.end(); ++i)
		uniq.insert(i->second);
	if (uniq.size() == 1)
		return *uniq.begin();
	return "";
	}

/*** CIF token stream (metadata categories only) ***/

static void SplitCifTag(const string &tag, string &cat, string &item)
	{
	cat.clear();
	item.clear();
	string t = tag;
	if (!t.empty() && t[0] == '_')
		t = t.substr(1);
	size_t dot = t.find('.');
	if (dot == string::npos)
		{
		cat = t;
		return;
		}
	cat = t.substr(0, dot);
	item = t.substr(dot + 1);
	}

static bool IsCifTag(const string &tok)
	{
	return tok.size() >= 2 && tok[0] == '_' && tok.find('.') != string::npos;
	}

static bool IsCifStopTok(const string &tok, bool eof)
	{
	if (eof)
		return true;
	if (tok == "loop_" || tok == "stop_")
		return true;
	if (StartsWith(tok, "data_") || StartsWith(tok, "save_") ||
	  StartsWith(tok, "global_"))
		return true;
	if (IsCifTag(tok))
		return true;
	return false;
	}

static void TokenizeCifLine(const string &line, vector<string> &toks)
	{
	toks.clear();
	const size_t n = line.size();
	size_t i = 0;
	while (i < n)
		{
		unsigned char c = (unsigned char) line[i];
		if (isspace(c))
			{
			++i;
			continue;
			}
		if (c == '#')
			return;
		if (c == '\'' || c == '"')
			{
			char quote = (char) c;
			++i;
			size_t start = i;
			while (i < n && line[i] != quote)
				++i;
			toks.push_back(line.substr(start, i - start));
			if (i < n)
				++i;
			continue;
			}
		size_t start = i;
		while (i < n && !isspace((unsigned char) line[i]) && line[i] != '#')
			++i;
		toks.push_back(line.substr(start, i - start));
		if (i < n && line[i] == '#')
			return;
		}
	}

class CifLineTokStream
	{
public:
	const vector<string> &m_Lines;
	uint m_LineIdx = 0;
	vector<string> m_Buf;
	uint m_BufIdx = 0;
	string m_Peek;
	bool m_HasPeek = false;
	bool m_Eof = false;

	explicit CifLineTokStream(const vector<string> &Lines)
		: m_Lines(Lines)
		{
		}

	bool Next(string &tok)
		{
		if (m_HasPeek)
			{
			tok = m_Peek;
			m_HasPeek = false;
			return true;
			}
		for (;;)
			{
			if (m_BufIdx < SIZE(m_Buf))
				{
				tok = m_Buf[m_BufIdx++];
				return true;
				}
			if (m_LineIdx >= SIZE(m_Lines))
				{
				m_Eof = true;
				tok.clear();
				return false;
				}
			const string &raw = m_Lines[m_LineIdx++];
			string line = raw;
			// strip CR
			if (!line.empty() && line[line.size() - 1] == '\r')
				line.resize(line.size() - 1);
			if (line.empty())
				continue;
			if (line[0] == ';')
				{
				vector<string> chunks;
				chunks.push_back(line.substr(1));
				for (;;)
					{
					if (m_LineIdx >= SIZE(m_Lines))
						{
						tok.clear();
						for (size_t i = 0; i < chunks.size(); ++i)
							{
							if (i)
								tok.push_back('\n');
							tok += chunks[i];
							}
						m_Buf.clear();
						m_BufIdx = 0;
						return true;
						}
					string more = m_Lines[m_LineIdx++];
					if (!more.empty() && more[more.size() - 1] == '\r')
						more.resize(more.size() - 1);
					if (!more.empty() && more[0] == ';')
						{
						tok.clear();
						for (size_t i = 0; i < chunks.size(); ++i)
							{
							if (i)
								tok.push_back('\n');
							tok += chunks[i];
							}
						m_Buf.clear();
						m_BufIdx = 0;
						string rest = more.substr(1);
						StripWhiteSpace(rest);
						if (!rest.empty())
							TokenizeCifLine(rest, m_Buf);
						return true;
						}
					chunks.push_back(more);
					}
				}
			string stripped = line;
			size_t p = 0;
			while (p < stripped.size() && isspace((unsigned char) stripped[p]))
				++p;
			if (p < stripped.size() && stripped[p] == '#')
				continue;
			TokenizeCifLine(line, m_Buf);
			m_BufIdx = 0;
			}
		}

	bool Peek(string &tok)
		{
		if (m_HasPeek)
			{
			tok = m_Peek;
			return true;
			}
		if (!Next(tok))
			return false;
		m_Peek = tok;
		m_HasPeek = true;
		return true;
		}
	};

static bool KeepCifCategory(const string &cat)
	{
	return cat == "entry" || cat == "struct" || cat == "entity" ||
	  cat == "struct_ref" || cat == "ma_target_ref_db_details";
	}

static void ParseCifLoop(CifLineTokStream &ts,
	map<string, vector<map<string, string> > > &categories)
	{
	vector<string> fields;
	string tok;
	for (;;)
		{
		if (!ts.Peek(tok) || !IsCifTag(tok))
			break;
		ts.Next(tok);
		fields.push_back(tok);
		}
	if (fields.empty())
		return;

	string cat0, item0;
	SplitCifTag(fields[0], cat0, item0);
	vector<string> items;
	items.reserve(fields.size());
	for (size_t i = 0; i < fields.size(); ++i)
		{
		string cat, item;
		SplitCifTag(fields[i], cat, item);
		items.push_back(item);
		}
	const uint n = SIZE(fields);
	const bool keep = KeepCifCategory(cat0);
	const bool skip_atom = (cat0 == "atom_site");

	vector<map<string, string> > rows;
	for (;;)
		{
		string peek;
		bool have = ts.Peek(peek);
		if (IsCifStopTok(peek, !have))
			break;

		vector<string> row_vals;
		row_vals.reserve(n);
		for (uint i = 0; i < n; ++i)
			{
			string peek2;
			bool have2 = ts.Peek(peek2);
			if (IsCifStopTok(peek2, !have2) && row_vals.empty())
				{
				if (keep && !rows.empty())
					categories[cat0].insert(categories[cat0].end(),
					  rows.begin(), rows.end());
				return;
				}
			string val;
			if (!ts.Next(val))
				{
				if (keep && !rows.empty())
					categories[cat0].insert(categories[cat0].end(),
					  rows.begin(), rows.end());
				return;
				}
			row_vals.push_back(val);
			}
		if (skip_atom)
			continue;
		if (keep)
			{
			map<string, string> row;
			for (uint i = 0; i < n; ++i)
				row[items[i]] = row_vals[i];
			rows.push_back(row);
			}
		}
	if (keep && !rows.empty())
		categories[cat0].insert(categories[cat0].end(),
		  rows.begin(), rows.end());
	}

static void RefsFromRows(const vector<map<string, string> > &rows,
	map<string, vector<string> > &by_entity)
	{
	by_entity.clear();
	vector<string> all_refs;
	for (size_t i = 0; i < rows.size(); ++i)
		{
		const map<string, string> &row = rows[i];
		string db, acc;
		map<string, string>::const_iterator it;
		it = row.find("db_name");
		if (it != row.end())
			db = it->second;
		it = row.find("pdbx_db_accession");
		if (it != row.end() && CifValuePresent(it->second))
			acc = it->second;
		else
			{
			it = row.find("db_accession");
			if (it != row.end() && CifValuePresent(it->second))
				acc = it->second;
			else
				{
				it = row.find("db_code");
				if (it != row.end())
					acc = it->second;
				}
			}
		string ref = FormatDbRef(db, acc);
		if (ref.empty())
			continue;
		all_refs.push_back(ref);
		string eid;
		it = row.find("entity_id");
		if (it != row.end())
			{
			eid = it->second;
			StripWhiteSpace(eid);
			}
		if (CifValuePresent(eid))
			by_entity[eid].push_back(ref);
		}
	if (by_entity.empty() && !all_refs.empty())
		by_entity[""] = all_refs;
	}

void ExtractCifMeta(const vector<string> &Lines, const string &FallbackLabel,
	string &Entry, string &Title,
	map<string, string> &MolByEntity,
	map<string, vector<string> > &RefsByEntity)
	{
	Entry = FallbackLabel;
	Title.clear();
	MolByEntity.clear();
	RefsByEntity.clear();

	map<string, vector<map<string, string> > > categories;
	map<string, map<string, string> > kv_rows;
	string data_name;

	CifLineTokStream ts(Lines);
	string tok;
	while (ts.Next(tok))
		{
		if (StartsWith(tok, "data_"))
			{
			if (data_name.empty())
				data_name = tok.substr(5);
			continue;
			}
		if (tok == "global_" || tok == "stop_" || StartsWith(tok, "save_"))
			continue;
		if (tok == "loop_")
			{
			ParseCifLoop(ts, categories);
			continue;
			}
		if (IsCifTag(tok))
			{
			string val;
			if (!ts.Next(val))
				break;
			string cat, item;
			SplitCifTag(tok, cat, item);
			if (cat == "atom_site")
				continue;
			if (KeepCifCategory(cat) && !item.empty())
				kv_rows[cat][item] = val;
			continue;
			}
		}

	for (map<string, map<string, string> >::const_iterator it = kv_rows.begin();
	  it != kv_rows.end(); ++it)
		{
		if (categories.find(it->first) == categories.end())
			categories[it->first].push_back(it->second);
		}

	Entry.clear();
	{
	map<string, vector<map<string, string> > >::const_iterator it =
		categories.find("entry");
	if (it != categories.end())
		{
		for (size_t i = 0; i < it->second.size(); ++i)
			{
			map<string, string>::const_iterator id =
				it->second[i].find("id");
			if (id != it->second[i].end() && CifValuePresent(id->second))
				{
				Entry = id->second;
				StripWhiteSpace(Entry);
				break;
				}
			}
		}
	}
	if (Entry.empty() && CifValuePresent(data_name))
		Entry = data_name;
	if (Entry.empty())
		Entry = FallbackLabel;

	Title.clear();
	{
	map<string, vector<map<string, string> > >::const_iterator it =
		categories.find("struct");
	if (it != categories.end())
		{
		for (size_t i = 0; i < it->second.size(); ++i)
			{
			map<string, string>::const_iterator tit =
				it->second[i].find("title");
			if (tit != it->second[i].end())
				{
				Title = CleanStructTextCopy(tit->second);
				if (!Title.empty())
					break;
				}
			}
		}
	}

	{
	map<string, vector<map<string, string> > >::const_iterator it =
		categories.find("entity");
	if (it != categories.end())
		{
		for (size_t i = 0; i < it->second.size(); ++i)
			{
			const map<string, string> &row = it->second[i];
			string eid;
			map<string, string>::const_iterator id = row.find("id");
			if (id != row.end())
				{
				eid = id->second;
				StripWhiteSpace(eid);
				}
			string desc;
			map<string, string>::const_iterator d = row.find("pdbx_description");
			if (d != row.end())
				desc = CleanStructTextCopy(d->second);
			if (CifValuePresent(eid) && !desc.empty())
				MolByEntity[eid] = desc;
			}
		}
	}

	{
	map<string, vector<map<string, string> > >::const_iterator it =
		categories.find("struct_ref");
	if (it != categories.end())
		RefsFromRows(it->second, RefsByEntity);
	}
	if (RefsByEntity.empty())
		{
		map<string, vector<map<string, string> > >::const_iterator it =
			categories.find("ma_target_ref_db_details");
		if (it != categories.end())
			RefsFromRows(it->second, RefsByEntity);
		}
	}

/*** PDB header metadata ***/

static string PdbField(const string &line, int start1, int end1)
	{
	string padded = line;
	if ((int) padded.size() < end1)
		padded.resize(end1, ' ');
	return padded.substr(size_t(start1 - 1), size_t(end1 - start1 + 1));
	}

static void PdbRecordText(const string &line, string &out)
	{
	if (line.size() >= 11)
		{
		string mid = line.substr(6, 4);
		StripWhiteSpace(mid);
		bool mid_ok = mid.empty();
		if (!mid_ok)
			{
			mid_ok = true;
			for (size_t i = 0; i < mid.size(); ++i)
				{
				if (!isdigit((unsigned char) mid[i]))
					{
					mid_ok = false;
					break;
					}
				}
			}
		if (mid_ok)
			{
			if (line.size() >= 80)
				out = line.substr(10, 70);
			else if (line.size() > 10)
				out = line.substr(10);
			else
				out.clear();
			return;
			}
		}
	if (line.size() > 6)
		{
		out = line.substr(6);
		size_t p = 0;
		while (p < out.size() && isspace((unsigned char) out[p]))
			++p;
		out = out.substr(p);
		}
	else
		out.clear();
	}

static void ParseCompnd(const string &text, map<string, string> &mol_by_chain)
	{
	mol_by_chain.clear();
	string current_mol;
	vector<string> current_chains;

	auto flush = [&]()
		{
		if (!current_mol.empty())
			{
			for (size_t i = 0; i < current_chains.size(); ++i)
				{
				const string &ch = current_chains[i];
				if (!ch.empty() && mol_by_chain.find(ch) == mol_by_chain.end())
					mol_by_chain[ch] = current_mol;
				}
			}
		};

	vector<string> tokens;
	{
	size_t start = 0;
	for (size_t i = 0; i <= text.size(); ++i)
		{
		if (i == text.size() || text[i] == ';')
			{
			string token = text.substr(start, i - start);
			StripWhiteSpace(token);
			if (!token.empty())
				tokens.push_back(token);
			start = i + 1;
			}
		}
	}

	for (size_t ti = 0; ti < tokens.size(); ++ti)
		{
		const string &token = tokens[ti];
		size_t colon = token.find(':');
		if (colon == string::npos)
			continue;
		string key = token.substr(0, colon);
		string val = token.substr(colon + 1);
		StripWhiteSpace(key);
		StripWhiteSpace(val);
		ToUpper(key);
		if (key == "MOL_ID")
			{
			flush();
			current_mol.clear();
			current_chains.clear();
			}
		else if (key == "MOLECULE")
			current_mol = val;
		else if (key == "CHAIN")
			{
			current_chains.clear();
			size_t s = 0;
			for (size_t i = 0; i <= val.size(); ++i)
				{
				if (i == val.size() || val[i] == ',')
					{
					string c = val.substr(s, i - s);
					StripWhiteSpace(c);
					if (!c.empty())
						current_chains.push_back(c);
					s = i + 1;
					}
				}
			}
		}
	flush();
	}

static bool ParseDbrefLine(const string &line, string &chain, string &db,
	string &acc)
	{
	string rec = line.size() >= 6 ? line.substr(0, 6) : line;
	StripWhiteSpace(rec);
	if (rec == "DBREF")
		{
		chain = PdbField(line, 13, 13);
		StripWhiteSpace(chain);
		vector<string> fields;
		if (line.size() > 26)
			SplitWhite(line.substr(26), fields);
		db = fields.empty() ? "" : fields[0];
		acc = fields.size() > 1 ? fields[1] : "";
		return true;
		}
	if (rec == "DBREF2")
		{
		chain = PdbField(line, 13, 13);
		StripWhiteSpace(chain);
		vector<string> fields;
		if (line.size() > 18)
			SplitWhite(line.substr(18), fields);
		db.clear();
		acc.clear();
		if (!fields.empty() && IsUnpDb(fields[0]))
			{
			db = fields[0];
			acc = fields.size() > 1 ? fields[1] : "";
			}
		else if (!fields.empty())
			acc = fields[0];
		return true;
		}
	return false;
	}

void ExtractPdbMeta(const vector<string> &Lines, string &Entry,
	string &Title,
	map<string, string> &MolByChain,
	map<string, vector<string> > &RefsByChain)
	{
	Title.clear();
	MolByChain.clear();
	RefsByChain.clear();

	string title_parts;
	string compnd_parts;
	map<string, string> dbref1_db;
	bool entry_from_header = false;

	for (size_t i = 0; i < Lines.size(); ++i)
		{
		const string &line = Lines[i];
		string rec = line.size() >= 6 ? line.substr(0, 6) : line;
		StripWhiteSpace(rec);
		if (rec == "HEADER" && !entry_from_header)
			{
			string code = PdbField(line, 63, 66);
			StripWhiteSpace(code);
			if (code.size() == 4)
				{
				string up = code;
				ToUpper(up);
				if (up != "XXXX")
					{
					Entry = code;
					entry_from_header = true;
					}
				}
			}
		else if (rec == "TITLE")
			{
			string part;
			PdbRecordText(line, part);
			title_parts += part;
			}
		else if (rec == "COMPND")
			{
			string part;
			PdbRecordText(line, part);
			compnd_parts += part;
			}
		else if (rec == "DBREF")
			{
			string chain, db, acc;
			if (ParseDbrefLine(line, chain, db, acc))
				{
				string ref = FormatDbRef(db, acc);
				if (!ref.empty() && !chain.empty())
					RefsByChain[chain].push_back(ref);
				}
			}
		else if (rec == "DBREF1")
			{
			string chain = PdbField(line, 13, 13);
			StripWhiteSpace(chain);
			string db = PdbField(line, 27, 32);
			StripWhiteSpace(db);
			if (!chain.empty())
				dbref1_db[chain] = db;
			}
		else if (rec == "DBREF2")
			{
			string chain, db, acc;
			if (ParseDbrefLine(line, chain, db, acc))
				{
				if (!chain.empty() && !CifValuePresent(db))
					{
					map<string, string>::const_iterator it =
						dbref1_db.find(chain);
					if (it != dbref1_db.end())
						db = it->second;
					}
				string ref = FormatDbRef(db, acc);
				if (!ref.empty() && !chain.empty())
					RefsByChain[chain].push_back(ref);
				}
			}
		else if (rec == "ATOM" || rec == "HETATM" || rec == "MODEL")
			break;
		}

	Title = CleanStructTextCopy(title_parts);
	ParseCompnd(compnd_parts, MolByChain);
	}
