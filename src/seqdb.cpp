#include "myutils.h"
#include "seqdb.h"
#include "flat_helpers.h"

void SeqDB::SetLabelToIndex()
	{
	m_LabelToIndex.clear();
	const uint N = SIZE(m_Labels);
	for (uint i = 0; i < N; ++i)
		m_LabelToIndex[m_Labels[i]] = i;
	}

uint SeqDB::GetSeqIndex(const string &Label, bool FailOnError) const
	{
	map<string, uint>::const_iterator p = m_LabelToIndex.find(Label);
	if (p == m_LabelToIndex.end())
		{
		if (FailOnError)
			Die("Not found >%s", Label.c_str());
		return UINT_MAX;
		}
	return p->second;
	}

unsigned SeqDB::AddSeq(const string &Label, const string &Seq)
	{
	unsigned SeqIndex = SIZE(m_Seqs);
	unsigned L = SIZE(Seq);
	if (SeqIndex == 0)
		{
		m_ColCount = L;
		m_IsAligned = true;
		}
	else if (L != m_ColCount)
		m_IsAligned = false;
	m_Labels.push_back(Label);
	m_Seqs.push_back(Seq);
	return SeqIndex;
	}

const byte *SeqDB::GetByteSeq(unsigned SeqIndex) const
	{
	asserta(SeqIndex < SIZE(m_Seqs));
	return (const byte *) m_Seqs[SeqIndex].c_str();
	}

const string &SeqDB::GetSeq(unsigned SeqIndex) const
	{
	asserta(SeqIndex < SIZE(m_Seqs));
	return m_Seqs[SeqIndex];
	}

const string &SeqDB::GetLabel(unsigned SeqIndex) const
	{
	asserta(SeqIndex < SIZE(m_Labels));
	return m_Labels[SeqIndex];
	}

unsigned SeqDB::GetSeqLength(unsigned SeqIndex) const
	{
	asserta(SeqIndex < SIZE(m_Seqs));
	return SIZE(m_Seqs[SeqIndex]);
	}

void SeqDB::FromFasta(const string &FileName, bool AllowGaps)
	{
	Clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	string Label;
	string Seq;
	bool HaveLabel = false;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			break;
		if (Line.empty())
			continue;
		if (Line[0] == '>')
			{
			if (HaveLabel && !Seq.empty())
				AddSeq(Label, Seq);
			Label = Line.substr(1);
			trunc_label(Label);
			Seq.clear();
			HaveLabel = true;
			}
		else
			{
			if (!AllowGaps)
				{
				string s;
				for (uint i = 0; i < SIZE(Line); ++i)
					{
					char c = Line[i];
					if (c != '-' && c != '.')
						s.push_back(c);
					}
				Seq += s;
				}
			else
				Seq += Line;
			}
		}
	if (HaveLabel && !Seq.empty())
		AddSeq(Label, Seq);
	CloseStdioFile(f);
	}
