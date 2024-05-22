#pragma once

// Reference is set of pairwise alignments which
//  are not necessarily consistent, calculates Q only.
class QScorer2
	{
public:
	const MSA *m_Test;
	const MultiSequence *m_Ref;

public:
	double Run(const MSA &Test, const MultiSequence &Ref);
	double Run(const MultiSequence &Test, const MultiSequence &Ref);
	double GetQ(const string &T1, const string &T2,
	  const string &R1, const string &R2) const;
	void GetPosToCol(const string &Row,
	  vector<uint> &PosToCol) const;
	void StripGaps(const string &Row, string &Seq) const;
	};
