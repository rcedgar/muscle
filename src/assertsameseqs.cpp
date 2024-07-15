#include "muscle.h"

static uint g_AssertOkCount = 0;

void _AssertSeqsEq(const char *FileName, uint LineNr,
  const MultiSequence &MSA1, const MultiSequence &MSA2)
	{
	const uint SeqCount1 = MSA1.GetSeqCount();
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount1; ++SeqIndex1)
		{
		const Sequence *Seq1 = MSA1.GetSequence((int) SeqIndex1);
		const string &Label = Seq1->m_Label;
		uint SeqIndex2 = MSA2.GetSeqIndex(Label);
		const Sequence *Seq2 = MSA2.GetSequence((int) SeqIndex2);

		//uint GSI1 = Seq1->GetGSI();
		//uint GSI2 = Seq2->GetGSI();

		Sequence *uSeq1 = Seq1->CopyDeleteGaps();
		Sequence *uSeq2 = Seq2->CopyDeleteGaps();
		int Length1 = uSeq1->GetLength();
		int Length2 = uSeq2->GetLength();

		const vector<char> &v1 = uSeq1->m_CharVec;
		const vector<char> &v2 = uSeq2->m_CharVec;
		if (v1 != v2)
			{
			Log("\n");
			Log("AssertSeqsEq >%s\n", Label.c_str());
			Log("Seq1[%d]  ", Length1);
			for (int i = 1; i < Length1; ++i)
				Log("%c", v1[i]);
			Log("\n");
			Log("Seq2[%d]  ", Length2);
			for (int i = 1; i < Length2; ++i)
				Log("%c", v2[i]);
			Log("\n");
			Die("AssertSeqsEq %s:%u", FileName, LineNr);
			}

		DeleteSequence(uSeq1);
		DeleteSequence(uSeq2);
		}
	}

void _AssertSeqsEqInput(const char *File, uint Line, const MultiSequence &MS)
	{
	const MultiSequence &GlobalMS = GetGlobalInputMS();
	const uint GN = GetGlobalMSSeqCount();

	const uint SeqCount = MS.GetSeqCount();

	set<uint> GSIs;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence *Seq = MS.GetSequence(i);
		const Sequence *InputSeq = &GetGlobalInputSeqByLabel(Seq->m_Label);
		const string &Label = string(MS.GetLabel(i));
		const string &GlobalLabel = InputSeq->m_Label;
		uint GSI = GetGSIByLabel(Label);
		if (GlobalLabel != Label)
			{
			Die("%s:%u AssertSeqsEqInput Seq(%u) GSI %u label '%s' != '%s'",
			  File, Line, i, GSI, Label.c_str(), GlobalLabel.c_str());
			}

		GSIs.insert(GSI);

		const Sequence *UngappedInputSeq = InputSeq->CopyDeleteGaps();
		const uint L = UngappedInputSeq->GetLength();
		const Sequence *UngappedSeq = Seq->CopyDeleteGaps();
		const uint MSL = UngappedSeq->GetLength();
		if (L != MSL)
			Die("%s:%u AssertSeqsEqInput Seq(%u) GSI=%u L=%u, MSL=%u, label=%s",
			  File, Line, i, GSI, L, MSL, Label.c_str());

		for (uint Pos = 0; Pos < L; ++Pos)
			{
			char InputChar = UngappedInputSeq->GetChar(Pos);
			char Char = UngappedSeq->GetChar(Pos);
			if (toupper(InputChar) != toupper(Char))
				Die("%s:%u AssertSeqsEqInput Seq(%u) GSI=%u Pos[%u]=%c,%c label=%s",
				  File, Line, i, GSI, Pos, Char, InputChar, Label.c_str());
			}

		DeleteSequence(UngappedInputSeq);
		DeleteSequence(UngappedSeq);
		}
	}

void _AssertSameSeqsVec(const char *File, uint Line, 
  const MultiSequence &MS, vector<const MultiSequence *> &v)
	{
	MultiSequence *CombinedMS = new MultiSequence;
	const uint N = SIZE(v);
	for (uint i = 0; i < N; ++i)
		{
		const MultiSequence *MS = v[i];
		const uint n = MS->GetSeqCount();
		for (uint j = 0; j < n; ++j)
			{
			const Sequence *Seq = MS->GetSequence(j);
			CombinedMS->AddSequence(Seq, false);
			}
		}
	_AssertSameSeqs(File, Line, MS, *CombinedMS);
	++g_AssertOkCount;
	delete CombinedMS;
	}

void _AssertSameSeqsVec(const char *File, uint Line, 
  const MultiSequence &MS, vector<MultiSequence *> &v)
	{
	MultiSequence *CombinedMS = new MultiSequence;
	const uint N = SIZE(v);
	for (uint i = 0; i < N; ++i)
		{
		const MultiSequence *MS = v[i];
		const uint n = MS->GetSeqCount();
		for (uint j = 0; j < n; ++j)
			{
			const Sequence *Seq = MS->GetSequence(j);
			CombinedMS->AddSequence(Seq, false);
			}
		}
	_AssertSameSeqs(File, Line, MS, *CombinedMS);
	++g_AssertOkCount;
	delete CombinedMS;
	}

void _AssertSameSeqsJoin(const char *File, uint Line, 
  const MultiSequence &MS1, const MultiSequence &MS2, const MultiSequence &MS12)
	{
	vector<const MultiSequence *> v;
	v.push_back(&MS1);
	v.push_back(&MS2);
	_AssertSameSeqsVec(File, Line, MS12, v);
	}

uint GetAssertSameSeqsOkCount()
	{
	return g_AssertOkCount;
	}

void _AssertSameLabels(const char *File, uint Line, const MultiSequence &MS)
	{
	const MultiSequence &GlobalMS = GetGlobalInputMS();
	const uint GN = GetGlobalMSSeqCount();

	const uint SeqCount = MS.GetSeqCount();

	set<uint> GSIs;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence *Seq = MS.GetSequence(i);
		uint GSI = GetGSIByLabel(Seq->m_Label);
		if (GSI >= GN)
			Die("%s:%u AssertSameLabels GSI1=%u > GN=%u",
			  File, Line, GSI, GN);
		if (GSIs.find(GSI) != GSIs.end())
			Die("%s:%u AssertSameLabels dupe GSI=%u",
			  File, Line, GSI);

		const string &Label = string(MS.GetLabel(i));
		const string &GlobalLabel = string(GlobalMS.GetLabel(GSI));
		if (GlobalLabel != Label)
			Die("%s:%u AssertSameLabels Seq(%u) GSI %u label '%s' != '%s'",
			  File, Line, i, GSI, Label.c_str(), GlobalLabel.c_str());

		GSIs.insert(GSI);
		}
	}

void _AssertSameSeqs(const char *File, uint Line, 
  const MultiSequence &MS1, const MultiSequence &MS2)
	{
	const MultiSequence &GlobalMS = GetGlobalInputMS();
	const uint GN = GetGlobalMSSeqCount();

	const uint SeqCount = MS1.GetSeqCount();
	const uint SeqCount2 = MS2.GetSeqCount();
	if (SeqCount2 != SeqCount)
		Die("%s:%u AssertSameSeqs N1=%u, N22=%u",
		  File, Line, SeqCount, SeqCount2);

	set<uint> GSIs1;
	set<uint> GSIs2;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence *Seq1 = MS1.GetSequence(i);
		const Sequence *Seq2 = MS2.GetSequence(i);
		uint GSI1 = GetGSIByLabel(Seq1->m_Label);
		uint GSI2 = GetGSIByLabel(Seq2->m_Label);
		if (GSI1 >= GN)
			Die("%s:%u AssertSameSeqs GSI1=%u > GN=%u",
			  File, Line, GSI1, GN);
		if (GSI2 >= GN)
			Die("%s:%u AssertSameSeqs GSI2=%u > GN=%u",
			  File, Line, GSI2, GN);
		if (GSIs1.find(GSI1) != GSIs1.end())
			Die("%s:%u AssertSameSeqs dupe GSI1=%u",
			  File, Line, GSI1);
		if (GSIs2.find(GSI2) != GSIs2.end())
			Die("%s:%u AssertSameSeqs dupe GSI2=%u",
			  File, Line, GSI2, GN);

		const string &Label1 = string(MS1.GetLabel(i));
		const string &Label2 = string(MS2.GetLabel(i));

		const string &GlobalLabel1 = string(GlobalMS.GetLabel(GSI1));
		const string &GlobalLabel2 = string(GlobalMS.GetLabel(GSI2));

		if (GlobalLabel1 != Label1)
			Die("%s:%u AssertSameSeqs Seq1(%u) GI %u label '%s' != '%s'",
			  File, Line, i, GSI1, Label1.c_str(), GlobalLabel1.c_str());

		if (GlobalLabel2 != Label2)
			Die("%s:%u AssertSameSeqs Seq2(%u) GI %u label '%s' != '%s'",
			  File, Line, i, GSI2, Label2.c_str(), GlobalLabel2.c_str());

		GSIs1.insert(GSI1);
		GSIs2.insert(GSI2);
		}

	for (set<uint>::const_iterator p = GSIs1.begin();
	  p != GSIs1.end(); ++p)
		{
		uint GSI1 = *p;
		if (GSIs2.find(GSI1) == GSIs2.end())
			Die("%s:%u AssertSameSeqs GSI1=%u missing in MS2",
			  File, Line, GSI1);
		}

	for (set<uint>::const_iterator p = GSIs2.begin();
	  p != GSIs2.end(); ++p)
		{
		uint GSI2 = *p;
		if (GSIs1.find(GSI2) == GSIs1.end())
			Die("%s:%u AssertSameSeqs GSI2=%u missing in MS1",
			  File, Line, GSI2);
		}

	++g_AssertOkCount;
	}
