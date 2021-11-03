#include "muscle.h"

/////////////////////////////////////////////////////////////////
// MultiSequence::Project()
//
// Given a set of indices, extract all sequences from the current
// MultiSequence object whose index is included in the set.
// Then, project the multiple alignments down to the desired
// subset, and return the projection as a new MultiSequence
// object.
/////////////////////////////////////////////////////////////////
MultiSequence* MultiSequence::Project(const set<int>& indices)
	{
	vector<vector<char> *> newPtrs(indices.size());

	assert(indices.size() != 0);

	// grab old data
	//vector<vector<char>::iterator> oldPtrs(indices.size());
	//for (set<int>::const_iterator iter = indices.begin();
	//  iter != indices.end(); ++iter)
	//	oldPtrs[i++] = GetSequence(*iter)->GetDataPtr();

	int i = 0;
	vector<const char *> oldPtrs(indices.size());
	for (set<int>::const_iterator iter = indices.begin();
	  iter != indices.end(); ++iter)
		oldPtrs[i++] = GetSequence(*iter)->GetCharPtr1();

	// compute new length
	int oldLength = GetSequence(*indices.begin())->GetLength();
	int newLength = 0;
	for (i = 1; i <= oldLength; i++)
		{
		// check to see if there is a gap in every sequence of the set
		bool found = false;
		for (int j = 0; !found && j < (int)indices.size(); j++)
			found = (oldPtrs[j][i] != '-');

		// if not, then this column counts towards the sequence length
		if (found) newLength++;
		}

	// build new alignments
	for (i = 0; i < (int)indices.size(); i++)
		{
		newPtrs[i] = new vector<char>; assert(newPtrs[i]);
		newPtrs[i]->push_back('@');
		}

	// add all needed columns
	for (i = 1; i <= oldLength; i++)
		{
		// make sure column is not gapped in all sequences in the set
		bool found = false;
		for (int j = 0; !found && j < (int)indices.size(); j++)
			found = (oldPtrs[j][i] != '-');

		// if not, then add it
		if (found)
			{
			for (int j = 0; j < (int)indices.size(); j++)
				newPtrs[j]->push_back(oldPtrs[j][i]);
			}
		}

	// wrap sequences in MultiSequence object
	MultiSequence* ret = new MultiSequence();
	i = 0;
	for (set<int>::const_iterator iter = indices.begin();
	  iter != indices.end(); ++iter)
		{
		const Sequence *OldSeq = GetSequence(*iter);
		Sequence *NewSeq = NewSequence();
		asserta(NewSeq != 0);
		vector<char> *DataPtr = newPtrs[i++];

		const string &Label = OldSeq->m_Label;
		uint GSI = OldSeq->GetGSI();
		uint SMI = OldSeq->GetSMI();

		NewSeq->Create(DataPtr, Label, GSI, SMI);
		ret->AddSequence(NewSeq, true);
		}

	return ret;
	}
