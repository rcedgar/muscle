#include "muscle.h"
#include "msa.h"
#include "tree.h"
#include "seqvect.h"

#if	DEBUG
static SeqVect *g_ptrMuscleSeqVect = 0;
static MSA MuscleInputMSA;

void SetMuscleInputMSA(MSA &msa)
	{
	MuscleInputMSA.Copy(msa);
	}

void SetMuscleSeqVect(SeqVect &v)
	{
	g_ptrMuscleSeqVect = &v;
	}

void ValidateMuscleIdsSeqVect(const MSA &msa)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const unsigned uId = msa.GetSeqId(uSeqIndex);
		const char *ptrNameMSA = msa.GetSeqName(uSeqIndex);
		const char *ptrName = g_ptrMuscleSeqVect->GetSeqName(uId);
		if (0 != strcmp(ptrNameMSA, ptrName))
			Quit("ValidateMuscleIdsSeqVect, names don't match");
		}
	}

void ValidateMuscleIdsMSA(const MSA &msa)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const unsigned uId = msa.GetSeqId(uSeqIndex);
		const char *ptrNameMSA = msa.GetSeqName(uSeqIndex);
		const char *ptrName = MuscleInputMSA.GetSeqName(uId);
		if (0 != strcmp(ptrNameMSA, ptrName))
			{
			Log("Input MSA:\n");
			MuscleInputMSA.LogMe();
			Log("MSA being tested:\n");
			msa.LogMe();
			Log("Id=%u\n", uId);
			Log("Input name=%s\n", ptrName);
			Log("Test name=%s\n", ptrNameMSA);
			Quit("ValidateMuscleIdsMSA, names don't match");
			}
		}
	}

void ValidateMuscleIds(const MSA &msa)
	{
	if (0 != g_ptrMuscleSeqVect)
		ValidateMuscleIdsSeqVect(msa);
	else if (0 != MuscleInputMSA.GetSeqCount())
		ValidateMuscleIdsMSA(msa);
	else
		Quit("ValidateMuscleIds, ptrMuscleSeqVect=0 && 0 == MuscleInputMSA.SeqCount()");

	}

void ValidateMuscleIdsSeqVect(const Tree &tree)
	{
	const unsigned uSeqCount = g_ptrMuscleSeqVect->GetSeqCount();
	const unsigned uNodeCount = tree.GetNodeCount();
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (!tree.IsLeaf(uNodeIndex))
			continue;
		const unsigned uId = tree.GetLeafId(uNodeIndex);
		if (uId >= uSeqCount)
			{
			tree.LogMe();
			Quit("Leaf with node index %u has id=%u, there are %u seqs",
			  uNodeIndex, uId, uSeqCount);
			}
		const char *ptrNameTree = tree.GetLeafName(uNodeIndex);
		const char *ptrName = g_ptrMuscleSeqVect->GetSeqName(uId);
		if (0 != strcmp(ptrNameTree, ptrName))
			Quit("ValidateMuscleIds: names don't match");
		}
	}

void ValidateMuscleIdsMSA(const Tree &tree)
	{
	const unsigned uNodeCount = tree.GetNodeCount();
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (!tree.IsLeaf(uNodeIndex))
			continue;
		const unsigned uId = tree.GetLeafId(uNodeIndex);
		const char *ptrNameTree = tree.GetLeafName(uNodeIndex);
		const char *ptrName = MuscleInputMSA.GetSeqName(uId);
		if (0 != strcmp(ptrNameTree, ptrName))
			Quit("ValidateMuscleIds: names don't match");
		}
	}

void ValidateMuscleIds(const Tree &tree)
	{
	if (0 != g_ptrMuscleSeqVect)
		ValidateMuscleIdsSeqVect(tree);
	else if (0 != MuscleInputMSA.GetSeqCount())
		ValidateMuscleIdsMSA(tree);
	else
		Quit("ValidateMuscleIds, ptrMuscleSeqVect=0 && 0 == MuscleInputMSA.SeqCount");
	}
#endif
