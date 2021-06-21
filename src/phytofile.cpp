#include "muscle.h"
#include "tree.h"
#include "textfile.h"

unsigned Tree::GetAnyNonLeafNode() const
	{
	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		if (!IsLeaf(uNodeIndex))
			return uNodeIndex;
	return NULL_NEIGHBOR;
	}

void Tree::ToFile(const string &FileName) const
	{
	if (FileName.empty())
		return;
	TextFile TF(FileName.c_str(), true);
	ToFile(TF);
	TF.Close();
	}

void Tree::ToFile(TextFile &File) const
	{
	if (IsRooted())
		{
		ToFileNodeRooted(File, m_uRootNodeIndex);
		File.PutString(";\n");
		return;
		}

// Unrooted.
	unsigned uNodeIndex = GetAnyNonLeafNode();

	File.PutString("(\n");
	ToFileNodeUnrooted(File, m_uNeighbor1[uNodeIndex], uNodeIndex);
	File.PutString(",\n");
	ToFileNodeUnrooted(File, m_uNeighbor2[uNodeIndex], uNodeIndex);
	File.PutString(",\n");
	ToFileNodeUnrooted(File, m_uNeighbor3[uNodeIndex], uNodeIndex);
	File.PutString(");\n");
	}

void Tree::ToFileNodeUnrooted(TextFile &File, unsigned uNodeIndex, unsigned uParent) const
	{
	assert(!IsRooted());

	bool bGroup = !IsLeaf(uNodeIndex);
	if (bGroup)
		File.PutString("(\n");

	if (IsLeaf(uNodeIndex))
		File.PutString(GetName(uNodeIndex));
	else
		{
		ToFileNodeUnrooted(File, GetFirstNeighbor(uNodeIndex, uParent), uNodeIndex);
		File.PutString(",\n");
		ToFileNodeUnrooted(File, GetSecondNeighbor(uNodeIndex, uParent), uNodeIndex);
		}

	if (bGroup)
		File.PutString(")");

	if (HasEdgeLength(uNodeIndex, uParent))
		File.PutFormat(":%g", GetEdgeLength(uNodeIndex, uParent));
	File.PutString("\n");
	}

void Tree::ToFileNodeRooted(TextFile &File, unsigned uNodeIndex) const
	{
	assert(IsRooted());

	bool bGroup = !IsLeaf(uNodeIndex) || IsRoot(uNodeIndex);
	if (bGroup)
		File.PutString("(\n");

	if (IsLeaf(uNodeIndex))
		File.PutString(GetName(uNodeIndex));
	else
		{
		ToFileNodeRooted(File, GetLeft(uNodeIndex));
		File.PutString(",\n");
		ToFileNodeRooted(File, GetRight(uNodeIndex));
		}

	if (bGroup)
		File.PutString(")");

	if (!IsRoot(uNodeIndex))
		{
		unsigned uParent = GetParent(uNodeIndex);
		if (HasEdgeLength(uNodeIndex, uParent))
			File.PutFormat(":%g", GetEdgeLength(uNodeIndex, uParent));
		}
	File.PutString("\n");
	}
