#pragma once

#include "upgma5.h"
#include "xdpmem.h"

class GTBuilder;

class GTBNode
	{
public:
	GTBuilder *m_Builder = 0;
	vector<uint> m_SeqIndexes;
	vector<uint> m_SeedSeqIndexes;
	GTBNode *m_Parent = 0;
	vector<GTBNode *> m_Children;
	UPGMA5 m_UPGMA;
	Tree m_Tree;

public:
	void Run();
	double GetProtDist(uint SeqIndexi, uint SeqIndexj);
	const byte *GetByteSeq(uint SeqIndex, uint &L) const;
	const char *GetLabel(uint SeqIndex) const;
	void DoAll();
	};
