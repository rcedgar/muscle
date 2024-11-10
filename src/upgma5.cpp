#include "muscle.h"
#include "upgma5.h"
#include "textfile.h"
#include "tree.h"

// UPGMA clustering in O(N^2) time and space.

#define	TRACE	0

static inline float AVG(float x, float y)
	{
	return (x + y)/2;
	}

void UPGMA5::LogMe() const
	{
	Log("Dist matrix\n");
	Log("     ");
	for (uint i = 0; i < m_LeafCount; ++i)
		{
		if (UINT_MAX == m_NodeIndex[i])
			continue;
		Log("  %5u", m_NodeIndex[i]);
		}
	Log("\n");

	for (uint i = 0; i < m_LeafCount; ++i)
		{
		if (UINT_MAX == m_NodeIndex[i])
			continue;
		Log("%5u  ", m_NodeIndex[i]);
		for (uint j = 0; j < m_LeafCount; ++j)
			{
			if (UINT_MAX == m_NodeIndex[j])
				continue;
			if (i == j)
				Log("       ");
			else
				{
				uint v = TriangleSubscript(i, j);
				Log("%5.2g  ", m_Dist[v]);
				}
			}
		Log("  %s", m_Labels[i].c_str());
		Log("\n");
		}

	Log("\n");
	Log("    i   Node   NrNb      Dist\n");
	Log("-----  -----  -----  --------\n");
	for (uint i = 0; i < m_LeafCount; ++i)
		{
		if (UINT_MAX == m_NodeIndex[i])
			continue;
		Log("%5u  %5u  %5u  %8.3f\n",
		  i,
		  m_NodeIndex[i],
		  m_NearestNeighbor[i],
		  m_MinDist[i]);
		}

	Log("\n");
	Log(" Node      L      R  Height  LLength  RLength\n");
	Log("-----  -----  -----  ------  -------  -------\n");
	for (uint i = 0; i <= m_InternalNodeIndex; ++i)
		Log("%5u  %5u  %5u  %6.2g  %6.2g  %6.2g\n",
		  i,
		  m_Left[i],
		  m_Right[i],
		  m_Height[i],
		  m_LeftLength[i],
		  m_RightLength[i]);
	}

void UPGMA5::Run(const string &sLinkage, Tree &tree)
	{
	LINKAGE Link = LINKAGE_Undefined;
	if      (sLinkage == "min") Link = LINKAGE_Min;
	else if (sLinkage == "max") Link = LINKAGE_Max;
	else if (sLinkage == "avg") Link = LINKAGE_Avg;
	else if (sLinkage == "biased") Link = LINKAGE_Biased;
	else
		Die("UPGMA5::Run(Linkage=%s)", sLinkage.c_str());
	Run(Link, tree);
	}

void UPGMA5::Run(LINKAGE Linkage, Tree &tree)
	{
	m_LeafCount = SIZE(m_Labels);
	asserta(SIZE(m_DistMx) == m_LeafCount);
	for (uint i = 0; i < m_LeafCount; ++i)
		asserta(SIZE(m_DistMx[i]) == m_LeafCount);

	m_TriangleSize = (m_LeafCount*(m_LeafCount - 1))/2;
	m_InternalNodeCount = m_LeafCount - 1;

	m_Dist = myalloc(float, m_TriangleSize);

	m_NodeIndex = myalloc(uint, m_LeafCount);
	m_NearestNeighbor = myalloc(uint, m_LeafCount);
	m_MinDist = myalloc(float, m_LeafCount);
	uint *Ids = myalloc(uint, m_LeafCount);
	char **Names = myalloc(char *, m_LeafCount);

	m_Left = myalloc(uint, m_InternalNodeCount);
	m_Right = myalloc(uint, m_InternalNodeCount);
	m_Height = myalloc(float, m_InternalNodeCount);
	m_LeftLength = myalloc(float, m_InternalNodeCount);
	m_RightLength = myalloc(float, m_InternalNodeCount);

	for (uint i = 0; i < m_LeafCount; ++i)
		{
		m_MinDist[i] = FLT_MAX;
		m_NodeIndex[i] = i;
		m_NearestNeighbor[i] = UINT_MAX;
		Ids[i] = i;
		Names[i] = mystrsave(m_Labels[i].c_str());
		}

	for (uint i = 0; i < m_InternalNodeCount; ++i)
		{
		m_Left[i] = UINT_MAX;
		m_Right[i] = UINT_MAX;
		m_LeftLength[i] = FLT_MAX;
		m_RightLength[i] = FLT_MAX;
		m_Height[i] = FLT_MAX;
		}

// Compute initial NxN triangular distance matrix.
// Store minimum distance for each full (not triangular) row.
// Loop from 1, not 0, because "row" is 0, 1 ... i-1,
// so nothing to do when i=0.
	for (uint i = 1; i < m_LeafCount; ++i)
		{
		uint Base = TriangleSubscript(i, 0);
		//float *Row = m_Dist + Base;
		for (uint j = 0; j < i; ++j)
			{
			float d = m_DistMx[i][j];
			if (d < 0)
				{
				d = 0;
				m_DistMx[i][j] = 0;
				m_DistMx[j][i] = 0;
				}
			m_Dist[Base++] = d;
			if (d < m_MinDist[i])
				{
				m_MinDist[i] = d;
				m_NearestNeighbor[i] = j;
				}
			if (d < m_MinDist[j])
				{
				m_MinDist[j] = d;
				m_NearestNeighbor[j] = i;
				}
			}
		asserta(Base <= m_TriangleSize);
		}

#if	TRACE
	Log("Initial state:\n");
	LogMe();
#endif

	const uint JoinCount = m_LeafCount - 1;
	for (m_InternalNodeIndex = 0; m_InternalNodeIndex < JoinCount;
	  ++m_InternalNodeIndex)
		{
		ProgressStep(m_InternalNodeIndex, JoinCount, "UPGMA5");
#if	TRACE
		Log("\n");
		Log("Internal node index %5u\n", m_InternalNodeIndex);
		Log("-------------------------\n");
#endif

	// Find nearest neighbors
		uint Lmin = UINT_MAX;
		uint Rmin = UINT_MAX;
		float dtMinDist = FLT_MAX;
		for (uint j = 0; j < m_LeafCount; ++j)
			{
			if (UINT_MAX == m_NodeIndex[j])
				continue;

			float d = m_MinDist[j];
			if (d < dtMinDist)
				{
				dtMinDist = d;
				Lmin = j;
				Rmin = m_NearestNeighbor[j];
				assert(UINT_MAX != Rmin);
				assert(UINT_MAX != m_NodeIndex[Rmin]);
				}
			}

		assert(Lmin != UINT_MAX);
		assert(Rmin != UINT_MAX);
		assert(dtMinDist != FLT_MAX);

#if	TRACE
		Log("Nearest neighbors Lmin %u[=%u] Rmin %u[=%u] dist %.3g\n",
		  Lmin,
		  m_NodeIndex[Lmin],
		  Rmin,
		  m_NodeIndex[Rmin],
		  dtMinDist);
#endif

	// Compute distances to new node
	// New node overwrites row currently assigned to Lmin
		float dtNewMinDist = FLT_MAX;
		uint uNewNearestNeighbor = UINT_MAX;
		for (uint j = 0; j < m_LeafCount; ++j)
			{
			if (j == Lmin || j == Rmin)
				continue;
			if (UINT_MAX == m_NodeIndex[j])
				continue;

			const uint vL = TriangleSubscript(Lmin, j);
			const uint vR = TriangleSubscript(Rmin, j);
			const float dL = m_Dist[vL];
			const float dR = m_Dist[vR];
			float dtNewDist = 0;

			switch (Linkage)
				{
			case LINKAGE_Avg:
				dtNewDist = AVG(dL, dR);
				break;

			case LINKAGE_Min:
				dtNewDist = min(dL, dR);
				break;

			case LINKAGE_Max:
				dtNewDist = max(dL, dR);
				break;

			case LINKAGE_Biased:
				dtNewDist = 0.1f*AVG(dL, dR) + (1 - 0.1f)*min(dL, dR);
				break;

			default:
				Die("UPGMA5: Invalid LINKAGE_%u", Linkage);
				}

		// Nasty special case.
		// If nearest neighbor of j is Lmin or Rmin, then make the new
		// node (which overwrites the row currently occupied by Lmin)
		// the nearest neighbor. This situation can occur when there are
		// equal distances in the matrix. If we don't make this fix,
		// the nearest neighbor pointer for j would become invalid.
		// (We don't need to test for == Lmin, because in that case
		// the net change needed is zero due to the change in row
		// numbering).
			if (m_NearestNeighbor[j] == Rmin)
				m_NearestNeighbor[j] = Lmin;

#if	TRACE
			Log("New dist to %u = (%u/%.3g + %u/%.3g)/2 = %.3g\n",
			  j, Lmin, dL, Rmin, dR, dtNewDist);
#endif
			m_Dist[vL] = dtNewDist;
			if (dtNewDist < dtNewMinDist)
				{
				dtNewMinDist = dtNewDist;
				uNewNearestNeighbor = j;
				}
			}

		assert(m_InternalNodeIndex < m_LeafCount - 1 || FLT_MAX != dtNewMinDist);
		assert(m_InternalNodeIndex < m_LeafCount - 1 || UINT_MAX != uNewNearestNeighbor);

		const uint v = TriangleSubscript(Lmin, Rmin);
		const float dLR = m_Dist[v];
		const float dHeightNew = dLR/2;
		const uint uLeft = m_NodeIndex[Lmin];
		const uint uRight = m_NodeIndex[Rmin];
		const float HeightLeft =
		  uLeft < m_LeafCount ? 0 : m_Height[uLeft - m_LeafCount];
		const float HeightRight =
		  uRight < m_LeafCount ? 0 : m_Height[uRight - m_LeafCount];

		m_Left[m_InternalNodeIndex] = uLeft;
		m_Right[m_InternalNodeIndex] = uRight;
		m_LeftLength[m_InternalNodeIndex] = dHeightNew - HeightLeft;
		m_RightLength[m_InternalNodeIndex] = dHeightNew - HeightRight;
		m_Height[m_InternalNodeIndex] = dHeightNew;

	// Row for left child overwritten by row for new node
		m_NodeIndex[Lmin] = m_LeafCount + m_InternalNodeIndex;
		m_NearestNeighbor[Lmin] = uNewNearestNeighbor;
		m_MinDist[Lmin] = dtNewMinDist;

	// Delete row for right child
		m_NodeIndex[Rmin] = UINT_MAX;

#if	TRACE
		Log("\nInternalNodeIndex=%u Lmin=%u Rmin=%u\n",
		  m_InternalNodeIndex, Lmin, Rmin);
		LogMe();
#endif
		}

	uint uRoot = m_LeafCount - 2;
	tree.Create(m_LeafCount, uRoot, m_Left, m_Right, m_LeftLength, m_RightLength,
	  Ids, Names);

#if	TRACE
	tree.LogMe();
#endif

	myfree(m_Dist);

	myfree(m_NodeIndex);
	myfree(m_NearestNeighbor);
	myfree(m_MinDist);
	myfree(m_Height);

	myfree(m_Left);
	myfree(m_Right);
	myfree(m_LeftLength);
	myfree(m_RightLength);
	
	for (uint i = 0; i < m_LeafCount; ++i)
		myfree(Names[i]);
	myfree(Names);
	myfree(Ids);

	m_Dist = 0;
	m_NodeIndex = 0;
	m_NearestNeighbor = 0;
	m_MinDist = 0;
	m_Height = 0;
	m_Left = 0;
	m_LeftLength = 0;
	m_RightLength = 0;
	Names = 0;
	Ids = 0;
	}

void UPGMA5::Clear()
	{
	m_Labels.clear();
	m_DistMx.clear();
	m_LabelToIndex.clear();
	}

void UPGMA5::Init(const vector<string> &Labels,
  const vector<vector<float> > &DistMx)
	{
	Clear();

	m_DistMx = DistMx;
	m_Labels = Labels;
	for (uint i = 0; i < SIZE(Labels); ++i)
		{
		const string &Label = Labels[i];
		if (m_LabelToIndex.find(Label) != m_LabelToIndex.end())
			Die("UPGMA5::Init(), duplicate label >%s", Label.c_str());
		m_LabelToIndex[Label] = i;
		}
	m_LeafCount = SIZE(m_Labels);
	}

void UPGMA5::AddLabel(const string &Label)
	{
	asserta(!Label.empty());
	if (m_LabelToIndex.find(Label) != m_LabelToIndex.end())
		return;
	uint Index = SIZE(m_Labels);
	m_Labels.push_back(Label);
	m_LabelToIndex[Label] = Index;
	}

uint UPGMA5::GetLabelIndex(const string &Label) const
	{
	map<string, uint>::const_iterator p = m_LabelToIndex.find(Label);
	asserta(p != m_LabelToIndex.end());
	uint Index = p->second;
	return Index;
	}

void UPGMA5::ReadDistMx(const string &FileName)
	{
	Progress("Reading dist mx...");
// Pass 1, labels
	FILE *f = OpenStdioFile(FileName);

	m_Labels.clear();
	m_LabelToIndex.clear();
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		AddLabel(Label1);
		AddLabel(Label2);
		}

	m_LeafCount = SIZE(m_Labels);
	m_DistMx.clear();
	m_DistMx.resize(m_LeafCount);
	for (uint i = 0; i < m_LeafCount; ++i)
		m_DistMx[i].resize(m_LeafCount, FLT_MAX);

// Pass 2, distances
	SetStdioFilePos(f, 0);
	uint LineNr = 0;
	while (ReadLineStdioFile(f, Line))
		{
		++LineNr;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		uint Index1 = GetLabelIndex(Label1);
		uint Index2 = GetLabelIndex(Label2);
		float Dist = (float) StrToFloat(Fields[2]);
		if (Index1 == Index2)
			Die("Line %u Index1=%u Index2=%u Label1='%s' Label2='%s'",
			  LineNr, Index1, Index2, Label1.c_str(), Label2.c_str());
		m_DistMx[Index1][Index2] = Dist;
		m_DistMx[Index2][Index1] = Dist;
		}

	CloseStdioFile(f);
	Progress(" done.\n");
	}

// Reseek distmx format.
// first line is distmx\tN
// next N lines are i\tLabel
// then distances are idx1\tidx2\tTS
void UPGMA5::ReadDistMx2(const string &FileName)
	{
	Progress("Reading dist mx (reseek format)...");
// Pass 1, labels
	FILE *f = OpenStdioFile(FileName);


	string Hdr;
	vector<string> HdrFields;
	bool Ok = ReadLineStdioFile(f, Hdr);
	asserta(Ok);
	Split(Hdr, HdrFields, '\t');
	asserta(SIZE(HdrFields) == 2);
	asserta(HdrFields[0] == "distmx");
	m_LeafCount = StrToUint(HdrFields[1]);
	asserta(m_LeafCount > 2);

	string Line;
	vector<string> Fields;
	m_Labels.clear();
	m_LabelToIndex.clear();
	for (uint Idx = 0; Idx < m_LeafCount; ++Idx)
		{
		Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		
		asserta(StrToUint(Fields[0]) == Idx);
		const string &Label = Fields[1];
		AddLabel(Label);
		}
	asserta(SIZE(m_Labels) == m_LeafCount);

	m_DistMx.clear();
	m_DistMx.resize(m_LeafCount);
	for (uint i = 0; i < m_LeafCount; ++i)
		m_DistMx[i].resize(m_LeafCount, 0);

// Pass 2, distances
	uint DistCount = 0;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		uint Index1 = StrToUint(Fields[0]);
		uint Index2 = StrToUint(Fields[1]);
		asserta(Index1 < m_LeafCount);
		asserta(Index2 < m_LeafCount);
		if (Index1 == Index2)
			continue;
		float Dist = (float) StrToFloat(Fields[2]);
		m_DistMx[Index1][Index2] = Dist;
		m_DistMx[Index2][Index1] = Dist;
		++DistCount;
		}
	ProgressLog("%u pair-wise distances\n", DistCount);
	if (DistCount < m_LeafCount)
		Die("Distance matrix too sparse");

	CloseStdioFile(f);
	Progress(" done.\n");
	}

void UPGMA5::FixEADistMx()
	{
	for (uint i = 0; i < m_LeafCount; ++i)
		{
		m_DistMx[i][i] = 0;
		for (uint j = 0; j < i; ++j)
			{
			float d = m_DistMx[i][j];
			asserta(d >= 0 && d <= 1);
			float NewDist = 1 - d;

			m_DistMx[i][j] = NewDist;
			m_DistMx[j][i] = NewDist;
			}
		}
	}

void UPGMA5::ScaleDistMx(bool InputIsSimilarity)
	{
	const float SCALE = 10.0f;
	float MinDist = m_DistMx[0][1];
	float MaxDist = m_DistMx[0][1];
	for (uint i = 0; i < m_LeafCount; ++i)
		{
		for (uint j = 0; j < i; ++j)
			{
			float d = m_DistMx[i][j];
			asserta(m_DistMx[j][i] == d);
			MinDist = min(MinDist, d);
			MaxDist = max(MaxDist, d);
			}
		}
	ProgressLog("Re-scaling, min %.4g, max %.4g\n", MinDist, MaxDist);

	float MinNewDist = FLT_MAX;
	float MaxNewDist = FLT_MAX;
	for (uint i = 0; i < m_LeafCount; ++i)
		{
		for (uint j = 0; j < i; ++j)
			{
			float d = m_DistMx[i][j];
			//float NewDist = SCALE*(MaxDist - d)/(MaxDist - MinDist);
			float NewDist = FLT_MAX;
			if (InputIsSimilarity)
				NewDist = SCALE*(MaxDist - d)/(MaxDist - MinDist);
			else
				NewDist = SCALE*(d - MinDist)/(MaxDist - MinDist);

			if (MinNewDist == FLT_MAX || NewDist < MinNewDist)
				MinNewDist = NewDist;
			if (MaxNewDist == FLT_MAX || NewDist > MaxNewDist)
				MaxNewDist = NewDist;

			m_DistMx[i][j] = NewDist;
			m_DistMx[j][i] = NewDist;
			}
		}
	ProgressLog("Scaled min dist %.3g, max %.3g. scale\n",
	  MinNewDist, MaxNewDist, SCALE);
	}

void cmd_upgma5()
	{
	const string &InputFileName = opt(upgma5);
	const string &OutputFileName = opt(output);

	UPGMA5 U;
	if (opt(reseek))
		{
		U.ReadDistMx2(InputFileName);
		U.ScaleDistMx();
		}
	else
		{
		U.ReadDistMx(InputFileName);
		if (opt(scaledist))
			U.ScaleDistMx();
		else if (opt(eadist))
			U.FixEADistMx();
		}

	LINKAGE Linkage = LINKAGE_Avg;
	string sLink = "avg";
	if (optset_linkage)
		{
		sLink = opt(linkage);
		if (sLink == "avg")
			Linkage = LINKAGE_Avg;
		else if (sLink == "min")
			Linkage = LINKAGE_Min;
		else if (sLink == "max")
			Linkage = LINKAGE_Max;
		else if (sLink == "biased")
			Linkage = LINKAGE_Biased;
		else
			Die("Invalid -linkage %s", sLink.c_str());
		}
	ProgressLog("UPGMA5(%s)\n", sLink.c_str());

	Tree t;
	U.Run(Linkage, t);
	TextFile TF(OutputFileName.c_str(), true);
	t.ToFile(TF);
	TF.Close();

	ProgressLog("All done.\n");
	}
