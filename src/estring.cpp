//#include "muscle.h"
//#include "estring.h"
//
//void EString::FromPathA(const string &Path)
//	{
//	Clear();
//	const uint ColCount = SIZE(Path);
//	if (ColCount == 0)
//		return;
//
//	char c0 = Path[0];
//	if (c0 == 'M' || c0 == 'D')
//		m_ES.push_back(1);
//	else if (c0 == 'I')
//		m_ES.push_back(-1);
//	else
//		asserta(false);
//
//	uint k = 0;
//	char cPrev = c0;
//	for (uint Col = 1; Col < ColCount; ++Col)
//		{
//		char c = Path[Col];
//
//		switch (c2(cPrev, c))
//			{
//		case c2('M', 'M'):
//		case c2('D', 'D'):
//		case c2('D', 'M'):
//		case c2('M', 'D'):
//			assert(SIZE(m_ES) == k + 1);
//			++(m_ES[k]);
//			break;
//
//		case c2('I', 'D'):
//		case c2('I', 'M'):
//			assert(SIZE(m_ES) == k && m_ES[k] > 0);
//			m_ES.push_back(1);
//			++k;
//			break;
//
//		case c2('M', 'I'):
//		case c2('D', 'I'):
//			m_ES.push_back(-1);
//			++k;
//			break;
//
//		case c2('I', 'I'):
//			assert(SIZE(m_ES) == k && m_ES[k] < 0);
//			--(m_ES[k]);
//			break;
//
//		default:
//			asserta(false);
//			}
//
//		cPrev = c;
//		}
//	}
//
//void EString::FromPathB(const string &Path)
//	{
//	Clear();
//	const uint ColCount = SIZE(Path);
//	if (ColCount == 0)
//		return;
//
//	char c0 = Path[0];
//	if (c0 == 'M' || c0 == 'I')
//		m_ES.push_back(1);
//	else if (c0 == 'D')
//		m_ES.push_back(-1);
//	else
//		asserta(false);
//
//	uint k = 0;
//	char cPrev = c0;
//	for (uint Col = 1; Col < ColCount; ++Col)
//		{
//		char c = Path[Col];
//
//		switch (c2(cPrev, c))
//			{
//		case c2('M', 'M'):
//		case c2('I', 'I'):
//		case c2('I', 'M'):
//		case c2('M', 'I'):
//			assert(SIZE(m_ES) == k + 1);
//			++(m_ES[k]);
//			break;
//
//		case c2('D', 'I'):
//		case c2('D', 'M'):
//			assert(SIZE(m_ES) == k && k > 0 && m_ES[k] > 0);
//			m_ES.push_back(1);
//			++k;
//			break;
//
//		case c2('M', 'D'):
//		case c2('I', 'D'):
//			m_ES.push_back(-1);
//			++k;
//			break;
//
//		case c2('D', 'D'):
//			assert(SIZE(m_ES) == k && k > 0 && m_ES[k] < 0);
//			--(m_ES[k]);
//			break;
//
//		default:
//			asserta(false);
//			}
//
//		cPrev = c;
//		}
//	}
//
//void EString::EStringsToPath(const EString &ESA,
//  const EString &ESB, string &Path)
//	{
//	Path.clear();
//	uint iA = 0;
//	uint iB = 0;
//	const vector<int> &vA = ESA.m_ES;
//	const vector<int> &vB = ESB.m_ES;
//	uint NA = SIZE(vA);
//	uint NB = SIZE(vB);
//	int nA = vA[iA++];
//	int nB = vB[iB++];
//	for (;;)
//		{
//		char cType = 0;
//		if (nA > 0)
//			{
//			if (nB > 0)
//				{
//				cType = 'M';
//				--nA;
//				--nB;
//				}
//			else if (nB < 0)
//				{
//				cType = 'D';
//				--nA;
//				++nB;
//				}
//			else
//				asserta(false);
//			}
//		else if (nA < 0)
//			{
//			if (nB > 0)
//				{
//				cType = 'I';
//				++nA;
//				--nB;
//				}
//			else
//				asserta(false);
//			}
//		else
//			asserta(false);
//
//		Path += cType;
//
//		if (nA == 0)
//			{
//			if (iA == NA)
//				{
//				asserta(iB == NB);
//				break;
//				}
//			assert(iA < SIZE(vA));
//			nA = vA[iA++];
//			}
//		if (nB == 0)
//			{
//			assert(iB < SIZE(vB));
//			nB = vB[iB++];
//			}
//		}
//	}
//
//void EString::FromMul(const EString &ES1, const EString &ES2)
//	{
//	Clear();
//	}
