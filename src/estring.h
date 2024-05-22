//#pragma once
//
///***
//An "estring" is an edit string that operates on a sequence.
//An estring is represented as a vector of integers.
//It is interpreted in order of increasing suffix.
//A positive value n means copy n letters.
//A negative value -n means insert n gaps.
//Consecutive entries must have opposite sign, i.e. the
//shortest possible representation must be used.
//
//A "tpair" is a traceback path for a pairwise alignment
//represented as two estrings, one for each sequence.
//***/
//
//#define c2(c,d)	(((unsigned char) (c)) << 8 | (unsigned char) (d))
//
//class EString
//	{
//public:
//	vector<int> m_ES;
//
//public:
//	void Clear()
//		{
//		m_ES.resize(0);
//		}
//
//	void FromPathA(const string &Path);
//	void FromPathB(const string &Path);
//	void FromMul(const EString &ES1, const EString &ES2);
//
//public:
//	static void EStringsToPath(const EString &ESA,
//	  const EString &ESB, string &Path);
//	};
