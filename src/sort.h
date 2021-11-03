#ifndef sort_h
#define sort_h

#include "myutils.h"
#include "countsort.h"
#include <map>

#define StartTimer(x)	/* empty */
#define EndTimer(x)	/* empty */

inline void Range(vector<unsigned> &v, unsigned N)
	{
	v.clear();
	v.reserve(N);
	for (unsigned i = 0; i < N; ++i)
		v.push_back(i);
	}

inline void Range(unsigned *v, unsigned N)
	{
	for (unsigned i = 0; i < N; ++i)
		v[i] = i;
	}

void Range(unsigned *v, unsigned n);

template<class T, bool Desc> void QuickSortInPlaceRecurse(T *Values, int left, int right)
	{
	int i = left;
	int j = right;
	int Mid = (left + right)/2;
	T pivot = Values[Mid];

	while (i <= j)
		{
		if (Desc)
			{
			while (Values[i] > pivot)
				i++;
			while (Values[j] < pivot)
				j--;
			}
		else
			{
			while (Values[i] < pivot)
				i++;
			while (Values[j] > pivot)
				j--;
			}

		if (i <= j)
			{
			swap(Values[i], Values[j]);
			i++;
			j--;
			}
		}

	if (left < j)
		QuickSortInPlaceRecurse<T, Desc>(Values, left, j);

	if (i < right)
		QuickSortInPlaceRecurse<T, Desc>(Values, i, right);
	}

template<class T, bool Desc> void QuickSortOrderRecurse(const T *Values, int left, int right, unsigned *Order)
	{
	int i = left;
	int j = right;
	int Mid = (left + right)/2;
	T pivot = Values[Order[Mid]];

	while (i <= j)
		{
		if (Desc)
			{
			while (Values[Order[i]] > pivot)
				i++;
			while (Values[Order[j]] < pivot)
				j--;
			}
		else
			{
			while (Values[Order[i]] < pivot)
				i++;
			while (Values[Order[j]] > pivot)
				j--;
			}

		if (i <= j)
			{
			swap(Order[i], Order[j]);
			i++;
			j--;
			}
		}

	if (left < j)
		QuickSortOrderRecurse<T, Desc>(Values, left, j, Order);

	if (i < right)
		QuickSortOrderRecurse<T, Desc>(Values, i, right, Order);
	}

template<class T> void QuickSortInPlace(T *Values, unsigned N)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortInPlace);
	QuickSortInPlaceRecurse<T, false>(Values, 0, int(N-1));
	EndTimer(QuickSortInPlace);
	}

template<class T> void QuickSortInPlaceDesc(T *Values, unsigned N)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortInPlaceDesc);
	QuickSortInPlaceRecurse<T, true>(Values, 0, int(N-1));
	EndTimer(QuickSortInPlaceDesc);
	}

template<class T> void QuickSortOrder(const T *Values, unsigned N, unsigned *Order)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortOrder);
	Range(Order, N);
	QuickSortOrderRecurse<T, false>(Values, 0, int(N-1), Order);
	EndTimer(QuickSortOrder);
	}

template<class T> void QuickSortOrderDesc(const T *Values, unsigned N, unsigned *Order)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortOrderDesc);
	Range(Order, N);
	QuickSortOrderRecurse<T, true>(Values, 0, int(N-1), Order);
	EndTimer(QuickSortOrderDesc);
	}

template<class T> void QuickSortSubset(const T *Values, unsigned N, unsigned *Subset)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortSubset);
	QuickSortOrderRecurse<T, false>(Values, 0, int(N-1), Subset);
	EndTimer(QuickSortSubset);
	}

template<class T> void QuickSortSubsetDesc(const T *Values, unsigned N, unsigned *Subset)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortSubsetDesc);
	QuickSortOrderRecurse<T, true>(Values, 0, int(N-1), Subset);
	EndTimer(QuickSortSubsetDesc);
	}

template<class t> float GetCountFromMapFloat(map<t, float> &Map, const t &Key, bool Fail = true)
	{
	typename map<t, float>::const_iterator p = Map.find(Key);
	if (p == Map.end())
		{
		if (Fail)
			Die("GetCountFromMapFloat(), key not found");
		return 0.0f;
		}
	return p->second;
	}

template<class t> unsigned GetCountFromMap(map<t, unsigned> &Map, const t &Key, bool Fail = true)
	{
	typename map<t, unsigned>::const_iterator p = Map.find(Key);
	if (p == Map.end())
		{
		if (Fail)
			Die("GetCountFromMap(), key not found");
		return 0;
		}
	return p->second;
	}

template<class t> void IncCountMapFloat(map<t, float> &Map, const t &Key, float n)
	{
	if (Map.find(Key) == Map.end())
		Map[Key] = n;
	else
		Map[Key] += n;
	}

template<class t> void IncCountMap(map<t, unsigned> &Map, const t &Key, unsigned n = 1)
	{
	if (Map.find(Key) == Map.end())
		Map[Key] = n;
	else
		Map[Key] += n;
	}

template<class t> void VecToCountMap(const vector<t> &Values, map<t, unsigned> &Map)
	{
	Map.clear();
	const unsigned N = SIZE(Values);
	for (unsigned i = 0; i < N; ++i)
		{
		t Value = Values[i];
		IncCountMap(Map, Value, 1);
		}
	}

inline void CountMapToVecs(const map<string, unsigned> &Map,
  vector<string> &Keys, vector<unsigned> &Counts)
	{
	Keys.clear();
	Counts.clear();
	vector<string> Keys1;
	vector<unsigned> Counts1;
	for (map<string, unsigned>::const_iterator p = Map.begin(); p != Map.end(); ++p)
		{
		Keys1.push_back(p->first);
		Counts1.push_back(p->second);
		}
	const unsigned N = SIZE(Keys1);
	if (N == 0)
		return;
	unsigned *Order = myalloc(unsigned, N);
	QuickSortOrderDesc(Counts1.data(), N, Order);
	for (unsigned k = 0; k < N; ++k)
		{
		unsigned i = Order[k];
		Keys.push_back(Keys1[i]);
		Counts.push_back(Counts1[i]);
		}
	myfree(Order);
	}

//template<class t> void CountMapToVecsTpl(const map<t, unsigned> &Map,
//  vector<t> &Keys, vector<unsigned> &Counts)
//	{
//	Keys.clear();
//	Counts.clear();
//	vector<t> Keys1;
//	vector<unsigned> Counts1;
//	for (map<t, unsigned>::const_iterator p = Map.begin(); p != Map.end(); ++p)
//		{
//		Keys1.push_back(p->first);
//		Counts1.push_back(p->second);
//		}
//	const unsigned N = SIZE(Keys1);
//	if (N == 0)
//		return;
//	unsigned *Order = myalloc(unsigned, N);
//	QuickSortOrderDesc(Counts1.data(), N, Order);
//	for (unsigned k = 0; k < N; ++k)
//		{
//		unsigned i = Order[k];
//		Keys.push_back(Keys1[i]);
//		Counts.push_back(Counts1[i]);
//		}
//	myfree(Order);
//	}

template<class T, bool Desc> void QuickSortIndexesRecurse(const T *Values, int left, int right, unsigned *Indexes)
	{
	int i = left;
	int j = right;
	int Mid = (left + right)/2;
	T pivot = Values[Indexes[Mid]];

	while (i <= j)
		{
		if (Desc)
			{
			while (Values[Indexes[i]] > pivot)
				i++;
			while (Values[Indexes[j]] < pivot)
				j--;
			}
		else
			{
			while (Values[Indexes[i]] < pivot)
				i++;
			while (Values[Indexes[j]] > pivot)
				j--;
			}

		if (i <= j)
			{
			swap(Indexes[i], Indexes[j]);
			i++;
			j--;
			}
		}

	if (left < j)
		QuickSortIndexesRecurse<T, Desc>(Values, left, j, Indexes);

	if (i < right)
		QuickSortIndexesRecurse<T, Desc>(Values, i, right, Indexes);
	}

template<class T> void QuickSortIndexesInPlaceDesc(const T *Values, unsigned N, unsigned *Indexes)
	{
	if (N == 0)
		return;

	QuickSortIndexesRecurse<T, true>(Values, 0, int(N-1), Indexes);
	}

#endif // sort_h
