#ifndef DistFunc_h
#define DistFunc_h

class DistFunc
	{
public:
	const char *m_DistType = "DistFunc::m_DistType";

public:
	DistFunc();
	virtual ~DistFunc();

public:
	virtual void SetCount(unsigned uCount);
	virtual void SetDist(unsigned uIndex1, unsigned uIndex2, float dDist);

	void SetName(unsigned uIndex, const char szName[]);
	void SetId(unsigned uIndex, unsigned uId);
	const char *GetName(unsigned uIndex) const;
	unsigned GetId(unsigned uIndex) const;

	virtual float GetDist(unsigned uIndex1, unsigned uIndex2) const;
	virtual unsigned GetCount() const;

	void LogMe() const;

protected:
	unsigned VectorIndex(unsigned uIndex, unsigned uIndex2) const;
	unsigned VectorLength() const;

private:
	unsigned m_uCount;
	unsigned m_uCacheCount;
	float *m_Dists;
	char **m_Names;
	unsigned *m_Ids;
	};

#endif	// DistFunc_h
