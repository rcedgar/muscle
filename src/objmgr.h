#ifndef objmgr_h
#define objmgr_h

#include "objtype.h"
#include "obj.h"

class Obj;

#define T(x)	class x;
#include "objtypes.h"

const char *ObjTypeToStr(ObjType Type);

class ObjMgr
	{
	friend class Obj;

private:
	static ObjMgr **m_OMs;
	static unsigned m_ThreadCount;

private:
	Obj *m_Free[OTCount];
	Obj *m_Busy[OTCount];

#if	DEBUG
	bool m_Validate;
	unsigned m_BusyCounts[OTCount];
	unsigned m_GetCallCounts[OTCount];
	unsigned m_FreeCallCounts[OTCount];
	unsigned m_AllocCallCounts[OTCount];
#endif

private:
	ObjMgr();

public:
	static void Up(Obj *pObj);
	static void Down(Obj *pObj);
	static void LogGlobalStats();

#define T(x)	\
		static x *Get##x() { return (x *) StaticGetObj(OT_##x); } \
		static x *Get##x(const char *FileName, unsigned LineNr) { return (x *) StaticGetObj(OT_##x, FileName, LineNr); }
#include "objtypes.h"

	static ObjMgr *GetObjMgr();
	static Obj *StaticGetObj(ObjType Type);
	static Obj *StaticGetObj(ObjType Type, const char *FileName, unsigned LineNr);

public:
	static void UpdateGlobalStats();
	void ThreadUpdateGlobalStats();
	Obj *ThreadGetObj(ObjType Type);
	Obj *ThreadGetObj(ObjType Type, const char *FileName, unsigned LineNr);

#define T(x)	\
		x *Get__##x() { return (x *) ThreadGetObj(OT_##x); }
#include "objtypes.h"

	void ThreadUp(Obj *pObj)
		{
#if	DEBUG
		if (m_Validate)
			Validate();
#endif

		assert(pObj->m_RefCount != 0);
		++pObj->m_RefCount;

#if	DEBUG
		if (m_Validate)
			Validate();
#endif
		}

	void ThreadDown(Obj *pObj)
		{
#if	DEBUG
		if (m_Validate)
			Validate();
#endif
#if	TRACK_OBJ_THREAD
		asserta(pObj->m_OMPThreadIndex == omp_get_thread_num());
#endif
		assert(pObj->m_RefCount > 0);
		--pObj->m_RefCount;
		if (pObj->m_RefCount == 0)
			{
			ObjType Type = pObj->m_Type;
#if	DEBUG
			assert(m_BusyCounts[Type] > 0);
			--(m_BusyCounts[Type]);
#endif
			FreeObj(pObj);
			pObj->OnZeroRefCount();
			}

#if	DEBUG
		if (m_Validate)
			Validate();
#endif
		}

#if	DEBUG
	void Validate() const;
#endif

	Obj *AllocNew(ObjType Type);
	void FreeObj(Obj *obj);
	unsigned GetFreeCount(ObjType Type) const;
	unsigned GetBusyCount(ObjType Type) const;
	unsigned GetMaxRefCount(ObjType Type) const;
	float GetTotalMem(ObjType Type) const;
#if	DEBUG
	void ValidateType(ObjType Type) const;
#endif
	};

#endif // objmgr_h
