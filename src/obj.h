#ifndef obj_h
#define obj_h

#define TRACK_OBJ_THREAD	0

#include "objtype.h"

/***
Objects are created and deleted only by ObjMgr.

Objects can contain other objects.

Getting new object:
	ptr = ObjMgr::GetXxx()
	New object has ref count = 1.

Adding reference:
	Convention: one variable per counted reference.
    assert(var == 0)
	var = ptr
	OM->Up(var)
	Up() is similar to IUnknown::AddRef().

Releasing reference:
	ObjMgr::Downvar)
	var = 0
	Down() is similar to IUknown::Release().

When objects contain references to other objects:
	In Init/Create()-like function, m_xxx = OM->GetXxx().
	In OnZeroRefCount(), call ObjMgr::Downm_xxx).

There is one ObjMgr per thread, avoids need for mutexes etc.
No Clear() member, this is confusing.

Memory allocation convention (not enforced or expressed):
	Objects use grow-only memory strategy.
	D'tor frees memory.
    D'tor is the only place memory is freed (except growing).
***/

class ObjMgr;

class Obj
	{
	friend class ObjMgr;
#if	TRACK_OBJ_THREAD
public:
	int m_OMPThreadIndex;
#endif

protected:
	ObjType m_Type;
	unsigned m_RefCount;
	Obj *m_Fwd;
	Obj *m_Bwd;

protected:
	Obj(ObjType Type)
		{
		m_Type = Type;
		m_RefCount = 0;
		m_Fwd = 0;
		m_Bwd = 0;
		}

	virtual ~Obj()
		{
		}

public:
	unsigned GetRefCount() const
		{
		return m_RefCount;
		}

// Override if contains objects, must Down() them.
	virtual void OnZeroRefCount() {}
	virtual unsigned GetMemBytes() const = 0;
	};

#endif // obj_h
