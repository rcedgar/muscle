#ifndef gobuff_h
#define gobuff_h

#include "myutils.h"

template<class T, unsigned SizeInc = 32, bool CopyOnGrow = false,
  bool ZeroOnGrow = false> class GoBuff
	{
public:
	unsigned MaxSize;
	unsigned Size;
	T *Data;

public:
	GoBuff()
		{
		MaxSize = 0;
		Size = 0;
		Data = 0;
		}

	~GoBuff() { Free(); }

	void Free()
		{
		myfree(Data);
		Size = 0;
		Data = 0;
		}

	void Alloc(unsigned n)
		{
		if (n <= MaxSize)
			return;

		unsigned NewMaxSize = n + SizeInc;
		T *NewBuffer = myalloc(T, NewMaxSize);
		if (Size > 0)
			{
			if (CopyOnGrow)
				memcpy(NewBuffer, Data, Size*sizeof(T));
			myfree(Data);
			}
		if (ZeroOnGrow)
			memset(NewBuffer, 0, NewMaxSize*sizeof(T));
		Data = NewBuffer;
		MaxSize = NewMaxSize;
		}

	unsigned GetMemUseBytes() const
		{
		return (MaxSize*sizeof(T));
		}
	};

const unsigned GROW64K = 0x10000;

#endif // gobuff_h
