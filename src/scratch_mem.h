#pragma once

class scratch_mem
	{
public:
	uint8_t *mem = 0;
	size_t size = 0;
	uint8_t *ptr = 0;
	bool owner = false;
#if DEBUG
	size_t used = 0;
#endif

public:
	scratch_mem(size_t n)
		{
		asserta(n > 0);
		mem = myalloc64(uint8_t, n);
		size = n;
		ptr = mem;
		owner = true;
#if DEBUG
		memset(mem, 0xff, n);
		used = 0;
#endif
		}

	scratch_mem(uint8_t *buffer, size_t n)
		{
		asserta(n > 0);
		mem = buffer;
		size = n;
		ptr = buffer;
		owner = false;
#if DEBUG
		used = 0;
#endif
		}

	~scratch_mem()
		{
		if (owner)
			myfree(mem);
		}

	void reset()
		{
		ptr = mem;
#if DEBUG
		used = 0;
#endif
		}

	template<class t> t *get(uint n)
		{
		uint8_t *tmp_ptr = ptr;
		ptr += n*sizeof(t);
#if DEBUG
		used += n*sizeof(t);
#endif
		asserta(size_t(ptr - mem) <= size);
		return (t *) tmp_ptr;
		}
	};
