#pragma once

#include "flat_enum.h"
#include "flat_dist_types.h"

void log_flat_stats(const string &msg = "");

#define TRACK_ACTIVE	0
#define TRACK_SRC		0

// Global atomics shared by all threads
// Simpler, faster and smaller compared to ObjMgr
// where one object per thread.
#if TRACK_ACTIVE
extern atomic<int64_t> g_flat_creates[FE_N];
extern atomic<int64_t> g_flat_destroys[FE_N];
extern atomic<int64_t> g_flat_bytes[FE_N];
#endif

#if TRACK_SRC
extern std::list<void *> g_flat_obj_list;
extern mutex g_flat_obj_list_lock;
#endif

template<typename T, FE fe>
class flat_base
	{
private:
	atomic<int> m_refcount;

public:
	T* m_data;
	uint32_t m_size;
#if TRACK_SRC
	FE m_fe = fe;
	const char *m_srcfile = 0;
	int m_srcline = 0;
	list<void *>::iterator m_list_iter;
#endif

protected:
	flat_base() = delete;

	flat_base(uint32_t n)
		{
		m_refcount = 0;
		m_size = n;
		m_data = n == 0 ? 0 : (T*) aligned_malloc(n*sizeof(T));
#if TRACK_ACTIVE
		++g_flat_creates[fe];
		g_flat_bytes[fe] += n*sizeof(T);
		m_refcount = 1;
#endif
#if TRACK_SRC
		m_srcfile = 0;
		m_srcline = 0;
		g_flat_obj_list_lock.lock();
		m_list_iter = g_flat_obj_list.insert(
			g_flat_obj_list.begin(), (void*) this);
		g_flat_obj_list_lock.unlock();
#endif
		}

	~flat_base()
		{
		asserta(m_refcount == 0);
#if TRACK_ACTIVE
		++g_flat_destroys[fe];
		g_flat_bytes[fe] -= m_size*sizeof(T);
#endif
		if (m_data) aligned_free(m_data);
#if TRACK_SRC
		if (m_srcfile)
			{
			g_flat_obj_list_lock.lock();
			g_flat_obj_list.erase(m_list_iter);
			g_flat_obj_list_lock.unlock();
			}
#endif
		}

private:
	void release_ref()
		{
		assert(m_refcount > 0);
		--m_refcount;
		if (m_refcount == 0)
			delete this;
		}

public:
	void add_ref()
		{
		++m_refcount;
		}

	int get_refcount() const
		{
		return m_refcount;
		}

	void truncate(uint L)
		{
		if (m_size <= L) return;
		m_size = L;
		}

	void falloc(uint32_t n)
		{
		asserta(m_size == 0);
		m_size = n;
		m_data = (T*) aligned_malloc(m_size*sizeof(T));
#if TRACK_ACTIVE
		g_flat_bytes[fe] += n*sizeof(T);
#endif
		}

public:
	template<class U>
	static void release(U*& ptr)
		{
		if (ptr)
			{
			ptr->release_ref();
			ptr = 0;
			}
		}
	};

template<typename T, FE fe>
class flat_vec : public flat_base<T, fe>
	{
protected:
	flat_vec() : flat_base<T, fe>() {} // explicit do-nothing
	flat_vec(uint32_t n) : flat_base<T, fe>(n) {}
	};

template<typename T, FE fe>
class flat_mx : public flat_base<T, fe>
	{
public:
	uint32_t m_rows = 0;
	uint32_t m_cols = 0;

protected:
	flat_mx() = delete;

	flat_mx(uint32_t rows, uint32_t cols) : flat_base<T, fe>(0)
		{
		falloc2(rows, cols);
		}

public:
	void falloc2(uint32_t rows, uint32_t cols)
		{
		flat_base<T, fe>::falloc(rows*cols);
		m_rows = rows;
		m_cols = cols;
		}

	void set(uint32_t i, uint32_t j, T value)
		{
		assert(i < m_rows);
		assert(j < m_cols);
		this->m_data[i*m_cols + j] = value;
		}

	T get(uint32_t i, uint32_t j) const
		{
		assert(i < m_rows);
		assert(j < m_cols);
		return this->m_data[i*m_cols + j];
		}
	};

#include "chainxyz.h"
#include "chaindistmx.h"
#include "chainaa.h"
#include "nnvec.h"
#include "sidvec.h"

#if TRACK_SRC
#define newflat(...) newflat_src(__VA_ARGS__, __FILE__, __LINE__)
#else
#define newflat(...) newflat(__VA_ARGS__)
#endif