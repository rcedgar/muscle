#include "rng.h"
#include "muscle.h"

void RNG::srand_opt()
  {
  uint32_t Seed = optd(randseed, 1);
  if (Seed == 0) Seed = (uint32_t) (time(0)*getpid());
  srand(Seed);
  }


uint32_t SLCG::randu32()
  {
  m_state = m_state*214013 + 2531011;
  return m_state;
  }

void SLCG::srand(uint32_t Seed)
  {
  m_state = Seed;
  for (int i = 0; i < 10; ++i)
    randu32();
  }


uint32_t MWCG::randu32()
  {
  uint64_t sum = 2111111111*(uint64_t) m_state[3] + 1492*(uint64_t) m_state[2] + 1776*(uint64_t) m_state[1] + 5115*(uint64_t) m_state[0] + m_state[4];
  m_state[3] = m_state[2];
  m_state[2] = m_state[1];
  m_state[1] = m_state[0];
  m_state[4] = (uint32_t) (sum >> 32);
  m_state[0] = (uint32_t) sum;
  return m_state[0];
  }

void MWCG::srand(uint32_t Seed)
  {
  SLCG rng = SLCG();
  rng.srand(Seed);
  for (unsigned i = 0; i < 5; i++)
    m_state[i] = rng.randu32();
  for (unsigned i = 0; i < 100; i++)
    randu32();
  }
