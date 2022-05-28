#pragma once

#include <stdint.h>
#include <time.h>

#ifdef	WIN32
#include <process.h>
#else
#include <unistd.h>
#endif


// Abstract class for Random Number Generators.
class RNG
  {
public:
  virtual uint32_t randu32() = 0;
  void srand_opt();
  virtual void srand(uint32_t Seed) = 0;
  };

// Simple Linear Congruential Generator
// Bad properties; used just to initialize the better generator.
// Numerical values used by Microsoft C, according to wikipedia:
// http://en.wikipedia.org/wiki/Linear_congruential_generator
class SLCG: public RNG
  {
private:
  uint32_t m_state = 1;

public:
  SLCG() { srand(1); };
  SLCG(uint32_t Seed) { srand(Seed); };
  uint32_t randu32();
  void srand(uint32_t Seed);
  };

// A Multiply-With-Carry random number generator, see:
// http://en.wikipedia.org/wiki/Multiply-with-carry
// The particular multipliers used here were found on
// the web where they are attributed to George Marsaglia.
class MWCG: public RNG
  {
private:
  uint32_t m_state[5];
public:
  MWCG() { srand(1); };
  MWCG(uint32_t Seed) { srand(Seed); };
  uint32_t randu32();
  void srand(uint32_t Seed);
  };
