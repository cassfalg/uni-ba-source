/* small utility to measure cache and memory read access times
 afb 10/2008, 10/2015
 */

// The memory buffer is organized as an array of pointers where
//    - all pointers point into the very same buffer and where
//    - beginning from any pointer all other pointers are
//      referenced directly or indirectly, and where
//    - the pointer chain is randomized
//
// Once such a memory buffer has been set up, we measure the time of
//
//    void** p = (void**) memory[0];
//    while (count-- > 0) {
//       p = (void**) *p;
//    }
//
// The "p = (void**) *p" construct causes all memory accesses to
// be serialized, i.e. the next access can only be started whenever
// the previous is finished.
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include "utl_walltime.h"

volatile void* global; // to defeat optimizations


/* create a cyclic pointer chain that covers all words
 in a memory section of the given size in a randomized order */
void**
create_random_chain (std::size_t size)
{
  std::size_t len = size / sizeof(void*);
  void** memory = new void*[len];

  // shuffle indices
  size_t* indices = new std::size_t[len];
  for (std::size_t i = 0; i < len; ++i)
    {
      indices[i] = i;
    }
  for (std::size_t i = 0; i < len - 1; ++i)
    {
      std::size_t j = i + lrand48 () % (len - i);
      if (i != j)
	{
	  std::swap (indices[i], indices[j]);
	}
    }
  // fill memory with pointer references
  for (std::size_t i = 1; i < len; ++i)
    {
      memory[indices[i - 1]] = (void*) &memory[indices[i]];
    }
  memory[indices[len - 1]] = (void*) &memory[indices[0]];
  delete[] indices;
  return memory;
}

/* follow a pointer chain the given number of times and
 return the measured time */
double
chase_pointers (void** memory, std::size_t count)
{
  double t0 = walltime ();
  // chase the pointers count times
  void** p = (void**) memory;
  while (count-- > 0)
    {
      p = (void**) *p;
    }
  double t1 = walltime ();
  global = *p;
  return t1 - t0;
}

unsigned int
log2 (std::size_t val)
{
  unsigned int count = 0;
  while (val >>= 1)
    {
      ++count;
    }
  return count;
}

#ifndef MIN_SIZE
#define MIN_SIZE 1024
#endif
#ifndef MAX_SIZE
#define MAX_SIZE 1024 * 1024 * 32
#endif
#ifndef GRANULARITY
#define GRANULARITY 1u
#endif

int
main ()
{
  double t=0;
  std::cout << "memsize[b]      [kb]     time in ns   speed in GiB/s" << std::endl;
  for (std::size_t memsize = MIN_SIZE; memsize <= MAX_SIZE;
      memsize += (1 << (std::max (GRANULARITY, log2 (memsize/8)) - GRANULARITY)))
    {
      void** memory = create_random_chain (memsize);
      std::size_t count = std::max (memsize * 16, (std::size_t) 1 << 30);
      t = chase_pointers (memory, count);
      delete[] memory;
      double ns = t * 1000000000 / count;
      std::cout << " " << std::setw (10) << memsize;
      std::cout << " " << std::setw (8) << memsize / 1024;
      std::cout << "  " << std::fixed << std::setw (10) << std::setprecision (5)  << ns;

	  double volume = (double) sizeof(void*) * count;
	  double speed = (double) volume / t / (1<<30); /* in GiB/s */
	  std::cout << "     " << std::setw(8) <<
	  std::fixed << std::setprecision(5) << speed;

      std::cout << std::endl;
      std::cout.flush ();
    }
}
