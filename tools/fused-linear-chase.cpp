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
//    void** p = memory;
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
#include <sys/times.h>

volatile void* global; // to defeat optimizations

/* return real time in seconds since start of the process */
double
walltime ()
{
#ifndef CLK_TCK
  /* CLK_TCK specifies the number of ticks per second */
  static int CLK_TCK = 0;
  if (!CLK_TCK)
    {
      CLK_TCK = sysconf (_SC_CLK_TCK);
    }
#endif
  struct tms timebuf;
  /* times returns the number of real time ticks passed since start */
  return (double) times (&timebuf) / CLK_TCK;
}

unsigned int
log2 (unsigned int val)
{
  unsigned int count = 0;
  while (val >>= 1)
    {
      ++count;
    }
  return count;
}

/* reverse the given number of bits within val */
unsigned int
bit_reverse (unsigned int val, unsigned int bits)
{
  unsigned int result = 0;
  while (bits > 0)
    {
      result = (result << 1) | (val & 1);
      val >>= 1;
      --bits;
    }
  return result;
}

/* generate a bit-reversal permutation; see https://oeis.org/A030109 */
void
gen_bit_reversal_permutation (unsigned int* seq, unsigned int bits,
			      unsigned int count)
{
  /* generate a bit-reversal permutation for integers from 0 to (2^bits)-1 */
  unsigned int maxval = 1 << bits;
  for (unsigned int val = 0; val < maxval; ++val)
    {
      seq[val] = bit_reverse (val, bits);
    }
  /* cut sequence short if count is not a power of 2, i.e. count < 2^bits */
  unsigned int current = maxval;
  unsigned int index = 0;
  while (current > count)
    {
      while (seq[index] < count)
	++index;
      --current;
      seq[index] = seq[current];
    }
}

/* create a cyclic pointer chain where the individual locations
 are stride bytes apart */
void**
create_linear_chain (std::size_t size, std::size_t stride)
{
  std::size_t len = size / sizeof(void*);
  void** memory = new void*[len];

  /* if we have multiple runs through the same buffer
   make sure that we operate with offsets where it appears
   more likely that the associated lines are not yet in
   one of the caches;
   to achieve this we operate with bit reversal permutations,
   if runs == 8 we would get following sequence

   0 4 2 6 1 5 3 7
   */
  unsigned int runs = stride / sizeof(void*);
  unsigned int bits = log2 (runs);
  if ((1 << bits) != runs)
    ++bits;
  unsigned int* offset = new unsigned int[1 << bits];
  gen_bit_reversal_permutation (offset, bits, runs);

  /* generate the actual pointer chain */
  unsigned int run = 0;
  void** last = nullptr;
  for (unsigned int run = 0; run < runs; ++run)
    {
      char* next = (char*) memory + offset[run] * sizeof(void*);
      if (last)
	{
	  *last = (void*) next;
	}
      last = (void**) next;
      for (;;)
	{
	  char* next = (char*) last + stride;
	  if (next >= (char*) memory + size)
	    break;
	  *last = (void*) next;
	  last = (void**) next;
	}
    }
  *last = (void*) memory; /* close the cycle */
  return memory;
}

template<typename Body, typename Object>
  inline void
  fused_action (Body body, Object& object)
  {
    body (object);
  }
template<typename Body, typename Object, typename ... Objects>
  inline void
  fused_action (Body body, Object& object, Objects&... objects)
  {
    body (object);
    fused_action (body, objects...);
  }

template<typename ... Pointers>
  double
  fused_chase (std::size_t count, Pointers&... ptrs)
  {
    double t0 = walltime ();
    // chase the pointers count times
    while (count-- > 0)
      {
	fused_action ([](void**& p)
	  { p = (void**) *p;},
		      ptrs...);
      }
    double t1 = walltime ();
    // defeat the optimization that removes the chasing
    fused_action ([](void**& p)
      { global = *p;},
		  ptrs...);
    return t1 - t0;
  }

#ifndef MIN_SIZE
#define MIN_SIZE (sizeof(void*))
#endif
#ifndef MAX_SIZE
#define MAX_SIZE 120
#endif

int
main ()
{
  void **p1, **p2, **p3, **p4, **p5, **p6, **p7, **p8;
  void **m1, **m2, **m3, **m4, **m5, **m6, **m7, **m8;
  std::cout << "                                          "
      << "data access speeds in GiB/s" << std::endl;
  std::cout << "     fuse";
  for (int i = 1; i <= 8; ++i)
    std::cout << std::setw (12) << i;
  std::cout << std::endl;
  std::cout << "   stride" << std::endl;
  for (std::size_t stride = MIN_SIZE; stride <= MAX_SIZE; stride +=
      sizeof(void*))
    {
      size_t memsize = std::min ((size_t) 1 << 26,
				 stride * 1024 * sizeof(void*));
      std::size_t count = (std::size_t) 1 << 30;
      std::cout << " " << std::setw (8) << stride;

      auto print_result = [=](int fuse, double t)
	{
	  double volume = (double) sizeof(void*) * count * fuse;
	  double speed = (double) volume / t / (1<<30); /* in GiB/s */
	  std::cout << "  " << std::setw(10) <<
	  std::fixed << std::setprecision(5) << speed;
	  std::cout.flush();
	};

      fused_action ([=](void**& p)
	{
	  p = create_linear_chain(memsize, stride);
	},
		    m1, m2, m3, m4, m5, m6, m7, m8);
      p1 = m1;
      p2 = m2;
      p3 = m3;
      p4 = m4;
      p5 = m5;
      p6 = m6;
      p7 = m7;
      p8 = m8;

      print_result (1, fused_chase (count, p1));
      print_result (2, fused_chase (count, p1, p2));
      print_result (3, fused_chase (count, p1, p2, p3));
      print_result (4, fused_chase (count, p1, p2, p3, p4));
      print_result (5, fused_chase (count, p1, p2, p3, p4, p5));
      print_result (6, fused_chase (count, p1, p2, p3, p4, p5, p6));
      print_result (7, fused_chase (count, p1, p2, p3, p4, p5, p6, p7));
      print_result (8, fused_chase (count, p1, p2, p3, p4, p5, p6, p7, p8));

      fused_action ([](void**& p)
	{
	  delete[] p;
	},
		    m1, m2, m3, m4, m5, m6, m7, m8);

      std::cout << std::endl;
      std::cout.flush ();
    }
}
