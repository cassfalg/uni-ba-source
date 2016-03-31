#ifndef HPC_UTIL_WALLTIME_H
#define HPC_UTIL_WALLTIME_H 1

#include <chrono>

namespace hpc { namespace util {


  //
  // Timer
  //
  template <typename T>
  struct WallTime
  {
      void
      tic()
      {
          t0 = std::chrono::high_resolution_clock::now();
      }

      T
      toc()
      {
          using namespace std::chrono;

          elapsed = high_resolution_clock::now() - t0;
          return duration<T,seconds::period>(elapsed).count();
      }

      std::chrono::high_resolution_clock::time_point t0;
      std::chrono::high_resolution_clock::duration   elapsed;
  };

} } // namespace util, hpc

#endif // HPC_UTIL_WALLTIME_H
