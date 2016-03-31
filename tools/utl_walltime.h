/*
 * utl_walltime.h
 *
 *  Created on: 26.03.2016
 *      Author: Andreas Borchert
 */

#ifndef TOOLS_UTL_WALLTIME_H_
#define TOOLS_UTL_WALLTIME_H_

#include <sys/times.h>
#include <unistd.h>

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




#endif /* TOOLS_UTL_WALLTIME_H_ */
