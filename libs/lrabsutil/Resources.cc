/**************************************************************
*       
*       LRABS Package // Resources.cc
*
*   Author:
*         Florian Pigorsch
*         University of Freiburg
*         pigorsch@informatik.uni-freiburg.de
*
*   Last revision:
*         $Revision: 716 $
*         $Author: pigorsch $
*         $Date$
*
***************************************************************/

#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
    
#include "Resources.hh"

double lrabs::cpuTime() 
{
    struct rusage u;
    getrusage( RUSAGE_SELF, &u );
    
    return static_cast<double>( u.ru_utime.tv_sec ) + 0.000001 * static_cast<double>( u.ru_utime.tv_usec );
}

double lrabs::getMaxMem() 
{
    struct rusage u;
    getrusage( RUSAGE_SELF, &u );
    
    return static_cast<double>( u.ru_maxrss ) * 1024.0 / static_cast<double>( getpagesize() );
}

lrabs::Timer::
Timer()
{
    reset();
}

void 
lrabs::Timer::
reset()
{
    theStartTime = currentTime();
}

double 
lrabs::Timer::
startTime() const
{
    return theStartTime;
}

double 
lrabs::Timer::
currentTime() const
{
    return lrabs::cpuTime();
}

double 
lrabs::Timer::
elapsedTime() const
{
    double t = currentTime() - startTime();
    if( t < 0 ) t = 0;
    return t;
}
