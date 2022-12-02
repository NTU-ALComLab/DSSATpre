/**********************************************************************
 * time_check.hpp
 *
 * Stefan Kupferschmid <skupfers@informatik.uni-freiburg.de>
 *
 **********************************************************************/
#ifndef MINICRAIG_TIME_CHECK_HPP
#define MINICRAIG_TIME_CHECK_HPP

//#define MAX_SOLVING_TIME 900000000ll // Overall CPU limit: 900 seconds.
#define MAX_SOLVING_TIME 9000000ll // Overall CPU limit: 9 seconds.

class TimeOut{
public:
    TimeOut();
    ~TimeOut();
    void init();
    void check();
private:
    long long start_time_;

    // Return the current "world time".
    long long getTime(void);
};

#endif /* TIME_CHECK_HPP */
