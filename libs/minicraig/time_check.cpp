#include "time_check.hpp"
#include <sys/time.h>
#include <cstdlib>
#include <iostream>


TimeOut::TimeOut():start_time_(0){}

TimeOut::~TimeOut(){}
    
void TimeOut::init(){
    start_time_ = getTime();
}
    
void TimeOut::check(void){
    
    if( getTime() - start_time_ > MAX_SOLVING_TIME ){
        std::cout << "\nTimeout!!!" << std::endl;
        exit(0);
    }
}
    
// Return the current "world time".
long long TimeOut::getTime(void){
    timeval time;
    gettimeofday(&time, NULL);
    long long wtime  = time.tv_sec;
    wtime            = wtime * 1000000l;
    wtime           += time.tv_usec;
    return wtime;
}


