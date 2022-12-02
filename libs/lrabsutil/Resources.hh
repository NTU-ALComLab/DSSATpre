/**************************************************************
*       
*       LRABS Package // Resources.hh
*
*   Author:
*         Florian Pigorsch
*         University of Freiburg
*         pigorsch@informatik.uni-freiburg.de
*
*   Last revision:
*         $Revision: 689 $
*         $Author: pigorsch $
*         $Date$
*
***************************************************************/

#ifndef LRABS_RESOURCES_HH
#define LRABS_RESOURCES_HH

#include <iostream>
#include <string>

namespace lrabs
{
    /*!
     * return the process' current cpu time in seconds.
     */
    double cpuTime();

    /*!
     * return the process' maximum memory usage in kilo bytes.
     */
    double getMaxMem();

    class Timer
    {
    public:
        Timer();

        void reset();

        double startTime() const;
        double currentTime() const;
        double elapsedTime() const;

    private:
        double theStartTime;
    };

    class PrintScopeTime
    {
    public:
        explicit PrintScopeTime( const std::string& scopeName, bool enabled = true, std::ostream& os = std::cout )
        :enabled_( enabled ),
         os_( os )
        {
            if( enabled_ )
            {
                scopeName_ = scopeName;
                startTime_ = cpuTime();
                os << "BEGIN " << scopeName_ << " @ " << startTime_ << std::endl;
            }
        }

        ~PrintScopeTime()
        {
            if( enabled_ )
            {
                double endTime = cpuTime();
                double timeRange = endTime - startTime_;

                os_ << "END   " << scopeName_ << " @ " << endTime << " " << timeRange << std::endl;
            }
        }

    private:
        bool enabled_;
        std::string scopeName_;
        double startTime_;
        std::ostream& os_;

        PrintScopeTime();
        PrintScopeTime( const PrintScopeTime& );
        PrintScopeTime& operator=( const PrintScopeTime& );
    };

    class AddScopeTime
    {
    public:
        explicit AddScopeTime( double& t )
            :timeSum_( t ),
             startTime_( cpuTime() )
            {}

        ~AddScopeTime()
            {
                double timeRange = cpuTime() - startTime_;
                if( timeRange > 0 )
                {
                    timeSum_ += timeRange;
                }
            }

    private:
        double& timeSum_;
        double startTime_;

        AddScopeTime();
        AddScopeTime( const AddScopeTime& );
        AddScopeTime& operator=( const AddScopeTime& );
    };
}

#endif /* LRABS_RESOURCES_HH */
