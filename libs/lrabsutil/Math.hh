/**************************************************************
*       
*       LRABS // Math.hh
*
*       Copyright (C) 2007 Florian Pigorsch
*
*       Author:
*         Florian Pigorsch
*         University of Freiburg
*         pigorsch@informatik.uni-freiburg.de
*
*       Last revision:
*         $Revision: 156 $
*         $Author: pigorsch $
*         $Date$
*
***************************************************************/

#ifndef LRABS_MATH_HH
#define LRABS_MATH_HH

namespace lrabs
{
    template<typename T>
    T abs( const T& t );
    
    template<typename T>
    const T& max( const T& t1, const T& t2 );
    
    template<typename T>
    const T& max( const T& t1, const T& t2, const T& t3 );

    template<typename T>
    const T& min( const T& t1, const T& t2 );
    
    template<typename T>
    const T& min( const T& t1, const T& t2, const T& t3 );

    template<typename T>
    unsigned int countOnes( T v );
}

#include "Math.icc"

#endif /* LRABS_MATH_HH */
