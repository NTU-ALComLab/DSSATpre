/**************************************************************
*       
*       LRABS // Random.cc
*
*       Copyright (C) 2006 Florian Pigorsch
*
*       Author:
*         Florian Pigorsch
*         University of Freiburg
*         pigorsch@informatik.uni-freiburg.de
*
*       Last revision:
*         $Revision: 711 $
*         $Author: pigorsch $
*         $Date$
*
***************************************************************/

#include "Random.hh"

#include <climits>

lrabs::Random::
Random( unsigned int seed )
    :_seed( seed ),
     _boolPool(0),
     _boolPoolRemaining( 0 )
{}

bool
lrabs::Random::
getBool()
{
    if( _boolPoolRemaining == 0 )
    {
        _boolPool = getRawRandom();
        _boolPoolRemaining = ( sizeof( unsigned int ) * CHAR_BIT ) - 1;
    }

    bool result = ( _boolPool & 1 );
    --_boolPoolRemaining;
    _boolPool >>= 1;

    return result;
}
