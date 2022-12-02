/**************************************************************
 *
 *       LRABS // BoolVector.cc
 *
 *       Copyright (C) 2007 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 236 $
 *         $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#include "BoolVector.hh"

#include <cstring>
#include <cassert>
#include <ostream>
#include <limits>

lrabs::BoolVector::
BoolVector()
    :_size( 0 ),
     _binCount( 0 ),
     _lastBinMask( 0 ),
     _bins( 0 )
{}

lrabs::BoolVector::
BoolVector( std::size_t size, bool initial )
    :_size( 0 ),
     _binCount( 0 ),
     _lastBinMask( 0 ),
     _bins( 0 )
{
    if( size > 0 )
    {
        initialize( size, initial );
    }
}

lrabs::BoolVector::
BoolVector( const BoolVector& v )
    :_size( v._size ),
     _binCount( v._binCount ),
     _lastBinMask( v._lastBinMask ),
     _bins( 0 )
{
    if( _binCount != 0 )
    {
        _bins = new BinType[ _binCount ];

        memcpy( _bins, v._bins, _binCount * sizeof( BinType ) );
    }
}

lrabs::BoolVector::
~BoolVector()
{
    if( _bins != 0 ) 
    {
        delete[] _bins;
    }
}

lrabs::BoolVector& 
lrabs::BoolVector::
operator=( const lrabs::BoolVector& v )
{
    if( this == &v ) return *this;

    if( _binCount != v._binCount )
    {
        _binCount = v._binCount;
        delete[] _bins;

        if( v._bins != nullptr )
        {
            _bins = new BinType[ _binCount ];
        }
        else
        {
            _bins = nullptr;
        }
    }

    _size = v._size;
    _lastBinMask = v._lastBinMask;

    if( _binCount != 0 && _bins != nullptr)
    {
        memcpy( _bins, v._bins, _binCount * sizeof( BinType ) );
    }

    return *this;
}

std::size_t
lrabs::BoolVector::
countTrue() const
{
    /*
     * Use Kerninghan's method for counting 1-bits in a word. Published in 1988,
     * the C Programming Language 2nd Ed. (by Brian W. Kernighan and Dennis M. Ritchie)
     * mentions this in exercise 2-9. See also:
     * http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetKernighan
     */

    std::size_t count = 0;

    std::size_t fullBins = _binCount;
    if( _lastBinMask != 0 ) --fullBins;

    /* use the fast method for the first "fullBins" bins */
    for( std::size_t i = 0; i != fullBins; ++i )
    {
        BinType b = _bins[ i ];

        while( b != 0 )
        {
            /* clear the least significant bit set */
            b &= b - 1;
            ++count;
        }
    }

    /* count 1-bits  in the last (partial) bin */
    if( _lastBinMask != 0 )
    {
        BinType b = _bins[ _binCount - 1 ] & _lastBinMask;
        while( b != 0 )
        {
            /* clear the least significant bit set */
            b &= b - 1;
            ++count;
        }
    }

    return count;
}

long
lrabs::BoolVector::
getFirstTrue() const
{
    std::size_t fullBins = _binCount;
    if( _lastBinMask != 0 ) --fullBins;

    long index = 0;

    for( std::size_t i = 0; i != fullBins; ++i )
    {
        BinType b = _bins[ i ];

        if( b != 0 )
        {
            while( true )
            {
                if( b & 1 )
                {
                    return index;
                }
                else 
                {
                    b >>= 1;
                    ++index;
                }
            }
        }
        else
        {
            index += BinSize;
        }
    }

    if( _lastBinMask != 0 )
    {
        BinType b = _bins[ _binCount - 1 ] & _lastBinMask;
        if( b != 0 )
        {
            while( true )
            {
                if( b & 1 )
                {
                    return index;
                }
                else 
                {
                    b >>= 1;
                    ++index;
                }
            }
        }
    }

    return -1;
}



void
lrabs::BoolVector::
initialize( std::size_t size, bool initial )
{
    assert( uninitialized() );

    assert( size > 0 );

    _size = size;
    _binCount = _size / BinSize;
    _lastBinMask = 0;

    if( _size % BinSize != 0 ) 
    {
        ++_binCount;

        std::size_t additionalBits = _size - ( ( _binCount - 1 ) * BinSize );
        assert( additionalBits > 0 );

        for( std::size_t i = 0; i != additionalBits; ++i )
        {
            _lastBinMask |= (BinType)1ul << i;
        }
    }

    _bins = new BinType[ _binCount ];

    memset( _bins, ( initial ? std::numeric_limits<unsigned char>::max() : 0 ), _binCount * sizeof( BinType ) );
}

void
lrabs::BoolVector::
set( std::size_t index, bool value )
{
    const std::size_t bin = index / BinSize;
    const std::size_t bit = index % BinSize;

    if( value )
    {
        _bins[ bin ] |= (BinType)1ul << bit;
    }
    else
    {
        _bins[ bin ] &= ~( (BinType)1ul << bit );
    }
}

void
lrabs::BoolVector::
flip( std::size_t index )
{
    const std::size_t bin = index / BinSize;
    const std::size_t bit = index % BinSize;

    _bins[ bin ] ^= (BinType)1ul << bit;
}

void
lrabs::BoolVector::
set( bool value )
{
    assert( _bins != nullptr );
    assert( _binCount != 0 );
    memset( _bins, ( value ? std::numeric_limits<unsigned char>::max() : 0 ), _binCount * sizeof( BinType ) );
}

bool
lrabs::BoolVector::
allTrue() const
{
    if( _lastBinMask == 0 )
    {
        for( std::size_t i = 0; i != _binCount; ++i )
        {
            if( _bins[ i ] != ~(BinType)0ul )
            {
                return false;
            }
        }
    }
    /* last bin is not full */
    else
    {
        assert( _binCount >= 1 );
        for( std::size_t i = 0; i != _binCount - 1; ++i )
        {
            if( _bins[ i ] != ~(BinType)0ul )
            {
                return false;
            }
        }
        if( ( _bins[ _binCount - 1 ] & _lastBinMask ) != _lastBinMask ) 
        {
            return false;
        }
    }

    return true;
}

bool
lrabs::BoolVector::
allFalse() const
{
    if( _lastBinMask == 0 )
    {
        for( std::size_t i = 0; i != _binCount; ++i )
        {
            if( _bins[ i ] != (BinType)0ul )
            {
                return false;
            }
        }
    }
    /* last bin is not full */
    else
    {
        assert( _binCount >= 1 );
        for( std::size_t i = 0; i != _binCount - 1; ++i )
        {
            if( _bins[ i ] != (BinType)0ul )
            {
                return false;
            }
        }

        if( ( _bins[ _binCount - 1 ] & _lastBinMask ) != (BinType)0ul ) 
        {
            return false;
        }
    }

    return true;
}

bool
lrabs::BoolVector::
operator==( const lrabs::BoolVector& b ) const
{
    if( _size != b._size ) return false;

    if( _lastBinMask == 0 )
    {
        for( std::size_t i = 0; i != _binCount; ++i )
        {
            if( _bins[ i ] != b._bins[ i ] )
            {
                return false;
            }
        }
    }
    else
    {
        assert( _binCount >= 1 );

        for( std::size_t i = 0; i != _binCount - 1; ++i )
        {
            if( _bins[ i ] != b._bins[ i ] )
            {
                return false;
            }
        }

        if( ( _bins[ _binCount - 1 ] & _lastBinMask ) != ( b._bins[ _binCount - 1 ] & _lastBinMask ) )
        {
            return false;
        }
    }

    return true;
}

lrabs::BoolVector
lrabs::BoolVector::
operator~() const
{
    BoolVector result( *this );
    result.flip();
    return result;
}

void
lrabs::BoolVector::
flip()
{
    for( std::size_t i = 0; i != _binCount; ++i )
    {
        _bins[ i ] = ~_bins[ i ];
    }
}

lrabs::BoolVector& 
lrabs::BoolVector::
operator &=( const lrabs::BoolVector& v )
{
    assert( _size == v._size );

    for( std::size_t i = 0; i != _binCount; ++i )
    {
        _bins[ i ] &= v._bins[ i ];
    }

    return *this;
}

lrabs::BoolVector& 
lrabs::BoolVector::
operator |=( const lrabs::BoolVector& v )
{
    assert( _size == v._size );

    for( std::size_t i = 0; i != _binCount; ++i )
    {
        _bins[ i ] |= v._bins[ i ];
    }

    return *this;
}

lrabs::BoolVector& 
lrabs::BoolVector::
operator ^=( const lrabs::BoolVector& v )
{
    assert( _size == v._size );

    for( std::size_t i = 0; i != _binCount; ++i )
    {
        _bins[ i ] ^= v._bins[ i ];
    }

    return *this;
}

bool
lrabs::BoolVector::
intersectAndCheckEmpty( lrabs::BoolVector& v1, const lrabs::BoolVector& v2 )
{
    assert( v1._size == v2._size );

    bool empty = true;

    if( v1._lastBinMask == 0 )
    {
        for( std::size_t i = 0; i != v1._binCount; ++i )
        {
            v1._bins[ i ] &= v2._bins[ i ];
            if( v1._bins[ i ] != (BinType)0ul ) empty = false;
        }
    }
    else
    {
        for( unsigned int i = 0; i != v1._binCount - 1; ++i )
        {
            v1._bins[ i ] &= v2._bins[ i ];
            if( v1._bins[ i ] != (BinType)0ul ) empty = false;
        }
        v1._bins[ v1._binCount - 1 ] &= v2._bins[ v1._binCount - 1 ];
        if( ( v1._bins[ v1._binCount - 1 ] & v1._lastBinMask ) != (BinType)0 )
        {
            empty = false;
        }
    }

    return empty;
}

std::ostream& 
lrabs::operator<<( std::ostream& os, const lrabs::BoolVector& b )
{
    for( std::size_t i = 0; i != b.size(); ++i )
    {
        if( b.get( i ) ) os << "1";
        else os << "0";
    }

    return os;
}
