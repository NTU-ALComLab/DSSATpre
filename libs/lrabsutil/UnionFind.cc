/**************************************************************
 *
 *       LRABS // UnionFind.hh
 *
 *       Copyright (C) 2009 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 749 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#include "UnionFind.hh"

#include <algorithm>
#include <iostream>

#include "Assert.hh"

lrabs::UnionFind::
UnionFind( std::size_t initialCapacity )
    : _capacity( ( initialCapacity >= 1024 ) ? initialCapacity : 1024 ),
      _size( 0 )
{
    _parent = new std::size_t[_capacity];
    _rank = new std::size_t[_capacity];
    _valid = new bool[_capacity];
}

lrabs::UnionFind::
UnionFind( const lrabs::UnionFind& other )
    : _capacity( other._capacity ),
      _size( other._size )
{
    _parent = new std::size_t[_capacity];
    _rank = new std::size_t[_capacity];
    _valid = new bool[_capacity];

    std::copy( other._parent, other._parent + _size, _parent );
    std::copy( other._rank, other._rank + _size, _rank );
    std::copy( other._valid, other._valid + _size, _valid );
}

lrabs::UnionFind::
~UnionFind()
{
    delete[] _parent;
    delete[] _rank;
    delete[] _valid;
}

lrabs::UnionFind& 
lrabs::UnionFind::
operator=( const lrabs::UnionFind& other )
{
    if( &other == this ) return *this;

    if( _capacity < other._capacity )
    {
        delete[] _parent;
        delete[] _rank;
        delete[] _valid;

        _capacity = other._capacity;

        _parent = new std::size_t[_capacity];
        _rank = new std::size_t[_capacity];
        _valid = new bool[_capacity];
    }

    _size = other._size;

    std::copy( other._parent, other._parent+_size, _parent );
    std::copy( other._rank, other._rank+_size, _rank );
    std::copy( other._valid, other._valid+_size, _valid );

    return *this;
}

bool
lrabs::UnionFind::
isValid( std::size_t t ) const noexcept
{
    if( t >= _size )
    {
        return false;
    }

    return _valid[t];
}

void
lrabs::UnionFind::
makeSet( std::size_t t )
{
    assert( !isValid( t ) );

    grow( t );

    assert( t < _capacity );

    if( t >= _size )
    {
        std::fill( _valid + _size, _valid + t, false );
        _size = t + 1;
    }

    _rank[t] = 0;
    _parent[t] = t;
    _valid[t] = true;
}

std::size_t
lrabs::UnionFind::
find( std::size_t t ) const
{
    if( _parent[t] == t )
    {
        return t;
    }
    else
    {
        _parent[t] = find( _parent[t] );
        return _parent[t];
    }
}

std::size_t
lrabs::UnionFind::
unite( std::size_t t1, std::size_t t2 )
{
    const std::size_t root1 = find( t1 );
    const std::size_t root2 = find( t2 );

    if( _rank[root1] > _rank[root2] )
    {
        _parent[root2] = root1;
        return root1;
    }
    else if( _rank[root1] < _rank[root2] )
    {
        _parent[root1] = root2;
        return root2;
    }
    else if( root1 != root2 )
    {
        _parent[root2] = root1;
        ++_rank[root1];
        return root1;
    }
    else
    {
        return root1;
    }
}

void 
lrabs::UnionFind::
grow( std::size_t maxitem )
{
    if( maxitem < _capacity ) return;

    while( _capacity <= maxitem )
    {
        _capacity *= 2;
    }

    std::size_t* newParent = new std::size_t[_capacity];
    std::copy( _parent, _parent + _size, newParent );
    delete[] _parent;
    _parent = newParent;

    std::size_t* newRank = new std::size_t[_capacity];
    std::copy( _rank, _rank + _size, newRank );
    delete[] _rank;
    _rank = newRank;

    bool* newValid = new bool[_capacity];
    std::copy( _valid, _valid + _size, newValid );
    delete[] _valid;
    _valid = newValid;
}
