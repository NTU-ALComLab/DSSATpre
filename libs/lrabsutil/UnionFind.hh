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
 *         $Revision: 752 $
 *         $Author: pigorsch $
 *       
 ***************************************************************/

#ifndef LRABS_UNIONFIND_HH
#define LRABS_UNIONFIND_HH

#include <cstddef>

namespace lrabs
{
    class UnionFind
    {
    public:
        explicit UnionFind( std::size_t initialCapacity = 1024 );
        UnionFind( const UnionFind& other );
        ~UnionFind();
        UnionFind& operator=( const UnionFind& other );

        bool isValid( std::size_t t ) const noexcept;
        void makeSet( std::size_t t );
        std::size_t find( std::size_t t ) const;
        std::size_t unite( std::size_t t1, std::size_t t2 );

    private:
        void grow( std::size_t maxitem );

        std::size_t _capacity;
        std::size_t _size;
        mutable std::size_t* _parent;
        std::size_t* _rank;
        bool* _valid;
    };
}

#endif /* LRABS_UNIONFIND_HH */
