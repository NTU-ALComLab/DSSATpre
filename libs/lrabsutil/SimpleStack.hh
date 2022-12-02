/**************************************************************
 *       
 *       LRABS // SimpleStack.hh
 *       
 *       Copyright (C) 2008 Florian Pigorsch
 *       
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *       
 *       Last revision:
 *         $Revision: 715 $
 *         $Author: pigorsch $
 *       
 ***************************************************************/

#ifndef LRABS_SIMPLESTACK_HH
#define LRABS_SIMPLESTACK_HH

#include <algorithm>

#include "Assert.hh"

namespace lrabs
{
    /* only works with non class types, e.g. int, pointers, ... */
    template<typename T>
    class SimpleStack
    {
    public:
        explicit SimpleStack( std::size_t initialcapacity = 1024 );
        SimpleStack( const SimpleStack& s );
        ~SimpleStack();

        SimpleStack& operator=( const SimpleStack& s );

        bool empty() const;
        const T& top() const;
        void pop();
        void push( const T& t );
        void clear();
        std::size_t size() const;
        const T& operator[]( std::size_t index ) const;

        std::size_t capacity() const;
        
    private:
        long _top;
        std::size_t _capacity;
        T* _data;
    };
}

#include "SimpleStack.icc"

#endif /* LRABS_SIMPLESTACK_HH */
