#ifndef LRABS_VECTOR_HH
#define LRABS_VECTOR_HH

#include <vector>

namespace lrabs
{
    template<typename T>
    class Vector
    {
    public:
        Vector& operator+( const T& t );
        Vector& operator<<( const T& t );

        Vector& operator+( const Vector<T>& t );
        Vector& operator<<( const Vector<T>& t );

        std::vector<T> stdvec() const;

    private:
        std::vector<T> _data;
    };
}

#include "Vector.icc"

#endif
