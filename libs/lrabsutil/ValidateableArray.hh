#ifndef LRABS_VALIDATEABLEARRAY_HH
#define LRABS_VALIDATEABLEARRAY_HH

#include "Assert.hh"

namespace lrabs
{
    template<typename T>
    class ValidateableArray
    {
    public:
        ValidateableArray();
        explicit ValidateableArray( unsigned int size );
        ValidateableArray( const ValidateableArray& a );

        ~ValidateableArray();

        ValidateableArray& operator=( const ValidateableArray& a );
        
        unsigned int size() const;
        
        void clear();
        void invalidate();
        void resize( unsigned int size );
        
        bool isValid( unsigned int index ) const;
        const T& getItem( unsigned int index ) const;
        void setItem( unsigned int index, const T& value );
        void invalidateItem( unsigned int index );
        
        bool operator==( const ValidateableArray& a ) const;
        
    private:       
        unsigned int currentTimeStamp_;
        unsigned int capacity_;
        unsigned int size_;
        T* data_;
        unsigned int* timeStamps_;
    };

    std::ostream& operator<<( std::ostream& os, const ValidateableArray& a );
}

#include "ValidateableArray.icc"

#endif /* LRABS_VALIDATEABLEARRAY_HH */
