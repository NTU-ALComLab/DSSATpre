#ifndef _SET_ITERATORS_HPP
#define _SET_ITERATORS_HPP

#include <type_traits>

namespace lrabs {

// forward declarations
template <typename Container1, typename Container2>
class union_iterator;

template <typename Container1, typename Container2>
class intersection_iterator;

template <typename Container1, typename Container2>
class difference_iterator;

template <typename Container1, typename Container2>
void swap(union_iterator<Container1, Container2>& lhs, union_iterator<Container1, Container2>& rhs);

template <typename Container1, typename Container2>
void swap(intersection_iterator<Container1, Container2>& lhs, intersection_iterator<Container1, Container2>& rhs);

template <typename Container1, typename Container2>
void swap(difference_iterator<Container1, Container2>& lhs, difference_iterator<Container1, Container2>& rhs);


/**
 * \brief Iterator to iterate over the union of two sorted sequences without duplicates.
 *
 * The iterator is a forward iterator, providing only read access.
 * The iterator visits the elements in increasing order.
 */
template <typename Container1, typename Container2>
class union_iterator
{
public:
    /**
     * \brief Creates a union iterator from iterators for the two sequences
     *
     * \param _a iterator pointing to the first element of the first sequence
     * \param _a_end iterator pointing after the last element of the first sequence
     * \param _b iterator pointing to the first element of the second sequence
     * \param _b_end iterator pointing after the last element of the second sequence
     */
    union_iterator(typename Container1::const_iterator _a,
                   typename Container1::const_iterator _a_end,
                   typename Container2::const_iterator _b,
                   typename Container2::const_iterator _b_end):
        a(_a),
        a_end(_a_end),
        b(_b),
        b_end(_b_end)
    {}

    const typename Container1::value_type& operator*() const {
        if (b == b_end || (a != a_end && *a < *b)) return *a; else return *b;
    }

    union_iterator<Container1, Container2>& operator++() // prefix increment
    {
        if (a == a_end && b != b_end) ++b;
        else if (a != a_end && b == b_end) ++a;
        else if (*a < *b) ++a;
        else if (*b < *a) ++b;
        else { ++a; ++b; }

        return *this;
    }

    union_iterator<Container1, Container2> operator++(int) // postfix increment
    {
        union_iterator<Container1, Container2> tmp(*this);
        operator++(); // prefix-increment this instance
        return tmp;
    }

    bool operator==(const union_iterator<Container1, Container2>& other) const { return a == other.a && b == other.b; }
    bool operator!=(const union_iterator<Container1, Container2>& other) const { return a != other.a || b != other.b; }


private:
    typename Container1::const_iterator a;
    typename Container1::const_iterator a_end;
    typename Container2::const_iterator b;
    typename Container2::const_iterator b_end;

    friend void swap<Container1, Container2>(union_iterator<Container1, Container2>& lhs, union_iterator<Container1, Container2>& rhs);

    static_assert(std::is_same<typename Container1::value_type, typename Container2::value_type>::value,
                  "The elements of both containers must have the same type");

};


/**
 * \brief Iterator to iterate over the intersection of two sorted sequences without duplicates.
 *
 * The iterator is a forward iterator, providing only read access.
 * The iterator visits the elements in increasing order.
 */
template <typename Container1, typename Container2>
class intersection_iterator
{
public:
    intersection_iterator(typename Container1::const_iterator _a,
                          typename Container1::const_iterator _a_end,
                          typename Container2::const_iterator _b,
                          typename Container2::const_iterator _b_end):
        a(_a),
        a_end(_a_end),
        b(_b),
        b_end(_b_end)
    {
        while (a != a_end && b != b_end && *a != *b) {
            if (*a < *b) ++a;
            else ++b;
        }
        if (a == a_end) b = b_end;
        else if (b == b_end) a = a_end;
    }

    const typename Container1::value_type& operator*() const {
        return *a;
    }

    intersection_iterator<Container1, Container2>& operator++() // prefix increment
    {
        if (a==a_end) return *this;

        ++a;
        ++b;

        while (a != a_end && b != b_end && *a != *b) {
            if (*a < *b) ++a;
            else ++b;
        }
        if (a == a_end) b = b_end;
        else if (b == b_end) a = a_end;

        return *this;
    }

    intersection_iterator<Container1, Container2> operator++(int) // postfix increment
    {
        intersection_iterator<Container1, Container2> tmp(*this);
        operator++(); // prefix-increment this instance
        return tmp;
    }

    bool operator==(const intersection_iterator<Container1, Container2>& other) const { return a == other.a && b == other.b; }
    bool operator!=(const intersection_iterator<Container1, Container2>& other) const { return a != other.a || b != other.b; }


private:
    typename Container1::const_iterator a;
    typename Container1::const_iterator a_end;
    typename Container2::const_iterator b;
    typename Container2::const_iterator b_end;

    friend void swap<Container1, Container2>(intersection_iterator<Container1, Container2>& lhs, intersection_iterator<Container1, Container2>& rhs);

    static_assert(std::is_same<typename Container1::value_type, typename Container2::value_type>::value,
                  "The elements of both containers must have the same type");

};


/**
 * \brief Iterator to iterate over the difference set of two sorted sequences without duplicates.
 *
 * The iterator is a forward iterator, providing only read access.
 * The iterator visits the elements in increasing order.
 */
template <typename Container1, typename Container2>
class difference_iterator
{
public:
    difference_iterator(typename Container1::const_iterator _a,
                          typename Container1::const_iterator _a_end,
                          typename Container2::const_iterator _b,
                          typename Container2::const_iterator _b_end):
        a(_a),
        a_end(_a_end),
        b(_b),
        b_end(_b_end)
    {
        while (a != a_end && b != b_end) {
            if (*a == *b) { ++a; ++b; }
            else if (*a > *b) { ++b; }
            else break;
        }
        if (a == a_end) b = b_end;
    }

    const typename Container1::value_type& operator*() const {
        return *a;
    }

    difference_iterator<Container1, Container2>& operator++() // prefix increment
    {
        ++a;

        while (a != a_end && b != b_end) {
            if (*a == *b) { ++a; ++b; }
            else if (*a > *b) { ++b; }
            else break;
        }
        if (a == a_end) b = b_end;

        return *this;
    }

    difference_iterator<Container1, Container2> operator++(int) // postfix increment
    {
        difference_iterator<Container1, Container2> tmp(*this);
        operator++(); // prefix-increment this instance
        return tmp;
    }

    bool operator==(const difference_iterator<Container1, Container2>& other) const { return a == other.a && b == other.b; }
    bool operator!=(const difference_iterator<Container1, Container2>& other) const { return a != other.a || b != other.b; }


private:
    typename Container1::const_iterator a;
    typename Container1::const_iterator a_end;
    typename Container2::const_iterator b;
    typename Container2::const_iterator b_end;

    friend void swap<Container1, Container2>(difference_iterator<Container1, Container2>& lhs, difference_iterator<Container1, Container2>& rhs);

    static_assert(std::is_same<typename Container1::value_type, typename Container2::value_type>::value,
                  "The elements of both containers must have the same type");

};



template <typename Container1, typename Container2>
inline union_iterator<Container1, Container2> union_iterator_begin(const Container1& a, const Container2& b)
{
    return union_iterator<Container1, Container2>(a.cbegin(), a.cend(), b.cbegin(), b.cend());
}


template <typename Container1, typename Container2>
inline difference_iterator<Container1, Container2> difference_iterator_begin(const Container1& a, const Container2& b)
{
    return difference_iterator<Container1, Container2>(a.cbegin(), a.cend(), b.cbegin(), b.cend());
}


template <typename Container1, typename Container2>
inline intersection_iterator<Container1, Container2> intersection_iterator_begin(const Container1& a, const Container2& b)
{
    return intersection_iterator<Container1, Container2>(a.cbegin(), a.cend(), b.cbegin(), b.cend());
}


template <typename Container1, typename Container2>
inline union_iterator<Container1, Container2> union_iterator_end(const Container1& a, const Container2& b)
{
    return union_iterator<Container1, Container2>(a.cend(), a.cend(), b.cend(), b.cend());
}


template <typename Container1, typename Container2>
inline intersection_iterator<Container1, Container2> intersection_iterator_end(const Container1& a, const Container2& b)
{
    return intersection_iterator<Container1, Container2>(a.cend(), a.cend(), b.cend(), b.cend());
}


template <typename Container1, typename Container2>
inline difference_iterator<Container1, Container2> difference_iterator_end(const Container1& a, const Container2& b)
{
    return difference_iterator<Container1, Container2>(a.cend(), a.cend(), b.cend(), b.cend());
}


template <typename Container1, typename Container2>
inline void swap(union_iterator<Container1, Container2>& lhs, union_iterator<Container1, Container2>& rhs)
{
    std::swap(lhs.a, rhs.a);
    std::swap(lhs.a_end, rhs.a_end);
    std::swap(lhs.b, rhs.b);
    std::swap(lhs.b_end, rhs.b_end);
}


template <typename Container1, typename Container2>
inline void swap(intersection_iterator<Container1, Container2>& lhs, intersection_iterator<Container1, Container2>& rhs)
{
    std::swap(lhs.a, rhs.a);
    std::swap(lhs.a_end, rhs.a_end);
    std::swap(lhs.b, rhs.b);
    std::swap(lhs.b_end, rhs.b_end);
}


template <typename Container1, typename Container2>
inline void swap(difference_iterator<Container1, Container2>& lhs, difference_iterator<Container1, Container2>& rhs)
{
    std::swap(lhs.a, rhs.a);
    std::swap(lhs.a_end, rhs.a_end);
    std::swap(lhs.b, rhs.b);
    std::swap(lhs.b_end, rhs.b_end);
}

} // end namespace lrabs

#endif
