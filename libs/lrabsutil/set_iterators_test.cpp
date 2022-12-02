#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <list>
#include "SetIterators.hh"

int main()
{
    std::set<unsigned int> a({2,4,1,5,9});
    std::set<unsigned int> b({3,0,4,8,7,2,1});

    std::cout << "union: ";
    const auto u_end = union_iterator_end(a, b);
    for (auto it = union_iterator_begin(a, b); it != u_end; ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;

    std::cout << "intersection: ";
    const auto i_end = intersection_iterator_end(a, b);
    for (auto it = intersection_iterator_begin(a, b); it != i_end; ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;

    std::cout << "difference: ";
    const auto d_end = difference_iterator_end(b, a);
    for (auto it = difference_iterator_begin(b, a); it != d_end; ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;

    std::cout << "applying std::for_each: ";
    std::for_each(difference_iterator_begin(b, a), difference_iterator_end(b, a),
        [](unsigned int a) { std::cout << a << " "; });
    std::cout << std::endl;

    return 0;
}
