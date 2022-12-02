/*
 * This file is part of HQSpre.
 *
 * Copyright 2013-21 Tobias Schubert, Sven Reimer, Ralf Wimmer
 * Albert-Ludwigs-Universitaet Freiburg, Freiburg im Breisgau, Germany
 *
 * HQSpre is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * HQSpre is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with HQSpre. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HQSPRE_VARHEAP_HPP_
#define HQSPRE_VARHEAP_HPP_

#include <cstdint>
#include <vector>

#include "aux.hpp"


/**
 * \file varheap.hpp
 * \author Tobias Schubert
 * \author Sven Reimer
 * \author Ralf Wimmer
 * \date 2013-2021
 * \brief Heap data structure to be used as a priority queue
 */

namespace hqspre {

/**
 * \brief Used to construct a min-heap which returns the element with the lowest value
 */
template <typename Key, typename Value>
class DescendingOrder
{
   public:
    explicit DescendingOrder(const std::vector<Value>& a) : _act(a) {}

    // Do not copy or move a DescendingOrder!
    DescendingOrder(const DescendingOrder<Key, Value>& other) = delete;
    DescendingOrder(DescendingOrder<Key, Value>&& other) = delete;
    DescendingOrder<Key, Value>& operator=(const DescendingOrder<Key, Value>& other) = delete;
    DescendingOrder<Key, Value>& operator=(DescendingOrder<Key, Value>&& other) = delete;

    bool operator()(const Key x, const Key y) const noexcept { return (_act[x] > _act[y]); }

   private:
    const std::vector<Value>& _act;
};

/**
 * \brief Used to construct a max-heap which returns the element with the
 * highest value.
 */
template <typename Key, typename Value>
class AscendingOrder
{
   public:
    explicit AscendingOrder(const std::vector<Value>& a) : _act(a) {}

    // Do not copy or move an AscendingOrder!
    AscendingOrder(const AscendingOrder<Key, Value>& other) = delete;
    AscendingOrder(AscendingOrder<Key, Value>&& other) = delete;
    AscendingOrder<Key, Value>& operator=(const AscendingOrder<Key, Value>& other) = delete;
    AscendingOrder<Key, Value>& operator=(AscendingOrder<Key, Value>&& other) = delete;

    bool operator()(const Key x, const Key y) const noexcept { return (_act[x] < _act[y]); }

   private:
    const std::vector<Value>& _act;
};

template <template <typename, typename> class Order, typename Key, typename Value>
class VarHeap
{
   public:
    /**
     * \brief Constructs an empty heap with a given order.
     */
    explicit VarHeap(const std::vector<Value>& weights) :
        _heap(),
        _position(),
        _size(0),
        _variables(0),
        _comp(new Order<Key, Value>(weights))
    {}

    /**
     * \brief Copy constructor, which needs the correct weight vector!
     */
    VarHeap(const std::vector<Value>& weights, const VarHeap<Order, Key, Value>& other) :
        _heap(other._heap),
        _position(other._position),
        _size(other._size),
        _variables(other._variables),
        _comp(new Order<Key, Value>(weights))
    {}

    /**
     * \brief Move constructor, which needs the correct weight vector!
     */
    VarHeap(const std::vector<Value>& weights, VarHeap<Order, Key, Value>&& other) :
        _heap(std::move(other._heap)),
        _position(std::move(other._position)),
        _size(other._size),
        _variables(other._variables),
        _comp(new Order<Key, Value>(weights))
    {}

    /**
     * Replacement for the assignment operator, which does not work because it needs
     * the correct weight vector.
     */
    void copy(const std::vector<Value>& weights, const VarHeap<Order, Key, Value>& other)
    {
        delete _comp;
        _comp = new Order<Key, Value>(weights);
        _heap = other._heap;
        _position = other._position;
        _size = other._size;
        _variables = other._variables;
    }

    /**
     * Replacement for the move assignment operator, which does not work because it needs
     * the correct weight vector.
     */
    void move(const std::vector<Value>& weights, VarHeap<Order, Key, Value>&& other)
    {
        delete _comp;
        _comp = new Order<Key, Value>(weights);
        _heap = std::move(other._heap);
        _position = std::move(other._position);
        _size = other._size;
        _variables = other._variables;
    }

    // WARNING: Do not copy or move VarHeap objects directly using the default
    // copy or move constructors/operators. _comp contains a
    // reference to a weight vector that is used for sorting the objects
    // on the heap. If the copy constructor is used, then both the original
    // VarHeap and its copy use the same weight vector, which is typically
    // not sound!
    VarHeap(VarHeap<Order, Key, Value>&& other) noexcept = delete;
    VarHeap(const VarHeap<Order, Key, Value>& other) = delete;
    VarHeap& operator=(VarHeap<Order, Key, Value>&& other) noexcept = delete;
    VarHeap& operator=(const VarHeap& other) = delete;

    ~VarHeap() noexcept
    {
        delete _comp;
    }

    /**
     * \brief Updates all data structures to be able to handle "var" variables.
     */
    void resize(const std::size_t var)
    {
        _position.resize(var + 1, -1);
        _heap.resize(var, 0);
        _variables = var;
    }

    /**
     * \brief Clears the entire heap.
     */
    void clear() noexcept
    {
        for (std::size_t m = 0; m < _size; ++m) {
            _position[_heap[m]] = -1;
        }
        _size = 0;
    }

    /**
     * \brief Returns whether the heap is empty or not.
     */
    bool empty() const noexcept { return (_size == 0); }

    /**
     * \brief Returns the current size of `_heap`.
     */
    std::size_t size() const noexcept { return _size; }

    /**
     * \brief Returns whether "var" is an element of "_heap".
     */
    bool inHeap(const Key var) const noexcept
    {
        // "var" has to be less or equal "_variables".
        val_assert(var <= _variables);

        // Return whether "var" is part of the heap.
        return (_position[var] > -1);
    }

    /**
     * \brief Updates the position of "var" within the heap.
     */
    void update(const Key var) noexcept
    {
        // "var" has to be less or equal "_variables".
        val_assert(var <= _variables);

        // "var" has to be an element of "_heap".
        val_assert(_position[var] > -1);

        // Ensure that the heap property holds. Assumes, that "var's"
        // activity has been incremented.
        shiftUpwards(_position[var]);
    }

    /**
     * \brief Inserts variable "var" into "_heap".
     */
    void insert(const Key var) noexcept
    {
        // "var" has to be less or equal "_variables".
        val_assert(var <= _variables);

        // If "var" is an element of "_heap", we've got a problem.
        val_assert(_position[var] == -1);

        // Update "_position".
        _position[var] = static_cast<std::int32_t>(_size);

        // Add "var" to "_heap".
        _heap[_size] = var;

        // Increment "_size".
        ++_size;

        // Ensure that the heap property holds.
        shiftUpwards(_position[var]);
    }

    /**
     * \brief Returns the root variable.
     */
    Key top() noexcept
    {
        // If "_heap" is empty, we've got a problem.
        val_assert(_size > 0);

        // Get the root variable.
        const Key var = _heap[0];

        // Decrement "_size".
        --_size;

        // Overwrite "_heap[0]" with the last element of "_heap".
        _heap[0] = _heap[_size];

        // Update "_position".
        _position[var] = -1;

        // If we removed the last element from the heap, we can skip the following
        // operations.
        if (_size > 0) {
            // Update "_position".
            _position[_heap[0]] = 0;

            // Ensure that the heap property holds.
            shiftDownwards(0);
        }

        // Return "var".
        return var;
    }

    /**
     * \brief Removes "var" from the heap.
     */
    void remove(const Key var) noexcept
    {
        // If "_heap" is empty, we've got a problem.
        val_assert(_size > 0);

        // Initialization.
        const int pos = _position[var];

        // If "var" is not part of the heap, we've got a problem.
        val_assert(pos > -1);

        // Decrement "_size".
        --_size;

        // Overwrite "_heap[pos]" with the last element of "_heap".
        _heap[pos] = _heap[_size];

        // Update "_position".
        _position[var] = -1;

        // If we removed the right-most element from the heap, we can skip the
        // following operations.
        if (static_cast<std::size_t>(pos) != _size) {
            // Update "_position".
            _position[_heap[pos]] = static_cast<int>(pos);

            // Ensure that the heap property holds.
            shiftDownwards(pos);
        }

        // Consistency check.
        val_assert(_position[var] == -1);
    }


   private:
    /**
     * \brief Returns the position of the "father" of the element stored on
     * position "pos".
     */
    static Key father(const Key pos) noexcept { return ((pos - 1) >> 1); }

    /**
     * \brief Returns the position of the left "son" of the element stored on
     * position "pos".
     */
    static Key left(const Key pos) noexcept { return ((pos << 1) + 1); }

    /**
     * \brief Returns the position of the right "son" of the element stored on
     * position "pos".
     */
    static Key right(const Key pos) noexcept { return ((pos + 1) << 1); }

    /**
     * \brief Ensures the heap property by shifting the element on position "pos"
     * of "_heap" upwards.
     */
    void shiftUpwards(Key pos) noexcept
    {
        // Get the variable stored on position "pos".
        const Key var = _heap[pos];

        // Determine the correct position of "var" within "_heap".
        while (pos > 0 && (*_comp)(var, _heap[father(pos)])) {
            _heap[pos]            = _heap[father(pos)];
            _position[_heap[pos]] = static_cast<int>(pos);
            pos                   = father(pos);
        }

        // Store "var" at position "pos".
        _heap[pos]     = var;
        _position[var] = static_cast<int>(pos);
    }

    /**
     * \brief Ensures the heap property by shifting the element on position "pos"
     * of "_heap" downwards.
     */
    void shiftDownwards(Key pos)
    {
        // Get the variable stored on position "pos".
        const Key var = _heap[pos];

        // Determine the correct position of "var" within "_heap".
        while (left(pos) < _size) {
            const Key r     = right(pos);
            const Key child = r < _size && (*_comp)(_heap[r], _heap[left(pos)]) ? r : left(pos);

            if (!(*_comp)(_heap[child], var)) {
                break;
            }

            _heap[pos]            = _heap[child];
            _position[_heap[pos]] = static_cast<int>(pos);
            pos                   = child;
        }

        // Store "var" at position "pos".
        _heap[pos]     = var;
        _position[var] = static_cast<int>(pos);
    }

    /// The heap.
    std::vector<Key> _heap;

    /// The position of a particular variable within "_heap".
    std::vector<std::int32_t> _position;

    /// The current size of "_heap".
    std::size_t _size;

    /// The maximum number of variables for which memory has been reserved.
    std::size_t _variables;

    /// Ordering class
    Order<Key, Value>* _comp;
};

}  // end namespace hqspre

#endif
