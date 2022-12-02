/*
 * SimpleAIG.hpp
 *
 *  Created on: Feb 17, 2009
 *      Author: Florian Pigorsch
 */

#ifndef MINICRAIG_SIMPLE_AIG_HPP_
#define MINICRAIG_SIMPLE_AIG_HPP_

#include <cassert>
#include <iostream>
#include <vector>
#include <tr1/unordered_map>

#define SIMPLEAIG_HASHLOOKUP

namespace minicraig {

class SimpleAIGEdge {
public:
	SimpleAIGEdge();

	explicit SimpleAIGEdge( int other );

	SimpleAIGEdge( const SimpleAIGEdge& other );

	SimpleAIGEdge& operator=( const SimpleAIGEdge& other );

	bool operator==( const SimpleAIGEdge& other ) const;

	bool operator<( const SimpleAIGEdge& other ) const;
        bool operator>( const SimpleAIGEdge& other ) const;

	SimpleAIGEdge operator!() const;
	
	bool isNegated() const;
	
	int getNodeIndex() const;

	unsigned hash() const;
	
	bool isValid() const;
	
	void invalidate();
	
	bool isConstant() const;
		
	int getRawIndex() const;
	
private:
	int index_;
};
	
class SimpleAIGNode
{
public:
    SimpleAIGNode();
    explicit SimpleAIGNode( int variableIndex );
    SimpleAIGNode( const SimpleAIGEdge& parent1, const SimpleAIGEdge& parent2 );

    bool isVar() const;
    int varIndex() const;

    SimpleAIGEdge p1() const;
    SimpleAIGEdge p2() const;

private:
    SimpleAIGEdge p1_, p2_;
};

class SimpleAIG
{
public:
    SimpleAIG();
    ~SimpleAIG();

    static SimpleAIGEdge getFalse();
    static SimpleAIGEdge getTrue();

    unsigned int size() const;
    unsigned int varNodes() const;
    unsigned int andNodes() const;

    void clear();

    SimpleAIGEdge addVar( int variableIndex );

    static SimpleAIGEdge Negate( const SimpleAIGEdge& e );
    SimpleAIGEdge And( SimpleAIGEdge e1, SimpleAIGEdge e2 );
    SimpleAIGEdge Or( const SimpleAIGEdge& e1, const SimpleAIGEdge& e2 );
    SimpleAIGEdge MultiAnd( const std::vector<SimpleAIGEdge>& edges );
    SimpleAIGEdge MultiOr( const std::vector<SimpleAIGEdge>& edges );

    void print( const SimpleAIGEdge& root, std::ostream& os = std::cout ) const;
    SimpleAIGEdge parse( std::istream& is );

    static int getIndex( const SimpleAIGEdge& e );
    static bool isNegated( const SimpleAIGEdge& e );
    static bool isConstant( const SimpleAIGEdge& e );
    //bool isValid( const SimpleAIGEdge& e ) const;

    const SimpleAIGNode& getNode( int index ) const;
    unsigned int getConeSize( const SimpleAIGEdge& e ) const;

private:
    std::vector<SimpleAIGNode> nodes_;
    unsigned int _varNodes;
    unsigned int _andNodes;

    struct hashEdgePair
    {
        size_t operator()( const std::pair<SimpleAIGEdge, SimpleAIGEdge>& p ) const
        {
            return ( p.first.hash() << 16 ) + p.second.hash();
        }
    };
    typedef std::tr1::unordered_map<std::pair<SimpleAIGEdge, SimpleAIGEdge>, unsigned int, hashEdgePair> NodeHashMap;
    NodeHashMap hashmap_;
};

}

#include "simple_aig_inline.hpp"

#endif /* SIMPLE_AIG_HPP_ */
