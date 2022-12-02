/*
 * simple_aig.cpp
 *
 *  Created on: Feb 17, 2009
 *      Author: Florian Pigorsch
 */

#include "simple_aig.hpp"

#include <map>
#include <set>
#include <sstream>
#include <stack>
#include <utility>
#include <vector>

namespace minicraig {
    
SimpleAIG::
SimpleAIG()
:_varNodes( 0 ),
 _andNodes( 0 )
{
	nodes_.push_back( SimpleAIGNode() );
}

SimpleAIG::
~SimpleAIG()
{}

SimpleAIGEdge
SimpleAIG::
addVar( int variableIndex )
{
    ++_varNodes;

	nodes_.push_back( SimpleAIGNode( variableIndex ) );
	return SimpleAIGEdge( ( nodes_.size() - 1 ) << 1 );
}

SimpleAIGEdge
SimpleAIG::
And(SimpleAIGEdge e1, SimpleAIGEdge e2 )
{
    assert( e1.isValid() && e2.isValid() );

	if( e1 == getFalse() || e2 == getFalse() )
	{
		return getFalse();
	}
	else if( e1 == getTrue() )
	{
		return e2;
	}
	else if( e2 == getTrue() )
	{
		return e1;
	}
	else if( e1 == e2 )
	{
		return e1;
	}
	else if( e1 == Negate( e2 ) )
	{
		return getFalse();
	}

#ifdef SIMPLEAIG_HASHLOOKUP
	if( e1 > e2 ) std::swap( e1, e2 );
	std::pair<SimpleAIGEdge, SimpleAIGEdge> ep( e1, e2 );
	NodeHashMap::const_iterator p = hashmap_.find( ep );
	if( p != hashmap_.end() )
	{
	    return SimpleAIGEdge( p->second << 1 );
	}
	else
	{
	    hashmap_[ep] = nodes_.size();
	}
#endif

	nodes_.push_back( SimpleAIGNode( e1, e2 ) );
    SimpleAIGEdge result = SimpleAIGEdge( ( nodes_.size() - 1 ) << 1 );

    ++_andNodes;

    return result;
}

SimpleAIGEdge
SimpleAIG::
MultiAnd( const std::vector<SimpleAIGEdge>& edges )
{
    for( unsigned int i = 0; i != edges.size(); ++i )
    {
        assert( edges[i].isValid() );
    }
    
	if( edges.empty() )
	{
		return getTrue();
	}
	else if( edges.size() == 1 )
	{
		return edges[0];
	}
	else
	{
		SimpleAIGEdge currentEdge = getTrue();

		for( unsigned int i = 0; i != edges.size(); ++i )
		{
			currentEdge = And( currentEdge, edges[i] );
		}

		return currentEdge;
	}
}

SimpleAIGEdge
SimpleAIG::
MultiOr( const std::vector<SimpleAIGEdge>& edges )
{
    for( unsigned int i = 0; i != edges.size(); ++i )
    {
        assert( edges[i].isValid() );
    }

	if( edges.empty() )
	{
		return getFalse();
	}
	else if( edges.size() == 1 )
	{
		return edges[0];
	}
	else
	{
		SimpleAIGEdge currentEdge = getFalse();

		for( unsigned int i = 0; i != edges.size(); ++i )
		{
			currentEdge = Or( currentEdge, edges[i] );
		}

		return currentEdge;
	}
}

void
SimpleAIG::
print( const SimpleAIGEdge& root, std::ostream& os ) const
{
    if( !root.isValid() )
    {
        std::cout << "root INVALID!\n" << std::endl;
        return;
    }
    
	if( root == getFalse() )
	{
		os << "root FALSE\n";
		return;
	}
	else if( root == getTrue() )
	{
		os << "root TRUE\n";
		return;
	}

	std::stack<int> pending;
	std::set<int> processed;

	os << "root " << ( isNegated( root ) ? "!" : "" ) << getIndex( root ) << "\n";

	pending.push( getIndex( root ) );

	while( !pending.empty() )
	{
		if( processed.find( pending.top() ) != processed.end() )
		{
			pending.pop();
			continue;
		}

		int index = pending.top();
		pending.pop();
		const SimpleAIGNode& n = getNode( index );

		processed.insert( index );

		if( n.isVar() )
		{
			os << index << " var( " << n.varIndex() << " )\n";
		}
		else
		{
			os
			<< index << " and( "
			<< ( isNegated( n.p1() ) ? "!" : "" ) << getIndex( n.p1() )
			<< " "
			<< ( isNegated( n.p2() ) ? "!" : "" ) << getIndex( n.p2() )
			<< " )\n";

			pending.push( getIndex( n.p1() ) );
			pending.push( getIndex( n.p2() ) );
		}
	}
}

int string2int( const std::string& s )
{
    std::istringstream iss;
    iss.str( s );
        
    int value;
    iss >> value;
    
    // s was not parsed completely
    if( !( iss.eof() ) )
    {
        std::cout << "ERROR: cannot parse integer '" << s << "'" << std::endl;
        throw;
    }
    
    return value;
}

std::vector<std::string>
splitstring( const std::string& input )
{
    std::vector<std::string> result;

    int lastpos = 0;
    while( true )
    {
        int pos = input.find( ' ', lastpos );
        if( pos < 0 )
        {
            result.push_back( input.substr( lastpos ) );
            break;
        }
        else
        {
            result.push_back( input.substr( lastpos, pos - lastpos ) );
            lastpos = pos + 1;
        }
    }
    
    return result;
}

int negatedstring2int(  const std::string& input )
{
    assert( input.size() > 0 );

    if( input[0] == '!' )
    {
        assert( input.size() > 1 );
        return 1 + 2 * string2int( input.substr( 1 ) );
    }
    else
    {
        return 2 * string2int( input );
    }
}

SimpleAIGEdge
SimpleAIG::
parse( std::istream& is )
{
    clear();
    
    std::string s;
    std::vector<std::string> sp;

    
    std::getline( is, s );
    sp = splitstring( s );
    assert( sp.size() == 2 );
    assert( sp[0] == "root" );
    
    if( sp[1] == "TRUE" )
    {
        return getTrue();
    }
    else if( sp[1] == "FALSE" )
    {
        return getFalse();
    }
    
    int root = negatedstring2int( sp[1] );
    
    std::map<int, std::pair<int, int> > ands;
    std::map<int, SimpleAIGEdge> mapping;

    while( std::getline( is, s ) )
    {
        /* X var( Y ) */
        /* X and( Y, Z ) */
        
        sp = splitstring( s );

        if( sp.size() != 4 && sp.size() != 5 )
        {
            std::cout << "ERROR: bad number of tokens in line '" << s << "'" << std::endl;
        }
        assert( sp.size() == 4 || sp.size() == 5 );
        assert( sp[0].size() > 0 );
        
        int index = string2int( sp[0] );
        assert( mapping.find( index ) == mapping.end() );
        assert( ands.find( index ) == ands.end() );
        
        if( sp.size() == 4 )
        {
            assert( sp[1] == "var(" && sp[2].size() > 0 && sp[3] == ")" );
            int varindex = string2int( sp[2] );
            mapping[index] = addVar( varindex );
        }
        else
        {
            assert( sp[1] == "and(" && sp[2].size() > 0 && sp[3].size() > 0 && sp[4] == ")" );

            /* sp[2] may have a trailing ',' -> remove it */
            if( sp[2][sp[2].size()-1] == ',' )
            {
                sp[2].erase( sp[2].size() - 1 );
            }
            
            int p1 = negatedstring2int( sp[2] );
            int p2 = negatedstring2int( sp[3] );
            
            ands[index] = std::pair<int, int>( p1, p2 );
        }
    }

    std::stack<int> pending;
    pending.push( root >> 1 );
    while( !pending.empty() )
    {
        if( mapping.find( pending.top() ) != mapping.end() )
        {
            pending.pop();
        }
        else
        {
            std::map<int, std::pair<int, int> >::const_iterator a = ands.find( pending.top() );
            assert( a != ands.end() );

            
            std::map<int, SimpleAIGEdge>::const_iterator p1 = mapping.find( a->second.first >> 1 );
            if( p1 == mapping.end() )
            {
                pending.push( a->second.first >> 1 );
                continue;
            }

            std::map<int, SimpleAIGEdge>::const_iterator p2 = mapping.find( a->second.second >> 1 );
            if( p2 == mapping.end() )
            {
                pending.push( a->second.second >> 1 );
                continue;
            }

            mapping[pending.top()] = And( ( a->second.first & 1 ) ? Negate( p1->second ) : p1->second,
                                          ( a->second.second & 1 ) ? Negate( p2->second ) : p2->second );
            
            pending.pop();
        }
    }

    std::map<int, SimpleAIGEdge>::const_iterator proot = mapping.find(root >> 1);
    assert( proot != mapping.end() );

    SimpleAIGEdge result = ( ( root & 1 ) ? Negate( proot->second ) : proot->second );
    
    return result;
}

unsigned int
SimpleAIG::
getConeSize( const SimpleAIGEdge& root ) const
{
    assert( root.isValid() );
    
	if( root == getFalse() || root == getTrue() )
	{
		return 0;
	}

	std::stack<int> pending;
	std::set<int> processed;

	pending.push( getIndex( root ) );

	while( !pending.empty() )
	{
		if( processed.find( pending.top() ) != processed.end() )
		{
			pending.pop();
			continue;
		}

		int index = pending.top();
		pending.pop();
		const SimpleAIGNode& n = getNode( index );

		processed.insert( index );

		if( !n.isVar() )
		{
			pending.push( getIndex( n.p1() ) );
			pending.push( getIndex( n.p2() ) );
		}
	}

    return processed.size();
}

}
