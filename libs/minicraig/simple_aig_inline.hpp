/*
 * simple_aig_inline.hpp
 *
 *  Created on: Feb 17, 2009
 *      Author: Florian Pigorsch
 */

inline
minicraig::SimpleAIGEdge::
SimpleAIGEdge()
: index_( -1 )
{}

inline
minicraig::SimpleAIGEdge::
SimpleAIGEdge( int other )
: index_( other )
{
    assert( isValid() );
}


inline
minicraig::SimpleAIGEdge::
SimpleAIGEdge( const minicraig::SimpleAIGEdge& other )
: index_( other.index_ )
{
    /*assert( isValid() );*/
}

inline
minicraig::SimpleAIGEdge&
minicraig::SimpleAIGEdge::
operator=( const minicraig::SimpleAIGEdge& other )
{
    index_ = other.getRawIndex();
    /*assert( isValid() );*/
    return *this;
}

inline
bool
minicraig::SimpleAIGEdge::
operator==( const minicraig::SimpleAIGEdge& other ) const
{
    return getRawIndex() == other.getRawIndex();
}

inline
bool
minicraig::SimpleAIGEdge::
operator<( const minicraig::SimpleAIGEdge& other ) const
{
    return getRawIndex() < other.getRawIndex();
}

inline
bool
minicraig::SimpleAIGEdge::
operator>( const minicraig::SimpleAIGEdge& other ) const
{
    return getRawIndex() > other.getRawIndex();
}

inline
minicraig::SimpleAIGEdge
minicraig::SimpleAIGEdge::
operator!() const
{
    assert( isValid() );
    return SimpleAIGEdge( getRawIndex() ^ 0x1 );
}

inline
bool
minicraig::SimpleAIGEdge::
isNegated() const
{
    assert( isValid() );
    return getRawIndex() & 0x1;
}

inline
int
minicraig::SimpleAIGEdge::
getNodeIndex() const
{
    assert( isValid() );
    return getRawIndex() >> 1;
}

inline
unsigned
minicraig::SimpleAIGEdge::
hash() const
{
    return (unsigned)index_;
}

inline
bool
minicraig::SimpleAIGEdge::
isValid() const
{
    return getRawIndex() >= 0;
}

inline
void
minicraig::SimpleAIGEdge::
invalidate()
{
    index_ = -1;
}

inline
bool
minicraig::SimpleAIGEdge::
isConstant() const
{
    return getNodeIndex() == 0;
}

inline
int
minicraig::SimpleAIGEdge::
getRawIndex() const
{
    return index_;
}


inline
minicraig::SimpleAIGNode::
SimpleAIGNode()
: p1_(),
  p2_()
{}

inline
minicraig::SimpleAIGNode::
SimpleAIGNode( int variableIndex )
: p1_( variableIndex ),
  p2_( 0 )
{}

inline
minicraig::SimpleAIGNode::
SimpleAIGNode( const minicraig::SimpleAIGEdge& parent1, const minicraig::SimpleAIGEdge& parent2 )
: p1_( parent1 ),
  p2_( parent2 )
{}

bool
inline
minicraig::SimpleAIGNode::
isVar() const
{
    return ( p2().getRawIndex() == 0 );
}

inline
int
minicraig::SimpleAIGNode::
varIndex() const
{
    return p1().getRawIndex();
}

inline
minicraig::SimpleAIGEdge
minicraig::SimpleAIGNode::
p1() const
{
    return p1_;
}

inline
minicraig::SimpleAIGEdge
minicraig::SimpleAIGNode::
p2() const
{
    return p2_;
}

inline
minicraig::SimpleAIGEdge
minicraig::SimpleAIG::
getFalse()
{
    return SimpleAIGEdge( 0 );
}

inline
minicraig::SimpleAIGEdge
minicraig::SimpleAIG::
getTrue()
{
    return SimpleAIGEdge( 1 );
}

unsigned int
inline
minicraig::SimpleAIG::
size() const
{
    return nodes_.size();
}

unsigned int
inline
minicraig::SimpleAIG::
varNodes() const
{
    return _varNodes;
}

unsigned int
inline
minicraig::SimpleAIG::
andNodes() const
{
    return _andNodes;
}

inline
void
minicraig::SimpleAIG::
clear()
{
    if( nodes_.size() != 1 )
    {
        nodes_.clear();
#ifdef SIMPLEAIG_HASHLOOKUP
        hashmap_.clear();
#endif
        nodes_.push_back( SimpleAIGNode() );

        _varNodes = 0;
        _andNodes = 0;
    }
}

inline
minicraig::SimpleAIGEdge
minicraig::SimpleAIG::
Negate( const SimpleAIGEdge& e )
{
    return !e;
}

inline
minicraig::SimpleAIGEdge
minicraig::SimpleAIG::
Or( const SimpleAIGEdge& e1, const SimpleAIGEdge& e2 )
{
    return Negate( And( Negate( e1 ), Negate( e2 ) ) );
}

inline
int
minicraig::SimpleAIG::
getIndex( const minicraig::SimpleAIGEdge& e )
{
    return e.getNodeIndex();
}

inline
bool
minicraig::SimpleAIG::
isNegated( const minicraig::SimpleAIGEdge& e )
{
    return e.isNegated();
}

inline
bool
minicraig::SimpleAIG::
isConstant( const minicraig::SimpleAIGEdge& e )
{
    return e.isConstant();
}

/*
inline
bool
minicraig::SimpleAIG::
isValid( const minicraig::SimpleAIGEdge& e ) const
{
    return (
        e.isValid() &&
        getIndex( e ) < (int)nodes_.size() );
}
 */

inline
const minicraig::SimpleAIGNode&
minicraig::SimpleAIG::
getNode( int index ) const
{
    return nodes_[index];
}
