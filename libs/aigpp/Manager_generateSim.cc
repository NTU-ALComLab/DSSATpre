/**************************************************************
 *
 *       AIGPP Package // Manager_generateSim.cc
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 373 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#include "Manager.hh"

#include <vector>

#include <lrabsutil/Math.hh>
#include <lrabsutil/Random.hh>
#include <lrabsutil/Resources.hh>

#if 0
namespace aigpp
{
    enum TriVal
    {
        C0 = 0, /* 00 */
        C1 = 3, /* 11 */
        CX = 1  /* 01 */
    };

    inline
    void setTriVec( unsigned int& b, int index, TriVal t )
    {
        b &= ~( 3 << ( 2 * index ) );
        b |= t << ( 2 * index );
    }

    inline
    void set01( unsigned int& b, int index, bool t )
    {
        b &= ~( 1 << index );

        if( t )
        {
            b |= ( 1 << index );
        }
    }

    inline
    void setTriVec( aigpp::Node* n, int index, TriVal t )
    {
        unsigned int a = n->tempSim();
        setTriVec( a, index, t );
        n->setTempSim( a );
    }

    inline
    void set01( aigpp::Node* n, int index, bool t )
    {
        unsigned int a = n->tempSim();
        set01( a, index, t );
        n->setTempSim( a );
    }

    inline
    void setTriVec( aigpp::Node* n, TriVal t )
    {
        if( t == C0 )
        {
            n->setTempSim( 0x00000000 );
        }
        else if( t == C1 )
        {
            n->setTempSim( 0xFFFFFFFF );
        }
        else
        {
            n->setTempSim( 0x55555555 );
        }
    }

    inline
    void set01( aigpp::Node* n, bool t )
    {
        if( t )
        {
            n->setTempSim( 0xFFFFFFFF );
        }
        else
        {
            n->setTempSim( 0x00000000 );
        }
    }

    inline
    TriVal getTriVec( aigpp::Node* n, int index )
    {
        return TriVal( ( n->tempSim() >> ( 2 * index ) ) & 3 );
    }

    inline
    bool get01( aigpp::Node* n, int index )
    {
        return ( ( n->tempSim() >> index ) & 1 );
    }

    inline
    void getTriVec( aigpp::Node* n, TriVal* val )
    {
        unsigned int v = n->tempSim();

        for( int i = 0; i != 16; ++i )
        {
            val[ i ] = TriVal( v & 3 );
            v >>= 2;
        }
    }

    inline
    unsigned int negateTriVec( unsigned int v )
    {

        return ~( ( ( v & 0xAAAAAAAA ) >> 1 ) |
                  ( ( v & 0x55555555 ) << 1 ) );
    }

    inline
    void updateTriVec( aigpp::Node* n )
    {
        unsigned int a = n->parent1().node()->tempSim();
        unsigned int b = n->parent2().node()->tempSim();

        if( a != 0x55555555 ||
            b != 0x55555555 )
        {
            if( n->parent1().isInverted() && a != 0x55555555 ) a = negateTriVec( a );
            if( n->parent2().isInverted() && b != 0x55555555 ) b = negateTriVec( b );

            a = a & b;
            if( n->isNAND() && a != 0x55555555 ) a = negateTriVec( a );
        }

        n->setTempSim( a );
    }

    inline
    void update01( aigpp::Node* n )
    {
        unsigned int a = n->parent1().node()->tempSim();
        unsigned int b = n->parent2().node()->tempSim();

        if( n->parent1().isInverted() ) a = ~a;
        if( n->parent2().isInverted() ) b = ~b;

        a = a & b;
        if( n->isNAND() ) a = ~a;

        n->setTempSim( a );
    }

    void propagate01X( aigpp::Node* nodes )
    {
        // std::cout << __func__ << std::endl;

        for( Node* n = nodes; n != 0; n = n->next() )
        {
            if( n->isVar() ) continue;

            updateTriVec( n );
        }
    }

    void propagate01( aigpp::Node* nodes )
    {
        // std::cout << __func__ << std::endl;

        for( Node* n = nodes; n != 0; n = n->next() )
        {
            if( n->isVar() ) continue;

            update01( n );
        }
    }

    struct ListItem
    {
        Node* n;
        ListItem* next;
    };

    struct Class
    {
        Node** nodes;
        int size;
    };

    struct Classes
    {
        Class* classes;
        int size;
    };

    void propagate01X( ListItem* unassignedNodes )
    {
        // std::cout << __func__ << std::endl;
        for( ListItem* n = unassignedNodes; n != 0; n = n->next )
        {
            if( n->n->isVar() ) continue;

            updateTriVec( n->n );
        }
    }

    ListItem* propagate01XAndUpdateList( ListItem* unassignedNodes )
    {
        // std::cout << __func__ << std::endl;

        ListItem* newList = 0;
        ListItem* lastItem = 0;
        ListItem* next;

        for( ListItem* n = unassignedNodes; n != 0; )
        {
            next = n->next;
            n->next = 0;

            if( !( n->n->isVar() ) )
            {
                updateTriVec( n->n );
            }

            if( n->n->tempSim() == 0x55555555 )
            {
                if( lastItem != 0 )
                {
                    lastItem->next = n;
                }
                else
                {
                    newList = n;
                }

                lastItem = n;
            }
            else
            {
                delete n;
            }

            n = next;
        }

        return newList;
    }

    ListItem* initializeUnassigned( Node* nodes )
    {
        /* unary nodes are marked with flag #256 */
        /* mark non toplevel nodes with flag #512 */

        while( true )
        {
            bool again = false;

            for( Node* n = nodes; n != 0; n = n->next() )
            {
                if( n->isVar() ) continue;
                if( n->flag<1024>() ) continue;

                n->parent1().node()->setFlag<512>();
                n->parent2().node()->setFlag<512>();
            }

            for( Node* n = nodes; n != 0; n = n->next() )
            {
                if( n->flag<256>() && !( n->flag<512>() ) && !( n->flag<1024>() ) )
                {
                    again = true;
                    n->setFlag<1024>();
                }

                n->unsetFlag<512>();
            }

            if( !again ) break;
        }

        ListItem* unassignedNodes = 0;
        ListItem* lastItem = 0;
        for( Node* n = nodes; n != 0; n = n->next() )
        {
            if( n->flag<1024>() ) continue;

            ListItem* newItem = new ListItem;
            newItem->n = n;
            newItem->next = 0;

            if( lastItem != 0 )
            {
                lastItem->next = newItem;
            }
            else
            {
                unassignedNodes = newItem;
            }

            lastItem = newItem;
        }

        return unassignedNodes;
    }


    ListItem* updateUnassigned( ListItem* nodes )
    {
        // std::cout << __func__ << std::endl;
        ListItem* newList = 0;
        ListItem* lastItem = 0;
        ListItem* next;

        for( ListItem* n = nodes; n != 0; )
        {
            next = n->next;
            n->next = 0;


            if( n->n->tempSim() == 0x55555555 )
            {
                if( lastItem != 0 )
                {
                    lastItem->next = n;
                }
                else
                {
                    newList = n;
                }

                lastItem = n;
            }
            else
            {
                delete n;
            }

            n = next;
        }

        return newList;
    }

    void deleteList( ListItem* n )
    {
        // std::cout << __func__ << std::endl;
        ListItem* next;

        while( n != 0 )
        {
            next = n->next;
            delete n;
            n = next;
        }
    }

    void calcScoreByNewSimParallel01_32( const Classes& classes, double* score )
    {
        // std::cout << __func__ << std::endl;

        const int combinations = 32;

        for( int i = 0; i != combinations; ++i )
        {
            score[ i ] = 0.0;
        }


        Class* classit = classes.classes;
        for( int i = classes.size; i != 0; --i, ++classit )
        {
            int count0[combinations] = {0,0,0,0,0,0,0,0,
                                        0,0,0,0,0,0,0,0,
                                        0,0,0,0,0,0,0,0,
                                        0,0,0,0,0,0,0,0};
            int count1[combinations] = {0,0,0,0,0,0,0,0,
                                        0,0,0,0,0,0,0,0,
                                        0,0,0,0,0,0,0,0,
                                        0,0,0,0,0,0,0,0};

            Node** n = classit->nodes;
            for( int j = classit->size; j != 0; --j, ++n )
            {
                const unsigned int v = ( *n )->tempSim();

                for( int c = 0; c != combinations; ++c )
                {
                    if( ( v & ( 1 << c ) ) != 0 )
                    {
                        ++count1[c];
                    }
                    else
                    {
                        ++count0[c];
                    }
                }
            }

            for( int c = 0; c != combinations; ++c )
            {
                score[c] += double( count0[c]*(count0[c]-1) + count1[c]*(count1[c]-1) );
            }
        }

        for( int c = 0; c != combinations; ++c )
        {
            score[ c ] *= 0.5;
        }
    }

#    define INCCOUNT(sim, mask, count) \
        switch (sim & mask) {          \
            case 0:                    \
                ++count##_0;           \
                break;                 \
            case mask:                 \
                ++count##_1;           \
        }
#    define INCCOUNT2(sim, mask, count0, count1) \
        switch (sim & mask) {                    \
            case 0:                              \
                ++count0;                        \
                break;                           \
            case mask:                           \
                ++count1;                        \
        }
#    define UPDSCORE(score, count) (score += double(count##_0 * (count##_0 - 1) + count##_1 * (count##_1 - 1)))
#    define UPDSCORE2(score, count, all)                                        \
        int count##_x2 = (all - count##_0 - count##_1) / 4;                     \
        score += double((count##_0 + count##_x2) * (count##_0 + count##_x2 - 1) \
                        + (count##_1 + count##_x2) * (count##_1 + count##_x2 - 1))
#    define UPDSCORE3(score, count0, count1, all)                                                                   \
        {                                                                                                           \
            const int countX = (all - count0 - count1) / 4;                                                         \
            score += double((count0 + countX) * (count0 + countX - 1) + (count1 + countX) * (count1 + countX - 1)); \
        }

//#define UPDSCORE2(score,count,all) UPDSCORE(score,count)

    void calcScoreByNewSimParallel4( const Classes& classes, double* score )
    {
        // std::cout << __func__ << std::endl;

        const int combinations = 16;

        for( int i = 0; i != combinations; ++i )
        {
            score[ i ] = 0.0;
        }


        Class* classit = classes.classes;
        for( int i = classes.size; i != 0; --i, ++classit )
        {
            int count0_0 = 0, count0_1 = 0;
            int count1_0 = 0, count1_1 = 0;
            int count2_0 = 0, count2_1 = 0;
            int count3_0 = 0, count3_1 = 0;
            int count4_0 = 0, count4_1 = 0;
            int count5_0 = 0, count5_1 = 0;
            int count6_0 = 0, count6_1 = 0;
            int count7_0 = 0, count7_1 = 0;
            int count8_0 = 0, count8_1 = 0;
            int count9_0 = 0, count9_1 = 0;
            int countA_0 = 0, countA_1 = 0;
            int countB_0 = 0, countB_1 = 0;
            int countC_0 = 0, countC_1 = 0;
            int countD_0 = 0, countD_1 = 0;
            int countE_0 = 0, countE_1 = 0;
            int countF_0 = 0, countF_1 = 0;
            /*
            int count0[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            int count1[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            */
            Node** n = classit->nodes;
            for( int j = classit->size; j != 0; --j, ++n )
            {
                const unsigned int v = ( *n )->tempSim();

                INCCOUNT( v,        0x3, count0 );
                INCCOUNT( v,        0xC, count1 );
                INCCOUNT( v,       0x30, count2 );
                INCCOUNT( v,       0xC0, count3 );
                INCCOUNT( v,      0x300, count4 );
                INCCOUNT( v,      0xC00, count5 );
                INCCOUNT( v,     0x3000, count6 );
                INCCOUNT( v,     0xC000, count7 );
                INCCOUNT( v,    0x30000, count8 );
                INCCOUNT( v,    0xC0000, count9 );
                INCCOUNT( v,   0x300000, countA );
                INCCOUNT( v,   0xC00000, countB );
                INCCOUNT( v,  0x3000000, countC );
                INCCOUNT( v,  0xC000000, countD );
                INCCOUNT( v, 0x30000000, countE );
                INCCOUNT( v, 0xC0000000, countF );

                /*
                INCCOUNT2( v,        0x3, count0[0x0], count1[0x0] );
                INCCOUNT2( v,        0xC, count0[0x1], count1[0x1] );
                INCCOUNT2( v,       0x30, count0[0x2], count1[0x2] );
                INCCOUNT2( v,       0xC0, count0[0x3], count1[0x3] );
                INCCOUNT2( v,      0x300, count0[0x4], count1[0x4] );
                INCCOUNT2( v,      0xC00, count0[0x5], count1[0x5] );
                INCCOUNT2( v,     0x3000, count0[0x6], count1[0x6] );
                INCCOUNT2( v,     0xC000, count0[0x7], count1[0x7] );
                INCCOUNT2( v,    0x30000, count0[0x8], count1[0x8] );
                INCCOUNT2( v,    0xC0000, count0[0x9], count1[0x9] );
                INCCOUNT2( v,   0x300000, count0[0xA], count1[0xA] );
                INCCOUNT2( v,   0xC00000, count0[0xB], count1[0xB] );
                INCCOUNT2( v,  0x3000000, count0[0xC], count1[0xC] );
                INCCOUNT2( v,  0xC000000, count0[0xD], count1[0xD] );
                INCCOUNT2( v, 0x30000000, count0[0xE], count1[0xE] );
                INCCOUNT2( v, 0xC0000000, count0[0xF], count1[0xF] );
                */
            }

            const int count = classit->size;

            UPDSCORE2( score[ 0x0 ], count0, count );
            UPDSCORE2( score[ 0x1 ], count1, count );
            UPDSCORE2( score[ 0x2 ], count2, count );
            UPDSCORE2( score[ 0x3 ], count3, count );
            UPDSCORE2( score[ 0x4 ], count4, count );
            UPDSCORE2( score[ 0x5 ], count5, count );
            UPDSCORE2( score[ 0x6 ], count6, count );
            UPDSCORE2( score[ 0x7 ], count7, count );
            UPDSCORE2( score[ 0x8 ], count8, count );
            UPDSCORE2( score[ 0x9 ], count9, count );
            UPDSCORE2( score[ 0xA ], countA, count );
            UPDSCORE2( score[ 0xB ], countB, count );
            UPDSCORE2( score[ 0xC ], countC, count );
            UPDSCORE2( score[ 0xD ], countD, count );
            UPDSCORE2( score[ 0xE ], countE, count );
            UPDSCORE2( score[ 0xF ], countF, count );

            /*
            UPDSCORE3( score[ 0x0 ], count0[0x0], count1[0x0], count );
            UPDSCORE3( score[ 0x1 ], count0[0x1], count1[0x1], count );
            UPDSCORE3( score[ 0x2 ], count0[0x2], count1[0x2], count );
            UPDSCORE3( score[ 0x3 ], count0[0x3], count1[0x3], count );
            UPDSCORE3( score[ 0x4 ], count0[0x4], count1[0x4], count );
            UPDSCORE3( score[ 0x5 ], count0[0x5], count1[0x5], count );
            UPDSCORE3( score[ 0x6 ], count0[0x6], count1[0x6], count );
            UPDSCORE3( score[ 0x7 ], count0[0x7], count1[0x7], count );
            UPDSCORE3( score[ 0x8 ], count0[0x8], count1[0x8], count );
            UPDSCORE3( score[ 0x9 ], count0[0x9], count1[0x9], count );
            UPDSCORE3( score[ 0xA ], count0[0xA], count1[0xA], count );
            UPDSCORE3( score[ 0xB ], count0[0xB], count1[0xB], count );
            UPDSCORE3( score[ 0xC ], count0[0xC], count1[0xC], count );
            UPDSCORE3( score[ 0xD ], count0[0xD], count1[0xD], count );
            UPDSCORE3( score[ 0xE ], count0[0xE], count1[0xE], count );
            UPDSCORE3( score[ 0xF ], count0[0xF], count1[0xF], count );
            */
        }

        for( int i = 0; i != combinations; ++i )
        {
            score[ i ] *= 0.5;
        }
    }

    void calcScoreByNewSimParallel3( const Classes& classes, double* score )
    {
        // std::cout << __func__ << std::endl;

        const int combinations = 8;

        for( int i = 0; i != combinations; ++i )
        {
            score[ i ] = 0.0;
        }


        Class* classit = classes.classes;
        for( int i = classes.size; i != 0; --i, ++classit )
        {
            int count0_0 = 0, count0_1 = 0;
            int count1_0 = 0, count1_1 = 0;
            int count2_0 = 0, count2_1 = 0;
            int count3_0 = 0, count3_1 = 0;
            int count4_0 = 0, count4_1 = 0;
            int count5_0 = 0, count5_1 = 0;
            int count6_0 = 0, count6_1 = 0;
            int count7_0 = 0, count7_1 = 0;

            Node** n = classit->nodes;
            for( int j = classit->size; j != 0; --j, ++n )
            {
                const unsigned int v = ( *n )->tempSim();

                INCCOUNT( v,        0x3, count0 );
                INCCOUNT( v,        0xC, count1 );
                INCCOUNT( v,       0x30, count2 );
                INCCOUNT( v,       0xC0, count3 );
                INCCOUNT( v,      0x300, count4 );
                INCCOUNT( v,      0xC00, count5 );
                INCCOUNT( v,     0x3000, count6 );
                INCCOUNT( v,     0xC000, count7 );
            }

            const int count = classit->size;
            UPDSCORE2( score[ 0x0 ], count0, count );
            UPDSCORE2( score[ 0x1 ], count1, count );
            UPDSCORE2( score[ 0x2 ], count2, count );
            UPDSCORE2( score[ 0x3 ], count3, count );
            UPDSCORE2( score[ 0x4 ], count4, count );
            UPDSCORE2( score[ 0x5 ], count5, count );
            UPDSCORE2( score[ 0x6 ], count6, count );
            UPDSCORE2( score[ 0x7 ], count7, count );
        }

        for( int i = 0; i != combinations; ++i )
        {
            score[ i ] *= 0.5;
        }
    }

    void calcScoreByNewSimParallel2( const Classes& classes, double* score )
    {
        // std::cout << __func__ << std::endl;

        const int combinations = 4;

        for( int i = 0; i != combinations; ++i )
        {
            score[ i ] = 0.0;
        }


        Class* classit = classes.classes;
        for( int i = classes.size; i != 0; --i, ++classit )
        {
            int count0_0 = 0, count0_1 = 0;
            int count1_0 = 0, count1_1 = 0;
            int count2_0 = 0, count2_1 = 0;
            int count3_0 = 0, count3_1 = 0;

            Node** n = classit->nodes;
            for( int j = classit->size; j != 0; --j, ++n )
            {
                const unsigned int v = ( *n )->tempSim();

                INCCOUNT( v,        0x3, count0 );
                INCCOUNT( v,        0xC, count1 );
                INCCOUNT( v,       0x30, count2 );
                INCCOUNT( v,       0xC0, count3 );
            }

            const int count = classit->size;
            UPDSCORE2( score[ 0x0 ], count0, count );
            UPDSCORE2( score[ 0x1 ], count1, count );
            UPDSCORE2( score[ 0x2 ], count2, count );
            UPDSCORE2( score[ 0x3 ], count3, count );
        }

        for( int i = 0; i != combinations; ++i )
        {
            score[ i ] *= 0.5;
        }
    }

    void calcScoreByNewSimParallel1( const Classes& classes, double* score )
    {
        // std::cout << __func__ << std::endl;

        const int combinations = 2;

        for( int i = 0; i != combinations; ++i )
        {
            score[ i ] = 0.0;
        }


        Class* classit = classes.classes;
        for( int i = classes.size; i != 0; --i, ++classit )
        {
            int count0_0 = 0, count0_1 = 0;
            int count1_0 = 0, count1_1 = 0;

            Node** n = classit->nodes;
            for( int j = classit->size; j != 0; --j, ++n )
            {
                const unsigned int v = ( *n )->tempSim();

                INCCOUNT( v,        0x3, count0 );
                INCCOUNT( v,        0xC, count1 );
            }

            const int count = classit->size;
            UPDSCORE2( score[ 0x0 ], count0, count );
            UPDSCORE2( score[ 0x1 ], count1, count );
        }

        for( int i = 0; i != combinations; ++i )
        {
            score[ i ] *= 0.5;
        }
    }

    double calcScore( const Classes& classes )
    {
        // std::cout << __func__ << std::endl;
        double score = 0;

        Class* classit = classes.classes;
        for( int i = classes.size; i != 0; --i, ++classit )
        {
            const int count = classit->size;
            score += double( count * ( count - 1 ) );
        }

        return score * 0.5;
    }

    void initializeClasses( Classes& classes, Node* nodes, int nodecount )
    {
        // std::cout << __func__ << std::endl;
        classes.size = 1;
        classes.classes = new Class[1];

        Class& c = classes.classes[0];
        c.size = 0;
        c.nodes = new Node*[nodecount];

        for( Node* n = nodes; n != 0; n = n->next() )
        {
            c.nodes[ c.size++ ] = n;
        }
    }

    void deleteClasses( Classes& classes )
    {
        // std::cout << __func__ << std::endl;

        Class* classit = classes.classes;
        for( int i = classes.size; i != 0; --i, ++classit )
        {
            delete[] classit->nodes;
        }
        delete[] classes.classes;
    }

    void updateClasses( Classes& classes, int nodecount )
    {
        //std::cout << __func__ << std::endl;

        Class* newClasses = new Class[ classes.size * 2 ];
        int newSize = 0;

        Class c0;
        c0.nodes = new Node*[ nodecount ];
        c0.size = 0;
        Class c1;
        c1.nodes = new Node*[ nodecount ];
        c1.size = 0;

        Class* classit = classes.classes;
        for( int i = classes.size; i != 0; --i, ++classit )
        {
            c0.size = 0;
            c1.size = 0;

            Node** n = classit->nodes;
            for( int j = classit->size; j != 0; --j, ++n )
            {
                if( ( *n )->tempSim() == 0x00000000 )
                {
                    c0.nodes[ c0.size++ ] = *n;
                }
                else
                {
                    c1.nodes[ c1.size++ ] = *n;
                }
            }

            if( c0.size >= 2 )
            {
                Class& c = newClasses[ newSize++ ];
                c.size = c0.size;
                c.nodes = new Node*[ c.size ];
                memcpy( c.nodes, c0.nodes, c.size * sizeof( Node* ) );
            }
            else if( c0.size == 1 )
            {
                c0.nodes[ 0 ]->setFlag<256>();
            }

            if( c1.size >= 2 )
            {
                Class& c = newClasses[ newSize++ ];
                c.size = c1.size;
                c.nodes = new Node*[ c.size ];
                memcpy( c.nodes, c1.nodes, c.size * sizeof( Node* ) );
            }
            else if( c1.size == 1 )
            {
                c1.nodes[ 0 ]->setFlag<256>();
            }
        }

        delete[] c0.nodes;
        delete[] c1.nodes;

        deleteClasses( classes );

        classes.size = newSize;
        classes.classes = newClasses;
    }

    void updateClasses01( Classes& classes, int nodecount )
    {
        //std::cout << __func__ << std::endl;

        Class* newClasses = new Class[ classes.size * 2 ];
        int newSize = 0;

        Class c0;
        c0.nodes = new Node*[ nodecount ];
        c0.size = 0;
        Class c1;
        c1.nodes = new Node*[ nodecount ];
        c1.size = 0;

        Class* classit = classes.classes;
        for( int i = classes.size; i != 0; --i, ++classit )
        {
            c0.size = 0;
            c1.size = 0;

            Node** n = classit->nodes;
            for( int j = classit->size; j != 0; --j, ++n )
            {
                if( ( *n )->tempSim() == 0x00000000 )
                {
                    c0.nodes[ c0.size++ ] = *n;
                }
                else
                {
                    assert( ( *n )->tempSim() == 0xFFFFFFFF );

                    c1.nodes[ c1.size++ ] = *n;
                }
            }

            if( c0.size >= 2 )
            {
                Class& c = newClasses[ newSize++ ];
                c.size = c0.size;
                c.nodes = new Node*[ c.size ];
                memcpy( c.nodes, c0.nodes, c.size * sizeof( Node* ) );
            }

            if( c1.size >= 2 )
            {
                Class& c = newClasses[ newSize++ ];
                c.size = c1.size;
                c.nodes = new Node*[ c.size ];
                memcpy( c.nodes, c1.nodes, c.size * sizeof( Node* ) );
            }
        }

        delete[] c0.nodes;
        delete[] c1.nodes;

        deleteClasses( classes );

        classes.size = newSize;
        classes.classes = newClasses;
    }

    void duplicate01X( ListItem* nodes, int index )
    {
        // std::cout << __func__ << std::endl;
        for( ListItem* n = nodes; n != 0; n = n->next )
        {
            if( !( n->n->isVar() ) )
            {
                setTriVec( n->n, getTriVec( n->n, index ) );
            }
        }
    }

    void duplicate01X( Node* nodes, int index )
    {
        // std::cout << __func__ << std::endl;
        for( Node* n = nodes; n != 0; n = n->next() )
        {
            if( !( n->isVar() ) )
            {
                setTriVec( n, getTriVec( n, index ) );
            }
        }
    }

    void duplicate01( Node* nodes, int index )
    {
        // std::cout << __func__ << std::endl;
        for( Node* n = nodes; n != 0; n = n->next() )
        {
            set01( n, get01( n, index ) );
        }
    }

    ListItem* duplicate01XAndUpdateList( ListItem* nodes, int index )
    {
        // std::cout << __func__ << std::endl;

        ListItem* newList = 0;
        ListItem* lastItem = 0;
        ListItem* next;

        for( ListItem* n = nodes; n != 0; )
        {
            next = n->next;
            n->next = 0;

            if( !( n->n->isVar() ) )
            {
                setTriVec( n->n, getTriVec( n->n, index ) );
            }

            if( n->n->tempSim() == 0x55555555 )
            {
                if( lastItem != 0 )
                {
                    lastItem->next = n;
                }
                else
                {
                    newList = n;
                }

                lastItem = n;
            }
            else
            {
                delete n;
            }

            n = next;
        }

        return newList;
    }

    void clearTempSim( Node* nodes )
    {
        // std::cout << __func__ << std::endl;
        for( Node* n = nodes; n != 0; n = n->next() )
        {
            if( n->isVar() )
            {
                n->setTempSim( 0 );
            }
            else
            {
                n->updateTempSim();
            }
        }
    }

    void setAssignment01( VariableAssignment& va, Node* var, bool value )
    {
        switch( value )
        {
        case false:
            va.setNegative( var->varIndex() );
            break;
        case true:
            va.setPositive( var->varIndex() );
            break;
        }
    }

    void setAssignment( VariableAssignment& va, Node* var, TriVal value )
    {
        switch( value )
        {
        case C0:
            va.setNegative( var->varIndex() );
            break;
        case C1:
            va.setPositive( var->varIndex() );
            break;
        case CX:
            va.setUnassigned( var->varIndex() );
            break;
        default:
            NEVER_GET_HERE;
        }
    }

}

double
aigpp::Manager::
generateBalancedSim( int numberOfPatterns, double minRandomRatio, double maxRandomRatio )
{
    // std::cout << __func__ << std::endl;
    double startTime = lrabs::cpuTime();

    /* temp sim is not valid after this function! */
    _tempSimUpToDate = false;


    lockFlags<256|512|1024>();

    int randomUpdates = 0;
    int parallelUpdates = 0;
    int singleUpdates = 0;
    double timeForRandomUpdates = 0;
    double timeForParallelUpdates = 0;
    double timeForSingleUpdates = 0;


    lrabs::Random rnd( /*seed = */42 );


    TriVal value0[ 16 ] = { C0, C1, C0, C1, C0, C1, C0, C1, C0, C1, C0, C1, C0, C1, C0, C1 };
    TriVal value1[ 16 ] = { C0, C0, C1, C1, C0, C0, C1, C1, C0, C0, C1, C1, C0, C0, C1, C1 };
    TriVal value2[ 16 ] = { C0, C0, C0, C0, C1, C1, C1, C1, C0, C0, C0, C0, C1, C1, C1, C1 };
    TriVal value3[ 16 ] = { C0, C0, C0, C0, C0, C0, C0, C0, C1, C1, C1, C1, C1, C1, C1, C1 };

    unsigned int
        bstring0 = 0,
        bstring1 = 0,
        bstring2 = 0,
        bstring3 = 0;

    for( int i = 0; i != 16; ++i )
    {
        setTriVec( bstring0, i, value0[ i ] );
        setTriVec( bstring1, i, value1[ i ] );
        setTriVec( bstring2, i, value2[ i ] );
        setTriVec( bstring3, i, value3[ i ] );
    }


    Classes classes;
    initializeClasses( classes, _nodes, nodeCount() );

    double score = calcScore( classes );
    std::cout << "bscore_" << minRandomRatio << "_" << maxRandomRatio << " (initial) = " << score << std::endl;

    double scores[16];

    std::vector<Node*> variableOrder( _variables.size(), 0 );
    VariableAssignment va( _variables.size() );

    for( int p = 0; p != numberOfPatterns; ++p )
    {
        /* create random variable order */
        variableOrder.assign( _variables.size(), 0 );
        for( std::vector<Node*>::iterator v = _variables.begin();
             v != _variables.end(); ++v )
        {
            unsigned int pos;
            do { pos = rnd.getUInt( _variables.size() ); } while( variableOrder[ pos ] != 0 );

            variableOrder[ pos ] = *v;
        }

        /* clear tempsim ( set X ) */
        for( std::vector<Node*>::iterator v = _variables.begin();
             v != _variables.end(); ++v )
        {
            setTriVec( *v, CX );
        }

        /* build list of unassigned nodes */
        ListItem* unassignedNodes = initializeUnassigned( _nodes );


        int numberOfRandom = int( ( ( ( maxRandomRatio - minRandomRatio ) * p ) / numberOfPatterns + minRandomRatio ) *_variables.size() );
        int numberOfQuads = ( ( _variables.size() - numberOfRandom ) / 4 ) * 4;

        randomUpdates += numberOfRandom;
        parallelUpdates += numberOfQuads;
        singleUpdates += variableOrder.size() - numberOfQuads - numberOfRandom;


        std::vector<Node*>::iterator firstNonRandom = variableOrder.begin() + numberOfRandom;
        std::vector<Node*>::iterator firstNonQuad = firstNonRandom + numberOfQuads;

        std::vector<Node*>::iterator v = variableOrder.begin();

        /* assign the first "numberOfRandom" variables with random values */
        double startTimeRandom = lrabs::cpuTime();
        for( /**/; v != firstNonRandom; ++v )
        {
            if( rnd.getBool() )
            {
                setTriVec( *v, C0 );
            }
            else
            {
                setTriVec( *v, C1 );
            }
        }
        unassignedNodes = propagate01XAndUpdateList( unassignedNodes );
        timeForRandomUpdates += ( lrabs::cpuTime() - startTimeRandom );

        /* greedily select assignments (4 variables in parallel) */
        double startTimeParallel = lrabs::cpuTime();
        for( /**/; v != firstNonQuad; v += 4 )
        {
            v[0]->setTempSim( bstring0 );
            v[1]->setTempSim( bstring1 );
            v[2]->setTempSim( bstring2 );
            v[3]->setTempSim( bstring3 );

            propagate01X( unassignedNodes );
            calcScoreByNewSimParallel4( classes, scores );

            score = -1.0;
            int minIndex = -1;
            for( int i = 0; i != 16; ++i )
            {
                if( minIndex == -1 || scores[ i ] < score )
                {
                    score = scores[ i ];
                    minIndex = i;
                }
            }
            assert( minIndex != -1 );

            const TriVal t0 = value0[ minIndex ];
            const TriVal t1 = value1[ minIndex ];
            const TriVal t2 = value2[ minIndex ];
            const TriVal t3 = value3[ minIndex ];

            setTriVec( v[0], t0 );
            setTriVec( v[1], t1 );
            setTriVec( v[2], t2 );
            setTriVec( v[3], t3 );

            setAssignment( va, v[0], t0 );
            setAssignment( va, v[1], t1 );
            setAssignment( va, v[2], t2 );
            setAssignment( va, v[3], t3 );

            unassignedNodes = duplicate01XAndUpdateList( unassignedNodes, minIndex );
        }
        timeForParallelUpdates += ( lrabs::cpuTime() - startTimeParallel );

        /* greedily select assignments (remaining <=3 variables) */
        double startTimeSingle = lrabs::cpuTime();
        if( v != variableOrder.end() )
        {
            int remaining = variableOrder.size() - numberOfRandom - numberOfQuads;
            int combinations = ( 1 << remaining );

            if( remaining == 1 )
            {
                v[0]->setTempSim( bstring3 );

                propagate01X( unassignedNodes );
                calcScoreByNewSimParallel1( classes, scores );

                score = -1.0;
                int minIndex = -1;
                for( int i = 0; i != combinations; ++i )
                {
                    if( minIndex == -1 || scores[ i ] < score )
                    {
                        score = scores[ i ];
                        minIndex = i;
                    }
                }
                assert( minIndex != -1 );

                setAssignment( va, v[0], value3[ minIndex ] );
            }
            else if( remaining == 2 )
            {
                v[0]->setTempSim( bstring2 );
                v[1]->setTempSim( bstring3 );

                propagate01X( unassignedNodes );
                calcScoreByNewSimParallel2( classes, scores );

                score = -1.0;
                int minIndex = -1;
                for( int i = 0; i != combinations; ++i )
                {
                    if( minIndex == -1 || scores[ i ] < score )
                    {
                        score = scores[ i ];
                        minIndex = i;
                    }
                }
                assert( minIndex != -1 );

                setAssignment( va, v[0], value2[ minIndex ] );
                setAssignment( va, v[1], value3[ minIndex ] );
            }
            else if( remaining == 3 )
            {
                v[0]->setTempSim( bstring1 );
                v[1]->setTempSim( bstring2 );
                v[2]->setTempSim( bstring3 );

                propagate01X( unassignedNodes );
                calcScoreByNewSimParallel3( classes, scores );

                score = -1.0;
                int minIndex = -1;
                for( int i = 0; i != combinations; ++i )
                {
                    if( minIndex == -1 || scores[ i ] < score )
                    {
                        score = scores[ i ];
                        minIndex = i;
                    }
                }
                assert( minIndex != -1 );

                setAssignment( va, v[0], value1[ minIndex ] );
                setAssignment( va, v[1], value2[ minIndex ] );
                setAssignment( va, v[2], value3[ minIndex ] );
            }
            else
            {
                NEVER_GET_HERE;
            }
        }

        timeForSingleUpdates += ( lrabs::cpuTime() - startTimeSingle );

        deleteList( unassignedNodes );
        updateClasses( classes, nodeCount() );


        /* print pattern statistics */
        std::cout << "bscore_" << minRandomRatio << "_" << maxRandomRatio << " (" << p << ") = " << score << std::endl;
    }

    deleteClasses( classes );
    clearTempSim( _nodes );

    for( Node* n = _nodes; n != 0; n = n->next() )
    {
        n->unsetFlag<256|512|1024>();
    }
    unlockFlags<256|512|1024>();

    std::cout << "random_time   = " << timeForRandomUpdates
              << " / " << randomUpdates
              << " = " << timeForRandomUpdates/randomUpdates << std::endl
              << "parallel_time = " << timeForParallelUpdates
              << " / " << parallelUpdates
              << " = " << timeForParallelUpdates/parallelUpdates << std::endl
              << "single_time   = " << timeForSingleUpdates
              << " / " << singleUpdates
              << " = " << timeForSingleUpdates/singleUpdates << std::endl;

    std::cout << "time = " << lrabs::cpuTime() - startTime << std::endl;

    return score;
}

double
aigpp::Manager::
generateRandomSim( int numberOfPatterns )
{
    // std::cout << __func__ << std::endl;
    double startTime = lrabs::cpuTime();

    /* temp sim is not valid after this function! */
    _tempSimUpToDate = false;

    lockFlags<256|512|1024>();

    lrabs::Random rnd( /*seed = */42 );

    Classes classes;
    initializeClasses(  classes, _nodes, nodeCount() );

    double score = calcScore( classes );
    std::cout << "rscore (initial) = " << score << std::endl;

    for( int p = 0; p != numberOfPatterns; ++p )
    {
        VariableAssignment va( _variables.size() );

        for( std::vector<Node*>::iterator v = _variables.begin();
             v != _variables.end(); ++v )
        {
            if( rnd.getBool() )
            {
                setTriVec( *v, C0 );
                va.setNegative( ( *v )->varIndex() );
            }
            else
            {
                setTriVec( *v, C1 );
                va.setPositive( ( *v )->varIndex() );
            }
        }

        propagate01X( _nodes );
        updateClasses( classes, nodeCount() );
        score = calcScore( classes );

        std::cout << "rscore (" << p << ") = " << score << std::endl;
    }

    deleteClasses( classes );
    clearTempSim( _nodes );

    for( Node* n = _nodes; n != 0; n = n->next() )
    {
        n->unsetFlag<256|512|1024>();
    }
    unlockFlags<256|512|1024>();

    std::cout << "time = " << lrabs::cpuTime() - startTime << std::endl;

    return score;
}

double
aigpp::Manager::
generateGreedyRandomSim( int numberOfPatterns )
{
    // std::cout << __func__ << std::endl;
    double startTime = lrabs::cpuTime();

    /* temp sim is not valid after this function! */
    _tempSimUpToDate = false;

    lrabs::Random rnd( /*seed = */42 );

    Classes classes;
    initializeClasses( classes, _nodes, nodeCount() );

    double score = calcScore( classes );
    std::cout << "grscore (initial) = " << score << std::endl;

    const int combinations = 32;

    double scores[combinations];

    VariableAssignment va( _variables.size() );

    for( int p = 0; p != numberOfPatterns; ++p )
    {
        /* create "#combinations" random patterns */
        for( std::vector<Node*>::iterator v = _variables.begin();
             v != _variables.end(); ++v )
        {
            for( int i = 0; i != combinations; ++i )
            {
                set01( *v, i, rnd.getBool() );
            }
        }

        propagate01( _nodes );

        /* greedily select best random pattern */

        calcScoreByNewSimParallel01_32( classes, scores );

        double oldScore = score;

        score = -1.0;
        int minIndex = -1;
        for( int i = 0; i != combinations; ++i )
        {
            if( minIndex == -1 || scores[ i ] < score )
            {
                score = scores[ i ];
                minIndex = i;
            }
        }
        assert( minIndex != -1 );



        for( std::vector<Node*>::iterator v = _variables.begin();
             v != _variables.end(); ++v )
        {
            setAssignment01( va, *v, get01( *v, minIndex ) );
        }

        duplicate01( _nodes, minIndex );
        updateClasses01( classes, nodeCount() );

        assert( score == calcScore( classes ) );
        assert( score <= oldScore );

        /* print pattern statistics */
        std::cout << "grscore (" << p << ") = " << score << std::endl;
    }

    deleteClasses( classes );
    clearTempSim( _nodes );

    std::cout << "time = " << lrabs::cpuTime() - startTime << std::endl;

    return score;
}

#endif
