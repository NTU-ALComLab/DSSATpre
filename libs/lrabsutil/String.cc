/**************************************************************
*       
*       LRABS // String.cc
*
*       Copyright (C) 2009 Florian Pigorsch
*
*       Author:
*         Florian Pigorsch
*         University of Freiburg
*         pigorsch@informatik.uni-freiburg.de
*
*       Last revision:
*         $Revision: 716 $
*         $Author: pigorsch $
*         $Date$
*
***************************************************************/

#include "String.hh"

#include <cctype>

void
lrabs::
trim( std::string& s )
{
    std::size_t firstNonWS = 0;
    std::size_t lastNonWS = std::string::npos;

    bool justWS = true;
    
    for( std::size_t i = 0 ; i != s.size(); ++i )
    {
        if( isspace( s[i] ) )
        {
            if( justWS )
            {
                ++firstNonWS;
            }
        }
        else
        {
            justWS = false;
            lastNonWS = i;
        }
    }
    
    if( firstNonWS == s.size() )
    {
        s.clear();
        return;
    }

    if( ( lastNonWS + 1 ) != s.size() )
    {
        /* there is trailing space */
        s.erase( lastNonWS + 1 );
    }

    if( firstNonWS != 0 )
    {
        s.erase( 0, firstNonWS );
    }
}


    
        
std::string
lrabs::
trimCopy( const std::string& s )
{
    std::size_t firstNonWS = 0;
    std::size_t lastNonWS = std::string::npos;

    bool justWS = true;
    
    for( std::size_t i = 0 ; i != s.size(); ++i )
    {
        if( isspace( s[i] ) )
        {
            if( justWS )
            {
                ++firstNonWS;
            }
        }
        else
        {
            justWS = false;
            lastNonWS = i;
        }
    }
    
    if( firstNonWS == s.size() )
    {
        return std::string();
    }

    return s.substr( firstNonWS, 1 + lastNonWS - firstNonWS );
}

bool 
lrabs::
equalCaseInsensitive( std::string s1, std::string s2 )
{
    if( s1.size() != s2.size() ) return false;

    for( std::size_t i = 0; i != s1.size(); ++i )
    {
        char c1 = ( ( s1[i] >= 'A' && s1[i] <= 'Z' ) ? ( s1[i] - 'A' + 'a'  ) : s1[i] );
        char c2 = ( ( s2[i] >= 'A' && s2[i] <= 'Z' ) ? ( s2[i] - 'A' + 'a'  ) : s2[i] );

        if( c1 != c2 ) return false;
    }

    return true;
}
