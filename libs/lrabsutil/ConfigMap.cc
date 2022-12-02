/**************************************************************
*       
*       LRABS // ConfigMap.cc
*
*       Copyright (C) 2006 Florian Pigorsch
*
*       Author:
*         Florian Pigorsch
*         University of Freiburg
*         pigorsch@informatik.uni-freiburg.de
*
*       Last revision:
*         $Revision: 585 $
*         $Author: pigorsch $
*         $Date$
*
***************************************************************/

#include "ConfigMap.hh"

#include <cctype>
#include <cstdio>

void
lrabs::ConfigMap::
readMap( std::istream& is )
{
  std::string s;
  
  while( is )
  {
    if( is.peek() == EOF )
    {
      break;
    }
    else if( isspace( is.peek() ) ) 
    {
      is.ignore();
    }
    else if( is.peek() == '#' ) 
    {
      getline( is, s );
    }
    else
    {
      getline( is, s );

      std::string key, value;
      
      unsigned int pos = 0;
      
      // skip initial spaces
      while( pos < s.size() && isspace( s[ pos ] ) ) ++pos;
      
      // read key
      while( pos < s.size() && !isspace( s[ pos ] ) && !( s[ pos ] == '=' ) ) 
      {
        key += s[ pos ];
        ++pos;
      }
      if( key.size() == 0 )
      {
        Exception e;
        e << "ERROR: readMap, error parsing key of line \"" << s << "\"";
        throw e;
      }
      
      // skip spaces after key
      while( pos < s.size() && isspace( s[ pos ] ) ) ++pos;
      
      // skip '='
      if( pos >= s.size() || s[ pos ] != '=' )
      {
        Exception e;
        e << "ERROR: readMap, error parsing '=' of line \"" << s << "\"";
        throw e;
      }
      else
      {
        ++pos;
      }
      
      
      // skip spaces after '='
      while( pos < s.size() && isspace( s[ pos ] ) ) ++pos;
      
      // read value
      while( pos < s.size() && !isspace( s[ pos ] ) ) 
      {
        value += s[ pos ];
        ++pos;
      }
      
      // skip spaces after value
      while( pos < s.size() && isspace( s[ pos ] ) ) ++pos;

      // there are some nonspace characters after the value
      if( pos < s.size() )
      {
        Exception e;
        e << "ERROR: readMap, error parsing value of line \"" << s << "\"";
        throw e;
      }
        
      if( hasKey( key ) )
      {
        Exception e;
        e << "WARNING: readMap, duplicate key \"" << key << "\" in map. overwriting...";
        throw e;
      }

      setValue( key, value );
    }
  }
}

bool 
lrabs::ConfigMap::
allKeysAreValid( const std::set<std::string>& validKeys ) const
{
  bool valid = true;
  
  for( const_iterator p = begin(); p != end(); ++p )
  {
    if( validKeys.find( p->first ) == validKeys.end() )
    {
      std::cout << "INFO: allKeysAreValid, key \"" << p->first << "\" is not valid" << std::endl;
      
      valid = false;
    }
  }
  
  return valid;
}

bool 
lrabs::ConfigMap::
containsAllKeys( const std::set<std::string>& keys ) const
{
  bool all = true;
  
  for( std::set<std::string>::const_iterator p = keys.begin(); p != keys.end(); ++p )
  {
    if( !hasKey( *p ) )
    {
      std::cout << "INFO: containsAllKeys, key \"" << *p << "\" is not in the map" << std::endl;
      all = false;
    }
  }
  
  return all;
}

void 
lrabs::ConfigMap::
writeMap( std::ostream& os ) const
{
  for( const_iterator p = begin(); p != end(); ++p )
  {
    os << p->first << "=" << p->second << std::endl;
  }
}

bool 
lrabs::ConfigMap::
hasKey( const std::string& key ) const
{
  return( find( key ) != end() );
}

