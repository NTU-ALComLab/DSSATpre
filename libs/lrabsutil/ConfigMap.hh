/**************************************************************
 *       
 *       LRABS // ConfigMap.hh
 *
 *       Copyright (C) 2006 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 236 $
 *         $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#ifndef LRABS_CONFIGMAP_HH
#define LRABS_CONFIGMAP_HH

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include "Exception.hh"
#include "String.hh"

namespace lrabs
{
  /*!
   * The ConfigMap class is basically a std::map<std::string, std::string>
   * with additional methods to load/store the maps contents from/to a
   * stream (e.g. file stream), and to manipulate key/value pairs.
   */
  class ConfigMap: public std::map<std::string, std::string>
  {
  public:
    
    /*! 
     *   Reads a ConfigMap from an input stream.
     *   Comment lines (starting with a '#') are ignored.
     *   Valid lines are of the form 'key=value' where neither key nor value 
     *   contain any whitespace characters. Arbitrary amounts of whitespace 
     *   are allowed before and after both key and value.
     *   If there is a format error, an error mnessage is printed and a 
     *   lrabs::Exception is thrown.
     *   If there are duplicate keys, a warning message is printed.
     *
     *   An example of a valid input file:
     *   \verbatim
# this is a comment

# a string
name = Bob

# phone
phone = 1234567

# some boolean flag
flag = true
\endverbatim
     *
     *   \param is input stream
     *
     *   \exception lrabs::Excetion on format error
     */
    void readMap( std::istream& is );

    /*! 
     *   Writes the ConfigMap to an output stream.
     *
     *   \param os output stream
     */
    void writeMap( std::ostream& os ) const ;

    /*! 
     *   Inserts a key/value pair to the ConfigMap. The value is converted
     *   to a string by using an std::ostringstream (and T's operator<<).
     *   The conversion result string may not contain any whitespace!
     *   If the map already contains the key, a warning message is printed.
     *
     *   \param key 
     *   \param value
     */
    template<typename T>
    void setValue( const std::string& key, const T& value );

    /*!
     *   Reads a value pair from the ConfigMap. Internally all values are stored as
     *   std::string object. Using this method you can explicitly convert the
     *   value to an arbitrary (but matching) target type T. 
     *   The target type T must have an input operator (operator>>).
     *
     *   If the map does not contain the key, an error message is printed and
     *   an empty exception is thrown.
     *   If the conversion  does not succedd, an error message is printed and
     *   a lrabs::Excetion is thrown. 
     *
     *   \param key 
     *
     *   \exception lrabs::Excetion if key is not found and if the conversion fails
     *   
     *   \return value
     */
    template<typename T>
    T getValue( const std::string& key ) const;
    
    /*!
     *   Reads a value pair from the ConfigMap. Internally all values are stored as
     *   std::string object. Using this method you can explicitly convert the
     *   value to an arbitrary (but matching) target type T. 
     *   The target type T must have an input operator (operator>>).
     *
     *   If the map does not contain the key, an error message is printed and
     *   an empty exception is thrown.
     *   If the conversion  does not succedd, an error message is printed and
     *   a lrabs::Excetion is thrown. 
     *
     *   \param key 
     *   \param mapping
     *
     *   \exception lrabs::Excetion if key is not found or if the value is not found in the mapping
     *   
     *   \return value
     */
    template<typename T>
    T getValue( const std::string& key, const std::map<std::string, T>& mapping ) const;
    
    /*! 
     *   Reads a value pair from the ConfigMap and stores it into target. 
     *   Internally all values are stored as std::string object. Using this 
     *   method you can explicitly convert the
     *   value to an arbitrary (but matching) target type T. 
     *   The target type T must have an input operator (operator>>).
     *
     *   If the map does not contain the key, no value is loaded into the 
     *   target and false is returned.
     *   If the conversion  does not succedd, an error message is printed and
     *   a lrabs::Excetion is thrown. 
     *
     *   \param key 
     *   \param target
     *
     *   \exception lrabs::Excetion if the conversion fails
     *   
     *   \return true if the key was found
     */
    template<typename T>
    bool getValue( const std::string& key, T& target ) const;
    
    template<typename T>
    bool getValue( const std::string& key, const std::map<std::string, T>& mapping, T& target ) const;
    
    /*! 
     *   Checks if the ConfigMap contains the key.
     *
     *   \param key 
     *
     *   \return true if the key is in the ConfigMap
     */
    bool hasKey( const std::string& key ) const;
    
    /*! 
     *   Checks if all keys in the ConfigMap are also in the validKeys set.
     *
     *   \param validKeys
     *
     *   \return true if all keys in the ConfigMap are also in the validKeys set
     */
    bool allKeysAreValid( const std::set<std::string>& validKeys ) const;
    
    /*! 
     *   Checks if all keys from the keys set are also in the ConfigMap.
     *
     *   \param keys
     *
     *   \return true if all keys from the keys set are also in the ConfigMap
     */
    bool containsAllKeys( const std::set<std::string>& keys ) const;
  };
  
}

#include "ConfigMap.icc"

#endif /* LRABS_CONFIGMAP_HH */
