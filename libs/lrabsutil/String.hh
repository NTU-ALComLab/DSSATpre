/**************************************************************
*       
*       LRABS Utilities // String.hh
*
*       Copyright (C) 2007 Florian Pigorsch
*
*       Author:
*         Florian Pigorsch
*         University of Freiburg
*         pigorsch@informatik.uni-freiburg.de
*
*       Last revision:
*         $Revision: 605 $
*         $Author: pigorsch $
*
***************************************************************/

#ifndef LRABS_STRING_HH
#define LRABS_STRING_HH

#include <sstream>
#include <string>

#include "Exception.hh"

namespace lrabs
{
    /*! 
     *   Converts a value of type T to a std::string. The conversion is done 
     *   by a std::ostringstream.
     *
     *   \param value
     *   \return string representation
     */
    template<typename T>
    std::string toString( const T& value );
    
    /*! 
     *   Converts a std::string to a value of type T. The conversion is done 
     *   by a std::istringstream.
     * 
     *   If the conversion is not successful, an Exception is thrown.
     *
     *   \param s
     * 
     *   \return T representation
     */
    template<typename T>
    T fromString( const std::string& s );


    /*!
     *  Removes leading and trailing whitespace from the std::string s.
     *  Whitespace is determined by the function 'isspace'.
     *
     *  \param s
     */
    void trim( std::string& s );
    
    /*!
     *  Removes leading and trailing whitespace from the std::string s and returnes
     *  the trimmed string.
     *  Whitespace is determined by the function 'isspace'.
     *
     *  \param s
     *  \return trimmed string
     */
    std::string trimCopy( const std::string& s );

    /*!
     *  Split the string 'input' into several substrings using the non-empty delimiter 'delimiter'.
     *  The individual substrings are stored in the container 'result'. The function returns the
     *  number of found substrings. The container 'result' is cleared (using its method 'clear()')
     *  before the substrings are inserted using its method 'push_back(...)'.
     *
     *  \param input input string
     *  \param delimiter delimiter string (must be non-empty)
     *  \param result
     *
     *  \return number of substrings
     */
    template<typename STRING_CONTAINER>
    int split( const std::string& input, const std::string& delimiter, STRING_CONTAINER& result );

    /*!
     *  Returns true if s1 and s2 are equal modulo case. Special characters (german umlauts, etc.)
     *  are not respected, only a-z
     */
    bool equalCaseInsensitive( std::string s1, std::string s2 );
}

#include "String.icc"

#endif /* LRABS_STRING_HH */
