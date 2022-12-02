/**************************************************************
*       
*       LRABS // Exception.hh
*
*       Copyright (C) 2006 Florian Pigorsch
*
*       Author:
*         Florian Pigorsch
*         University of Freiburg
*         pigorsch@informatik.uni-freiburg.de
*
*       Last revision:
*         $Revision: 141 $
*         $Author: pigorsch $
*         $Date$
*
***************************************************************/

#ifndef LRABS_EXCEPTION_HH
#define LRABS_EXCEPTION_HH

// std includes
#include <sstream>
#include <string>
#include <exception>

namespace lrabs
{
  /*!
   * An advanced exception class subclassing std::exception.
   * The lrabs::Exception class offers an operator<< which can be used to
   * specify an error message in a convenient way:
   * \code
   * ...
   * if( x <= 0.0 )
   * {
   *   lrabs::Exception e;
   *   e << "An error occured: x is negative (x=" << x << ")!";
   *   throw e;
   * }
   * ...
   * \endcode
   *
   * The error message of an lrabs::Exception can be accessed by its what()
   * method:
   * \code
   * try
   * {
   *   ...
   * }
   * catch( const lrabs::Exception& e )
   * {
   *   std::cerr << "Caught an exception! Its error message is" << std::endl
   *             << e.what() << std::endl;
   * }
   * \endcode
   */
  class Exception: public std::exception
  {
  public:
    /*!
     * Constructor. Constructs an empty Exception.
     */
    Exception();

    /*!
     * Copy constructor. Constructs a copy of e.
     *
     * \param e
     */
    Exception( const Exception& e );

    /*!
     * Move constructor.
     */

    Exception(Exception&& e) noexcept;

    /*!
     * Destructor.
     */
    ~Exception() noexcept override;

    Exception& operator=(const Exception& other) = delete;
    Exception& operator=(Exception&& other) = delete;

    /*!
     * Output operator for Exception objects. Use this operator to modify the Exception's
     * error message. Any object of type T can be given to this operator, as long
     * as there exists an std::ostream output operator for T.
     *
     * \param data object to be appended to the error message
     * 
     * \return reference to the modified object
     */
    template<typename T>
    Exception& operator<<( const T& data );

    /*!
     * Retrieve the Exception's error message.
     *
     * \return error message
     */
    const char* what() const noexcept override;

  protected:
    std::string _message;
  };

}

#include "Exception.icc"

#endif /* LRABS_EXCEPTION_HH */
