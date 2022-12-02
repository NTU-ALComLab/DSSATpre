/**************************************************************
*       
*       LRABS // Exception.cc
*
*       Copyright (C) 2006 Florian Pigorsch
*
*       Author:
*         Florian Pigorsch
*         University of Freiburg
*         pigorsch@informatik.uni-freiburg.de
*
*       Last revision:
*         $Revision: 128 $
*         $Author: pigorsch $
*         $Date$
*
***************************************************************/

#include "Exception.hh"

lrabs::Exception::
Exception()
        :_message()
{}

lrabs::Exception::
Exception( const lrabs::Exception& e )
        :std::exception( e ),
         _message( e._message )
{}

lrabs::Exception::
Exception( lrabs::Exception&& e ) noexcept
        :std::exception( std::move(e) ),
         _message( std::move(e._message) )
{}


lrabs::Exception::
~Exception() noexcept
{}
