/**************************************************************
 *       
 *       LRABS // Hashes.hh
 *       
 *       Copyright (C) 2007 Florian Pigorsch
 *       
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *       
 *       Last revision:
 *         $Revision: 442 $
 *         $Author: pigorsch $
 *       
 ***************************************************************/

#ifndef LRABS_HASHES_HH
#define LRABS_HASHES_HH

#include <cstddef>

namespace lrabs
{
    template<typename T>
    struct hashPointer
    {
        std::size_t operator()( T* p ) const
            {
                return (std::size_t)(p) >> 3u;
            }
    };
}

#endif /* LRABS_HASHES_HH */
