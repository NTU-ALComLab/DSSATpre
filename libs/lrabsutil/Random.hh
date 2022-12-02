/**************************************************************
 *       
 *       LRABS // Random.hh
 *
 *       Copyright (C) 2007 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 711 $
 *         $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#ifndef LRABS_RANDOM_HH
#define LRABS_RANDOM_HH

#include <cstdlib>

namespace lrabs
{
    /*!
     * A random number generator class with internal seed.
     *
     * Random unsigned integers can be created using Random's getUInt() 
     * method, booleans are created using getBool().
     *
     * An example:
     * \code
     * 
     * #include <iostream>
     * #include <lrabsutil/Random.hh>
     * 
     * int main()
     * {
     *   // create a random number generator with seed 42
     *   lrabs::Random rng( 42 );
     * 
     *   // generate ans output 100 random numbers
     *   for( int i = 0; i < 100; ++i )
     *   {
     *     unsigned int r = rng.getUInt();
     *     std::cout << "random number #" << i << ": " << r << std::endl;
     *   }
     * 
     *   return 0;
     * }
     * \endcode
     *
     */
    class Random
    {
    public:
        /*!
         * Constructor. Creates a random number generator with specified seed.
         * The default seed is 1.
         *
         * \param seed initial seed of the random number generator
         */
        explicit Random( unsigned int seed = 1 );

        /*!
         * Get a random bool.
         *
         * \return random bool
         */
        bool getBool();

        /*!
         * Get a random unsigned int.
         *
         * \return random unsigned int
         */
        inline unsigned int getUInt();

        /*!
         * Get a random unsigned int \e u with 0 <= \e u < \e max.
         *
         * \param max
         *
         * \return random unsigned int
         */
        inline unsigned int getUInt( unsigned int max );

        /*!
         * Get a random unsigned long. 
         * This method is 64-bit aware!
         *
         * \return random unsigned long
         */
        inline unsigned long getULong();

        /*!
         * Get a random unsigned int \e u with 0 <= \e u < \e max. 
         * This method is 64-bit aware!
         *
         * \param max
         *
         * \return random unsigned long
         */
        inline unsigned long getULong( unsigned long max );

    private:
        unsigned int getRawRandom();

        unsigned int _seed;
        unsigned int _boolPool;
        unsigned int _boolPoolRemaining;
    };
}

#include "Random.icc"

#endif /* LRABS_RANDOM_HH */
