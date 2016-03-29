/*
 * defs.h
 *
 *  Created on: Jan 30, 2015
 *      Author: dmarce1
 */

#ifndef DEFS_H_
#define DEFS_H_

#include <array>

//#define CORRECTION_ON
//#define CORRECTION_OPTIMIZE

using real = double;
using integer = long long;

constexpr integer NDIM = 3;
constexpr integer ncrit = 8;
constexpr integer nparts = 100000;
constexpr integer nchild = 8;
constexpr integer LP = 20;
//constexpr real soft_len = 1.0e-6;

#ifdef CORRECTION_ON
constexpr integer MP = 17;
#else
constexpr integer MP = 7;
#endif

using space_vector = std::array<real,3>;

#endif /* DEFS_H_ */
