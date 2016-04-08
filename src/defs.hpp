

/*  
    Copyright (c) 2016 Dominic C. Marcello

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DEFS_H_
#define DEFS_H_

#include <array>

#define CORRECTION_ON
#define CORRECTION_OPTIMIZE

using real = double;
using integer = long long;

constexpr integer NDIM = 3;
constexpr integer ncrit = 8;
constexpr integer nparts = 100000;
constexpr integer nchild = 8;
constexpr integer LP = 20;

#ifdef CORRECTION_ON
constexpr integer MP = 17;
#else
constexpr integer MP = 7;
#endif

using space_vector = std::array<real,3>;

#endif /* DEFS_H_ */
