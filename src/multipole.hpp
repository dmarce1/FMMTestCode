

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


#ifndef multipole_H_
#define multipole_H_

#include "defs.hpp"
#include <array>

class multipole: public std::array<real,MP> {
private:
	const real delta[3][3] = { { real(1.0), real(0.0), real(0.0) }, { real(0.0), real(1.0), real(0.0) }, { real(0.0),
			real(0.0), real(1.0) } };
	const size_t map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
#ifdef CORRECTION_ON
	const size_t map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } },
			{ { 2, 4, 5 }, { 4, 7, 8 }, { 5, 8, 9 } } };
#endif
public:
	multipole();
	real operator ()() const;
	real& operator ()();
	real operator ()(integer i, integer j) const;
	real& operator ()(integer i, integer j);
#ifdef CORRECTION_ON
	real operator ()(integer i, integer j, integer k) const;
	real& operator ()(integer i, integer j, integer k);
#endif
	multipole& operator =(const multipole& expansion);
	multipole& operator =(real expansion);
	multipole operator>>(const space_vector& dX) const;
	multipole& operator>>=(const space_vector& Y);
	std::array<real, MP>& operator +=(const std::array<real, MP>& vec);
	~multipole();
};

#endif /* multipole_H_ */
