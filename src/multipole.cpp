


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


#include "multipole.hpp"

multipole::multipole() {
}

real multipole::operator ()() const {
	return (*this)[0];
}
real& multipole::operator ()() {
	return (*this)[0];
}

real multipole::operator ()(integer i, integer j) const {
	return (*this)[1 + map2[i][j]];
}

real& multipole::operator ()(integer i, integer j) {
	return (*this)[1 + map2[i][j]];
}

#ifdef CORRECTION_ON
real multipole::operator ()(integer i, integer j, integer k) const {
	return (*this)[7 + map3[i][j][k]];
}

real& multipole::operator ()(integer i, integer j, integer k) {
	return (*this)[7 + map3[i][j][k]];
}
#endif

multipole& multipole::operator =(const multipole& expansion) {
	for (integer i = 0; i < MP; i++) {
		(*this)[i] = expansion[i];
	}
	return *this;
}

multipole& multipole::operator =(real expansion) {
	for (integer i = 0; i < MP; i++) {
		(*this)[i] = expansion;
	}
	return *this;
}

multipole multipole::operator>>(const space_vector& dX) const {
	multipole you = *this;
	you >>= dX;
	return you;
}

multipole& multipole::operator>>=(const space_vector& Y) {
	multipole& me = *this;
#ifdef CORRECTION_ON
	for (integer p = 0; p < 3; p++) {
		for (integer q = p; q < 3; q++) {
			for (integer r = q; r < 3; r++) {
				me(p, q, r) += me() * Y[p] * Y[q] * Y[r];
				me(p, q, r) += Y[p] * me(r, q);
				me(p, q, r) += Y[q] * me(p, r);
				me(p, q, r) += Y[r] * me(q, p);
			}
		}
	}
#endif
	for (integer p = 0; p < 3; p++) {
		for (integer q = p; q < 3; q++) {
			me(p, q) += me() * Y[p] * Y[q];
		}
	}
	return me;
}

std::array<real, MP>& multipole::operator +=(const std::array<real, MP>& vec) {
	for (integer i = 0; i < MP; i++) {
		(*this)[i] += vec[i];
	}
	return *this;
}

multipole::~multipole() {
}
