

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

#ifndef EXPAN222SION_H_
#define EXPAN222SION_H_

#include "defs.hpp"
#include "multipole.hpp"

class expansion: public std::array<real, LP> {
private:
	const real delta[3][3] = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
	const size_t map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	const size_t map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } },
			{ { 2, 4, 5 }, { 4, 7, 8 }, { 5, 8, 9 } } };
public:
	expansion& operator*=(real r) {
		for (integer i = 0; i != LP; ++i) {
			(*this)[i] *= r;
		}
		return *this;
	}
	expansion();
	real operator ()() const;
	real& operator ()();
	real operator ()(int i) const;
	real& operator ()(int i);
	real operator ()(int i, int j) const;
	real& operator ()(int i, int j);
	real operator ()(int i, int j, int k) const;
	real& operator ()(int i, int j, int k);
	expansion& operator =(const expansion& expansion);
	expansion& operator =(real expansion);
	expansion operator<<(const space_vector& dX) const;
	void translate_to_particle(const space_vector& dX, real& phi, space_vector& g) const;
	real translate_to_particle(const space_vector& dX) const;
	expansion& operator<<=(const space_vector& dX);
	std::array<real, LP>& operator +=(const std::array<real, LP>& vec);
	std::array<real, LP>& operator -=(const std::array<real, LP>& vec);
	void compute_D(const space_vector& Y);
	void invert();
	~expansion();
	std::array<expansion, NDIM> get_derivatives() const;
};

void multipole_interaction(expansion&, expansion&, const multipole&, const multipole&, space_vector dX);

/* namespace fmmx */
#endif /* expansion_H_ */
