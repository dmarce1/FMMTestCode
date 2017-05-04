

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


template<class Type>
class multipole: public std::array<Type,MP> {
private:
public:
	multipole();
	Type operator ()() const;
	Type& operator ()();
	Type operator ()(integer i, integer j) const;
	Type& operator ()(integer i, integer j);
#ifdef CORRECTION_ON
	Type operator ()(integer i, integer j, integer k) const;
	Type& operator ()(integer i, integer j, integer k);
#endif
	multipole<Type>& operator =(const multipole<Type>& expansion);
	multipole<Type>& operator =(Type expansion);
	multipole operator>>(const space_vector<Type>& dX) const;
	multipole<Type>& operator>>=(const space_vector<Type>& Y);
	std::array<Type, MP>& operator +=(const std::array<Type, MP>& vec);
	~multipole();
};


template<class Type>
inline multipole<Type>::multipole() {
}

template<class Type>
inline Type multipole<Type>::operator ()() const {
	return (*this)[0];
}

template<class Type>
inline Type& multipole<Type>::operator ()() {
	return (*this)[0];
}

template<class Type>
inline Type multipole<Type>::operator ()(integer i, integer j) const {
	static constexpr size_t map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	return (*this)[1 + map2[i][j]];
}

template<class Type>
inline Type& multipole<Type>::operator ()(integer i, integer j) {
	static constexpr size_t map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	return (*this)[1 + map2[i][j]];
}

#ifdef CORRECTION_ON
template<class Type>
inline Type multipole<Type>::operator ()(integer i, integer j, integer k) const {
	static constexpr size_t map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } },
			{ { 2, 4, 5 }, { 4, 7, 8 }, { 5, 8, 9 } } };
	return (*this)[7 + map3[i][j][k]];
}

template<class Type>
inline Type& multipole<Type>::operator ()(integer i, integer j, integer k) {
	static constexpr size_t map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } },
			{ { 2, 4, 5 }, { 4, 7, 8 }, { 5, 8, 9 } } };
	return (*this)[7 + map3[i][j][k]];
}
#endif

template<class Type>
inline multipole<Type>& multipole<Type>::operator =(const multipole<Type>& expansion) {
	for (integer i = 0; i < MP; i++) {
		(*this)[i] = expansion[i];
	}
	return *this;
}

template<class Type>
inline multipole<Type>& multipole<Type>::operator =(Type expansion) {
	for (integer i = 0; i < MP; i++) {
		(*this)[i] = expansion;
	}
	return *this;
}

template<class Type>
inline multipole<Type> multipole<Type>::operator>>(const space_vector<Type>& dX) const {
	multipole you = *this;
	you >>= dX;
	return you;
}

template<class Type>
inline std::array<Type, MP>& multipole<Type>::operator +=(const std::array<Type, MP>& vec) {
	for (integer i = 0; i < MP; i++) {
		(*this)[i] += vec[i];
	}
	return *this;
}

template<class Type>
inline multipole<Type>::~multipole() {
}


template<class Type>
inline multipole<Type>& multipole<Type>::operator>>=(const space_vector<Type>& Y) {
	multipole<Type>& me = *this;
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



#endif /* multipole_H_ */
