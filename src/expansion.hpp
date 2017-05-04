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
#include <cmath>

template<class Type>
class expansion: public std::array<Type, LP> {

public:
	expansion<Type>& operator*=(Type r) {
		for (integer i = 0; i != LP; ++i) {
			(*this)[i] *= r;
		}
		return *this;
	}
	expansion();
	Type operator ()() const;
	Type& operator ()();
	Type operator ()(int i) const;
	Type& operator ()(int i);
	Type operator ()(int i, int j) const;
	Type& operator ()(int i, int j);
	Type operator ()(int i, int j, int k) const;
	Type& operator ()(int i, int j, int k);
	expansion<Type>& operator =(const expansion<Type>& expansion);
	expansion<Type>& operator =(Type expansion);
	expansion<Type> operator<<(const space_vector<Type>& dX) const;
	void translate_to_particle(const space_vector<Type>& dX, Type& phi,
			space_vector<Type>& g) const;
	Type translate_to_particle(const space_vector<Type>& dX) const;
	expansion<Type>& operator<<=(const space_vector<Type>& dX);
	std::array<Type, LP>& operator +=(const std::array<Type, LP>& vec);
	std::array<Type, LP>& operator -=(const std::array<Type, LP>& vec);
	void compute_D(const space_vector<Type>& Y);
	void invert();
	~expansion();
	std::array<expansion<Type>, NDIM> get_derivatives() const;
};

template<class Type>
inline expansion<Type>::expansion() {
}

template<class Type>
inline Type expansion<Type>::operator ()() const {
	return (*this)[0];
}
template<class Type>
inline Type& expansion<Type>::operator ()() {
	return (*this)[0];
}

template<class Type>
inline Type expansion<Type>::operator ()(int i) const {
	return (*this)[1 + i];
}
template<class Type>
inline Type& expansion<Type>::operator ()(int i) {
	return (*this)[1 + i];
}

template<class Type>
inline Type expansion<Type>::operator ()(int i, int j) const {
	static constexpr size_t
	map2[3][3] = { {0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
	return (*this)[4 + map2[i][j]];
}
template<class Type>
inline Type& expansion<Type>::operator ()(int i, int j) {
	static constexpr size_t
	map2[3][3] = { {0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
	return (*this)[4 + map2[i][j]];
}

template<class Type>
inline Type expansion<Type>::operator ()(int i, int j, int k) const {
	static constexpr size_t
	map3[3][3][3] = { { {0, 1, 2}, {1, 3, 4}, {2, 4, 5}}, { {1, 3, 4}, {3, 6, 7}, {4, 7, 8}},
		{	{	2, 4, 5}, {4, 7, 8}, {5, 8, 9}}};

	return (*this)[10 + map3[i][j][k]];
}
template<class Type>
inline Type& expansion<Type>::operator ()(int i, int j, int k) {
	static constexpr size_t
	map3[3][3][3] = { { {0, 1, 2}, {1, 3, 4}, {2, 4, 5}}, { {1, 3, 4}, {3, 6, 7}, {4, 7, 8}},
		{	{	2, 4, 5}, {4, 7, 8}, {5, 8, 9}}};

	return (*this)[10 + map3[i][j][k]];
}

template<class Type>
inline expansion<Type>& expansion<Type>::operator =(
		const expansion<Type>& expansion) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] = expansion[i];
	}
	return *this;
}

template<class Type>
inline expansion<Type>& expansion<Type>::operator =(Type expansion) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] = expansion;
	}
	return *this;
}

template<class Type>
inline expansion<Type> expansion<Type>::operator<<(
		const space_vector<Type>& dX) const {
	expansion you = *this;
	you <<= dX;
	return you;
}

template<class Type>
inline expansion<Type>& expansion<Type>::operator<<=(
		const space_vector<Type>& dX) {
	expansion<Type>& me = *this;
	for (integer a = 0; a < 3; a++) {
		me() += me(a) * dX[a];
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			me() += me(a, b) * dX[a] * dX[b] * real(0.5);
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			me(a) += me(a, b) * dX[b];
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = 0; c < 3; c++) {
				me() += me(a, b, c) * dX[a] * dX[b] * dX[c] * (1.0 / 6.0);
			}
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = 0; c < 3; c++) {
				me(a) += me(a, b, c) * dX[b] * dX[c] * real(0.5);
			}
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = a; c < 3; c++) {
				me(a, c) += me(a, b, c) * dX[b];
			}
		}
	}
	return me;
}

template<class Type>
inline Type expansion<Type>::translate_to_particle(
		const space_vector<Type>& dX) const {
	const auto& L = *this;
	Type this_phi = L();
	for (integer a = 0; a < 3; a++) {
		this_phi += L(a) * dX[a];
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			this_phi += L(a, b) * dX[a] * dX[b] * real(0.5);
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = 0; c < 3; c++) {
				this_phi += L(a, b, c) * dX[a] * dX[b] * dX[c] * (1.0 / 6.0);
			}
		}
	}
	return this_phi;
}

template<class Type>
inline std::array<Type, LP>& expansion<Type>::operator +=(
		const std::array<Type, LP>& vec) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] += vec[i];
	}
	return *this;
}

template<class Type>
inline std::array<Type, LP>& expansion<Type>::operator -=(
		const std::array<Type, LP>& vec) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] -= vec[i];
	}
	return *this;
}

//void expansion::compute_D(const space_vector<Type>& Y) {
//}

template<class Type>
inline void expansion<Type>::invert() {
	expansion<Type>& me = *this;
	for (integer a = 0; a < 3; a++) {
		me(a) = -me(a);
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = a; b < 3; b++) {
			for (integer c = b; c < 3; c++) {
				me(a, b, c) = -me(a, b, c);
			}
		}
	}
}

template<class Type>
inline expansion<Type>::~expansion() {
}

static expansion<real> factor;

__attribute((constructor))
static void init_factors() {
	factor = 0.0;
	factor() += 1.0;
	for (integer a = 0; a < NDIM; ++a) {
		factor(a) += 1.0;
		for (integer b = 0; b < NDIM; ++b) {
			factor(a, b) += 1.0;
			for (integer c = 0; c < NDIM; ++c) {
				factor(a, b, c) += 1.0;
			}
		}
	}
}

template<class Type, class Type2>
inline void multipole_interaction(expansion<Type2>& L1,
		const multipole<Type>& M1, const multipole<Type2>& M2,
		space_vector<Type2> dX) {
	Type2 y0 = 0.0;
	for (integer d = 0; d != NDIM; ++d) {
		y0 += dX[d] * dX[d];
	}
	expansion<Type2> D;
	const Type2 r2inv = 1.0 / y0;
	const Type2 d0 = -sqrt(r2inv);
	const Type2 d1 = -d0 * r2inv;
	const Type2 d2 = -3.0 * d1 * r2inv;
	const Type2 d3 = -5.0 * d2 * r2inv;
#ifdef CORRECTION_ON
#ifndef CORRECTION_OPTIMIZE
	const Type2 d4 = -7.0 * d3 * r2inv;
#endif
#endif

	D() = 0.0;
	for (integer a = 0; a < 3; a++) {
		D(a) = 0.0;
		for (integer b = a; b < 3; b++) {
			D(a, b) = 0.0;
			for (integer c = b; c < 3; c++) {
				D(a, b, c) = 0.0;
			}
		}
	}

	D() += d0;
	for (integer a = 0; a < 3; a++) {
		D(a) += dX[a] * d1;
		D(a, a) += d1;
		D(a, a, a) += dX[a] * d2;
		for (integer b = a; b < 3; b++) {
			D(a, b) += dX[a] * dX[b] * d2;
			D(a, a, b) += dX[b] * d2;
			D(a, b, b) += dX[a] * d2;
			for (integer c = b; c < 3; c++) {
				D(a, b, c) += dX[a] * dX[b] * dX[c] * d3;
			}
		}
	}
#ifdef CORRECTION_ON
	std::array<expansion<Type2>, NDIM> D4;
	for (integer i = 0; i != NDIM; ++i) {
		D4[i] = 0.0;
	}
	for (integer j = 0; j != NDIM; ++j) {
		D4[j](j, j, j) += d2;
		for (integer k = j; k != NDIM; ++k) {
			D4[j](j, k, k) += d2;
			D4[k](j, j, k) += d2;
			for (integer l = k; l != NDIM; ++l) {
				D4[j](j, k, l) += dX[k] * dX[l] * d3;
				D4[l](j, k, l) += dX[j] * dX[k] * d3;
				D4[k](j, k, l) += dX[l] * dX[j] * d3;
			}
		}
	}
#ifndef CORRECTION_OPTIMIZE
	for (integer i = 0; i != NDIM; ++i) {
		for (integer j = 0; j != NDIM; ++j) {
			D4[i](j, j, j) += dX[i] * dX[j] * d3;
			for (integer k = j; k != NDIM; ++k) {
				D4[i](j, k, k) += dX[i] * dX[j] * d3;
				D4[i](j, j, k) += dX[i] * dX[k] * d3;
				for (integer l = k; l != NDIM; ++l) {
					D4[i](j, k, l) += dX[i] * dX[j] * dX[k] * dX[l] * d4;
				}
			}
		}
	}
#endif
#endif

	L1() += M2() * D();
	for (integer a = 0; a < 3; a++) {
		for (integer b = a; b < 3; b++) {
			L1() += M2(a, b) * D(a, b) * real(0.5) * factor(a, b);
		}
	}

	for (integer a = 0; a < 3; a++) {
		L1(a) += M2() * D(a);
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = b; c < 3; c++) {
				L1(a) += M2(c, b) * D(a, b, c) * real(0.5) * factor(c, b);
			}
		}
	}

	for (integer a = 0; a < 3; a++) {
		for (integer b = a; b < 3; b++) {
			L1(a, b) += M2() * D(a, b);
		}
	}

	for (integer a = 0; a < 3; a++) {
		for (integer b = a; b < 3; b++) {
			for (integer c = b; c < 3; c++) {
				L1(a, b, c) += M2() * D(a, b, c);
			}
		}
	}

#ifdef CORRECTION_ON
	for (integer i = 0; i != NDIM; ++i) {
		for (integer j = 0; j != NDIM; ++j) {
			for (integer k = j; k != NDIM; ++k) {
				for (integer l = k; l != NDIM; ++l) {
					L1(i) -= D4[i](j, k, l)
							* (M2(j, k, l) - M1(j, k, l) * M2() / M1())
							* (Type(1) / Type(6)) * factor(j, k, l);
				}
			}
		}
	}
#endif
}
/* namespace fmmx */
#endif /* expansion_H_ */
