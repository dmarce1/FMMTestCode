

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


#include "expansion.hpp"
#include <cmath>

expansion::expansion() {
}

real expansion::operator ()() const {
	return (*this)[0];
}
real& expansion::operator ()() {
	return (*this)[0];
}

real expansion::operator ()(int i) const {
	return (*this)[1 + i];
}
real& expansion::operator ()(int i) {
	return (*this)[1 + i];
}

real expansion::operator ()(int i, int j) const {
	return (*this)[4 + map2[i][j]];
}
real& expansion::operator ()(int i, int j) {
	return (*this)[4 + map2[i][j]];
}

real expansion::operator ()(int i, int j, int k) const {
	return (*this)[10 + map3[i][j][k]];
}
real& expansion::operator ()(int i, int j, int k) {
	return (*this)[10 + map3[i][j][k]];
}

expansion& expansion::operator =(const expansion& expansion) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] = expansion[i];
	}
	return *this;
}

expansion& expansion::operator =(real expansion) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] = expansion;
	}
	return *this;
}

expansion expansion::operator<<(const space_vector& dX) const {
	expansion you = *this;
	you <<= dX;
	return you;
}

expansion& expansion::operator<<=(const space_vector& dX) {
	expansion& me = *this;
	for (integer a = 0; a < 3; a++) {
		me() += me(a) * dX[a];
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			me() += me(a, b) * dX[a] * dX[b] * 0.5;
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
				me(a) += me(a, b, c) * dX[b] * dX[c] * 0.5;
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

real expansion::translate_to_particle(const space_vector& dX) const {
	const auto& L = *this;
	real this_phi = L();
	for (integer a = 0; a < 3; a++) {
		this_phi += L(a) * dX[a];
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			this_phi += L(a, b) * dX[a] * dX[b] * 0.5;
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

std::array<real, LP>& expansion::operator +=(const std::array<real, LP>& vec) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] += vec[i];
	}
	return *this;
}

std::array<real, LP>& expansion::operator -=(const std::array<real, LP>& vec) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] -= vec[i];
	}
	return *this;
}

//void expansion::compute_D(const space_vector& Y) {
//}

void expansion::invert() {
	expansion& me = *this;
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

expansion::~expansion() {
}

void multipole_interaction(expansion& L1, expansion& L2, const multipole& M1, const multipole& M2, space_vector dX) {

	const real delta[3][3] = { { real(1.0), real(0.0), real(0.0) }, { real(0.0), real(1.0), real(0.0) }, { real(0.0),
			real(0.0), real(1.0) } };

	real y0 = 0.0;
	for (integer d = 0; d != NDIM; ++d) {
		y0 += dX[d] * dX[d];
	}
	expansion D;
	const real r2inv = 1.0 / y0;
	const real d0 = -std::sqrt(r2inv);
	const real d1 = -d0 * r2inv;
	const real d2 = -3.0 * d1 * r2inv;
	const real d3 = -5.0 * d2 * r2inv;
#ifdef CORRECTION_ON
#ifndef CORRECTION_OPTIMIZE
	const real d4 = -7.0 * d3 * r2inv;
#endif
#endif
	D() = d0;
	for (integer a = 0; a < 3; a++) {
		D(a) = dX[a] * d1;
		for (integer b = a; b < 3; b++) {
			D(a, b) = dX[a] * dX[b] * d2 + delta[a][b] * d1;
			for (integer c = b; c < 3; c++) {
				D(a, b, c) = dX[a] * dX[b] * dX[c] * d3;
				D(a, b, c) += (delta[a][b] * dX[c] + delta[b][c] * dX[a] + delta[c][a] * dX[b]) * d2;
			}
		}
	}
#ifdef CORRECTION_ON
	std::array<expansion, NDIM> D4;
	for (integer i = 0; i != NDIM; ++i) {
		D4[i] = 0.0;
		for (integer j = 0; j != NDIM; ++j) {
			for (integer k = j; k != NDIM; ++k) {
				for (integer l = k; l != NDIM; ++l) {
					D4[i](j, k, l) += delta[i][j] * delta[k][l] * d2;
					D4[i](j, k, l) += delta[i][l] * delta[j][k] * d2;
					D4[i](j, k, l) += delta[i][k] * delta[l][j] * d2;

					D4[i](j, k, l) += delta[i][j] * dX[k] * dX[l] * d3;
					D4[i](j, k, l) += delta[i][l] * dX[j] * dX[k] * d3;
					D4[i](j, k, l) += delta[i][k] * dX[l] * dX[j] * d3;
#ifndef CORRECTION_OPTIMIZE
					D4[i](j, k, l) += dX[i] * dX[j] * delta[k][l] * d3;
					D4[i](j, k, l) += dX[i] * dX[l] * delta[j][k] * d3;
					D4[i](j, k, l) += dX[i] * dX[k] * delta[l][j] * d3;

					D4[i](j, k, l) += dX[i] * dX[j] * dX[k] * dX[l] * d4;
#endif
				}
			}
		}
	}
#endif

	L1() += M2() * D();
	L2() += M1() * D();
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			L1() += M2(a, b) * D(a, b) * 0.5;
			L2() += M1(a, b) * D(a, b) * 0.5;
		}
	}

	for (integer a = 0; a < 3; a++) {
		L1(a) += M2() * D(a);
		L2(a) -= M1() * D(a);
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = 0; c < 3; c++) {
				L1(a) += M2(c, b) * D(a, b, c) * 0.5;
				L2(a) -= M1(c, b) * D(a, b, c) * 0.5;
			}
		}
	}

	for (integer a = 0; a < 3; a++) {
		for (integer b = a; b < 3; b++) {
			L1(a, b) += M2() * D(a, b);
			L2(a, b) += M1() * D(a, b);
		}
	}

	for (integer a = 0; a < 3; a++) {
		for (integer b = a; b < 3; b++) {
			for (integer c = b; c < 3; c++) {
				L1(a, b, c) += M2() * D(a, b, c);
				L2(a, b, c) -= M1() * D(a, b, c);
			}
		}
	}



#ifdef CORRECTION_ON
	for (integer i = 0; i != NDIM; ++i) {
		for (integer j = 0; j != NDIM; ++j) {
			for (integer k = 0; k != NDIM; ++k) {
				for (integer l = 0; l != NDIM; ++l) {
					L1(i) -= D4[i](j, k, l) * (M2(j, k, l) - M1(j, k, l) * M2() / M1()) * (real(1) / real(6));
					L2(i) -= D4[i](j, k, l) * (M1(j, k, l) - M2(j, k, l) * M1() / M2()) * (real(1) / real(6));
				}
			}
		}
	}
#endif
}

/*
 std::array<expansion, NDIM> expansion::get_derivatives() const {
 std::array<expansion, NDIM> D;
 for (integer i = 0; i != NDIM; ++i) {
 D[i]() = (*this)(i);
 for (integer a = 0; a < NDIM; ++a) {
 D[i](a) = (*this)(i, a);
 for (integer b = a; b < NDIM; ++b) {
 D[i](a, b) = (*this)(i, a, b);
 for (integer c = b; c != NDIM; ++c) {
 D[i](a, b, c) = real(0.0);
 }
 }
 }
 }
 #ifdef CORRECTION_ON
 const real rinv = real(1) / r;
 const real D2 = -real(3.0) * std::pow(rinv, real(5));
 const real D3 = +real(15.) * std::pow(rinv, real(7));
 const real D4 = -real(105) * std::pow(rinv, real(9));
 std::array<expansion, NDIM> E1;
 for (integer i = 0; i != NDIM; ++i) {
 E1[i] = 0.0;
 for (integer j = 0; j != NDIM; ++j) {
 for (integer k = j; k != NDIM; ++k) {
 for (integer l = k; l != NDIM; ++l) {
 E1[i](j, k, l) += delta[i][j] * delta[k][l] * D2;
 E1[i](j, k, l) += delta[i][l] * delta[j][k] * D2;
 E1[i](j, k, l) += delta[i][k] * delta[l][j] * D2;

 E1[i](j, k, l) += delta[i][j] * R[k] * R[l] * D3;
 E1[i](j, k, l) += delta[i][l] * R[j] * R[k] * D3;
 E1[i](j, k, l) += delta[i][k] * R[l] * R[j] * D3;

 E1[i](j, k, l) += R[i] * R[j] * delta[k][l] * D3;
 E1[i](j, k, l) += R[i] * R[l] * delta[j][k] * D3;
 E1[i](j, k, l) += R[i] * R[k] * delta[l][j] * D3;

 E1[i](j, k, l) += R[i] * R[j] * R[k] * R[l] * D4;


 }
 }
 }
 }
 for (integer i = 0; i != NDIM; ++i) {
 for (integer j = 0; j != NDIM; ++j) {
 for (integer k = 0; k != NDIM; ++k) {
 for (integer l = 0; l != NDIM; ++l) {
 D[i]() -= E1[i](j, k, l) * M(j, k, l) * (real(1) / real(6));
 D[i]() += E1[i](j, k, l) * N(j, k, l) * M() / N() * (real(1) / real(6));
 }
 }
 }
 }
 #endif
 return D;
 }
 */

