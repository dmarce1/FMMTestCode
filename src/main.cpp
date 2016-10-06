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

#include "defs.hpp"
#include "multipole.hpp"
#include "expansion.hpp"
#include "simd.hpp"

#include <cassert>
#include <memory>
#include <unistd.h>
#include <set>
#include <cmath>
#include <list>
#include <vector>

real theta = 1.0;
integer ncells = 0;
bool analytic_computed = false;
static std::multiset<real> top10;

struct particle {
	space_vector<real> X;
	space_vector<real> g;
	space_vector<real> g_analytic;
	real M;
	real phi;
	real phi_analytic;
	bool operator<(const particle& other) {
		if (X[0] < other.X[0]) {
			return true;
		} else if (X[0] == other.X[0]) {
			if (X[1] < other.X[1]) {
				return true;
			} else if (X[1] == other.X[1]) {
				if (X[2] < other.X[2]) {
					return true;
				}
			}
		}
		return false;
	}

};

class cell {
	space_vector<real> Xcom;
	space_vector<real> Xmax;
	space_vector<real> Xmin;
	real Rmax;
	expansion<real> Lphi;
	multipole<real> M;
	integer level;
	bool is_leaf;
	std::vector<particle> particles;
	cell* children[nchild];

public:

	cell() {
		ncells++;
	}
	~cell() {
		if (!is_leaf) {
			for (integer i = 0; i != nchild; ++i) {
				delete children[i];
			}
		}
		--ncells;
	}
	bool operator<(const cell& other) {
		if (Xcom[0] < other.Xcom[0]) {
			return true;
		} else if (Xcom[0] == other.Xcom[0]) {
			if (Xcom[1] < other.Xcom[1]) {
				return true;
			} else if (Xcom[1] == other.Xcom[1]) {
				if (Xcom[2] < other.Xcom[2]) {
					return true;
				}
			}
		}
		return false;
	}

	space_vector<real> force_sum(bool norm) const {
		space_vector < real > sum;
		for (integer d = 0; d != NDIM; ++d) {
			sum[d] = 0.0;
		}
		if (is_leaf) {
			for (auto l = particles.begin(); l != particles.end(); ++l) {
				for (integer d = 0; d != NDIM; ++d) {
					if (norm) {
						sum[d] += std::abs(l->g[d] * l->M);
					} else {
						sum[d] += l->g[d] * l->M;
					}
				}
			}
		} else {
			for (integer c = 0; c != nchild; ++c) {
				auto v = children[c]->force_sum(norm);
				for (integer d = 0; d != NDIM; ++d) {
					sum[d] += v[d];
				}
			}
		}
		return sum;
	}

	space_vector<real> torque_sum(const bool norm) const {
		space_vector < real > sum;
		for (integer d = 0; d != NDIM; ++d) {
			sum[d] = 0.0;
		}
		auto f = [=](real a) {return norm ? std::abs(a) : a;};
		if (is_leaf) {
			for (auto l = particles.begin(); l != particles.end(); ++l) {
				sum[0] += f(
						+l->X[1] * l->g[2] * l->M - l->X[2] * l->g[1] * l->M);
				sum[1] += f(
						-l->X[0] * l->g[2] * l->M + l->X[2] * l->g[0] * l->M);
				sum[2] += f(
						+l->X[0] * l->g[1] * l->M - l->X[1] * l->g[0] * l->M);
			}
		} else {
			for (integer c = 0; c != nchild; ++c) {
				auto v = children[c]->torque_sum(norm);
				for (integer d = 0; d != NDIM; ++d) {
					sum[d] += v[d];
				}
			}
		}
		return sum;
	}
	void output(const char* filename) const {
		if (is_leaf) {
			FILE* fp = fopen(filename, "at");
			for (const auto& p : particles) {
				fprintf(fp, "%e %e %e\n", p.X[0], p.X[1], p.X[2]);
			}
			fclose(fp);
		} else {
			for (integer c = 0; c != nchild; ++c) {
				children[c]->output(filename);
			}
		}
	}

	bool is_well_separated_from(const cell& other) const {
		real Z = 0.0;
		for (integer d = 0; d != NDIM; ++d) {
			Z += std::pow(other.Xcom[d] - Xcom[d], real(2));
		}
		Z = sqrt(Z);
		assert(Rmax >= 0.0);
		assert(other.Rmax >= 0.0);
		return bool(Z > (Rmax + other.Rmax) / theta);
	}

	void compute_expansions(const expansion<real>& Lphi0,
			const space_vector<real>& X0) {
		space_vector < real > dX;
		for (integer i = 0; i != NDIM; ++i) {
			dX[i] = Xcom[i] - X0[i];
		}
		Lphi += Lphi0 << dX;
		if (is_leaf) {
			for (auto l = particles.begin(); l != particles.end(); ++l) {
				for (integer i = 0; i != NDIM; ++i) {
					dX[i] = l->X[i] - Xcom[i];
				}
				auto this_L = Lphi << dX;
				l->phi += this_L();
				for (integer d = 0; d != NDIM; ++d) {
					l->g[d] -= this_L(d);
				}
			}
		} else {
			for (integer i = 0; i != nchild; ++i) {
				children[i]->compute_expansions(Lphi, Xcom);
			}
		}
	}

	void compute_interactions(std::vector<cell*> cells) {
		static std::vector<std::vector<cell*>> __child_cells(100);

		std::vector<cell*>& child_cells = __child_cells[level];
		child_cells.resize(0);
		auto l = cells.begin();

		multipole<simd_vector> Mvec;
		space_vector<simd_vector> Rvec;

		int index = 0;

		for (auto l = cells.begin(); l != cells.end(); ++l) {
			if (this->is_well_separated_from(**l)) {
				if ((*l)->M() != 0.0 && M() != 0.0) {
					for (integer d = 0; d != MP; ++d) {
						Mvec[d][index] = ((*l)->M)[d];
					}
					for (integer d = 0; d != NDIM; ++d) {
						Rvec[d][index] = Xcom[d] - (*l)->Xcom[d];
					}
					++index;
					if (index == simd_len) {
						expansion<simd_vector> Lvec;
						for (integer l = 0; l != LP; ++l) {
							Lvec[l] = 0.0;
						}
						multipole_interaction(Lvec, M, Mvec, Rvec);
						for (integer l = 0; l != LP; ++l) {
							Lphi[l] += Lvec[l].sum();
						}
						index = 0;
					}
				}
			} else {
				if ((*l)->is_leaf) {
					if (this->is_leaf) {
						if (this != *l) {
							P2P(std::move(*l));
						} else {
							P2Pself(std::move(*l));

						}
					} else {
						child_cells.push_back(std::move(*l));
					}
				} else {
					for (integer c = 0; c != nchild; ++c) {
						child_cells.push_back((*l)->children[c]);
					}
				}
			}
		}
		if (index > 0) {
			for (; index < simd_len; ++index) {
				Rvec[0][index] = 1.0;
				Rvec[1][index] = 1.0;
				Rvec[2][index] = 1.0;
				for (integer l = 0; l != MP; ++l) {
					Mvec[l][index] = 0.0;
				}
			}
			expansion<simd_vector> Lvec;
			for (integer l = 0; l != LP; ++l) {
				Lvec[l] = 0.0;
			}
			multipole_interaction(Lvec, M, Mvec, Rvec);
			for (integer l = 0; l != LP; ++l) {
				Lphi[l] += Lvec[l].sum();
			}
		}

		if (is_leaf) {
			if (!child_cells.empty()) {
				this->compute_interactions(std::move(child_cells));
			}
		} else {
			for (integer c = 0; c != nchild; ++c) {
				children[c]->compute_interactions(child_cells);
			}
		}

	}

	void P2P(cell* other) {
		if (!(*this < *other)) {
			return;
		}

		for (auto n = particles.begin(); n != particles.end(); ++n) {
			int i = 0;
			simd_vector mX;
			simd_vector mY;
			simd_vector mZ;
			simd_vector mM;
			simd_vector phi, gx, gy, gz;
			const int sz = other->particles.size();
			while (i < sz) {
				for (int j = 0; j < simd_len; ++j) {
					if (i < sz) {
						const auto& p = other->particles;
						mX[j] = p[i].X[0];
						mY[j] = p[i].X[1];
						mZ[j] = p[i].X[2];
						mM[j] = p[i].M;
					} else {
						mX[j] = 0.0;
						mY[j] = 0.0;
						mZ[j] = 0.0;
						mM[j] = 0.0;
					}
					++i;
				}
				i -= simd_len;

				const simd_vector dX = simd_vector(n->X[0]) - mX;
				const simd_vector dY = simd_vector(n->X[1]) - mY;
				const simd_vector dZ = simd_vector(n->X[2]) - mZ;
				const simd_vector r2 = dX * dX + dY * dY + dZ * dZ;
				const simd_vector rinv = simd_vector(1.) / sqrt(r2);
				auto tmp1 = mM * rinv;
				auto tmp2 = n->M * rinv;
				n->phi -= tmp1.sum();
				phi = -tmp2;
				tmp1 /= r2;
				tmp2 /= r2;
				n->g[0] -= (dX * tmp1).sum();
				n->g[1] -= (dY * tmp1).sum();
				n->g[2] -= (dZ * tmp1).sum();
				gx = dX * tmp2;
				gy = dY * tmp2;
				gz = dZ * tmp2;
				for (int j = 0; j < simd_len && i < sz; ++j) {
					auto& p = other->particles;
					if (i < sz) {
						p[i].phi += phi[j];
						p[i].g[0] += gx[j];
						p[i].g[1] += gy[j];
						p[i].g[2] += gz[j];
					}
					++i;
				}
			}
		}
	}

	void P2Pself(const cell* other) {
		expansion<real> D;
		space_vector < real > dX;
		real r2, r3inv, rinv;
		for (auto n = particles.begin(); n != particles.end(); ++n) {
			auto m = n;
			m++;
			for (; m != other->particles.end(); ++m) {
				r2 = 0.0;
				for (integer d = 0; d != NDIM; ++d) {
					dX[d] = n->X[d] - m->X[d];
					r2 += dX[d] * dX[d];
				}
				rinv = real(1.0) / std::sqrt(r2);
				r3inv = rinv / r2;
				n->phi -= m->M * rinv;
				m->phi -= n->M * rinv;
				for (integer d = 0; d != NDIM; ++d) {
					n->g[d] -= dX[d] * m->M * r3inv;
					m->g[d] += dX[d] * n->M * r3inv;
				}
			}
		}
	}

	void compute_multipoles() {
		Lphi = 0.0;
		for (integer d = 0; d != NDIM; ++d) {
			Xcom[d] = 0.0;
		}
		for (auto l = particles.begin(); l != particles.end(); ++l) {
			(l)->phi = 0.0;
			for (integer d = 0; d != NDIM; ++d) {
				(l)->g[d] = 0.0;
			}
		}
		real m = 0.0;
		if (is_leaf) {
			for (auto l = particles.begin(); l != particles.end(); ++l) {
				for (integer d = 0; d != NDIM; ++d) {
					Xcom[d] += l->X[d] * l->M;
				}
				m += l->M;
			}
		} else {
			for (integer l = 0; l != nchild; ++l) {
				auto& c = children[l];
				c->compute_multipoles();
				for (integer d = 0; d != NDIM; ++d) {
					Xcom[d] += c->Xcom[d] * c->M();
				}
				m += c->M();
			}
		}
		if (m != real(0.0)) {
			for (integer d = 0; d != NDIM; ++d) {
				Xcom[d] /= m;
			}
		}

		real rmax_sal, rmax_ben;

		rmax_ben = 0.0;
		rmax_sal = 0.0;
		for (integer d = 0; d != NDIM; ++d) {
			rmax_sal += std::pow(std::max(Xmax[d] - Xcom[d], Xcom[d] - Xmin[d]),
					real(2));
		}
		rmax_sal = std::sqrt(rmax_sal);
		assert(rmax_sal > 0.0);

		M = 0.0;
		if (is_leaf) {
			rmax_ben = 0.0;
			for (auto l = particles.begin(); l != particles.end(); ++l) {
				multipole<real> this_M;
				this_M = 0.0;
				this_M() = l->M;
				space_vector < real > dX;
				for (integer i = 0; i != NDIM; ++i) {
					dX[i] = l->X[i] - Xcom[i];
				}
				this_M >>= dX;
				M += this_M;
				real this_r = 0.0;
				for (integer d = 0; d != NDIM; ++d) {
					this_r += std::pow(Xcom[d] - l->X[d], real(2));
				}
				this_r = std::sqrt(this_r);
				rmax_ben = std::max(rmax_ben, this_r);
			}
		} else {
			for (integer l = 0; l != nchild; ++l) {
				auto& c = children[l];
				space_vector < real > dX;
				for (integer d = 0; d != NDIM; ++d) {
					dX[d] = c->Xcom[d] - Xcom[d];
				}
				M += c->M >> dX;
				real this_r = 0.0;
				for (integer d = 0; d != NDIM; ++d) {
					this_r += std::pow(Xcom[d] - c->Xcom[d], real(2));
				}
				this_r = std::sqrt(this_r);
				rmax_ben = std::max(rmax_ben, this_r + c->Rmax);
			}
		}
		assert(rmax_ben >= 0.0);
		Rmax = std::min(rmax_sal, rmax_ben);
	}

	void direct_interaction_at(const space_vector<real>& X, real& phi,
			space_vector<real>& g) const {
		if (level == 0) {
			for (integer d = 0; d != NDIM; ++d) {
				g[d] = 0.0;
			}
			phi = 0.0;
		}
		if (is_leaf) {
			real r2, rinv, r3inv;
			space_vector < real > dX;
			for (auto m = particles.begin(); m != particles.end(); ++m) {
				r2 = 0.0;
				for (integer d = 0; d != NDIM; ++d) {
					dX[d] = X[d] - m->X[d];
					r2 += dX[d] * dX[d];
				}
				if (r2 > real(0.0)) {
					rinv = real(1.0) / std::sqrt(r2);
					r3inv = rinv / r2;
				} else {
					rinv = r3inv = 0.0;
				}
				phi -= m->M * rinv;
				for (integer d = 0; d != NDIM; ++d) {
					g[d] -= dX[d] * m->M * r3inv;
				}
			}
		} else {
			for (integer i = 0; i != nchild; ++i) {
				children[i]->direct_interaction_at(X, phi, g);
			}
		}
	}

	real compute_error(const cell& root) {
		if (level == 0) {
			top10.clear();
		}
		static int counter = 0;
		static real start_time;

		real err = 0.0;
		if (is_leaf) {
			std::vector<particle*> parts(particles.size());
			auto iter = particles.begin();
			for (std::size_t i = 0; i != particles.size(); ++i) {
				parts[i] = &(*iter);
				++iter;
			}
			if (!analytic_computed) {
				const std::size_t psize = parts.size();
#pragma omp parallel for
				for (std::size_t i = 0; i < psize; i++) {
					auto l = parts[i];
					root.direct_interaction_at(l->X, l->phi_analytic,
							l->g_analytic);
				}
			}

			for (std::size_t i = 0; i != parts.size(); i++) {
				auto l = parts[i];
				space_vector < real > &g = l->g_analytic;
				const real abs_g_exact = std::sqrt(
						g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
				const real abs_g_approx = std::sqrt(
						l->g[0] * l->g[0] + l->g[1] * l->g[1]
								+ l->g[2] * l->g[2]);
				const real eps = std::abs(abs_g_approx - abs_g_exact)
						/ abs_g_exact;
				{
					if (top10.size() < nparts / 100) {
						top10.insert(eps);
					} else {
						if (eps > *(top10.begin())) {
							top10.erase(top10.begin());
							top10.insert(eps);
						}
					}
				}
				err += eps;
			}
		} else {
			for (integer c = 0; c != nchild; ++c) {
				err += children[c]->compute_error(root);
			}
		}
		if (level == 0) {
			err /= nparts;
			//	printf("Time for direction solution = %e\n", (clock() - t0) / real(CLOCKS_PER_SEC));
		}

		if (counter == 0) {
			printf("Starting direct calculation\n");
			start_time = real(clock()) / real(CLOCKS_PER_SEC);
		}
		if (counter < ncells) {
			printf("%.2f\r", 100.0 * counter / real(ncells));
			fflush (stdout);
			++counter;
			if (counter == ncells) {
				analytic_computed = true;
				printf("100.00\n took %e seconds\n",
						real(clock()) / real(CLOCKS_PER_SEC - start_time));
				++counter;
				sleep(2);
			}
		}
		return err;
	}

	void initialize() {
		std::list<particle> these_parts;
		auto randx =
				[]() {
					return (real(rand())+real(0.5)) / (real(RAND_MAX) + real(1.0));
				};
		const integer ngalaxies = 5;
		Xmax[0] = Xmax[1] = Xmax[2] = 0.0;
		Xmin[0] = Xmin[1] = Xmin[2] = 0.0;
		for (integer gi = 0; gi != ngalaxies; ++gi) {
			const real this_m = randx();
			const real a = randx();
			const real x0 = real(1000.0) * (real(2.0) * randx() - real(1.0));
			const real y0 = real(1000.0) * (real(2.0) * randx() - real(1.0));
			const real z0 = real(1000.0) * (real(2.0) * randx() - real(1.0));
			for (integer i = 0; i != nparts / ngalaxies; ++i) {
				particle p;
				real x, y, z;
				real r0, p0;
				do {
					p0 = randx();
					r0 = (a / (real(1.0) - p0)) * (p0 + std::sqrt(p0));
				} while (r0 > real(1000.0) * a);
				real r;
				do {
					x = r0 * (real(2.0) * randx() - real(1.0));
					y = r0 * (real(2.0) * randx() - real(1.0));
					z = r0 * (real(2.0) * randx() - real(1.0));
					r = std::sqrt(x * x + y * y + z * z);
				} while (r > r0);
				x += x0;
				y += y0;
				z += z0;
				Xmax[0] = std::max(x, Xmax[0]);
				Xmax[1] = std::max(y, Xmax[1]);
				Xmax[2] = std::max(z, Xmax[2]);
				Xmin[0] = std::min(x, Xmin[0]);
				Xmin[1] = std::min(y, Xmin[1]);
				Xmin[2] = std::min(z, Xmin[2]);
				p.X[0] = x;
				p.X[1] = y;
				p.X[2] = z;
				p.M = randx() * this_m;
				these_parts.push_back(std::move(p));
			}
		}
		create_tree(Xmin, Xmax, these_parts, 0);

		//	printf("Created a tree with maximum level %i\n", int(maxlev));
	}

	integer create_tree(const space_vector<real>& xmin,
			const space_vector<real>& xmax, std::list<particle> parts,
			integer lev) {
		integer maxlev = lev;
		level = lev;
		Xmax = xmax;
		is_leaf = true;
		Xmin = xmin;
		space_vector<real> this_xmax, this_xmin;
		std::list<particle> these_parts;
		std::array < std::array<std::array<std::list<particle>, 2>, 2>, 2
				> child_parts;
		auto tmp = std::move(parts);
		if (tmp.size() > ncrit) {
			is_leaf = false;
			for (integer i = 0; i < nchild; ++i) {
				children[i] = new cell();
			}
			auto l = tmp.begin();
			while (l != tmp.end()) {
				const auto i = std::max(integer(0),
						std::min(integer(1),
								integer(
										real(2) * (l->X[0] - Xmin[0])
												/ (Xmax[0] - Xmin[0]))));
				const auto j = std::max(integer(0),
						std::min(integer(1),
								integer(
										real(2) * (l->X[1] - Xmin[1])
												/ (Xmax[1] - Xmin[1]))));
				const auto k = std::max(integer(0),
						std::min(integer(1),
								integer(
										real(2) * (l->X[2] - Xmin[2])
												/ (Xmax[2] - Xmin[2]))));
				child_parts[i][j][k].push_back(std::move(*l));
				l = tmp.erase(l);
			}
			for (integer i = 0; i < 2; ++i) {
				this_xmin[0] = Xmin[0]
						+ real(0.5) * real(i) * (Xmax[0] - Xmin[0]);
				this_xmax[0] = this_xmin[0] + real(0.5) * (Xmax[0] - Xmin[0]);
				for (integer j = 0; j < 2; ++j) {
					this_xmin[1] = Xmin[1]
							+ real(0.5) * real(j) * (Xmax[1] - Xmin[1]);
					this_xmax[1] = this_xmin[1]
							+ real(0.5) * (Xmax[1] - Xmin[1]);
					for (integer k = 0; k < 2; ++k) {
						this_xmin[2] = Xmin[2]
								+ real(0.5) * real(k) * (Xmax[2] - Xmin[2]);
						this_xmax[2] = this_xmin[2]
								+ real(0.5) * (Xmax[2] - Xmin[2]);
						maxlev = std::max(
								children[i * 4 + j * 2 + k]->create_tree(
										this_xmin, this_xmax,
										std::move(child_parts[i][j][k]),
										lev + 1), maxlev);
					}
				}
			}
		}
		//	printf( "%i\n", maxlev);
		particles = std::vector < particle > (tmp.begin(), tmp.end());
		return maxlev;
	}

	real solve() {
		real t00 = clock();
		compute_multipoles();
		std::vector < cell*> cells;
		for (integer i = 0; i != nchild; ++i) {
			cells.push_back(children[i]);
		}
		Lphi = 0.0;
		for (integer d = 0; d != NDIM; ++d) {
			Xcom[d] = 0.0;
		}
		for (integer i = 0; i != nchild; ++i) {
			children[i]->compute_interactions(cells);
		}
		for (integer i = 0; i != nchild; ++i) {
			children[i]->compute_expansions(Lphi, Xcom);
		}
		real t1 = clock();
		real total_time = (t1 - t00) / real(CLOCKS_PER_SEC);
		return total_time;
	}
};

int main() {
	printf("%16s %16s %16s %16s %16s %16s %16s\n", "theta", "nparts",
			"solve_time", "err", "err99", "gsum", "tsum");
	cell root;
	root.initialize();
	root.output("output.txt");
	for (theta = 1.0; theta > 0.15; theta -= 0.1) {
		space_vector < real > g_err;
		for (integer d = 0; d != NDIM; ++d) {
			g_err[d] = 0.0;
		}
		real solve_time = root.solve();
		//	srand(1);
		const real this_err = root.compute_error(root);
		//const real this_err = 0.0;
		auto g = root.force_sum(false);
		auto gnorm = root.force_sum(true);
		auto t = root.torque_sum(false);
		auto tnorm = root.torque_sum(true);
		real gsum = 0.0;
		real tsum = 0.0;
		real gsum_norm = 0.0;
		real tsum_norm = 0.0;
		real err99 = *(top10.begin());
		for (integer d = 0; d != NDIM; ++d) {
			gsum += g[d] * g[d];
			gsum_norm += gnorm[d] * gnorm[d];
			tsum += t[d] * t[d];
			tsum_norm += tnorm[d] * tnorm[d];
		}
		gsum /= gsum_norm;
		tsum /= tsum_norm;
		gsum = std::sqrt(gsum);
		tsum = std::sqrt(tsum);
		//	printf("Force  sum = %e \n", gsum);
		//	printf("Torque sum = %e \n", tsum);
		FILE* fp = fopen("data.dat", "at");
		printf("%16e %16lli %16e %16e %16e %16e %16e\n", theta, nparts,
				solve_time, this_err, err99, gsum, tsum);
		fprintf(fp, "%e %lli %e %e %e %e %e\n", theta, nparts, solve_time,
				this_err, err99, gsum, tsum);
		fclose(fp);
	}
	return 0;
}
