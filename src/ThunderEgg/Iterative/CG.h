/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef THUNDEREGG_ITERATIVE_CG_H
#define THUNDEREGG_ITERATIVE_CG_H
/**
 * @file
 *
 * @brief CG class
 */

#include <ThunderEgg/Iterative/BreakdownError.h>
#include <ThunderEgg/Iterative/Solver.h>
#include <ThunderEgg/Operator.h>
#include <ThunderEgg/Timer.h>

namespace ThunderEgg::Iterative {
/**
 * @brief CG iterative solver.
 *
 * @tparam D the number of Cartesian dimensions
 */
template<int D>
class CG : public Solver<D>
{
private:
  /**
   * @brief The maximum number of iterations
   */
  int max_iterations = 1000;
  /**
   * @brief The maximum number of iterations
   */
  double tolerance = 1e-12;
  /**
   * @brief The timer
   */
  std::shared_ptr<Timer> timer = nullptr;

  void applyWithPreconditioner(const Operator<D>* M_l,
                               const Operator<D>& A,
                               const Operator<D>* M_r,
                               const Vector<D>& x,
                               Vector<D>& b) const
  {
    if (M_l == nullptr && M_r == nullptr) {
      A.apply(x, b);
    } else if (M_l == nullptr && M_r != nullptr) {
      Vector<D> tmp = b.getZeroClone();
      M_r->apply(x, tmp);
      A.apply(tmp, b);
    }
  }

public:
  /**
   * @brief Clone this solver
   *
   * @return CG<D>* a newly allocated copy of this solver
   */
  CG<D>* clone() const override { return new CG<D>(*this); }
  /**
   * @brief Set the maximum number of iterations.
   *
   * Default is 1000
   *
   * @param max_iterations_in the maximum number of iterations
   */
  void setMaxIterations(int max_iterations_in) { max_iterations = max_iterations_in; };
  /**
   * @brief Get the maximum number of iterations
   *
   * Default is 1000
   *
   * @return int the maximum number of iterations
   */
  int getMaxIterations() const { return max_iterations; }
  /**
   * @brief Set the stopping tolerance
   *
   * Default is 1e-12
   *
   * @param tolerance_in the stopping tolerance
   */
  void setTolerance(double tolerance_in) { tolerance = tolerance_in; };
  /**
   * @brief Get the stopping tolerance
   *
   * Default is 1e-12
   *
   * @return double the stopping tolerance
   */
  double getTolerance() const { return tolerance; }
  /**
   * @brief Set the Timer object
   *
   * @param timer_in the Timer
   */
  void setTimer(std::shared_ptr<Timer> timer_in) { timer = timer_in; }

  /**
   * @brief Get the Timer object
   *
   * @return std::shared_ptr<Timer> the Timer
   */
  std::shared_ptr<Timer> getTimer() const { return timer; }

public:
  int solve(const Operator<D>& A,
            Vector<D>& x,
            const Vector<D>& b,
            const Operator<D>* Mr = nullptr,
            bool output = false,
            std::ostream& os = std::cout) const override
  {
    Vector<D> resid = b.getZeroClone();

    A.apply(x, resid);
    resid.scaleThenAdd(-1, b);

    Vector<D> initial_guess = x;
    x.set(0);

    double r0_norm = b.twoNorm();
    Vector<D> p = resid;
    Vector<D> ap = b.getZeroClone();

    double rho = resid.dot(resid);

    int num_its = 0;
    if (r0_norm == 0) {
      return num_its;
    }
    double residual = resid.twoNorm() / r0_norm;
    if (output) {
      char buf[100];
      sprintf(buf, "%5d %16.8e\n", num_its, residual);
      os << std::string(buf);
    }
    while (residual > tolerance && num_its < max_iterations) {
      if (timer) {
        timer->start("Iteration");
      }

      if (rho == 0) {
        throw BreakdownError("CG broke down, rho was 0 on iteration " + std::to_string(num_its));
      }

      applyWithPreconditioner(nullptr, A, Mr, p, ap);
      double alpha = rho / p.dot(ap);
      x.addScaled(alpha, p);
      resid.addScaled(-alpha, ap);

      double rho_new = resid.dot(resid);
      double beta = rho_new / rho;
      p.scaleThenAdd(beta, resid);

      num_its++;
      rho = rho_new;
      residual = resid.twoNorm() / r0_norm;

      if (output) {
        char buf[100];
        sprintf(buf, "%5d %16.8e\n", num_its, residual);
        os << std::string(buf);
      }
      if (timer) {
        timer->stop("Iteration");
      }
    }
    if (Mr != nullptr) {
      Mr->apply(x, resid);
      x.copy(resid);
    }
    x.add(initial_guess);
    return num_its;
  }
};
} // namespace ThunderEgg::Iterative
extern template class ThunderEgg::Iterative::CG<2>;
extern template class ThunderEgg::Iterative::CG<3>;
#endif