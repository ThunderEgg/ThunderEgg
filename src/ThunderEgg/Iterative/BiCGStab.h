/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2019-2021 Scott Aiton
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

#ifndef THUNDEREGG_ITERATIVE_BICGSTAB_H
#define THUNDEREGG_ITERATIVE_BICGSTAB_H
/**
 * @file
 *
 * @brief BiCGStab class
 */

#include <ThunderEgg/Iterative/BreakdownError.h>
#include <ThunderEgg/Iterative/Solver.h>
#include <ThunderEgg/Operator.h>
#include <ThunderEgg/Timer.h>

namespace ThunderEgg::Iterative {
/**
 * @brief BiCGStab iterative solver.
 *
 * @tparam D the number of Cartesian dimensions
 */
template<int D>
class BiCGStab : public Solver<D>
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

public:
  /**
   * @brief Clone this solver
   *
   * @return BiCGStab<D>* a newly allocated copy of this solver
   */
  BiCGStab<D>* clone() const override { return new BiCGStab<D>(*this); }
  /**
   * @brief Set the maximum number of iterations
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
  int solve(const Operator<D>& A, Vector<D>& x, const Vector<D>& b, const Operator<D>* Mr = nullptr, bool output = false, std::ostream& os = std::cout) const override
  {
    Vector<D> resid = b.getZeroClone();

    Vector<D> ms;
    Vector<D> mp;
    if (Mr != nullptr) {
      ms = b.getZeroClone();
      mp = b.getZeroClone();
    }

    A.apply(x, resid);
    resid.scaleThenAdd(-1, b);

    double r0_norm = b.twoNorm();
    Vector<D> rhat = resid;
    Vector<D> p = resid;
    Vector<D> ap = b.getZeroClone();
    Vector<D> as = b.getZeroClone();

    Vector<D> s = x.getZeroClone();
    double rho = rhat.dot(resid);

    int num_its = 0;
    double residual = resid.twoNorm() / r0_norm;
    if (output) {
      char buf[100];
      sprintf(buf, "%5d %16.8e\n", num_its, residual);
      os << std::string(buf);
    }
    if (r0_norm == 0) {
      return num_its;
    }
    while (residual > tolerance && num_its < max_iterations) {
      if (timer) {
        timer->start("Iteration");
      }

      if (rho == 0) {
        throw BreakdownError("BiCGStab broke down, rho was 0 on iteration " + std::to_string(num_its));
      }

      if (Mr != nullptr) {
        Mr->apply(p, mp);
        A.apply(mp, ap);
      } else {
        A.apply(p, ap);
      }
      double alpha = rho / rhat.dot(ap);
      s.copy(resid);
      s.addScaled(-alpha, ap);
      if (s.twoNorm() / r0_norm <= tolerance) {
        x.addScaled(alpha, p);
        if (timer) {
          timer->stop("Iteration");
        }
        break;
      }
      if (Mr != nullptr) {
        Mr->apply(s, ms);
        A.apply(ms, as);
      } else {
        A.apply(s, as);
      }
      double omega = as.dot(s) / as.dot(as);
      // update x and residual
      if (Mr != nullptr) {
        x.addScaled(alpha, mp, omega, ms);
      } else {
        x.addScaled(alpha, p, omega, s);
      }
      resid.addScaled(-alpha, ap, -omega, as);

      double rho_new = resid.dot(rhat);
      double beta = rho_new * alpha / (rho * omega);
      p.addScaled(-omega, ap);
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
    return num_its;
  }
};
} // namespace ThunderEgg::Iterative
extern template class ThunderEgg::Iterative::BiCGStab<2>;
extern template class ThunderEgg::Iterative::BiCGStab<3>;
#endif