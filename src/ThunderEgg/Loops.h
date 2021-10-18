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

#ifndef THUNDEREGG_LOOP_H
#define THUNDEREGG_LOOP_H
/**
 * @file
 *
 * @brief Loop templates
 */

#include <array>

namespace ThunderEgg {
/**
 * @brief Dimension templated loops
 */
class Loop
{
private:
  template<int start, int stop, typename T>
  class ULoop
  {
  public:
    static void inline loop_loop(T lambda)
    {
      lambda(start);
      ULoop<start + 1, stop, T>::loop_loop(lambda);
    }
  };

  template<int start, typename T>
  class ULoop<start, start, T>
  {
  public:
    static void inline loop_loop(T lambda) { lambda(start); }
  };
  template<typename T>
  class ULoop<0, -1, T>
  {
  public:
    static void inline loop_loop(T lambda) {}
  };

public:
  /**
   * @brief Unroll a fixed-length loop
   *
   * @tparam start the starting index
   * @tparam stop the stopping index, inclusive
   * @tparam T lambda type
   * @param lambda
   */
  template<int start, int stop, typename T>
  static inline void Unroll(T lambda)
  {
    ULoop<start, stop, T>::loop_loop(lambda);
  }

private:
  /**
   * @tparam Dir
   * @tparam T
   * @tparam A
   */
  template<int D, int Dir, typename T, typename A>
  class NestedLoop
  {
  public:
    static void inline nested_loop_loop(A coord, A start, A end, T lambda)
    {
      for (coord[Dir] = start[Dir]; coord[Dir] <= end[Dir]; coord[Dir]++) {
        NestedLoop<D, Dir - 1, T, A>::nested_loop_loop(coord, start, end, lambda);
      }
    }
  };

  template<int D, typename T, typename A>
  class NestedLoop<D, 0, T, A>
  {
  public:
    static void inline nested_loop_loop(A coord, A start, A end, T lambda)
    {
      for (coord[0] = start[0]; coord[0] <= end[0]; coord[0]++) {
        lambda(coord);
      }
    }
  };
  template<typename T, typename A>
  class NestedLoop<0, -1, T, A>
  {
  public:
    static void inline nested_loop_loop(A coord, A, A, T lambda) { lambda(coord); }
  };

public:
  /**
   * @brief Loop over a range of integer coordinates
   *
   * @tparam D the dimensions of the coordinates
   * @tparam T lambda type
   * @tparam A coordinate type
   * @param start initial coordinate
   * @param end final coordinate, inclusive
   * @param lambda the lambda function to call each time
   */
  template<int D, typename T, typename A>
  static inline void Nested(A start, A end, T lambda)
  {
    A coord = start;
    NestedLoop<D, D - 1, T, A>::nested_loop_loop(coord, start, end, lambda);
  }
  /**
   * @brief Loop over all the interior coordinates of a view
   *
   * @tparam D the dimension of the view
   * @tparam V the view type
   * @tparam T the lambda type
   * @param view the view to loop over
   * @param lambda the lambda function to call each time
   */
  template<int D, typename V, typename T>
  static inline void OverInteriorIndexes(const V& view, T lambda)
  {
    Nested<D>(view.getStart(), view.getEnd(), lambda);
  }
  /**
   * @brief Loop over all the coordinates of a view, including ghost cells
   *
   * @tparam D the dimension of the view
   * @tparam V the view type
   * @tparam T the lambda type
   * @param view the view to loop over
   * @param lambda the lambda function to call each time
   */
  template<int D, typename V, typename T>
  static inline void OverAllIndexes(const V& view, T lambda)
  {
    Nested<D>(view.getGhostStart(), view.getGhostEnd(), lambda);
  }
};
} // namespace ThunderEgg
#endif