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

#ifndef THUNDEREGG_DOMAIN_TOOLS_H
#define THUNDEREGG_DOMAIN_TOOLS_H
/**
 * @file
 *
 * @brief DomainTools class
 */

#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Vector.h>

namespace ThunderEgg {
/**
 * @brief Various tools for filling in values in a domain
 */
class DomainTools
{
private:
  /**
   * @brief This is not intended to be called directly. Set the values for a vector with the given
   * function.
   *
   * @tparam D the number of cartesian dimensions
   * @tparam T function type
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   */
  template<int D, typename T>
  static void
  _SetValues(const Domain<D>& domain, Vector<D>& vec, int component_index, T func)
  {
    SetValues(domain, vec, component_index, func);
  }

  /**
   * @brief This is not intended to be called directly. Set the values for a vector with the given
   * function.
   *
   * @tparam D the number of cartesian dimensions
   * @tparam T function type
   * @tparam Args additional functions
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   * @param args additional functions for additional components
   */
  template<int D, typename T, typename... Args>
  static void
  _SetValues(const Domain<D>& domain, Vector<D>& vec, int component_index, T func, Args... args)
  {
    SetValues(domain, vec, component_index, func);
    _SetValues(domain, vec, component_index + 1, args...);
  }

  /**
   * @brief This is not intended to be called directly. Set the values (including ghost values) for
   * a vector with the given function.
   *
   * @tparam D the number of cartesian dimensions
   * @tparam T function type
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   */
  template<int D, typename T>
  static void
  _SetValuesWithGhost(const Domain<D>& domain, Vector<D>& vec, int component_index, T func)
  {
    SetValuesWithGhost(domain, vec, component_index, func);
  }

  /**
   * @brief This is not intended to be called directly. Set the values (including ghost values) for
   * a vector with the given function.
   *
   * @tparam D the number of cartesian dimensions
   * @tparam T function type
   * @tparam Args additional functions
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   * @param args additional functions for additional components
   */
  template<int D, typename T, typename... Args>
  static void
  _SetValuesWithGhost(const Domain<D>& domain,
                      Vector<D>& vec,
                      int component_index,
                      T func,
                      Args... args)
  {
    SetValuesWithGhost(domain, vec, component_index, func);
    _SetValuesWithGhost(domain, vec, component_index + 1, args...);
  }

public:
  /**
   * @brief Given a path info object, get the coordinate from a given index into the patch.
   *
   * If one of the values is -1 or ns[axis] it will give the coordinate of the cooresponding
   * interface
   *
   * @tparam D the number of cartesian dimensions
   * @param pinfo the patch to get the coordinate for
   * @param coord the index in the patch
   * @param real_coord (output) the coordnitate of the index
   */
  template<int D>
    requires is_supported_dimension<D>
  static void
  GetRealCoord(const PatchInfo<D>& pinfo,
               const std::array<int, D>& coord,
               std::array<double, D>& real_coord);

  /**
   * @brief Given a path info object, get the coordinate from a given index into the patch.
   *
   * @tparam D the number of cartesian dimensions
   * @param pinfo the patch to get the coordinate for
   * @param coord the index in the patch
   * @param real_coord (output) the coordnitate of the index
   */
  template<int D>
    requires is_supported_dimension<D>
  static void
  GetRealCoordGhost(const PatchInfo<D>& pinfo,
                    const std::array<int, D>& coord,
                    std::array<double, D>& real_coord);

  /**
   * @brief Given a path info object and a side of the patch, get the coordinate from a given
   * index into the interface of the patch.
   *
   * @tparam D the number of cartesian dimensions
   * @param pinfo the patch to get the coordinate for
   * @param coord the index in the patch
   * @param s the side of the patch that the boundary is on
   * @param real_coord (output) the coordnitate of the index
   */
  template<int D>
    requires is_supported_dimension<D>
  static void
  GetRealCoordBound(const PatchInfo<D>& pinfo,
                    const std::array<int, D - 1>& coord,
                    Side<D> s,
                    std::array<double, D>& real_coord);

  /**
   * @brief Set the values for a vector with the given function.
   *
   * @tparam D the number of cartesian dimensions
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   */
  template<int D>
    requires is_supported_dimension<D>
  static void
  SetValues(const Domain<D>& domain,
            Vector<D>& vec,
            int component_index,
            std::function<double(const std::array<double, (int)D>&)> func);

  /**
   * @brief Set the values for a vector with the given function.
   *
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   */
  static void
  SetValues(const Domain<3>& domain,
            Vector<3>& vec,
            int component_index,
            std::function<double(double, double, double)> func);

  /**
   * @brief Set the values for a vector with the given function.
   *
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   */
  static void
  SetValues(const Domain<2>& domain,
            Vector<2>& vec,
            int component_index,
            std::function<double(double, double)> func);

  /**
   * @brief Set the values for a vector with the given functions
   *
   * @tparam D the number of cartesian dimensions
   * @tparam Args additional functions
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   * @param args additional functions for additional components
   */
  template<int D, typename... Args>
  static void
  SetValues(const Domain<D>& domain,
            Vector<D>& vec,
            std::function<double(const std::array<double, D>&)> func,
            Args... args)
  {
    _SetValues(domain, vec, 0, func, args...);
  }

  /**
   * @brief Set the values (including ghost values) for a vector with the given function.
   *
   * @tparam D the number of cartesian dimensions
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   */
  template<int D>
    requires is_supported_dimension<D>
  static void
  SetValuesWithGhost(const Domain<D>& domain,
                     Vector<D>& vec,
                     int component_index,
                     std::function<double(const std::array<double, (int)D>&)> func);

  /**
   * @brief Set the values (including ghost values) for a vector with the given function.
   *
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   */
  static void
  SetValuesWithGhost(const Domain<3>& domain,
                     Vector<3>& vec,
                     int component_index,
                     std::function<double(double, double, double)> func);

  /**
   * @brief Set the values (including ghost values) for a vector with the given function.
   *
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   */
  static void
  SetValuesWithGhost(const Domain<2>& domain,
                     Vector<2>& vec,
                     int component_index,
                     std::function<double(double, double)> func);

  /**
   * @brief Set the values (including ghost values) for a vector with the given functions
   *
   * @tparam D the number of cartesian dimensions
   * @tparam T function type
   * @tparam Args additional functions
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   * @param args additional functions for additional components
   */
  template<int D, typename... Args>
  static void
  SetValuesWithGhost(const Domain<D>& domain,
                     Vector<D>& vec,
                     std::function<double(const std::array<double, D>&)> func,
                     Args... args)
  {
    _SetValuesWithGhost(domain, vec, 0, func, args...);
  }

  /**
   * @brief Set the values (including ghost values) for a vector with the given functions
   *
   * @tparam Args additional functions
   * @param domain the domain that we are setting values for
   * @param vec the vector to set values in
   * @param component_index the component to set
   * @param func the function
   * @param args additional functions for additional components
   */
  template<typename... Args>
  static void
  SetValuesWithGhost(const Domain<3>& domain,
                     Vector<3>& vec,
                     std::function<double(double, double, double)> func,
                     Args... args)
  {
    _SetValuesWithGhost(domain, vec, 0, func, args...);
  }

  /**
   * @brief Integrate a vector over the domain.
   *
   * @param u the vector
   * @return double the result of the integral
   */
  template<int D>
    requires is_supported_dimension<D>
  static double
  Integrate(const Domain<D>& domain, const Vector<D>& u);
};

extern template void
DomainTools::GetRealCoord<2>(const PatchInfo<2>& pinfo,
                             const std::array<int, 2>& coord,
                             std::array<double, 2>& real_coord);

extern template void
DomainTools::GetRealCoord<3>(const PatchInfo<3>& pinfo,
                             const std::array<int, 3>& coord,
                             std::array<double, 3>& real_coord);

extern template void
DomainTools::GetRealCoordGhost<2>(const PatchInfo<2>& pinfo,
                                  const std::array<int, 2>& coord,
                                  std::array<double, 2>& real_coord);

extern template void
DomainTools::GetRealCoordGhost<3>(const PatchInfo<3>& pinfo,
                                  const std::array<int, 3>& coord,
                                  std::array<double, 3>& real_coord);

extern template void
DomainTools::GetRealCoordBound<2>(const PatchInfo<2>& pinfo,
                                  const std::array<int, 1>& coord,
                                  Side<2> s,
                                  std::array<double, 2>& real_coord);

extern template void
DomainTools::GetRealCoordBound<3>(const PatchInfo<3>& pinfo,
                                  const std::array<int, 2>& coord,
                                  Side<3> s,
                                  std::array<double, 3>& real_coord);

extern template void
DomainTools::SetValues<2>(const Domain<2>& domain,
                          Vector<2>& vec,
                          int component_index,
                          std::function<double(const std::array<double, (int)2>&)> func);

extern template void
DomainTools::SetValues<3>(const Domain<3>& domain,
                          Vector<3>& vec,
                          int component_index,
                          std::function<double(const std::array<double, (int)3>&)> func);

extern template void
DomainTools::SetValues<2>(const Domain<2>& domain,
                          Vector<2>& vec,
                          int component_index,
                          std::function<double(const std::array<double, (int)2>&)> func);

extern template void
DomainTools::SetValuesWithGhost<3>(const Domain<3>& domain,
                                   Vector<3>& vec,
                                   int component_index,
                                   std::function<double(const std::array<double, (int)3>&)> func);

extern template double
DomainTools::Integrate<2>(const Domain<2>& domain, const Vector<2>& u);
extern template double
DomainTools::Integrate<3>(const Domain<3>& domain, const Vector<3>& u);

} // namespace ThunderEgg
#endif