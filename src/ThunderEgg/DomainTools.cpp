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

#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/RuntimeError.h>

namespace ThunderEgg {
template<int D>
  requires is_supported_dimension<D>
void
DomainTools::GetRealCoord(const PatchInfo<D>& pinfo,
                          const std::array<int, D>& coord,
                          std::array<double, D>& real_coord)
{
  Loop::Unroll<0, D - 1>([&](int dir) {
    if (coord[dir] == -1) {
      real_coord[dir] = pinfo.starts[dir];
    } else if (coord[dir] == pinfo.ns[dir]) {
      real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] * pinfo.ns[dir];
    } else {
      real_coord[dir] =
        pinfo.starts[dir] + pinfo.spacings[dir] / 2.0 + pinfo.spacings[dir] * coord[dir];
    }
  });
}

template<int D>
  requires is_supported_dimension<D>
void
DomainTools::GetRealCoordGhost(const PatchInfo<D>& pinfo,
                               const std::array<int, D>& coord,
                               std::array<double, D>& real_coord)
{
  Loop::Unroll<0, D - 1>([&](int dir) {
    real_coord[dir] =
      pinfo.starts[dir] + pinfo.spacings[dir] / 2.0 + pinfo.spacings[dir] * coord[dir];
  });
}

template<int D>
  requires is_supported_dimension<D>
void
DomainTools::GetRealCoordBound(const PatchInfo<D>& pinfo,
                               const std::array<int, D - 1>& coord,
                               Side<D> s,
                               std::array<double, D>& real_coord)
{
  for (size_t dir = 0; dir < s.getAxisIndex(); dir++) {
    if (coord[dir] == -1) {
      real_coord[dir] = pinfo.starts[dir];
    } else if (coord[dir] == pinfo.ns[dir]) {
      real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] * pinfo.ns[dir];
    } else {
      real_coord[dir] =
        pinfo.starts[dir] + pinfo.spacings[dir] / 2.0 + pinfo.spacings[dir] * coord[dir];
    }
  }
  if (s.isLowerOnAxis()) {
    real_coord[s.getAxisIndex()] = pinfo.starts[s.getAxisIndex()];
  } else {
    real_coord[s.getAxisIndex()] = pinfo.starts[s.getAxisIndex()] +
                                   pinfo.spacings[s.getAxisIndex()] * pinfo.ns[s.getAxisIndex()];
  }
  for (size_t dir = s.getAxisIndex() + 1; dir < D; dir++) {
    if (coord[dir - 1] == -1) {
      real_coord[dir] = pinfo.starts[dir];
    } else if (coord[dir - 1] == pinfo.ns[dir]) {
      real_coord[dir] = pinfo.starts[dir] + pinfo.spacings[dir] * pinfo.ns[dir];
    } else {
      real_coord[dir] =
        pinfo.starts[dir] + pinfo.spacings[dir] / 2.0 + pinfo.spacings[dir] * coord[dir - 1];
    }
  }
}

template<int D>
  requires is_supported_dimension<D>
void
DomainTools::SetValues(const Domain<D>& domain,
                       Vector<D>& vec,
                       int component_index,
                       std::function<double(const std::array<double, (int)D>&)> func)
{
  if (component_index >= vec.getNumComponents()) {
    throw RuntimeError("Invalid component to set");
  }
  std::array<double, D> real_coord;
  for (int i = 0; i < vec.getNumLocalPatches(); i++) {
    ComponentView<double, D> ld = vec.getComponentView(component_index, i);
    auto pinfo = domain.getPatchInfoVector()[i];
    Loop::Nested<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D>& coord) {
      GetRealCoord<D>(pinfo, coord, real_coord);
      ld[coord] = func(real_coord);
    });
  }
}

void
DomainTools::SetValues(const Domain<3>& domain,
                       Vector<3>& vec,
                       int component_index,
                       std::function<double(double, double, double)> func)
{
  if (component_index >= vec.getNumComponents()) {
    throw RuntimeError("Invalid component to set");
  }
  for (int i = 0; i < vec.getNumLocalPatches(); i++) {
    ComponentView<double, 3> ld = vec.getComponentView(component_index, i);
    const PatchInfo<3>& pinfo = domain.getPatchInfoVector()[i];
    double dx = pinfo.spacings[0];
    double dy = pinfo.spacings[1];
    double dz = pinfo.spacings[2];
    for (int zi = 0; zi < pinfo.ns[2]; zi++) {
      double z = pinfo.starts[2] + 0.5 * dz + zi * dz;
      for (int yi = 0; yi < pinfo.ns[1]; yi++) {
        double y = pinfo.starts[1] + 0.5 * dy + yi * dy;
        for (int xi = 0; xi < pinfo.ns[0]; xi++) {
          double x = pinfo.starts[0] + 0.5 * dx + xi * dx;
          ld(xi, yi, zi) = func(x, y, z);
        }
      }
    }
  }
}

void
SetValues(const Domain<2>& domain,
          Vector<2>& vec,
          int component_index,
          std::function<double(double, double)> func)
{
  if (component_index >= vec.getNumComponents()) {
    throw RuntimeError("Invalid component to set");
  }
  for (int i = 0; i < vec.getNumLocalPatches(); i++) {
    ComponentView<double, 2> ld = vec.getComponentView(component_index, i);
    const PatchInfo<2>& pinfo = domain.getPatchInfoVector()[i];
    double dx = pinfo.spacings[0];
    double dy = pinfo.spacings[1];
    for (int yi = 0; yi < pinfo.ns[1]; yi++) {
      double y = pinfo.starts[1] + 0.5 * dy + yi * dy;
      for (int xi = 0; xi < pinfo.ns[0]; xi++) {
        double x = pinfo.starts[0] + 0.5 * dx + xi * dx;
        ld(xi, yi) = func(x, y);
      }
    }
  }
}

template<int D>
  requires is_supported_dimension<D>
void
DomainTools::SetValuesWithGhost(const Domain<D>& domain,
                                Vector<D>& vec,
                                int component_index,
                                std::function<double(const std::array<double, (int)D>&)> func)
{
  if (component_index >= vec.getNumComponents()) {
    throw RuntimeError("Invalid component to set");
  }
  std::array<double, D> real_coord;
  for (int i = 0; i < vec.getNumLocalPatches(); i++) {
    ComponentView<double, D> ld = vec.getComponentView(component_index, i);
    auto pinfo = domain.getPatchInfoVector()[i];
    Loop::Nested<D>(ld.getGhostStart(), ld.getGhostEnd(), [&](const std::array<int, D>& coord) {
      GetRealCoordGhost<D>(pinfo, coord, real_coord);
      ld[coord] = func(real_coord);
    });
  }
}

void
DomainTools::SetValuesWithGhost(const Domain<3>& domain,
                                Vector<3>& vec,
                                int component_index,
                                std::function<double(double, double, double)> func)
{
  if (component_index >= vec.getNumComponents()) {
    throw RuntimeError("Invalid component to set");
  }
  for (int i = 0; i < vec.getNumLocalPatches(); i++) {
    ComponentView<double, 3> ld = vec.getComponentView(component_index, i);
    const PatchInfo<3>& pinfo = domain.getPatchInfoVector()[i];
    int num_ghost = pinfo.num_ghost_cells;
    double dx = pinfo.spacings[0];
    double dy = pinfo.spacings[1];
    double dz = pinfo.spacings[2];
    for (int zi = -num_ghost; zi < pinfo.ns[2] + num_ghost; zi++) {
      double z = pinfo.starts[2] + 0.5 * dz + zi * dz;
      for (int yi = -num_ghost; yi < pinfo.ns[1] + num_ghost; yi++) {
        double y = pinfo.starts[1] + 0.5 * dy + yi * dy;
        for (int xi = -num_ghost; xi < pinfo.ns[0] + num_ghost; xi++) {
          double x = pinfo.starts[0] + 0.5 * dx + xi * dx;
          ld(xi, yi, zi) = func(x, y, z);
        }
      }
    }
  }
}

void
DomainTools::SetValuesWithGhost(const Domain<2>& domain,
                                Vector<2>& vec,
                                int component_index,
                                std::function<double(double, double)> func)
{
  if (component_index >= vec.getNumComponents()) {
    throw RuntimeError("Invalid component to set");
  }
  for (int i = 0; i < vec.getNumLocalPatches(); i++) {
    ComponentView<double, 2> ld = vec.getComponentView(component_index, i);
    const PatchInfo<2>& pinfo = domain.getPatchInfoVector()[i];
    int num_ghost = pinfo.num_ghost_cells;
    double dx = pinfo.spacings[0];
    double dy = pinfo.spacings[1];
    for (int yi = -num_ghost; yi < pinfo.ns[1] + num_ghost; yi++) {
      double y = pinfo.starts[1] + 0.5 * dy + yi * dy;
      for (int xi = -num_ghost; xi < pinfo.ns[0] + num_ghost; xi++) {
        double x = pinfo.starts[0] + 0.5 * dx + xi * dx;
        ld(xi, yi) = func(x, y);
      }
    }
  }
}

template<int D>
  requires is_supported_dimension<D>
double
DomainTools::Integrate(const Domain<D>& domain, const Vector<D>& u)
{
  double sum = 0;

  for (const auto& pinfo : domain.getPatchInfoVector()) {
    for (int c = 0; c < u.getNumComponents(); c++) {
      ComponentView<const double, D> u_data = u.getComponentView(c, pinfo.local_index);

      double patch_sum = 0;
      Loop::Nested<D>(u_data.getStart(), u_data.getEnd(), [&](std::array<int, D> coord) {
        patch_sum += u_data[coord];
      });

      for (size_t i = 0; i < D; i++) {
        patch_sum *= pinfo.spacings[i];
      }
      sum += patch_sum;
    }
  }
  double retval;
  MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return retval;
}

template void
DomainTools::GetRealCoord<2>(const PatchInfo<2>& pinfo,
                             const std::array<int, 2>& coord,
                             std::array<double, 2>& real_coord);

template void
DomainTools::GetRealCoord<3>(const PatchInfo<3>& pinfo,
                             const std::array<int, 3>& coord,
                             std::array<double, 3>& real_coord);

template void
DomainTools::GetRealCoordGhost<2>(const PatchInfo<2>& pinfo,
                                  const std::array<int, 2>& coord,
                                  std::array<double, 2>& real_coord);

template void
DomainTools::GetRealCoordGhost<3>(const PatchInfo<3>& pinfo,
                                  const std::array<int, 3>& coord,
                                  std::array<double, 3>& real_coord);

template void
DomainTools::GetRealCoordBound<2>(const PatchInfo<2>& pinfo,
                                  const std::array<int, 1>& coord,
                                  Side<2> s,
                                  std::array<double, 2>& real_coord);

template void
DomainTools::GetRealCoordBound<3>(const PatchInfo<3>& pinfo,
                                  const std::array<int, 2>& coord,
                                  Side<3> s,
                                  std::array<double, 3>& real_coord);

template void
DomainTools::SetValues<2>(const Domain<2>& domain,
                          Vector<2>& vec,
                          int component_index,
                          std::function<double(const std::array<double, (int)2>&)> func);

template void
DomainTools::SetValues<3>(const Domain<3>& domain,
                          Vector<3>& vec,
                          int component_index,
                          std::function<double(const std::array<double, (int)3>&)> func);

template void
DomainTools::SetValuesWithGhost<2>(const Domain<2>& domain,
                                   Vector<2>& vec,
                                   int component_index,
                                   std::function<double(const std::array<double, (int)2>&)> func);

template void
DomainTools::SetValuesWithGhost<3>(const Domain<3>& domain,
                                   Vector<3>& vec,
                                   int component_index,
                                   std::function<double(const std::array<double, (int)3>&)> func);

template double
DomainTools::Integrate<2>(const Domain<2>& domain, const Vector<2>& u);
template double
DomainTools::Integrate<3>(const Domain<3>& domain, const Vector<3>& u);

} // namespace ThunderEgg