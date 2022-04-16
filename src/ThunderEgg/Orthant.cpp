/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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

#include <ThunderEgg/Orthant.h>
#include <ThunderEgg/tpl/json.hpp>

namespace ThunderEgg {

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<D>::Orthant(const unsigned char val_in)
  : val(val_in)
{
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<D>::Orthant()
{
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<D>
Orthant<D>::null()
{
  return Orthant<D>(num_orthants);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<1>
Orthant<D>::lower()
{
  return Orthant<1>(0b0);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<1>
Orthant<D>::upper()
{
  return Orthant<1>(0b1);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<2>
Orthant<D>::sw()
{
  return Orthant<2>(0b00);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<2>
Orthant<D>::se()
{
  return Orthant<2>(0b01);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<2>
Orthant<D>::nw()
{
  return Orthant<2>(0b10);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<2>
Orthant<D>::ne()
{
  return Orthant<2>(0b11);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<3>
Orthant<D>::bsw()
{
  return Orthant<3>(0b000);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<3>
Orthant<D>::bse()
{
  return Orthant<3>(0b001);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<3>
Orthant<D>::bnw()
{
  return Orthant<3>(0b010);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<3>
Orthant<D>::bne()
{
  return Orthant<3>(0b011);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<3>
Orthant<D>::tsw()
{
  return Orthant<3>(0b100);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<3>
Orthant<D>::tse()
{
  return Orthant<3>(0b101);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<3>
Orthant<D>::tnw()
{
  return Orthant<3>(0b110);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<3>
Orthant<D>::tne()
{
  return Orthant<3>(0b111);
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<D>::Range::Iterator::Iterator(Orthant<D> o_in)
  : o(o_in)
{
}

template<int D>
  requires is_supported_orthant_dimension<D>
const Orthant<D>&
Orthant<D>::Range::Iterator::operator++()
{
  ++o.val;
  return o;
}

template<int D>
  requires is_supported_orthant_dimension<D>
const Orthant<D>&
Orthant<D>::Range::Iterator::operator*() const
{
  return o;
}

template<int D>
  requires is_supported_orthant_dimension<D>
const Orthant<D>*
Orthant<D>::Range::Iterator::operator->() const
{
  return &o;
}

template<int D>
  requires is_supported_orthant_dimension<D> bool
Orthant<D>::Range::Iterator::operator==(const Iterator& b) const
{
  return o.val == b.o.val;
}

template<int D>
  requires is_supported_orthant_dimension<D> bool
Orthant<D>::Range::Iterator::operator!=(const Iterator& b) const
{
  return o.val != b.o.val;
}

template<int D>
  requires is_supported_orthant_dimension<D>
typename Orthant<D>::Range::Iterator
Orthant<D>::Range::begin()
{
  return Iterator(Orthant<D>(0));
}

template<int D>
  requires is_supported_orthant_dimension<D>
typename Orthant<D>::Range::Iterator
Orthant<D>::Range::end()
{
  return Iterator(null());
}

template<int D>
  requires is_supported_orthant_dimension<D>
typename Orthant<D>::Range
Orthant<D>::getValues()
{
  return Range();
}

template<int D>
  requires is_supported_orthant_dimension<D>
size_t
Orthant<D>::getIndex() const
{
  return val;
}

template<int D>
  requires is_supported_orthant_dimension<D>
Orthant<D>
Orthant<D>::getNbrOnSide(Side<D> s) const
{
  Orthant<D> retval = *this;
  // flip the bit for that side
  retval.val ^= (0x1 << s.getAxisIndex());
  return retval;
}

template<int D>
  requires is_supported_orthant_dimension<D>
std::array<Side<D>, D>
Orthant<D>::getInteriorSides() const
{
  std::array<Side<D>, D> retval;
  for (size_t i = 0; i < D; i++) {
    size_t side = 2 * i;
    if (!((1 << i) & val)) {
      side |= 1;
    }
    retval[i] = Side<D>(side);
  }
  return retval;
}

template<int D>
  requires is_supported_orthant_dimension<D>
std::array<Side<D>, D>
Orthant<D>::getExteriorSides() const
{
  std::array<Side<D>, D> retval;
  for (size_t i = 0; i < D; i++) {
    size_t side = 2 * i;
    if ((1 << i) & val) {
      side |= 1;
    }
    retval[i] = Side<D>(side);
  }
  return retval;
}

template<int D>
  requires is_supported_orthant_dimension<D> bool
Orthant<D>::isOnSide(Side<D> s) const
{
  int idx = s.getIndex() / 2;
  int remainder = s.getIndex() % 2;
  bool is_bit_set = val & (0x1 << idx);
  return is_bit_set == remainder;
}

template<int D>
  requires is_supported_orthant_dimension<D> bool
Orthant<D>::isHigherOnAxis(size_t axis) const
{
  return val & (0b1 << axis);
}

template<int D>
  requires is_supported_orthant_dimension<D> bool
Orthant<D>::isLowerOnAxis(size_t axis) const
{
  return !(val & (0b1 << axis));
}

template<int D>
  requires is_supported_orthant_dimension<D>
template<int Dm1>
  requires(Dm1 == D - 1) && (D - 1 >= 0)
Orthant<Dm1> Orthant<D>::collapseOnAxis(size_t axis) const
{
  size_t upper_mask = (~0x0U) << axis;
  return Orthant<D - 1>(((val >> 1) & upper_mask) | (val & ~upper_mask));
}

template<int D>
  requires is_supported_orthant_dimension<D>
std::array<Orthant<D>, Orthant<D>::num_orthants / 2>
Orthant<D>::getValuesOnSide(Side<D> s)
{
  unsigned int bit_to_insert = s.getAxisIndex();
  unsigned int set_bit = s.isLowerOnAxis() ? 0 : 1;
  unsigned int lower_mask = ~((~0x0U) << bit_to_insert);
  unsigned int upper_mask = (~0x0U) << (bit_to_insert + 1);

  std::array<Orthant<D>, Orthant<D>::num_orthants / 2> retval;
  for (size_t i = 0; i < Orthant<D>::num_orthants / 2; i++) {
    size_t value = (i << 1) & upper_mask;
    value |= i & lower_mask;
    value |= set_bit << bit_to_insert;
    retval[i] = Orthant<D>(value);
  }
  return retval;
}

template<int D>
  requires is_supported_orthant_dimension<D> bool
Orthant<D>::operator==(const Orthant<D>& other) const
{
  return val == other.val;
}

template<int D>
  requires is_supported_orthant_dimension<D> bool
Orthant<D>::operator!=(const Orthant<D>& other) const
{
  return val != other.val;
}

template<int D>
  requires is_supported_orthant_dimension<D> bool
Orthant<D>::operator<(const Orthant<D>& other) const
{
  return val < other.val;
}

std::ostream&
operator<<(std::ostream& os, const Orthant<0>& o)
{
  if (o == Orthant<0>::null()) {
    os << "Orthant<0>::null()";
  } else {
    os << "Orthant<0> invalid value: " << o.getIndex();
  }
  return os;
}

std::ostream&
operator<<(std::ostream& os, const Orthant<1>& o)
{
  if (o == Orthant<1>::lower()) {
    os << "Orthant<1>::lower()";
  } else if (o == Orthant<1>::upper()) {
    os << "Orthant<1>::upper()";
  } else if (o == Orthant<1>::null()) {
    os << "Orthant<1>::null()";
  } else {
    os << "Orthant<1> invalid value: " << o.getIndex();
  }
  return os;
}

std::ostream&
operator<<(std::ostream& os, const Orthant<2>& o)
{
  if (o == Orthant<2>::sw()) {
    os << "Orthant<2>::sw()";
  } else if (o == Orthant<2>::se()) {
    os << "Orthant<2>::se()";
  } else if (o == Orthant<2>::nw()) {
    os << "Orthant<2>::nw()";
  } else if (o == Orthant<2>::ne()) {
    os << "Orthant<2>::ne()";
  } else if (o == Orthant<2>::null()) {
    os << "Orthant<2>::null()";
  } else {
    os << "Orthant<2> invalid value: " << o.getIndex();
  }
  return os;
}

std::ostream&
operator<<(std::ostream& os, const Orthant<3>& o)
{
  if (o == Orthant<3>::bsw()) {
    os << "Orthant<3>::bsw()";
  } else if (o == Orthant<3>::bse()) {
    os << "Orthant<3>::bse()";
  } else if (o == Orthant<3>::bnw()) {
    os << "Orthant<3>::bnw()";
  } else if (o == Orthant<3>::bne()) {
    os << "Orthant<3>::bne()";
  } else if (o == Orthant<3>::tsw()) {
    os << "Orthant<3>::tsw()";
  } else if (o == Orthant<3>::tse()) {
    os << "Orthant<3>::tse()";
  } else if (o == Orthant<3>::tnw()) {
    os << "Orthant<3>::tnw()";
  } else if (o == Orthant<3>::tne()) {
    os << "Orthant<3>::tne()";
  } else if (o == Orthant<3>::null()) {
    os << "Orthant<3>::null()";
  } else {
    os << "Orthant<3> invalid value: " << o.getIndex();
  }
  return os;
}

void
to_json(tpl::nlohmann::json& j, const Orthant<0>&)
{
  j = nullptr;
}

void
to_json(tpl::nlohmann::json& j, const Orthant<1>& o)
{
  if (o == Orthant<1>::lower()) {
    j = "LOWER";
  } else if (o == Orthant<1>::upper()) {
    j = "UPPER";
  } else {
    j = nullptr;
  }
}

void
to_json(tpl::nlohmann::json& j, const Orthant<2>& o)
{
  if (o == Orthant<2>::sw()) {
    j = "SW";
  } else if (o == Orthant<2>::se()) {
    j = "SE";
  } else if (o == Orthant<2>::nw()) {
    j = "NW";
  } else if (o == Orthant<2>::ne()) {
    j = "NE";
  } else {
    j = nullptr;
  }
}

void
to_json(tpl::nlohmann::json& j, const Orthant<3>& o)
{
  if (o == Orthant<3>::bsw()) {
    j = "BSW";
  } else if (o == Orthant<3>::bse()) {
    j = "BSE";
  } else if (o == Orthant<3>::bnw()) {
    j = "BNW";
  } else if (o == Orthant<3>::bne()) {
    j = "BNE";
  } else if (o == Orthant<3>::tsw()) {
    j = "TSW";
  } else if (o == Orthant<3>::tse()) {
    j = "TSE";
  } else if (o == Orthant<3>::tnw()) {
    j = "TNW";
  } else if (o == Orthant<3>::tne()) {
    j = "TNE";
  } else {
    j = nullptr;
  }
}

void
from_json(const tpl::nlohmann::json&, Orthant<0>& o)
{
  o = Orthant<0>::null();
}

void
from_json(const tpl::nlohmann::json& j, Orthant<1>& o)
{
  if (j == "LOWER") {
    o = Orthant<1>::lower();
  } else if (j == "UPPER") {
    o = Orthant<1>::upper();
  } else {
    o = Orthant<1>::null();
  }
}

void
from_json(const tpl::nlohmann::json& j, Orthant<2>& o)
{
  if (j == "SW") {
    o = Orthant<2>::sw();
  } else if (j == "SE") {
    o = Orthant<2>::se();
  } else if (j == "NW") {
    o = Orthant<2>::nw();
  } else if (j == "NE") {
    o = Orthant<2>::ne();
  } else {
    o = Orthant<2>::null();
  }
}

void
from_json(const tpl::nlohmann::json& j, Orthant<3>& o)
{
  if (j == "BSW") {
    o = Orthant<3>::bsw();
  } else if (j == "BSE") {
    o = Orthant<3>::bse();
  } else if (j == "BNW") {
    o = Orthant<3>::bnw();
  } else if (j == "BNE") {
    o = Orthant<3>::bne();
  } else if (j == "TSW") {
    o = Orthant<3>::tsw();
  } else if (j == "TSE") {
    o = Orthant<3>::tse();
  } else if (j == "TNW") {
    o = Orthant<3>::tnw();
  } else if (j == "TNE") {
    o = Orthant<3>::tne();
  } else {
    o = Orthant<3>::null();
  }
}

// EXPLICIT INSTANTIATIONS

template class Orthant<0>;

template class Orthant<1>;

template Orthant<0>
Orthant<1>::collapseOnAxis(size_t axis) const;

template class Orthant<2>;

template Orthant<1>
Orthant<2>::collapseOnAxis(size_t axis) const;

template class Orthant<3>;

template Orthant<2>
Orthant<3>::collapseOnAxis(size_t axis) const;

} // namespace ThunderEgg
