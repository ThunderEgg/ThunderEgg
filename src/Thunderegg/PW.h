/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
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

#ifndef THUNDEREGG_PW_H
#define THUNDEREGG_PW_H
#include <petscsys.h>
namespace Thunderegg
{
/**
 * @brief This is a shared pointer class for petsc objects
 *
 * @tparam X the petsc object type
 */
template <typename X> class PW
{
	protected:
	X       obj       = nullptr;
	size_t *ref_count = nullptr;
	bool    owns      = true;
	void    decrement()
	{
		if (ref_count != nullptr) {
			--(*ref_count);
			if (*ref_count == 0) {
				if (obj != nullptr && owns) { PetscObjectDestroy((PetscObject *) &obj); }
				delete ref_count;
			}
		}
	}

	public:
	/**
	 * @brief Construct a new PW object
	 *
	 * @param obj the petsc object
	 * @param owns whether this should deallocate the petsc object
	 */
	explicit PW(X obj = nullptr, bool owns = true)
	{
		this->obj  = obj;
		this->owns = owns;
		ref_count  = new size_t;
		*ref_count = 1;
	}
	/**
	 * @brief Copy constructor
	 */
	PW(const PW<X> &other)
	{
		ref_count = other.ref_count;
		obj       = other.obj;
		owns      = other.owns;
		++(*ref_count);
	}
	/**
	 * @brief Destroy the PW object
	 */
	~PW()
	{
		decrement();
	}
	/**
	 * @brief Assignment operator
	 */
	PW &operator=(const PW<X> &other)
	{
		decrement();
		ref_count = other.ref_count;
		obj       = other.obj;
		++(*ref_count);
		return *this;
	}
	/**
	 * @brief reset the pointer with a new petsc object
	 *
	 * @param obj the petsc object
	 * @param owns whether this pointer is responsible for destroying the petsc object
	 */
	void reset(const X &obj, bool owns)
	{
		decrement();
		this->obj  = obj;
		this->owns = owns;
		ref_count  = new size_t;
		*ref_count = 1;
	}
	/**
	 * @brief return wether this object contains a pointer
	 */
	bool isSet()
	{
		return obj != nullptr;
	}
	/*
	PW &operator=(const X &other)
	{
	    decrement();
	    ref_count    = new size_t;
	    obj          = other;
	    owns         = false;
	    (*ref_count) = 1;
	    return *this;
	} */
	/**
	 * @brief cast to petsc object
	 *
	 * @return X the petsc object
	 */
	operator X() const
	{
		return obj;
	}
	/**
	 * @brief the pointer to the petsc object
	 */
	X *operator&()
	{
		return &obj;
	}
};
/**
 * @brief This is PW object that does not have implicit cast to the petsc object.
 *
 * This is useful for return types from functions
 *
 * @tparam X the petsc object type
 */
template <typename X> class PW_explicit : public PW<X>
{
	public:
	PW_explicit() = delete;
	PW_explicit(const PW_explicit<X> &other) : PW<X>(other) {}
	PW_explicit(const PW<X> &other) : PW<X>(other) {}
	operator X() = delete;
};
} // namespace Thunderegg
#endif
