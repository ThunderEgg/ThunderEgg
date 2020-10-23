/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_TIMER_H
#define THUNDEREGG_TIMER_H
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/tpl/json.hpp>
#include <list>
#include <map>
#include <mpi.h>
#include <ostream>
#include <string>
namespace ThunderEgg
{
/**
 * @brief Class for keeping track of parallel timings
 *
 * This class will keep track of timings and present them in a nice format.
 *
 * This class also keeps track of nested timings. For following
 *
 * 		timer.start("A");
 * 		timer.start("B");
 * 		timer.stop("B");
 * 		timer.stop("A");
 *
 * outputs
 *
 * 		A
 * 		-
 * 		time (sec): 1.0
 *
 *
 * 		A -> B
 * 		------
 * 		time(sec):
 */
class Timer
{
	private:
	/**
	 * @brief The mpi communicator
	 */
	MPI_Comm comm;
	/**
	 * @brief Simple structure for keeping track of a timing
	 */
	class Timing;
	/**
	 * @brief the root timing, this is not really a timing itself, it just contains other timings
	 */
	std::unique_ptr<Timing> root;
	/**
	 * @brief The stack that keeps track of what timing we are one. Each sequential timing is nested
	 * in the other.
	 */
	std::list<std::reference_wrapper<Timing>> stack;
	std::map<int, nlohmann::json>             domains;

	public:
	/**
	 * @brief Construct a new empty Timer object
	 */
	Timer(MPI_Comm comm);
	/**
	 * @brief Destruct a Timer object
	 */
	~Timer();
	/**
	 * @brief Start a new timing
	 *
	 * @param name the name of the timing
	 */
	void start(const std::string &name);
	/**
	 * @brief Stop a timing
	 *
	 * @param name the name of the timing
	 *
	 * @exception TimerException if the name does not match the name of the last started timing.
	 */
	void stop(const std::string &name);
	/**
	 * @brief add a domain to to timer
	 *
	 * @param domain_id the id of the domain
	 * @param domain the domain
	 * @exception RuntimerError if domain with same id was already added
	 */
	void addDomain(int domain_id, nlohmann::json domain);
	/**
	 * @brief Start a new Domain associated timing
	 *
	 * @param domain_id the id of the Domain
	 * @param name the name of the timing
	 * @exception RuntimerError if domain was not added with addDomain
	 */
	void startDomainTiming(int domain_id, const std::string &name);
	/**
	 * @brief Stop a Domain associated timing
	 *
	 * @param domain_id the id of the Domain
	 * @param name the name of the timing
	 *
	 * @exception RuntimeError if the domain id and name does not match the name of the last
	 * started timing.
	 */
	void stopDomainTiming(int domain_id, const std::string &name);
	/**
	 * @brief ostream operator for Timer, this is collective for all ranks, will only output on rank
	 * 0
	 *
	 * @param os the stream
	 * @param timer the timer
	 * @return std::ostream& the stream
	 */
	friend std::ostream &operator<<(std::ostream &os, const Timer &timer);
	/**
	 * @brief Convert a timer to a json serialization, this is collective over all processes will
	 * only result in a json object for rank 0, will be null for all other ranks
	 *
	 * @param j resulting json
	 * @param timer the timer
	 */
	friend void to_json(nlohmann::json &j, const Timer &timer);
	/**
	 * @brief Save a json representation of the timer to the file. This is collective over all
	 * processes.
	 *
	 * @param filename the file to save to
	 * @exception RuntimeError on rank 0 if the file cannot be opened for writing.
	 */
	void saveToFile(const std::string &filename) const;
};
} // namespace ThunderEgg
#endif
