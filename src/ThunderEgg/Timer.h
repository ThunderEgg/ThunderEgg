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
#include <chrono>
#include <deque>
#include <functional>
#include <list>
#include <map>
#include <ostream>
#include <string>
#include <vector>
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
	 * @brief Simple structure for keeping track of a timing
	 */
	class Timing
	{
		public:
		/**
		 * @brief Pointer to parent timer
		 */
		const Timing *parent = nullptr;
		/**
		 * @brief The domain id of the timing
		 *
		 * If no domain is associated, it is set to the max value of int
		 */
		int domain_id = std::numeric_limits<int>::max();
		/**
		 * @brief The name of the timing
		 */
		std::string name;
		/**
		 * @brief The number of calls for this timing
		 */
		size_t num_calls = 0;
		/**
		 * @brief Minimum time
		 */
		double max = std::numeric_limits<double>::min();
		/**
		 * @brief Maximum time
		 */
		double min = std::numeric_limits<double>::max();
		/**
		 * @brief Sum of all timings
		 */
		double sum = 0;
		/**
		 * @brief A list of timings that are nested in this timing
		 */
		std::list<Timing> timings;
		/**
		 * @brief A map from timing domain id and name to a reference of the nested timing
		 */
		std::map<std::tuple<int, std::string>, std::reference_wrapper<Timing>> timing_map;
		/**
		 * @brief The starting time of the latest timing
		 */
		std::chrono::steady_clock::time_point start_time;
		/**
		 * @brief Construct a new Timing object
		 */
		Timing() = default;
		/**
		 * @brief Construct a new Timing object
		 *
		 * @param parent pointer to parent timing
		 * @param domain_id the id of the domain associated with the timing
		 * @param name the name of the timing
		 */
		Timing(const Timing *parent, int domain_id, const std::string &name);
		/**
		 * @brief get a Timing that is nested in this timing
		 *
		 * @param domain_id the id of the domain associated with the timing
		 * @param name  the name of the timing
		 * @return Timing& A reference to the timing
		 */
		Timing &getTiming(int domain_id, const std::string &name);
		/**
		 * @brief Start a timing
		 */
		void start();
		/**
		 * @brief stop a timing
		 */
		void stop();
		/**
		 * @brief Print a the results of this timing and the results of the timings that are nested
		 * in this timing.
		 *
		 * @param parent_string the string of timings that this timing is nested in
		 * @param os the stream
		 */
		void        print(const std::string &parent_string, std::ostream &os) const;
		friend void to_json(nlohmann::json &j, const Timing &timing);
	};
	friend void to_json(nlohmann::json &j, const Timing &timing);
	/**
	 * @brief the root timing, this is not really a timing itself, it just contains other timings
	 */
	Timing root;
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
	Timer();
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
	 * @brief ostream operator for Timer
	 *
	 * @param os the stream
	 * @param timer the timer
	 * @return std::ostream& the stream
	 */
	friend std::ostream &operator<<(std::ostream &os, const Timer &timer);
	friend void          to_json(nlohmann::json &j, const Timer &timer);
};
} // namespace ThunderEgg
#endif
