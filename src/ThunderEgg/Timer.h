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
#include <chrono>
#include <deque>
#include <functional>
#include <list>
#include <map>
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
	struct Timing {
		/**
		 * @brief The name of the timing
		 */
		std::string name;
		/**
		 * @brief A list of times
		 */
		std::deque<double> times;
		/**
		 * @brief A list of timings that are neseted in this timing
		 */
		std::list<Timing> timings;
		/**
		 * @brief A map from timing name to a reference of the nested timing
		 */
		std::map<std::string, std::reference_wrapper<Timing>> timing_map;
		/**
		 * @brief The starting time of the latest timing
		 */
		std::chrono::steady_clock::time_point start;
		/**
		 * @brief Construct a new Timing object
		 */
		Timing() = default;
		/**
		 * @brief Construct a new Timing object
		 *
		 * @param name the name of the timing
		 */
		Timing(std::string name)
		{
			this->name = name;
		}
		/**
		 * @brief get a Timing that is nested in this timing
		 *
		 * @param name  the name of the timing
		 * @return Timing& A reference to the timing
		 */
		Timing &getTiming(std::string name)
		{
			auto timing_map_iter = timing_map.find(name);
			if (timing_map_iter == timing_map.end()) {
				timings.push_back(Timing(name));
				timing_map.emplace(name, timings.back());
				return timings.back();
			} else {
				return timing_map_iter->second;
			}
		}
		/**
		 * @brief Print a the results of this timing and the results of the timings that are nested
		 * in this timing.
		 *
		 * @param parent_string the string of timings that this timing is nested in
		 * @param os the stream
		 */
		void print(std::string parent_string, std::ostream &os) const
		{
			os << parent_string + name << std::endl;
			os << std::string((parent_string + name).size(), '-') << std::endl;

			if (times.size() == 1) {
				os << "   time (sec): " << times[0] << std::endl << std::endl;
			} else {
				os << "  total calls: " << times.size() << std::endl;
				double average = 0;
				for (double t : times) {
					average += t;
				}
				average /= times.size();

				os << "average (sec): " << average << std::endl;
				if (times.size() < 10) {
					os << "  times (sec):";
					for (double t : times) {
						os << " " << t;
					}
					os << std::endl;
				}
				os << std::endl;
			}
			for (const Timing &timing : timings) {
				timing.print(parent_string + name + " -> ", os);
			}
		}
	};
	/**
	 * @brief the root timing, this is not really a timing itself, it just contains other timings
	 */
	Timing root;
	/**
	 * @brief The stack that keeps track of what timing we are one. Each sequential timing is nested
	 * in the other.
	 */
	std::list<std::reference_wrapper<Timing>> stack;

	public:
	/**
	 * @brief Construct a new empty Timer object
	 */
	Timer()
	{
		stack.push_back(root);
	}
	/**
	 * @brief Start a new timing
	 *
	 * @param name the name of the timing
	 */
	void start(std::string name)
	{
		Timing &curr_timing = stack.back();
		Timing &next_timing = curr_timing.getTiming(name);
		next_timing.start   = std::chrono::steady_clock::now();
		stack.push_back(next_timing);
	}
	/**
	 * @brief Stop a timing
	 *
	 * @param name the name of the timing
	 *
	 * @exception TimerException if the name does not match the name of the last started timing.
	 */
	void stop(std::string name)
	{
		Timing &curr_timing = stack.back();
		if (curr_timing.name == name) {
			std::chrono::duration<double> time
			= std::chrono::steady_clock::now() - curr_timing.start;
			curr_timing.times.push_back(time.count());
			stack.pop_back();
		} else {
			throw RuntimeError("Timer was expecting to end \"" + curr_timing.name
			                   + "\", instead got \"" + name + "\"");
		}
	}
	/**
	 * @brief ostream operator for Timer
	 *
	 * @param os the stream
	 * @param timer the timer
	 * @return std::ostream& the stream
	 */
	friend std::ostream &operator<<(std::ostream &os, const Timer &timer)
	{
		timer.stack.empty();
		os << std::endl;
		os << "TIMING RESULTS" << std::endl;
		os << "==============" << std::endl << std::endl;

		for (const Timing &timing : timer.root.timings) {
			timing.print("", os);
		}
		return os;
	}
};
} // namespace ThunderEgg
#endif
