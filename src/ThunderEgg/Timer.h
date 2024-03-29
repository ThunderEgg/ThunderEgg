/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2017-2021 Scott Aiton
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
/**
 * @file
 *
 * @brief Timer class
 */

#include <ThunderEgg/Communicator.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/tpl/json_fwd.hpp>
#include <list>
#include <map>
#include <ostream>
#include <string>
namespace ThunderEgg {
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
  Communicator comm;
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
  std::map<int, tpl::nlohmann::json> domains;

public:
  /**
   * @brief Construct a new empty Timer object
   */
  explicit Timer(const Communicator& comm);
  /**
   * @brief Destruct a Timer object
   */
  ~Timer();
  /**
   * @brief Start a new timing
   *
   * @param name the name of the timing
   */
  void start(const std::string& name);
  /**
   * @brief Stop a timing
   *
   * @param name the name of the timing
   *
   * @exception TimerException if the name does not match the name of the last started timing.
   */
  void stop(const std::string& name);
  /**
   * @brief add a domain to to timer
   *
   * @param domain_id the id of the domain
   * @param domain the domain
   * @exception RuntimerError if domain with same id was already added
   */
  void addDomain(int domain_id, tpl::nlohmann::json domain);
  /**
   * @brief Start a new Domain associated timing
   *
   * @param domain_id the id of the Domain
   * @param name the name of the timing
   * @exception RuntimerError if domain was not added with addDomain
   */
  void startDomainTiming(int domain_id, const std::string& name);
  /**
   * @brief Stop a Domain associated timing
   *
   * @param domain_id the id of the Domain
   * @param name the name of the timing
   *
   * @exception RuntimeError if the domain id and name does not match the name of the last
   * started timing.
   */
  void stopDomainTiming(int domain_id, const std::string& name);
  /**
   * @brief Start a new Domain associated timing
   *
   * @param patch_id the id of the PatchInfo
   * @param domain_id the id of the Domain
   * @param name the name of the timing
   * @exception RuntimerError if domain was not added with addDomain
   */
  void startPatchTiming(int patch_id, int domain_id, const std::string& name);
  /**
   * @brief Stop a Domain associated timing
   *
   * @param patch_id the id of the PatchInfo
   * @param domain_id the id of the Domain
   * @param name the name of the timing
   *
   * @exception RuntimeError if the patch id, domain id and name does not match the name of the
   * last started timing.
   */
  void stopPatchTiming(int patch_id, int domain_id, const std::string& name);
  /**
   * @brief Add information to a timing
   *
   * Has to be called after start is called for the timing and before stop is called for the
   * timing
   *
   * @param name the name of the information
   * @param info the value of the information
   *
   * @exception RuntimeError there is no timing to add information to, or if adding int
   * information to existing double information
   */
  void addIntInfo(const std::string& name, int info);
  /**
   * @brief Add information to a timing
   *
   * Has to be called after start is called for the timing and before stop is called for the
   * timing
   *
   * @param name the name of the information
   * @param info the value of the information
   *
   * @exception RuntimeError there is no timing to add information to, or if adding double
   * information to existing int information
   */
  void addDoubleInfo(const std::string& name, double info);
  /**
   * @brief ostream operator for Timer, this is collective for all ranks, will only output on rank
   * 0
   *
   * @param os the stream
   * @param timer the timer
   * @return std::ostream& the stream
   */
  friend std::ostream& operator<<(std::ostream& os, const Timer& timer);
  /**
   * @brief Convert a timer to a json serialization, this is collective over all processes will
   * only result in a json object for rank 0, will be null for all other ranks
   *
   * @param j resulting json
   * @param timer the timer
   */
  friend void to_json(tpl::nlohmann::json& j, const Timer& timer);
  /**
   * @brief Save a json representation of the timer to the file. This is collective over all
   * processes.
   *
   * @param filename the file to save to
   * @exception RuntimeError on rank 0 if the file cannot be opened for writing.
   */
  void saveToFile(const std::string& filename) const;
};
} // namespace ThunderEgg
#endif
