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

#include <ThunderEgg/Timer.h>
#include <chrono>
#include <fstream>
#include <iomanip>

using namespace ThunderEgg::tpl;
namespace {
/**
 * @brief Stores scalar information associated with a timing
 */
class Info
{
public:
  /**
   * @brief Destroy the Info object
   */
  virtual ~Info() {}
  /**
   * @brief serialize Info as json object
   *
   * @param j json object to add to
   */
  virtual void to_json(nlohmann::json& j) = 0;
};
/**
 * @brief Stores integer information
 */
class IntInfo : public Info
{
private:
  /**
   * @brief name of the information
   */
  std::string name;
  /**
   * @brief Sum of all informations
   */
  long int sum = 0;
  /**
   * @brief Min value
   */
  int min = std::numeric_limits<int>::max();
  /**
   * @brief max value
   */
  int max = std::numeric_limits<int>::min();
  /**
   * @brief number of informations stored
   */
  int num_calls = 0;

public:
  /**
   * @brief Construct a new Int Info object
   *
   * @param name the name of the information
   */
  explicit IntInfo(const std::string& name)
    : name(name)
  {}
  /**
   * @brief Add information
   *
   * @param info the value to add
   */
  void addInfo(int info)
  {
    sum += info;
    min = std::min(min, info);
    max = std::max(max, info);
    num_calls++;
  }
  void to_json(nlohmann::json& j) override
  {
    j["name"] = name;
    j["sum"] = sum;
    j["min"] = min;
    j["max"] = max;
    j["num_calls"] = num_calls;
  }
};
/**
 * @brief Stores double information
 */
class DoubleInfo : public Info
{
private:
  /**
   * @brief name of the information
   */
  std::string name;
  /**
   * @brief sum of all informations
   */
  double sum = 0;
  /**
   * @brief min value
   */
  double min = std::numeric_limits<double>::max();
  /**
   * @brief max value
   */
  double max = std::numeric_limits<double>::lowest();
  /**
   * @brief number of informations stored
   */
  int num_calls = 0;

public:
  /**
   * @brief Construct a new Double Info object
   *
   * @param name the name of the information
   */
  explicit DoubleInfo(const std::string& name)
    : name(name)
  {}
  /**
   * @brief add information
   *
   * @param info the value to add
   */
  void addInfo(double info)
  {
    sum += info;
    min = std::min(min, info);
    max = std::max(max, info);
    num_calls++;
  }
  void to_json(nlohmann::json& j) override
  {
    j["name"] = name;
    j["sum"] = sum;
    j["min"] = min;
    j["max"] = max;
    j["num_calls"] = num_calls;
  }
};
void
to_json(nlohmann::json& j, const std::unique_ptr<Info>& info)
{
  info->to_json(j);
}
} // namespace
namespace ThunderEgg {
class Timer::Timing
{
public:
  /**
   * @brief Pointer to parent timer
   */
  const Timing* parent = nullptr;
  /**
   * @brief The name of the timing
   */
  std::string name;
  /**
   * @brief The domain id of the timing
   *
   * If no domain is associated, it is set to the max value of int
   */
  int domain_id = std::numeric_limits<int>::max();
  /**
   * @brief The patch id of the timing
   *
   * If no patch is associated, it is set to the max value of int
   */
  int patch_id = std::numeric_limits<int>::max();
  /**
   * @brief The number of calls for this timing
   */
  size_t num_calls = 0;
  /**
   * @brief Minimum time
   */
  double max = std::numeric_limits<double>::lowest();
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
   * @brief A list of informaiton associated with this timing
   */
  std::list<std::unique_ptr<Info>> infos;
  /**
   * @brief A map from timing name,domain id, and patch id to a reference of the nested timing
   */
  std::map<std::tuple<std::string, int, int>, std::reference_wrapper<Timing>> timing_map;
  /**
   * @brief A map from name to Info object
   */
  std::map<std::string, Info*> info_map;
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
   * @param patch_id the id of the patch associated with the timing
   * @param name the name of the timing
   */
  Timing(const Timing* parent, int patch_id, int domain_id, const std::string& name)
    : parent(parent)
    , name(name)
    , domain_id(domain_id)
    , patch_id(patch_id)
  {}
  /**
   * @brief get a Timing that is nested in this timing
   *
   * @param child_patch_id the id of the patch associated with the timing
   * @param child_domain_id the id of the domain associated with the timing
   * @param child_name  the name of the timing
   * @return Timing& A reference to the timing
   */
  Timing& getTiming(int child_patch_id, int child_domain_id, const std::string& child_name)
  {
    auto key = std::make_tuple(child_name, child_domain_id, child_patch_id);
    auto timing_map_iter = timing_map.find(key);
    if (timing_map_iter == timing_map.end()) {
      timings.push_back(Timing(this, child_patch_id, child_domain_id, child_name));
      timing_map.emplace(key, timings.back());
      return timings.back();
    } else {
      return timing_map_iter->second;
    }
  }
  /**
   * @brief Start a timing
   */
  void start() { start_time = std::chrono::steady_clock::now(); }
  /**
   * @brief stop a timing
   */
  void stop()
  {
    std::chrono::duration<double> duration = std::chrono::steady_clock::now() - start_time;
    double time = duration.count();
    sum += time;
    max = std::max(max, time);
    min = std::min(min, time);
    num_calls++;
  }
  /**
   * @brief Add information to this timing
   *
   * Has to be called after start is called for the timing and before stop is called for the
   * timing
   *
   * @param info_name the name of the information
   * @param info the value of the information
   *
   * @exception RuntimeError if adding int
   * information to existing double information
   */
  void addIntInfo(const std::string& info_name, int info)
  {
    auto pair = info_map.emplace(info_name, nullptr);
    IntInfo* info_ptr = dynamic_cast<IntInfo*>(pair.first->second);
    if (pair.second == true) {
      info_ptr = new IntInfo(info_name);
      infos.emplace_back(info_ptr);
      pair.first->second = info_ptr;
    } else if (info_ptr == nullptr) {
      throw RuntimeError("Adding int info to existing double info timing" + name + " for info " +
                         info_name);
    }
    info_ptr->addInfo(info);
  }
  /**
   * @brief Add information to this timing
   *
   * Has to be called after start is called for the timing and before stop is called for the
   * timing
   *
   * @param info_name the name of the information
   * @param info the value of the information
   *
   * @exception RuntimeError if adding double
   * information to existing int information
   */
  void addDoubleInfo(const std::string& info_name, double info)
  {
    auto pair = info_map.emplace(info_name, nullptr);
    DoubleInfo* info_ptr = dynamic_cast<DoubleInfo*>(pair.first->second);
    if (pair.second == true) {
      info_ptr = new DoubleInfo(info_name);
      infos.emplace_back(info_ptr);
      pair.first->second = info_ptr;
    } else if (info_ptr == nullptr) {
      throw RuntimeError("Adding int info to existing int info timing" + name + " for info " +
                         info_name);
    }
    info_ptr->addInfo(info);
  }
  friend void to_json(nlohmann::json& j, const Timing& timing)
  {
    if (timing.name != "") {
      j["name"] = timing.name;
    }
    if (timing.domain_id != std::numeric_limits<int>::max()) {
      j["domain_id"] = timing.domain_id;
    }
    if (timing.patch_id != std::numeric_limits<int>::max()) {
      j["patch_id"] = timing.patch_id;
    }
    if (timing.num_calls != 0) {
      j["num_calls"] = timing.num_calls;
      j["sum"] = timing.sum;
      j["max"] = timing.max;
      j["min"] = timing.min;
    }
    if (timing.timings.size() > 0) {
      j["timings"] = timing.timings;
    }
    if (timing.infos.size() > 0) {
      j["infos"] = timing.infos;
    }
  }
};
Timer::Timer(const Communicator& comm)
  : comm(comm)
  , root(new Timing())
{
  stack.push_back(*root);
}
Timer::~Timer() = default;
void
Timer::start(const std::string& name)
{
  startDomainTiming(std::numeric_limits<int>::max(), name);
}
void
Timer::stop(const std::string& name)
{
  stopDomainTiming(std::numeric_limits<int>::max(), name);
}
void
Timer::addDomain(int domain_id, nlohmann::json domain)
{
  auto pair = domains.emplace(domain_id, domain);
  if (!pair.second) {
    throw RuntimeError("Domain with id " + std::to_string(domain_id) +
                       " was already added to timer");
  }
}
void
Timer::startDomainTiming(int domain_id, const std::string& name)
{
  startPatchTiming(std::numeric_limits<int>::max(), domain_id, name);
}
void
Timer::stopDomainTiming(int domain_id, const std::string& name)
{
  stopPatchTiming(std::numeric_limits<int>::max(), domain_id, name);
}
void
Timer::startPatchTiming(int patch_id, int domain_id, const std::string& name)
{
  if (domain_id != std::numeric_limits<int>::max() && domains.find(domain_id) == domains.end()) {
    throw RuntimeError("Domain with id " + std::to_string(domain_id) + " was not added to timer");
  }
  Timing& curr_timing = stack.back();
  Timing& next_timing = curr_timing.getTiming(patch_id, domain_id, name);
  next_timing.start();
  stack.push_back(next_timing);
}
void
Timer::stopPatchTiming(int patch_id, int domain_id, const std::string& name)
{
  Timing& curr_timing = stack.back();
  if (curr_timing.patch_id == patch_id && curr_timing.domain_id == domain_id &&
      curr_timing.name == name && stack.size() > 1) {
    curr_timing.stop();
    stack.pop_back();
  } else {
    throw RuntimeError("Timer was expecting to end \"" + curr_timing.name + "\", instead got \"" +
                       name + "\"");
  }
}

void
Timer::addIntInfo(const std::string& name, int info)
{
  if (stack.size() == 1) {
    throw RuntimeError("No timing to add information to");
  }
  Timing& curr_timing = stack.back();
  curr_timing.addIntInfo(name, info);
}

void
Timer::addDoubleInfo(const std::string& name, double info)
{
  if (stack.size() == 1) {
    throw RuntimeError("No timing to add information to");
  }
  Timing& curr_timing = stack.back();
  curr_timing.addDoubleInfo(name, info);
}
static void
PrintMergedTimings(const Communicator& comm,
                   size_t max_name_size,
                   const std::string& parent_string,
                   std::ostream& os,
                   nlohmann::json& timings);
/**
 * @brief Output an information line
 *
 * @param os the ostream
 * @param max_name_size the max size of an information name
 * @param name the name
 * @param value the value
 */
static void
OutputLine(std::ostream& os, size_t max_name_size, const std::string& name, double value)
{
  os.width(max_name_size + 1);
  os << std::right << name << std::left;
  os.width(0);
  os << ": " << value << std::endl;
}

/**
 * @brief Print a timing
 *
 * @param comm the mpi comm
 * @param max_name_size the max size of a name
 * @param parent_string the parent string
 * @param os the output stream
 * @param timing the timing to print
 */
static void
PrintTiming(const Communicator& comm,
            size_t max_name_size,
            const std::string& parent_string,
            std::ostream& os,
            nlohmann::json& timing)
{
  std::string my_string = parent_string + timing["name"].get<std::string>();
  os << my_string << std::endl;
  os << std::string(my_string.size(), '-') << std::endl;

  int size = comm.getSize();

  if (timing["num_calls"].get<size_t>() == 1) {
    OutputLine(os, max_name_size, "time (sec)", timing["sum"].get<double>());
  } else {
    OutputLine(
      os, max_name_size, "average calls per rank", timing["num_calls"].get<double>() / size);
    OutputLine(os,
               max_name_size,
               "average (sec)",
               timing["sum"].get<double>() / timing["num_calls"].get<double>());
    OutputLine(os, max_name_size, "min (sec)", timing["min"].get<double>());
    OutputLine(os, max_name_size, "max (sec)", timing["max"].get<double>());
  }
  for (const auto& info : timing["infos"]) {
    if (info["num_calls"] == 1) {
      OutputLine(os, max_name_size, info["name"].get<std::string>(), info["sum"].get<double>());
    } else {
      OutputLine(os,
                 max_name_size,
                 info["name"].get<std::string>() + " avg",
                 info["sum"].get<double>() / info["num_calls"].get<double>());
      OutputLine(
        os, max_name_size, info["name"].get<std::string>() + " min", info["min"].get<double>());
      OutputLine(
        os, max_name_size, info["name"].get<std::string>() + " max", info["max"].get<double>());
    }
  }
  os << std::endl;
  PrintMergedTimings(comm, max_name_size, my_string + " -> ", os, timing["timings"]);
}
/**
 * @brief Merge one info object with another
 *
 * @param a the info object to merge to
 * @param b the info object to merge
 */
static void
MergeInfo(nlohmann::json& a, const nlohmann::json& b)
{
  a["num_calls"] = a["num_calls"].get<size_t>() + b["num_calls"].get<size_t>();
  a["sum"] = a["sum"].get<double>() + b["sum"].get<double>();
  a["min"] = std::min(a["min"].get<double>(), b["min"].get<double>());
  a["max"] = std::max(a["max"].get<double>(), b["max"].get<double>());
}
/**
 * @brief Merge together infos with same name
 *
 * @param a_infos the info array to merge to
 * @param b_infos the info array to merge
 */
static void
MergeInfos(nlohmann::json& a_infos, const nlohmann::json& b_infos)
{
  std::map<std::string, size_t> inserted_names;
  for (size_t i = 0; i < a_infos.size(); i++) {
    inserted_names[a_infos[i]["name"]] = i;
  }
  for (const nlohmann::json& info : b_infos) {
    auto pair = inserted_names.emplace(info["name"], a_infos.size());
    if (pair.second == true) {
      a_infos.push_back(info);
    } else {
      MergeInfo(a_infos[pair.first->second], info);
    }
  }
}
/**
 * @brief Merge one timing with another
 *
 * @param a the timing to merge to
 * @param b the timing to merge
 */
static void
MergeTiming(nlohmann::json& a, nlohmann::json& b)
{
  a["num_calls"] = a["num_calls"].get<size_t>() + b["num_calls"].get<size_t>();
  a["sum"] = a["sum"].get<double>() + b["sum"].get<double>();
  a["min"] = std::min(a["min"].get<double>(), b["min"].get<double>());
  a["max"] = std::max(a["max"].get<double>(), b["max"].get<double>());
  MergeInfos(a["infos"], b["infos"]);
  nlohmann::json& a_timings = a["timings"];
  for (const nlohmann::json& b_timing : b["timings"]) {
    a_timings.push_back(b_timing);
  }
}
/**
 * @brief Merge together timings with same name
 *
 * @param timings the timings array
 * @return nlohmann::json the resulting merged array
 */
static nlohmann::json
MergeTimings(nlohmann::json& timings)
{
  nlohmann::json merged_timings;
  std::map<std::string, size_t> inserted_names;
  for (nlohmann::json& timing : timings) {
    auto pair = inserted_names.emplace(timing["name"], merged_timings.size());
    if (pair.second == true) {
      merged_timings.push_back(timing);
    } else {
      MergeTiming(merged_timings[pair.first->second], timing);
    }
  }
  return merged_timings;
}
/**
 * @brief merge timings with same name together
 *
 * @param comm the mpi comm
 * @param max_name_size the max size of a name
 * @param parent_string the string of the parent timing
 * @param os the output stream
 * @param timings the timings to print
 */
static void
PrintMergedTimings(const Communicator& comm,
                   size_t max_name_size,
                   const std::string& parent_string,
                   std::ostream& os,
                   nlohmann::json& timings)
{
  nlohmann::json merged_timings = MergeTimings(timings);
  for (nlohmann::json& timing : merged_timings) {
    PrintTiming(comm, max_name_size, parent_string, os, timing);
  }
}
void
GetMaxInfoNameLength(const nlohmann::json& timings, size_t& max_name_size)
{
  for (const nlohmann::json& timing : timings) {
    if (timing.contains("infos")) {
      for (const nlohmann::json& info : timing["infos"]) {
        max_name_size = std::max(max_name_size, info["name"].size() + 4);
      }
    }
    if (timing.contains("timings")) {
      GetMaxInfoNameLength(timing["timings"], max_name_size);
    }
  }
}
std::ostream&
operator<<(std::ostream& os, const Timer& timer)
{
  if (timer.stack.size() > 1) {
    const Timer::Timing& curr_timing = timer.stack.back();
    throw RuntimeError("Cannot output Timer results with unfinished timings, check that all "
                       "timings have been stopped. Currently waiting on timing \"" +
                       curr_timing.name + "\"");
  }
  nlohmann::json timer_j = timer;

  if (timer.comm.getRank() == 0) {
    os << std::endl;
    os << "TIMING RESULTS" << std::endl;
    os << "==============" << std::endl << std::endl;

    if (timer_j["timings"] == nullptr) {
      os << "No timings to report." << std::endl << std::endl;
    } else {
      size_t max_name_size = std::strlen("average calls per rank");
      GetMaxInfoNameLength(timer_j["timings"], max_name_size);
      PrintMergedTimings(timer.comm, max_name_size, "", os, timer_j["timings"]);
    }
  }
  return os;
}
/**
 * @brief add rank value to timings array
 *
 * @param j the timing to array
 * @param rank the rank value to add
 */

static void
DecorateTimingsWithRank(nlohmann::json& timings_j, int rank)
{
  for (auto& timing : timings_j) {
    timing["rank"] = rank;
    if (timing.contains("timings")) {
      DecorateTimingsWithRank(timing["timings"], rank);
    }
  }
}
/**
 * @brief add rank value to timer object
 *
 * @param j the timing to array
 * @param rank the rank value to add
 */
static void
DecorateWithRank(nlohmann::json& j, int rank)
{
  if (j.contains("timings")) {
    DecorateTimingsWithRank(j["timings"], rank);
  }
}
/**
 * @brief merge domains
 *
 * @param j the domain array to merge to
 * @param incoming_j the domain array to merge
 */
static void
MergeIncomingDomains(nlohmann::json& j, const nlohmann::json& incoming_j)
{
  for (size_t i = 0; i < j.size(); i++) {
    nlohmann::json& patches = j[i][1];
    const nlohmann::json& incoming_patches = incoming_j[i][1];
    for (const auto& patch : incoming_patches) {
      patches.push_back(patch);
    }
  }
}
/**
 * @brief merge timing arrays
 *
 * @param j the timing array to merge to
 * @param incoming_j the timing array to merge
 */
static void
MergeIncomingTimings(nlohmann::json& j, const nlohmann::json& incoming_j)
{
  for (const auto& timing : incoming_j) {
    j.push_back(timing);
  }
}
/**
 * @brief Merge json objects from other ranks
 *
 * @param j the json object for this rank
 * @param incoming_j the json object from the other rank
 */
static void
MergeIncomingJson(nlohmann::json& j, nlohmann::json& incoming_j)
{
  if (j.contains("domains")) {
    MergeIncomingDomains(j["domains"], incoming_j["domains"]);
  }
  if (j.contains("timings") || incoming_j.contains("timings")) {
    MergeIncomingTimings(j["timings"], incoming_j["timings"]);
  }
}
void
to_json(nlohmann::json& output_j, const Timer& timer)
{
  int rank = timer.comm.getRank();
  int size = timer.comm.getSize();
  nlohmann::json j = *timer.root;

  if (timer.domains.size() > 0) {
    j["domains"] = timer.domains;
  }
  DecorateWithRank(j, rank);
  if (rank == 0) {
    for (int incoming_rank = 1; incoming_rank < size; incoming_rank++) {
      MPI_Status status;
      MPI_Probe(incoming_rank, 0, timer.comm.getMPIComm(), &status);

      int buffer_size;
      MPI_Get_count(&status, MPI_CHAR, &buffer_size);

      char incoming_j_string[buffer_size];
      MPI_Recv(incoming_j_string,
               buffer_size,
               MPI_CHAR,
               incoming_rank,
               0,
               timer.comm.getMPIComm(),
               &status);

      nlohmann::json incoming_j = nlohmann::json::parse(incoming_j_string);
      MergeIncomingJson(j, incoming_j);
    }
    if (j != nullptr) {
      j["comm_size"] = size;
    }
    for (const auto& el : j.items()) {
      output_j[el.key()] = el.value();
    }
  } else {
    std::string j_string = j.dump();
    MPI_Send(j_string.data(), (int)j_string.size() + 1, MPI_CHAR, 0, 0, timer.comm.getMPIComm());
  }
}
void
Timer::saveToFile(const std::string& filename) const
{
  nlohmann::json j = *this;
  if (comm.getRank() == 0) {
    std::ofstream out(filename, std::ofstream::out | std::ofstream::trunc);
    if (out.fail()) {
      throw RuntimeError("Failed to open file " + filename);
    }
    out << std::setw(4) << j;
    out.close();
  }
}
} // namespace ThunderEgg