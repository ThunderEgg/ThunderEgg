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
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/Timer.h>
#include <ThunderEgg/tpl/json.hpp>

#include <fstream>
#include <sstream>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

static Domain<2>
GetDomain(const Communicator& comm)
{
  vector<PatchInfo<2>> pinfos(1);

  int n = 10;
  double spacing = 0.01;
  int num_ghost = 1;

  pinfos[0].id = 0;
  pinfos[0].ns.fill(n);
  pinfos[0].spacings.fill(spacing);
  pinfos[0].num_ghost_cells = num_ghost;

  pinfos[0].rank = comm.getRank();

  Domain<2> d(comm, 1, { n, n }, num_ghost, pinfos.begin(), pinfos.end());
  return d;
}
static int
occurrences(const std::string& s, const std::string& target)
{
  int occurrences = 0;
  std::string::size_type pos = 0;
  while ((pos = s.find(target, pos)) != std::string::npos) {
    ++occurrences;
    pos += target.length();
  }
  return occurrences;
}
TEST_CASE("Timer to_json empty timer")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  nlohmann::json j = timer;
  REQUIRE_EQ(j, nullptr);
}
TEST_CASE("Timer to_json unassociated timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);

  timer.start("A");
  timer.stop("A");

  const nlohmann::json j = timer;

  if (comm.getRank() == 0) {
    REQUIRE_NE(j, nullptr);
    CHECK_EQ(j.size(), 2);
    CHECK_EQ(j["comm_size"], 2);

    CHECK_UNARY(j["timings"].is_array());
    CHECK_EQ(j["timings"].size(), 2);
    CHECK_EQ(j["timings"][0]["rank"], 0);
    CHECK_UNARY(j["timings"][0]["min"].is_number());
    CHECK_UNARY(j["timings"][0]["max"].is_number());
    CHECK_UNARY(j["timings"][0]["sum"].is_number());
    CHECK_UNARY(j["timings"][0]["num_calls"].is_number());
    CHECK_EQ(j["timings"][0]["name"], "A");
    CHECK_EQ(j["timings"][1]["rank"], 1);
    CHECK_UNARY(j["timings"][1]["min"].is_number());
    CHECK_UNARY(j["timings"][1]["max"].is_number());
    CHECK_UNARY(j["timings"][1]["sum"].is_number());
    CHECK_UNARY(j["timings"][1]["num_calls"].is_number());
    CHECK_EQ(j["timings"][1]["name"], "A");
  } else {
    CHECK_EQ(j, nullptr);
  }
}
TEST_CASE("Timer to_json unassociated timing with int info")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);

  timer.start("A");
  if (comm.getRank() == 0) {
    timer.addIntInfo("Example", 0);
  } else {
    timer.addIntInfo("Example", 1);
  }
  timer.stop("A");

  const nlohmann::json j = timer;

  if (comm.getRank() == 0) {
    REQUIRE_NE(j, nullptr);
    CHECK_EQ(j.size(), 2);
    CHECK_EQ(j["comm_size"], 2);

    CHECK_UNARY(j["timings"].is_array());
    CHECK_EQ(j["timings"].size(), 2);
    CHECK_EQ(j["timings"][0]["rank"], 0);
    CHECK_UNARY(j["timings"][0]["min"].is_number());
    CHECK_UNARY(j["timings"][0]["max"].is_number());
    CHECK_UNARY(j["timings"][0]["sum"].is_number());
    CHECK_UNARY(j["timings"][0]["num_calls"].is_number());
    CHECK_EQ(j["timings"][0]["name"], "A");

    CHECK_UNARY(j["timings"][0]["infos"].is_array());
    CHECK_EQ(j["timings"][0]["infos"].size(), 1);
    CHECK_EQ(j["timings"][0]["infos"][0]["name"], "Example");
    CHECK_EQ(j["timings"][0]["infos"][0]["min"], 0);
    CHECK_EQ(j["timings"][0]["infos"][0]["max"], 0);
    CHECK_EQ(j["timings"][0]["infos"][0]["num_calls"], 1);
    CHECK_EQ(j["timings"][0]["infos"][0]["sum"], 0);

    CHECK_EQ(j["timings"][1]["rank"], 1);
    CHECK_UNARY(j["timings"][1]["min"].is_number());
    CHECK_UNARY(j["timings"][1]["max"].is_number());
    CHECK_UNARY(j["timings"][1]["sum"].is_number());
    CHECK_UNARY(j["timings"][1]["num_calls"].is_number());
    CHECK_EQ(j["timings"][1]["name"], "A");

    CHECK_UNARY(j["timings"][1]["infos"].is_array());
    CHECK_EQ(j["timings"][1]["infos"].size(), 1);
    CHECK_EQ(j["timings"][1]["infos"][0]["name"], "Example");
    CHECK_EQ(j["timings"][1]["infos"][0]["min"], 1);
    CHECK_EQ(j["timings"][1]["infos"][0]["max"], 1);
    CHECK_EQ(j["timings"][1]["infos"][0]["num_calls"], 1);
    CHECK_EQ(j["timings"][1]["infos"][0]["sum"], 1);
  } else {
    CHECK_EQ(j, nullptr);
  }
}
TEST_CASE("Timer to_json unassociated timing with double info")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);

  timer.start("A");
  if (comm.getRank() == 0) {
    timer.addDoubleInfo("Example", 0);
  } else {
    timer.addDoubleInfo("Example", 1);
  }
  timer.stop("A");

  const nlohmann::json j = timer;

  if (comm.getRank() == 0) {
    REQUIRE_NE(j, nullptr);
    CHECK_EQ(j.size(), 2);
    CHECK_EQ(j["comm_size"], 2);

    CHECK_UNARY(j["timings"].is_array());
    CHECK_EQ(j["timings"].size(), 2);
    CHECK_EQ(j["timings"][0]["rank"], 0);
    CHECK_UNARY(j["timings"][0]["min"].is_number());
    CHECK_UNARY(j["timings"][0]["max"].is_number());
    CHECK_UNARY(j["timings"][0]["sum"].is_number());
    CHECK_UNARY(j["timings"][0]["num_calls"].is_number());
    CHECK_EQ(j["timings"][0]["name"], "A");

    CHECK_UNARY(j["timings"][0]["infos"].is_array());
    CHECK_EQ(j["timings"][0]["infos"].size(), 1);
    CHECK_EQ(j["timings"][0]["infos"][0]["name"], "Example");
    CHECK_EQ(j["timings"][0]["infos"][0]["min"], 0);
    CHECK_EQ(j["timings"][0]["infos"][0]["max"], 0);
    CHECK_EQ(j["timings"][0]["infos"][0]["num_calls"], 1);
    CHECK_EQ(j["timings"][0]["infos"][0]["sum"], 0);

    CHECK_EQ(j["timings"][1]["rank"], 1);
    CHECK_UNARY(j["timings"][1]["min"].is_number());
    CHECK_UNARY(j["timings"][1]["max"].is_number());
    CHECK_UNARY(j["timings"][1]["sum"].is_number());
    CHECK_UNARY(j["timings"][1]["num_calls"].is_number());
    CHECK_EQ(j["timings"][1]["name"], "A");

    CHECK_UNARY(j["timings"][1]["infos"].is_array());
    CHECK_EQ(j["timings"][1]["infos"].size(), 1);
    CHECK_EQ(j["timings"][1]["infos"][0]["name"], "Example");
    CHECK_EQ(j["timings"][1]["infos"][0]["min"], 1);
    CHECK_EQ(j["timings"][1]["infos"][0]["max"], 1);
    CHECK_EQ(j["timings"][1]["infos"][0]["num_calls"], 1);
    CHECK_EQ(j["timings"][1]["infos"][0]["sum"], 1);
  } else {
    CHECK_EQ(j, nullptr);
  }
}
TEST_CASE("Timer to_json two unassociated timings sequential")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.start("A");
  timer.stop("A");
  timer.start("B");
  timer.stop("B");
  const nlohmann::json j = timer;

  if (comm.getRank() == 0) {
    REQUIRE_NE(j, nullptr);
    CHECK_EQ(j.size(), 2);
    CHECK_EQ(j["comm_size"], 2);

    CHECK_UNARY(j["timings"].is_array());
    CHECK_EQ(j["timings"].size(), 4);

    CHECK_EQ(j["timings"][0]["rank"], 0);
    CHECK_UNARY(j["timings"][0]["min"].is_number());
    CHECK_UNARY(j["timings"][0]["max"].is_number());
    CHECK_UNARY(j["timings"][0]["sum"].is_number());
    CHECK_UNARY(j["timings"][0]["num_calls"].is_number());
    CHECK_EQ(j["timings"][0]["name"], "A");

    CHECK_EQ(j["timings"][1]["rank"], 0);
    CHECK_UNARY(j["timings"][1]["min"].is_number());
    CHECK_UNARY(j["timings"][1]["max"].is_number());
    CHECK_UNARY(j["timings"][1]["sum"].is_number());
    CHECK_UNARY(j["timings"][1]["num_calls"].is_number());
    CHECK_EQ(j["timings"][1]["name"], "B");

    CHECK_EQ(j["timings"][2]["rank"], 1);
    CHECK_UNARY(j["timings"][2]["min"].is_number());
    CHECK_UNARY(j["timings"][2]["max"].is_number());
    CHECK_UNARY(j["timings"][2]["sum"].is_number());
    CHECK_UNARY(j["timings"][2]["num_calls"].is_number());
    CHECK_EQ(j["timings"][2]["name"], "A");

    CHECK_EQ(j["timings"][3]["rank"], 1);
    CHECK_UNARY(j["timings"][3]["min"].is_number());
    CHECK_UNARY(j["timings"][3]["max"].is_number());
    CHECK_UNARY(j["timings"][3]["sum"].is_number());
    CHECK_UNARY(j["timings"][3]["num_calls"].is_number());
    CHECK_EQ(j["timings"][3]["name"], "B");
  } else {
    CHECK_EQ(j, nullptr);
  }
}
TEST_CASE("Timer to_json nested timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.start("A");
  timer.start("B");
  timer.stop("B");
  timer.stop("A");
  const nlohmann::json j = timer;
  if (comm.getRank() == 0) {
    REQUIRE_NE(j, nullptr);
    CHECK_EQ(j.size(), 2);
    CHECK_EQ(j["comm_size"], 2);

    CHECK_UNARY(j["timings"].is_array());
    CHECK_EQ(j["timings"].size(), 2);

    CHECK_EQ(j["timings"][0]["rank"], 0);
    CHECK_UNARY(j["timings"][0]["min"].is_number());
    CHECK_UNARY(j["timings"][0]["max"].is_number());
    CHECK_UNARY(j["timings"][0]["sum"].is_number());
    CHECK_UNARY(j["timings"][0]["num_calls"].is_number());
    CHECK_EQ(j["timings"][0]["name"], "A");
    CHECK_UNARY(j["timings"][0]["timings"].is_array());
    CHECK_EQ(j["timings"][0]["timings"].size(), 1);

    CHECK_EQ(j["timings"][1]["rank"], 1);
    CHECK_UNARY(j["timings"][1]["min"].is_number());
    CHECK_UNARY(j["timings"][1]["max"].is_number());
    CHECK_UNARY(j["timings"][1]["sum"].is_number());
    CHECK_UNARY(j["timings"][1]["num_calls"].is_number());
    CHECK_EQ(j["timings"][1]["name"], "A");
    CHECK_UNARY(j["timings"][1]["timings"].is_array());
    CHECK_EQ(j["timings"][1]["timings"].size(), 1);
  } else {
    CHECK_EQ(j, nullptr);
  }
}
TEST_CASE("Timer to_json domain timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.stopDomainTiming(0, "A");
  const nlohmann::json j = timer;
  if (comm.getRank() == 0) {
    REQUIRE_NE(j, nullptr);
    CHECK_EQ(j.size(), 3);
    CHECK_EQ(j["comm_size"], 2);

    CHECK_UNARY(j["domains"].is_array());
    CHECK_EQ(j["domains"].size(), 1);
    CHECK_UNARY(j["domains"][0].is_array());
    CHECK_EQ(j["domains"][0].size(), 2);
    CHECK_EQ(j["domains"][0][1].size(), 2);
    CHECK_UNARY(j["domains"][0][1].is_array());
    CHECK_EQ(j["domains"][0][1].size(), 2);
    CHECK_EQ(j["domains"][0][1][0]["rank"], 0);
    CHECK_EQ(j["domains"][0][1][1]["rank"], 1);

    CHECK_UNARY(j["timings"].is_array());
    CHECK_EQ(j["timings"].size(), 2);

    CHECK_EQ(j["timings"][0]["rank"], 0);
    CHECK_UNARY(j["timings"][0]["min"].is_number());
    CHECK_UNARY(j["timings"][0]["max"].is_number());
    CHECK_UNARY(j["timings"][0]["sum"].is_number());
    CHECK_UNARY(j["timings"][0]["num_calls"].is_number());
    CHECK_EQ(j["timings"][0]["name"], "A");
    CHECK_EQ(j["timings"][0]["domain_id"], 0);

    CHECK_EQ(j["timings"][1]["rank"], 1);
    CHECK_UNARY(j["timings"][1]["min"].is_number());
    CHECK_UNARY(j["timings"][1]["max"].is_number());
    CHECK_UNARY(j["timings"][1]["sum"].is_number());
    CHECK_UNARY(j["timings"][1]["num_calls"].is_number());
    CHECK_EQ(j["timings"][1]["name"], "A");
    CHECK_EQ(j["timings"][1]["domain_id"], 0);
  } else {
    REQUIRE_EQ(j, nullptr);
  }
}
TEST_CASE("Timer to_json domain timing only on rank 0")
{
  Communicator comm(MPI_COMM_WORLD);

  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  if (comm.getRank() == 0) {
    timer.startDomainTiming(0, "A");
    timer.stopDomainTiming(0, "A");
  }
  const nlohmann::json j = timer;
  if (comm.getRank() == 0) {
    REQUIRE_NE(j, nullptr);
    CHECK_EQ(j.size(), 3);
    CHECK_EQ(j["comm_size"], 2);

    CHECK_UNARY(j["domains"].is_array());
    CHECK_EQ(j["domains"].size(), 1);
    CHECK_UNARY(j["domains"][0].is_array());
    CHECK_EQ(j["domains"][0].size(), 2);
    CHECK_EQ(j["domains"][0][1].size(), 2);
    CHECK_UNARY(j["domains"][0][1].is_array());
    CHECK_EQ(j["domains"][0][1].size(), 2);
    CHECK_EQ(j["domains"][0][1][0]["rank"], 0);
    CHECK_EQ(j["domains"][0][1][1]["rank"], 1);

    CHECK_UNARY(j["timings"].is_array());
    CHECK_EQ(j["timings"].size(), 1);

    CHECK_EQ(j["timings"][0]["rank"], 0);
    CHECK_UNARY(j["timings"][0]["min"].is_number());
    CHECK_UNARY(j["timings"][0]["max"].is_number());
    CHECK_UNARY(j["timings"][0]["sum"].is_number());
    CHECK_UNARY(j["timings"][0]["num_calls"].is_number());
    CHECK_EQ(j["timings"][0]["name"], "A");
    CHECK_EQ(j["timings"][0]["domain_id"], 0);
  } else {
    REQUIRE_EQ(j, nullptr);
  }
}
TEST_CASE("Timer to_json domain timing only on rank 1")
{
  Communicator comm(MPI_COMM_WORLD);

  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  if (comm.getRank() == 1) {
    timer.startDomainTiming(0, "A");
    timer.stopDomainTiming(0, "A");
  }
  const nlohmann::json j = timer;
  if (comm.getRank() == 0) {
    REQUIRE_NE(j, nullptr);
    CHECK_EQ(j.size(), 3);
    CHECK_EQ(j["comm_size"], 2);

    CHECK_UNARY(j["domains"].is_array());
    CHECK_EQ(j["domains"].size(), 1);
    CHECK_UNARY(j["domains"][0].is_array());
    CHECK_EQ(j["domains"][0].size(), 2);
    CHECK_EQ(j["domains"][0][1].size(), 2);
    CHECK_UNARY(j["domains"][0][1].is_array());
    CHECK_EQ(j["domains"][0][1].size(), 2);
    CHECK_EQ(j["domains"][0][1][0]["rank"], 0);
    CHECK_EQ(j["domains"][0][1][1]["rank"], 1);

    CHECK_UNARY(j["timings"].is_array());
    CHECK_EQ(j["timings"].size(), 1);

    CHECK_EQ(j["timings"][0]["rank"], 1);
    CHECK_UNARY(j["timings"][0]["min"].is_number());
    CHECK_UNARY(j["timings"][0]["max"].is_number());
    CHECK_UNARY(j["timings"][0]["sum"].is_number());
    CHECK_UNARY(j["timings"][0]["num_calls"].is_number());
    CHECK_EQ(j["timings"][0]["name"], "A");
    CHECK_EQ(j["timings"][0]["domain_id"], 0);
  } else {
    REQUIRE_EQ(j, nullptr);
  }
}
TEST_CASE("Timer ostream empty timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "No timings to report"), 1);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer ostream unassociated timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.start("A");
  timer.stop("A");
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "A"), 1);
    CHECK_EQ(occurrences(s, "time (sec)"), 0);
    CHECK_EQ(occurrences(s, "average (sec)"), 1);
    CHECK_EQ(occurrences(s, "min (sec)"), 1);
    CHECK_EQ(occurrences(s, "max (sec)"), 1);
    CHECK_EQ(occurrences(s, "average calls per rank"), 1);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer ostream nested unassociated timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.start("A");
  timer.start("B");
  timer.stop("B");
  timer.stop("A");
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "A"), 2);
    CHECK_EQ(occurrences(s, "B"), 1);
    CHECK_EQ(occurrences(s, "time (sec)"), 0);
    CHECK_EQ(occurrences(s, "average (sec)"), 2);
    CHECK_EQ(occurrences(s, "min (sec)"), 2);
    CHECK_EQ(occurrences(s, "max (sec)"), 2);
    CHECK_EQ(occurrences(s, "average calls per rank"), 2);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer ostream sequential unassociated timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.start("A");
  timer.stop("A");
  timer.start("A");
  timer.stop("A");
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "A"), 1);
    CHECK_EQ(occurrences(s, "time (sec)"), 0);
    CHECK_EQ(occurrences(s, "average (sec)"), 1);
    CHECK_EQ(occurrences(s, "min (sec)"), 1);
    CHECK_EQ(occurrences(s, "max (sec)"), 1);
    CHECK_EQ(occurrences(s, "average calls per rank"), 1);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer ostream domain timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.stopDomainTiming(0, "A");
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "A"), 1);
    CHECK_EQ(occurrences(s, "time (sec)"), 0);
    CHECK_EQ(occurrences(s, "average (sec)"), 1);
    CHECK_EQ(occurrences(s, "min (sec)"), 1);
    CHECK_EQ(occurrences(s, "max (sec)"), 1);
    CHECK_EQ(occurrences(s, "average calls per rank"), 1);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer ostream domain timing two different domains sequential")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.addDomain(1, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.stopDomainTiming(0, "A");
  timer.startDomainTiming(1, "A");
  timer.stopDomainTiming(1, "A");
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "A"), 1);
    CHECK_EQ(occurrences(s, "time (sec)"), 0);
    CHECK_EQ(occurrences(s, "average (sec)"), 1);
    CHECK_EQ(occurrences(s, "min (sec)"), 1);
    CHECK_EQ(occurrences(s, "max (sec)"), 1);
    CHECK_EQ(occurrences(s, "average calls per rank"), 1);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer ostream domain timing two different domains sequential rank 0 has unrelated timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.addDomain(1, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.stopDomainTiming(0, "A");
  if (comm.getRank() == 0) {
    timer.startDomainTiming(0, "B");
    timer.stopDomainTiming(0, "B");
  }
  timer.startDomainTiming(1, "A");
  timer.stopDomainTiming(1, "A");
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "A"), 1);
    CHECK_EQ(occurrences(s, "B"), 1);
    CHECK_EQ(occurrences(s, "time (sec)"), 1);
    CHECK_EQ(occurrences(s, "average (sec)"), 1);
    CHECK_EQ(occurrences(s, "min (sec)"), 1);
    CHECK_EQ(occurrences(s, "max (sec)"), 1);
    CHECK_EQ(occurrences(s, "average calls per rank"), 1);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer ostream domain timing two different domains sequential nested")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.addDomain(1, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.startDomainTiming(0, "B");
  timer.stopDomainTiming(0, "B");
  timer.stopDomainTiming(0, "A");
  timer.startDomainTiming(1, "A");
  timer.startDomainTiming(1, "B");
  timer.stopDomainTiming(1, "B");
  timer.stopDomainTiming(1, "A");
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "A"), 2);
    CHECK_EQ(occurrences(s, "B"), 1);
    CHECK_EQ(occurrences(s, "time (sec)"), 0);
    CHECK_EQ(occurrences(s, "average (sec)"), 2);
    CHECK_EQ(occurrences(s, "min (sec)"), 2);
    CHECK_EQ(occurrences(s, "max (sec)"), 2);
    CHECK_EQ(occurrences(s, "average calls per rank"), 2);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer ostream domain timing two different domains sequential nested first domain has "
          "extra timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.addDomain(1, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.startDomainTiming(0, "C");
  timer.stopDomainTiming(0, "C");
  timer.startDomainTiming(0, "B");
  timer.stopDomainTiming(0, "B");
  timer.stopDomainTiming(0, "A");
  timer.startDomainTiming(1, "A");
  timer.startDomainTiming(1, "B");
  timer.stopDomainTiming(1, "B");
  timer.stopDomainTiming(1, "A");
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "A"), 3);
    CHECK_EQ(occurrences(s, "B"), 1);
    CHECK_EQ(occurrences(s, "C"), 1);
    CHECK_EQ(occurrences(s, "time (sec)"), 0);
    CHECK_EQ(occurrences(s, "average (sec)"), 3);
    CHECK_EQ(occurrences(s, "min (sec)"), 3);
    CHECK_EQ(occurrences(s, "max (sec)"), 3);
    CHECK_EQ(occurrences(s, "average calls per rank"), 3);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer ostream domain timing two different domains sequential nested second domain has "
          "extra timing")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.addDomain(1, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.startDomainTiming(0, "B");
  timer.stopDomainTiming(0, "B");
  timer.stopDomainTiming(0, "A");
  timer.startDomainTiming(1, "A");
  timer.startDomainTiming(0, "C");
  timer.stopDomainTiming(0, "C");
  timer.startDomainTiming(1, "B");
  timer.stopDomainTiming(1, "B");
  timer.stopDomainTiming(1, "A");
  stringstream ss;
  ss << timer;
  std::string s = ss.str();
  if (comm.getRank() == 0) {
    CHECK_EQ(occurrences(s, "A"), 3);
    CHECK_EQ(occurrences(s, "B"), 1);
    CHECK_EQ(occurrences(s, "C"), 1);
    CHECK_EQ(occurrences(s, "time (sec)"), 0);
    CHECK_EQ(occurrences(s, "average (sec)"), 3);
    CHECK_EQ(occurrences(s, "min (sec)"), 3);
    CHECK_EQ(occurrences(s, "max (sec)"), 3);
    CHECK_EQ(occurrences(s, "average calls per rank"), 3);
  } else {
    CHECK_EQ(s.size(), 0);
  }
}
TEST_CASE("Timer saveToFile new empty file")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.addDomain(1, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.startDomainTiming(0, "B");
  timer.stopDomainTiming(0, "B");
  timer.stopDomainTiming(0, "A");
  timer.startDomainTiming(1, "A");
  timer.startDomainTiming(0, "C");
  timer.stopDomainTiming(0, "C");
  timer.startDomainTiming(1, "B");
  timer.stopDomainTiming(1, "B");
  timer.stopDomainTiming(1, "A");

  if (comm.getRank() == 0) {
    std::remove("timerne2.json");
  }
  timer.saveToFile("timerne2.json");
  nlohmann::json j = timer;

  if (comm.getRank() == 0) {
    ifstream input("timerne2.json");
    nlohmann::json file_j;
    input >> file_j;
    nlohmann::json extra_j;
    CHECK_THROWS(input >> extra_j);
    CHECK_EQ(file_j.dump(), j.dump());
    input.close();
    std::remove("timerne2.json");
  }
}
TEST_CASE("Timer saveToFile overwrites file")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.addDomain(1, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.startDomainTiming(0, "B");
  timer.stopDomainTiming(0, "B");
  timer.stopDomainTiming(0, "A");
  timer.startDomainTiming(1, "A");
  timer.startDomainTiming(0, "C");
  timer.stopDomainTiming(0, "C");
  timer.startDomainTiming(1, "B");
  timer.stopDomainTiming(1, "B");
  timer.stopDomainTiming(1, "A");

  if (comm.getRank() == 0) {
    std::remove("timerow2.json");
  }
  timer.saveToFile("timerow2.json");
  timer.saveToFile("timerow2.json");
  nlohmann::json j = timer;

  if (comm.getRank() == 0) {
    ifstream input("timerow2.json");
    nlohmann::json file_j;
    input >> file_j;
    nlohmann::json extra_j;
    CHECK_THROWS(input >> extra_j);
    CHECK_EQ(file_j.dump(), j.dump());
    input.close();
    std::remove("timerow2.json");
  }
}
TEST_CASE("Timer saveToFile throws with nonexistant directory")
{
  Communicator comm(MPI_COMM_WORLD);
  Timer timer(comm);
  timer.addDomain(0, GetDomain(comm));
  timer.addDomain(1, GetDomain(comm));
  timer.startDomainTiming(0, "A");
  timer.startDomainTiming(0, "B");
  timer.stopDomainTiming(0, "B");
  timer.stopDomainTiming(0, "A");
  timer.startDomainTiming(1, "A");
  timer.startDomainTiming(0, "C");
  timer.stopDomainTiming(0, "C");
  timer.startDomainTiming(1, "B");
  timer.stopDomainTiming(1, "B");
  timer.stopDomainTiming(1, "A");

  if (comm.getRank() == 0) {
    CHECK_THROWS_AS(timer.saveToFile("surely/this/directory/does/not/exist/timer.json"), RuntimeError);
  } else {
    timer.saveToFile("surely/this/directory/does/not/exist/timer.json");
  }
}
