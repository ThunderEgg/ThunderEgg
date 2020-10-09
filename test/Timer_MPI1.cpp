#include "catch.hpp"
#include <ThunderEgg/Timer.h>
#include <sstream>
using namespace std;
using namespace ThunderEgg;
TEST_CASE("No Timing", "[Timer]")
{
	Timer timer;

	// check output
	stringstream writer;
	writer << timer;
	stringstream reader(writer.str());

	auto getNextLine = [&]() {
		string line;
		getline(reader, line);
		return line;
	};

	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "TIMING RESULTS");
	CHECK(getNextLine() == "==============");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Single Timing", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");

	// check output
	stringstream writer;
	writer << timer;
	stringstream reader(writer.str());

	auto getNextLine = [&]() {
		string line;
		getline(reader, line);
		return line;
	};

	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "TIMING RESULTS");
	CHECK(getNextLine() == "==============");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Two Timings Sequential", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");
	timer.start("B");
	timer.stop("B");

	// check output
	stringstream writer;
	writer << timer;
	stringstream reader(writer.str());

	auto getNextLine = [&]() {
		string line;
		getline(reader, line);
		return line;
	};

	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "TIMING RESULTS");
	CHECK(getNextLine() == "==============");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "B");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Two Timings Sequential Stop second before started", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");
	REQUIRE_THROWS_AS(timer.stop("B"), RuntimeError);
}
TEST_CASE("Two Timings Nested", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.start("B");
	timer.stop("B");
	timer.stop("A");

	// check output
	stringstream writer;
	writer << timer;
	stringstream reader(writer.str());

	auto getNextLine = [&]() {
		string line;
		getline(reader, line);
		return line;
	};

	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "TIMING RESULTS");
	CHECK(getNextLine() == "==============");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A -> B");
	CHECK(getNextLine() == "------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Two Timings Nested Wrong Order", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.start("B");
	REQUIRE_THROWS_AS(timer.stop("A"), RuntimeError);
	timer.stop("B");
	timer.stop("A");
}
TEST_CASE("Timer deconstructor throws with unfinished timing", "[Timer]")
{
	Timer *timer = new Timer();
	timer->start("A");
	timer->start("B");
	REQUIRE_THROWS_AS(delete timer, RuntimeError);
}
TEST_CASE("Timer setDomainId", "[Timer]")
{
	Timer timer;
	timer.setDomainId(0);
	timer.start("A");
	timer.stop("A");
	// check output
	stringstream writer;
	writer << timer;
	stringstream reader(writer.str());

	auto getNextLine = [&]() {
		string line;
		getline(reader, line);
		return line;
	};

	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "TIMING RESULTS");
	CHECK(getNextLine() == "==============");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "Domain 0");
	CHECK(getNextLine() == "========");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Timer setDomainId with unassociated timing before", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");
	timer.setDomainId(0);
	timer.start("A");
	timer.stop("A");
	// check output
	stringstream writer;
	writer << timer;
	stringstream reader(writer.str());

	auto getNextLine = [&]() {
		string line;
		getline(reader, line);
		return line;
	};

	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "TIMING RESULTS");
	CHECK(getNextLine() == "==============");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "Domain 0");
	CHECK(getNextLine() == "========");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Timer setDomainId with unassociated timing nested", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.setDomainId(0);
	timer.start("A");
	timer.stop("A");
	timer.stop("A");
	// check output
	stringstream writer;
	writer << timer;
	stringstream reader(writer.str());

	auto getNextLine = [&]() {
		string line;
		getline(reader, line);
		return line;
	};

	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "TIMING RESULTS");
	CHECK(getNextLine() == "==============");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "Domain 0");
	CHECK(getNextLine() == "========");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A -> A");
	CHECK(getNextLine() == "------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Timer setDomainId called again with calling unsetDomainId throws", "[Timer]")
{
	Timer timer;
	timer.setDomainId(0);
	CHECK_THROWS_AS(timer.setDomainId(2), RuntimeError);
}
TEST_CASE("Timer setDomainId two different ids", "[Timer]")
{
	Timer timer;
	timer.setDomainId(0);
	timer.start("A");
	timer.stop("A");
	timer.unsetDomainId();
	timer.setDomainId(1);
	timer.start("A");
	timer.stop("A");
	// check output
	stringstream writer;
	writer << timer;
	stringstream reader(writer.str());

	auto getNextLine = [&]() {
		string line;
		getline(reader, line);
		return line;
	};

	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "TIMING RESULTS");
	CHECK(getNextLine() == "==============");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "Domain 1");
	CHECK(getNextLine() == "========");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "Domain 0");
	CHECK(getNextLine() == "========");
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "A");
	CHECK(getNextLine() == "-");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}