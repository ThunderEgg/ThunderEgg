#include "catch.hpp"
#include <Thunderegg/Timer.h>
#include <sstream>
using namespace std;
using namespace Thunderegg;
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
	REQUIRE_THROWS_AS(timer.stop("B"), TimerException);
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
	REQUIRE_THROWS_AS(timer.stop("A"), TimerException);
}