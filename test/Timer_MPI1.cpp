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
TEST_CASE("Two Timings Sequential Stop with empty string second before started", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");
	REQUIRE_THROWS_AS(timer.stop(""), RuntimeError);
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
}
TEST_CASE("Timer ostream operator throws with unfinished timing", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.start("B");
	stringstream writer;
	REQUIRE_THROWS_AS(writer << timer, RuntimeError);
}
TEST_CASE("Timer DomainTiming", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	timer.stopDomainTiming(0, "A");
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
TEST_CASE("Timer  two sequential DomainTimings", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	timer.stopDomainTiming(0, "A");
	timer.startDomainTiming(0, "B");
	timer.stopDomainTiming(0, "B");
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
	CHECK(getNextLine() == "(Domain 0) A");
	CHECK(getNextLine() == "------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "(Domain 0) B");
	CHECK(getNextLine() == "------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Timer DomainTiming with unassociated timing before")
{
	Timer timer;
	timer.start("A");
	timer.stop("A");
	timer.startDomainTiming(0, "A");
	timer.stopDomainTiming(0, "A");
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
	CHECK(getNextLine() == "(Domain 0) A");
	CHECK(getNextLine() == "------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Timer DomainTiming nested in unassociated timing", "[Timer]")
{
	Timer timer;
	timer.start("A");
	timer.startDomainTiming(0, "A");
	timer.stopDomainTiming(0, "A");
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
	CHECK(getNextLine() == "A -> (Domain 0) A");
	CHECK(getNextLine() == "------------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Timer nested DomainTiming same domain", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	timer.startDomainTiming(0, "B");
	timer.stopDomainTiming(0, "B");
	timer.stopDomainTiming(0, "A");
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
	CHECK(getNextLine() == "(Domain 0) A");
	CHECK(getNextLine() == "------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "(Domain 0) A -> B");
	CHECK(getNextLine() == "-----------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Timer nested DomainTiming different domain", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	timer.startDomainTiming(1, "B");
	timer.stopDomainTiming(1, "B");
	timer.stopDomainTiming(0, "A");
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
	CHECK(getNextLine() == "(Domain 0) A");
	CHECK(getNextLine() == "------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "(Domain 0) A -> (Domain 1) B");
	CHECK(getNextLine() == "----------------------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Timer consecutive DomainTiming two different ids", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	timer.stopDomainTiming(0, "A");
	timer.startDomainTiming(1, "A");
	timer.stopDomainTiming(1, "A");
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
	CHECK(getNextLine() == "(Domain 0) A");
	CHECK(getNextLine() == "------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(getNextLine() == "(Domain 1) A");
	CHECK(getNextLine() == "------------");
	CHECK(getNextLine().find("time (sec):") != string::npos);
	CHECK(getNextLine() == "");
	CHECK(reader.eof());
}
TEST_CASE("Timer Two DomainTimings Sequential Stop second before started", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	timer.stopDomainTiming(0, "A");
	int          id   = GENERATE(0, 1);
	const string name = GENERATE("A", "B", "");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(id, name), RuntimeError);
}
TEST_CASE("Timer DomainTimings Nested Wrong id on stop", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(1, "A"), RuntimeError);
}
TEST_CASE("Timer DomainTimings Nested Wrong name on stop", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(0, "blah"), RuntimeError);
}
TEST_CASE("Timer DomainTimings Nested Wrong name and id on stop", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(1, "blah"), RuntimeError);
}
TEST_CASE("Timer Two DomainTimings Nested Wrong Order", "[Timer]")
{
	Timer timer;
	timer.startDomainTiming(0, "A");
	timer.startDomainTiming(1, "B");
	REQUIRE_THROWS_AS(timer.stopDomainTiming(0, "A"), RuntimeError);
}