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

#include <ThunderEgg/Timer.h>
namespace ThunderEgg
{
Timer::Timing::Timing(const Timing *parent, int domain_id, const std::string &name)
: parent(parent), domain_id(domain_id), name(name)
{
}
Timer::Timing &Timer::Timing::getTiming(int domain_id, const std::string &name)
{
	auto key             = std::make_tuple(domain_id, name);
	auto timing_map_iter = timing_map.find(key);
	if (timing_map_iter == timing_map.end()) {
		timings.push_back(Timing(this, domain_id, name));
		timing_map.emplace(key, timings.back());
		return timings.back();
	} else {
		return timing_map_iter->second;
	}
}
void Timer::Timing::start()
{
	start_time = std::chrono::steady_clock::now();
}
void Timer::Timing::stop()
{
	std::chrono::duration<double> duration = std::chrono::steady_clock::now() - start_time;
	double                        time     = duration.count();
	sum += time;
	max = std::max(max, time);
	min = std::min(min, time);
	num_calls++;
}
void Timer::Timing::print(const std::string &parent_string, std::ostream &os) const
{
	std::string my_string = parent_string;
	if (domain_id != std::numeric_limits<int>::max() && parent != nullptr
	    && domain_id != parent->domain_id) {
		my_string += "(Domain " + std::to_string(domain_id) + ") ";
	}
	my_string += name;
	os << my_string << std::endl;
	os << std::string(my_string.size(), '-') << std::endl;

	if (num_calls == 1) {
		os << "   time (sec): " << sum << std::endl << std::endl;
	} else {
		os << "  total calls: " << num_calls << std::endl;
		os << "average (sec): " << sum / num_calls << std::endl;
	}
	for (const Timing &timing : timings) {
		timing.print(my_string + " -> ", os);
	}
}
void to_json(nlohmann::json &j, const Timer::Timing &timing)
{
	if (timing.name != "") {
		j["name"] = timing.name;
	}
	if (timing.domain_id != std::numeric_limits<int>::max()) {
		j["domain_id"] = timing.domain_id;
	}
	if (timing.sum != 0) {
		j["sum"] = timing.sum;
	}
	if (timing.num_calls != 0) {
		j["num_calls"] = timing.num_calls;
	}
	if (timing.max != std::numeric_limits<double>::min()) {
		j["max"] = timing.max;
	}
	if (timing.min != std::numeric_limits<double>::max()) {
		j["min"] = timing.min;
	}
	if (timing.timings.size() > 0) {
		j["timings"] = timing.timings;
	}
}
Timer::Timer()
{
	stack.push_back(root);
}
void Timer::start(const std::string &name)
{
	startDomainTiming(std::numeric_limits<int>::max(), name);
}
void Timer::stop(const std::string &name)
{
	stopDomainTiming(std::numeric_limits<int>::max(), name);
}
void Timer::addDomain(int domain_id, nlohmann::json domain)
{
	auto pair = domains.emplace(domain_id, domain);
	if (!pair.second) {
		throw RuntimeError("Domain with id " + std::to_string(domain_id)
		                   + " was already added to timer");
	}
}
void Timer::startDomainTiming(int domain_id, const std::string &name)
{
	if (domain_id != std::numeric_limits<int>::max() && domains.find(domain_id) == domains.end()) {
		throw RuntimeError("Domain with id " + std::to_string(domain_id)
		                   + " was not added to timer");
	}
	Timing &curr_timing = stack.back();
	Timing &next_timing = curr_timing.getTiming(domain_id, name);
	next_timing.start();
	stack.push_back(next_timing);
}
void Timer::stopDomainTiming(int domain_id, const std::string &name)
{
	Timing &curr_timing = stack.back();
	if (curr_timing.domain_id == domain_id && curr_timing.name == name && stack.size() > 1) {
		curr_timing.stop();
		stack.pop_back();
	} else {
		throw RuntimeError("Timer was expecting to end \"" + curr_timing.name + "\", instead got \""
		                   + name + "\"");
	}
}
std::ostream &operator<<(std::ostream &os, const Timer &timer)
{
	if (timer.stack.size() > 1) {
		const Timer::Timing &curr_timing = timer.stack.back();
		throw RuntimeError(
		"Cannot output Timer results with unfinished timings, check that all timings have been stopped. Currently waiting on timing \""
		+ curr_timing.name + "\"");
	}
	os << std::endl;
	os << "TIMING RESULTS" << std::endl;
	os << "==============" << std::endl << std::endl;

	for (const Timer::Timing &timing : timer.root.timings) {
		timing.print("", os);
	}
	return os;
}
void to_json(nlohmann::json &j, const Timer &timer)
{
	if (timer.root.timings.size() > 0) {
		j = timer.root;
		if (timer.domains.size() > 0) {
			j["domains"] = timer.domains;
		}
	}
}
} // namespace ThunderEgg