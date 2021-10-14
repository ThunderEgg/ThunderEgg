#!/usr/bin/python3

liscense_header = """/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
"""
liscense_footer = """ *
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
"""

import os
import subprocess
from datetime import datetime
from datetime import timedelta

class authorship:
	def __init__(self, author, year):
		self.author = author
		self.min = year
		self.max = year

	def addyear(self, year):
		self.min = min(self.min,year)
		self.max = max(self.max,year)
	
	def range(self):
		if self.min == self.max:
			return "%i     "%(self.min)
		else:
			return "%i-%i"%(self.min,self.max)


	def __repr__(self):
		return "<author %s, min %i, max %i>"%(self.author,self.min,self.max)

	def __str__(self):
		return "Copyright (c) %s %s"%(self.range(),self.author)


def parse_author(line):
	return line[7:]

def parse_time(line):
	return int(line[12:])

def parse_time_zone(line):
	return int(line[10:])

def get_year(time, time_zone):
	dt = datetime.fromtimestamp(time) + timedelta(minutes=time_zone)
	return dt.year


def parse_authors(output):
	authorships = {}
	lines = output.splitlines()
	current_line = 0
	while current_line < len(lines):
		if lines[current_line].startswith("author"):
			author = parse_author(lines[current_line])
			current_line += 2
			time = parse_time(lines[current_line])
			current_line += 1
			time_zone = parse_time_zone(lines[current_line])
			current_line += 1
			year = get_year(time, time_zone)
			if "updated copyright headers" not in lines[current_line+4]:
				if authorships.get(author):
					authorships[author].addyear(year)
				else:
					authorships[author] = authorship(author,year)
		else:
			current_line += 1
		
	return authorships.values()

def remove_copyright_notice(lines):
	current_line = 0
	start = 0
	end = 0
	found = False
	while current_line < len(lines):
		if lines[current_line].startswith("/******"):
			tmp_start = current_line
			current_line += 1
			while "*/" not in lines[current_line]:
				if "Copyright" in lines[current_line]:
					found = True
				current_line += 1
			if found:
				start = tmp_start
				end = current_line
				break
			current_line += 1
		else:
			current_line += 1
	if found:
		del lines[start:end+1]

def add_copyright_notice(lines,filename):
	authorships = parse_authors(subprocess.run(["git","blame","--incremental",filename],capture_output=True, text=True).stdout)
	liscense = liscense_header
	for authorship in authorships:
		liscense += " *  " + str(authorship) + "\n"
	liscense += liscense_footer
	return liscense.splitlines(True) + lines


def shortlog(filename):
	lines = []
	with open(filename, "r+") as file:
		lines = file.readlines()
		remove_copyright_notice(lines)
		lines = add_copyright_notice(lines,filename)
		file.seek(0)
		file.writelines(lines)
		file.truncate()

for root in ['./test','./src','./apps']:
	for directory, subdirectories, files in os.walk(root):
		for file in files:
			if file.endswith(".cpp") or file.endswith(".h") or file.endswith(".h.in"):
				shortlog(directory+"/"+file)
