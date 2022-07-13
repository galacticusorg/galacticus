// Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
//           2019
//    Andrew Benson <abenson@carnegiescience.edu>
//
// This file is part of Galacticus.
//
//    Galacticus is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Galacticus is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

//% Implements an interface for getting the current virtual size of the running process.
#ifdef PROCPS
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

// Declare our function to be interoperable with Fortran.
extern "C" 
{
  long Memory_Usage_Get_C();
}

long Memory_Usage_Get_C() {
  using std::ios_base;
  using std::ifstream;
  using std::string;

  // Open a stream to our own stat file in the procfs file system.
  ifstream stat_stream("/proc/self/stat",ios_base::in);

  // Variables to read all of the entries in "stat" - most of these we will ignore.
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;

  // Variables to hold the memory stat data.
  unsigned long vsize;
  long rss;

  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	      >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	      >> utime >> stime >> cutime >> cstime >> priority >> nice
	      >> O >> itrealvalue >> starttime >> vsize >> rss; // No need to read beyond this as we've now got what we want.

  stat_stream.close();

  long page_size_b = sysconf(_SC_PAGE_SIZE); // Get page size so we can convert result to bytes.
  rss *= page_size_b;
  return rss;
}
#endif
