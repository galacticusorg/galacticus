// Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
//           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

#include <git2.h>

// Reference every libgit2 function that git2.c actually calls. Taking the
// address of an undeclared identifier is a hard compile error in C, so this
// probe only succeeds when the installed git2.h is new enough to declare all
// of them. That makes the Makefile select GIT2AVAIL only when git2.c will in
// fact compile, and otherwise fall back cleanly to GIT2UNAVAIL (using the
// `git` command line) instead of hard-failing later in git2.c -- which is what
// happened when an older libgit2 lacked git_graph_descendant_of.
// Keep this list in sync with the git_* calls in git2.c.
void *libgit2RequiredSymbols[] = {
  (void *) &git_repository_open_ext,
  (void *) &git_repository_head,
  (void *) &git_reference_type,
  (void *) &git_reference_target,
  (void *) &git_oid_tostr,
  (void *) &git_reference_free,
  (void *) &git_repository_free,
  (void *) &git_oid_fromstr,
  (void *) &git_oid_equal,
  (void *) &git_graph_descendant_of
};

int main() {
}
