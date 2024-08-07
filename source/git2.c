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

//% Implements Fortran-callable wrappers around the libgit2 library.

#ifdef GIT2AVAIL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <git2.h>
 
void repoHeadHash(char *repoPath, char hash[41]) {
    int err;
    git_repository *repo;
    git_reference *head;

    /* Initialize hash to "unknown" */
    strcpy(hash,"unknown");
    
    /* Initialize libgit2 */
    git_libgit2_init();
    
    /* Open the repository. */
    err = git_repository_open_ext(&repo, repoPath, 0, NULL);
    if (err == 0) {
      /* Find the HEAD of the repo. */
      err = git_repository_head(&head, repo);
      if (err == 0) {
	/* Extract the hash from HEAD */
	int referenceType = git_reference_type(head);
	if (referenceType == GIT_REF_OID) {
	  /* Get the hash. */
	  const git_oid *oid = git_reference_target(head);
	  git_oid_tostr(hash, 41, oid);
	}
      }
      git_reference_free(head);
    }
    git_repository_free(repo);
    return;
}

#endif
