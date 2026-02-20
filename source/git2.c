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
  /* Return the hash of the repo HEAD. */
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

int gitDescendantOf(char *repoPath, const char commitHash[41], const char ancestorHash[41]) {
  /* Determine if a hash is an ancestor of another hash. */

  /* Initialize return value to failure status. */
  int isDescendantOf = 2;
  /* Initialize libgit2 */
  git_libgit2_init();
  /* Open the repository. */
  git_repository *repo;
  int err = git_repository_open_ext(&repo, repoPath, 0, NULL);
  if (err == 0) {
    /* Get the OIDs for the two commits. */
    git_oid commit, ancestor;
    int errCommit   = git_oid_fromstr(&commit  ,commitHash  );
    int errAncestor = git_oid_fromstr(&ancestor,ancestorHash);
    /* Check if is a descendant. */
    if (errCommit == 0 && errAncestor == 0) {
      if (git_oid_equal(&commit,&ancestor)) {
	/* libgit2 doesn't consider a commit to be its own descendant  */
	isDescendantOf = 1;
      } else {
	isDescendantOf = git_graph_descendant_of(repo,&commit,&ancestor);
      }
    }
  }
  git_repository_free(repo);
  return isDescendantOf;
}

#endif
