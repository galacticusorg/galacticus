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

//% Implements Fortran-callable wrappers around the ANN (Approximate Nearest Neighbors) library.

#ifdef ANNAVAIL
#include <ANN/ANN.h>
#include <stdio.h>

// Declare our functions to be interoperable with Fortran.
extern "C" 
{
  ANNkd_tree * nearestNeighborsConstructorC(int n, int d, double *pa);
  void nearestNeighborsDestructorC(ANNkd_tree * ANN);
  void nearestNeighborsCloseC();
  void nearestNeighborsSearchC(ANNkd_tree * ANN, double *point, int neighborCount, double tolerance, int *neighborIndex, double *neighborDistance);
  int nearestNeighborsSearchFixedRadiusC(ANNkd_tree * ANN, double *point, double radiusSquared, int neighborCount, int *neighborIndex, double *neighborDistance, double tolerance);
}
  
ANNkd_tree * nearestNeighborsConstructorC(int n, int d, double *pa) {
  //% Fortran-callable wrapper around the ANN constructor.
  ANNpoint      ANNp;
  ANNpointArray ANNpa;
  int           i, j;

  // Allocate points array.
  ANNpa = annAllocPts(n,d);
  // Load points into ANN data structure.
  for(i=0;i<n;++i) {
    ANNp = ANNpa[i];
    for(j=0;j<d;++j) {
      ANNp[j] = pa[i+j*n];
    }
  }
  // Create the KD-tree.
  ANNkd_tree *ANN;
  ANN = new ANNkd_tree(ANNpa,n,d);
  return ANN;
}

void nearestNeighborsDestructorC(ANNkd_tree * ANN) {
  //% Fortran-callable wrapper around the ANN destructor.
  ANNpointArray ANNpa;

  if ( ANN != NULL ) {
    // Get a pointer to the array of points, then deallocate it.
    ANNpa = ANN->thePoints();
    annDeallocPts(ANNpa);
    // Explicitly destruct the ANN KD-tree object.
    ANN->~ANNkd_tree();
  }
  return;
}

void nearestNeighborsCloseC() {
  //% Fortran-callable wrapper around the ANN library close function.
  annClose();
  return;
}

void nearestNeighborsSearchC(ANNkd_tree * ANN, double *point, int neighborCount, double tolerance, int *neighborIndex, double *neighborDistance) {
  //% Fortran-callable wrapper around the ANN search function.
  ANNpoint     ANNp;
  ANNidxArray  ANNindices;
  ANNdistArray ANNdistances;
  int          j;

  // Create an ANN point.
  ANNp = annAllocPt(ANN->theDim());
  for(j=0;j<ANN->theDim();++j) {
    ANNp[j] = point[j];
  }
  // Allocate index and distance arrays.
  ANNindices   = new ANNidx [neighborCount];
  ANNdistances = new ANNdist[neighborCount];
  // Perform the search.
  ANN->annkSearch(ANNp,neighborCount,ANNindices,ANNdistances,tolerance);
  // Deallocate the point.
  annDeallocPt(ANNp);
  // Copy results to output arrays.
  for(j=0;j<neighborCount;++j) {
    neighborIndex   [j] = ANNindices  [j];
    neighborDistance[j] = ANNdistances[j];
  }
  // Clean up.
  delete [] ANNindices;
  delete [] ANNdistances;
  return;
}

int nearestNeighborsSearchFixedRadiusC(ANNkd_tree * ANN, double *point, double radiusSquared, int neighborCount, int *neighborIndex, double *neighborDistance, double tolerance) {
  //% Fortran-callable wrapper around the ANN search function.
  ANNpoint     ANNp;
  ANNidxArray  ANNindices;
  ANNdistArray ANNdistances;
  int          j;
  int          neighborsFound;

  // Create an ANN point.
  ANNp = annAllocPt(ANN->theDim());
  for(j=0;j<ANN->theDim();++j) {
    ANNp[j] = point[j];
  }
  // Allocate index and distance arrays.
  ANNindices   = new ANNidx [neighborCount];
  ANNdistances = new ANNdist[neighborCount];
  // Perform the search.
  neighborsFound = ANN->annkFRSearch(ANNp,radiusSquared,neighborCount,ANNindices,ANNdistances,tolerance);
  // Deallocate the point.
  annDeallocPt(ANNp);
  // Copy results to output arrays.
  for(j=0;j<neighborCount;++j) {
    neighborIndex   [j] = ANNindices  [j];
    neighborDistance[j] = ANNdistances[j];
  }
  // Clean up.
  delete [] ANNindices;
  delete [] ANNdistances;
  return neighborsFound;
}
#endif
