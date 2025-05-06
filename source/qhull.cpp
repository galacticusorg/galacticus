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

//% Implements Fortran-callable wrappers around the qhull library (http://www.qhull.org).

#ifdef QHULLAVAIL
#include <string>
#include "Qhull.h"

extern "C" 
{
  orgQhull::Qhull * convexHullConstructorC(int n, double *points[], int *status);
  void convexHullDestructorC(orgQhull::Qhull *qhull);
  double convexHullVolumeC(orgQhull::Qhull *qhull);
  bool convexHullPointIsInHullC(orgQhull::Qhull *qhull, double point[]);
}

orgQhull::Qhull * convexHullConstructorC(int n, double *points[], int *status)
{
  //% Constructor for convex hull objects.
   
    std::string comment       = ""; // rbox commands, see http://www.qhull.org/html/rbox.htm
    std::string qhull_command = ""; // For qhull commands, see http://www.qhull.org/html/qhull.htm

    try
      {
	int ndim = 3;
	orgQhull::Qhull *qhull = new orgQhull::Qhull(comment.c_str(), ndim, n, *points, qhull_command.c_str());
	*status = 0;
	return qhull;
      }
    catch (orgQhull::QhullError &e)
      {
        std::cerr << e.what() << std::endl;
	*status = e.errorCode();
        return NULL;
      }
}

void convexHullDestructorC(orgQhull::Qhull *qhull)
{
  //% Destroy a convex hull object.
  delete(qhull);
}

double convexHullVolumeC(orgQhull::Qhull *qhull)
{
  //% Return the volume of the convex hull.
  return qhull->volume();
}

bool convexHullPointIsInHullC(orgQhull::Qhull *qhull, double point[])
{
  //% Return whether or not a point is inside the hull.
   unsigned int isoutside_;
   double bestdist;
   facetT *facet = qh_findbestfacet(qhull->qh(),point,qh_ALL, &bestdist, &isoutside_);
   bool isInside = isoutside_ == 0;
   return isInside;
}
#endif
