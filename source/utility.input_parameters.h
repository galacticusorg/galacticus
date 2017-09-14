// Copyright 2009, 2010, 2011 Andrew Benson <abenson@obs.carnegiescience.edu>
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

//% Contains definitions of a C++ class to provide input parameters.

class inputParameters
{
  //% A C++ class which provides input parameters.
  void *selfParameters;
 public:
  inputParameters              (void *parameters    );
  ~inputParameters             (                    );
  inputParameters subParameters(char parameterName[]);
  void            setParameters(void *parameters    );
  void            value        (char parameterName[], double *parameterValue, double *defaultValue = NULL, int *errorStatus = NULL, bool *writeOutput = NULL);
  void            value        (char parameterName[], long   *parameterValue, long   *defaultValue = NULL, int *errorStatus = NULL, bool *writeOutput = NULL);
};

namespace globalParameters
{
  extern inputParameters parameters;
}
