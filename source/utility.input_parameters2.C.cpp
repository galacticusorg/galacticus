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

//% Contains the implementation of a C++ class which provides input parameters.

#include <cstddef>
#include <string.h>
#include <utility.input_parameters2.h>

extern "C"
{
  void inputParametersFinalize       (void *selfParameters);
  void inputParametersSetGlobalC     (void *parameters);
  void inputParametersSubParameters  (void *selfParameters,int parameterNameLength, char *parameterName, void   **subParameters);
  void inputParametersValueNameDouble(void *selfParameters,int parameterNameLength, char *parameterName, double *parameterValue, double *defaultValue, int *errorStatus, bool *writeOutput);
  void inputParametersValueNameLong  (void *selfParameters,int parameterNameLength, char *parameterName, long   *parameterValue, long   *defaultValue, int *errorStatus, bool *writeOutput);
}

namespace globalParameters
{
  //% Globally-visible input parameters object.
  inputParameters parameters = inputParameters(NULL);
}

inputParameters::inputParameters(void *parameters)
{
  //% Constructor for C++ {\normalfont \ttfamily inputParameters} class.
  selfParameters = parameters;
}

inputParameters::~inputParameters()
{
  //% Destructor for C++ {\normalfont \ttfamily inputParameters} class. 
  if ( &selfParameters != &globalParameters::parameters.selfParameters )
    {
      // Don't destruct the global parameters - this will be taken care of from the Fortran side.
      inputParametersFinalize(selfParameters);
    }
}

inputParameters inputParameters::subParameters(char parameterName[])
{
  //% Retrieve subparameters from a C++ {\normalfont \ttfamily inputParameters} class.
  void *subParameters;
  inputParametersSubParameters(selfParameters,strlen(parameterName),parameterName,&subParameters);
  return inputParameters(subParameters);
}

void inputParameters::value(char parameterName[], double *parameterValue, double *defaultValue, int *errorStatus, bool *writeOutput)
{
  //% Retrieve a {\normalfont \ttfamily double} value from a C++ {\normalfont \ttfamily inputParameters} class.
  inputParametersValueNameDouble(selfParameters,strlen(parameterName),parameterName,parameterValue,defaultValue,errorStatus,writeOutput);
}

void inputParameters::value(char parameterName[], long   *parameterValue, long *defaultValue  , int *errorStatus, bool *writeOutput)
{
  //% Retrieve a {\normalfont \ttfamily long} value from a C++ {\normalfont \ttfamily inputParameters} class.
  inputParametersValueNameLong(selfParameters,strlen(parameterName),parameterName,parameterValue,defaultValue,errorStatus,writeOutput);
}

void inputParameters::setParameters(void *parameters)
{
  //% Set the parameters accessible via the C++ interface.
  selfParameters=parameters;
}

void inputParametersSetGlobalC(void *parameters)
{
  //% Set the global input parameters accessible via the C++ interface.
  globalParameters::parameters.setParameters(parameters);
}
