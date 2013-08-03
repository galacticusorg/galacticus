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

//% Contains a set of templates and overloaded wrapper functions to provide a C interface the the input parameters module.

extern "C"
{
  void Get_Input_Parameter_Double (int parameterNameLength, char *parameterName, double *parameterValue, double *defaultValue);
  void Get_Input_Parameter_Integer(int parameterNameLength, char *parameterName, int    *parameterValue, int    *defaultValue);
}

void Get_Input_Parameter(int parameterNameLength, char *parameterName,double *parameterValue, double *defaultValue = NULL)
{
  //% Get the value of a {\tt double} parameter from C.
  Get_Input_Parameter_Double(parameterNameLength,parameterName,parameterValue,defaultValue);
}

void Get_Input_Parameter(int parameterNameLength, char *parameterName,int    *parameterValue, int    *defaultValue = NULL)
{
  //% Get the value of an {\tt int} parameter from C.
  Get_Input_Parameter_Integer(parameterNameLength,parameterName,parameterValue,defaultValue);
}
