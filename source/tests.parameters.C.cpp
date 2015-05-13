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

//% Contains test code for testing the C interface to input parameter functions.

#include <iostream>
#include <utility.input_parameters2.h>

extern "C"
{
  void testParametersC(bool results[3]);
}

void testParametersC(bool results[3])
{
  //% Test the C interface to input parameter functions.

  // Get parameters from the root-level of global parameters.
  double myDbl ;
  long   myLong;
  long   myLongDefault        = 34567;
  char    dblParameterName [] = "dblFromC" ;
  char   longParameterName [] = "longFromC";
  //@ <inputParameter>
  //@   <name>dblFromC</name>
  //@   <description>Test parameter</description>
  //@   <attachedTo>file</attachedTo>
  //@   <type>real</type>
  //@   <cardinality>1</cardinality>
  //@ </inputParameter>
  globalParameters::parameters.value( dblParameterName,&myDbl );
  //@ <inputParameter>
  //@   <name>longFromC</name>
  //@   <description>Test parameter</description>
  //@   <attachedTo>file</attachedTo>
  //@   <type>integer</type>
  //@   <cardinality>1</cardinality>
  //@ </inputParameter>
  globalParameters::parameters.value(longParameterName,&myLong,&myLongDefault);
  results[0] = (myDbl  == 1.2434);
  results[1] = (myLong == 12434 );

  // Get a set of subparameters and extract a value.
  char            parentParameterName    [] = "cosmologyParametersMethod";
  char            cosmologyParameterName [] = "OmegaMatter"              ;
  inputParameters subParameters = globalParameters::parameters.subParameters(parentParameterName);
  subParameters.value(cosmologyParameterName,&myDbl);
  results[2] = (myDbl == 0.2725);
}
