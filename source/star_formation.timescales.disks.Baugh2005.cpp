// Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
//
//
//    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
//
//    The California Institute of Technology shall allow RECIPIENT to use and
//    distribute this software subject to the terms of the included license
//    agreement with the understanding that:
//
//    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
//    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
//    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
//    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
//    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
//    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
//    USED.
//
//    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
//    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
//    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
//    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
//    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
//
//    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
//    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
//    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
//    USE OF THE SOFTWARE.
//
//    In addition, RECIPIENT also agrees that Caltech is under no obligation
//    to provide technical support for the Software.
//
//    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
//    of Derivative Works, public display or redistribution of the Software
//    other than those specified in the included license and the requirement
//    that all copies of the Software released be marked with the language
//    provided in this notice.
//
//    This software is separately available under negotiable license terms
//    from:
//    California Institute of Technology
//    Office of Technology Transfer
//    1200 E. California Blvd.
//    Pasadena, California 91125
//    http://www.ott.caltech.edu

//% Contains an implementation of the \cite{baugh_can_2005} star formation timescale rule for galactic disks.

namespace Star_Formation_Timescale_Disks_Baugh2005
{
  
//% Implements the \cite{baugh_can_2005} star formation timescale rule for galactic disks.
#include <string.h>
#include <math.h>
#include <stdio.h>
  
//: ./work/build/objects.tree_node.cWrappers.o
#include <objects.tree_node.cWrappers.h>
  
//: ./work/build/cosmology.functions.o
#include <cosmology.functions.h>
  
//: ./work/build/utility.input_parameters.o
#include <utility.input_parameters.h>
  
  // Define a type for our function.
  typedef double (func)(void *thisNode);
  
  // Declare our initialization function to be interoperable with Fortran.
  extern "C" 
  {
    void Star_Formation_Timescale_Disk_Baugh2005_Initialize(char *starFormationTimescaleDisksMethod,func **Star_Formation_Timescale_Disk_Get);
  }
  
  // Parameters of the timescale prescription.
  double starFormationDiskTimescale, starFormationDiskVelocityExponent, starFormationExpansionExponent;

  double Star_Formation_Timescale_Disk_Baugh2005(void *thisNode)
  {
    //% Return the star formation timescale for disks using the \cite{baugh_can_2005} prescription.
    double time, expansionFactor, timeScale, diskVelocity;
    double velocityNormalization  =200.0;

    diskVelocity=Tree_Node_Disk_Velocity(thisNode);
    if (diskVelocity <= 0.0) {
      timeScale=0.0;
    } else {
      time           =Tree_Node_Time(thisNode);
      expansionFactor=Expansion_Factor(time);
      timeScale      =starFormationDiskTimescale
	*pow(diskVelocity/velocityNormalization,starFormationDiskVelocityExponent)
	*pow(expansionFactor,starFormationExpansionExponent);
    }
    return timeScale;
  }
  
  //# <starFormationTimescaleDisksMethod>
  //#  <unitName>Star_Formation_Timescale_Disk_Baugh2005_Initialize</unitName>
  //# </starFormationTimescaleDisksMethod>
  void  Star_Formation_Timescale_Disk_Baugh2005_Initialize(char *starFormationTimescaleDisksMethod,func **Star_Formation_Timescale_Disk_Get)
  {
    //% Initialize the \cite{baugh_can_2005} star formation timescale in disks module.

    // Name of our method.
    char ourName           [] = "Baugh2005";

    // Names of our parameters.
    char timescaleParameter[] = "starFormationDiskTimescale";
    char velocityParameter [] = "starFormationDiskVelocityExponent";
    char expansionParameter[] = "starFormationExpansionExponent";

    // Default values for parameters.
    double starFormationDiskTimescaleDefault        =  8.0;
    double starFormationDiskVelocityExponentDefault = -3.0;
    double starFormationExpansionExponentDefault    =  0.0;

    // Test whether our method is selected.
    if (strcmp(starFormationTimescaleDisksMethod,ourName) == 0) {
      // It is so return a pointer to our timescale function.
      *Star_Formation_Timescale_Disk_Get=&Star_Formation_Timescale_Disk_Baugh2005;
      // Also read parameters.
      //@ <inputParameter>
      //@   <name>starFormationDiskTimescale</name>
      //@   <defaultValue>8.0</defaultValue>
      //@   <attachedTo>namespace</attachedTo>
      //@   <description>
      //@     The timescale (in Gyr) for star formation in the \cite{baugh_can_2005} prescription.
      //@   </description>
      //@ </inputParameter>
      Get_Input_Parameter(strlen(timescaleParameter),timescaleParameter,&starFormationDiskTimescale       ,&starFormationDiskTimescaleDefault       );
      //@ <inputParameter>
      //@   <name>starFormationDiskVelocityExponent</name>
      //@   <defaultValue>-3.0</defaultValue>
      //@   <attachedTo>namespace</attachedTo>
      //@   <description>
      //@     The exponent for disk velocity in the \cite{baugh_can_2005} prescription for star formation in galactic disks.
      //@   </description>
      //@ </inputParameter>
      Get_Input_Parameter(strlen(velocityParameter ),velocityParameter ,&starFormationDiskVelocityExponent,&starFormationDiskVelocityExponentDefault);
      //@ <inputParameter>
      //@   <name>starFormationExpansionExponent</name>
      //@   <defaultValue>0.0</defaultValue>
      //@   <attachedTo>namespace</attachedTo>
      //@   <description>
      //@     The exponent for expansion factor in the \cite{baugh_can_2005} prescription for star formation in galactic disks.
      //@   </description>
      //@ </inputParameter>
      Get_Input_Parameter(strlen(expansionParameter),expansionParameter,&starFormationExpansionExponent   ,&starFormationExpansionExponentDefault   );
    }

  }
  
}
