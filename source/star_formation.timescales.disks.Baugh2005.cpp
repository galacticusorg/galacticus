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

//% Contains an implementation of the \cite{baugh_can_2005} star formation timescale rule for galactic disks.

namespace Star_Formation_Timescale_Disks_Baugh2005
{
  
//% Implements the \cite{baugh_can_2005} star formation timescale rule for galactic disks.
#include <string.h>
#include <math.h>
#include <stdio.h>
  
//: ./work/build/objects.nodes.bindings.C.o
#include <objects.nodes.bindings.C.h>
  
//: ./work/build/cosmology.functions.o
#include <cosmologyFunctions.h>
  
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
    double                  time, expansionFactor, timeScale, diskVelocity;
    double                  velocityNormalization    =200.0    ;
    cosmologyFunctionsClass cosmologyFunctionsDefault          ;
    nodeComponentBasic      thisBasicComponent       (thisNode);
    nodeComponentDisk       thisDiskComponent        (thisNode);

    diskVelocity=thisDiskComponent.velocity();
    if (diskVelocity <= 0.0) {
      timeScale=0.0;
    } else {
      time           =thisBasicComponent.time();
      expansionFactor=cosmologyFunctionsDefault.expansionFactor(time);
      timeScale      =starFormationDiskTimescale
     	*pow(diskVelocity/velocityNormalization,starFormationDiskVelocityExponent)
     	*pow(expansionFactor                   ,starFormationExpansionExponent   );
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
      //@   <attachedTo>file</attachedTo>
      //@   <description>
      //@     The timescale (in Gyr) for star formation in the \cite{baugh_can_2005} prescription.
      //@   </description>
      //@   <type>real</type>
      //@   <cardinality>1</cardinality>
      //@   <group>starFormation</group>
      //@ </inputParameter>
      Get_Input_Parameter(strlen(timescaleParameter),timescaleParameter,&starFormationDiskTimescale       ,&starFormationDiskTimescaleDefault       );
      //@ <inputParameter>
      //@   <name>starFormationDiskVelocityExponent</name>
      //@   <defaultValue>-3.0</defaultValue>
      //@   <attachedTo>file</attachedTo>
      //@   <description>
      //@     The exponent for disk velocity in the \cite{baugh_can_2005} prescription for star formation in galactic disks.
      //@   </description>
      //@   <type>real</type>
      //@   <cardinality>1</cardinality>
      //@   <group>starFormation</group>
      //@ </inputParameter>
      Get_Input_Parameter(strlen(velocityParameter ),velocityParameter ,&starFormationDiskVelocityExponent,&starFormationDiskVelocityExponentDefault);
      //@ <inputParameter>
      //@   <name>starFormationExpansionExponent</name>
      //@   <defaultValue>0.0</defaultValue>
      //@   <attachedTo>file</attachedTo>
      //@   <description>
      //@     The exponent for expansion factor in the \cite{baugh_can_2005} prescription for star formation in galactic disks.
      //@   </description>
      //@   <type>real</type>
      //@   <cardinality>1</cardinality>
      //@   <group>starFormation</group>
      //@ </inputParameter>
      Get_Input_Parameter(strlen(expansionParameter),expansionParameter,&starFormationExpansionExponent   ,&starFormationExpansionExponentDefault   );
    }

  }
  
}
