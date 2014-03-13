!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

  !% An implementation of the intergalactic medium state class for a simplistic model of instantaneous and full reionization.

  !# <intergalacticMediumState name="intergalacticMediumStateSimple">
  !#  <description>The intergalactic medium is assumed to be instantaneously and fully reionized at a fixed redshift, and heated to a fixed temperature.</description>
  !# </intergalacticMediumState>

  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateSimple
     !% An \gls{igm} state class for a simple model in which the \gls{igm} is assumed to be instantaneously and fully reionized at
     !% a fixed redshift, and heated to a fixed temperature.
     private
     double precision :: reionizationTime, reionizationTemperature, preReionizationTemperature
   contains
     procedure :: electronFraction            => simpleElectronFraction
     procedure :: temperature                 => simpleTemperature
     procedure :: neutralHydrogenFraction     => simpleNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => simpleNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => simpleSinglyIonizedHeliumFraction
  end type intergalacticMediumStateSimple

  interface intergalacticMediumStateSimple
     !% Constructors for the simple intergalactic medium state class.
     module procedure simpleDefaultConstructor
     module procedure simpleConstructor
  end interface intergalacticMediumStateSimple

contains

  function simpleDefaultConstructor()
    !% Default constructor for the simple \gls{igm} state class.
    use Input_Parameters
    implicit none
    type            (intergalacticMediumStateSimple), target  :: simpleDefaultConstructor
    double precision                                          :: reionizationRedshift      , reionizationTemperature, &
         &                                                       preReionizationTemperature
    !@ <inputParameter>
    !@   <name>igmStateSimpleReionizationRedshift</name>
    !@   <defaultValue>9.97 (\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The redshift of reionization in the simple \gls{igm} state model.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('igmStateSimpleReionizationRedshift',reionizationRedshift,defaultValue=9.97d0)
    !@ <inputParameter>
    !@   <name>igmStateSimpleReionizationTemperature</name>
    !@   <defaultValue>$10^4$K</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The post-reionization temperature (in units of Kelvin) in the simple \gls{igm} state model.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('igmStateSimpleReionizationTemperature',reionizationTemperature,defaultValue=1.0d4)
    !@ <inputParameter>
    !@   <name>igmStateSimplePreReionizationTemperature</name>
    !@   <defaultValue>$10$K</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The pre-reionization temperature (in units of Kelvin) in the simple \gls{igm} state model.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('igmStateSimplePreReionizationTemperature',preReionizationTemperature,defaultValue=1.0d1)
    ! Construct the object.
    simpleDefaultConstructor=simpleConstructor(reionizationRedshift,reionizationTemperature,preReionizationTemperature)
    return
  end function simpleDefaultConstructor

  function simpleConstructor(reionizationRedshift,reionizationTemperature,preReionizationTemperature)
    !% Constructor for the simple \gls{igm} state class.
    use Cosmology_Functions
    implicit none
    type            (intergalacticMediumStateSimple), target        :: simpleConstructor
    double precision                                , intent(in   ) :: reionizationRedshift      , reionizationTemperature, &
         &                                                             preReionizationTemperature
    class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_

    cosmologyFunctions_                          => cosmologyFunctions()
    simpleConstructor%   reionizationTime        =  cosmologyFunctions_ %cosmicTime                 (                      &
         &                                           cosmologyFunctions_%expansionFactorFromRedshift (                     &
         &                                                                                            reionizationRedshift &
         &                                                                                           )                     &
         &                                                                                          )
    simpleConstructor%   reionizationTemperature =     reionizationTemperature
    simpleConstructor%preReionizationTemperature =  preReionizationTemperature
    return
  end function simpleConstructor
  
  double precision function simpleElectronFraction(self,time)
    !% Return the electron fraction of the \gls{igm} in the simple model.
    use Numerical_Constants_Astronomical
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    if (time > self%reionizationTime) then
       simpleElectronFraction=+      hydrogenByMassPrimordial                                  /atomicMassHydrogen &
            &                 +2.0d0*  heliumByMassPrimordial                                  /atomicMassHelium
    else
       simpleElectronFraction=0.0d0
    end if
    return
  end function simpleElectronFraction

  double precision function simpleNeutralHydrogenFraction(self,time)
    !% Return the neutral hydrogen fraction of the \gls{igm} in the simple model.
    use Numerical_Constants_Astronomical
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    if (time > self%reionizationTime) then
       simpleNeutralHydrogenFraction=0.0d0
    else
       simpleNeutralHydrogenFraction=1.0d0
    end if
    return
  end function simpleNeutralHydrogenFraction

  double precision function simpleNeutralHeliumFraction(self,time)
    !% Return the neutral helium fraction of the \gls{igm} in the simple model.
    use Numerical_Constants_Astronomical
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    if (time > self%reionizationTime) then
       simpleNeutralHeliumFraction=0.0d0
    else
       simpleNeutralHeliumFraction=1.0d0
    end if
    return
  end function simpleNeutralHeliumFraction

  double precision function simpleSinglyIonizedHeliumFraction(self,time)
    !% Return the singly-ionized helium fraction of the \gls{igm} in the simple model.
    use Numerical_Constants_Astronomical
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    if (time > self%reionizationTime) then
       simpleSinglyIonizedHeliumFraction=0.0d0
    else
       simpleSinglyIonizedHeliumFraction=1.0d0
    end if
    return
  end function simpleSinglyIonizedHeliumFraction

  double precision function simpleTemperature(self,time)
    !% Return the temperature of the \gls{igm} in the simple model.
    implicit none
    class           (intergalacticMediumStateSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: time
 
    if (time > self%reionizationTime) then
       simpleTemperature=self%   reionizationTemperature
    else
       simpleTemperature=self%preReionizationTemperature
    end if
    return
  end function simpleTemperature
