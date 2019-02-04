!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !% Implementation of the \cite{baugh_can_2005} timescale for star formation in galactic disks

  use Cosmology_Functions

  !# <starFormationTimescaleDisks name="starFormationTimescaleDisksBaugh2005">
  !#  <description>The \cite{baugh_can_2005} timescale for star formation in galactic disks.</description>
  !# </starFormationTimescaleDisks>
  type, extends(starFormationTimescaleDisksClass) :: starFormationTimescaleDisksBaugh2005
     !% Implementation of the \cite{baugh_can_2005} timescale for star formation in galactic disks.
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
     double precision                                    :: timescaleValue         , exponentVelocity, &
          &                                                 exponentExpansionFactor
   contains
     final     ::              baugh2005Destructor
     procedure :: timescale => baugh2005Timescale
  end type starFormationTimescaleDisksBaugh2005

  interface starFormationTimescaleDisksBaugh2005
     !% Constructors for the {\normalfont \ttfamily baugh2005} timescale for star formation in disks class.
     module procedure baugh2005ConstructorParameters
     module procedure baugh2005ConstructorInternal
  end interface starFormationTimescaleDisksBaugh2005

  double precision, parameter :: baugh2005VelocityVirialNormalization=200.0d0

contains

  function baugh2005ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily baugh2005} timescale for star formation in disks class which takes a
    !% parameter set as input.
    use Input_Parameters
    implicit none
    type            (starFormationTimescaleDisksBaugh2005)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass             ), pointer       :: cosmologyFunctions_
    double precision                                                      :: timescale              , exponentVelocity, &
         &                                                                   exponentExpansionFactor

    !# <inputParameter>
    !#   <name>timescale</name>
    !#   <defaultValue>8.0d0</defaultValue>
    !#   <description>The timescale (in Gyr) for star formation in the \cite{baugh_can_2005} prescription.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !#   <group>starFormation</group>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentVelocity</name>
    !#   <defaultValue>-3.0d0</defaultValue>
    !#   <description>The exponent for disk velocity in the \cite{baugh_can_2005} prescription for star formation in galactic disks.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !#   <group>starFormation</group>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentExpansionFactor</name>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The exponent for expansion factor in the \cite{baugh_can_2005} prescription for star formation in galactic disks.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !#   <group>starFormation</group>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    self=starFormationTimescaleDisksBaugh2005(timescale,exponentVelocity,exponentExpansionFactor,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    return
  end function baugh2005ConstructorParameters

  function baugh2005ConstructorInternal(timescale,exponentVelocity,exponentExpansionFactor,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily baugh2005} timescale for star formation in disks class.
    implicit none
    type            (starFormationTimescaleDisksBaugh2005)                          :: self
    double precision                                        , intent(in   )         :: timescale              , exponentVelocity, &
         &                                                                             exponentExpansionFactor
    class           (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="exponentVelocity, exponentExpansionFactor, *cosmologyFunctions_"/>

    self%timescaleValue=timescale
    return
  end function baugh2005ConstructorInternal

  subroutine baugh2005Destructor(self)
    !% Destructor for the {\normalfont \ttfamily baugh2005} timescale for star formation in disks class.
    implicit none
    type(starFormationTimescaleDisksBaugh2005), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end subroutine baugh2005Destructor

  double precision function baugh2005Timescale(self,node)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\normalfont \ttfamily node} in the halo scaling timescale model.
    use Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk
    implicit none
    class           (starFormationTimescaleDisksBaugh2005), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node 
    class           (nodeComponentBasic                  ), pointer               :: basic
    class           (nodeComponentDisk                   ), pointer               :: disk
    double precision                                                              :: expansionFactor, velocityDisk, &
         &                                                                           time

    basic        => node%basic   ()
    disk         => node%disk    ()
    velocityDisk =  disk%velocity()
    if (velocityDisk <= 0.0) then
       baugh2005Timescale=0.0d0
    else
       time              =+basic%time                               (    )
       expansionFactor   =+self %cosmologyFunctions_%expansionFactor(time)
       baugh2005Timescale=+self %timescaleValue                                                                 &
            &             *(velocityDisk   /baugh2005VelocityVirialNormalization)**self%exponentVelocity        &
            &             * expansionFactor                                      **self%exponentExpansionFactor
    end if
    return
  end function baugh2005Timescale
