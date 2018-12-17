!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implementation of a timescale for star formation in galactic disks which scales with the circular velocity of the host halo.

  use Kind_Numbers
  use Math_Exponentiation
  use Cosmology_Functions
  use Dark_Matter_Profiles

  !# <starFormationTimescaleDisks name="starFormationTimescaleDisksVelocityMaxScaling" defaultThreadPrivate="yes">
  !#  <description>A velocityMaxScaling timescale for star formation in galactic disks.</description>
  !# </starFormationTimescaleDisks>
  type, extends(starFormationTimescaleDisksClass) :: starFormationTimescaleDisksVelocityMaxScaling
     !% Implementation of a velocityMaxScaling timescale for star formation in galactic disks.
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_
     class           (darkMatterProfileClass ), pointer :: darkMatterProfile_
     double precision                                   :: expansionFactorFactorPrevious, exponentVelocity            , &
          &                                                exponentRedshift             , timescaleNormalization      , &
          &                                                timescaleStored              , velocityMaximumPrevious     , &
          &                                                velocityFactorPrevious       , expansionFactorPrevious
     logical                                            :: timescaleComputed
     integer         (kind_int8                       ) :: lastUniqueID
     type            (fastExponentiator               ) :: velocityExponentiator        , expansionFactorExponentiator
   contains
     final     ::                     velocityMaxScalingDestructor
     procedure :: timescale        => velocityMaxScalingTimescale
     procedure :: calculationReset => velocityMaxScalingCalculationReset
  end type starFormationTimescaleDisksVelocityMaxScaling

  interface starFormationTimescaleDisksVelocityMaxScaling
     !% Constructors for the {\normalfont \ttfamily velocityMaxScaling} timescale for star formation in disks class.
     module procedure velocityMaxScalingConstructorParameters
     module procedure velocityMaxScalingConstructorInternal
  end interface starFormationTimescaleDisksVelocityMaxScaling

  double precision, parameter :: velocityMaxScalingVelocityNormalization=200.0d0

contains

  function velocityMaxScalingConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily velocityMaxScaling} timescale for star formation in disks class which takes a
    !% parameter set as input.
    use Input_Parameters
    implicit none
    type            (starFormationTimescaleDisksVelocityMaxScaling)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                      ), pointer       :: cosmologyFunctions_
    class           (darkMatterProfileClass                       ), pointer       :: darkMatterProfile_
    double precision                                                               :: timescale           , exponentVelocity, &
         &                                                                            exponentRedshift
    
    ! Get parameters of for the timescale calculation.
    !# <inputParameter>
    !#   <name>timescale</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The timescale for star formation in the velocity maximum scaling timescale model for disks.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentVelocity</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The exponent of virial velocity in the timescale for star formation in the velocity maximum scaling timescale model for disks.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentRedshift</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The exponent of redshift in the timescale for star formation in the velocity maximum scaling timescale model for disks.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !# <objectBuilder class="darkMatterProfile"  name="darkMatterProfile_"  source="parameters"/>
    self=starFormationTimescaleDisksVelocityMaxScaling(timescale,exponentVelocity,exponentRedshift,cosmologyFunctions_,darkMatterProfile_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function velocityMaxScalingConstructorParameters

  function velocityMaxScalingConstructorInternal(timescale,exponentVelocity,exponentRedshift,cosmologyFunctions_,darkMatterProfile_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily velocityMaxScaling} timescale for star formation in disks class.
    implicit none
    type            (starFormationTimescaleDisksVelocityMaxScaling)                        :: self
    double precision                                               , intent(in   )         :: timescale          , exponentVelocity, &
         &                                                                                    exponentRedshift
    class           (cosmologyFunctionsClass                      ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterProfileClass                       ), intent(in   ), target :: darkMatterProfile_
    !# <constructorAssign variables="exponentVelocity, exponentRedshift, *cosmologyFunctions_, *darkMatterProfile_"/>

    self%lastUniqueID                 =-1_kind_int8
    self%timescaleComputed            =.false.
    self%velocityMaximumPrevious      =-1.0d0
    self%velocityFactorPrevious       =-1.0d0
    self%expansionFactorPrevious      =-1.0d0
    self%expansionFactorFactorPrevious=-1.0d0
    ! Compute the normalization of the timescale.
    self%timeScaleNormalization=+timescale                                                     &
         &                      /velocityMaxScalingVelocityNormalization**self%exponentVelocity
    ! Initialize exponentiators.
    self%velocityExponentiator       =fastExponentiator(1.0d+0,1.0d+3,self%exponentVelocity,1.0d+1,abortOutsideRange=.false.)
    self%expansionFactorExponentiator=fastExponentiator(1.0d-3,1.0d+0,self%exponentRedshift,1.0d+3,abortOutsideRange=.false.)
    return
  end function velocityMaxScalingConstructorInternal

  subroutine velocityMaxScalingDestructor(self)
    !% Destructor for the {\normalfont \ttfamily velocityMaxScaling} timescale for star formation in disks class.
    implicit none
    type(starFormationTimescaleDisksVelocityMaxScaling), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_" />
    !# <objectDestructor name="self%darkMatterProfile_"/>
    return
  end subroutine velocityMaxScalingDestructor

  subroutine velocityMaxScalingCalculationReset(self,node)
    !% Reset the halo scaling disk star formation timescale calculation.
    implicit none
    class(starFormationTimescaleDisksVelocityMaxScaling), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node

    self%timescaleComputed=.false.
    self%lastUniqueID     =node%uniqueID()
    return
  end subroutine velocityMaxScalingCalculationReset

  double precision function velocityMaxScalingTimescale(self,node)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\normalfont \ttfamily node} in the halo scaling timescale model.
    use Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (starFormationTimescaleDisksVelocityMaxScaling), intent(inout), target :: self
    type            (treeNode                                     ), intent(inout), target :: node
    class           (nodeComponentBasic                           ), pointer               :: basic
    double precision                                                                       :: expansionFactor, velocityMaximum
    
    
    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Compute the timescale if necessary.
    if (.not.self%timescaleComputed) then
       ! Get virial velocity and expansion factor.
       basic           => node%basic                                      (            )
       velocityMaximum =  self%darkMatterProfile_ %circularVelocityMaximum(node        )
       expansionFactor =  self%cosmologyFunctions_%expansionFactor        (basic%time())
       ! Compute the velocity factor.
       if (velocityMaximum /= self%velocityMaximumPrevious) then
           self%velocityMaximumPrevious=velocityMaximum
           self%velocityFactorPrevious =velocityMaximum**self%exponentVelocity
       end if
       ! Compute the expansion-factor factor.
       if (expansionFactor /= self%expansionFactorPrevious) then
          self%expansionFactorPrevious      =      expansionFactor
          self%expansionFactorFactorPrevious=1.0d0/expansionFactor**self%exponentRedshift
       end if
       ! Computed the timescale.
       self%timescaleStored=+self%timeScaleNormalization        &
            &               *self%velocityFactorPrevious        &
            &               *self%expansionFactorFactorPrevious
       ! Record that the timescale is now computed.
       self%timescaleComputed=.true.
    end if
    ! Return the stored timescale.
    velocityMaxScalingTimescale=self%timescaleStored
    return
  end function velocityMaxScalingTimescale
