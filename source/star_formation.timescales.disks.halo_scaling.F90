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

  !% Implementation of a timescale for star formation feedback in galactic disks which scales with the circular velocity of the host halo.

  use Cosmology_Functions
  use Dark_Matter_Halo_Scales

  !# <starFormationTimescaleDisks name="starFormationTimescaleDisksHaloScaling" defaultThreadPrivate="yes">
  !#  <description>A haloScaling timescale for star formation feedback in galactic disks.</description>
  !# </starFormationTimescaleDisks>
  type, extends(starFormationTimescaleDisksClass) :: starFormationTimescaleDisksHaloScaling
     !% Implementation of a haloScaling timescale for star formation feedback in galactic disks.
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
     double precision                                    :: expansionFactorFactorPrevious, exponentVelocityVirial , &
          &                                                 exponentRedshift             , timescaleNormalization , &
          &                                                 timescaleStored              , velocityPrevious       , &
          &                                                 velocityFactorPrevious       , expansionFactorPrevious
     logical                                             :: timescaleComputed
     integer         (kind_int8                        ) :: lastUniqueID
   contains
     final     ::                     haloScalingDestructor
     procedure :: timescale        => haloScalingTimescale
     procedure :: calculationReset => haloScalingCalculationReset
  end type starFormationTimescaleDisksHaloScaling

  interface starFormationTimescaleDisksHaloScaling
     !% Constructors for the {\normalfont \ttfamily haloScaling} timescale for star formation in disks class.
     module procedure haloScalingConstructorParameters
     module procedure haloScalingConstructorInternal
  end interface starFormationTimescaleDisksHaloScaling

  double precision, parameter :: haloScalingVelocityVirialNormalization=200.0d0

contains

  function haloScalingConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily haloScaling} timescale for star formation feedback in disks class which takes a
    !% parameter set as input.
    use Input_Parameters
    implicit none
    type            (starFormationTimescaleDisksHaloScaling)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    double precision                                                        :: timescale           , exponentVelocityVirial, &
         &                                                                     exponentRedshift

    !# <inputParameter>
    !#   <name>timescale</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The timescale for star formation in the halo scaling timescale model for disks.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentVelocityVirial</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The exponent of virial velocity in the timescale for star formation in the halo scaling timescale model for disks.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentRedshift</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The exponent of redshift in the timescale for star formation in the halo scaling timescale model for disks.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=starFormationTimescaleDisksHaloScaling(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function haloScalingConstructorParameters

  function haloScalingConstructorInternal(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily haloScaling} timescale for star formation in disks class.
    implicit none
    type            (starFormationTimescaleDisksHaloScaling)                        :: self
    double precision                                        , intent(in   )         :: timescale           , exponentVelocityVirial, &
         &                                                                             exponentRedshift
    class           (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_
    !# <constructorAssign variables="exponentVelocityVirial, exponentRedshift, *cosmologyFunctions_, *darkMatterHaloScale_"/>
    
    self%lastUniqueID                 =-1_kind_int8
    self%timescaleComputed            =.false.
    self%velocityPrevious             =-1.0d0
    self%velocityFactorPrevious       =-1.0d0
    self%expansionFactorPrevious      =-1.0d0
    self%expansionFactorFactorPrevious=-1.0d0
    ! Compute the normalization of the timescale.
    self%timeScaleNormalization=+timescale                                                           &
         &                      /haloScalingVelocityVirialNormalization**self%exponentVelocityVirial
    return
  end function haloScalingConstructorInternal

  subroutine haloScalingDestructor(self)
    !% Destructor for the {\normalfont \ttfamily haloScaling} timescale for star formation in disks class.
    implicit none
    type(starFormationTimescaleDisksHaloScaling), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_" />
    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine haloScalingDestructor

  subroutine haloScalingCalculationReset(self,node)
    !% Reset the halo scaling disk star formation timescale calculation.
    implicit none
    class(starFormationTimescaleDisksHaloScaling), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node

    self%timescaleComputed=.false.
    self%lastUniqueID     =node%uniqueID()
    return
  end subroutine haloScalingCalculationReset

  double precision function haloScalingTimescale(self,node)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\normalfont \ttfamily node} in the halo scaling timescale model.
    use Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (starFormationTimescaleDisksHaloScaling), intent(inout), target :: self
    type            (treeNode                              ), intent(inout), target :: node
    class           (nodeComponentBasic                    ), pointer               :: basic
    double precision                                                                :: expansionFactor, velocityVirial
    
    
    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Compute the timescale if necessary.
    if (.not.self%timescaleComputed) then
       ! Get virial velocity and expansion factor.
       basic           => node%basic                               (            )
       velocityVirial  =  self%darkMatterHaloScale_%virialVelocity (node        )
       expansionFactor =  self%cosmologyFunctions_ %expansionFactor(basic%time())
       ! Compute the velocity factor.
       if (velocityVirial /= self%velocityPrevious) then
           self%velocityPrevious      =velocityVirial
           self%velocityFactorPrevious=velocityVirial**self%exponentVelocityVirial
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
    haloScalingTimescale=self%timescaleStored
    return
  end function haloScalingTimescale
