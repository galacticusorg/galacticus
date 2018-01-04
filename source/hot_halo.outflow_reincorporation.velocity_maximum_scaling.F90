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

!% An implementation of the hot halo outflow reincorporation class which uses simple scalings based on the halo maximum circular
!% velocity.

  use Math_Exponentiation

  !# <hotHaloOutflowReincorporation name="hotHaloOutflowReincorporationVelocityMaximumScaling">
  !#  <description>An implementation of the hot halo outflow reincorporation class which uses simple scalings based on the halo maximum circular velocity.</description>
  !# </hotHaloOutflowReincorporation>
  type, extends(hotHaloOutflowReincorporationClass) :: hotHaloOutflowReincorporationVelocityMaximumScaling
     !% An implementation of the hot halo outflow reincorporation class which uses simple scalings based on the halo maximum circular velocity.
     private
     double precision                    :: timeScaleNormalization , velocityExponent            , &
          &                                 redshiftExponent       , velocityMaximumFactor       , &
          &                                 expansionFactorFactor  , rateStored                  , &
          &                                 timeScaleMinimum
     logical                             :: velocityMaximumComputed, expansionFactorComputed     , &
          &                                 rateComputed
     integer         (kind=kind_int8   ) :: lastUniqueID
     ! Fast exponentiation tables for rapid computation of the outflow rate.
     type            (fastExponentiator) :: velocityExponentiator  , expansionFactorExponentiator
   contains
     procedure :: calculationReset => velocityMaximumScalingCalculationReset
     procedure :: rate             => velocityMaximumScalingRate
  end type hotHaloOutflowReincorporationVelocityMaximumScaling

  interface hotHaloOutflowReincorporationVelocityMaximumScaling
     !% Constructors for the {\normalfont \ttfamily velocityMaximumScaling} hot halo outflow reincorporation class.
     module procedure velocityMaximumScalingDefaultConstructor
     module procedure velocityMaximumScalingConstructor
  end interface hotHaloOutflowReincorporationVelocityMaximumScaling

  ! Initialization state.
  logical                     :: velocityMaximumScalingInitialized       =.false.
  logical                     :: velocityMaximumScalingDefaultInitialized=.false.

  ! Normalization parameter.
  double precision, parameter :: velocityMaximumScalingVelocityNormalization=200.0d0

  ! Default parameters.
  double precision            :: velocityMaximumScalingVelocityExponent, velocityMaximumScalingRedshiftExponent, &
       &                         velocityMaximumScalingTimeScale       , velocityMaximumScalingTimescaleMinimum

contains

  function velocityMaximumScalingDefaultConstructor() result(self)
    !% Default constructor for the velocityMaximumScaling hot halo outflow reincorporation class.
    use Input_Parameters
    implicit none
    type(hotHaloOutflowReincorporationVelocityMaximumScaling) :: self
    
    if (.not.velocityMaximumScalingDefaultInitialized) then
       !$omp critical(hotHaloOutflowReincorporationVMSDefaultInitialize)
       if (.not.velocityMaximumScalingDefaultInitialized) then
          !# <inputParameter>
          !#   <name>hotHaloOutflowReincorporationVlctyMxSclngTimeScale</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>1.0d0</defaultValue>
          !#   <description>The timescale in the velocity maximum scaling model for outflow reincorporation.</description>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !#   <variable>velocityMaximumScalingTimeScale</variable>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>hotHaloOutflowReincorporationVlctyMxSclngVelocityExponent</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>0.0d0</defaultValue>
          !#   <description>The exponent of maximum circular velocity in the velocity maximum scaling model for outflow reincorporation.</description>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !#   <variable>velocityMaximumScalingVelocityExponent</variable>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>hotHaloOutflowReincorporationVlctyMxSclngRedshiftExponent</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>-1.5d0</defaultValue>
          !#   <description>The exponent of redshift in the velocity maximum scaling model for outflow reincorporation.</description>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !#   <variable>velocityMaximumScalingRedshiftExponent</variable>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>hotHaloOutflowReincorporationVlctyMxSclngTimescaleMinimum</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>1.0d-3</defaultValue>
          !#   <description>The minimum timescale for outflow reincorporation in the velocity maximum scaling model.</description>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !#   <variable>velocityMaximumScalingTimescaleMinimum</variable>
          !# </inputParameter>
          ! Record that class is initialized.
          velocityMaximumScalingDefaultInitialized=.true.
       end if
       !$omp end critical(hotHaloOutflowReincorporationVMSDefaultInitialize)
    end if
    self=hotHaloOutflowReincorporationVelocityMaximumScaling(                                        &
         &                                                   velocityMaximumScalingTimeScale       , &
         &                                                   velocityMaximumScalingTimescaleMinimum, &
         &                                                   velocityMaximumScalingVelocityExponent, &
         &                                                   velocityMaximumScalingRedshiftExponent  &
         &                                                  )
    return
  end function velocityMaximumScalingDefaultConstructor
  
  function velocityMaximumScalingConstructor(timeScale,timeScaleMinimum,velocityExponent,redshiftExponent) result(self)
    !% Default constructor for the velocityMaximumScaling hot halo outflow reincorporation class.
    implicit none
    type            (hotHaloOutflowReincorporationVelocityMaximumScaling)                :: self
    double precision                                                     , intent(in   ) :: timeScale       , velocityExponent, &
         &                                                                                  redshiftExponent, timeScaleMinimum
    
    ! Initialize.
    call velocityMaximumScalingInitalize()
    ! Construct the object.
    self%timeScaleNormalization =timeScale/velocityMaximumScalingVelocityNormalization**velocityExponent
    self%timeScaleMinimum       =timeScaleMinimum
    self%velocityExponent       =velocityExponent
    self%redshiftExponent       =redshiftExponent
    self%lastUniqueID           =-1_kind_int8
    self%velocityMaximumFactor  =-1.0d0
    self%expansionFactorFactor  =-1.0d0
    self%rateStored             =-1.0d0
    self%velocityMaximumComputed=.false.
    self%expansionFactorComputed=.false.
    self%rateComputed           =.false.
    ! Initialize exponentiators.
    self%velocityExponentiator       =fastExponentiator(1.0d+0,1.0d+3,self%velocityExponent,1.0d+1,abortOutsideRange=.false.)
    self%expansionFactorExponentiator=fastExponentiator(1.0d-3,1.0d+0,self%redshiftExponent,1.0d+3,abortOutsideRange=.false.)
    return
  end function velocityMaximumScalingConstructor

  subroutine velocityMaximumScalingInitalize()
    !% Initialize the {\normalfont \ttfamily velocityMaximumScaling} hot halo outflow reincorporation class.
    use Galacticus_Error
    implicit none
    
    if (.not.velocityMaximumScalingInitialized) then
       !$omp critical(hotHaloOutflowReincorporationVelocityMaximumScalingInitialize)
       if (.not.velocityMaximumScalingInitialized) then
          if (.not.defaultHotHaloComponent%outflowedMassIsGettable())                                                       &
               & call Galacticus_Error_Report                                                                               &
               &   (                                                                                                        &
               &    'the "outflowedMass" properties of the hotHalo component must be gettable.'                          // &
               &    Galacticus_Component_List(                                                                              &
               &                              'hotHalo'                                                                   , &
               &                              defaultHotHaloComponent%outflowedMassAttributeMatch(requireGettable=.true.)   &
               &                             )                                                                           // &
               &    {introspection:location}                                                                                &
               &   )
          ! Record that class is initialized.
          velocityMaximumScalingInitialized=.true.
       end if
       !$omp end critical(hotHaloOutflowReincorporationVelocityMaximumScalingInitialize)
    end if
    return
  end subroutine velocityMaximumScalingInitalize
  
  subroutine velocityMaximumScalingCalculationReset(self,node)
    !% Reset the halo scales calculation.
    implicit none
    class(hotHaloOutflowReincorporationVelocityMaximumScaling), intent(inout) :: self
    type (treeNode                                           ), intent(inout) :: node

    self%velocityMaximumComputed=.false.
    self%expansionFactorComputed=.false.
    self%rateComputed           =.false.
    self%lastUniqueID           =node%uniqueID()
    return
  end subroutine velocityMaximumScalingCalculationReset

  double precision function velocityMaximumScalingRate(self,node)
    !% Return the rate of mass reincorporation for outflowed gas in the hot halo.
    use Cosmology_Functions
    use Dark_Matter_Profiles
    implicit none
    class           (hotHaloOutflowReincorporationVelocityMaximumScaling), intent(inout) :: self
    type            (treeNode                                           ), intent(inout) :: node
    class           (cosmologyFunctionsClass                            ), pointer       :: cosmologyFunctions_
    class           (darkMatterProfileClass                             ), pointer       :: darkMatterProfile_
    class           (nodeComponentBasic                                 ), pointer       :: basic
    class           (nodeComponentHotHalo                               ), pointer       :: hotHalo
    double precision                                                                     :: timeScale

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Get required components.
    ! Compute velocity maximum factor.
    if (.not.self%velocityMaximumComputed) then
       darkMatterProfile_           =>                                         darkMatterProfile                          (    )
       self%velocityMaximumFactor   =  self%velocityExponentiator%exponentiate(darkMatterProfile_ %circularVelocityMaximum(node))
       self%velocityMaximumComputed = .true.
    end if
    ! Compute expansion factor factor.
    if (.not.self%expansionFactorComputed) then
       basic                        =>                                                      node               %basic          (            )
       cosmologyFunctions_          =>                                                      cosmologyFunctions                 (            )
       self%expansionFactorFactor   =  1.0d0/self%expansionFactorExponentiator%exponentiate(cosmologyFunctions_%expansionFactor(basic%time()))
       self%expansionFactorComputed = .true.
    end if
    ! Compute the rate.
    if (.not.self%rateComputed) then
       hotHalo          =>  node    %hotHalo               ()
       timeScale        =  max(                                &
            &                  +self%timeScaleNormalization    &
            &                  *self%velocityMaximumFactor     &
            &                  *self%expansionFactorFactor   , &
            &                   self%timeScaleMinimum          &
            &                 )
       self%rateStored  =  +hotHalo %outflowedMass         ()  &
            &              /timeScale
       self%rateComputed=  .true.
    end if
    ! Return the rate.
    velocityMaximumScalingRate=self%rateStored
    return
  end function velocityMaximumScalingRate
