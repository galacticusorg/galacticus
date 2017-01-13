!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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
  
  !# <hotHaloOutflowReincorporation name="hotHaloOutflowReincorporationVelocityMaximumScaling">
  !#  <description>An implementation of the hot halo outflow reincorporation class which uses simple scalings based on the halo maximum circular velocity.</description>
  !# </hotHaloOutflowReincorporation>
  type, extends(hotHaloOutflowReincorporationClass) :: hotHaloOutflowReincorporationVelocityMaximumScaling
     !% An implementation of the hot halo outflow reincorporation class which uses simple scalings based on the halo maximum circular velocity.
     private
     double precision                 :: timeScaleNormalization , velocityExponent       , &
          &                              redshiftExponent       , velocityMaximumFactor  , &
          &                              expansionFactorFactor  , rateStored             , &
          &                              timeScaleMinimum
     logical                          :: velocityMaximumComputed, expansionFactorComputed, &
          &                              rateComputed
     integer         (kind=kind_int8) :: lastUniqueID
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

  function velocityMaximumScalingDefaultConstructor()
    !% Default constructor for the velocityMaximumScaling hot halo outflow reincorporation class.
    use Input_Parameters
    implicit none
    type(hotHaloOutflowReincorporationVelocityMaximumScaling) :: velocityMaximumScalingDefaultConstructor
    
    if (.not.velocityMaximumScalingDefaultInitialized) then
       !$omp critical(hotHaloOutflowReincorporationVMSDefaultInitialize)
       if (.not.velocityMaximumScalingDefaultInitialized) then
          !@ <inputParameter>
          !@   <name>hotHaloOutflowReincorporationVlctyMxSclngTimeScale</name>
          !@   <defaultValue>$1$~Gyr</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The timescale in the velocity maximum scaling model for outflow reincorporation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloOutflowReincorporationVlctyMxSclngTimeScale',velocityMaximumScalingTimeScale,defaultValue=1.0d0)
          !@ <inputParameter>
          !@   <name>hotHaloOutflowReincorporationVlctyMxSclngVelocityExponent</name>
          !@   <defaultValue>$0$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The exponent of maximum circular velocity in the velocity maximum scaling model for outflow reincorporation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloOutflowReincorporationVlctyMxSclngVelocityExponent',velocityMaximumScalingVelocityExponent,defaultValue=0.0d0)
          !@ <inputParameter>
          !@   <name>hotHaloOutflowReincorporationVlctyMxSclngRedshiftExponent</name>
          !@   <defaultValue>$-1.5$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The exponent of redshift in the velocity maximum scaling model for outflow reincorporation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloOutflowReincorporationVlctyMxSclngRedshiftExponent',velocityMaximumScalingRedshiftExponent,defaultValue=-1.5d0)
          !@ <inputParameter>
          !@   <name>hotHaloOutflowReincorporationVlctyMxSclngTimescaleMinimum</name>
          !@   <defaultValue>$0.001$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The minimum timescale for outflow reincorporation in the velocity maximum scaling model.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloOutflowReincorporationVlctyMxSclngTimescaleMinimum',velocityMaximumScalingTimescaleMinimum,defaultValue=1.0d-3)
          ! Record that class is initialized.
          velocityMaximumScalingDefaultInitialized=.true.
       end if
       !$omp end critical(hotHaloOutflowReincorporationVMSDefaultInitialize)
    end if
    velocityMaximumScalingDefaultConstructor=velocityMaximumScalingConstructor(                                        &
         &                                                                     velocityMaximumScalingTimeScale       , &
         &                                                                     velocityMaximumScalingTimescaleMinimum, &
         &                                                                     velocityMaximumScalingVelocityExponent, &
         &                                                                     velocityMaximumScalingRedshiftExponent  &
         &                                                                    )
    return
  end function velocityMaximumScalingDefaultConstructor
  
  function velocityMaximumScalingConstructor(timeScale,timeScaleMinimum,velocityExponent,redshiftExponent)
    !% Default constructor for the velocityMaximumScaling hot halo outflow reincorporation class.
    implicit none
    type            (hotHaloOutflowReincorporationVelocityMaximumScaling)                :: velocityMaximumScalingConstructor
    double precision                                                     , intent(in   ) :: timeScale                        , velocityExponent, &
         &                                                                                  redshiftExponent                 , timeScaleMinimum
    
    ! Initialize.
    call velocityMaximumScalingInitalize()
    ! Construct the object.
    velocityMaximumScalingConstructor%timeScaleNormalization =timeScale/velocityMaximumScalingVelocityNormalization**velocityExponent
    velocityMaximumScalingConstructor%timeScaleMinimum       =timeScaleMinimum
    velocityMaximumScalingConstructor%velocityExponent       =velocityExponent
    velocityMaximumScalingConstructor%redshiftExponent       =redshiftExponent
    velocityMaximumScalingConstructor%lastUniqueID           =-1_kind_int8
    velocityMaximumScalingConstructor%velocityMaximumFactor  =-1.0d0
    velocityMaximumScalingConstructor%expansionFactorFactor  =-1.0d0
    velocityMaximumScalingConstructor%rateStored             =-1.0d0
    velocityMaximumScalingConstructor%velocityMaximumComputed=.false.
    velocityMaximumScalingConstructor%expansionFactorComputed=.false.
    velocityMaximumScalingConstructor%rateComputed           =.false.
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
               &    'velocityMaximumScalingInitalize'                                                                     , &
               &    'the "outflowedMass" properties of the hotHalo component must be gettable.'                          // &
               &    Galacticus_Component_List(                                                                              &
               &                              'hotHalo'                                                                   , &
               &                              defaultHotHaloComponent%outflowedMassAttributeMatch(requireGettable=.true.)   &
               &                             )                                                                              &
               &   )
          ! Record that class is initialized.
          velocityMaximumScalingInitialized=.true.
       end if
       !$omp end critical(hotHaloOutflowReincorporationVelocityMaximumScalingInitialize)
    end if
    return
  end subroutine velocityMaximumScalingInitalize
  
  subroutine velocityMaximumScalingCalculationReset(self,thisNode)
    !% Reset the halo scales calculation.
    implicit none
    class(hotHaloOutflowReincorporationVelocityMaximumScaling), intent(inout) :: self
    type (treeNode                                           ), intent(inout) :: thisNode

    self%velocityMaximumComputed=.false.
    self%expansionFactorComputed=.false.
    self%rateComputed           =.false.
    self%lastUniqueID           =thisNode%uniqueID()
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
       darkMatterProfile_           => darkMatterProfile                          (    )
       self%velocityMaximumFactor   =  darkMatterProfile_ %circularVelocityMaximum(node)**self%velocityExponent
       self%velocityMaximumComputed = .true.
    end if
    ! Compute expansion factor factor.
    if (.not.self%expansionFactorComputed) then
       basic                        =>       node               %basic          (            )
       cosmologyFunctions_          =>       cosmologyFunctions                 (            )
       self%expansionFactorFactor   =  1.0d0/cosmologyFunctions_%expansionFactor(basic%time())**self%redshiftExponent
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
