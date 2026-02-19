!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
An implementation of the hot halo outflow reincorporation class which uses simple scalings based on the halo maximum circular
velocity.
!!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Kind_Numbers            , only : kind_int8
  use :: Math_Exponentiation     , only : fastExponentiator

  !![
  <hotHaloOutflowReincorporation name="hotHaloOutflowReincorporationVelocityMaximumScaling">
   <description>
    A hot halo outflow reincorporation class which uses simple scalings based on the halo maximum circular
    velocity. Specifically,
    \begin{equation}
     \dot{M}_\mathrm{reincorporation} = M_\mathrm{outflowed} / t_\mathrm{reincorporation},
    \end{equation}
    where
    \begin{equation}
     t_\mathrm{reincorporation} = \tau_\mathrm{reincorporation} \left( { V_\mathrm{max} \over 200 \hbox{km/s}}
     \right)^{\alpha_\mathrm{reincorporation}} (1+z)^{\beta_\mathrm{reincorporation}},
    \end{equation}
    where $\tau_\mathrm{reincorporation}=${\normalfont \ttfamily [timeScale]}, $\alpha_\mathrm{reincorporation}=${\normalfont
    \ttfamily [velocityExponent]}, and $\beta_\mathrm{reincorporation}=${\normalfont \ttfamily [redshiftExponent]}.
   </description>
  </hotHaloOutflowReincorporation>
  !!]
  type, extends(hotHaloOutflowReincorporationClass) :: hotHaloOutflowReincorporationVelocityMaximumScaling
     !!{
     An implementation of the hot halo outflow reincorporation class which uses simple scalings based on the halo maximum circular velocity.
     !!}
     private
     class           (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_     => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_   => null()
     double precision                                     :: timeScaleNormalization           , velocityExponent            , &
          &                                                  redshiftExponent                 , velocityMaximumFactor       , &
          &                                                  expansionFactorFactor            , rateStored                  , &
          &                                                  timeScaleMinimum                 , timeScale
     logical                                              :: velocityMaximumComputed          , expansionFactorComputed     , &
          &                                                  rateComputed
     integer         (kind=kind_int8            )         :: lastUniqueID
     ! Fast exponentiation tables for rapid computation of the outflow rate.
     type            (fastExponentiator         )         :: velocityExponentiator            , expansionFactorExponentiator
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                     velocityMaximumScalingDestructor
     procedure :: autoHook         => velocityMaximumScalingAutoHook
     procedure :: calculationReset => velocityMaximumScalingCalculationReset
     procedure :: rate             => velocityMaximumScalingRate
  end type hotHaloOutflowReincorporationVelocityMaximumScaling

  interface hotHaloOutflowReincorporationVelocityMaximumScaling
     !!{
     Constructors for the \refClass{hotHaloOutflowReincorporationVelocityMaximumScaling} hot halo outflow reincorporation class.
     !!}
     module procedure velocityMaximumScalingConstructorParameters
     module procedure velocityMaximumScalingConstructorInternal
  end interface hotHaloOutflowReincorporationVelocityMaximumScaling

  ! Normalization parameter.
  double precision, parameter :: velocityNormalization=200.0d0

contains

  function velocityMaximumScalingConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily velocityMaximumScaling} hot halo outflow reincorporation class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hotHaloOutflowReincorporationVelocityMaximumScaling)                :: self
    type            (inputParameters                                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                            ), pointer       :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass                          ), pointer       :: darkMatterProfileDMO_
    double precision                                                                     :: timeScale            , velocityExponent, &
         &                                                                                  redshiftExponent     , timeScaleMinimum

    !![
    <inputParameter>
      <name>timeScale</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The timescale in the velocity maximum scaling model for outflow reincorporation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>velocityExponent</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent of maximum circular velocity in the velocity maximum scaling model for outflow reincorporation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftExponent</name>
      <defaultValue>-1.5d0</defaultValue>
      <description>The exponent of redshift in the velocity maximum scaling model for outflow reincorporation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timescaleMinimum</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The minimum timescale for outflow reincorporation in the velocity maximum scaling model.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=hotHaloOutflowReincorporationVelocityMaximumScaling(timeScale,timescaleMinimum,velocityExponent,redshiftExponent,cosmologyFunctions_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function velocityMaximumScalingConstructorParameters

  function velocityMaximumScalingConstructorInternal(timeScale,timeScaleMinimum,velocityExponent,redshiftExponent,cosmologyFunctions_,darkMatterProfileDMO_) result(self)
    !!{
    Default constructor for the velocityMaximumScaling hot halo outflow reincorporation class.
    !!}
    use :: Error           , only : Component_List         , Error_Report
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none
    type            (hotHaloOutflowReincorporationVelocityMaximumScaling)                        :: self
    double precision                                                     , intent(in   )         :: timeScale            , velocityExponent, &
         &                                                                                          redshiftExponent     , timeScaleMinimum
    class           (cosmologyFunctionsClass                            ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass                          ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="timeScale, velocityExponent, redshiftExponent, timeScaleMinimum, *cosmologyFunctions_, *darkMatterProfileDMO_"/>
    !!]

    ! Validate.
    if (.not.defaultHotHaloComponent%outflowedMassIsGettable())                                            &
         & call Error_Report                                                                               &
         &   (                                                                                             &
         &    'the "outflowedMass" properties of the hotHalo component must be gettable.'               // &
         &    Component_List(                                                                              &
         &                   'hotHalo'                                                                   , &
         &                   defaultHotHaloComponent%outflowedMassAttributeMatch(requireGettable=.true.)   &
         &                  )                                                                           // &
         &    {introspection:location}                                                                     &
         &   )
    ! Construct the object.
    self%timeScaleNormalization =timeScale/velocityNormalization**velocityExponent
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
  end function velocityMaximumScalingConstructorInternal

  subroutine velocityMaximumScalingAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(hotHaloOutflowReincorporationVelocityMaximumScaling), intent(inout) :: self

    call calculationResetEvent%attach(self,velocityMaximumScalingCalculationReset,openMPThreadBindingAllLevels,label='hotHaloOutflowReincorporationVelocityMaximumScaling')
    return
  end subroutine velocityMaximumScalingAutoHook

  subroutine velocityMaximumScalingDestructor(self)
    !!{
    Destructor for the \glc\ format merger tree importer class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(hotHaloOutflowReincorporationVelocityMaximumScaling), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    if (calculationResetEvent%isAttached(self,velocityMaximumScalingCalculationReset)) call calculationResetEvent%detach(self,velocityMaximumScalingCalculationReset)
    return
  end subroutine velocityMaximumScalingDestructor

  subroutine velocityMaximumScalingCalculationReset(self,node,uniqueID)
    !!{
    Reset the halo scales calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (hotHaloOutflowReincorporationVelocityMaximumScaling), intent(inout) :: self
    type   (treeNode                                           ), intent(inout) :: node
    integer(kind_int8                                          ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%velocityMaximumComputed=.false.
    self%expansionFactorComputed=.false.
    self%rateComputed           =.false.
    self%lastUniqueID           =uniqueID
    return
  end subroutine velocityMaximumScalingCalculationReset

  double precision function velocityMaximumScalingRate(self,node)
    !!{
    Return the rate of mass reincorporation for outflowed gas in the hot halo.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic   , nodeComponentHotHalo, treeNode
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (hotHaloOutflowReincorporationVelocityMaximumScaling), intent(inout) :: self
    type            (treeNode                                           ), intent(inout) :: node
    class           (nodeComponentBasic                                 ), pointer       :: basic
    class           (nodeComponentHotHalo                               ), pointer       :: hotHalo
    class           (massDistributionClass                              ), pointer       :: massDistribution_
    double precision                                                                     :: timeScale

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Get required components.
    ! Compute velocity maximum factor.
    if (.not.self%velocityMaximumComputed) then
       massDistribution_            =>       self%darkMatterProfileDMO_       %get         (node                                                                        )
       self%velocityMaximumFactor   =        self%velocityExponentiator       %exponentiate(massDistribution_%velocityRotationCurveMaximum                (            ))
       self%velocityMaximumComputed = .true.
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    end if
    ! Compute expansion factor factor.
    if (.not.self%expansionFactorComputed) then
       basic                        =>                                                      node                                          %basic          (            )
       self%expansionFactorFactor   =  1.0d0/self%expansionFactorExponentiator%exponentiate(self             %cosmologyFunctions_         %expansionFactor(basic%time()))
       self%expansionFactorComputed = .true.
    end if
    ! Compute the rate.
    if (.not.self%rateComputed) then
       hotHalo          =>      node%hotHalo               ()
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
