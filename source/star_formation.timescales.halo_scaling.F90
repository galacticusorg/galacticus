!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implementation of a timescale for star formation which scales with the circular velocity of the host halo.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <starFormationTimescale name="starFormationTimescaleHaloScaling">
   <description>
    A star formation timescale class in which the timescale scales with halo properties. Specifically,
    \begin{equation}
     \tau_\star = \tau_\mathrm{\star,0} \left( {V_\mathrm{vir} \over 200\hbox{km/s}} \right)^{\alpha_\star} (1+z)^{\beta_\star},
    \end{equation}
    where $\tau_\mathrm{\star,0}=${\normalfont \ttfamily [timescale]}, $\alpha_\star=${\normalfont \ttfamily
    [exponentVelocityVirial]}, and $\beta_\star=${\normalfont \ttfamily [exponentRedshift]}.
   </description>
  </starFormationTimescale>
  !!]
  type, extends(starFormationTimescaleClass) :: starFormationTimescaleHaloScaling
     !!{
     Implementation of a haloScaling timescale for star formation.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_           => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_          => null()
     double precision                                    :: expansionFactorFactorPrevious          , exponentVelocityVirial , &
          &                                                 exponentRedshift                       , timescaleNormalization , &
          &                                                 timescaleStored                        , velocityPrevious       , &
          &                                                 velocityFactorPrevious                 , expansionFactorPrevious, &
          &                                                 timeScale_
     logical                                             :: timescaleComputed
     integer         (kind_int8                        ) :: lastUniqueID
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                     haloScalingDestructor
     procedure :: autoHook         => haloScalingAutoHook
     procedure :: timescale        => haloScalingTimescale
     procedure :: calculationReset => haloScalingCalculationReset
  end type starFormationTimescaleHaloScaling

  interface starFormationTimescaleHaloScaling
     !!{
     Constructors for the \refClass{starFormationTimescaleHaloScaling} timescale for star formation class.
     !!}
     module procedure haloScalingConstructorParameters
     module procedure haloScalingConstructorInternal
  end interface starFormationTimescaleHaloScaling

  double precision, parameter :: velocityVirialNormalization=200.0d0

contains

  function haloScalingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationTimescaleHaloScaling} timescale for star formation feedback class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationTimescaleHaloScaling)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
    double precision                                                   :: timescale           , exponentVelocityVirial, &
         &                                                                exponentRedshift

    !![
    <inputParameter>
      <name>timescale</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The timescale for star formation in the halo scaling timescale model.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponentVelocityVirial</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent of virial velocity in the timescale for star formation in the halo scaling timescale model.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponentRedshift</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent of redshift in the timescale for star formation in the halo scaling timescale model.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=starFormationTimescaleHaloScaling(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function haloScalingConstructorParameters

  function haloScalingConstructorInternal(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{starFormationTimescaleHaloScaling} timescale for star formation class.
    !!}
    implicit none
    type            (starFormationTimescaleHaloScaling)                        :: self
    double precision                                   , intent(in   )         :: timescale           , exponentVelocityVirial, &
         &                                                                        exponentRedshift
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="exponentVelocityVirial, exponentRedshift, *cosmologyFunctions_, *darkMatterHaloScale_"/>
    !!]

    self%lastUniqueID                 =-1_kind_int8
    self%timescaleComputed            =.false.
    self%velocityPrevious             =-1.0d0
    self%velocityFactorPrevious       =-1.0d0
    self%expansionFactorPrevious      =-1.0d0
    self%expansionFactorFactorPrevious=-1.0d0
    ! Compute the normalization of the timescale.
    self%timeScale_            =+timescale
    self%timeScaleNormalization=+timescale                                                &
         &                      /velocityVirialNormalization**self%exponentVelocityVirial
    return
  end function haloScalingConstructorInternal

  subroutine haloScalingAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(starFormationTimescaleHaloScaling), intent(inout) :: self

    call calculationResetEvent%attach(self,haloScalingCalculationReset,openMPThreadBindingAllLevels,label='starFormationTimescaleHaloScaling')
    return
  end subroutine haloScalingAutoHook

  subroutine haloScalingDestructor(self)
    !!{
    Destructor for the \refClass{starFormationTimescaleHaloScaling} timescale for star formation class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(starFormationTimescaleHaloScaling), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    if (calculationResetEvent%isAttached(self,haloScalingCalculationReset)) call calculationResetEvent%detach(self,haloScalingCalculationReset)
    return
  end subroutine haloScalingDestructor

  subroutine haloScalingCalculationReset(self,node,uniqueID)
    !!{
    Reset the halo scaling star formation timescale calculation.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    use :: Kind_Numbers    , only : kind_int8
    implicit none
    class  (starFormationTimescaleHaloScaling), intent(inout) :: self
    type   (treeNode                         ), intent(inout) :: node
    integer(kind_int8                        ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%timescaleComputed=.false.
    self%lastUniqueID     =uniqueID
    return
  end subroutine haloScalingCalculationReset

  double precision function haloScalingTimescale(self,component) result(timescale)
    !!{
    Returns the timescale (in Gyr) for star formation in the given {\normalfont \ttfamily component} in the halo scaling
    timescale model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (starFormationTimescaleHaloScaling), intent(inout) :: self
    class           (nodeComponent                    ), intent(inout) :: component
    class           (nodeComponentBasic               ), pointer       :: basic
    double precision                                                   :: expansionFactor, velocityVirial

    ! Check if node differs from previous one for which we performed calculations.
    if (component%hostNode%uniqueID() /= self%lastUniqueID) call self%calculationReset(component%hostNode,component%hostNode%uniqueID())
    ! Compute the timescale if necessary.
    if (.not.self%timescaleComputed) then
       ! Get virial velocity and expansion factor.
       basic           => component%hostNode%basic                               (                    )
       velocityVirial  =  self              %darkMatterHaloScale_%velocityVirial (component%hostNode  )
       expansionFactor =  self              %cosmologyFunctions_ %expansionFactor(basic    %time    ())
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
    timescale=self%timescaleStored
    return
  end function haloScalingTimescale
