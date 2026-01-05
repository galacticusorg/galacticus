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
  Implementation of a cooling rate class in which the cooling rate scales with the mass of the halo.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  !![
  <coolingRate name="coolingRateSimpleScaling">
   <description>
    A cooling rate class in which the cooling rate scales with the mass of the halo. Specifically, the cooling rate is given by
    \begin{equation}
    \dot{M}_\mathrm{cool} = M_\mathrm{hot}/\tau_\mathrm{cool}(M_\mathrm{halo},z) ,
    \end{equation}
    where 
    \begin{equation}
    \tau_\mathrm{cool}=\tau_\mathrm{cool,0} (1+z)^{\beta_\mathrm{cool}} \exp \left(\left[{M_\mathrm{halo} \over
    M_\mathrm{transition}}\right]^{\gamma_\mathrm{cool}}\right),
    \end{equation}
    $\tau_\mathrm{cool,0}=${\normalfont \ttfamily [timescale]}, $\beta_\mathrm{cool}=${\normalfont \ttfamily
    [exponentRedshift]}, $M_\mathrm{transition}=${\normalfont \ttfamily [massCutOff]}, and $\gamma_\mathrm{cool}=${\normalfont
    \ttfamily [exponentCutOff]}.
   </description>
  </coolingRate>
  !!]
  type, extends(coolingRateClass) :: coolingRateSimpleScaling
     !!{
     Implementation of cooling rate class in which the cooling rate scales with the mass of the halo.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     ! Parameters controlling the cooling rate.
     double precision                                   :: timescale                         , exponentRedshift , &
          &                                                widthCutOff                       , massCutOff       , &
          &                                                exponentCutOff                    , normalization
     ! Stored values for rapid re-use.
     double precision                                   :: expansionFactorPrevious           , massBasicPrevious, &
          &                                                coolingRateStored
   contains
     final     ::         simpleScalingDestructor
     procedure :: rate => simpleScalingRate
  end type coolingRateSimpleScaling

  interface coolingRateSimpleScaling
     !!{
     Constructors for the simple scaling cooling rate class.
     !!}
     module procedure simpleScalingConstructorParameters
     module procedure simpleScalingConstructorInternal
  end interface coolingRateSimpleScaling

contains

  function simpleScalingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the simple scaling cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coolingRateSimpleScaling)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    double precision                                          :: timeScale          , exponentCutOff, &
         &                                                       exponentRedshift   , widthCutOff   , &
         &                                                       massCutOff

    !![
    <inputParameter>
      <name>timescale</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The timescale (in Gyr) for cooling in low mass halos at $z=0$ in the simple scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentRedshift</name>
      <source>parameters</source>
      <defaultValue>-1.5d0</defaultValue>
      <description>The exponent of $(1+z)$ in the cooling timescale for low mass halos in the simple scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>massCutOff</name>
      <source>parameters</source>
      <defaultValue>200.0d0</defaultValue>
      <description>The halo mass scale appearing in the exponential term for cooling timescale in the simple cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>widthCutOff</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The width appearing in the exponential term for cooling timescale in the simple scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentCutOff</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The exponent appearing in the exponential term for cooling timescale in the simple scaling cooling rate model.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=coolingRateSimpleScaling(timeScale,exponentRedshift,massCutOff,widthCutOff,exponentCutOff,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function simpleScalingConstructorParameters

  function simpleScalingConstructorInternal(timeScale,exponentRedshift,massCutOff,widthCutOff,exponentCutOff,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the simple scaling cooling rate class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List          , Error_Report
    use :: Galacticus_Nodes, only : defaultBasicComponent   , defaultHotHaloComponent
    implicit none
    type            (coolingRateSimpleScaling)                        :: self
    double precision                          , intent(in   )         :: timeScale          , exponentRedshift, &
         &                                                               massCutOff         , widthCutOff     , &
         &                                                               exponentCutOff
    class           (cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_

    !![
    <constructorAssign variables="timeScale, exponentRedshift, massCutOff, widthCutOff, exponentCutOff, *cosmologyFunctions_"/>
    !!]
    ! Check that the properties we need are gettable.
    if (.not.defaultHotHaloComponent%massIsGettable())                                                            &
         & call Error_Report(                                                                                     &
         &                   'Hot halo component must have gettable mass.'                                     // &
         &                   Component_List(                                                                      &
         &                                  'hotHalo'                                                          ,  &
         &                                   defaultHotHaloComponent%massAttributeMatch(requireGettable=.true.)   &
         &                                 )                                                                   // &
         &                   {introspection:location}                                                             &
         &                  )
    if     (                                                                                                      &
         &  .not.(                                                                                                &
         &         defaultBasicComponent%massIsGettable()                                                         &
         &        .and.                                                                                           &
         &         defaultBasicComponent%timeIsGettable()                                                         &
         &       )                                                                                                &
         & ) call Error_Report(                                                                                   &
         &                     'Basic component must have gettable mass and time.'//                              &
         &                     Component_List(                                                                    &
         &                                    'basic'                                                          ,  &
         &                                     defaultBasicComponent%massAttributeMatch(requireGettable=.true.)   &
         &                                    .intersection.                                                      &
         &                                     defaultBasicComponent%timeAttributeMatch(requireGettable=.true.)   &
         &                                   )                                                                 // &
         &                     {introspection:location}                                                           &
         &                    )
    ! Initialize stored solutions.
    self%expansionFactorPrevious=-1.0d0
    self%massBasicPrevious      =-1.0d0
    self%coolingRateStored      =-1.0d0
    return
  end function simpleScalingConstructorInternal

  subroutine simpleScalingDestructor(self)
    !!{
    Destructor for the simple scaling cooling rate class.
    !!}
    implicit none
    type(coolingRateSimpleScaling), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine simpleScalingDestructor

  double precision function simpleScalingRate(self,node)
    !!{
    Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for a model in which this rate scales with the mass of the halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentHotHalo, treeNode
    implicit none
    class           (coolingRateSimpleScaling), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , parameter     :: expArgumentMaximum=100.0d0
    class           (nodeComponentBasic      ), pointer       :: basic
    class           (nodeComponentHotHalo    ), pointer       :: hotHalo
    double precision                                          :: expFactor                 , expansionFactor, &
         &                                                       expArgument

    ! Compute expansion factor.
    basic               => node                    %basic          (            )
    hotHalo             => node                    %hotHalo        (            )
    expansionFactor     =  self%cosmologyFunctions_%expansionFactor(basic%time())
    if (expansionFactor /= self%expansionFactorPrevious .or. basic%mass() /= self%massBasicPrevious) then
       expArgument=+log10(                 &
            &              basic%mass()    &
            &             /self%massCutOff &
            &            )                 &
            &      /self%widthCutOff
       if (expArgument < expArgumentMaximum) then
          expFactor=1.0d0/(1.0d0+exp(+expArgument))**self%exponentCutOff
       else
          expFactor=             exp(-expArgument   *self%exponentCutOff)
       end if
       self%coolingRateStored      =+expansionFactor**self%exponentRedshift &
            &                       /self%timescale                         &
            &                       *expFactor
       self%expansionFactorPrevious= expansionFactor
       self%massBasicPrevious      = basic          %mass()
    end if
    simpleScalingRate=hotHalo%mass()*self%coolingRateStored
    return
  end function simpleScalingRate
