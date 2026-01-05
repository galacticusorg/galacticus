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
  Implementation of a calculation of cold mode infall rates assuming infall on a dynamical timescale.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <coldModeInfallRate name="coldModeInfallRateDynamicalTime">
   <description>
    A cold mode infall rate class in which the infall rate from the cold mode is given by
    \begin{equation}
    \dot{M}_\mathrm{infall, cold mode} = f_\mathrm{infall, cold mode}{M_\mathrm{cold mode} \over \tau_\mathrm{dyn}},
    \end{equation}
    where $f_\mathrm{infall, cold mode}=${\normalfont \ttfamily [dynamicalRateFraction]}.
   </description>
  </coldModeInfallRate>
  !!]
  type, extends(coldModeInfallRateClass) :: coldModeInfallRateDynamicalTime
     !!{
     Implementation of a calculation of cold mode infall rates assuming infall on a dynamical timescale.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: dynamicalRateFraction
   contains
     final     ::               dynamicalTimeDestructor
     procedure :: infallRate => dynamicalTimeInfallRate
  end type coldModeInfallRateDynamicalTime

  interface coldModeInfallRateDynamicalTime
     !!{
     Constructors for the dynamicalTime cooling time class.
     !!}
     module procedure dynamicalTimeConstructorParameters
     module procedure dynamicalTimeConstructorInternal
  end interface coldModeInfallRateDynamicalTime

contains

  function dynamicalTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the dynamical time cooling time class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coldModeInfallRateDynamicalTime)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    double precision                                                 :: dynamicalRateFraction

    !![
    <inputParameter>
      <name>dynamicalRateFraction</name>
      <defaultValue>2.0d0</defaultValue>
      <source>parameters</source>
      <description>The fraction of the inverse dynamical time to use as the rate for infall of the cold mode component.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=coldModeInfallRateDynamicalTime(dynamicalRateFraction,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function dynamicalTimeConstructorParameters

  function dynamicalTimeConstructorInternal(dynamicalRateFraction,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the dynamical time cooling time class.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none
    type            (coldModeInfallRateDynamicalTime)                        :: self
    class           (darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                 , intent(in   )         :: dynamicalRateFraction
    !![
    <constructorAssign variables="dynamicalRateFraction, *darkMatterHaloScale_"/>
    !!]

    ! Check that the properties we need are gettable.
    if (.not.defaultHotHaloComponent%massColdIsGettable())                       &
         & call Error_Report(                                                    &
         &                   'hot halo component must have gettable cold mass'// &
         &                   {introspection:location}                            &
         &                  )
    return
  end function dynamicalTimeConstructorInternal

  subroutine dynamicalTimeDestructor(self)
    !!{
    Destructor for the dynamical time cold mode infall rate class.
    !!}
    implicit none
    type(coldModeInfallRateDynamicalTime), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine dynamicalTimeDestructor

  double precision function dynamicalTimeInfallRate(self,node)
    !!{
    Computes the cold mode infall rate as a fraction of the halo dynamical time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class(coldModeInfallRateDynamicalTime), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node
    class(nodeComponentHotHalo           ), pointer       :: hotHalo

    hotHalo                 =>  node                        %hotHalo              (    )
    dynamicalTimeInfallRate =  +self                        %dynamicalRateFraction       &
         &                     /self   %darkMatterHaloScale_%timescaleDynamical   (node) &
         &                     *hotHalo                     %massCold             (    )
    return
  end function dynamicalTimeInfallRate
