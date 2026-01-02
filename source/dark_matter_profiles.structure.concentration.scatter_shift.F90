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
  An implementation of dark matter halo profile concentrations that shifts the mean concentration--mass relation up/down by a specified number of $\sigma$.
  !!}

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationScatterShift">
   <description>
    A dark matter profile concentration class in which the concentration is computed by shifting the results of another
    concentration class up/down by a specified number of $\sigma$.
   </description>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationScatterShift
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{gao_redshift_2008}.
     !!}
     private
     class           (darkMatterProfileConcentrationClass), pointer :: darkMatterProfileConcentration_ => null()
     double precision                                               :: scatter                                  , sigmaShift
   contains
     final     ::                                   scatterShiftDestructor
     procedure :: concentration                  => scatterShiftConcentration
     procedure :: densityContrastDefinition      => scatterShiftDensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => scatterShiftDarkMatterProfileDefinition
  end type darkMatterProfileConcentrationScatterShift

  interface darkMatterProfileConcentrationScatterShift
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationScatterShift} dark matter halo profile concentration class.
     !!}
     module procedure scatterShiftConstructorParameters
     module procedure scatterShiftConstructorInternal
  end interface darkMatterProfileConcentrationScatterShift

contains

  function scatterShiftConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationScatterShift} dark matter halo profile concentration class which takes a parameter
    list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileConcentrationScatterShift)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (darkMatterProfileConcentrationClass       ), pointer       :: darkMatterProfileConcentration_
    double precision                                                            :: scatter                        , sigmaShift

    !![
    <objectBuilder class="darkMatterProfileConcentration" name="darkMatterProfileConcentration_" source="parameters"/>
    <inputParameter>
      <name>scatter</name>
      <source>parameters</source>
      <description>The scatter (in dex) to assume in the halo concentration distribution at fixed mass.</description>
    </inputParameter>
    <inputParameter>
      <name>sigmaShift</name>
      <source>parameters</source>
      <description>The number of $\sigma$ by which to shift the concentration.</description>
    </inputParameter>
    !!]
    self=darkMatterProfileConcentrationScatterShift(scatter,sigmaShift,darkMatterProfileConcentration_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileConcentration_"/>
    !!]
    return
  end function scatterShiftConstructorParameters

  function scatterShiftConstructorInternal(scatter,sigmaShift,darkMatterProfileConcentration_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileConcentrationScatterShift} dark matter halo profile concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Virial_Density_Contrast, only : fixedDensityTypeCritical
    implicit none
    type            (darkMatterProfileConcentrationScatterShift)                         :: self
    class           (darkMatterProfileConcentrationClass       ), intent(in   ), target  :: darkMatterProfileConcentration_
    double precision                                            , intent(in   )          :: scatter                        , sigmaShift
    !![
    <constructorAssign variables="scatter, sigmaShift, *darkMatterProfileConcentration_"/>
    !!]
 
    return
  end function scatterShiftConstructorInternal

  subroutine scatterShiftDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationScatterShift} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationScatterShift), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileConcentration_" />
    !!]
    return
  end subroutine scatterShiftDestructor

  double precision function scatterShiftConcentration(self,node) result(concentration)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} by shifting relative to another
    concentration calculation.
    !!}
    implicit none
    class(darkMatterProfileConcentrationScatterShift), intent(inout), target :: self
    type (treeNode                                  ), intent(inout), target :: node

    ! Get the mean concentration from the reference class, and apply a shift.
    concentration=+          self%darkMatterProfileConcentration_%concentrationMean(node) &
         &        *10.0d0**(                                                              &
         &                  +self                                %scatter                 &
         &                  *self                                %sigmaShift              &
         &                 )
    return
  end function scatterShiftConcentration

  function scatterShiftDensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    reference class.
    !!}
    implicit none
    class(virialDensityContrastClass                ), pointer       :: scatterShiftDensityContrastDefinition
    class(darkMatterProfileConcentrationScatterShift), intent(inout) :: self

    scatterShiftDensityContrastDefinition => self%darkMatterProfileConcentration_%densityContrastDefinition()
    return
  end function scatterShiftDensityContrastDefinition

  function scatterShiftDarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    reference class.
    !!}
    implicit none
    class(darkMatterProfileDMOClass                 ), pointer       :: scatterShiftDarkMatterProfileDefinition
    class(darkMatterProfileConcentrationScatterShift), intent(inout) :: self

    scatterShiftDarkMatterProfileDefinition => self%darkMatterProfileConcentration_%darkMatterProfileDMODefinition()
    return
  end function scatterShiftDarkMatterProfileDefinition

