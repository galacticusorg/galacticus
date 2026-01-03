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
  An implementation of dark matter halo profile scale radii in which radii from another class are limited to enforce bounds on concentration.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusConcentrationLimiter">
   <description>Dark matter halo scale radii from another class are limited to enforce bounds on concentration.</description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusConcentrationLimiter
     !!{
     A dark matter halo profile scale radius class in which radii from another class are limited to enforce bounds on concentration.
     !!}
     private
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_          => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     double precision                                             :: concentrationMinimum                   , concentrationMaximum
   contains
     final     ::           concentrationLimiterDestructor
     procedure :: radius => concentrationLimiterRadius
  end type darkMatterProfileScaleRadiusConcentrationLimiter

  interface darkMatterProfileScaleRadiusConcentrationLimiter
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusConcentrationLimiter} dark matter halo profile scale radius class.
     !!}
     module procedure concentrationLimiterConstructorParameters
     module procedure concentrationLimiterConstructorInternal
  end interface darkMatterProfileScaleRadiusConcentrationLimiter

contains

  function concentrationLimiterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusConcentrationLimiter} dark matter halo profile scale radius class which takes a
    parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileScaleRadiusConcentrationLimiter)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                        ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass               ), pointer       :: darkMatterProfileScaleRadius_
    double precision                                                                  :: concentrationMinimum         , concentrationMaximum

    !![
    <inputParameter>
      <name>concentrationMinimum</name>
      <description>The minimum concentration to allow for halos.</description>
      <source>parameters</source>
    </inputParameter>   
    <inputParameter>
      <name>concentrationMaximum</name>
      <description>The maximum concentration to allow for halos.</description>
      <source>parameters</source>
    </inputParameter>   
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    self=darkMatterProfileScaleRadiusConcentrationLimiter(concentrationMinimum,concentrationMaximum,darkMatterHaloScale_,darkMatterProfileScaleRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function concentrationLimiterConstructorParameters

  function concentrationLimiterConstructorInternal(concentrationMinimum,concentrationMaximum,darkMatterHaloScale_,darkMatterProfileScaleRadius_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileScaleRadiusConcentrationLimiter} dark matter halo profile scale radius class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (darkMatterProfileScaleRadiusConcentrationLimiter)                        :: self
    class           (darkMatterHaloScaleClass                        ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass               ), intent(in   ), target :: darkMatterProfileScaleRadius_
    double precision                                                  , intent(in   )         :: concentrationMinimum         , concentrationMaximum
    !![
    <constructorAssign variables="concentrationMinimum, concentrationMaximum, *darkMatterHaloScale_, *darkMatterProfileScaleRadius_"/>
    !!]

    if (concentrationMaximum <= concentrationMinimum) call Error_Report('cₘₐₓ > cₘᵢₙ is required'//{introspection:location})
    return
  end function concentrationLimiterConstructorInternal

  subroutine concentrationLimiterDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileScaleRadiusConcentrationLimiter} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(darkMatterProfileScaleRadiusConcentrationLimiter), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"         />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    !!]
    return
  end subroutine concentrationLimiterDestructor

  double precision function concentrationLimiterRadius(self,node)
    !!{
    Compute the scale radius of the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class           (darkMatterProfileScaleRadiusConcentrationLimiter), intent(inout), target :: self
    type            (treeNode                                        ), intent(inout), target :: node
    double precision                                                                          :: radiusVirial

    radiusVirial              =        self%darkMatterHaloScale_         %radiusVirial(node)
    concentrationLimiterRadius=min(                                                                                     &
         &                                                                radiusVirial      /self%concentrationMinimum, &
         &                         max(                                                                                 &
         &                                                                radiusVirial      /self%concentrationMaximum, &
         &                             self%darkMatterProfileScaleRadius_%radius      (node)                            &
         &                            )                                                                                 &
         &                        )
    return
  end function concentrationLimiterRadius
