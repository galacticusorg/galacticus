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
  Implements calculations of incompleteness assuming a normal distribution of surface brightnesses.
  !!}

  !![
  <massFunctionIncompleteness name="massFunctionIncompletenessSurfaceBrightness">
   <description>
    A mass function incompleteness class which models the surface brightness distribution of galaxies as a normal distribution
    with mean $\langle \mu \rangle (M) = \alpha \log_{10}(M/M_0)+\beta$, with root-variance $\sigma$, where
    $\alpha=${\normalfont \ttfamily [slope]}, $\beta=${\normalfont \ttfamily [zeroPoint]}, and $\sigma=${\normalfont \ttfamily
    [scatter]}. The completeness is the fraction of this distribution above the surface brightness limit given by {\normalfont
    \ttfamily [limit]}.
   </description>
  </massFunctionIncompleteness>
  !!]
  type, extends(massFunctionIncompletenessClass) :: massFunctionIncompletenessSurfaceBrightness
     !!{
     A class implementing incompleteness calculations assuming a normal distribution of surface brightnesses.
     !!}
     private
     double precision :: limit  , zeroPoint, &
          &              slope  , offset   , &
          &              scatter
   contains
     procedure :: completeness => surfaceBrightnessCompleteness
  end type massFunctionIncompletenessSurfaceBrightness

  interface massFunctionIncompletenessSurfaceBrightness
     !!{
     Constructors for the \refClass{massFunctionIncompletenessSurfaceBrightness} incompleteness class.
     !!}
     module procedure surfaceBrightnessConstructorParameters
     module procedure surfaceBrightnessConstructorInternal
  end interface massFunctionIncompletenessSurfaceBrightness

contains

  function surfaceBrightnessConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massFunctionIncompletenessSurfaceBrightness} incompleteness class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massFunctionIncompletenessSurfaceBrightness)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    double precision                                                             :: limit     , zeroPoint, &
         &                                                                          slope     , offset   , &
         &                                                                          scatter

    !![
    <inputParameter>
      <name>limit</name>
      <description>Limiting surface brightness for mass function incompleteness calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>zeroPoint</name>
      <description>Mass zero point for the mass function incompleteness surface brightness model, i.e. $M_0$ in $\mu(M) = \alpha \log_{10}(M/M_0)+\beta$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>slope</name>
      <description>Slope of mass function incompleteness surface brightness model, i.e. $\alpha$ in $\mu(M) = \alpha \log_{10}(M/M_0)+\beta$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>offset</name>
      <description>Offset in the mass function incompleteness surface brightness model, i.e. $beta$ in $\mu(M) = \alpha \log_{10}(M/M_0)+\beta$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scatter</name>
      <description>Scatter in the mass function incompleteness surface brightness model.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=massFunctionIncompletenessSurfaceBrightness(limit,zeroPoint,slope,offset,scatter)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function surfaceBrightnessConstructorParameters

  function surfaceBrightnessConstructorInternal(limit,zeroPoint,slope,offset,scatter) result(self)
    !!{
    Internal constructor for the \refClass{massFunctionIncompletenessSurfaceBrightness} incompleteness class.
    !!}
    implicit none
    type            (massFunctionIncompletenessSurfaceBrightness)                :: self
    double precision                                             , intent(in   ) :: limit  , zeroPoint, &
         &                                                                          slope  , offset   , &
         &                                                                          scatter
    !![
    <constructorAssign variables="limit, zeroPoint, slope, offset, scatter"/>
    !!]

    return
  end function surfaceBrightnessConstructorInternal

  double precision function surfaceBrightnessCompleteness(self,mass)
    !!{
    Return the completeness.
    !!}
    use :: Error_Functions, only : Error_Function
    implicit none
    class           (massFunctionIncompletenessSurfaceBrightness), intent(inout) :: self
    double precision                                             , intent(in   ) :: mass
    double precision                                                             :: surfaceBrightnessMean, limitNormalized

    surfaceBrightnessMean        =+      self%slope                  &
         &                        *log10(                            &
         &                               +     mass                  &
         &                               /self%zeroPoint             &
         &                              )                            &
         &                        +       self%offset
    limitNormalized              =+(                                 &
         &                          +     self%limit                 &
         &                          -          surfaceBrightnessMean &
         &                         )                                 &
         &                        /       self%scatter
    surfaceBrightnessCompleteness=+0.5d0                             &
         &                        *(                                 &
         &                          +1.0d0                           &
         &                          +Error_Function(limitNormalized) &
         &                        )
    return
  end function surfaceBrightnessCompleteness
