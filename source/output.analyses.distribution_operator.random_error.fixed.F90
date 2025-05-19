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
Implements a random error output analysis distribution operator class.
!!}
  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRandomErrorFixed">
   <description>A random error output analysis distribution operator class.</description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorRandomError) :: outputAnalysisDistributionOperatorRandomErrorFixed
     !!{
     A random error output distribution operator class which has a fixed error magnitude.
     !!}
     private
     double precision :: rootVariance_
   contains
     procedure :: rootVariance => randomErrorFixedRootVariance
  end type outputAnalysisDistributionOperatorRandomErrorFixed

  interface outputAnalysisDistributionOperatorRandomErrorFixed
     !!{
     Constructors for the {\normalfont \ttfamily randomErrorFixed} output analysis distribution operator class.
     !!}
     module procedure randomErrorFixedConstructorParameters
     module procedure randomErrorFixedConstructorInternal
  end interface outputAnalysisDistributionOperatorRandomErrorFixed

contains

  function randomErrorFixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily randomErrorFixed} output analysis distribution operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisDistributionOperatorRandomErrorFixed)                :: self
    type            (inputParameters                                   ), intent(inout) :: parameters
    double precision                                                                    :: rootVariance_

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>rootVariance</name>
      <source>parameters</source>
      <variable>rootVariance_</variable>
      <description>The root variance of the random error distribution.</description>
    </inputParameter>
    !!]
    ! Construct the object.
    self=outputAnalysisDistributionOperatorRandomErrorFixed(rootVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function randomErrorFixedConstructorParameters

  function randomErrorFixedConstructorInternal(rootVariance_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily randomErrorFixed} output analysis distribution operator class.
    !!}
    implicit none
    type            (outputAnalysisDistributionOperatorRandomErrorFixed)                :: self
    double precision                                                    , intent(in   ) :: rootVariance_
    !![
    <constructorAssign variables="rootVariance_"/>
    !!]

    return
  end function randomErrorFixedConstructorInternal

  double precision function randomErrorFixedRootVariance(self,propertyValue,node)
    !!{
    Return the root-variance in the fixed random error distribution operator.
    !!}
    implicit none
    class           (outputAnalysisDistributionOperatorRandomErrorFixed), intent(inout) :: self
    double precision                                                    , intent(in   ) :: propertyValue
    type            (treeNode                                          ), intent(inout) :: node
    !$GLC attributes unused :: propertyValue, node

    randomErrorFixedRootVariance=self%rootVariance_
    return
  end function randomErrorFixedRootVariance
