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
  Implements a identity output analysis distribution operator class.
  !!}

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorIdentity">
   <description>A identity output analysis distribution operator class.</description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorIdentity
     !!{
     A identity output distribution operator class.
     !!}
     private
   contains
     procedure :: operateScalar       => identityOperateScalar
     procedure :: operateDistribution => identityOperateDistribution
  end type outputAnalysisDistributionOperatorIdentity

  interface outputAnalysisDistributionOperatorIdentity
     !!{
     Constructors for the \refClass{outputAnalysisDistributionOperatorIdentity} output analysis class.
     !!}
     module procedure identityConstructorParameters
  end interface outputAnalysisDistributionOperatorIdentity

contains

  function identityConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionOperatorIdentity} output analysis distribution operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisDistributionOperatorIdentity)                :: self
    type(inputParameters                           ), intent(inout) :: parameters

    self=outputAnalysisDistributionOperatorIdentity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identityConstructorParameters

  function identityOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement a identity output analysis distribution operator.
    !!}
    use :: Arrays_Search, only : searchArray
    implicit none
    class           (outputAnalysisDistributionOperatorIdentity), intent(inout)                                        :: self
    double precision                                            , intent(in   )                                        :: propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType ), intent(in   )                                        :: propertyType
    double precision                                            , intent(in   ), dimension(:)                          :: propertyValueMinimum    , propertyValueMaximum
    integer         (c_size_t                                  ), intent(in   )                                        :: outputIndex
    type            (treeNode                                  ), intent(inout)                                        :: node
    double precision                                                           , dimension(size(propertyValueMinimum)) :: identityOperateScalar
    integer         (c_size_t                                  )                                                       :: binIndex
    !$GLC attributes unused :: self, outputIndex, propertyType, node

    ! Initialize distribution to zero.
    identityOperateScalar=0.0d0
    ! Find the corresponding bin in the array.
    binIndex=searchArray(propertyValueMinimum,propertyValue)
    ! Check if value lies within range.
    if (binIndex <= 0) return
    ! Capture the final bin.
    if (binIndex == size(propertyValueMinimum)-1_c_size_t .and. propertyValue > propertyValueMaximum(binIndex)) binIndex=size(propertyValueMinimum)
    ! Add weight to distribution if within the bin.
    if     (                                                 &
         &   propertyValue >= propertyValueMinimum(binIndex) &
         &  .and.                                            &
         &   propertyValue <  propertyValueMaximum(binIndex) &
         & ) identityOperateScalar(binIndex)=1.0d0
    return
  end function identityOperateScalar

  function identityOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement a identity output analysis distribution operator.
    !!}
    implicit none
    class           (outputAnalysisDistributionOperatorIdentity), intent(inout)                                        :: self
    double precision                                            , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType ), intent(in   )                                        :: propertyType
    double precision                                            , intent(in   ), dimension(:)                          :: propertyValueMinimum          , propertyValueMaximum
    integer         (c_size_t                                  ), intent(in   )                                        :: outputIndex
    type            (treeNode                                  ), intent(inout)                                        :: node
    double precision                                                           , dimension(size(propertyValueMinimum)) :: identityOperateDistribution
    !$GLC attributes unused :: self, propertyValueMinimum, propertyValueMaximum, outputIndex, propertyType, node

    identityOperateDistribution=distribution
    return
  end function identityOperateDistribution
