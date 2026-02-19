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
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRandomError" abstract="yes">
   <description>
    A random error output analysis distribution operator class. The weight of each galaxy is integrated over every bin of the
    histogram using a Gaussian kernel. This is an abstract class---the width of the Gaussian kernel must be provided by a
    concrete class.
   </description>
  </outputAnalysisDistributionOperator>
  !!]
  type, abstract, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorRandomError
     !!{
     A random error output distribution operator class.
     !!}
     private
   contains
     !![
     <methods>
       <method description="Return the root-variance to apply to the distribution." method="rootVariance" />
     </methods>
     !!]
     procedure                                           :: operateScalar       => randomErrorOperateScalar
     procedure                                           :: operateDistribution => randomErrorOperateDistribution
     procedure(randomErrorOperateRootVariance), deferred :: rootVariance
  end type outputAnalysisDistributionOperatorRandomError

  abstract interface
     double precision function randomErrorOperateRootVariance(self,propertyValue,node)
       !!{
       Abstract interface for the root variance method of random error output analysis distribution operators.
       !!}
       import outputAnalysisDistributionOperatorRandomError, treeNode
       class           (outputAnalysisDistributionOperatorRandomError), intent(inout) :: self
       double precision                                               , intent(in   ) :: propertyValue
       type            (treeNode                                     ), intent(inout) :: node
     end function randomErrorOperateRootVariance
  end interface

contains

  function randomErrorOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement a random error output analysis distribution operator.
    !!}
    implicit none
    class           (outputAnalysisDistributionOperatorRandomError), intent(inout)                                        :: self
    double precision                                               , intent(in   )                                        :: propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   )                                        :: propertyType
    double precision                                               , intent(in   ), dimension(:)                          :: propertyValueMinimum    , propertyValueMaximum
    integer         (c_size_t                                     ), intent(in   )                                        :: outputIndex
    type            (treeNode                                     ), intent(inout)                                        :: node
    double precision                                                              , dimension(size(propertyValueMinimum)) :: randomErrorOperateScalar
    double precision                                                                                                      :: rootVariance
    !$GLC attributes unused :: outputIndex, propertyType

    rootVariance=self%rootVariance(propertyValue,node)
    if     (                               &
         &   propertyValue == +huge(0.0d0) &
         &  .or.                           &
         &   propertyValue == -huge(0.0d0) &
         & ) then
       randomErrorOperateScalar=+0.0d0
    else
       randomErrorOperateScalar=+0.5d0                                                                &
            &                   *(                                                                    &
            &                     +erf((propertyValueMaximum-propertyValue)/rootVariance/sqrt(2.0d0)) &
            &                     -erf((propertyValueMinimum-propertyValue)/rootVariance/sqrt(2.0d0)) &
            &                    )
    end if
    return
  end function randomErrorOperateScalar

  function randomErrorOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement a random error output analysis distribution operator.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (outputAnalysisDistributionOperatorRandomError), intent(inout)                                        :: self
    double precision                                               , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   )                                        :: propertyType
    double precision                                               , intent(in   ), dimension(:)                          :: propertyValueMinimum          , propertyValueMaximum
    integer         (c_size_t                                     ), intent(in   )                                        :: outputIndex
    type            (treeNode                                     ), intent(inout)                                        :: node
    double precision                                                              , dimension(size(propertyValueMinimum)) :: randomErrorOperateDistribution
    !$GLC attributes unused :: self, distribution, propertyValueMinimum, propertyValueMaximum, outputIndex, propertyType, node

    randomErrorOperateDistribution=0.0d0
    call Error_Report('not implemented'//{introspection:location})
    return
  end function randomErrorOperateDistribution
