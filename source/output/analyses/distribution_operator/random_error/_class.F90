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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Implements a random error output analysis distribution operator class.
!!}
  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRandomError" abstract="yes" docformat="rst">
   <description>
   A random error output analysis distribution operator class. The weight of each galaxy is integrated over every bin of the histogram using a Gaussian kernel. This is an abstract class---the width of the Gaussian kernel must be provided by a concrete class.
   </description>
  </outputAnalysisDistributionOperator>
  !!]
  type, abstract, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorRandomError
     !!{RST
     A random error output distribution operator class.
     !!}
     private
     type(outputAnalysisDistributionOperatorIdentity) :: identity=outputAnalysisDistributionOperatorIdentity()
   contains
     !![
     <methods docformat="rst">
       <method description="Return the root-variance to apply to the distribution." method="rootVariance" />
     </methods>
     !!]
     procedure                                           :: operateScalar       => randomErrorOperateScalar
     procedure                                           :: operateDistribution => randomErrorOperateDistribution
     procedure(randomErrorOperateRootVariance), deferred :: rootVariance
  end type outputAnalysisDistributionOperatorRandomError

  abstract interface
     double precision function randomErrorOperateRootVariance(self,propertyValue,node)
       !!{RST
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
    !!{RST
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
    if (rootVariance > 0.0d0) then
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
    else
       ! Zero variance - treat as an identity operator.
       randomErrorOperateScalar=self%identity%operateScalar(propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    end if
    return
  end function randomErrorOperateScalar

  function randomErrorOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{RST
    Implement a random error output analysis distribution operator acting on a distribution. As for
    :galacticus-class:`outputAnalysisDistributionOperatorGravitationalLensing`, the weight in each bin of the input distribution
    is taken to be spread uniformly across that bin. It is then convolved with a normal distribution whose root variance is that
    at the center of the bin, which gives, for a source bin :math:`[a,b]` and target bin :math:`[A,B]`, a fraction

    .. math::

       f = \frac{\sigma}{b-a} \left[ G\left(\frac{B-a}{\sigma}\right) - G\left(\frac{B-b}{\sigma}\right) - G\left(\frac{A-a}{\sigma}\right) + G\left(\frac{A-b}{\sigma}\right) \right],

    where :math:`G(u) = u \Phi(u) + \phi(u)` is the integral of the cumulative normal distribution, :math:`\Phi(u)`, and
    :math:`\phi(u)` is the normal distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (outputAnalysisDistributionOperatorRandomError), intent(inout)                                        :: self
    double precision                                               , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   )                                        :: propertyType
    double precision                                               , intent(in   ), dimension(:)                          :: propertyValueMinimum          , propertyValueMaximum
    integer         (c_size_t                                     ), intent(in   )                                        :: outputIndex
    type            (treeNode                                     ), intent(inout)                                        :: node
    double precision                                                              , dimension(size(propertyValueMinimum)) :: randomErrorOperateDistribution
    double precision                                                                                                      :: rootVariance                  , widthBin
    integer                                                                                                               :: i                             , j
    !$GLC attributes unused :: propertyType, outputIndex

    randomErrorOperateDistribution=0.0d0
    do j=1,size(distribution)
       if (distribution(j) == 0.0d0) cycle
       rootVariance=self%rootVariance(0.5d0*(propertyValueMinimum(j)+propertyValueMaximum(j)),node)
       widthBin    =propertyValueMaximum(j)-propertyValueMinimum(j)
       if (rootVariance <= 0.0d0 .or. widthBin <= 0.0d0) then
          ! Zero variance (or a degenerate bin) - the weight stays in its own bin.
          randomErrorOperateDistribution(j)=randomErrorOperateDistribution(j)+distribution(j)
          cycle
       end if
       do i=1,size(distribution)
          randomErrorOperateDistribution(i)=+randomErrorOperateDistribution(i)                                                          &
               &                            +distribution(j)                                                                            &
               &                            *rootVariance                                                                               &
               &                            /widthBin                                                                                   &
               &                            *(                                                                                          &
               &                              +normalCumulativeIntegral((propertyValueMaximum(i)-propertyValueMinimum(j))/rootVariance) &
               &                              -normalCumulativeIntegral((propertyValueMaximum(i)-propertyValueMaximum(j))/rootVariance) &
               &                              -normalCumulativeIntegral((propertyValueMinimum(i)-propertyValueMinimum(j))/rootVariance) &
               &                              +normalCumulativeIntegral((propertyValueMinimum(i)-propertyValueMaximum(j))/rootVariance) &
               &                             )
       end do
    end do
    return

  contains

    double precision function normalCumulativeIntegral(u)
      !!{RST
      The integral of the cumulative normal distribution, :math:`G(u) = u \Phi(u) + \phi(u)`.
      !!}
      implicit none
      double precision, intent(in   ) :: u

      normalCumulativeIntegral=+u                       &
           &                   *0.5d0                   &
           &                   *erfc(-u   /sqrt(2.0d0)) &
           &                   +exp (-u**2/     2.0d0)  &
           &                   /sqrt(2.0d0*Pi)
      return
    end function normalCumulativeIntegral

  end function randomErrorOperateDistribution
