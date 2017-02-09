!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a random error output analysis distribution operator class. 
  !# <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRandomError">
  !#  <description>A random error output analysis distribution operator class.</description>
  !# </outputAnalysisDistributionOperator>
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorRandomError
     !% A random error output distribution operator class.
     private
     double precision :: rootVariance
   contains
     procedure :: operateScalar       => randomErrorOperateScalar
     procedure :: operateDistribution => randomErrorOperateDistribution
  end type outputAnalysisDistributionOperatorRandomError

  interface outputAnalysisDistributionOperatorRandomError
     !% Constructors for the ``randomError'' output analysis distribution operator class.
     module procedure randomErrorConstructorParameters
     module procedure randomErrorConstructorInternal
  end interface outputAnalysisDistributionOperatorRandomError

contains

  function randomErrorConstructorParameters(parameters) result(self)
    !% Constructor for the ``randomError'' output analysis distribution operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (outputAnalysisDistributionOperatorRandomError)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    double precision                                                               :: rootVariance
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>rootVariance</name>
    !#   <source>parameters</source>
    !#   <variable>rootVariance</variable>
    !#   <description>The root variance of the random error distribution.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    ! Construct the object.
    self=outputAnalysisDistributionOperatorRandomError(rootVariance)
    return
  end function randomErrorConstructorParameters

  function randomErrorConstructorInternal(rootVariance) result(self)
    !% Constructor for the ``randomError'' output analysis distribution operator class which takes a parameter set as input.
    implicit none
    type            (outputAnalysisDistributionOperatorRandomError)                :: self
    double precision                                               , intent(in   ) :: rootVariance
    !# <constructorAssign variables="rootVariance"/>

    return
  end function randomErrorConstructorInternal

  function randomErrorOperateScalar(self,propertyValue,propertyValueMinimum,propertyValueMaximum)
    !% Implement a random error output analysis distribution operator.
    implicit none
    class           (outputAnalysisDistributionOperatorRandomError), intent(inout)                                        :: self
    double precision                                               , intent(in   )                                        :: propertyValue
    double precision                                               , intent(in   ), dimension(:)                          :: propertyValueMinimum    , propertyValueMaximum
    double precision                                                              , dimension(size(propertyValueMinimum)) :: randomErrorOperateScalar
    
    randomErrorOperateScalar=+0.5d0                                                                     &
         &                   *(                                                                         &
         &                     +erf((propertyValueMaximum-propertyValue)/self%rootVariance/sqrt(2.0d0)) &
         &                     -erf((propertyValueMinimum-propertyValue)/self%rootVariance/sqrt(2.0d0)) &
         &                    )
    return
  end function randomErrorOperateScalar

  function randomErrorOperateDistribution(self,distribution,propertyValueMinimum,propertyValueMaximum)
    !% Implement a random error output analysis distribution operator.
    use Galacticus_Error
    implicit none
    class           (outputAnalysisDistributionOperatorRandomError), intent(inout)                                        :: self
    double precision                                               , intent(in   ), dimension(:)                          :: distribution
    double precision                                               , intent(in   ), dimension(:)                          :: propertyValueMinimum          , propertyValueMaximum
    double precision                                                              , dimension(size(propertyValueMinimum)) :: randomErrorOperateDistribution
    !GCC$ attributes unused :: self, distribution, propertyValueMinimum, propertyValueMaximum

    randomErrorOperateDistribution=0.0d0
    call Galacticus_Error_Report('randomErrorOperateDistribution','not implemented')
    return
  end function randomErrorOperateDistribution
