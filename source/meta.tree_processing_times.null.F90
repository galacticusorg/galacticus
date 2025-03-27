!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
Implements a merger tree processing time estimator that provides no estimates.
!!}

  !![
  <metaTreeProcessingTime name="metaTreeProcessingTimeNull">
   <description>
    A merger tree processing time class provides no estimates.
   </description>
  </metaTreeProcessingTime>
  !!]
  type, extends(metaTreeProcessingTimeClass) :: metaTreeProcessingTimeNull
     !!{
     A merger tree processing time estimator that provides no estimates.
     !!}
     private
   contains
     procedure :: time          => nullTime
     procedure :: timeRemaining => nullTimeRemaining
  end type metaTreeProcessingTimeNull

  interface metaTreeProcessingTimeNull
     !!{
     Constructors for the {\normalfont \ttfamily null} merger tree processing time estimator.
     !!}
     module procedure nullConstructorParameters
  end interface metaTreeProcessingTimeNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily null} merger tree processing time estimator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(metaTreeProcessingTimeNull)                :: self
    type(inputParameters           ), intent(inout) :: parameters

    self=metaTreeProcessingTimeNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  double precision function nullTime(self,massTree)
    !!{
    Return a null estimate of the time to process a merger tree.
    !!}
    implicit none
    class           (metaTreeProcessingTimeNull), intent(inout) :: self
    double precision                            , intent(in   ) :: massTree
    !$GLC attributes unused :: self, massTree
    
    nullTime=-1.0d0
    return
  end function nullTime

  double precision function nullTimeRemaining(self,tree,timeFinal)
    !!{
    Return a null estimate of the time to process a merger tree.
    !!}
    implicit none
    class           (metaTreeProcessingTimeNull), intent(inout) :: self
    type            (mergerTree                ), intent(inout) :: tree
    double precision                            , intent(in   ) :: timeFinal
    !$GLC attributes unused :: self, tree, timeFinal
    
    nullTimeRemaining=-1.0d0
    return
  end function nullTimeRemaining

