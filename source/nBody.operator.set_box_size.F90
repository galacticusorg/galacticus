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
  Implements an N-body data operator which sets the box size of the data set.
  !!}
  
  !![
  <nbodyOperator name="nbodyOperatorSetBoxSize">
   <description>An N-body data operator which sets the box size of the data set.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSetBoxSize
     !!{
     An N-body data operator which shifts sets the box size of the data set.
     !!}
     private
     double precision :: boxSize
   contains
     procedure :: operate => setBoxSizeOperate
  end type nbodyOperatorSetBoxSize

  interface nbodyOperatorSetBoxSize
     !!{
     Constructors for the {\normalfont \ttfamily setBoxSize} N-body operator class.
     !!}
     module procedure setBoxSizeConstructorParameters
     module procedure setBoxSizeConstructorInternal
  end interface nbodyOperatorSetBoxSize

contains

  function setBoxSizeConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily setBoxSize} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type           (nbodyOperatorSetBoxSize)                :: self
    type           (inputParameters        ), intent(inout) :: parameters
    double precision                                        :: boxSize
    
    !![
    <inputParameter>
      <name>boxSize</name>
      <source>parameters</source>
      <description>The box size to set.</description>
    </inputParameter>
    !!]
    self=nbodyOperatorSetBoxSize(boxSize)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function setBoxSizeConstructorParameters

  function setBoxSizeConstructorInternal(boxSize) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily setBoxSize} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorSetBoxSize)                :: self
    double precision                         , intent(in   ) :: boxSize
    !![
    <constructorAssign variables="boxSize"/>
    !!]

    return
  end function setBoxSizeConstructorInternal

  subroutine setBoxSizeOperate(self,simulations)
    !!{
    Set the box size for the data set.
    !!}
    use :: Display, only : displayIndent, displayMessage, displayUnindent, verbosityLevelStandard
    implicit none
    class  (nbodyOperatorSetBoxSize), intent(inout)               :: self
    type   (nBodyData              ), intent(inout), dimension(:) :: simulations
    integer                                                       :: i
    
    call displayIndent('setting box size',verbosityLevelStandard)
    do i=1,size(simulations)
       call simulations(i)%attributesReal%set('boxSize',self%boxSize)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine setBoxSizeOperate
