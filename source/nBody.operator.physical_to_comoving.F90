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
Implements an N-body data operator which converts physical to comoving coordinates.
!!}

  !![
  <nbodyOperator name="nbodyOperatorPhysicalToComoving">
   <description>An N-body data operator which converts physical to comoving coordinates.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorPhysicalToComoving
     !!{
     An N-body data operator which converts physical to comoving coordinates.
     !!}
     private
   contains
     procedure :: operate => physicalToComovingOperate
  end type nbodyOperatorPhysicalToComoving

  interface nbodyOperatorPhysicalToComoving
     !!{
     Constructors for the {\normalfont \ttfamily physicalToComoving} N-body operator class.
     !!}
     module procedure physicalToComovingConstructorParameters
  end interface nbodyOperatorPhysicalToComoving

contains

  function physicalToComovingConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily physicalToComoving} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nbodyOperatorPhysicalToComoving)                :: self
    type(inputParameters                ), intent(inout) :: parameters

 
    self=nbodyOperatorPhysicalToComoving()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function physicalToComovingConstructorParameters

  subroutine physicalToComovingOperate(self,simulations)
    !!{
    Convert positions from physical to comoving coordinates of N-body particles.
    !!}
    implicit none
    class           (nbodyOperatorPhysicalToComoving), intent(inout)                 :: self
    type            (nBodyData                      ), intent(inout), dimension(:  ) :: simulations
    double precision                                 , pointer      , dimension(:,:) :: position
    double precision                                 , pointer      , dimension(:  ) :: expansionFactor
    integer                                                                          :: iSimulation    , i
    
    do iSimulation=1,size(simulations)
       expansionFactor => simulations(iSimulation)%propertiesReal     %value('expansionFactor')
       position        => simulations(iSimulation)%propertiesRealRank1%value('position'       )
       !$omp parallel workshare
       forall(i=1:3)
          position(i,:)=+position       (i,:) &
               &        /expansionFactor
       end forall
       !$omp end parallel workshare
       nullify(expansionFactor)
       nullify(position       )
    end do
    return
  end subroutine physicalToComovingOperate

