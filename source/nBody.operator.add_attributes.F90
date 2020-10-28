!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements an N-body data operator which adds attributes to the data.

  !# <nbodyOperator name="nbodyOperatorAddAttributes">
  !#  <description>An N-body data operator which adds attributes to the data.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorAddAttributes
     !% An N-body data operator which adds attributes.
     private
     double precision                , allocatable, dimension(:) :: values
     type            (varying_string), allocatable, dimension(:) :: names
   contains
     procedure :: operate => addAttributesOperate
  end type nbodyOperatorAddAttributes

  interface nbodyOperatorAddAttributes
     !% Constructors for the ``addAttributes'' N-body operator class.
     module procedure addAttributesConstructorParameters
     module procedure addAttributesConstructorInternal
  end interface nbodyOperatorAddAttributes

contains

  function addAttributesConstructorParameters(parameters) result (self)
    !% Constructor for the ``addAttributes'' N-body operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorAddAttributes)                              :: self
    type            (inputParameters           ), intent(inout)               :: parameters
    double precision                            , allocatable  , dimension(:) :: values
    type            (varying_string            ), allocatable  , dimension(:) :: names

    allocate(names (parameters%count('names' )))
    allocate(values(parameters%count('values')))
    !# <inputParameter>
    !#   <name>names</name>
    !#   <source>parameters</source>
    !#   <description>A list of attribute names.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>values</name>
    !#   <source>parameters</source>
    !#   <description>A list of attribute values.</description>
    !# </inputParameter>
    self=nbodyOperatorAddAttributes(names,values)
    !# <inputParametersValidate source="parameters"/>
    return
  end function addAttributesConstructorParameters

  function addAttributesConstructorInternal(names,values) result (self)
    !% Internal constructor for the ``addAttributes'' N-body operator class.
    implicit none
    type            (nbodyOperatorAddAttributes)                              :: self
    double precision                            , intent(in   ), dimension(:) :: values
    type            (varying_string            ), intent(in   ), dimension(:) :: names
    !# <constructorAssign variables="names, values"/>

    return
  end function addAttributesConstructorInternal

  subroutine addAttributesOperate(self,simulations)
    !% Add attributes to the simulations.
    use :: Galacticus_Display, only : Galacticus_Display_Indent , Galacticus_Display_Unindent, verbosityStandard
#ifdef USEMPI
    use :: MPI_Utilities     , only : mpiSelf
#endif
    implicit none
    class  (nbodyOperatorAddAttributes), intent(inout)               :: self
    type   (nBodyData                 ), intent(inout), dimension(:) :: simulations
    integer                                                          :: iSimulation , i

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call Galacticus_Display_Indent('add attributes',verbosityStandard)
#ifdef USEMPI
    end if
#endif
    do iSimulation=1_c_size_t,size(simulations)
       do i=1,size(self%names)
          call simulations(iSimulation)%attributesReal%set(self%names(i),self%values(i))
       end do
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call Galacticus_Display_Unindent('done',verbosityStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine addAttributesOperate

