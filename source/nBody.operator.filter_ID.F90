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
Implements an N-body data operator which filters particles by ID.
!!}
  
  !![
  <nbodyOperator name="nbodyOperatorFilterID">
   <description>An N-body data operator which filters particles by ID.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorFilterID
     !!{
     An N-body data operator which filters particles by ID.
     !!}
     private
     integer(c_size_t      ), allocatable, dimension(:) :: idSelection
     type   (varying_string)                            :: idSelectionFileName
   contains
     procedure :: operate => filterIDOperate
  end type nbodyOperatorFilterID

  interface nbodyOperatorFilterID
     !!{
     Constructors for the {\normalfont \ttfamily filterID} N-body operator class.
     !!}
     module procedure filterIDConstructorParameters
     module procedure filterIDConstructorInternal
  end interface nbodyOperatorFilterID

contains

  function filterIDConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily filterID} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: IO_HDF5         , only : hdf5Object
    use :: HDF5_Access     , only : hdf5Access
    implicit none
    type   (nbodyOperatorFilterID)                              :: self
    type   (inputParameters      ), intent(inout)               :: parameters
    type   (varying_string       )                              :: idSelectionFileName
    integer(c_size_t             ), allocatable  , dimension(:) :: idSelection
    type   (hdf5Object           )                              :: idFile

    if      (parameters%isPresent('idSelection'        )) then
       allocate(idSelection(parameters%count('idSelection')))
       !![
       <inputParameter>
	 <name>idSelection</name>
	 <source>parameters</source>
	 <description>The IDs of particles to retain.</description>
       </inputParameter>
       !!]
    else if (parameters%isPresent('idSelectionFileName')) then
       !![
       <inputParameter>
	 <name>idSelectionFileName</name>
	 <source>parameters</source>
	 <description>The name of a file containing the IDs of particles to retain.</description>
       </inputParameter>
       !!]
       !$ call hdf5Access%set()
       idFile=hdf5Object(char(idSelectionFileName))
       call idFile%readDataset('id',idSelection)
       !$ call hdf5Access%unset()
    else
       call Error_Report('either "idSelection" of "idSelectionFileName" must be provided'//{introspection:location})
    end if
    self=nbodyOperatorFilterID(idSelection)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    if (parameters%isPresent('idSelectionFileName')) self%idSelectionFileName=idSelectionFileName
    return
  end function filterIDConstructorParameters

  function filterIDConstructorInternal(idSelection) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily filterID} N-body operator class.
    !!}
    use :: Sorting, only : sort
    implicit none
    type   (nbodyOperatorFilterID)                              :: self
    integer(c_size_t             ), intent(in   ), dimension(:) :: idSelection
    !![
    <constructorAssign variables="idSelection"/>
    !!]

    self%idSelectionFileName=''
    call sort(self%idSelection)
    return
  end function filterIDConstructorInternal

  subroutine filterIDOperate(self,simulations)
    !!{
    Filter particles outside of a cuboid region.
    !!}
    use :: Arrays_Search, only : searchArray
    use :: Display      , only : displayIndent, displayMessage  , displayUnindent, verbosityLevelStandard
    implicit none
    class           (nbodyOperatorFilterID), intent(inout)                 :: self
    type            (nBodyData            ), intent(inout), dimension(  :) :: simulations
    logical                                , allocatable  , dimension(  :) :: mask
    integer         (c_size_t             ), pointer      , dimension(  :) :: propertyInteger     , propertyIntegerFiltered
    double precision                       , pointer      , dimension(  :) :: propertyReal        , propertyRealFiltered
    integer         (c_size_t             ), pointer      , dimension(:,:) :: propertyIntegerRank1, propertyIntegerRank1Filtered
    double precision                       , pointer      , dimension(:,:) :: propertyRealRank1   , propertyRealRank1Filtered
    integer                                                                :: i                   , j                           , &
         &                                                                    k
    integer         (c_size_t             )                                :: countFiltered       , idIndex
    
    call displayIndent('filter on particle ID',verbosityLevelStandard)
    do i=1,size(simulations)
       propertyInteger => simulations(i)%propertiesInteger%value('particleID')
       allocate(mask(size(propertyInteger)))
       !$omp parallel do private(idIndex)
       do j=1,size(mask)
          idIndex=searchArray(self%idSelection,propertyInteger(j))
          mask(j)=idIndex > 0_c_size_t .and. self%idSelection(idIndex) == propertyInteger(j)
       end do
       !$omp end parallel do
       nullify(propertyInteger)
       countFiltered=count(mask)
       ! Filter properties.
       !! Integer properties.
       do j=1,simulations(i)%propertiesInteger    %size()
          propertyInteger      => simulations(i)%propertiesInteger    %value(j)
          allocate(propertyIntegerFiltered    (                                  countFiltered))
          propertyIntegerFiltered             =pack(propertyInteger          ,mask)
          call simulations(i)%propertiesInteger    %set(simulations(i)%propertiesInteger        %key(j),propertyIntegerFiltered    )
          deallocate(propertyInteger            )
          nullify   (propertyIntegerFiltered    )
       end do
       !! Real properties.
       do j=1,simulations(i)%propertiesReal   %size()
          propertyReal         => simulations(i)%propertiesReal       %value(j)
          allocate(propertyRealFiltered   (                                      countFiltered))
          propertyRealFiltered                =pack(propertyReal             ,mask)
          call simulations(i)%propertiesReal       %set(simulations(i)%propertiesReal           %key(j),propertyRealFiltered       )
          deallocate(propertyReal               )
          nullify   (propertyRealFiltered       )
       end do
       !! Integer rank-1 properties.
       do j=1,simulations(i)%propertiesIntegerRank1%size()
          propertyIntegerRank1 => simulations(i)%propertiesIntegerRank1%value(j)
          allocate(propertyIntegerRank1Filtered(size(propertyIntegerRank1,dim=1),countFiltered))
          do k=1,size(propertyIntegerRank1,dim=1)
             propertyIntegerRank1Filtered(k,:)=pack(propertyIntegerRank1(k,:),mask)
          end do
          call simulations(i)%propertiesIntegerRank1%set(simulations(i)%propertiesIntegerRank1%key(j),propertyIntegerRank1Filtered)
          deallocate(propertyIntegerRank1        )
          nullify   (propertyIntegerRank1Filtered)
       end do
       !! Real rank-1 properties.
       do j=1,simulations(i)%propertiesRealRank1   %size()
          propertyRealRank1    => simulations(i)%propertiesRealRank1   %value(j)
          allocate(propertyRealRank1Filtered   (size(propertyRealRank1   ,dim=1),countFiltered))
          do k=1,size(propertyRealRank1   ,dim=1)
             propertyRealRank1Filtered  (k,:)=  pack(propertyRealRank1   (k,:),mask)
          end do
          call simulations(i)%propertiesRealRank1   %set(simulations(i)%propertiesRealRank1   %key(j),propertyRealRank1Filtered   )
          deallocate(propertyRealRank1           )
          nullify   (propertyRealRank1Filtered   )
       end do
       deallocate(mask)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine filterIDOperate
