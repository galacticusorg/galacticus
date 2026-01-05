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
Implements an N-body data operator which subsamples points at a given rate.
!!}

  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
 
  !![
  <nbodyOperator name="nbodyOperatorSubsample">
   <description>An N-body data operator which filters out particles based on a property range.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSubsample
     !!{
     An N-body data operator which filters out particles based on a property range.
     !!}
     private
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
     double precision                                      :: rate
  contains
    final :: subsampleDestructor
     procedure :: operate => subsampleOperate
  end type nbodyOperatorSubsample

  interface nbodyOperatorSubsample
     !!{
     Constructors for the \refClass{nbodyOperatorSubsample} N-body operator class.
     !!}
     module procedure subsampleConstructorParameters
     module procedure subsampleConstructorInternal
  end interface nbodyOperatorSubsample

contains
  
  function subsampleConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorSubsample} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nbodyOperatorSubsample    )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    double precision                                            :: rate

    !![
    <inputParameter>
      <name>rate</name>
      <source>parameters</source>
      <description>The rate at which to subsample points.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorSubsample(rate,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function subsampleConstructorParameters

  function subsampleConstructorInternal(rate,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorSubsample} N-body operator class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (nbodyOperatorSubsample    )                        :: self
    class           (randomNumberGeneratorClass), intent(in   ), target :: randomNumberGenerator_
    double precision                            , intent(in   )         :: rate
    !![
    <constructorAssign variables="rate, *randomNumberGenerator_"/>
    !!]

    if (rate <= 0.0d0 .or. rate > 1.0d0) call Error_Report('range must be in the interval [0,1)'//{introspection:location})
    return
  end function subsampleConstructorInternal

  subroutine subsampleDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorSubsample} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorSubsample), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine subsampleDestructor

  subroutine subsampleOperate(self,simulations)
    !!{
    Identify and flag particles which have been always isolated.
    !!}
    use :: Display              , only : displayIndent      , displayMessage  , displayUnindent, verbosityLevelStandard
    use :: Error                , only : Error_Report
    use :: NBody_Simulation_Data, only : propertyTypeInteger, propertyTypeReal
    implicit none
    class           (nbodyOperatorSubsample), intent(inout)                 :: self
    type            (nBodyData             ), intent(inout), dimension(  :) :: simulations
    logical                                 , allocatable  , dimension(  :) :: mask
    integer         (c_size_t              ), pointer      , dimension(  :) :: propertyInteger     , propertyIntegerFiltered
    double precision                        , pointer      , dimension(  :) :: propertyReal        , propertyRealFiltered
    integer         (c_size_t              ), pointer      , dimension(:,:) :: propertyIntegerRank1, propertyIntegerRank1Filtered
    double precision                        , pointer      , dimension(:,:) :: propertyRealRank1   , propertyRealRank1Filtered
    integer                                                                 :: i                   , j
    integer         (c_size_t              )                                :: k                   , countSubsampled             , &
         &                                                                     countOriginal
    
    call displayIndent('subsample points',verbosityLevelStandard)
    do i=1,size(simulations)
       countOriginal=-1_c_size_t
       if (countOriginal < 0_c_size_t .and. simulations(i)%propertiesInteger     %size() > 0) then
          propertyInteger      => simulations(i)%propertiesInteger     %value(1)
          countOriginal        =  size(propertyInteger     ,dim=1)
       end if
       if (countOriginal < 0_c_size_t .and. simulations(i)%propertiesReal        %size() > 0) then
          propertyReal         => simulations(i)%propertiesReal        %value(1)
          countOriginal        =  size(propertyReal        ,dim=1)
       end if
       if (countOriginal < 0_c_size_t .and. simulations(i)%propertiesIntegerRank1%size() > 0) then
          propertyIntegerRank1 => simulations(i)%propertiesIntegerRank1%value(1)
          countOriginal        =  size(propertyIntegerRank1,dim=2)
       end if
       if (countOriginal < 0_c_size_t .and. simulations(i)%propertiesRealRank1   %size() > 0) then
          propertyRealRank1    => simulations(i)%propertiesRealRank1   %value(1)
          countOriginal        =  size(propertyRealRank1   ,dim=2)
       end if
       if (countOriginal <= 0_c_size_t) call Error_Report('no properties available to subsample'//{introspection:location})
       allocate(mask(countOriginal))
       do k=1,size(mask)
          mask(k)=self%randomNumberGenerator_%uniformSample() < self%rate
       end do
       countSubsampled=count(mask)
       ! Subsample default properties.
       !! Integer properties.
       do j=1,simulations(i)%propertiesInteger    %size()
          propertyInteger      => simulations(i)%propertiesInteger    %value(j)
          allocate(propertyIntegerFiltered    (                                  countSubsampled))
          propertyIntegerFiltered             =pack(propertyInteger          ,mask)
          call simulations(i)%propertiesInteger    %set(simulations(i)%propertiesInteger        %key(j),propertyIntegerFiltered    )
          deallocate(propertyInteger            )
          nullify   (propertyIntegerFiltered    )
       end do
       !! Real properties.
       do j=1,simulations(i)%propertiesReal   %size()
          propertyReal         => simulations(i)%propertiesReal       %value(j)
          allocate(propertyRealFiltered   (                                      countSubsampled))
          propertyRealFiltered                =pack(propertyReal             ,mask)
          call simulations(i)%propertiesReal       %set(simulations(i)%propertiesReal           %key(j),propertyRealFiltered       )
          deallocate(propertyReal               )
          nullify   (propertyRealFiltered       )
       end do
       !! Integer rank-1 properties.
       do j=1,simulations(i)%propertiesIntegerRank1%size()
          propertyIntegerRank1 => simulations(i)%propertiesIntegerRank1%value(j)
          allocate(propertyIntegerRank1Filtered(size(propertyIntegerRank1,dim=1),countSubsampled))
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
          allocate(propertyRealRank1Filtered   (size(propertyRealRank1   ,dim=1),countSubsampled))
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
  end subroutine subsampleOperate
