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

!% Contains a module which implements an N-body data operator which subsamples points at a given rate.

  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
 
  !# <nbodyOperator name="nbodyOperatorSubsample">
  !#  <description>An N-body data operator which filters out particles based on a property range.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorSubsample
     !% An N-body data operator which filters out particles based on a property range.
     private
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
     double precision                                      :: rate
  contains
    final :: subsampleDestructor
     procedure :: operate => subsampleOperate
  end type nbodyOperatorSubsample

  interface nbodyOperatorSubsample
     !% Constructors for the {\normalfont \ttfamily subsample} N-body operator class.
     module procedure subsampleConstructorParameters
     module procedure subsampleConstructorInternal
  end interface nbodyOperatorSubsample

contains
  
  function subsampleConstructorParameters(parameters) result (self)
    !% Constructor for the {\normalfont \ttfamily subsample} N-body operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nbodyOperatorSubsample    )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    double precision                                            :: rate

    !# <inputParameter>
    !#   <name>rate</name>
    !#   <source>parameters</source>
    !#   <description>The rate at which to subsample points.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    self=nbodyOperatorSubsample(rate,randomNumberGenerator_)
    !# <objectDestructor name="randomNumberGenerator_"/>
    !# <inputParametersValidate source="parameters"/>
    return
  end function subsampleConstructorParameters

  function subsampleConstructorInternal(rate,randomNumberGenerator_) result (self)
    !% Internal constructor for the {\normalfont \ttfamily subsample} N-body operator class.
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    implicit none
    type            (nbodyOperatorSubsample    )                        :: self
    class           (randomNumberGeneratorClass), intent(in   ), target :: randomNumberGenerator_
    double precision                            , intent(in   )         :: rate
    !# <constructorAssign variables="rate, *randomNumberGenerator_"/>

    if (rate <= 0.0d0 .or. rate > 1.0d0) call Galacticus_Error_Report('range must be in the interval [0,1)'//{introspection:location})
    return
  end function subsampleConstructorInternal

  subroutine subsampleDestructor(self)
    !% Destructor for the {\normalfont \ttfamily subsample} N-body operator class.
    implicit none
    type(nbodyOperatorSubsample), intent(inout) :: self

    !# <objectDestructor name="self%randomNumberGenerator_"/>
    return
  end subroutine subsampleDestructor

  subroutine subsampleOperate(self,simulations)
    !% Identify and flag particles which have been always isolated.
    use :: Galacticus_Display   , only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, verbosityStandard
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    use :: NBody_Simulation_Data, only : propertyTypeInteger      , propertyTypeReal
    implicit none
    class           (nbodyOperatorSubsample), intent(inout)                 :: self
    type            (nBodyData             ), intent(inout), dimension(  :) :: simulations
    logical                                 , allocatable  , dimension(  :) :: mask
    integer         (c_size_t              ), allocatable  , dimension(  :) :: propertyInteger, propertyIntegerSubsampled
    double precision                        , allocatable  , dimension(  :) :: propertyReal   , propertyRealSubsampled
    double precision                        , allocatable  , dimension(:,:) :: position
    integer                                                                 :: i              , j
    integer         (c_size_t              )                                :: k              , countSubsampled
    
    call Galacticus_Display_Indent('subsample points',verbosityStandard)
    do i=1,size(simulations)
       allocate(mask(size(simulations(i)%particleIDs)))
       do k=1,size(mask)
          mask(k)=self%randomNumberGenerator_%uniformSample() < self%rate
       end do
       ! Subsample default properties.
       countSubsampled=count(mask)
       !! Position
       allocate(position(3,countSubsampled))
       do j=1,3
          position(j,:)=pack(simulations(i)%position(j,:),mask)
       end do
       deallocate(simulations(i)%position)
       call move_alloc(position,simulations(i)%position)
       !! Velocity
       allocate(position(3,countSubsampled))
       do j=1,3
          position(j,:)=pack(simulations(i)%velocity(j,:),mask)
       end do
       deallocate(simulations(i)%velocity)
       call move_alloc(position,simulations(i)%velocity)
       !! IDs
       allocate(propertyInteger(countSubsampled))
       propertyInteger=pack(simulations(i)%particleIDs,mask)
       deallocate(simulations(i)%particleIDs)
       call move_alloc(propertyInteger,simulations(i)%particleIDs)
       ! Subsample other properties.
       !! Integer properties.
       allocate(propertyIntegerSubsampled(countSubsampled))
       do j=1,simulations(i)%propertiesInteger%size()
          propertyInteger          =simulations(i)%propertiesInteger%value(j)
          propertyIntegerSubsampled=pack(propertyInteger,mask)
          call simulations(i)%propertiesInteger%set(simulations(i)%propertiesInteger%key(j),propertyIntegerSubsampled)
       end do
       deallocate(propertyIntegerSubsampled)
       !! Real properties.
       allocate(propertyRealSubsampled   (countSubsampled))
       do j=1,simulations(i)%propertiesReal   %size()
          propertyReal          =simulations(i)%propertiesReal   %value(j)
          propertyRealSubsampled=pack(propertyReal   ,mask)
          call simulations(i)%propertiesReal   %set(simulations(i)%propertiesReal   %key(j),propertyRealSubsampled   )
       end do
       deallocate(propertyRealSubsampled)
       deallocate(mask                  )
    end do
    call Galacticus_Display_Unindent('done',verbosityStandard)
    return
  end subroutine subsampleOperate
