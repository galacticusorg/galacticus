!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implementation of a composite mass distribution class.
  !!}
  
  !![
  <massDistribution name="massDistributionComposite">
    <description>A composite mass distribution.</description>
    <linkedList type="massDistributionList" variable="massDistributions" next="next" object="massDistribution_" objectType="massDistributionClass"/>
  </massDistribution>
  !!]

  type, public :: massDistributionList
     class(massDistributionClass), pointer :: massDistribution_ => null()
     type (massDistributionList ), pointer :: next              => null()
  end type massDistributionList


  type, public, extends(massDistributionClass) :: massDistributionComposite
     !!{
     A composite mass distribution class.
     !!}
     type(massDistributionList                   ), pointer :: massDistributions => null()
     type(enumerationMassDistributionSymmetryType)          :: symmetry_
   contains
     !![
     <methods>
       <method method="initialize" description="Initialize the mass distribution after construction."/>
     </methods>
     !!]
     final     ::                          compositeDestructor
     procedure :: initialize            => compositeInitialize
     procedure :: symmetry              => compositeSymmetry
     procedure :: isDimensionless       => compositeIsDimensionless
     procedure :: massTotal             => compositeMassTotal
     procedure :: acceleration          => compositeAcceleration
     procedure :: tidalTensor           => compositeTidalTensor
     procedure :: density               => compositeDensity
     procedure :: densityGradientRadial => compositeDensityGradientRadial
     procedure :: densityRadialMoment   => compositeDensityRadialMoment
     procedure :: potential             => compositePotential
     procedure :: massEnclosedBySphere  => compositeMassEnclosedBySphere
     procedure :: positionSample        => compositePositionSample
  end type massDistributionComposite

  interface massDistributionComposite
     !!{
     Constructors for the {\normalfont \ttfamily composite} mass distribution class.
     !!}
     module procedure compositeConstructorParameters
     module procedure compositeConstructorInternal
  end interface massDistributionComposite

contains

  function compositeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily composite} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (massDistributionComposite)                :: self
    type   (inputParameters          ), intent(inout) :: parameters
    type   (massDistributionList     ), pointer       :: massDistribution_
    integer                                           :: i

    massDistribution_ => null()
    do i=1,parameters%copiesCount('massDistribution',zeroIfNotPresent=.true.)
       if (associated(massDistribution_)) then
          allocate(massDistribution_%next)
          massDistribution_ => massDistribution_%next
       else
          allocate(self%massDistributions)
          massDistribution_ => self             %massDistributions
       end if
       !![
       <objectBuilder class="massDistribution" name="massDistribution_%massDistribution_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="mmassDistribution"/>
    !!]
    call self%initialize()
    return
  end function compositeConstructorParameters

  function compositeConstructorInternal(massDistributions) result(self)
    !!{
    Internal constructor for ``composite'' mass distribution class.
    !!}
    type(massDistributionComposite)                        :: self
    type(massDistributionList     ), target, intent(in   ) :: massDistributions
    type(massDistributionList     ), pointer               :: massDistribution_

    self             %massDistributions => massDistributions
    massDistribution_                   => massDistributions
    do while (associated(massDistribution_))
       !![
       <referenceCountIncrement owner="massDistribution_" object="massDistribution_"/>
       !!]
       massDistribution_ => massDistribution_%next
    end do
    call self%initialize()
    return
  end function compositeConstructorInternal

  subroutine compositeDestructor(self)
    !!{
    Destructor for composite mass distributions.
    !!}
    implicit none
    type(massDistributionComposite), intent(inout) :: self
    type(massDistributionList     ), pointer       :: massDistribution_, massDistributionNext

    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          massDistributionNext => massDistribution_%next
          !![
          <objectDestructor name="massDistribution_%massDistribution_"/>
          !!]
          deallocate(massDistribution_)
          massDistribution_ => massDistributionNext
       end do
    end if
    return
  end subroutine compositeDestructor

  subroutine compositeInitialize(self)
    !!{
    Destructor for composite mass distributions.
    !!}
    implicit none
    class(massDistributionComposite              ), intent(inout) :: self
    type (massDistributionList                   ), pointer       :: massDistribution_
    type (enumerationMassDistributionSymmetryType)                :: symmetry_

    ! Begin by assuming the highest degree of symmetry.
    self%symmetry_=massDistributionSymmetrySpherical
    ! Examine each distribution.
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          ! Dimensionless mass distributions are not allowed.
          if (massDistribution_%massDistribution_%isDimensionless())                                                                         &
               & call Error_Report('dimensionless mass distributions can not be part of a composite distribution'//{introspection:location})
          ! Determine the symmetry of the composite distribution.
          symmetry_=massDistribution_%massDistribution_%symmetry()
          select case (symmetry_%ID)
          case (massDistributionSymmetryNone       %ID)
             ! Distribution has no symmetry. Therefore, the composite distribution also has no symmetry.
             self%symmetry_=massDistributionSymmetryNone
          case (massDistributionSymmetryCylindrical%ID)
             ! Distribution has spherical symmetry. If the composite distribution so far is spherically symmetric, we can reduce
             ! the symmetry of the overall system to cylindrical.
             if (self%symmetry_ == massDistributionSymmetrySpherical)  &
                  & self%symmetry_=massDistributionSymmetryCylindrical
          case (massDistributionSymmetrySpherical  %ID)
             ! This is the highest symmetry possible - no need to change anything.
          case default
             call Error_Report('unknown symmetry'//{introspection:location})
          end select
          ! Move to the next mass distribution.
          massDistribution_ => massDistribution_%next
       end do
    end if
    return
  end subroutine compositeInitialize

  function compositeSymmetry(self) result(symmetry)
    !!{
    Return the symmetry of a composite mass distribution.
    !!}
    implicit none
    type (enumerationMassDistributionSymmetryType)                :: symmetry
    class(massDistributionComposite              ), intent(inout) :: self

    symmetry=self%symmetry_
    return
  end function compositeSymmetry

  logical function compositeIsDimensionless(self)
    !!{
    Return the dimensionless nature of a composite mass distribution.
    !!}
    implicit none
    class(massDistributionComposite), intent(inout) :: self

    ! Composite distributions are never dimensionless.
    compositeIsDimensionless=.false.
    return
  end function compositeIsDimensionless

  double precision function compositeMassTotal(self,componentType,massType)
    !!{
    Return the total mass of a composite mass distribution.
    !!}
    implicit none
    class(massDistributionComposite   ), intent(inout)           :: self
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type (massDistributionList        ), pointer                 :: massDistribution_

    compositeMassTotal=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          if (massDistribution_ %massDistribution_%matches(componentType,massType))     &
               & compositeMassTotal =  +compositeMassTotal                              &
               &                       +massDistribution_ %massDistribution_%massTotal()
          massDistribution_         =>  massDistribution_ %next
       end do
    end if
    return
  end function compositeMassTotal

  double precision function compositeDensity(self,coordinates,componentType,massType)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a composite mass distribution.
    !!}
    implicit none
    class(massDistributionComposite   ), intent(inout)           :: self
    class(coordinate                  ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type (massDistributionList        ), pointer                 :: massDistribution_

    compositeDensity=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          if (massDistribution_ %massDistribution_%matches(componentType,massType))              &
               & compositeDensity  =  +compositeDensity                                          &
               &                      +massDistribution_ %massDistribution_%density(coordinates)
          massDistribution_        =>  massDistribution_ %next
       end do
    end if
    return
  end function compositeDensity

  double precision function compositeDensityGradientRadial(self,coordinates,logarithmic,componentType,massType)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a composite mass distribution.
    !!}
    implicit none
    class           (massDistributionComposite  ), intent(inout)           :: self
    class           (coordinate                 ), intent(in   )           :: coordinates
    logical                                      , intent(in   ), optional :: logarithmic
    type           (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type           (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (massDistributionList       ), pointer                 :: massDistribution_
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]
    
    compositeDensityGradientRadial=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          if (massDistribution_ %massDistribution_%matches(componentType,massType))                                                                          &
               & compositeDensityGradientRadial  =  +compositeDensityGradientRadial                                                                          &
               &                                    +massDistribution_             %massDistribution_%densityGradientRadial(coordinates,logarithmic=.false.)
          massDistribution_                      =>  massDistribution_             %next
       end do
       if (logarithmic_) then
          compositeDensityGradientRadial=+compositeDensityGradientRadial                         &
               &                         *coordinates                   %rSpherical(           ) &
               &                         /self                          %density   (coordinates)
       end if
    end if
    return
  end function compositeDensityGradientRadial

  double precision function compositeDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType)
    !!{
    Return the given radial density moment in a composite mass distribution.
    !!}
    implicit none
    class           (massDistributionComposite  ), intent(inout)           :: self
    double precision                             , intent(in   )           :: moment
    double precision                             , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                      , intent(  out), optional :: isInfinite
    type           (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type           (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (massDistributionList       ), pointer                 :: massDistribution_

    compositeDensityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          if (massDistribution_ %massDistribution_%matches(componentType,massType))                                                                                 &
               & compositeDensityRadialMoment =  +compositeDensityRadialMoment                                                                                      &
               &                                 +massDistribution_           %massDistribution_%densityRadialMoment(moment,radiusMinimum,radiusMaximum,isInfinite)
          if (present(isInfinite)) then
             if (isInfinite) return
          end if
          massDistribution_            =>  massDistribution_           %next
       end do
    end if
    return
  end function compositeDensityRadialMoment

  double precision function compositeMassEnclosedBySphere(self,radius,componentType,massType)
    !!{
    Computes the mass enclosed within a sphere in a composite mass distribution.
    !!}
    implicit none
    class           (massDistributionComposite  ), intent(inout), target   :: self
    double precision                             , intent(in   )           :: radius
    type           (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type           (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (massDistributionList       )               , pointer  :: massDistribution_

    compositeMassEnclosedBySphere=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          if (massDistribution_ %massDistribution_%matches(componentType,massType))                                   &
               & compositeMassEnclosedBySphere  =  +compositeMassEnclosedBySphere                                     &
               &                                   +massDistribution_ %massDistribution_%massEnclosedBySphere(radius)
          massDistribution_                     =>  massDistribution_ %next
       end do
    end if
    return
  end function compositeMassEnclosedBySphere

  function compositeAcceleration(self,coordinates,componentType,massType)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for a composite mass distribution.
    !!}
    implicit none
    double precision                              , dimension(3  )          :: compositeAcceleration
    class           (massDistributionComposite   ), intent(inout)           :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (massDistributionList        ), pointer                 :: massDistribution_

    compositeAcceleration=[0.0d0,0.0d0,0.0d0]
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          if (massDistribution_ %massDistribution_%matches(componentType,massType))                           &
               & compositeAcceleration  =  +compositeAcceleration                                             &
               &                           +massDistribution_    %massDistribution_%acceleration(coordinates)
          massDistribution_             =>  massDistribution_    %next
       end do
    end if   
    return
  end function compositeAcceleration

  function compositeTidalTensor(self,coordinates,componentType,massType)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} for exponential disk mass distributions.
    !!}
    implicit none
    type (tensorRank2Dimension3Symmetric)                          :: compositeTidalTensor
    class(massDistributionComposite     ), intent(inout)           :: self
    class(coordinate                    ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType  ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType       ), intent(in   ), optional :: massType
    type (massDistributionList          ), pointer                 :: massDistribution_

    compositeTidalTensor=tensorRank2Dimension3Symmetric()
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          if (massDistribution_ %massDistribution_%matches(componentType,massType))                        &
               & compositeTidalTensor  =  +compositeTidalTensor                                            &
               &                          +massDistribution_   %massDistribution_%tidalTensor(coordinates)
          massDistribution_            =>  massDistribution_   %next
       end do
    end if
    return
  end function compositeTidalTensor
  
  double precision function compositePotential(self,coordinates,componentType,massType)
    !!{
    Return the gravitational potential for a composite mass distribution.
    !!}
    implicit none
    class(massDistributionComposite   ), intent(inout)           :: self
    class(coordinate                  ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type (massDistributionList        ), pointer                 :: massDistribution_

    compositePotential=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          if (massDistribution_ %massDistribution_%matches(componentType,massType))                  &
               & compositePotential  =  +compositePotential                                          &
               &                        +massDistribution_ %massDistribution_%potential(coordinates)
          massDistribution_          =>  massDistribution_ %next
       end do
    end if
    return
  end function compositePotential
  
  function compositePositionSample(self,randomNumberGenerator_,componentType,massType)
    !!{
    Sample a position from a composite distribution.
    !!}
    implicit none
    double precision                              , dimension(3)            :: compositePositionSample
    class           (massDistributionComposite   ), intent(inout)           :: self
    class           (randomNumberGeneratorClass  ), intent(inout)           :: randomNumberGenerator_
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (massDistributionList        ), pointer                 :: massDistribution_
    double precision                                                        :: massCumulative

    compositePositionSample=[0.0d0,0.0d0,0.0d0]
    if (associated(self%massDistributions)) then
       massCumulative    =  +self                  %massTotal    () &
            &               *randomNumberGenerator_%uniformSample()
       massDistribution_ =>  self%massDistributions
       do while (associated(massDistribution_))
          if (massDistribution_ %massDistribution_%matches(componentType,massType))   &
               & massCumulative=+                                    massCumulative   &
               &                -massDistribution_%massDistribution_%massTotal     ()
          if (massCumulative <= 0.0d0) then
             compositePositionSample=massDistribution_%massDistribution_%positionSample(randomNumberGenerator_)
             return
          end if
          massDistribution_   =>  massDistribution_ %next
       end do
    end if
    
    return
  end function compositePositionSample
