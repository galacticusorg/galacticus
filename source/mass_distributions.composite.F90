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
     type   (massDistributionList                   ), pointer :: massDistributions => null()
     type   (enumerationMassDistributionSymmetryType)          :: symmetry_
     logical                                                   :: isSingleComponent          , isCollisionless
   contains
     !![
     <methods>
       <method method="initialize" description="Initialize the mass distribution after construction."  />
       <method method="subset"     description="Return a subset of a composite mass distribution."     />
       <method method="describe"   description="Display a description of a composite mass distibution."/>
     </methods>
     !!]
     final     ::                            compositeDestructor
     procedure :: initialize              => compositeInitialize
     procedure :: subset                  => compositeSubset
     procedure :: describe                => compositeDescribe
     procedure :: matches                 => compositeMatches
     procedure :: symmetry                => compositeSymmetry
     procedure :: isSphericallySymmetric  => compositeIsSphericallySymmetric
     procedure :: isDimensionless         => compositeIsDimensionless
     procedure :: massTotal               => compositeMassTotal
     procedure :: acceleration            => compositeAcceleration
     procedure :: tidalTensor             => compositeTidalTensor
     procedure :: density                 => compositeDensity
     procedure :: surfaceDensity          => compositeSurfaceDensity
     procedure :: densityGradientRadial   => compositeDensityGradientRadial
     procedure :: densityRadialMoment     => compositeDensityRadialMoment
     procedure :: densitySphericalAverage => compositeDensitySphericalAverage
     procedure :: densitySquareIntegral   => compositeDensitySquareIntegral
     procedure :: potential               => compositePotential
     procedure :: energy                  => compositeEnergy
     procedure :: massEnclosedBySphere    => compositeMassEnclosedBySphere
     procedure :: radiusEnclosingMass     => compositeRadiusEnclosingMass
     procedure :: radiusEnclosingDensity  => compositeRadiusEnclosingDensity
     procedure :: rotationCurve           => compositeRotationCurve
     procedure :: rotationCurveGradient   => compositeRotationCurveGradient
     procedure :: chandrasekharIntegral   => compositeChandrasekharIntegral
     procedure :: positionSample          => compositePositionSample
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
    type(massDistributionComposite)                         :: self
    type(massDistributionList     ), pointer, intent(in   ) :: massDistributions
    type(massDistributionList     ), pointer                :: massDistribution_

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
    Initialize a composite mass distribution.
    !!}
    implicit none
    class  (massDistributionComposite              ), intent(inout) :: self
    type   (massDistributionList                   ), pointer       :: massDistribution_
    type   (enumerationMassDistributionSymmetryType)                :: symmetry_
    logical                                                         :: firstComponent   , haveKinematics
    
    ! Begin by assuming the highest degree of symmetry.
    self%symmetry_=massDistributionSymmetrySpherical
    ! Begin by assuming a single component.
    self%isSingleComponent=.true.
    firstComponent        =.true.
    ! Begin by assuming a collisionless distribution.
    self%isCollisionless  =.true.
    haveKinematics        =.true.
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
          ! Record if we have multiple components.
          if (firstComponent) then
             firstComponent=.false.
          else
             self%isSingleComponent=.false.
          end if
          ! Check for collisional components.
          if (associated(massDistribution_%massDistribution_%kinematicsDistribution_)) then
             if (massDistribution_%massDistribution_%kinematicsDistribution_%isCollisional()) self%isCollisionless=.false.
          else
             haveKinematics=.false.
          end if
          ! Move to the next mass distribution.
          massDistribution_ => massDistribution_%next
       end do
       ! Establish a kinematics distribution.
       if (haveKinematics) then
          if (self%isSingleComponent) then
             ! For a single component, simply use the kinematics distribution from that component
             call self%setKinematicsDistribution(self%massDistributions%massDistribution_%kinematicsDistribution_)
          else if (self%isCollisionless) then
             ! Construct a collisionless mass distribution.
             allocate(kinematicsDistributionCollisionless :: self%kinematicsDistribution_)
             select type (kinematicsDistribution_ => self%kinematicsDistribution_)
             type is (kinematicsDistributionCollisionless)
                !![
		<referenceConstruct owner="self" object="kinematicsDistribution_" nameAssociated="kinematicsDistribution_" constructor="kinematicsDistributionCollisionless()"/>
	        !!]
             end select
          end if
       end if
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

  logical function compositeIsSphericallySymmetric(self) result(isSphericallySymmetric)
    !!{
    Return true if the distribution is spherically symmetric.
    !!}
    implicit none
    class(massDistributionComposite), intent(inout) :: self

    isSphericallySymmetric=self%symmetry_ == massDistributionSymmetrySpherical
    return
  end function compositeIsSphericallySymmetric

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

  logical function compositeMatches(self,componentType,massType)
    !!{
    Return the total mass of a composite mass distribution.
    !!}
    implicit none
    class(massDistributionComposite   ), intent(inout)           :: self
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type (massDistributionList        ), pointer                 :: massDistribution_
    !![
    <optionalArgument name="componentType" defaultsTo="componentTypeAll"/>
    <optionalArgument name="massType"      defaultsTo="massTypeAll"     />
    !!]

    ! Assume no match by default.
    compositeMatches=.false.
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          ! If any component matches, report that as a match.
          if (massDistribution_ %massDistribution_%matches(componentType,massType)) then
             compositeMatches=.true.
             exit
          end if
          massDistribution_ => massDistribution_%next
       end do
    end if
    return
  end function compositeMatches

  subroutine compositeDescribe(self)
    !!{
    Display a description of a composite mass distribution.
    !!}
    use :: Display, only : displayMessage, displayIndent, displayUnindent
    implicit none
    class(massDistributionComposite), intent(inout) :: self
    type (massDistributionList     ), pointer       :: massDistribution_
    
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          select type (massDistribution__ => massDistribution_ %massDistribution_)
          class is (massDistributionComposite)
             call displayIndent (massDistribution_%massDistribution_%objectType())
             call massDistribution__%describe()
             call displayUnindent("")
          class default
             call displayMessage(massDistribution_%massDistribution_%objectType())
          end select
          massDistribution_ => massDistribution_%next
       end do
    end if
    return
  end subroutine compositeDescribe
  
  function compositeSubset(self,componentType,massType) result(subset)
    !!{
    Return the subset of the composite distribution that matches.
    !!}
    implicit none
    class(massDistributionClass       ), pointer                 :: subset
    class(massDistributionComposite   ), intent(inout)           :: self
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType
    class(massDistributionClass       ), pointer                 :: massDistribution___
    type (massDistributionList        ), pointer                 :: subsetHead         , subsetNext, &
         &                                                          massDistribution_
    !![
    <optionalArgument name="componentType" defaultsTo="componentTypeAll"/>
    <optionalArgument name="massType"      defaultsTo="massTypeAll"     />
    !!]

    subset => null()
    if (associated(self%massDistributions)) then
       subsetHead        => null()
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          select type (massDistribution__ => massDistribution_ %massDistribution_)
          class is (massDistributionComposite)
             massDistribution___ => massDistribution__%subset(componentType_,massType_)
             if (associated(massDistribution___)) then
                if (associated(subsetHead)) then
                   allocate(subsetNext%next)
                   subsetNext => subsetNext%next
                else
                   allocate(subsetHead     )
                   subsetNext => subsetHead
                end if
                subsetNext%massDistribution_ => massDistribution___
                ! !![
                ! <referenceCountIncrement owner="subsetNext" object="massDistribution_"/>
                ! !!]
             end if
          class default
             if (massDistribution__%matches(componentType,massType)) then
                if (associated(subsetHead)) then
                   allocate(subsetNext%next)
                   subsetNext => subsetNext%next
                else
                   allocate(subsetHead     )
                   subsetNext => subsetHead
                end if
                subsetNext%massDistribution_ => massDistribution__
                ! !![
                ! <referenceCountIncrement owner="subsetNext" object="massDistribution_"/>
                ! !!]
             end if
          end select
          massDistribution_ => massDistribution_%next
       end do
       if (associated(subsetHead)) then
          allocate(massDistributionComposite :: subset)
          select type(subset)
          type is (massDistributionComposite)
             !![
             <referenceConstruct object="subset" constructor="massDistributionComposite(subsetHead)"/>
             !!]
          end select
       end if
    end if
    return
  end function compositeSubset

  double precision function compositeMassTotal(self)
    !!{
    Return the total mass of a composite mass distribution.
    !!}
    implicit none
    class(massDistributionComposite), intent(inout) :: self
    type (massDistributionList     ), pointer       :: massDistribution_

    compositeMassTotal=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeMassTotal =  +compositeMassTotal                               &
               &                +massDistribution_ %massDistribution_%massTotal()
          massDistribution_  =>  massDistribution_ %next
       end do
    end if
    return
  end function compositeMassTotal

  double precision function compositeDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a composite mass distribution.
    !!}
    implicit none
    class(massDistributionComposite), intent(inout) :: self
    class(coordinate               ), intent(in   ) :: coordinates
    type (massDistributionList     ), pointer       :: massDistribution_

    compositeDensity=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeDensity  =  +compositeDensity                                          &
               &               +massDistribution_ %massDistribution_%density(coordinates)          
          massDistribution_ =>  massDistribution_ %next
       end do
    end if
    return
  end function compositeDensity

  double precision function compositeDensitySphericalAverage(self,radius)
    !!{
    Return the spherically-averaged density at the specified {\normalfont \ttfamily radius} in a composite mass distribution.
    !!}
    implicit none
    class           (massDistributionComposite), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
    type            (massDistributionList     ), pointer       :: massDistribution_

    compositeDensitySphericalAverage=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeDensitySphericalAverage  =  +compositeDensitySphericalAverage                                                   &
               &                               +massDistribution_               %massDistribution_%densitySphericalAverage(radius)          
          massDistribution_                 =>  massDistribution_               %next
       end do
    end if
    return
  end function compositeDensitySphericalAverage

  double precision function compositeDensitySquareIntegral(self,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Return the integral of the square of the density within the given radial interval.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionComposite), intent(inout)           :: self
    double precision                           , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                    , intent(  out), optional :: isInfinite

    if (self%isSingleComponent) then
       compositeDensitySquareIntegral=self%massDistributions%massDistribution_%densitySquareIntegral(radiusMinimum,radiusMaximum)
    else
       compositeDensitySquareIntegral=0.0d0
       call Error_Report('support for ∫ dr ρ²(r) of multiple components is not implemented'//{introspection:location})
    end if
    return
  end function compositeDensitySquareIntegral

  double precision function compositeSurfaceDensity(self,coordinates)
    !!{
    Return the surface density at the specified {\normalfont \ttfamily coordinates} in a composite mass distribution.
    !!}
    use :: Coordinates, only : coordinate
    implicit none
    class(massDistributionComposite), intent(inout) :: self
    class(coordinate               ), intent(in   ) :: coordinates
    type (massDistributionList     ), pointer       :: massDistribution_


    compositeSurfaceDensity=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeSurfaceDensity  =  +compositeSurfaceDensity                                          &
               &                      +massDistribution_ %massDistribution_%surfaceDensity(coordinates)          
          massDistribution_        =>  massDistribution_ %next
       end do
    end if
    return
  end function compositeSurfaceDensity

  double precision function compositeDensityGradientRadial(self,coordinates,logarithmic)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a composite mass distribution.
    !!}
    implicit none
    class  (massDistributionComposite), intent(inout), target   :: self
    class  (coordinate               ), intent(in   )           :: coordinates
    logical                           , intent(in   ), optional :: logarithmic
    type   (massDistributionList     ), pointer                 :: massDistribution_
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]
    
    compositeDensityGradientRadial=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeDensityGradientRadial  =  +compositeDensityGradientRadial                                                                          &
               &                             +massDistribution_             %massDistribution_%densityGradientRadial(coordinates,logarithmic=.false.)
          massDistribution_               =>  massDistribution_             %next
       end do
       if (logarithmic_) then
          compositeDensityGradientRadial=+compositeDensityGradientRadial                         &
               &                         *coordinates                   %rSpherical(           ) &
               &                         /self                          %density   (coordinates)
       end if
    end if
    return
  end function compositeDensityGradientRadial

  double precision function compositeDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Return the given radial density moment in a composite mass distribution.
    !!}
    implicit none
    class           (massDistributionComposite), intent(inout)           :: self
    double precision                           , intent(in   )           :: moment
    double precision                           , intent(in   ), optional :: radiusMinimum    , radiusMaximum
    logical                                    , intent(  out), optional :: isInfinite
    type            (massDistributionList     ), pointer                 :: massDistribution_

    compositeDensityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeDensityRadialMoment =  +compositeDensityRadialMoment                                                                                      &
               &                          +massDistribution_           %massDistribution_%densityRadialMoment(moment,radiusMinimum,radiusMaximum,isInfinite)
          if (present(isInfinite)) then
             if (isInfinite) return
          end if
          massDistribution_            =>  massDistribution_           %next
       end do
    end if
    return
  end function compositeDensityRadialMoment

  double precision function compositeMassEnclosedBySphere(self,radius)
    !!{
    Computes the mass enclosed within a sphere in a composite mass distribution.
    !!}
    implicit none
    class           (massDistributionComposite), intent(inout), target   :: self
    double precision                           , intent(in   )           :: radius
    type            (massDistributionList     )               , pointer  :: massDistribution_

    compositeMassEnclosedBySphere=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeMassEnclosedBySphere =  +compositeMassEnclosedBySphere                                     &
               &                           +massDistribution_ %massDistribution_%massEnclosedBySphere(radius)
          massDistribution_             =>  massDistribution_ %next
       end do
    end if
    return
  end function compositeMassEnclosedBySphere

  double precision function compositeRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for composite mass distributions.
    !!}    
    implicit none
    class           (massDistributionComposite), intent(inout), target   :: self
    double precision                           , intent(in   ), optional :: mass, massFractional

    if (self%isSingleComponent) then
       radius=self%massDistributions%massDistribution_%radiusEnclosingMass         (mass,massFractional)
    else
       radius=self                                    %radiusEnclosingMassNumerical(mass,massFractional)
    end if
    return
  end function compositeRadiusEnclosingMass

  double precision function compositeRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for composite mass distributions.
    !!}    
    implicit none
    class           (massDistributionComposite), intent(inout), target   :: self
    double precision                           , intent(in   )           :: density
    double precision                           , intent(in   ), optional :: radiusGuess

    if (self%isSingleComponent) then
       radius=self%massDistributions%massDistribution_%radiusEnclosingDensity         (density,radiusGuess)
    else
       radius=self                                    %radiusEnclosingDensityNumerical(density,radiusGuess)
    end if
    return
  end function compositeRadiusEnclosingDensity
  
  double precision function compositeRotationCurve(self,radius)
    !!{
    Return the rotation curve for a composite mass distribution.
    !!}
    implicit none
    class           (massDistributionComposite), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
    type            (massDistributionList     ), pointer       :: massDistribution_

    compositeRotationCurve=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeRotationCurve =  +compositeRotationCurve                                        &
               &                    +massDistribution_ %massDistribution_%rotationCurve(radius)**2
          massDistribution_      =>  massDistribution_ %next
       end do
       compositeRotationCurve=sqrt(compositeRotationCurve)
    end if
    return
  end function compositeRotationCurve

  double precision function compositeRotationCurveGradient(self,radius)
    !!{
    Return the gradient of the rotation curve for a composite mass distribution.
    !!}
    implicit none
    class           (massDistributionComposite), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
    type            (massDistributionList     ), pointer       :: massDistribution_
    
    compositeRotationCurveGradient=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeRotationCurveGradient =  +compositeRotationCurveGradient                                                 &
               &                            +massDistribution_             %massDistribution_%rotationCurveGradient(radius)
          massDistribution_              =>  massDistribution_             %next
       end do
    end if
    return
  end function compositeRotationCurveGradient

  function compositeAcceleration(self,coordinates)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for a composite mass distribution.
    !!}
    implicit none
    double precision                              , dimension(3  ) :: compositeAcceleration
    class           (massDistributionComposite   ), intent(inout)  :: self
    class           (coordinate                  ), intent(in   )  :: coordinates
    type            (massDistributionList        ), pointer        :: massDistribution_

    compositeAcceleration=[0.0d0,0.0d0,0.0d0]
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeAcceleration  =  +compositeAcceleration                                             &
               &                    +massDistribution_    %massDistribution_%acceleration(coordinates)
          massDistribution_      =>  massDistribution_    %next
       end do
    end if   
    return
  end function compositeAcceleration

  function compositeTidalTensor(self,coordinates)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} for exponential disk mass distributions.
    !!}
    implicit none
    type (tensorRank2Dimension3Symmetric)                :: compositeTidalTensor
    class(massDistributionComposite     ), intent(inout) :: self
    class(coordinate                    ), intent(in   ) :: coordinates
    type (massDistributionList          ), pointer       :: massDistribution_

    compositeTidalTensor=tensorRank2Dimension3Symmetric()
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeTidalTensor  =  +compositeTidalTensor                                            &
               &                   +massDistribution_   %massDistribution_%tidalTensor(coordinates)
          massDistribution_     =>  massDistribution_   %next
       end do
    end if
    return
  end function compositeTidalTensor
  
  double precision function compositePotential(self,coordinates,status)
    !!{
    Return the gravitational potential for a composite mass distribution.
    !!}
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    implicit none
    class(massDistributionComposite        ), intent(inout), target   :: self
    class(coordinate                       ), intent(in   )           :: coordinates
    type (enumerationStructureErrorCodeType), intent(  out), optional :: status
    type (massDistributionList             ), pointer                 :: massDistribution_

    if (present(status)) status=structureErrorCodeSuccess
    compositePotential=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositePotential  =  +compositePotential                                                 &
               &                 +massDistribution_ %massDistribution_%potential(coordinates,status)
          if (present(status).and.status /= structureErrorCodeSuccess) return
          massDistribution_   =>  massDistribution_ %next
       end do
    end if
    return
  end function compositePotential

  double precision function compositeEnergy(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the energy of the mass distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionComposite), intent(inout), target  :: self
    double precision                           , intent(in   )          :: radiusOuter
    class           (massDistributionClass    ), intent(inout), target  :: massDistributionEmbedding
    class           (massDistributionClass    )               , pointer :: self_

    self_ => self
    if (self%isSingleComponent .and. associated(self_,massDistributionEmbedding)) then
       energy=self%massDistributions%massDistribution_%energy(radiusOuter,self%massDistributions%massDistribution_)
    else
       energy=0.0d0
       call Error_Report('support for energy of multiple components is not implemented'//{introspection:location})
    end if
    return
  end function compositeEnergy
  
  function compositeChandrasekharIntegral(self,massDistributionEmbedding,massDistributionPerturber,massPerturber,coordinates,velocity)
    !!{
    Compute the Chandrasekhar integral at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    implicit none
    double precision                              , dimension(3)  :: compositeChandrasekharIntegral
    class           (massDistributionComposite   ), intent(inout) :: self
    class           (massDistributionClass       ), intent(inout) :: massDistributionEmbedding     , massDistributionPerturber
    double precision                              , intent(in   ) :: massPerturber
    class           (coordinate                  ), intent(in   ) :: coordinates                   , velocity
    type            (massDistributionList        ), pointer       :: massDistribution_

    compositeChandrasekharIntegral=0.0d0
    if (associated(self%massDistributions)) then
       massDistribution_ => self%massDistributions
       do while (associated(massDistribution_))
          compositeChandrasekharIntegral =  +compositeChandrasekharIntegral                                                                                                                              &
               &                            +massDistribution_%massDistribution_%chandrasekharIntegral(massDistribution_%massDistribution_,massDistributionPerturber,massPerturber,coordinates,velocity)
          massDistribution_              =>  massDistribution_%next
       end do
    end if
    return
  end function compositeChandrasekharIntegral

  function compositePositionSample(self,randomNumberGenerator_)
    !!{
    Sample a position from a composite distribution.
    !!}
    implicit none
    double precision                              , dimension(3)  :: compositePositionSample
    class           (massDistributionComposite   ), intent(inout) :: self
    class           (randomNumberGeneratorClass  ), intent(inout) :: randomNumberGenerator_
    type            (massDistributionList        ), pointer       :: massDistribution_
    double precision                                              :: massCumulative

    compositePositionSample=[0.0d0,0.0d0,0.0d0]
    if (associated(self%massDistributions)) then
       massCumulative    =  +self                  %massTotal    () &
            &               *randomNumberGenerator_%uniformSample()
       massDistribution_ =>  self%massDistributions
       do while (associated(massDistribution_))
          massCumulative=+                                    massCumulative   &
               &         -massDistribution_%massDistribution_%massTotal     ()
          if (massCumulative <= 0.0d0) then
             compositePositionSample=massDistribution_%massDistribution_%positionSample(randomNumberGenerator_)
             return
          end if
          massDistribution_ => massDistribution_%next
       end do
    end if
    
    return
  end function compositePositionSample
