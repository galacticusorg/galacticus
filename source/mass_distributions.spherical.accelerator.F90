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
  Implements an accelerator for spherical mass distributions.
  !!}

  use :: Binary_Search_Trees , only : binaryTree

  !![
  <massDistribution name="massDistributionSphericalAccelerator">
   <description>
     Accelerates spherical mass distribution classes by storing previous results for the enclosed mass and interpolating where
     possible.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalDecorator) :: massDistributionSphericalAccelerator
     !!{
     Implementation of a finite resolution spherical mass distribution.
     !!}
     private
     type            (binaryTree) :: treeMassEnclosed
     double precision             :: toleranceRelative, factorRadiusMaximum, factorRadiusLogarithmicMaximum
   contains
     final     ::                         sphericalAcceleratorDestructor
     procedure :: density              => sphericalAcceleratorDensity
     procedure :: massEnclosedBySphere => sphericalAcceleratorMassEnclosedBySphere
     procedure :: useUndecorated       => sphericalAcceleratorUseUndecorated
  end type massDistributionSphericalAccelerator

  interface massDistributionSphericalAccelerator
     !!{
     Constructors for the {\normalfont \ttfamily sphericalAccelerator} mass distribution class.
     !!}
     module procedure sphericalAcceleratorConstructorParameters
     module procedure sphericalAcceleratorConstructorInternal
  end interface massDistributionSphericalAccelerator

contains

  function sphericalAcceleratorConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sphericalAccelerator} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalAccelerator)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (massDistributionClass               ), pointer       :: massDistribution_
    double precision                                                      :: toleranceRelative, factorRadiusMaximum
    type            (varying_string                      )                :: componentType    , massType           , &
         &                                                                   nonAnalyticSolver

    !![
    <inputParameter>
      <name>toleranceRelative</name>
      <defaultValue>1.0d-2</defaultValue>
      <source>parameters</source>
      <description>The tolerance with which to accept accelerated estimates.</description>
    </inputParameter>
    <inputParameter>
      <name>factorRadiusMaximum</name>
      <defaultValue>3.0d0</defaultValue>
      <source>parameters</source>
      <description>The maximum factor by which to interpolate in radius.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="massDistribution" name="massDistribution_" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       self=massDistributionSphericalAccelerator(toleranceRelative,factorRadiusMaximum,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function sphericalAcceleratorConstructorParameters
  
  function sphericalAcceleratorConstructorInternal(toleranceRelative,factorRadiusMaximum,nonAnalyticSolver,massDistribution_,componentType,massType) result(self)
    !!{
    Constructor for {\normalfont \ttfamily sphericalAccelerator} mass distribution class.
    !!}
    implicit none
    type            (massDistributionSphericalAccelerator)                          :: self
    class           (massDistributionSpherical           ), intent(in   ), target   :: massDistribution_
    double precision                                      , intent(in   )           :: toleranceRelative, factorRadiusMaximum
    type            (enumerationNonAnalyticSolversType   ), intent(in   )           :: nonAnalyticSolver
    type            (enumerationComponentTypeType        ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType             ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="toleranceRelative, factorRadiusMaximum, *massDistribution_, nonAnalyticSolver, componentType, massType"/>
    !!]

    self%factorRadiusLogarithmicMaximum=+log(sqrt(factorRadiusMaximum))
    self%dimensionless                 =self%massDistribution_%isDimensionless()
    return
  end function sphericalAcceleratorConstructorInternal

  subroutine sphericalAcceleratorDestructor(self)
    !!{
    Destructor for the abstract {\normalfont \ttfamily massDistributionSphericalAccelerator} class.
    !!}
    implicit none
    type(massDistributionSphericalAccelerator), intent(inout) :: self
    
    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine sphericalAcceleratorDestructor

  logical function sphericalAcceleratorUseUndecorated(self) result(useUndecorated)
    !!{
    Determines whether to use the undecorated solution.
    !!}
    implicit none
    class(massDistributionSphericalAccelerator), intent(inout) :: self

    useUndecorated=.false.
    return
  end function sphericalAcceleratorUseUndecorated

  double precision function sphericalAcceleratorDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an accelerated mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalAccelerator), intent(inout) :: self
    class(coordinate                          ), intent(in   ) :: coordinates

    density=self%massDistribution_%density(coordinates)
    return
  end function sphericalAcceleratorDensity

  double precision function sphericalAcceleratorMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for accelerated mass distributions.
    !!}
    use :: Binary_Search_Trees , only : binaryTreeNode
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (massDistributionSphericalAccelerator), intent(inout), target :: self
    double precision                                      , intent(in   )         :: radius
    type            (binaryTreeNode                      ), pointer               :: left1            , left2        , &
         &                                                                           right1           , right2
    double precision                                                              :: massEnclosed1    , massEnclosed2, &
         &                                                                           radiusLogarithmic
    logical                                                                       :: found

    found            =.false.
    radiusLogarithmic=log(radius)
    call self%treeMassEnclosed%bracket(radiusLogarithmic,left1,right1)
    if (associated(left1).and.associated(right1)) then
       if (associated(left1,right1)) then
          mass=exp(left1%value)
          found                  =.true.
       else
          if     (                                                                     &
               &   +radiusLogarithmic- left1%key < self%factorRadiusLogarithmicMaximum &
               &  .and.                                                                &
               &   -radiusLogarithmic+right1%key < self%factorRadiusLogarithmicMaximum &
               & ) then
             left2  =>  left1%predecessor()
             right2 => right1%  successor()
             if (associated(left2).and.associated(right2)) then
                massEnclosed1=(radiusLogarithmic-left1%key)*(right1%value-left1%value)/(right1%key-left1%key)+left1%value
                massEnclosed2=(radiusLogarithmic-left2%key)*(right2%value-left2%value)/(right2%key-left2%key)+left2%value
                if (Values_Agree(massEnclosed1,massEnclosed2,relTol=self%toleranceRelative)) then
                   mass=exp(massEnclosed1)
                   found                  =.true.
                end if
             end if
          end if
       end if
    end if
    if (.not.found) then
       mass=self%massDistribution_%massEnclosedBySphere(radius)
       call self%treeMassEnclosed%insert(radiusLogarithmic,log(mass))
    end if
    return
  end function sphericalAcceleratorMassEnclosedBySphere
