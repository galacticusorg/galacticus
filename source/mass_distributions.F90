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
Contains a module which implements a class that provides mass distributions.
!!}

module Mass_Distributions
  !!{
  Implements a class that provides mass distributions.
  !!}
  use :: Coordinates               , only : coordinate
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType     , enumerationMassTypeType, massTypeAll           , massTypeDark                     , &
       &                                    massTypeBaryonic                 , massTypeGalactic       , massTypeGaseous       , massTypeStellar                  , &
       &                                    massTypeBlackHole                , componentTypeAll       , componentTypeUnknown  , massTypeUnknown                  , &
       &                                    componentTypeDisk                , componentTypeSpheroid  , componentTypeBlackHole, enumerationStructureErrorCodeType
  use :: Numerical_Random_Numbers  , only : randomNumberGeneratorClass
  use :: Tensors                   , only : tensorRank2Dimension3Symmetric
  use :: Numerical_Interpolation   , only : interpolator
  private
  public  :: massDistributionMatches_
  
  !![
  <functionClass>
   <name>massDistribution</name>
   <descriptiveName>Mass Distributions</descriptiveName>
   <description>Class providing mass distributions.</description>
   <destructor>
    <code>
     call kinematicsDistributionDestructor(self)
    </code>
   </destructor>
   <method name="setKinematicsDistribution">
     <description>Set the kinematics distribution for this mass distribution.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>class(kinematicsDistributionClass), intent(in   ) :: kinematicsDistribution_</argument>
     <code>
       call kinematicsDistributionAcquire(self,kinematicsDistribution_)
     </code>
   </method>
   <method name="kinematicsDistribution">
     <description>Get a pointer to the kinematics distribution for this mass distribution.</description>
     <type>class(kinematicsDistributionClass)</type>
     <pass>yes</pass>
     <code>
      massDistributionKinematicsDistribution => self%kinematicsDistribution_
      call kinematicsDistributionIncrement(self)
     </code>
   </method>
   <method name="setTypes">
     <description>Set the component and mass types of the mass distribution.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>type(enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
     <argument>type(enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
     <code>
       if (present(componentType)) self%componentType=componentType
       if (present(     massType)) self%     massType=     massType
     </code>
   </method>
   <method name="subset">
      <description>Return the subset of the mass distribution matching the given {\normalfont componentType} and {\normalfont \ttfamily massType}.</description>
      <type>class(massDistributionClass)</type>
      <pass>yes</pass>
      <argument>type(enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
      <argument>type(enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
      <code>
        if (self%matches(componentType,massType)) then
           call selfAcquire(self,massDistributionSubset)
        else
           massDistributionSubset => null()
        end if
      </code>
   </method>
   <method name="describe" >
    <description>Display a description of the mass distribution.</description>
    <type>void</type>
    <pass>yes</pass>
    <modules>Display</modules>
    <code>
     call displayMessage('unknown')
    </code>
   </method>
   <method name="matches" >
    <description>Return true if this mass distribution matches the specified component and mass type.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type(enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
    <code>
     massDistributionMatches=massDistributionMatches_(self%componentType,self%massType,componentType,massType)
    </code>
   </method>
   <method name="symmetry" >
    <description>Return the symmetry of the distribution.</description>
    <type>type(enumerationMassDistributionSymmetryType)</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      massDistributionSymmetry=massDistributionSymmetryNone
    </code>
   </method>
   <method name="isSphericallySymmetric" >
    <description>Return true if the distribution is spherically symmetric.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     massDistributionIsSphericallySymmetric=.false.
    </code>
   </method>
   <method name="assumeMonotonicDecreasingSurfaceDensity" >
    <description>Return true if the distribution can be assumed to have a monotonically decreasing surface density.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     massDistributionAssumeMonotonicDecreasingSurfaceDensity=.false.
    </code>
   </method>
   <method name="isDimensionless" >
    <description>Return true if the distribution is dimensionless.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     massDistributionIsDimensionless=self%dimensionless
    </code>
   </method>
   <method name="massTotal" >
    <description>Return the total mass of the distribution.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="acceleration" >
    <description>Return the gravitational acceleration due to the distribution at the given coordinates.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   ) :: coordinates  </argument>
   </method>
   <method name="tidalTensor" >
    <description>Return the gravitational tidal tensor due to the distribution at the given coordinates.</description>
    <type>type(tensorRank2Dimension3Symmetric)</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   ) :: coordinates  </argument>
   </method>
   <method name="density" >
    <description>Return the density of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   ) :: coordinates  </argument>
   </method>
   <method name="densitySphericalAverage" >
    <description>Return the average density on a spherical shell of the gievn radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: radius </argument>
   </method>
   <method name="densityGradientRadial" >
    <description>Return the radial gradient of density of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class  (coordinate), intent(in   )           :: coordinates</argument>
    <argument>logical            , intent(in   ), optional :: logarithmic</argument>
   </method>
   <method name="potential" >
    <description>Return the gravitational potential of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class(coordinate                       ), intent(in   )           :: coordinates  </argument>
    <argument>type (enumerationStructureErrorCodeType), intent(  out), optional :: status       </argument>
   </method>
   <method name="potentialIsAnalytic" >
     <description>Return true if the gravitational potential for this distribution has an analytic form.</description>
     <type>logical</type>
     <pass>yes</pass>
     <code>
       massDistributionPotentialIsAnalytic=.false.
     </code>
   </method>
   <method name="potentialDifference" >
    <description>Return the difference in that gravitational potential of the distribution between the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class(coordinate                       ), intent(in   )           :: coordinates1, coordinates2</argument>
    <argument>type (enumerationStructureErrorCodeType), intent(  out), optional :: status                    </argument>
    <modules>Galactic_Structure_Options</modules>
    <code>
      double precision :: potential1, potential2

      massDistributionPotentialDifference=0.0d0
      if (self%potentialIsAnalytic()) then
         potential1=self%potential(coordinates1,status)
         if (present(status)) then
            if (status /= structureErrorCodeSuccess) return
         end if
         potential2=self%potential(coordinates2,status)
         if (present(status)) then
            if (status /= structureErrorCodeSuccess) return
         end if
         massDistributionPotentialDifference=potential1-potential2
      else
         massDistributionPotentialDifference=self%potentialDifferenceNumerical(coordinates1,coordinates2,status)
      end if
    </code>
   </method>
   <method name="potentialDifferenceNumerical" >
    <description>Return the difference in that gravitational potential of the distribution between the given coordinates using a numerical calculation.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class(coordinate                       ), intent(in   )           :: coordinates1, coordinates2</argument>
    <argument>type (enumerationStructureErrorCodeType), intent(  out), optional :: status                    </argument>
    <modules>Galactic_Structure_Options</modules>
    <code>
      massDistributionPotentialDifferenceNumerical=massDistributionPotentialDifferenceNumerical_(self,coordinates1,coordinates2,status)
    </code>
   </method>
   <method name="massEnclosedBySphere" >
    <description>Return the mass enclosed in the distribution by a sphere of given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ) :: radius</argument>
   </method>
   <method name="massEnclosedByCylinder" >
    <description>Return the mass enclosed in the distribution by a cylinder of given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ) :: radius</argument>
   </method>
   <method name="radiusEnclosingMass" >
    <description>Return the radius enclosing a specified mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ), optional :: mass, massFractional</argument>
    <code>
      massDistributionRadiusEnclosingMass=self%radiusEnclosingMassNumerical(mass,massFractional)
    </code>
   </method>
   <method name="radiusEnclosingMassNumerical" >
    <description>Return the radius enclosing a specified mass using a numerical calculation.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ), optional :: mass, massFractional</argument>
    <modules>Root_Finder</modules>
    <code>
      type            (rootFinder)            :: finder
      double precision            , parameter :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6
      double precision                        :: massTarget

      if      (present(mass          )) then
       massTarget=     mass
      else if (present(massFractional)) then
       massTarget=self%massTotal()*massFractional
      else
       massTarget=0.0d0
       call Error_Report('either "mass" or "massFractional" must be provided'//{introspection:location})
      end if
      if (massTarget &lt;= 0.0d0 .or. self%massEnclosedBySphere(0.0d0) &gt;= massTarget) then
       massDistributionRadiusEnclosingMassNumerical=0.0d0
       return
      end if
      finder      =rootFinder(                                                             &amp;
            &amp;             rootFunction                 =massEnclosedRoot             , &amp;
            &amp;             toleranceAbsolute            =toleranceAbsolute            , &amp;
            &amp;             toleranceRelative            =toleranceRelative            , &amp;
            &amp;             solverType                   =GSL_Root_fSolver_Brent       , &amp;
            &amp;             rangeExpandUpward            =2.0d0                        , &amp;
            &amp;             rangeExpandDownward          =0.5d0                        , &amp;
            &amp;             rangeExpandType              =rangeExpandMultiplicative    , &amp;
            &amp;             rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &amp;
            &amp;             rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &amp;
            &amp;            )
      call self%solverSet  (massTarget=massTarget)
      massDistributionRadiusEnclosingMassNumerical=finder%find(rootGuess=1.0d0)
      call self%solverUnset(                     )
    </code>
   </method>
   <method name="radiusCylindricalEnclosingMass" >
    <description>Return the cylindrical radius enclosing a specified mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ), optional :: mass, massFractional</argument>
    <code>
      massDistributionRadiusCylindricalEnclosingMass=self%radiusCylindricalEnclosingMassNumerical(mass,massFractional)
    </code>
   </method>
   <method name="radiusCylindricalEnclosingMassNumerical" >
    <description>Return the cylindrical radius enclosing a specified mass using a numerical calculation.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ), optional :: mass, massFractional</argument>
    <modules>Root_Finder</modules>
    <code>
      type            (rootFinder)            :: finder
      double precision            , parameter :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6
      double precision                        :: massTarget

      if      (present(mass          )) then
       massTarget=     mass
      else if (present(massFractional)) then
       massTarget=self%massTotal()*massFractional
      else
       massTarget=0.0d0
       call Error_Report('either "mass" or "massFractional" must be provided'//{introspection:location})
      end if
      if (massTarget &lt;= 0.0d0 .or. self%massEnclosedByCylinder(0.0d0) &gt;= massTarget) then
       massDistributionRadiusCylindricalEnclosingMassNumerical=0.0d0
       return
      end if
      finder      =rootFinder(                                                             &amp;
            &amp;             rootFunction                 =massEnclosedCylindricalRoot  , &amp;
            &amp;             toleranceAbsolute            =toleranceAbsolute            , &amp;
            &amp;             toleranceRelative            =toleranceRelative            , &amp;
            &amp;             solverType                   =GSL_Root_fSolver_Brent       , &amp;
            &amp;             rangeExpandUpward            =2.0d0                        , &amp;
            &amp;             rangeExpandDownward          =0.5d0                        , &amp;
            &amp;             rangeExpandType              =rangeExpandMultiplicative    , &amp;
            &amp;             rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &amp;
            &amp;             rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &amp;
            &amp;            )
      call self%solverSet  (massTarget=massTarget)
      massDistributionRadiusCylindricalEnclosingMassNumerical=finder%find(rootGuess=1.0d0)
      call self%solverUnset(                     )
    </code>
   </method>
   <method name="radiusEnclosingDensity" >
    <description>Return the radius enclosing a specified density.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   )           :: density    </argument>
    <argument>double precision, intent(in   ), optional :: radiusGuess</argument>
    <code>
      massDistributionRadiusEnclosingDensity=self%radiusEnclosingDensityNumerical(density,radiusGuess)
    </code>
   </method>
   <method name="radiusEnclosingDensityNumerical" >
    <description>Return the radius enclosing a specified density using a numerical calculation.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   )           :: density    </argument>
    <argument>double precision, intent(in   ), optional :: radiusGuess</argument>
    <modules>Root_Finder</modules>
    <code>
      type            (rootFinder)            :: finder
      double precision            , parameter :: toleranceAbsolute=0.0d0  , toleranceRelative=1.0d-3
      double precision                        :: radiusGuess_

      finder    =rootFinder(                                                             &amp;
           &amp;            rootFunction                 =densityEnclosedRoot          , &amp;
           &amp;            toleranceAbsolute            =toleranceAbsolute            , &amp;
           &amp;            toleranceRelative            =toleranceRelative            , &amp;
           &amp;            solverType                   =GSL_Root_fSolver_Brent       , &amp;
           &amp;            rangeExpandUpward            =2.0d0                        , &amp;
           &amp;            rangeExpandDownward          =0.5d0                        , &amp;
           &amp;            rangeExpandType              =rangeExpandMultiplicative    , &amp;
           &amp;            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &amp;
           &amp;            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative  &amp;
           &amp;           )
      radiusGuess_                                    =     self%radiusEnclosingDensityPrevious__
      if (present(radiusGuess)) radiusGuess_=radiusGuess
      call self%solverSet  (densityTarget=density)
      massDistributionRadiusEnclosingDensityNumerical =     finder%find(rootGuess=radiusGuess_)
      call self%solverUnset(                     )
      self%radiusEnclosingDensityPrevious__           =     massDistributionRadiusEnclosingDensityNumerical
    </code>
   </method>
   <method name="radiusEnclosingSurfaceDensity" >
    <description>Return the radius enclosing a specified surface density.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   )           :: densitySurface</argument>
    <argument>double precision, intent(in   ), optional :: radiusGuess   </argument>
    <code>
      massDistributionRadiusEnclosingSurfaceDensity=self%radiusEnclosingSurfaceDensityNumerical(densitySurface,radiusGuess)
    </code>
   </method>
   <method name="radiusEnclosingSurfaceDensityNumerical" >
    <description>Return the radius enclosing a specified surface density using a numerical calculation.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   )           :: densitySurface</argument>
    <argument>double precision, intent(in   ), optional :: radiusGuess   </argument>
    <modules>Root_Finder</modules>
    <code>
      type            (rootFinder)            :: finder
      double precision            , parameter :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
      double precision                        :: radiusGuess_

      finder     =rootFinder(                                                             &amp;
           &amp;             rootFunction                 =densitySurfaceEnclosedRoot   , &amp;
           &amp;             toleranceAbsolute            =toleranceAbsolute            , &amp;
           &amp;             toleranceRelative            =toleranceRelative            , &amp;
           &amp;             solverType                   =GSL_Root_fSolver_Brent       , &amp;
           &amp;             rangeExpandUpward            =2.0d0                        , &amp;
           &amp;             rangeExpandDownward          =0.5d0                        , &amp;
           &amp;             rangeExpandType              =rangeExpandMultiplicative    , &amp;
           &amp;             rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &amp;
           &amp;             rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative  &amp;
           &amp;            )
      radiusGuess_=self%radiusEnclosingDensitySurfacePrevious__
      if (present(radiusGuess)) radiusGuess_=radiusGuess
      call self%solverSet  (densitySurfaceTarget=densitySurface)
      massDistributionRadiusEnclosingSurfaceDensityNumerical=finder%find(rootGuess=radiusGuess_)
      call self%solverUnset(                                   )
      self%radiusEnclosingDensitySurfacePrevious__          =massDistributionRadiusEnclosingSurfaceDensityNumerical
    </code>
   </method>
   <method name="radiusFromSpecificAngularMomentum" >
    <description>Return the radius corresponding to a given specific angular momentum.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ) :: angularMomentumSpecific</argument>
    <code>
      massDistributionRadiusFromSpecificAngularMomentum=self%radiusFromSpecificAngularMomentumNumerical(angularMomentumSpecific)
    </code>
   </method>
   <method name="radiusFromSpecificAngularMomentumNumerical" >
    <description>Return the radius corresponding to a given specific angular momentum using a numerical calculation.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ) :: angularMomentumSpecific</argument>
    <modules>Root_Finder</modules>
    <code>
      type            (rootFinder)            :: finder
      double precision            , parameter :: toleranceAbsolute=0.0d0  , toleranceRelative=1.0d-6

      if (angularMomentumSpecific &lt;= 0.0d0) then
         massDistributionRadiusFromSpecificAngularMomentumNumerical=+0.0d0
      else         
         finder     =rootFinder(                                                             &amp;
              &amp;             rootFunction                 =specificAngularMomentumRoot  , &amp;
              &amp;             toleranceAbsolute            =toleranceAbsolute            , &amp;
              &amp;             toleranceRelative            =toleranceRelative            , &amp;
              &amp;             solverType                   =GSL_Root_fSolver_Brent       , &amp;
              &amp;             rangeExpandUpward            =2.0d0                        , &amp;
              &amp;             rangeExpandDownward          =0.5d0                        , &amp;
              &amp;             rangeExpandType              =rangeExpandMultiplicative    , &amp;
              &amp;             rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &amp;
              &amp;             rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &amp;
              &amp;            )
         call self%solverSet  (angularMomentumSpecificTarget=angularMomentumSpecific)
         massDistributionRadiusFromSpecificAngularMomentumNumerical =     finder%find(rootGuess=1.0d0)
         call self%solverUnSet(                                                     )
      end if
    </code>
   </method>
   <method name="rotationCurve" >
    <description>Return the rotation curve at the given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: radius</argument>
   </method>
   <method name="rotationCurveGradient" >
    <description>Return the rotation curve gradient, $\mathrm{d}V^2/\mathrm{d}r$, at the given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: radius</argument>
   </method>
   <method name="velocityRotationCurveMaximum" >
    <description>Return the maximum velocity in the rotation curve.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <code>
      massDistributionVelocityRotationCurveMaximum=self%rotationCurve(self%radiusRotationCurveMaximum())
    </code>
   </method>
   <method name="radiusRotationCurveMaximum" >
    <description>Return the radius of the maximum velocity in the rotation curve.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <modules>Root_Finder</modules>
    <code>
      massDistributionRadiusRotationCurveMaximum=self%radiusRotationCurveMaximumNumerical()
    </code>
   </method>
   <method name="radiusRotationCurveMaximumNumerical" >
    <description>Return the radius of the maximum velocity in the rotation curve.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <modules>Root_Finder Error</modules>
    <code>
      type            (rootFinder)            :: finder
      double precision            , parameter :: toleranceAbsolute=0.0d0  , toleranceRelative=1.0d-06, &amp;
         &amp;                                   radiusTiny       =1.0d-9 , radiusHuge       =1.0d+30
      integer                                 :: status
      
      finder     =rootFinder(                                                              &amp;
           &amp;             rootFunction                 =rotationCurveMaximumRoot      , &amp;
           &amp;             toleranceAbsolute            =toleranceAbsolute             , &amp;
           &amp;             toleranceRelative            =toleranceRelative             , &amp;
           &amp;             solverType                   =GSL_Root_fSolver_Brent        , &amp;
           &amp;             rangeExpandUpward            =2.0d0                         , &amp;
           &amp;             rangeExpandDownward          =0.5d0                         , &amp;
           &amp;             rangeExpandType              =rangeExpandMultiplicative     , &amp;
           &amp;             rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive , &amp;
           &amp;             rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative , &amp;
           &amp;             rangeDownwardLimit           =radiusTiny                    , &amp;
           &amp;             rangeUpwardLimit             =radiusHuge                      &amp;
           &amp;            )
      call self%solverSet  ()
      massDistributionRadiusRotationCurveMaximumNumerical =     finder%find(rootGuess=1.0d0,status=status)
      call self%solverUnset()
      if (status /= errorStatusSuccess .and. .not.self%tolerateVelocityMaximumFailure) &amp;
            &amp; call Error_Report('failed to find radius of maximum circular velocity'//{introspection:location})
    </code>
   </method>
   <method name="surfaceDensity" >
     <description>Return the surface density at the given coordinates.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>class(coordinate), intent(in   ) :: coordinates  </argument>
   </method>
   <method name="surfaceDensityRadialMoment" >
     <description>Return the surface density at the given coordinates.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   )           :: moment                      </argument>
     <argument>double precision, intent(in   ), optional :: radiusMinimum, radiusMaximum</argument>
     <argument>logical         , intent(  out), optional :: isInfinite                  </argument>
   </method>
   <method name="densityRadialMoment" >
    <description>Return the radial moment of the distribution.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: moment                      </argument>
    <argument>double precision, intent(in   ), optional :: radiusMinimum, radiusMaximum</argument>
    <argument>logical         , intent(  out), optional :: isInfinite                  </argument>
   </method>
   <method name="densitySquareIntegral" >
    <description>Return the integral over the square of the density of the distribution.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), optional :: radiusMinimum, radiusMaximum</argument>
    <argument>logical         , intent(  out), optional :: isInfinite                  </argument>
   </method>
   <method name="chandrasekharIntegral" >
    <description>Return the Chandresekhar integral of the distribution.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>class           (massDistributionClass), intent(inout) :: massDistributionEmbedding, massDistributionPerturber</argument>
    <argument>double precision                       , intent(in   ) :: massPerturber                                       </argument>
    <argument>class           (coordinate           ), intent(in   ) :: coordinates              , velocity                 </argument>
   </method>
   <method name="radiusFreefall" >
    <description>Return the radius at which the freefall time to the center equals the given {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="radiusFreefallIncreaseRate" >
    <description>Return the rate of increase of the freefall radius corresponding to the given {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="fourierTransform" >
    <description>Return the spherically-symmetrized Fourier transform of the density profile at the given wavenumber.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: radiusOuter  , wavenumber</argument>
   </method>
   <method name="energy" >
    <description>Return the total energy of the distribution within the given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision                       , intent(in   )         :: radiusOuter              </argument>
    <argument>class           (massDistributionClass), intent(inout), target :: massDistributionEmbedding</argument>
   </method>
   <method name="positionSample" >
    <description>Return a position sampled from the distribution.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>class(randomNumberGeneratorClass  ), intent(inout) :: randomNumberGenerator_</argument>
   </method>
   <method name="solverSet" >
    <description>Set a sub-module scope pointers on a stack to allow recursive calls to functions.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   ), dimension(3), optional :: position1 , position2    , vectorUnit                                                     </argument>
    <argument>double precision, intent(in   )              , optional :: massTarget, densityTarget, angularMomentumSpecificTarget, densitySurfaceTarget, separation</argument>
    <code>
      integer                                        :: i
      type   (massSolver), allocatable, dimension(:) :: solvers_
      if (allocated(massSolvers)) then
         if (massSolversCount == size(massSolvers)) then
            call move_alloc(massSolvers,solvers_)
            allocate(massSolvers(size(solvers_)+massSolversIncrement))
            massSolvers(1:size(solvers_))=solvers_
            do i=1,size(solvers_)
               nullify(solvers_(i)%self)
            end do
            deallocate(solvers_)
         end if
      else
         allocate(massSolvers(massSolversIncrement))
      end if
      massSolversCount=massSolversCount+1
                                                  massSolvers(massSolversCount)%self                          => self
      if (present(separation                   )) massSolvers(massSolversCount)%separation                    =  separation
      if (present(position1                    )) massSolvers(massSolversCount)%position1                     =  position1
      if (present(position2                    )) massSolvers(massSolversCount)%position2                     =  position2
      if (present(vectorUnit                   )) massSolvers(massSolversCount)%vectorUnit                    =  vectorUnit
      if (present(massTarget                   )) massSolvers(massSolversCount)%massTarget                    =  massTarget
      if (present(densityTarget                )) massSolvers(massSolversCount)%densityTarget                 =  densityTarget
      if (present(angularMomentumSpecificTarget)) massSolvers(massSolversCount)%angularMomentumSpecificTarget =  angularMomentumSpecificTarget
      if (present(densitySurfaceTarget         )) massSolvers(massSolversCount)%densitySurfaceTarget          =  densitySurfaceTarget
    </code>
   </method>
   <method name="solverUnset" >
    <description>Unset a sub-module scope pointers on the stack.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      massSolvers(massSolversCount)%self => null()
      massSolversCount=massSolversCount-1
    </code>
   </method>
   <data>class           (kinematicsDistributionClass ), pointer :: kinematicsDistribution_                 => null()              </data>
   <data>logical                                                 :: dimensionless                           =  .false.             </data>
   <data>logical                                                 :: tolerateVelocityMaximumFailure          =  .false.             </data>
   <data>type            (enumerationComponentTypeType)          :: componentType                           =  componentTypeUnknown</data>
   <data>type            (enumerationMassTypeType     )          :: massType                                =  massTypeUnknown     </data>
   <data>double precision                                        :: radiusEnclosingDensityPrevious__        =  1.0d0               </data>
   <data>double precision                                        :: radiusEnclosingDensitySurfacePrevious__ =  1.0d0               </data>
   </functionClass>
  !!]

  !![
  <functionClass>
   <name>kinematicsDistribution</name>
   <descriptiveName>Kinematics Distributions</descriptiveName>
   <description>Class providing kinematics distributions.</description>
   <method name="isCollisional" >
    <description>Return true if the kinematics is collisional.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="temperature" >
    <description>Return the temperature of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   ) :: coordinates</argument>
    <code>
      !$GLC attributes unused :: self, coordinates
      kinematicsDistributionTemperature=0.0d0
    </code>
   </method>
   <method name="temperatureGradientLogarithmic" >
    <description>Return the logarithmic gradient of the temperature of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate), intent(in   ) :: coordinates</argument>
    <code>
      !$GLC attributes unused :: self, coordinates
      kinematicsDistributionTemperatureGradientLogarithmic=0.0d0
    </code>
   </method>
   <method name="velocityRadial" >
    <description>Return the mean radial velocity at the given coordinate.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate           ), intent(in   ) :: coordinates              </argument>
    <argument>class(massDistributionClass), intent(inout) :: massDistributionEmbedding</argument>
    <code>
      !$GLC attributes unused :: self, coordinates, massDistributionEmbedding
      kinematicsDistributionVelocityRadial=0.0d0
    </code>
   </method>
   <method name="velocityDispersion1D" >
    <description>Return the 1D velocity dispersion at the given coordinate.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate           ), intent(in   )         :: coordinates                                 </argument>
    <argument>class(massDistributionClass), intent(inout), target :: massDistribution_, massDistributionEmbedding</argument>
    <code>
      kinematicsDistributionVelocityDispersion1D=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    </code>
   </method>
   <method name="velocityDispersion1DNumerical" >
    <description>Return the 1D velocity dispersion at the given coordinate by numerically solving the Jeans equation.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate           ), intent(in   )         :: coordinates                                 </argument>
    <argument>class(massDistributionClass), intent(inout), target :: massDistribution_, massDistributionEmbedding</argument>
    <code>
      call jeansEquationSolver(self,coordinates%rSpherical(),massDistribution_,massDistributionEmbedding)
      kinematicsDistributionVelocityDispersion1DNumerical=self%velocityDispersion1D__%interpolate(log(coordinates%rSpherical()))
    </code>
   </method>
   <method name="jeansEquationRadius" >
    <description>Return the radius variable used in solving the Jeans equation that corresponds to a given physical radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                       , intent(in   ) :: radius                   </argument>
    <argument>class           (massDistributionClass), intent(inout) :: massDistributionEmbedding</argument>
    <code>
      !$GLC attributes unused :: massDistributionEmbedding
      kinematicsDistributionJeansEquationRadius=radius
    </code>
   </method>
   <method name="jeansEquationIntegrand" >
    <description>Integrand for Jeans equation.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                       , intent(in   ) :: radius                                      </argument>
    <argument>class           (massDistributionClass), intent(inout) :: massDistribution_, massDistributionEmbedding</argument>
    <modules>Numerical_Constants_Astronomical Coordinates</modules>
    <code>
      type(coordinateSpherical) :: coordinates
      if (radius > 0.0d0) then
        coordinates                                 = [radius,0.0d0,0.0d0]
        kinematicsDistributionJeansEquationIntegrand=+gravitationalConstant_internal                                 &amp;
             &amp;                                   *massDistributionEmbedding%massEnclosedBySphere(radius     )    &amp;
             &amp;                                   *massDistribution_        %density             (coordinates)    &amp;
             &amp;                                   /                                               radius      **2
      else
        kinematicsDistributionJeansEquationIntegrand=+0.0d0
      end if
    </code>
   </method>
   <method name="solverSet" >
    <description>Set a sub-module scope pointers on a stack to allow recursive calls to functions.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class           (massDistributionClass), intent(in   ), target :: massDistribution_, massDistributionEmbedding</argument>
    <code>
      integer                                              :: i
      type   (kinematicsSolver), allocatable, dimension(:) :: solvers_
      if (allocated(solvers)) then
         if (solversCount == size(solvers)) then
            call move_alloc(solvers,solvers_)
            allocate(solvers(size(solvers_)+solversIncrement))
            solvers(1:size(solvers_))=solvers_
            do i=1,size(solvers_)
               nullify(solvers_(i)%self                     )
               nullify(solvers_(i)%massDistribution_        )
               nullify(solvers_(i)%massDistributionEmbedding)
            end do
            deallocate(solvers_)
         end if
      else
         allocate(solvers(solversIncrement))
      end if
      solversCount=solversCount+1
      solvers(solversCount)%self                      => self
      solvers(solversCount)%massDistribution_         => massDistribution_
      solvers(solversCount)%massDistributionEmbedding => massDistributionEmbedding
    </code>
   </method>
   <method name="solverUnset" >
    <description>Unset a sub-module scope pointers on the stack.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      solvers(solversCount)%self                      => null()
      solvers(solversCount)%massDistribution_         => null()
      solvers(solversCount)%massDistributionEmbedding => null()
      solversCount=solversCount-1
    </code>
   </method>
   <data>type            (interpolator), allocatable               :: velocityDispersion1D__                                                                                       </data>
   <data>double precision              , allocatable, dimension(:) :: velocityDispersionRadialVelocity__                     , velocityDispersionRadialRadius__                    </data>
   <data>double precision                                          :: velocityDispersionRadialRadiusMinimum__   =+huge(0.0d0), velocityDispersionRadialRadiusMaximum__=-huge(0.0d0)</data>
   <data>double precision                                          :: velocityDispersionRadialRadiusOuter__     =+huge(0.0d0)                                                      </data>
   <data>double precision                                          :: toleranceRelativeVelocityDispersion       =1.0d-6                                                            </data>
   <data>double precision                                          :: toleranceRelativeVelocityDispersionMaximum=1.0d-3                                                            </data>
  </functionClass>
  !!]

  !![
  <functionClass>
   <name>massDistributionHeating</name>
   <descriptiveName>Heating of Mass Distributions</descriptiveName>
   <description>Class providing heating models for mass distributions.</description>
   <method name="specificEnergy" >
    <description>Return the specific energy at the given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                       , intent(in   ) :: radius           </argument>
    <argument>class           (massDistributionClass), intent(inout) :: massDistribution_</argument>
   </method>
   <method name="specificEnergyGradient" >
    <description>Return the radial gradient of the specific energy at the given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                       , intent(in   ) :: radius           </argument>
    <argument>class           (massDistributionClass), intent(inout) :: massDistribution_</argument>
   </method>
   <method name="specificEnergyIsEverywhereZero" >
    <description>Return true if the specific energy is zero everywhere (i.e. no heating).</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

  ! Enumeration of mass distribution symmetries.
  !![
  <enumeration>
   <name>massDistributionSymmetry</name>
   <description>Specifies the symmetry of {\normalfont \ttfamily massDistribution} objects.</description>
   <visibility>public</visibility>
   <entry label="none"        />
   <entry label="cylindrical" />
   <entry label="spherical"   />
  </enumeration>
  !!]

  ! Enumeration of non-analytic solver options.
  !![
  <enumeration>
   <name>nonAnalyticSolvers</name>
   <description>Used to specify the type of solution to use when no analytic solution is available.</description>
   <encodeFunction>yes</encodeFunction>
   <visibility>public</visibility>
   <validator>yes</validator>
   <entry label="fallThrough"/>
   <entry label="numerical"  />
  </enumeration>
  !!]

  ! Module-scope pointers used in integrand functions and root finding.
  type :: kinematicsSolver
     class(kinematicsDistributionClass), pointer :: self              => null()
     class(massDistributionClass      ), pointer :: massDistribution_ => null(), massDistributionEmbedding => null()
  end type kinematicsSolver
  type   (kinematicsSolver), allocatable, dimension(:) :: solvers
  integer                  , parameter                 :: solversIncrement=10
  integer                                              :: solversCount    = 0
  !$omp threadprivate(solvers,solversCount)
  
  ! Module-scope pointers used in integrand functions and root finding.
  type :: massSolver
     class           (massDistributionClass), pointer      :: self                          => null()
     double precision                       , dimension(3) :: position1                               , position2           , &
          &                                                   vectorUnit
     double precision                                      :: massTarget                              , densityTarget       , &
          &                                                   angularMomentumSpecificTarget           , densitySurfaceTarget, &
          &                                                   separation
  end type massSolver
  type   (massSolver), allocatable, dimension(:) :: massSolvers
  integer            , parameter                 :: massSolversIncrement=10
  integer                                        :: massSolversCount    = 0
  !$omp threadprivate(massSolvers,massSolversCount)
  
contains
  
  subroutine kinematicsDistributionDestructor(self)
    !!{
    Destroy a kinematics distribution.
    !!}
    implicit none
    type(massDistributionClass ), intent(inout) :: self

    !![
    <objectDestructor name="self%kinematicsDistribution_"/>
    !!]
    return
  end subroutine kinematicsDistributionDestructor
  
  subroutine selfAcquire(self,self_)
    !!{
    Acquire a reference to a mass distribution.
    !!}
    implicit none
    class(massDistributionClass), intent(inout), target  :: self
    class(massDistributionClass), intent(  out), pointer :: self_

    !![
    <referenceAcquire target="self_" source="self"/>
    !!]
    return
  end subroutine selfAcquire

  subroutine kinematicsDistributionAcquire(self,kinematicsDistribution_)
    !!{
    Acquire a reference to a kinematics distribution.
    !!}
    implicit none
    class(massDistributionClass      ), intent(inout)         :: self
    class(kinematicsDistributionClass), intent(in   ), target :: kinematicsDistribution_

    !![
    <objectDestructor name="self%kinematicsDistribution_"/>
    <referenceAcquire owner="self" target="kinematicsDistribution_" source="kinematicsDistribution_"/>
    !!]
    return
  end subroutine kinematicsDistributionAcquire
  
  subroutine kinematicsDistributionIncrement(self)
    !!{
    Increment the reference count to a kinematics distribution.
    !!}
    implicit none
    class(massDistributionClass), intent(inout) :: self

    !![
    <referenceCountIncrement owner="self" object="kinematicsDistribution_"/>
    !!]
    return
  end subroutine kinematicsDistributionIncrement

  logical function massDistributionMatches_(componentTypeTarget,massTypeTarget,componentType,massType)
    !!{
    Determine if the requested mass and component types match that of the target mass distribution.
    !!}
    implicit none
    type(enumerationComponentTypeType), intent(in   )           :: componentTypeTarget
    type(enumerationMassTypeType     ), intent(in   )           :: massTypeTarget
    type(enumerationComponentTypeType), intent(in   ), optional :: componentType
    type(enumerationMassTypeType     ), intent(in   ), optional :: massType
    
    if (present(componentType)) then
       massDistributionMatches_   = componentType == componentTypeAll                                                                                                                                                           &
            &                      .or.                                                                                                                                                                                         &
            &                       componentType == componentTypeTarget
       if (massDistributionMatches_.and.present(massType)) then
          massDistributionMatches_= (massType      == massTypeAll                                                                                                                                                             ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeDark       .and.  massTypeTarget      == massTypeDark                                                                                                           ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeGaseous    .and.  massTypeTarget      == massTypeGaseous                                                                                                        ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeStellar    .and.                                                massTypeTarget      == massTypeStellar                                                          ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeBlackHole  .and.                                                                                                  massTypeTarget      == massTypeBlackHole      ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeBaryonic   .and. (massTypeTarget      == massTypeGaseous   .or. massTypeTarget      == massTypeStellar                                                         )) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeGalactic   .and. (massTypeTarget      == massTypeGaseous   .or. massTypeTarget      == massTypeStellar       .or. massTypeTarget      == massTypeBlackHole      ) &
               &                                                         .and. (componentTypeTarget == componentTypeDisk .or. componentTypeTarget == componentTypeSpheroid .or. componentTypeTarget == componentTypeBlackHole)) 
       end if
    else
       if (present(massType)) then
          ! Match on mass distribution only.
          massDistributionMatches_= (massType      == massTypeAll                                                                                                                                                             ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeDark       .and.  massTypeTarget      == massTypeDark                                                                                                           ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeGaseous    .and.  massTypeTarget      == massTypeGaseous                                                                                                        ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeStellar    .and.                                                massTypeTarget      == massTypeStellar                                                          ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeBlackHole  .and.                                                                                                  massTypeTarget      == massTypeBlackHole      ) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeBaryonic   .and. (massTypeTarget      == massTypeGaseous   .or. massTypeTarget      == massTypeStellar                                                         )) &
               &                   .or.                                                                                                                                                                                         &
               &                    (massType      == massTypeGalactic   .and. (massTypeTarget      == massTypeGaseous   .or. massTypeTarget      == massTypeStellar       .or. massTypeTarget      == massTypeBlackHole      ) &
               &                                                         .and. (componentTypeTarget == componentTypeDisk .or. componentTypeTarget == componentTypeSpheroid .or. componentTypeTarget == componentTypeBlackHole))
       else
          ! No selection was requested, so we always match.
          massDistributionMatches_=.true.
       end if
    end if
    return
  end function massDistributionMatches_
  
  double precision function massEnclosedRoot(radius)
    !!{
    Root function used in finding radii enclosing a target mass.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    massEnclosedRoot=+massSolvers(massSolversCount)%self%massEnclosedBySphere(radius) &
         &           -massSolvers(massSolversCount)     %massTarget
    return
  end function massEnclosedRoot
  
  double precision function massEnclosedCylindricalRoot(radius)
    !!{
    Root function used in finding cylindrical radii enclosing a target mass.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    massEnclosedCylindricalRoot=+massSolvers(massSolversCount)%self%massEnclosedByCylinder(radius) &
         &                      -massSolvers(massSolversCount)     %massTarget
    return
  end function massEnclosedCylindricalRoot
  
  double precision function densityEnclosedRoot(radius)
    !!{
    Root function used in finding radii enclosing a target density.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    densityEnclosedRoot=+3.0d0                                                              &
         &              /4.0d0                                                              &
         &              /Pi                                                                 &
         &              *massSolvers(massSolversCount)%self%massEnclosedBySphere(radius)    &
         &              /                                   radius                      **3 &
         &              -massSolvers(massSolversCount)%     densityTarget
    return
  end function densityEnclosedRoot
  
  double precision function densitySurfaceEnclosedRoot(radius)
    !!{
    Root function used in finding radii enclosing a target surface density.
    !!}
    use :: Coordinates, only : coordinateCylindrical, assignment(=)
    implicit none
    double precision                       , intent(in   ) :: radius
    type            (coordinateCylindrical)                :: coordinates

    coordinates               =[radius,0.0d0,0.0d0]
    densitySurfaceEnclosedRoot=+massSolvers(massSolversCount)%self%surfaceDensity      (coordinates) &
         &                     -massSolvers(massSolversCount)     %densitySurfaceTarget
    return
  end function densitySurfaceEnclosedRoot
  
  double precision function specificAngularMomentumRoot(radius)
    !!{
    Root function used in finding radii corresponding to a target specific angular momentum.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    specificAngularMomentumRoot=+massSolvers(massSolversCount)%self%rotationCurve                (radius) &
         &                      *                                   radius                                &
         &                      -massSolvers(massSolversCount)     %angularMomentumSpecificTarget
    return
  end function specificAngularMomentumRoot
  
  double precision function rotationCurveMaximumRoot(radius)
    !!{
    Root function used in finding the radius corresponding to the peak of the rotation curve.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    rotationCurveMaximumRoot=massSolvers(massSolversCount)%self%rotationCurveGradient(radius)
    return
  end function rotationCurveMaximumRoot

  double precision function massDistributionPotentialDifferenceNumerical_(self,coordinates1,coordinates2,status) result(potentialDifference)
    !!{
    Numerically calculate the potential difference between the provided coordinates.
    !!}
    use :: Coordinates               , only : coordinateCartesian      , assignment(=)
    use :: Vectors                   , only : Vector_Magnitude
    use :: Numerical_Integration     , only : integrator
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess, structureErrorCodeIntegration
    use :: Error                     , only : Error_Report             , errorStatusSuccess
    implicit none
    class           (massDistributionClass            ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates1 , coordinates2
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    type            (coordinateCartesian              )                          :: coordinates1_, coordinates2_
    double precision                                   , dimension(3)            :: position1    , position2    , &
         &                                                                          vectorUnit
    double precision                                                             :: separation
    type            (integrator                       )                          :: integrator_
    integer                                                                      :: status_
    
    if (present(status)) status=structureErrorCodeSuccess
    coordinates1_=coordinates1
    coordinates2_=coordinates2
    position1    =coordinates1_
    position2    =coordinates2_
    separation   =Vector_Magnitude(position1-position2)
    if (separation == 0.0d0) then
       potentialDifference=0.0d0
    else
       vectorUnit   =+(            &
            &          +position1  &
            &          -position2  &
            &         )            &
            &        /  separation
       integrator_=integrator(potentialDifferenceIntegrand,toleranceRelative=1.0d-3)
       call self%solverSet  (position1=position1,position2=position2,vectorUnit=vectorUnit,separation=separation)
       potentialDifference=integrator_%integrate(0.0d0,separation,status=status_)
       call self%solverUnset(                                         )
       if (status_ /= errorStatusSuccess) then
          if (present(status)) then
             status=structureErrorCodeIntegration
          else
             call Error_Report("integration of potential difference failed"//{introspection:location})
          end if
       end if
    end if
    return
  end function massDistributionPotentialDifferenceNumerical_

  double precision function potentialDifferenceIntegrand(distance)
    !!{
    Integrand used in computing potential differences.
    !!}
    use :: Coordinates                     , only : coordinateCartesian, assignment(=)
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    double precision                     , intent(in   ) :: distance
    double precision                     , dimension(3)  :: position   , acceleration
    type            (coordinateCartesian)                :: coordinates
    
    position                    =+massSolvers(massSolversCount)    %position2                &
         &                       +massSolvers(massSolversCount)    %vectorUnit               &
         &                       *                                  distance
    coordinates                 =                                   position
    acceleration                =massSolvers(massSolversCount)%self%acceleration(coordinates)
    potentialDifferenceIntegrand=-Dot_Product(acceleration,massSolvers(massSolversCount)%vectorUnit) &
         &                       *MpcPerKmPerSToGyr
    return
  end function potentialDifferenceIntegrand
  
  subroutine jeansEquationSolver(self,radius,massDistribution_,massDistributionEmbedding)
    !!{
    Solve the Jeans equation numerically to find the 1D velocity dispersion.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Display                , only : displayIndent       , displayUnindent     , displayMessage
    use            :: Error                  , only : Error_Report        , errorStatusSuccess  , GSL_Error_Details
    use            :: Numerical_Integration  , only : integrator
    use            :: Numerical_Ranges       , only : Make_Range          , rangeTypeLogarithmic
    use            :: Table_Labels           , only : extrapolationTypeFix
    use            :: Numerical_Interpolation, only : gsl_interp_linear
    use            :: String_Handling        , only : operator(//)
    use            :: Coordinates            , only : coordinateSpherical , assignment(=)
    implicit none
    class           (kinematicsDistributionClass), intent(inout)              :: self
    double precision                             , intent(in   )              :: radius
    class           (massDistributionClass      ), intent(inout), target      :: massDistribution_               , massDistributionEmbedding
    double precision                                            , parameter   :: radiusTinyFactor        =1.0d-9 , factorDensityLarge       =1.0d+5
    double precision                                            , parameter   :: countPointsPerOctave    =2.0d0
    double precision                                            , parameter   :: toleranceFactor         =2.0d0
    double precision                             , dimension(:) , allocatable :: velocityDispersions             , radii
    double precision                                                          :: radiusMinimum                   , radiusMaximum                   , &
         &                                                                       toleranceRelative               , density                         , &
         &                                                                       jeansIntegral                   , radiusOuter_                    , &
         &                                                                       radiusLower                     , radiusUpper                     , &
         &                                                                       radiusLowerJeansEquation        , radiusUpperJeansEquation        , &
         &                                                                       densityMaximum                  , densityOuter_                   , &
         &                                                                       jeansIntegralPrevious           , toleranceRelativePrevious
    integer         (c_size_t                   )                             :: countRadii                      , iMinimum                        , &
         &                                                                       iMaximum                        , i
    integer                                                                   :: status
    type            (coordinateSpherical        )                             :: coordinates
    type            (integrator                 )                             :: integrator_
    logical                                                                   :: remakeTable

    ! Determine if the table must be rebuilt.
    remakeTable=.false.
    if (.not.allocated(self%velocityDispersionRadialVelocity__)) then
       remakeTable=.true.
    else
       remakeTable= radius < self%velocityDispersionRadialRadiusMinimum__ &
            &      .or.                                                   &
            &       radius > self%velocityDispersionRadialRadiusMaximum__
    end if
    if (remakeTable) then
       integrator_=integrator(jeansEquationIntegrand_,toleranceRelative=self%toleranceRelativeVelocityDispersion,intervalsMaximum=10000_c_size_t)
       toleranceRelativePrevious=self%toleranceRelativeVelocityDispersion
       ! Find the range of radii at which to compute the velocity dispersion, and construct the arrays.
       call self%solverSet(massDistribution_,massDistributionEmbedding)
       !! Set an initial range of radii that brackets the requested radius.
       radiusMinimum=0.5d0*radius
       radiusMaximum=2.0d0*radius
       !! Round to the nearest factor of 2.
       radiusMinimum=2.0d0**floor  (log(radiusMinimum)/log(2.0d0))
       radiusMaximum=2.0d0**ceiling(log(radiusMaximum)/log(2.0d0))
       !! Expand to encompass any pre-existing range.
       if (allocated(self%velocityDispersionRadialRadius__)) then
          radiusMinimum=min(radiusMinimum,self%velocityDispersionRadialRadiusMinimum__)
          radiusMaximum=max(radiusMaximum,self%velocityDispersionRadialRadiusMaximum__)
       end if
       !! Set a suitable outer radius for integration. We require that the density at the outer radius be much smaller than that
       !! at the maximum radius. This should ensure that any constribution to the Jeans integral from beyond this outer radius is
       !! negligible.
       !!! Start at the maximum radius and gradually increase the outer radius until the density is sufficiently small.
       coordinates    =[radiusMaximum,0.0d0,0.0d0]
       densityMaximum=massDistribution_%density(coordinates)
       radiusOuter_  =radiusMaximum
       densityOuter_ =densityMaximum
       do while (densityOuter_ > densityMaximum/factorDensityLarge)
          radiusOuter_ =radiusOuter_*2.0d0
          coordinates  =[radiusOuter_,0.0d0,0.0d0]
          densityOuter_=massDistribution_%density(coordinates)
       end do
       !! Construct arrays.
       countRadii=nint(log(radiusMaximum/radiusMinimum)/log(2.0d0)*countPointsPerOctave+1.0d0)
       allocate(radii              (countRadii))
       allocate(velocityDispersions(countRadii))
       radii=Make_Range(radiusMinimum,radiusMaximum,int(countRadii),rangeTypeLogarithmic)
       ! Copy in any usable results from any previous solution.
       !! Assume by default that no previous solutions are usable.
       iMinimum=+huge(0_c_size_t)
       iMaximum=-huge(0_c_size_t)
       !! Check that a pre-existing solution exists.
       if (allocated(self%velocityDispersionRadialRadius__)) then
          !! Check that the outer radius for integration has not changed - if it has we need to recompute the full solution for
          !! consistency.
          if (radiusOuter_ == self%velocityDispersionRadialRadiusOuter__) then
             iMinimum=nint(log(self%velocityDispersionRadialRadiusMinimum__/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
             iMaximum=nint(log(self%velocityDispersionRadialRadiusMaximum__/radiusMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
             velocityDispersions(iMinimum:iMaximum)=self%velocityDispersionRadialVelocity__
          end if
       end if
       ! Solve for the velocity dispersion where old results were unavailable.
       jeansIntegralPrevious=0.0d0
       do i=countRadii,1,-1
          ! Skip cases for which we have a pre-existing solution.
          if (i >= iMinimum .and. i <= iMaximum) cycle
          ! Find the limits for the integral.
          if (i == countRadii) then
             radiusUpper=radiusOuter_
          else
             radiusUpper=radii(i+1)
          end if
          radiusLower   =radii(i  )
          ! Reset the accumulated Jeans integral if necessary.
          if (i == iMinimum-1) then
             coordinates          = [radii(iMinimum),0.0d0,0.0d0]
             jeansIntegralPrevious=+                  velocityDispersions(iMinimum   )**2 &
                  &                *massDistribution_%density            (coordinates)
          end if
          ! If the interval is wholly outside of the outer radius, the integral is zero.
          if (radiusLower > radiusOuter_) then
             jeansIntegral         =0.0d0
             velocityDispersions(i)=0.0d0
          else
             ! Evaluate the integral.
             coordinates             =[radiusLower,0.0d0,0.0d0]
             density                 =massDistribution_%density            (coordinates                                                                )
             radiusLowerJeansEquation=self             %jeansEquationRadius(radiusLower                                      ,massDistributionEmbedding)
             radiusUpperJeansEquation=self             %jeansEquationRadius(radiusUpper                                      ,massDistributionEmbedding)
             jeansIntegral           =integrator_      %integrate          (radiusLowerJeansEquation,radiusUpperJeansEquation,status                   )
             if (status /= errorStatusSuccess) then
                ! Integration failed.
                toleranceRelative=+     toleranceFactor                     &
                     &            *self%toleranceRelativeVelocityDispersion
                do while (toleranceRelative <= self%toleranceRelativeVelocityDispersionMaximum)
                   call integrator_%toleranceSet(toleranceRelative=toleranceRelative)
                   toleranceRelativePrevious=toleranceRelative
                   jeansIntegral=integrator_%integrate(radiusLowerJeansEquation,radiusUpperJeansEquation,status)
                   if (status == errorStatusSuccess) then
                      exit
                   else
                      if (toleranceRelative >= self%toleranceRelativeVelocityDispersionMaximum) then
                         exit
                      else
                         toleranceRelative=min(                                                  &
                              &                +     toleranceFactor                             &
                              &                *     toleranceRelative                         , &
                              &                +self%toleranceRelativeVelocityDispersionMaximum  &
                              &               )
                      end if
                   end if
                end do
                if (status /= errorStatusSuccess) then
                   block
                     integer                   :: line
                     type     (varying_string) :: reason, file
                     character(len=24        ) :: label
                     call GSL_Error_Details(reason,file,line,status)
                     call integrator_%toleranceSet(toleranceRelative=toleranceRelativePrevious)
                     call displayIndent('Jeans equation integration failure report')
                     call displayMessage(var_str('radius ')//i//' of '//countRadii)
                     write (label,'(e24.16)') radiusLowerJeansEquation
                     call displayMessage('radiusLowerJeansEquation = '//trim(label)//' Mpc'             )
                     write (label,'(e24.16)') radiusUpperJeansEquation
                     call displayMessage('radiusUpperJeansEquation = '//trim(label)//' Mpc'             )
                     write (label,'(e24.16)') radiusOuter_
                     call displayMessage('radiusOuter              = '//trim(label)//' Mpc'             )
                     write (label,'(e24.16)') toleranceRelative
                     call displayMessage('toleranceRelative        = '//trim(label)                     )
                     write (label,'(e24.16)') jeansIntegral
                     call displayMessage('jeansIntegral            = '//trim(label)//' M☉ Mpc⁻³ km² s⁻²')
                     call displayIndent('Mass distribution properties:')
                     call displayMessage('Type: '//massDistribution_%objectType())
                     call massDistribution_        %describe()
                     call displayUnindent('done')
                     call displayIndent('Mass distribution embedding properties:')
                     call displayMessage('Type: '//massDistributionEmbedding%objectType())
                     call massDistributionEmbedding%describe()
                     call displayUnindent('done')
                     call displayUnindent('done')
                     call Error_Report(var_str('integration of Jeans equation failed - GSL Error ')//status//': "'//reason//'" at line '//line//' of file "'//file//'"'//{introspection:location})
                   end block
                end if
                call integrator_%toleranceSet(toleranceRelative=self%toleranceRelativeVelocityDispersion)
                toleranceRelativePrevious=self%toleranceRelativeVelocityDispersion
             end if
             if (density <= 0.0d0) then
                ! Density is zero - the velocity dispersion is undefined. If the Jeans integral is also zero this is acceptable - we've
                ! been asked for the velocity dispersion in a region of zero density, so we simply return zero dispersion as it should have
                ! no consequence. If the Jeans integral is non-zero however, then something has gone wrong.
                velocityDispersions(i)=0.0d0
                if (jeansIntegral+jeansIntegralPrevious > 0.0d0) call Error_Report('undefined velocity dispersion'//{introspection:location})
             else
                velocityDispersions(i)=sqrt(                         &
                     &                      +(                       &
                     &                        +jeansIntegral         &
                     &                        +jeansIntegralPrevious &
                     &                       )                       &
                     &                      /density                 &
                     &                     )
             end if
          end if
          jeansIntegralPrevious=+jeansIntegralPrevious &
               &                +jeansIntegral
       end do
       call self%solverUnset()
       ! Build the interpolator.
       if (allocated(self%velocityDispersion1D__)) deallocate(self%velocityDispersion1D__)
       allocate(self%velocityDispersion1D__)
       self%velocityDispersion1D__=interpolator(log(radii),velocityDispersions,interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeFix)
       ! Store the current results for future re-use.
       if (allocated(self%velocityDispersionRadialRadius__  )) deallocate(self%velocityDispersionRadialRadius__  )
       if (allocated(self%velocityDispersionRadialVelocity__)) deallocate(self%velocityDispersionRadialVelocity__)
       allocate(self%velocityDispersionRadialRadius__  (countRadii))
       allocate(self%velocityDispersionRadialVelocity__(countRadii))
       self%velocityDispersionRadialRadius__       =radii
       self%velocityDispersionRadialVelocity__     =velocityDispersions
       self%velocityDispersionRadialRadiusMinimum__=radiusMinimum
       self%velocityDispersionRadialRadiusMaximum__=radiusMaximum
       self%velocityDispersionRadialRadiusOuter__  =radiusOuter_
    end if
    return
  end subroutine jeansEquationSolver

  double precision function jeansEquationIntegrand_(radius)
    !!{
    Integrand for the Jeans equation.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    jeansEquationIntegrand_=solvers(solversCount)%self%jeansEquationIntegrand(radius,solvers(solversCount)%massDistribution_,solvers(solversCount)%massDistributionEmbedding)
    return
  end function jeansEquationIntegrand_
  
end module Mass_Distributions
