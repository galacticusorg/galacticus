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
Contains a module which implements a class that provides mass distributions.
!!}

module Mass_Distributions
  !!{
  Implements a class that provides mass distributions.
  !!}
  use :: Coordinates               , only : coordinate
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType     , enumerationMassTypeType, massTypeAll         , massTypeDark   , &
       &                                    massTypeBaryonic                 , massTypeGalactic       , massTypeGaseous     , massTypeStellar, &
       &                                    massTypeBlackHole                , componentTypeAll       , componentTypeUnknown, massTypeUnknown, &
       &                                    enumerationStructureErrorCodeType
  use :: Numerical_Random_Numbers  , only : randomNumberGeneratorClass
  use :: Tensors                   , only : tensorRank2Dimension3Symmetric
  private

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
   <method name="matches" >
    <description>Return true if this mass distribution matches the specified component and mass type.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type(enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
    <code>
      type(enumerationComponentTypeType) :: componentType_
      type(enumerationMassTypeType     ) :: massType_
      if (present(componentType)) then
         componentType_=componentType
      else
         componentType_=componentTypeAll
      end if
      if (present(massType     )) then
         massType_     =massType     
      else
         massType_     =massTypeAll
      end if
      massDistributionMatches= (                                                                                                                                                             &amp;
         &amp;                   massType_      == massTypeAll                                                                                                                               &amp;
         &amp;                  .or.                                                                                                                                                         &amp;
         &amp;                   massType_      == massTypeDark       .and.  self%massType == massTypeDark                                                                                   &amp;
         &amp;                  .or.                                                                                                                                                         &amp;
         &amp;                   massType_      == massTypeBaryonic   .and. (self%massType == massTypeGaseous .or. self%massType == massTypeStellar                                        ) &amp;
         &amp;                  .or.                                                                                                                                                         &amp;
         &amp;                   massType_      == massTypeGalactic   .and. (self%massType == massTypeGaseous .or. self%massType == massTypeStellar .or. self%massType == massTypeBlackHole) &amp;
         &amp;                  .or.                                                                                                                                                         &amp;
         &amp;                   massType_      == massTypeGaseous    .and.  self%massType == massTypeGaseous                                                                                &amp;
         &amp;                  .or.                                                                                                                                                         &amp;
         &amp;                   massType_      == massTypeStellar    .and.  self%massType ==                                       massTypeStellar                                          &amp;
         &amp;                  .or.                                                                                                                                                         &amp;
         &amp;                   massType_      == massTypeBlackHole  .and.  self%massType ==                                                                             massTypeBlackHole  &amp;
         &amp;                 )                                                                                                                                                             &amp;
         &amp;                .and.                                                                                                                                                          &amp;
         &amp;                 (                                                                                                                                                             &amp;
         &amp;                   componentType_ == componentTypeAll                                                                                                                          &amp;
         &amp;                  .or.                                                                                                                                                         &amp;
         &amp;                   componentType_ == self%componentType                                                                                                                        &amp;
         &amp;                 )      
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
    <argument>type(enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type(enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="acceleration" >
    <description>Return the gravitational acceleration due to the distribution at the given coordinates.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>class(coordinate                  ), intent(in   )           :: coordinates  </argument>
    <argument>type (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="tidalTensor" >
    <description>Return the gravitational tidal tensor due to the distribution at the given coordinates.</description>
    <type>type(tensorRank2Dimension3Symmetric)</type>
    <pass>yes</pass>
    <argument>class(coordinate                  ), intent(in   )           :: coordinates  </argument>
    <argument>type (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="density" >
    <description>Return the density of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate                  ), intent(in   )           :: coordinates  </argument>
    <argument>type (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="densitySphericalAverage" >
    <description>Return the average density on a spherical shell of the gievn radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                              , intent(in   )           :: radius       </argument>
    <argument>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="densityGradientRadial" >
    <description>Return the radial gradient of density of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class  (coordinate                  ), intent(in   )           :: coordinates</argument>
    <argument>logical                              , intent(in   ), optional :: logarithmic</argument>
    <argument>type   (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type   (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="potential" >
    <description>Return the gravitational potential of the distribution at the given coordinates.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate                       ), intent(in   )           :: coordinates  </argument>
    <argument>type (enumerationComponentTypeType     ), intent(in   ), optional :: componentType</argument>
    <argument>type (enumerationMassTypeType          ), intent(in   ), optional :: massType     </argument>
    <argument>type (enumerationStructureErrorCodeType), intent(  out), optional :: status       </argument>
   </method>
   <method name="massEnclosedBySphere" >
    <description>Return the mass enclosed in the distribution by a sphere of given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision                              , intent(in   )           :: radius       </argument>
    <argument>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="radiusEnclosingMass" >
    <description>Return the radius enclosing a specified mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision                              , intent(in   ), optional :: mass         , massFractional</argument>
    <argument>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType                </argument>
    <argument>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType                     </argument>
    <modules>Root_Finder</modules>
    <code>
      type            (rootFinder), save      :: finder
      logical                     , save      :: finderConstructed=.false.
      !$omp threadprivate(finder,finderConstructed)
      double precision            , parameter :: toleranceAbsolute=0.0d0  , toleranceRelative=1.0d-6

      if      (present(mass          )) then
       massTarget=     mass
      else if (present(massFractional)) then
       massTarget=self%massTotal(componentType,massType)*massFractional
      else
       call Error_Report('either "mass" or "massFractional" must be provided'//{introspection:location})
      end if
      if (massTarget &lt;= 0.0d0 .or. .not.self%matches(componentType,massType)) then
       massDistributionRadiusEnclosingMass=0.0d0
       return
      end if
      if (.not.finderConstructed) then
       finder           =rootFinder(                                                             &amp;
            &amp;                   rootFunction                 =massEnclosedRoot             , &amp;
            &amp;                   toleranceAbsolute            =toleranceAbsolute            , &amp;
            &amp;                   toleranceRelative            =toleranceRelative            , &amp;
            &amp;                   solverType                   =GSL_Root_fSolver_Brent       , &amp;
            &amp;                   rangeExpandUpward            =2.0d0                        , &amp;
            &amp;                   rangeExpandDownward          =0.5d0                        , &amp;
            &amp;                   rangeExpandType              =rangeExpandMultiplicative    , &amp;
            &amp;                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &amp;
            &amp;                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &amp;
            &amp;                  )
       finderConstructed=.true.
    end if
    self_                               =&gt; self
    massDistributionRadiusEnclosingMass =     finder%find(rootGuess=1.0d0)
    </code>
   </method>
   <method name="rotationCurve" >
    <description>Return the rotation curve at the given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                              , intent(in   )           :: radius       </argument>
    <argument>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="rotationCurveGradient" >
    <description>Return the rotation curve gradient at the given radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                              , intent(in   )           :: radius       </argument>
    <argument>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
    <argument>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="surfaceDensity" >
     <description>Return the surface density at the given coordinates.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>class(coordinate                  ), intent(in   )           :: coordinates  </argument>
     <argument>type (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
     <argument>type (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
   </method>
   <method name="surfaceDensityRadialMoment" >
     <description>Return the surface density at the given coordinates.</description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision                              , intent(in   )           :: moment                      </argument>
     <argument>double precision                              , intent(in   ), optional :: radiusMinimum, radiusMaximum</argument>
     <argument>logical                                       , intent(  out), optional :: isInfinite                  </argument>
     <argument>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType               </argument>
     <argument>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType                    </argument>
   </method>
   <method name="densityRadialMoment" >
    <description>Return the radial moment of the distribution.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                              , intent(in   )           :: moment                      </argument>
    <argument>double precision                              , intent(in   ), optional :: radiusMinimum, radiusMaximum</argument>
    <argument>logical                                       , intent(  out), optional :: isInfinite                  </argument>
    <argument>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType               </argument>
    <argument>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType                    </argument>
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
    <argument>class           (massDistributionClass       ), intent(inout)           :: massDistributionEmbedding, massDistributionPerturber</argument>
    <argument>class           (coordinate                  ), intent(in   )           :: coordinates              , velocity                 </argument>
    <argument>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType                                       </argument>
    <argument>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType                                            </argument>
    </method>
   <method name="positionSample" >
    <description>Return a position sampled from the distribution.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>class(randomNumberGeneratorClass  ), intent(inout)           :: randomNumberGenerator_</argument>
    <argument>type (enumerationComponentTypeType), intent(in   ), optional :: componentType         </argument>
    <argument>type (enumerationMassTypeType     ), intent(in   ), optional :: massType              </argument>
   </method>
   <data>class  (kinematicsDistributionClass ), pointer :: kinematicsDistribution_ => null()              </data>
   <data>logical                                        :: dimensionless                                  </data>
   <data>type   (enumerationComponentTypeType)          :: componentType           =  componentTypeUnknown</data>
   <data>type   (enumerationMassTypeType     )          :: massType                =  massTypeUnknown     </data>
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
   </method>
   <method name="velocityDispersion1D" >
    <description>Return the 1D velocity dispersion at the given coordinate.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(coordinate           ), intent(in   ) :: coordinates              </argument>
    <argument>class(massDistributionClass), intent(inout) :: massDistributionEmbedding</argument>
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

  ! Module-scope variables used in root finding.
  class           (massDistributionClass), pointer :: self_
  double precision                                 :: massTarget
  !$omp threadprivate(self_,massTarget)
  
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
  
  subroutine kinematicsDistributionAcquire(self,kinematicsDistribution_)
    !!{
    Acquire a reference to a kinematics distribution.
    !!}
    implicit none
    class(massDistributionClass      ), intent(inout)         :: self
    class(kinematicsDistributionClass), intent(in   ), target :: kinematicsDistribution_

    !![
    <referenceAcquire owner="self" target="kinematicsDistribution_" source="kinematicsDistribution_"/>
    !!]
    return
  end subroutine kinematicsDistributionAcquire
  
  double precision function massEnclosedRoot(radius)
    !!{
    Root function used in finding radii enclosing a target mass.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    massEnclosedRoot=+self_%massEnclosedBySphere(radius) &
         &           -      massTarget
    return
  end function massEnclosedRoot
  
end module Mass_Distributions
