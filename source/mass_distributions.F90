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
   <method name="positionSample" >
    <description>Return a position sampled from the distribution.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>class(randomNumberGeneratorClass  ), intent(inout)           :: randomNumberGenerator_</argument>
    <argument>type (enumerationComponentTypeType), intent(in   ), optional :: componentType         </argument>
    <argument>type (enumerationMassTypeType     ), intent(in   ), optional :: massType              </argument>
   </method>
   <data>logical                               :: dimensionless                     </data>
   <data>type   (enumerationComponentTypeType) :: componentType=componentTypeUnknown</data>
   <data>type   (enumerationMassTypeType     ) :: massType     =massTypeUnknown     </data>
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

end module Mass_Distributions
