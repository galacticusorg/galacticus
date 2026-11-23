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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a property extractor class for the density at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale, darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusDefinitions

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDensityProfile" docformat="rst">
   <description>
   A property extractor that returns the mass density profile (in :math:`\mathrm{M}_\odot \, \mathrm{Mpc}^{-3}`) of a galaxy or
   halo component at a user-specified set of radii. The ``radiusSpecifiers`` parameter provides a list of radius definitions
   (e.g. multiples of the virial radius, disk radius, or half-mass radius), supporting both galactic structural radii and fixed
   physical radii. If ``includeRadii`` is ``true``, the corresponding radii (in Mpc) are also written to the output file as a
   second column alongside the density values.

   A radius specifier can evaluate to zero---most often because the component on which it is based is absent or empty in the node
   in question. The density is then finite only if the mass distribution being evaluated has a finite central density; where it
   does not, the requested property does not exist. By default this is reported as a fatal error naming the offending specifier;
   setting ``zeroRadiusIsFatal`` to ``false`` instead writes the undefined-value sentinel in place of the density.

   By default the density is evaluated in the equatorial plane of the component (polar angle :math:`\theta=90^\circ`). If
   ``polarAngles`` is given, it is instead evaluated at each radius at each of the listed polar angles, measured from the
   component's symmetry axis, so that the profile of an axisymmetric component (such as a disk) is tabulated on a grid in
   :math:`(r,\theta)`. The output then has one column per (radius, angle) pair, with the radius varying fastest, and the polar
   angle of each column is written as its column value.
   </description>
   <stateStorable>
    <!-- Not stored: radii - radius definitions, set by the constructor and never changed afterwards. -->
    <exclude variables="radii"/>
   </stateStorable>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorDensityProfile
     !!{RST
     A property extractor class for the density at a set of radii.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_ => null()
     integer                                                               :: radiiCount                    , elementCount_    , &
          &                                                                   anglesCount
     logical                                                               :: includeRadii                  , zeroRadiusIsFatal
     type            (varying_string          ), allocatable, dimension(:) :: radiusSpecifiers
     type            (radiusDefinitions       )                            :: radii
     double precision                          , allocatable, dimension(:) :: polarAngles
   contains
     final     ::                       densityProfileDestructor
     procedure :: columnDescriptions => densityProfileColumnDescriptions
     procedure :: size               => densityProfileSize
     procedure :: elementCount       => densityProfileElementCount
     procedure :: extract            => densityProfileExtract
     procedure :: names              => densityProfileNames
     procedure :: descriptions       => densityProfileDescriptions
     procedure :: unitsInSI          => densityProfileUnitsInSI
     procedure :: units              => densityProfileUnits
  end type nodePropertyExtractorDensityProfile

  interface nodePropertyExtractorDensityProfile
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorDensityProfile` property extractor class.
     !!}
     module procedure densityProfileConstructorParameters
     module procedure densityProfileConstructorInternal
  end interface nodePropertyExtractorDensityProfile

contains

  function densityProfileConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorDensityProfile` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorDensityProfile)                              :: self
    type            (inputParameters                    ), intent(inout)               :: parameters
    type            (varying_string                     ), allocatable  , dimension(:) :: radiusSpecifiers
    class           (darkMatterHaloScaleClass           ), pointer                     :: darkMatterHaloScale_
    logical                                                                            :: includeRadii        , zeroRadiusIsFatal
    double precision                                     , allocatable  , dimension(:) :: polarAngles

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !![
    <inputParameter docformat="rst">
      <name>radiusSpecifiers</name>
      <description>
      A list of radius specifiers at which to output the density profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>includeRadii</name>
      <defaultValue>.false.</defaultValue>
      <description>
      Specifies whether or not the radii at which density data are output should also be included in the output file.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>zeroRadiusIsFatal</name>
      <defaultValue>.true.</defaultValue>
      <description>
      Specifies whether a radius specifier which evaluates to zero, in a node in which the density diverges at zero radius,
      should be reported as a fatal error. If ``false``, the undefined-value sentinel is written in place of the density
      instead.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    if (parameters%isPresent('polarAngles')) then
       allocate(polarAngles(parameters%count('polarAngles')))
       !![
       <inputParameter docformat="rst">
         <name>polarAngles</name>
         <description>
         A list of polar angles (in degrees, measured from the symmetry axis of the component, so that :math:`90^\circ` is the
         equatorial plane) at which to evaluate the density at each radius. If absent, the density is evaluated in the
         equatorial plane only.
         </description>
         <source>parameters</source>
       </inputParameter>
       !!]
    else
       allocate(polarAngles(1))
       polarAngles=90.0d0
    end if
    self=nodePropertyExtractorDensityProfile(radiusSpecifiers,polarAngles,includeRadii,zeroRadiusIsFatal,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function densityProfileConstructorParameters

  function densityProfileConstructorInternal(radiusSpecifiers,polarAngles,includeRadii,zeroRadiusIsFatal,darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorDensityProfile` property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (nodePropertyExtractorDensityProfile)                              :: self
    type            (varying_string                     ), intent(in   ), dimension(:) :: radiusSpecifiers
    double precision                                     , intent(in   ), dimension(:) :: polarAngles
    class           (darkMatterHaloScaleClass           ), intent(in   ), target       :: darkMatterHaloScale_
    logical                                              , intent(in   )               :: includeRadii        , zeroRadiusIsFatal
    !![
    <constructorAssign variables="radiusSpecifiers, polarAngles, includeRadii, zeroRadiusIsFatal, *darkMatterHaloScale_"/>
    !!]

    if (includeRadii) then
       self%elementCount_=2
    else
       self%elementCount_=1
    end if
    self%radiiCount =size(radiusSpecifiers)
    self%anglesCount=size(polarAngles     )
    if (self%anglesCount < 1) call Error_Report('at least one polar angle must be given'//{introspection:location})
    if (any(polarAngles < 0.0d0 .or. polarAngles > 180.0d0)) &
         & call Error_Report('polar angles must lie in the range [0,180] degrees'//{introspection:location})
    call self%radii%decode(radiusSpecifiers)
    return
  end function densityProfileConstructorInternal

  subroutine densityProfileDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorDensityProfile` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDensityProfile), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine densityProfileDestructor

  integer function densityProfileElementCount(self,time)
    !!{RST
    Return the number of elements in the ``densityProfile`` property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityProfile), intent(inout) :: self
    double precision                                     , intent(in   ) :: time
    !$GLC attributes unused :: time

    densityProfileElementCount=self%elementCount_
    return
  end function densityProfileElementCount

  function densityProfileSize(self,time)
    !!{RST
    Return the number of array elements in the ``densityProfile`` property extractors.
    !!}
    implicit none
    integer         (c_size_t                           )                :: densityProfileSize
    class           (nodePropertyExtractorDensityProfile), intent(inout) :: self
    double precision                                     , intent(in   ) :: time
    !$GLC attributes unused :: time

    densityProfileSize=self%radiiCount*self%anglesCount
    return
  end function densityProfileSize

  function densityProfileExtract(self,node,time,instance)
    !!{RST
    Implement a ``densityProfile`` property extractor.
    !!}
    use :: Galactic_Structure_Radii_Definitions, only : radiusResolver        , radiusUndefined
    use :: Galacticus_Nodes                    , only : treeNode
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Coordinates                         , only : coordinateSpherical  , assignment(=)
    use :: Numerical_Constants_Math            , only : Pi
    implicit none
    double precision                                     , dimension(:,:), allocatable :: densityProfileExtract
    class           (nodePropertyExtractorDensityProfile), intent(inout) , target      :: self
    type            (treeNode                           ), intent(inout) , target      :: node
    double precision                                     , intent(in   )               :: time
    type            (multiCounter                       ), intent(inout) , optional    :: instance
    class           (massDistributionClass              ), pointer                     :: massDistribution_
    type            (radiusResolver                     )                              :: resolver
    type            (coordinateSpherical                )                              :: coordinates
    integer                                                                            :: i                    , j, &
         &                                                                                k
    double precision                                                                   :: radius
    !$GLC attributes unused :: time, instance

    allocate(densityProfileExtract(self%radiiCount*self%anglesCount,self%elementCount_))
    resolver=radiusResolver(self%radii,node,self%darkMatterHaloScale_)
    do i=1,self%radiiCount
       ! Columns are ordered with the radius varying fastest, so the column for radius i at angle j is i+(j-1)*radiiCount.
       call resolver%evaluate(i,radius)
       if (radius < 0.0d0) then
          ! The radius is undefined in this node - report the sentinel rather than evaluating the density there.
          do j=1,self%anglesCount
             k=i+(j-1)*self%radiiCount
             densityProfileExtract       (k,1)=radiusUndefined
             if (self%includeRadii)                            &
                  & densityProfileExtract(k,2)=radiusUndefined
          end do
          cycle
       end if
       massDistribution_ => node%massDistribution(                                                  &
            &                                     componentType=self%radii%specifiers(i)%component, &
            &                                     massType     =self%radii%specifiers(i)%mass       &
            &                                    )
       if     (                                                             &
            &   radius                                             <= 0.0d0 &
            &  .and.                                                        &
            &   massDistribution_%densitySlopeLogarithmicCentral() <  0.0d0 &
            & ) then
          ! The radius is zero, but the density of this mass distribution diverges at the center (or its central slope is
          ! unknown), so the density can not be evaluated here.
          if (self%zeroRadiusIsFatal)                                                                                               &
               & call resolver%reportZeroRadius(i,'densityProfile','the density of this mass distribution diverges at zero radius')
          do j=1,self%anglesCount
             k=i+(j-1)*self%radiiCount
             densityProfileExtract       (k,1)=radiusUndefined
             if (self%includeRadii)                            &
                  & densityProfileExtract(k,2)=radius
          end do
          !![
          <objectDestructor name="massDistribution_"/>
          !!]
          cycle
       end if
       do j=1,self%anglesCount
          k                                 =  i+(j-1)*self%radiiCount
          ! Written as Pi*(angle/180) so that the default 90° gives exactly π/2, as before polarAngles existed.
          coordinates                       =  [radius,Pi*(self%polarAngles(j)/180.0d0),0.0d0]
          densityProfileExtract       (k,1) =  massDistribution_%density(                        &
               &                                                         coordinates=coordinates &
               &                                                        )
          if (self%includeRadii)                                                                 &
               & densityProfileExtract(k,2) =                                        radius
       end do
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    end do
    return
  end function densityProfileExtract

  subroutine densityProfileNames(self,names,time)
    !!{RST
    Return the names of the ``densityProfile`` properties.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityProfile), intent(inout)                             :: self
    double precision                                     , intent(in   ), optional                   :: time
    type            (varying_string                     ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names(1)="densityProfile"
    if (self%includeRadii) names(2)="densityProfileRadius"
    return
  end subroutine densityProfileNames

  subroutine densityProfileDescriptions(self,descriptions,time)
    !!{RST
    Return descriptions of the ``densityProfile`` property.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityProfile), intent(inout)                             :: self
    double precision                                     , intent(in   ), optional                   :: time
    type            (varying_string                     ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions       (1)="Density at a given radius [M☉/Mpc⁻³]."
    if (self%includeRadii)                                            &
         & descriptions(2)="Radius at which density is output [Mpc]."
    return
  end subroutine densityProfileDescriptions

  subroutine densityProfileColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnits,time)
    !!{RST
    Return column descriptions of the ``densityProfile`` property.
    !!}
    use            :: Units_MetaData          , only : unitType
    use            :: Numerical_Constants_Math, only : Pi
    use, intrinsic :: ISO_C_Binding           , only : c_int
    implicit none
    class           (nodePropertyExtractorDensityProfile), intent(inout)                            :: self
    double precision                                     , intent(in   ), optional                  :: time
    type            (varying_string                     ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                     , intent(inout), dimension(:), allocatable :: values
    type            (varying_string                     ), intent(  out)                            :: valuesDescription
    type            (unitType                           ), intent(  out)                            :: valuesUnits
    integer                                                                                         :: j
    character       (len=16                             )                                           :: label
    !$GLC attributes unused :: time

    if (self%anglesCount == 1 .and. self%polarAngles(1) == 90.0d0) then
       ! Equatorial plane only - columns are labeled by radius alone.
       allocate(descriptions(self%radiiCount))
       allocate(values      (              0))
       valuesDescription=var_str('')
       valuesUnits      =unitType(1.0d0)
       descriptions     =self%radii%specifiers%name
    else
       ! A grid in radius and polar angle - label each column by both, and give the polar angle as the column value.
       allocate(descriptions(self%radiiCount*self%anglesCount))
       allocate(values      (self%radiiCount*self%anglesCount))
       do j=1,self%anglesCount
          ! The label is built by concatenation, not a fixed-width edit descriptor: "θ" and "°" are multi-byte in UTF-8.
          write (label,'(f9.4)') self%polarAngles(j)
          descriptions((j-1)*self%radiiCount+1:j*self%radiiCount)=self%radii%specifiers%name//" θ="//trim(adjustl(label))//"°"
          values      ((j-1)*self%radiiCount+1:j*self%radiiCount)=self%polarAngles(j)
       end do
       valuesDescription=var_str('Polar angle, measured from the symmetry axis of the component, at which the density is evaluated [degrees].')
       valuesUnits      =unitType(Pi/180.0d0,"degrees","degree")
    end if
    return
  end subroutine densityProfileColumnDescriptions

  function densityProfileUnitsInSI(self,time)
    !!{RST
    Return the units of the ``densityProfile`` properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                     , allocatable  , dimension(:) :: densityProfileUnitsInSI
    class           (nodePropertyExtractorDensityProfile), intent(inout)               :: self
    double precision                                     , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(densityProfileUnitsInSI(self%elementCount_))
    densityProfileUnitsInSI       (1)=massSolar/megaParsec**3
    if (self%includeRadii)                                    &
         & densityProfileUnitsInSI(2)=          megaParsec
    return
  end function densityProfileUnitsInSI

  function densityProfileUnits(self,time) result(units)
    !!{RST
    Return the units of the densityProfile properties.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    use :: Units_MetaData                  , only : unitType
    implicit none
    type            (unitType                           ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorDensityProfile), intent(inout)             :: self
    double precision                                     , intent(in   ), optional   :: time

    allocate(units(self%elementCount_))
    units       (1)=unitType(massSolar/megaParsec**3,description='M☉/Mpc³',quantity='solMass/Mpc^3')
    if (self%includeRadii)                                                                           &
         & units(2)=unitType(          megaParsec   ,description='Mpc'    ,quantity='Mpc'          )
    return
  end function densityProfileUnits
