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

!+    Contributions to this file made by: Andrew Robertson, Andrew Benson
  
  use            :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
  use            :: Geometry_Lightcones, only : geometryLightcone , geometryLightconeClass
  use, intrinsic :: ISO_C_Binding      , only : c_size_t

  !![
  <nodePropertyExtractor name="nodePropertyExtractorLightcone" docformat="rst">
   <description>
   A lightcone output extractor property extractor class. The position (and velocity and redshift) of a galaxy within a lightcone will be extracted. Specifically, these properties are extracted as:

   ``lightconePositionX``
      Position of the galaxy (in comoving Mpc) along the radial direction of the lightcone;

   ``lightconePositionY``
      Position of the galaxy (in comoving Mpc) along the 1\ :math:`^\mathrm{st}` angular direction of the lightcone;

   ``lightconePositionZ``
      Position of the galaxy (in comoving Mpc) along the 2\ :math:`^\mathrm{nd}` angular direction of the lightcone;

   ``lightconePositionObservedX``
      Position of the galaxy (in comoving Mpc) along the radial direction of the lightcone, accounting for the effects of line-of-sight peculiar velocity (included only if ``[includeObservedPosition]``\ =\ ``true``);

   ``lightconePositionObservedY``
      Position of the galaxy (in comoving Mpc) along the 1\ :math:`^\mathrm{st}` angular direction of the lightcone, accounting for the effects of line-of-sight peculiar velocity (included only if ``[includeObservedPosition]``\ =\ ``true``);

   ``lightconePositionObservedZ``
      Position of the galaxy (in comoving Mpc) along the 2\ :math:`^\mathrm{nd}` angular direction of the lightcone, accounting for the effects of line-of-sight peculiar velocity (included only if ``[includeObservedPosition]``\ =\ ``true``);

   ``lightconeVelocityX``
      Velocity of the galaxy (in km/s) along the radial direction of the lightcone;

   ``lightconeVelocityY``
      Velocity of the galaxy (in km/s) along the 1\ :math:`^\mathrm{st}` angular direction of the lightcone;

   ``lightconeVelocityZ``
      Velocity of the galaxy (in km/s) along the 2\ :math:`^\mathrm{nd}` angular direction of the lightcone;

   ``lightconeRedshiftCosmological``
      Redshift of the galaxy in the lightcone\footnoteNote that this will not, in general, be precisely the same as the redshift corresponding to the output time.;

   ``lightconeRedshiftObserved``
      Observed redshift of the galaxy, accounting for the effects of line-of-sight peculiar velocity (included only if ``[includeObservedRedshift]``\ =\ ``true``);

   ``lightconeAngularTheta``
      Angular distance from pole of coordinate system (i.e. :math:`\theta` in a spherical coordinate system; included only if ``[includeAngularCoordinates]``\ =\ ``true``) [radians]

   ``lightconeAngularPhi``
      Angular distance around the pole of coordinate system system (i.e. :math:`\phi` in a spherical coordinate system; included only if ``[includeAngularCoordinates]``\ =\ ``true``) [radians]

   ``angularWeight``
      The mean number density of this galaxy per unit area on the sky (in degrees\ :math:`^{-2}`).

   In order to allow this output a lightcone geometry (see ``geometryLightcone``) must be specified.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorLightcone
     !!{RST
     A property extractor which extracts lightcone properties.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_      => null()
     class           (geometryLightconeClass ), pointer                   :: geometryLightcone_       => null()
     integer                                                              :: elementCount_                     , redshiftObservedOffset   , &
          &                                                                  angularCoordinatesOffset          , positionObservedOffset
     integer         (c_size_t               )                            :: instanceIndex
     type            (varying_string         ), allocatable, dimension(:) :: names_                            , descriptions_
     type            (unitType               ), allocatable, dimension(:) :: units_
     logical                                                              :: includeObservedRedshift           , includeAngularCoordinates, &
          &                                                                  atCrossing                        , failIfNotInLightcone     , &
          &                                                                  includeObservedPosition
   contains
     final     ::                 lightconeDestructor
     procedure :: elementCount => lightconeElementCount
     procedure :: extract      => lightconeExtract
     procedure :: names        => lightconeNames
     procedure :: descriptions => lightconeDescriptions
     procedure :: unitsInSI    => lightconeUnitsInSI
     procedure :: units        => lightconeUnits
     procedure :: addInstances => lightconeAddInstances
  end type nodePropertyExtractorLightcone

  interface nodePropertyExtractorLightcone
     !!{RST
     Constructors for the ``nodePropertyExtractorLightcone`` property extractor class.
     !!}
     module procedure lightconeConstructorParameters
     module procedure lightconeConstructorInternal
  end interface nodePropertyExtractorLightcone

contains

  function lightconeConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``nodePropertyExtractorLightcone`` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorLightcone)                :: self
    type   (inputParameters               ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    class  (geometryLightconeClass        ), pointer       :: geometryLightcone_
    logical                                                :: includeObservedRedshift, includeAngularCoordinates, &
         &                                                    includeObservedPosition, atCrossing               , &
         &                                                    failIfNotInLightcone

    !![
    <inputParameter docformat="rst">
      <name>includeObservedRedshift</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>
      If true output the observed redshift (i.e. including the effects of peculiar velocities).
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>includeAngularCoordinates</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>
      If true output angular coordinates in the lightcone.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>includeObservedPosition</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>
      If true output the observed position (i.e. including the effects of peculiar velocities).
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>atCrossing</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>
      If true output positions/velocities at the time of lightcone crossing. Otherwise, output positions at the output time.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>failIfNotInLightcone</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a node that is not in the lightcone will cause a fatal error. Otherwise, such nodes are simply assigned unphysical values for lightcone properties.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="geometryLightcone"  name="geometryLightcone_"  source="parameters"/>
    !!]
    self=nodePropertyExtractorLightcone(includeObservedRedshift,includeObservedPosition,includeAngularCoordinates,atCrossing,failIfNotInLightcone,cosmologyFunctions_,geometryLightcone_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="geometryLightcone_" />
    !!]
    return
  end function lightconeConstructorParameters

  function lightconeConstructorInternal(includeObservedRedshift,includeObservedPosition,includeAngularCoordinates,atCrossing,failIfNotInLightcone,cosmologyFunctions_,geometryLightcone_) result(self)
    !!{RST
    Internal constructor for the ``nodePropertyExtractorLightcone`` property extractor class.
    !!}
    use :: Numerical_Constants_Astronomical, only : degreesToRadians, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Units_MetaData                  , only : unitType
    implicit none
    type   (nodePropertyExtractorLightcone)                        :: self
    logical                                , intent(in   )         :: includeObservedRedshift, includeAngularCoordinates, &
         &                                                            includeObservedPosition, atCrossing               , &
         &                                                            failIfNotInLightcone
    class  (cosmologyFunctionsClass       ), intent(in   ), target :: cosmologyFunctions_
    class  (geometryLightconeClass        ), intent(in   ), target :: geometryLightcone_
    !![
    <constructorAssign variables="includeObservedRedshift, includeObservedPosition, includeAngularCoordinates, atCrossing, failIfNotInLightcone, *cosmologyFunctions_, *geometryLightcone_"/>
    !!]

    self%elementCount_=8
    if (includeObservedRedshift) then
       self%redshiftObservedOffset  =self%elementCount_
       self%elementCount_           =self%elementCount_+1
    end if
    if (includeObservedPosition) then
       self%positionObservedOffset  =self%elementCount_
       self%elementCount_           =self%elementCount_+3
    end if
    if (includeAngularCoordinates) then
       self%angularCoordinatesOffset=self%elementCount_
       self%elementCount_           =self%elementCount_+2
    end if
    allocate(self%names_       (self%elementCount_))
    allocate(self%descriptions_(self%elementCount_))
    allocate(self%units_       (self%elementCount_))
    self%names_       (1)='lightconePositionX'
    self%descriptions_(1)='Position of galaxy in lightcone (radial axis) [Mpc].'
    self%units_       (1)=unitType(megaParsec,'Mpc','Mpc',isComoving=.true.)
    self%names_       (2)='lightconePositionY'
    self%descriptions_(2)='Position of galaxy in lightcone (1st angular axis) [Mpc].'
    self%units_       (2)=unitType(megaParsec,'Mpc','Mpc',isComoving=.true.)
    self%names_       (3)='lightconePositionZ'
    self%descriptions_(3)='Position of galaxy in lightcone (2nd angular axis) [Mpc].'
    self%units_       (3)=unitType(megaParsec,'Mpc','Mpc',isComoving=.true.)
    self%names_       (4)='lightconeVelocityX'
    self%descriptions_(4)='Velocity of galaxy in lightcone (radial axis) [km/s].'
    self%units_       (4)=unitType(kilo,'km/s','km/s')
    self%names_       (5)='lightconeVelocityY'
    self%descriptions_(5)='Velocity of galaxy in lightcone (1st angular axis) [km/s].'
    self%units_       (5)=unitType(kilo,'km/s','km/s')
    self%names_       (6)='lightconeVelocityZ'
    self%descriptions_(6)='Velocity of galaxy in lightcone (2nd angular axis) [km/s].'
    self%units_       (6)=unitType(kilo,'km/s','km/s')
    self%names_       (7)='lightconeRedshiftCosmological'
    self%descriptions_(7)='Cosmological redshift of galaxy in lightcone (does not include the effects of line-of-sight peculiar velocity).'
    self%units_       (7)=unitType(1.0d0)
    self%names_       (8)='angularWeight'
    self%descriptions_(8)='Number of such galaxies per unit area [degrees⁻²].'
    self%units_       (8)=unitType(degreesToRadians**2,'degrees⁻²','deg^-2')
    if (includeObservedRedshift) then
       self%names_       (self%redshiftObservedOffset  +1)='lightconeRedshiftObserved'
       self%descriptions_(self%redshiftObservedOffset  +1)='Observed redshift of galaxy in lightcone (includes the effects of line-of-sight peculiar velocity).'
       self%units_       (self%redshiftObservedOffset  +1)=unitType(1.0d0)
    end if
    if (includeObservedPosition) then
       self%names_       (self%positionObservedOffset  +1)='lightconePositionObservedX'
       self%names_       (self%positionObservedOffset  +2)='lightconePositionObservedY'
       self%names_       (self%positionObservedOffset  +3)='lightconePositionObservedZ'
       self%descriptions_(self%positionObservedOffset  +1)='Observed position of galaxy in lightcone (includes the effects of line-of-sight peculiar velocity; radial axis) [Mpc].'
       self%descriptions_(self%positionObservedOffset  +2)='Observed position of galaxy in lightcone (includes the effects of line-of-sight peculiar velocity; 1st angular axis) [Mpc].'
       self%descriptions_(self%positionObservedOffset  +3)='Observed position of galaxy in lightcone (includes the effects of line-of-sight peculiar velocity; 2nd angular axis) [Mpc].'
       self%units_       (self%positionObservedOffset  +1)=unitType(megaParsec,'Mpc','Mpc',isComoving=.true.)
       self%units_       (self%positionObservedOffset  +2)=unitType(megaParsec,'Mpc','Mpc',isComoving=.true.)
       self%units_       (self%positionObservedOffset  +3)=unitType(megaParsec,'Mpc','Mpc',isComoving=.true.)
    end if
    if (includeAngularCoordinates) then
       self%names_       (self%angularCoordinatesOffset+1)='lightconeAngularTheta'
       self%descriptions_(self%angularCoordinatesOffset+1)='Angular distance from pole of coordinate system (i.e. θ in a spherical coordinate system) [radians]'
       self%units_       (self%angularCoordinatesOffset+1)=unitType(1.0d0,'radians','rad')
       self%names_       (self%angularCoordinatesOffset+2)='lightconeAngularPhi'
       self%descriptions_(self%angularCoordinatesOffset+2)='Angular distance around the pole of coordinate system (i.e. φ in a spherical coordinate system) [radians]'
       self%units_       (self%angularCoordinatesOffset+2)=unitType(1.0d0,'radians','rad')
    end if
    return
  end function lightconeConstructorInternal

  subroutine lightconeDestructor(self)
    !!{RST
    Destructor for the ``nodePropertyExtractorLightcone`` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLightcone), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%geometryLightcone_" />
    !!]
    return
  end subroutine lightconeDestructor

  integer function lightconeElementCount(self,time)
    !!{RST
    Return the number of elements in the lightcone property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorLightcone), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    !$GLC attributes unused :: time

    lightconeElementCount=self%elementCount_
    return
  end function lightconeElementCount

  function lightconeExtract(self,node,time,instance)
    !!{RST
    Implement a lightcone output extractor.
    !!}
    use :: Cosmology_Functions_Options     , only : distanceTypeComoving
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : degreesToRadians
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    double precision                                , dimension(:) , allocatable :: lightconeExtract
    class           (nodePropertyExtractorLightcone), intent(inout), target      :: self
    type            (treeNode                      ), intent(inout), target      :: node
    double precision                                , intent(in   )              :: time
    type            (multiCounter                  ), intent(inout), optional    :: instance
    double precision                                                             :: velocityBeta    , redshiftObserved      , &
         &                                                                          distanceRadial  , distanceRadialObserved
    !$GLC attributes unused :: time

    if (.not.present(instance).and..not.self%atCrossing) call Error_Report('instance is required'//{introspection:location})
    allocate(lightconeExtract(self%elementCount_))
    if     (                                                                        &
         &   .not.self                   %atCrossing                                &
         &  .and.                                                                   &
         &   .not.self%geometryLightcone_%isInLightcone(node,atPresentEpoch=.true.) &
         & ) then
       ! Node is not in the lightcone. If this is allowed, simply assign unphysical values.
       if (self%failIfNotInLightcone) call Error_Report('node is not in lightcone'//{introspection:location})
       lightconeExtract=-huge(0.0d0)
       return
    end if
    ! Node is in the lightcone - compute properties on the lightcone.
    if (self%atCrossing) then
       lightconeExtract(1:3)=self%geometryLightcone_%positionLightconeCrossing(node                                   )
       lightconeExtract(4:6)=self%geometryLightcone_%velocityLightconeCrossing(node                                   )
    else
       lightconeExtract(1:3)=self%geometryLightcone_%position                 (node,instance%state(self%instanceIndex))
       lightconeExtract(4:6)=self%geometryLightcone_%velocity                 (node,instance%state(self%instanceIndex))
    end if
    lightconeExtract   (7  )=self%cosmologyFunctions_   %redshiftFromExpansionFactor(                          &
         &                    self%cosmologyFunctions_  %expansionFactor             (                         &
         &                     self%cosmologyFunctions_ %timeAtDistanceComoving       (                        &
         &                      Vector_Magnitude                                       (lightconeExtract(1:3)) &
         &                                                                            )                        &
         &                                                                           )                         &
         &                                                                          )
    lightconeExtract   (8  )=self%geometryLightcone_%solidAngle()/degreesToRadians**2
    ! Compute observed redshift if needed.
    if (self%includeObservedRedshift .or. self%includeObservedPosition) then
       ! Compute the line-of-sight peculiar velocity divided by the speed of light: β_los=v_los/c.
       velocityBeta    =+Dot_Product     (lightconeExtract(4:6),lightconeExtract(1:3)) &
            &           /Vector_Magnitude(                      lightconeExtract(1:3)) &
            &           *kilo                                                          &
            &           /speedLight
       ! Compute the observed redshift. This is given by:
       !  1 + zₒ = (1 + zₕ) (1 + zₚ),
       ! where zₒ is observed redshift, zₕ is cosmological redshift (due
       ! to Hubble expansion), and zₚ is the peculiar redshift, given by
       ! (in the non-relativistic limit) zₚ = β_los (e.g. Davis et al.;
       ! 2011; ApJ; 741; 67; below eqn. 4;
       ! https://ui.adsabs.harvard.edu/abs/2011ApJ...741...67D).
       redshiftObserved=-1.0d0                        &
            &           +(+1.0d0+lightconeExtract(7)) &
            &           *(+1.0d0+velocityBeta       )
    else
       redshiftObserved=-huge(0.0d0)
    end if
    if (self%includeObservedRedshift) lightconeExtract(self%redshiftObservedOffset+1)=redshiftObserved
    if (self%includeObservedPosition) then
       distanceRadial                                                               =Vector_Magnitude(lightconeExtract(1:3)) 
       distanceRadialObserved                                                       =self%cosmologyFunctions_%distanceComovingConvert(output=distanceTypeComoving,redshift=redshiftObserved)
       lightconeExtract(self%positionObservedOffset+1:self%positionObservedOffset+3)=+lightconeExtract      (1:3) &
            &                                                                        *distanceRadialObserved      &
            &                                                                        /distanceRadial
    end if
    if (self%includeAngularCoordinates) then
       lightconeExtract(self%angularCoordinatesOffset+1)=atan2(                              &
            &                                                  sqrt(                         &
            &                                                       +lightconeExtract(2)**2  &
            &                                                       +lightconeExtract(3)**2  &
            &                                                      )                       , &
            &                                                        lightconeExtract(1)     &
            &                                                  )
       lightconeExtract(self%angularCoordinatesOffset+2)=atan2(                              &
            &                                                        lightconeExtract(3)   , &
            &                                                        lightconeExtract(2)     &
            &                                                 )
    end if
    return
  end function lightconeExtract

  subroutine lightconeAddInstances(self,node,instance)
    !!{RST
    Implement adding of instances to a lightcone property extractor.
    !!}
    use :: Error          , only : Error_Report
    use :: String_Handling, only : operator(//)
    implicit none
    class    (nodePropertyExtractorLightcone), intent(inout) :: self
    type     (treeNode                      ), intent(inout) :: node
    type     (multiCounter                  ), intent(inout) :: instance
    integer  (c_size_t                      )                :: replicationCount
    type     (varying_string                )                :: message
    character(len=5                         )                :: label

    if (self%atCrossing) then
       replicationCount=1_c_size_t
    else
       replicationCount=self%geometryLightcone_%replicationCount(node)
       if (replicationCount < 1_c_size_t) then
          if (self%failIfNotInLightcone) then
             if (self%geometryLightcone_%isInLightcone(node,atPresentEpoch=.true.)) then
                label="true"
             else
                label="false"
             end if
             message=var_str("Node ")//node%index()//" of tree "//node%hostTree%index//" appears in "//replicationCount//"(<1) replicants - this should not happen - lightcone intersection reports '"//trim(label)//"'"
             call Error_Report(message//{introspection:location})
          else
             ! Nodes not in the lightcone are to be accepted.
             replicationCount=1
          end if
       end if
       call instance%append(replicationCount)
       self%instanceIndex=instance%count()
    end if
    return
  end subroutine lightconeAddInstances

  subroutine lightconeNames(self,time,names)
    !!{RST
    Return the names of the lightcone properties.
    !!}
    implicit none
    class           (nodePropertyExtractorLightcone), intent(inout)                             :: self
    double precision                                , intent(in   )                             :: time
    type            (varying_string                ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names=self%names_
    return
  end subroutine lightconeNames

  subroutine lightconeDescriptions(self,time,descriptions)
    !!{RST
    Return the descriptions of the lightcone properties.
    !!}
    implicit none
    class           (nodePropertyExtractorLightcone), intent(inout)                             :: self
    double precision                                , intent(in   )                             :: time
    type            (varying_string                ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions=self%descriptions_
    return
  end subroutine lightconeDescriptions

  function lightconeUnitsInSI(self,time)
    !!{RST
    Return the units of the lightcone properties in the SI system.
    !!}
    implicit none
    double precision                                , dimension(:) , allocatable :: lightconeUnitsInSI
    class           (nodePropertyExtractorLightcone), intent(inout)              :: self
    double precision                                , intent(in   )              :: time
    integer                                                                      :: i
    !$GLC attributes unused :: time

    allocate(lightconeUnitsInSI(self%elementCount_))
    do i=1,self%elementCount_
       lightconeUnitsInSI(i)=self%units_(i)%unitsInSI
    end do
    return
  end function lightconeUnitsInSI

  function lightconeUnits(self,time) result(units)
    !!{RST
    Return the units of the lightcone properties.
    !!}
    implicit none
    type            (unitType                      ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorLightcone), intent(inout)             :: self
    double precision                                , intent(in   )             :: time
    !$GLC attributes unused :: time

    allocate(units(self%elementCount_))
    units=self%units_
    return
  end function lightconeUnits
