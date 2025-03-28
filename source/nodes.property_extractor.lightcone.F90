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

!+    Contributions to this file made by: Andrew Robertson, Andrew Benson
  
  use            :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
  use            :: Geometry_Lightcones, only : geometryLightcone , geometryLightconeClass
  use, intrinsic :: ISO_C_Binding      , only : c_size_t

  !![
  <nodePropertyExtractor name="nodePropertyExtractorLightcone">
   <description>
    A lightcone output extractor property extractor class. The position (and velocity and redshift) of a galaxy within a
    lightcone will be extracted. Specifically, these properties are extracted as:
    \begin{description}
     \item [{\normalfont \ttfamily lightconePositionX}] Position of the galaxy (in comoving Mpc) along the radial direction of
     the lightcone;
     \item [{\normalfont \ttfamily lightconePositionY}] Position of the galaxy (in comoving Mpc) along the 1$^\mathrm{st}$
     angular direction of the lightcone;
     \item [{\normalfont \ttfamily lightconePositionZ}] Position of the galaxy (in comoving Mpc) along the 2$^\mathrm{nd}$
     angular direction of the lightcone;
     \item [{\normalfont \ttfamily lightconeVelocityX}] Velocity of the galaxy (in km/s) along the radial direction of the
     lightcone;
     \item [{\normalfont \ttfamily lightconeVelocityY}] Velocity of the galaxy (in km/s) along the 1$^\mathrm{st}$ angular
     direction of the lightcone;
     \item [{\normalfont \ttfamily lightconeVelocityZ}] Velocity of the galaxy (in km/s) along the 2$^\mathrm{nd}$ angular
     direction of the lightcone;
     \item [{\normalfont \ttfamily lightconeRedshift}] Redshift of the galaxy in the lightcone\footnote{Note that this will
     not, in general, be precisely the same as the redshift corresponding to the output time.};
     \item [{\normalfont \ttfamily angularWeight}] The mean number density of this galaxy per unit area on the sky (in
     degrees$^{-2}$).
    \end{description}
    In order to allow this output a lightcone geometry (see \refPhysics{geometryLightcone}) must be specified.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorLightcone
     !!{
     A property extractor which extracts lightcone properties.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_      => null()
     class           (geometryLightconeClass ), pointer                   :: geometryLightcone_       => null()
     integer                                                              :: elementCount_                     , redshiftObservedOffset   , &
          &                                                                  angularCoordinatesOffset
     integer         (c_size_t               )                            :: instanceIndex
     type            (varying_string         ), allocatable, dimension(:) :: names_                            , descriptions_
     double precision                         , allocatable, dimension(:) :: unitsInSI_
     logical                                                              :: includeObservedRedshift           , includeAngularCoordinates, &
          &                                                                  atCrossing                        , failIfNotInLightcone
   contains
     final     ::                 lightconeDestructor
     procedure :: elementCount => lightconeElementCount
     procedure :: extract      => lightconeExtract
     procedure :: names        => lightconeNames
     procedure :: descriptions => lightconeDescriptions
     procedure :: unitsInSI    => lightconeUnitsInSI
     procedure :: addInstances => lightconeAddInstances
  end type nodePropertyExtractorLightcone

  interface nodePropertyExtractorLightcone
     !!{
     Constructors for the {\normalfont \ttfamily lightcone} output extractor class.
     !!}
     module procedure lightconeConstructorParameters
     module procedure lightconeConstructorInternal
  end interface nodePropertyExtractorLightcone

contains

  function lightconeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily lightcone} output extractor property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorLightcone)                :: self
    type   (inputParameters               ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    class  (geometryLightconeClass        ), pointer       :: geometryLightcone_
    logical                                                :: includeObservedRedshift, includeAngularCoordinates, &
         &                                                    atCrossing             , failIfNotInLightcone

    !![
    <inputParameter>
      <name>includeObservedRedshift</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true output the observed redshift (i.e. including the effects of peculiar velocities).</description>
    </inputParameter>
    <inputParameter>
      <name>includeAngularCoordinates</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true output angular coordinates in the lightcone.</description>
    </inputParameter>
    <inputParameter>
      <name>atCrossing</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true output positions/velocities at the time of lightcone crossing. Otherwise, output positions at the output time.</description>
    </inputParameter>
    <inputParameter>
      <name>failIfNotInLightcone</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, a node that is not in the lightcone will cause a fatal error. Otherwise, such nodes are simply assigned unphysical values for lightcone properties.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="geometryLightcone"  name="geometryLightcone_"  source="parameters"/>
    !!]
    self=nodePropertyExtractorLightcone(includeObservedRedshift,includeAngularCoordinates,atCrossing,failIfNotInLightcone,cosmologyFunctions_,geometryLightcone_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="geometryLightcone_" />
    !!]
    return
  end function lightconeConstructorParameters

  function lightconeConstructorInternal(includeObservedRedshift,includeAngularCoordinates,atCrossing,failIfNotInLightcone,cosmologyFunctions_,geometryLightcone_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily lightcone} output extractor property extractor class.
    !!}
    use :: Numerical_Constants_Astronomical, only : degreesToRadians, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    type   (nodePropertyExtractorLightcone)                        :: self
    logical                                , intent(in   )         :: includeObservedRedshift, includeAngularCoordinates, &
         &                                                            atCrossing             , failIfNotInLightcone
    class  (cosmologyFunctionsClass       ), intent(in   ), target :: cosmologyFunctions_
    class  (geometryLightconeClass        ), intent(in   ), target :: geometryLightcone_
    !![
    <constructorAssign variables="includeObservedRedshift, includeAngularCoordinates, atCrossing, failIfNotInLightcone, *cosmologyFunctions_, *geometryLightcone_"/>
    !!]

    self%elementCount_=8
    if (includeObservedRedshift) then
       self%redshiftObservedOffset  =self%elementCount_
       self%elementCount_           =self%elementCount_+1
    end if
    if (includeAngularCoordinates) then
       self%angularCoordinatesOffset=self%elementCount_
       self%elementCount_           =self%elementCount_+2
    end if
    allocate(self%names_       (self%elementCount_))
    allocate(self%descriptions_(self%elementCount_))
    allocate(self%unitsInSI_   (self%elementCount_))
    self%names_       (1)='lightconePositionX'
    self%descriptions_(1)='Position of galaxy in lightcone (radial axis) [Mpc].'
    self%unitsInSI_   (1)=megaParsec
    self%names_       (2)='lightconePositionY'
    self%descriptions_(2)='Position of galaxy in lightcone (1st angular axis) [Mpc].'
    self%unitsInSI_   (2)=megaParsec
    self%names_       (3)='lightconePositionZ'
    self%descriptions_(3)='Position of galaxy in lightcone (2nd angular axis) [Mpc].'
    self%unitsInSI_   (3)=megaParsec
    self%names_       (4)='lightconeVelocityX'
    self%descriptions_(4)='Velocity of galaxy in lightcone (radial axis) [km/s].'
    self%unitsInSI_   (4)=kilo
    self%names_       (5)='lightconeVelocityY'
    self%descriptions_(5)='Velocity of galaxy in lightcone (1st angular axis) [km/s].'
    self%unitsInSI_   (5)=kilo
    self%names_       (6)='lightconeVelocityZ'
    self%descriptions_(6)='Velocity of galaxy in lightcone (2nd angular axis) [km/s].'
    self%unitsInSI_   (6)=kilo
    self%names_       (7)='lightconeRedshiftCosmological'
    self%descriptions_(7)='Cosmological redshift of galaxy in lightcone (does not include the effects of line-of-sight peculiar velocity).'
    self%unitsInSI_   (7)=0.0d0
    self%names_       (8)='angularWeight'
    self%descriptions_(8)='Number of such galaxies per unit area [degrees⁻²].'
    self%unitsInSI_   (8)=degreesToRadians**2
    if (includeObservedRedshift) then
       self%names_       (self%redshiftObservedOffset  +1)='lightconeRedshiftObserved'
       self%descriptions_(self%redshiftObservedOffset  +1)='Observed redshift of galaxy in lightcone (includes the effects of line-of-sight peculiar velocity).'
       self%unitsInSI_   (self%redshiftObservedOffset  +1)=0.0d0
    end if
    if (includeAngularCoordinates) then
       self%names_       (self%angularCoordinatesOffset+1)='lightconeAngularTheta'
       self%descriptions_(self%angularCoordinatesOffset+1)='Angular distance from pole of coordinate system (i.e. θ in a spherical coordinate system) [radians]'
       self%unitsInSI_   (self%angularCoordinatesOffset+1)=1.0d0
       self%names_       (self%angularCoordinatesOffset+2)='lightconeAngularPhi'
       self%descriptions_(self%angularCoordinatesOffset+2)='Angular distance around the pole of coordinate system (i.e. φ in a spherical coordinate system) [radians]'
       self%unitsInSI_   (self%angularCoordinatesOffset+2)=1.0d0
    end if
    return
  end function lightconeConstructorInternal

  subroutine lightconeDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily lightcone} output extractor property extractor class.
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
    !!{
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
    !!{
    Implement a lightcone output extractor.
    !!}
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
    double precision                                                             :: velocityBeta
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
    if (self%includeObservedRedshift) then
       ! Compute the line-of-sight peculiar velocity divided by the speed of light: β_los=v_los/c.
       velocityBeta                                     =+Dot_Product     (lightconeExtract(4:6),lightconeExtract(1:3)) &
            &                                            /Vector_Magnitude(                      lightconeExtract(1:3)) &
            &                                            *kilo                                                          &
            &                                            /speedLight
       ! Compute the observed redshift. This is given by:
       !  1 + zₒ = (1 + zₕ) (1 + zₚ),
       ! where zₒ is observed redshift, zₕ is cosmological redshift (due
       ! to Hubble expansion), and zₚ is the peculiar redshift, given by
       ! (in the non-relativistic limit) zₚ = β_los (e.g. Davis et al.;
       ! 2011; ApJ; 741; 67; below eqn. 4;
       ! https://ui.adsabs.harvard.edu/abs/2011ApJ...741...67D).
       lightconeExtract(self%redshiftObservedOffset+1)=  -1.0d0                      &
            &                                          +(+1.0d0+lightconeExtract(7)) &
            &                                          *(+1.0d0+velocityBeta       )
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
    !!{
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
    !!{
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
    !!{
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
    !!{
    Return the units of the lightcone properties in the SI system.
    !!}
    implicit none
    double precision                                , dimension(:) , allocatable :: lightconeUnitsInSI
    class           (nodePropertyExtractorLightcone), intent(inout)              :: self
    double precision                                , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(lightconeUnitsInSI(self%elementCount_))
    lightconeUnitsInSI=self%unitsInSI_
    return
  end function lightconeUnitsInSI

