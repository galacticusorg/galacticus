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

  !+    Contributions to this file made by: Yu Zhao

  !!{
  Implementation of a heated soliton+NFW mass distribution class.
  !!}
  
  !![
  <massDistribution name="massDistributionSolitonNFWHeated">
    <description>
      A mass distribution class for fuzzy dark matter halos \citep{schive_understanding_2014} consisting of soliton core for small
      radii, transitioning to a heated \gls{nfw} profile at larger radii.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionSolitonNFWHeated
     !!{
     A mass distribution class for fuzzy dark matter halos \citep{schive_understanding_2014} consisting of soliton core for small
     radii, transitioning to a heated \gls{nfw} profile at larger radii.
     !!}
     private
     double precision                                 :: densitySolitonCentral             , radiusCore            , &
          &                                              radiusSoliton                     , massSolitionTransition, &
          &                                              massNFWTransition                 , densityTransition
     class           (massDistributionClass), pointer :: massDistributionSoliton_ => null()
     class           (massDistributionClass), pointer :: massDistributionHeated_  => null()
   contains
     final     ::                           solitonNFWHeatedDestructor
     procedure :: density                => solitonNFWHeatedDensity
     procedure :: densityGradientRadial  => solitonNFWHeatedDensityGradientRadial
     procedure :: massEnclosedBySphere   => solitonNFWHeatedMassEnclosedBySphere
     procedure :: radiusEnclosingMass    => solitonNFWHeatedRadiusEnclosingMass
     procedure :: radiusEnclosingDensity => solitonNFWHeatedRadiusEnclosingDensity
     procedure :: densityRadialMoment    => solitonNFWHeatedDensityRadialMoment
  end type massDistributionSolitonNFWHeated
  
  interface massDistributionSolitonNFWHeated
     !!{
     Constructors for the \refClass{massDistributionSolitonNFWHeated} mass distribution class.
     !!}
     module procedure massDistributionSolitonNFWHeatedConstructorParameters
     module procedure massDistributionSolitonNFWHeatedConstructorInternal
  end interface massDistributionSolitonNFWHeated

contains
  
  function massDistributionSolitonNFWHeatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSolitonNFWHeated} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSolitonNFWHeated)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (massDistributionClass           ), pointer       :: massDistributionHeated_
    double precision                                                  :: radiusCore             , radiusSoliton              , &
         &                                                               densitySolitonCentral  , toleranceRelativePotential
    logical                                                           :: dimensionless
    type            (varying_string                  )                :: componentType          , massType
    
    !![
    <inputParameter>
      <name>radiusCore</name>
      <description>The soliton core radius.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusSoliton</name>
      <description>The soliton radius.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densitySolitonCentral</name>
      <description>The central density of the soliton.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativePotential</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The relative tolerance for numerical solutions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the NFW profile is considered to be dimensionless.</description>
      <source>parameters</source>
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
    <objectBuilder class="massDistribution"        name="massDistributionHeated_"        source="parameters"/>
    !!]
    select type (massDistributionHeated_)
    class is (massDistributionSpherical)
       self=massDistributionSolitonNFWHeated(radiusCore,radiusSoliton,densitySolitonCentral,toleranceRelativePotential,dimensionless,massDistributionHeated_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistributionHeated_"/>
    !!]
    return
  end function massDistributionSolitonNFWHeatedConstructorParameters

  function massDistributionSolitonNFWHeatedConstructorInternal(radiusCore,radiusSoliton,densitySolitonCentral,toleranceRelativePotential,dimensionless,massDistributionHeated_,componentType,massType) result(self)
    !!{
    Internal constructor for the \refClass{massDistributionSolitonNFWHeated} mass distribution class.
    !!}
     use :: Error                     , only : Error_Report
     use :: Galactic_Structure_Options, only : componentTypeDarkHalo
     use :: Coordinates               , only : coordinateSpherical  , assignment(=)
     implicit none
     type            (massDistributionSolitonNFWHeated)                       :: self
     class           (massDistributionSpherical       ), intent(in), target   :: massDistributionHeated_
     double precision                                  , intent(in), optional :: radiusCore             , radiusSoliton             , &
          &                                                                      densitySolitonCentral  , toleranceRelativePotential
     logical                                           , intent(in), optional :: dimensionless
     type            (enumerationComponentTypeType    ), intent(in), optional :: componentType
     type            (enumerationMassTypeType         ), intent(in), optional :: massType
     type            (coordinateSpherical             )                       :: coordinates
     type            (kinematicsDistributionSoliton   ), pointer              :: kinematicsDistribution_
     !![
     <constructorAssign variables="radiusCore,radiusSoliton,densitySolitonCentral,toleranceRelativePotential,dimensionless,*massDistributionHeated_,componentType,massType"/>
     !!]

     if (present(radiusCore   )) then
        self%radiusCore             =+radiusCore
     else
        call Error_Report('core radius must be specified')
     end if
     if (present(radiusSoliton   )) then
        self%radiusSoliton          =+radiusSoliton
     else
        call Error_Report('soliton radius must be specified')
     end if
     if (present(densitySolitonCentral)) then
        self%densitySolitonCentral  =+densitySolitonCentral
     else
        call Error_Report('densitySolitonCentral must be specified')
     end if
     if (present(toleranceRelativePotential)) then
        self%toleranceRelativePotential = toleranceRelativePotential
     else
        self%toleranceRelativePotential = 1.0d-3
     end if
     if (present(dimensionless)) then
        self%dimensionless          =dimensionless
     else
        self%dimensionless          =.false.
     end if

     allocate(massDistributionSoliton :: self%massDistributionSoliton_)
     select type (massDistributionSoliton_ => self%massDistributionSoliton_)
     type is (massDistributionSoliton)
        !![
	<referenceConstruct isResult="yes" owner="self" nameAssociated="massDistributionSoliton_" object="massDistributionSoliton_">
	  <constructor>
	    massDistributionSoliton(                                               &amp;
	    &amp;                  radiusCore            = radiusCore           , &amp;
	    &amp;                  densitySolitonCentral = densitySolitonCentral, &amp;
            &amp;                  componentType         = componentTypeDarkHalo, &amp;
            &amp;                  massType              = massTypeDark           &amp;
	    &amp;                 )
	  </constructor>
	</referenceConstruct>
        !!]
        allocate(kinematicsDistribution_)
        !![
	<referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionSoliton()"/>
        !!]
        call massDistributionSoliton_%setKinematicsDistribution(kinematicsDistribution_)
        !![
	<objectDestructor name="kinematicsDistribution_"       />
        !!]
     end select
     ! Determine the masses of soliton profiles at the soliton radius.
     coordinates                =[self%radiusSoliton,0.0d0,0.0d0]
     self%massSolitionTransition=+self%massDistributionSoliton_%massEnclosedBySphere(self%radiusSoliton)
     self%massNFWTransition     =+self%massDistributionHeated_ %massEnclosedBySphere(self%radiusSoliton)
     self%densityTransition     =+self%massDistributionSoliton_%density             (coordinates       )
     return
  end function massDistributionSolitonNFWHeatedConstructorInternal

  subroutine solitonNFWHeatedDestructor(self)
    !!{
    Destructor for the \refClass{massDistributionSolitonNFWHeated} mass distribution class.
    !!}
    implicit none
    type(massDistributionSolitonNFWHeated), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistributionHeated_" />
    <objectDestructor name="self%massDistributionSoliton_"/>
    !!]
    return
  end subroutine solitonNFWHeatedDestructor

  double precision function solitonNFWHeatedDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a solitonNFWHeated mass distribution.
    !!}
    implicit none
    class(massDistributionSolitonNFWHeated), intent(inout) :: self
    class(coordinate                      ), intent(in   ) :: coordinates

    if (coordinates%rSpherical() < self%radiusSoliton) then
        solitonNFWHeatedDensity=+self%massDistributionSoliton_%density(coordinates)
    else
        solitonNFWHeatedDensity=+self%massDistributionHeated_ %density(coordinates)
    end if
    return
  end function solitonNFWHeatedDensity

  double precision function solitonNFWHeatedDensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a soliton+heated NFW mass distribution.
    !!}
    implicit none
    class  (massDistributionSolitonNFWHeated), intent(inout), target   :: self
    class  (coordinate                      ), intent(in   )           :: coordinates
    logical                                  , intent(in   ), optional :: logarithmic

    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    if (coordinates%rSpherical() < self%radiusSoliton) then
        densityGradientRadial=+self%massDistributionSoliton_%densityGradientRadial(coordinates,logarithmic)
    else
        densityGradientRadial=+self%massDistributionHeated_ %densityGradientRadial(coordinates,logarithmic)
    end if
    return
  end function solitonNFWHeatedDensityGradientRadial

  double precision function solitonNFWHeatedMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for soliton+heated NFW mass distributions.
    !!}
    implicit none
    class           (massDistributionSolitonNFWHeated), intent(inout), target :: self
    double precision                                  , intent(in   )         :: radius
    
    if (radius < self%radiusSoliton) then
        mass=+self%massDistributionSoliton_%massEnclosedBySphere(radius)
    else
        mass=+self%massDistributionHeated_ %massEnclosedBySphere(radius) &
           & -self%massNFWTransition                                     &
           & +self%massSolitionTransition
    end if
    return
  end function solitonNFWHeatedMassEnclosedBySphere
  
  double precision function solitonNFWHeatedRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for soliton+heated NFW mass distributions.
    !!}    
    implicit none
    class           (massDistributionSolitonNFWHeated), intent(inout), target   :: self
    double precision                                  , intent(in   ), optional :: mass      , massFractional
    double precision                                                            :: massHeated

    if (mass < +self%massSolitionTransition) then
        radius    =+self%massDistributionSoliton_%radiusEnclosingMass(mass      ,massFractional)
    else
        massHeated=+mass                              &
           &       -self%massSolitionTransition       &
           &       +self%massNFWTransition
        radius    =+self%massDistributionHeated_% radiusEnclosingMass(massHeated,massFractional)
    end if
    return
  end function solitonNFWHeatedRadiusEnclosingMass
  
  double precision function solitonNFWHeatedRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mean density for soliton+heated NFW mass distributions.
    !!}
    implicit none
    class           (massDistributionSolitonNFWHeated), intent(inout), target   :: self
    double precision                                  , intent(in   )           :: density
    double precision                                  , intent(in   ), optional :: radiusGuess

    if (density > +self%densityTransition) then
        radius = +self%massDistributionSoliton_%radiusEnclosingDensity(density,radiusGuess)
    else
        radius = +self%massDistributionHeated_%radiusEnclosingDensity(density,radiusGuess)
    end if
    return    
  end function solitonNFWHeatedRadiusEnclosingDensity

  double precision function solitonNFWHeatedDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Returns a radial density moment for the solitonNFWHeated mass distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionSolitonNFWHeated), intent(inout)           :: self
    double precision                                  , intent(in   )           :: moment
    double precision                                  , intent(in   ), optional :: radiusMinimum            , radiusMaximum
    logical                                           , intent(  out), optional :: isInfinite
    logical                                                                     :: isInfiniteSoliton=.false., isInfiniteHeated=.false.

    if (present(radiusMaximum)) then
       if (radiusMaximum <= self%radiusSoliton) then
          ! Only in the soliton region
          if (.not.present(radiusMinimum) .and. moment <= -1.0d0 ) then
             isInfiniteSoliton=.true.
             solitonNFWHeatedDensityRadialMoment=+0.0d0
          else
             solitonNFWHeatedDensityRadialMoment=+self%massDistributionSoliton_%densityRadialMoment(moment,     radiusMinimum,     radiusMaximum)
          end if
       else if (present(radiusMinimum) .and. radiusMinimum >= self%radiusSoliton) then
          ! Only in the heated region, no need to consider inifinite
          solitonNFWHeatedDensityRadialMoment   =+self%massDistributionHeated_ %densityRadialMoment(moment,     radiusMinimum,     radiusMaximum)
       else
          ! Consider both the soliton region and heated region
          if (.not.present(radiusMinimum) .and. moment <= -1.0d0 ) then
             isInfiniteSoliton=.true.
             solitonNFWHeatedDensityRadialMoment=+0.0d0
          else
             solitonNFWHeatedDensityRadialMoment=+self%massDistributionSoliton_%densityRadialMoment(moment,     radiusMinimum,self%radiusSoliton) &
                  &                              +self%massDistributionHeated_ %densityRadialMoment(moment,self%radiusSoliton,     radiusMaximum)
          end if
       end if
    else
       ! Consider both the soliton region and heated region
       if (                                  moment >= +2.0d0) isInfiniteHeated =.true.
       if (.not.present(radiusMinimum) .and. moment <= -1.0d0) isInfiniteSoliton=.true.
       if (.not.isInfiniteSoliton .and. .not.isInfiniteHeated) then
           solitonNFWHeatedDensityRadialMoment  =+self%massDistributionSoliton_%densityRadialMoment(moment,     radiusMinimum,self%radiusSoliton) &
              &                                  +self%massDistributionHeated_ %densityRadialMoment(moment,self%radiusSoliton,     radiusMaximum)
       else
           solitonNFWHeatedDensityRadialMoment  =+0.0d0
       end if
    end if
    ! Indicate if the moment is infinite.
    if (present(isInfinite)) then
       isInfinite=isInfiniteSoliton .or. isInfiniteHeated
    else if (isInfiniteSoliton .or. isInfiniteHeated) then
       call Error_Report('requested radial density moment is infinite'//{introspection:location})
    end if
    return
  end function solitonNFWHeatedDensityRadialMoment

