!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !% Contains a class which implements the tidal heating rate model of \cite{gnedin_tidal_1999}.

  use Cosmology_Parameters
  use Dark_Matter_Halo_Scales

  !# <satelliteTidalHeatingRate name="satelliteTidalHeatingRateGnedin1999">
  !#  <description>A satellite tidal heating rate class which implements the tidal heating rate model of \cite{gnedin_tidal_1999}.</description>
  !# </satelliteTidalHeatingRate>
  type, extends(satelliteTidalHeatingRateClass) :: satelliteTidalHeatingRateGnedin1999
     !% A satellite tidal heating rate class which implements the tidal heating rate model of \cite{gnedin_tidal_1999}.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: epsilon             , gamma
   contains
     final     ::                gnedin1999Destructor
     procedure :: heatingRate => gnedin1999HeatingRate
  end type satelliteTidalHeatingRateGnedin1999

  interface satelliteTidalHeatingRateGnedin1999
     !% Constructors for the gnedin1999 satellite tidal heating rate class.
     module procedure gnedin1999ConstructorParameters
     module procedure gnedin1999ConstructorInternal
  end interface satelliteTidalHeatingRateGnedin1999

contains
  
  function gnedin1999ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily gnedin1999} satellite tidal heating rate class which builds the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (satelliteTidalHeatingRateGnedin1999)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyParametersClass           ), pointer       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    double precision                                                     :: epsilon             , gamma

    !# <inputParameter>
    !#   <name>epsilon</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>3.0d0</defaultValue>
    !#   <description>Parameter, $\epsilon$, controlling the tidal heating rate of satellites in the {\normalfont \ttfamily Gnedin1999} method.</description>
    !#   <group></group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>gamma</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>2.5d0</defaultValue>
    !#   <description>Parameter, $\gamma$, controlling the tidal heating rate of satellites in the {\normalfont \ttfamily Gnedin1999} method.</description>
    !#   <group></group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=satelliteTidalHeatingRateGnedin1999(epsilon,gamma,cosmologyParameters_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function gnedin1999ConstructorParameters

  function gnedin1999ConstructorInternal(epsilon,gamma,cosmologyParameters_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily gnedin1999} satellite tidal heating rate class.
    implicit none
    type            (satelliteTidalHeatingRateGnedin1999)                        :: self
    class           (cosmologyParametersClass           ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                     , intent(in)            :: epsilon, gamma
    !# <constructorAssign variables="epsilon, gamma, *cosmologyParameters_, *darkMatterHaloScale_"/>

    return
  end function gnedin1999ConstructorInternal

  subroutine gnedin1999Destructor(self)
    !% Default constructor for the {\normalfont \ttfamily gnedin1999} satellite tidal heating rate class.
    implicit none
    type(satelliteTidalHeatingRateGnedin1999), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine gnedin1999Destructor

  double precision function gnedin1999HeatingRate(self,node)
    !% Return the tidal heating rate for satellite halos assuming the model of \cite{gnedin_tidal_1999}.
    use Galacticus_Nodes                  , only : nodeComponentBasic, nodeComponentSatellite
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    use Error_Functions
    use Galactic_Structure_Densities
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Rotation_Curves
    use Galactic_Structure_Options
    use Vectors
    use Tensors
    implicit none
    class           (satelliteTidalHeatingRateGnedin1999), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    class           (nodeComponentSatellite             ), pointer       :: satellite
    type            (treeNode                           ), pointer       :: nodeHost
    class           (nodeComponentBasic                 ), pointer       :: basic
    double precision                                     , dimension(3)  :: position                 , velocity
    double precision                                                     :: satelliteMass            , densityParent            , &
         &                                                                  massEnclosedHost         , velocityCircularSatellite, &
         &                                                                  radius                   , speed                    , &
         &                                                                  timescaleShock           , heatingRateNormalized    , &
         &                                                                  orbitalFrequencySatellite, radiusHalfMassSatellite  , &
         &                                                                  satelliteHalfMass        , fractionDarkMatter
    type            (tensorRank2Dimension3Symmetric     )                :: tidalTensor              , tidalTensorPathIntegrated, &
         &                                                                  positionTensor

    ! Construct required properties of satellite and host.
    nodeHost                  => node     %mergesWith               (                                                             )
    satellite                 => node     %satellite                (                                                             )
    satelliteMass             =  satellite%boundMass                (                                                             )
    position                  =  satellite%position                 (                                                             )
    velocity                  =  satellite%velocity                 (                                                             )
    tidalTensorPathIntegrated =  satellite%tidalTensorPathIntegrated(                                                             )
    radius                    =  Vector_Magnitude                   (         position                                            )
    speed                     =  Vector_Magnitude                   (         velocity                                            )
    densityParent             =  Galactic_Structure_Density         (nodeHost,position,coordinateSystemCartesian                  )
    massEnclosedHost          =  Galactic_Structure_Enclosed_Mass   (nodeHost,radius                                              )
    positionTensor            =  Vector_Outer_Product               (         position                          ,symmetrize=.true.)
    ! Find the universal dark matter fraction.
    fractionDarkMatter        =  +(                                        &
         &                         +self%cosmologyParameters_%OmegaMatter() &
         &                         -self%cosmologyParameters_%OmegaBaryon() &
         &                        )                                         &
         &                       /  self%cosmologyParameters_%OmegaMatter()
    ! Find the gravitational tidal tensor.
    tidalTensor=                                                                                        &
         & -(gravitationalConstantGalacticus*massEnclosedHost         /radius**3)*tensorIdentityR2D3Sym &
         & +(gravitationalConstantGalacticus*massEnclosedHost*3.0d0   /radius**5)*positionTensor        &
         & -(gravitationalConstantGalacticus*densityParent   *4.0d0*Pi/radius**2)*positionTensor
    ! Find the orbital frequency at the half mass radius of the satellite.
    basic                    => node%basic()
    satelliteHalfMass        =  +0.50d0             &
         &                      *fractionDarkMatter &
         &                      *min(               &
         &                           satelliteMass, &
         &                           basic%mass()   &
         &                      )    
    radiusHalfMassSatellite  =  Galactic_Structure_Radius_Enclosing_Mass(                                         &
         &                                                             node                                     , &
         &                                                             mass                   =satelliteHalfMass, &
         &                                                             componentType          =componentTypeAll , &
         &                                                             massType               =massTypeDark       &
         &                                                            )
    velocityCircularSatellite=Galactic_Structure_Rotation_Curve       (                                           &
         &                                                             node                                     , &
         &                                                             radiusHalfMassSatellite                  , &
         &                                                             componentType          =componentTypeAll , &
         &                                                             massType               =massTypeDark       &
         &                                                            )
    if (radiusHalfMassSatellite > 0.0d0) then
       ! Compute the orbital frequency.
       orbitalFrequencySatellite =  +velocityCircularSatellite &
            &                       /radiusHalfMassSatellite   &
            &                       *gigaYear                  &
            &                       *kilo                      &
            &                       /megaParsec
    else
       ! No well-defined half mass radius exists in the satellite. Use the virial orbital frequency instead.
       orbitalFrequencySatellite =  +self%darkMatterHaloScale_%virialVelocity(node) &
            &                       /self%darkMatterHaloScale_%virialRadius  (node) &
            &                       *gigaYear                                       &
            &                       *kilo                                           &
            &                       /megaParsec
    end if
    ! Find the shock timescale (i.e. crossing time in the radial direction).
    timescaleShock=+megaParsec &
         &         /kilo       &
         &         /gigaYear   &
         &         *radius     &
         &         /speed
    ! Compute the heating rate.
    heatingRateNormalized=+self%epsilon                                          &
         &                /(                                                     &
         &                  +1.0d0                                               &
         &                  +(                                                   &
         &                    +timescaleShock                                    &
         &                    *orbitalFrequencySatellite                         &
         &                   )**2                                                &
         &                 )**self%gamma                                         &
         &                /3.0d0                                                 &
         &                *tidalTensor%doubleContract(tidalTensorPathIntegrated) &
         &                *(kilo*gigaYear/megaParsec)**2
    ! Limit the heating rate to be non-negative.
    gnedin1999HeatingRate=max(heatingRateNormalized,0.0d0)
    return
  end function gnedin1999HeatingRate


