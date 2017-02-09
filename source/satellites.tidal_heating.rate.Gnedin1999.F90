!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module with a \cite{gnedin_tidal_1999} implementation of calculations of satellite tidal heating.

module Tidal_Heating_Rate_Gnedin
  !% Implements \cite{gnedin_tidal_1999} value of calculations of satellite tidal heating.
  implicit none
  private
  public :: Satellite_Tidal_Heating_Rate_Gnedin_Initialize

  ! Parameters controlling the heating rate.
  double precision :: satelliteTidalHeatingGnedinGamma, satelliteTidalHeatingGnedinEpsilon

contains

  !# <satelliteTidalHeatingMethod>
  !#  <unitName>Satellite_Tidal_Heating_Rate_Gnedin_Initialize</unitName>
  !# </satelliteTidalHeatingMethod>
  subroutine Satellite_Tidal_Heating_Rate_Gnedin_Initialize(satelliteTidalHeatingMethod,Satellite_Tidal_Heating_Rate)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                     ),          intent(in   ) :: satelliteTidalHeatingMethod
    procedure(Satellite_Tidal_Heating_Rate_Gnedin), pointer, intent(inout) :: Satellite_Tidal_Heating_Rate

    if (satelliteTidalHeatingMethod == 'Gnedin1999') then
       Satellite_Tidal_Heating_Rate => Satellite_Tidal_Heating_Rate_Gnedin
       !@ <inputParameter>
       !@   <name>satelliteTidalHeatingGnedinEpsilon</name>
       !@   <defaultValue>3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Parameter, $\epsilon$, controlling the tidal heating rate of satellites in the {\normalfont \ttfamily Gnedin1999} method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group></group>
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteTidalHeatingGnedinEpsilon',satelliteTidalHeatingGnedinEpsilon,defaultValue=3.0d0)
       !@ <inputParameter>
       !@   <name>satelliteTidalHeatingGnedinGamma</name>
       !@   <defaultValue>2.5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Parameter, $\gamma$, controlling the tidal heating rate of satellites in the {\normalfont \ttfamily Gnedin1999} method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group></group>
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteTidalHeatingGnedinGamma',satelliteTidalHeatingGnedinGamma,defaultValue=2.5d0)
     end if
    return
  end subroutine Satellite_Tidal_Heating_Rate_Gnedin_Initialize

  double precision function Satellite_Tidal_Heating_Rate_Gnedin(thisNode)
    !% Return the \cite{gnedin_tidal_1999} rate for satellite tidal heating.
    use Galacticus_Nodes
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    use Error_Functions
    use Cosmology_Parameters
    use Galactic_Structure_Densities
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Rotation_Curves
    use Galactic_Structure_Options
    use Vectors
    use Tensors
    use Dark_Matter_Halo_Scales
    implicit none
    class           (nodeComponentSatellite        ), pointer       :: thisSatellite
    type            (treeNode                      ), intent(inout) :: thisNode
    type            (treeNode                      ), pointer       :: hostNode
    class           (cosmologyParametersClass      ), pointer       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    class           (nodeComponentBasic            ), pointer       :: basic
    double precision                                , dimension(3)  :: position                 , velocity
    double precision                                                :: satelliteMass            , parentDensity            , &
         &                                                             parentEnclosedMass       , velocityCircularSatellite, &
         &                                                             radius                   , speed                    , &
         &                                                             timescaleShock           , heatingRateNormalized    , &
         &                                                             orbitalFrequencySatellite, radiusHalfMassSatellite  , &
         &                                                             satelliteHalfMass        , darkMatterFraction
    type            (tensorRank2Dimension3Symmetric)                :: tidalTensor              , tidalTensorPathIntegrated, &
         &                                                             positionTensor
    
    ! Construct required properties of satellite and host.
    hostNode                  => thisNode     %mergesWith               ()
    thisSatellite             => thisNode     %satellite                ()
    satelliteMass             =  thisSatellite%boundMass                ()
    position                  =  thisSatellite%position                 ()
    velocity                  =  thisSatellite%velocity                 ()
    tidalTensorPathIntegrated =  thisSatellite%tidalTensorPathIntegrated()
    radius                    =  Vector_Magnitude                (         position                                            )
    speed                     =  Vector_Magnitude                (         velocity                                            )
    parentDensity             =  Galactic_Structure_Density      (hostNode,position,coordinateSystemCartesian                  )
    parentEnclosedMass        =  Galactic_Structure_Enclosed_Mass(hostNode,radius                                              )
    positionTensor            =  Vector_Outer_Product            (         position                          ,symmetrize=.true.)
    ! Find the universal dark matter fraction.
    cosmologyParameters_      => cosmologyParameters()    
    darkMatterFraction        =  +(                                    &
         &                         +cosmologyParameters_%OmegaMatter() &
         &                         -cosmologyParameters_%OmegaBaryon() &
         &                        )                                    &
         &                       /  cosmologyParameters_%OmegaMatter()
    ! Find the gravitational tidal tensor.
    tidalTensor=                                                                                          &
         & -(gravitationalConstantGalacticus*parentEnclosedMass         /radius**3)*tensorIdentityR2D3Sym &
         & +(gravitationalConstantGalacticus*parentEnclosedMass*3.0d0   /radius**5)*positionTensor        &
         & -(gravitationalConstantGalacticus*parentDensity     *4.0d0*Pi/radius**2)*positionTensor
    ! Find the orbital frequency at the half mass radius of the satellite.
    basic             => thisNode%basic()
    satelliteHalfMass =  0.50d0*darkMatterFraction*min(satelliteMass,basic%mass())    
    radiusHalfMassSatellite  =Galactic_Structure_Radius_Enclosing_Mass(                                           &
         &                                                             thisNode                                 , &
         &                                                             mass                   =satelliteHalfMass, &
         &                                                             componentType          =componentTypeAll , &
         &                                                             massType               =massTypeDark     , &
         &                                                             haloLoaded             =.true.             &
         &                                                            )
    velocityCircularSatellite=Galactic_Structure_Rotation_Curve       (                                           &
         &                                                             thisNode                                 , &
         &                                                             radiusHalfMassSatellite                  , &
         &                                                             componentType          =componentTypeAll , &
         &                                                             massType               =massTypeDark     , &
         &                                                             haloLoaded             =.true.             &
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
       darkMatterHaloScale_      =>  darkMatterHaloScale                (        )
       orbitalFrequencySatellite =  +darkMatterHaloScale_%virialVelocity(thisNode) &
            &                       /darkMatterHaloScale_%virialRadius  (thisNode) &
            &                       *gigaYear                                      &
            &                       *kilo                                          &
            &                       /megaParsec
    end if
    ! Find the shock timescale (i.e. crossing time in the radial direction).
    timescaleShock=megaParsec/kilo/gigaYear*radius/speed
    ! Compute the heating rate.
    heatingRateNormalized=                                          &
         &    satelliteTidalHeatingGnedinEpsilon                    &
         &   /(                                                     &
         &      1.0d0                                               &
         &     +(                                                   &
         &       +timescaleShock                                    &
         &       *orbitalFrequencySatellite                         &
         &      )**2                                                &
         &    )**satelliteTidalHeatingGnedinGamma                   &
         &   /3.0d0                                                 &
         &   *tidalTensor%doubleContract(tidalTensorPathIntegrated) &
         &   *(kilo*gigaYear/megaParsec)**2
    ! Limit the heating rate to be non-negative.
    Satellite_Tidal_Heating_Rate_Gnedin=max(heatingRateNormalized,0.0d0)
    return
  end function Satellite_Tidal_Heating_Rate_Gnedin

end module Tidal_Heating_Rate_Gnedin
