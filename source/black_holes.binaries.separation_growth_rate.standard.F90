!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

  !% Implements a black hole binary separation growth class which follows a modified version of \cite{volonteri_assembly_2003},
  !% including terms for dynamical friction, hardening due to scattering of stars and emission of gravitational waves.

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !# <blackHoleBinarySeparationGrowthRate name="blackHoleBinarySeparationGrowthRateStandard">
  !#  <description>A black hole binary separation growth class which follows a modified version of \cite{volonteri_assembly_2003}, including terms for dynamical friction, hardening due to scattering of stars and emission of gravitational waves.</description>
  !# </blackHoleBinarySeparationGrowthRate>
  type, extends(blackHoleBinarySeparationGrowthRateClass) :: blackHoleBinarySeparationGrowthRateStandard
     !% A black hole binary separation growth class which follows a modified version of \cite{volonteri_assembly_2003}, including
     !% terms for dynamical friction, hardening due to scattering of stars and emission of gravitational waves.
     private
     logical                                    :: stellarDensityChangeBinaryMotion, computeVelocityDispersion
     class  (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::               standardDestructor
     procedure :: growthRate => standardGrowthRate
  end type blackHoleBinarySeparationGrowthRateStandard

  interface blackHoleBinarySeparationGrowthRateStandard
     !% Constructors for the {\normalfont \ttfamily standard} black hole binary recoil class.
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface blackHoleBinarySeparationGrowthRateStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily standard} black hole binary separation growth rate class which takes a parameter
    !% set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (blackHoleBinarySeparationGrowthRateStandard)                :: self
    type   (inputParameters                            ), intent(inout) :: parameters
    class  (darkMatterHaloScaleClass                   ), pointer       :: darkMatterHaloScale_
    logical                                                             :: stellarDensityChangeBinaryMotion, computeVelocityDispersion

    !# <inputParameter>
    !#   <name>stellarDensityChangeBinaryMotion</name>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>The change in density due to the black hole's motion.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>computeVelocityDispersion</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not the velocity dispersion of dark matter and stars should be computed using Jeans equation
    !#      in black hole binary hardening calculations. If {\normalfont \ttfamily false}, then the velocity dispersions are assumed to equal
    !#      the characteristic velocity of dark matter and spheroid.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=blackHoleBinarySeparationGrowthRateStandard(stellarDensityChangeBinaryMotion,computeVelocityDispersion,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function standardConstructorParameters

  function standardConstructorInternal(stellarDensityChangeBinaryMotion,computeVelocityDispersion,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily standard} black hole binary separation growth class.
    implicit none
    type   (blackHoleBinarySeparationGrowthRateStandard)                        :: self
    class  (darkMatterHaloScaleClass                   ), intent(in   ), target :: darkMatterHaloScale_
    logical                                             , intent(in   )         :: stellarDensityChangeBinaryMotion, computeVelocityDispersion
    !# <constructorAssign variables="stellarDensityChangeBinaryMotion,computeVelocityDispersion,*darkMatterHaloScale_"/>

    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !% Destructor for the {\normalfont \ttfamily standard} black hole binary separation growth class.
    implicit none
    type(blackHoleBinarySeparationGrowthRateStandard), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine standardDestructor

  double precision function standardGrowthRate(self,blackHole)
    !% Returns an initial separation growth rate for a binary black holes that follows a modified version of
    !% \cite{volonteri_assembly_2003}.
    use :: Galactic_Structure_Densities               , only : Galactic_Structure_Density
    use :: Galactic_Structure_Options                 , only : componentTypeDarkHalo                     , componentTypeSpheroid          , coordinateSystemCylindrical, massTypeDark, &
          &                                                    massTypeGalactic                          , massTypeStellar
    use :: Galactic_Structure_Rotation_Curve_Gradients, only : Galactic_Structure_Rotation_Curve_Gradient
    use :: Galactic_Structure_Rotation_Curves         , only : Galactic_Structure_Rotation_Curve
    use :: Galactic_Structure_Velocity_Dispersions    , only : Galactic_Structure_Velocity_Dispersion
    use :: Galacticus_Display                         , only : Galacticus_Display_Indent                 , Galacticus_Display_Message     , Galacticus_Display_Unindent
    use :: Galacticus_Error                           , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                           , only : nodeComponentBlackHole                    , nodeComponentSpheroid          , treeNode
    use :: Numerical_Constants_Astronomical           , only : Mpc_per_km_per_s_To_Gyr                   , gravitationalConstantGalacticus
    use :: Numerical_Constants_Math                   , only : Pi
    use :: Numerical_Constants_Physical               , only : speedLight
    use :: Numerical_Constants_Prefixes               , only : milli
    implicit none
    class           (blackHoleBinarySeparationGrowthRateStandard), intent(inout) :: self
    class           (nodeComponentBlackHole                     ), intent(inout) :: blackHole
    type            (treeNode                                   ), pointer       :: node
    class           (nodeComponentBlackHole                     ), pointer       :: blackHoleCentral
    class           (nodeComponentSpheroid                      ), pointer       :: spheroid
    double precision                                             , parameter     :: hardeningRateDimensionless     =15.0d0
    double precision                                             , parameter     :: outerRadiusMultiplier          =10.0d0
    double precision                                             , parameter     :: dynamicalFrictionMinimumRadius =0.1d0
    double precision                                                             :: coulombLogarithmDarkMatter            , coulombLogarithmSpheroid       , &
         &                                                                          densityDarkMatter                     , densitySpheroid                , &
         &                                                                          densityStellar                        , dynamicalFrictionAcceleration  , &
         &                                                                          dynamicalFrictionXDarkMatter          , dynamicalFrictionXSpheroid     , &
         &                                                                          radiusHardBinary                      , rateGravitationalWaves         , &
         &                                                                          rateScattering                        , rateScatteringDynamicalFriction, &
         &                                                                          rateScatteringStars                   , rotationCurveGradient          , &
         &                                                                          stellarDensityFractionRemaining       , velocityDispersionDarkMatter   , &
         &                                                                          velocityDispersionSpheroid
    character       (len=24                                     )                :: message

    ! Get the host node.
    node             => blackHole%host()
    ! Get the primary (central) black hole of this node.
    blackHoleCentral => node%blackHole(instance=1)
    ! Return a zero separation growth rate if the black hole has non-positive radial position or either black hole (active or
    ! central) has negative mass.
    if     (                                        &
         &   blackHole   %radialPosition() <= 0.0d0 &
         &  .or.                                    &
         &   blackHole   %mass          () <= 0.0d0 &
         &  .or.                                    &
         &   blackHoleCentral%mass      () <= 0.0d0 &
         & ) then
       standardGrowthRate=0.0d0
       return
    end if
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Compute the velocity dispersion of stars and dark matter.
    if (self%computeVelocityDispersion) then
       velocityDispersionSpheroid  =Galactic_Structure_Velocity_Dispersion(                                                                           &
            &                                                                                                            node,                        &
            &                                                              blackHole                     %radialPosition(    ),                       &
            &                                                              spheroid                      %radius        (    )*outerRadiusMultiplier, &
            &                                                              componentTypeSpheroid                                                    , &
            &                                                              massTypeStellar                                                            &
            &                                                             )
       velocityDispersionDarkMatter=Galactic_Structure_Velocity_Dispersion(                                                                           &
            &                                                                                                            node,                        &
            &                                                              blackHole                     %radialPosition(    ),                       &
            &                                                              self     %darkMatterHaloScale_%virialRadius  (node)*outerRadiusMultiplier, &
            &                                                              componentTypeDarkHalo                                                    , &
            &                                                              massTypeDark                                                               &
            &                                                             )
    else
       velocityDispersionSpheroid  =spheroid                     %velocity      (    )
       velocityDispersionDarkMatter=self    %darkMatterHaloScale_%virialVelocity(node)
    end if
    ! Compute the separation growth rate due to emission of gravitational waves.
    rateGravitationalWaves=-(                                        &
         &                   +256.0d0                                &
         &                   *gravitationalConstantGalacticus**3     &
         &                   *  blackHole       %mass          ()    &
         &                   *  blackHoleCentral%mass          ()    &
         &                   *(                                      &
         &                     +blackHole       %mass          ()    &
         &                     +blackHoleCentral%mass          ()    &
         &                    )                                      &
         &                  )                                        &
         &                 /(                                        &
         &                   +5.0d0                                  &
         &                   *(speedLight*milli)**5                  &
         &                   *  blackHole       %radialPosition()**3 &
         &                  )                                        &
         &                 /Mpc_per_km_per_s_To_Gyr
    ! Compute the hard binary radius, where shrinking of the binary switches from dynamical friction to hardening due to strong
    ! scattering of individual stars.
    if (velocityDispersionSpheroid > 0.0d0) then
       radiusHardBinary=+(                                 &
            &             +gravitationalConstantGalacticus &
            &             *(                               &
            &               +blackHole       %mass()       &
            &               +blackHoleCentral%mass()       &
            &              )                               &
            &            )                                 &
            &           /(                                 &
            &             +4.0d0                           &
            &             *velocityDispersionSpheroid**2   &
            &            )
    else
       radiusHardBinary=0.0d0
    end if
    ! First, check if the change in stellar density due to the binary's inward motion is to be computed.
    stellarDensityFractionRemaining=1.0d0
    ! If it does change, we first compute the fraction of that change according to Volonteri et al. (2003) else we set it as the
    ! normal density.
    if (self%stellarDensityChangeBinaryMotion .and. velocityDispersionSpheroid > 0.0d0)                                          &
         &  stellarDensityFractionRemaining=(+                                                 blackHole       %radialPosition() &
         &                                   *                                                 velocityDispersionSpheroid**2     &
         &                                   *(4.0d0/3.0d0)/(gravitationalConstantGalacticus*(+blackHoleCentral%mass          () &
         &                                                                                    +blackHole       %mass          () &
         &                                                                                   )                                   &
         &                                   *log                                            (                                   &
         &                                                   gravitationalConstantGalacticus*  blackHole       %mass          () &
         &                                                                                    /4.0d0                             &
         &                                                                                    /velocityDispersionSpheroid**2     &
         &                                                                                    /blackHole       %radialPosition() &
         &                                                                                   )                                   &
         &                                                  )                                                                    &
         &                                  )**2
    ! Limit the density fraction to unity.
    stellarDensityFractionRemaining=min(stellarDensityFractionRemaining,1.0d0)
    ! Compute the stellar density, accounting for any loss.
    densityStellar= Galactic_Structure_Density(node                                        , &
         &                                     [blackHole%radialPosition(),0.0d0,0.0d0]    , &
         &                                     coordinateSystem=coordinateSystemCylindrical, &
         &                                     componentType   =componentTypeSpheroid      , &
         &                                     massType        =massTypeStellar              &
         &                                    )                                              &
         &         *stellarDensityFractionRemaining
    ! Compute the hardening rate due to strong scattering of individual stars.
    if (velocityDispersionSpheroid > 0.0d0) then
       rateScatteringStars=-gravitationalConstantGalacticus &
            &              *densityStellar                  &
            &              *blackHole%radialPosition()**2   &
            &              *hardeningRateDimensionless      &
            &              /velocityDispersionSpheroid      &
            &              /Mpc_per_km_per_s_To_Gyr
    else
       rateScatteringStars=0.0d0
    end if
    ! Check if the binary has sufficiently large separation that we should compute the rate of hardening due to dynamical friction.
    if (blackHole%radialPosition() > dynamicalFrictionMinimumRadius*radiusHardBinary) then
       ! Compute the total density, including dark matter.
       densitySpheroid  =Galactic_Structure_Density(node                                        , &
            &                                       [blackHole%radialPosition(),0.0d0,0.0d0]    , &
            &                                       coordinateSystem=coordinateSystemCylindrical, &
            &                                       massType        =massTypeGalactic             &
            &                                      )
       densityDarkMatter=Galactic_Structure_Density(node                                        , &
            &                                       [blackHole%radialPosition(),0.0d0,0.0d0]    , &
            &                                       coordinateSystem=coordinateSystemCylindrical, &
            &                                       massType        =massTypeDark                 &
            &                                      )
       ! Compute the Coulomb logarithms for dynamical friction.
       coulombLogarithmSpheroid  = (                                     &
            &                          blackHole       %radialPosition() &
            &                       *  velocityDispersionSpheroid**2     &
            &                      )                                     &
            &                     /(                                     &
            &                        gravitationalConstantGalacticus     &
            &                       *(                                   &
            &                          blackHole       %mass          () &
            &                         +blackHoleCentral%mass          () &
            &                        )                                   &
            &                      )
       coulombLogarithmDarkMatter= (                                     &
            &                          blackHole       %radialPosition() &
            &                       *  velocityDispersionDarkMatter**2   &
            &                      )                                     &
            &                     /(                                     &
            &                        gravitationalConstantGalacticus     &
            &                       *(                                   &
            &                          blackHole       %mass          () &
            &                         +blackHoleCentral%mass          () &
            &                        )                                   &
            &                      )
       ! Compute the rotation curve of the galaxy and the additional contribution from the active black hole. Add them in
       ! quadrature to get an estimate of the actual orbital speed of the black hole binary.
       ! Precompute the "X" term appearing in the dynamical friction formula.
       if (velocityDispersionSpheroid > 0.0d0) then
          dynamicalFrictionXSpheroid  = Galactic_Structure_Rotation_Curve(                            &
               &                                                          node                      , &
               &                                                          blackHole%radialPosition()  &
               &                                                         )                            &
               &                       /sqrt(2.0d0)                                                   &
               &                       /velocityDispersionSpheroid
       else
          dynamicalFrictionXSpheroid=0.0d0
       end if
       dynamicalFrictionXDarkMatter= Galactic_Structure_Rotation_Curve(                            &
            &                                                          node                      , &
            &                                                          blackHole%radialPosition()  &
            &                                                         )                            &
            &                 /sqrt(2.0d0)                                                         &
            &                 /velocityDispersionDarkMatter
       ! Compute the acceleration due to dynamical friction.
       dynamicalFrictionAcceleration=-4.0d0                                                                 &
            &                        *Pi                                                                    &
            &                        *gravitationalConstantGalacticus**2                                    &
            &                        *blackHole%mass()                                                      &
            &                        *0.5d0                                                                 &
            &                        /Galactic_Structure_Rotation_Curve(node,blackHole%radialPosition())**2 &
            &                        /Mpc_per_km_per_s_To_Gyr                                               &
            &                        *(                                                                     &
            &                           densitySpheroid                                                     &
            &                          *log(1.0d0+coulombLogarithmSpheroid**2)                              &
            &                          *(                                                                   &
            &                             erf(dynamicalFrictionXSpheroid)                                   &
            &                            -(                                                                 &
            &                               2.0d0                                                           &
            &                              *dynamicalFrictionXSpheroid                                      &
            &                              /sqrt(Pi)                                                        &
            &                              *exp(-dynamicalFrictionXSpheroid**2)                             &
            &                             )                                                                 &
            &                           )                                                                   &
            &                          +densityDarkMatter                                                   &
            &                          *log(1.0d0+coulombLogarithmDarkMatter**2)                            &
            &                          *(                                                                   &
            &                             erf(dynamicalFrictionXDarkMatter)                                 &
            &                            -(                                                                 &
            &                               2.0d0                                                           &
            &                              *dynamicalFrictionXDarkMatter                                    &
            &                              /sqrt(Pi)                                                        &
            &                              *exp(-dynamicalFrictionXDarkMatter**2)                           &
            &                             )                                                                 &
            &                           )                                                                   &
            &                         )
       ! Compute the radial inflow velocity due to dynamical friction.
       rotationCurveGradient          =(                                                                             &
            &                           +Galactic_Structure_Rotation_Curve         (node,blackHole%radialPosition()) &
            &                           +                                                blackHole%radialPosition()  &
            &                           *Galactic_Structure_Rotation_Curve_Gradient(node,blackHole%radialPosition()) &
            &                          )
       if (rotationCurveGradient == 0.0d0) then
          call Galacticus_Display_Indent('dynamical friction calculation report')
          write (message,'(a,i12)  ') 'nodeIndex = ',node%index()
          write (message,'(a,e12.6)') '     V(r) = ',Galactic_Structure_Rotation_Curve         (node,blackHole%radialPosition())
          call Galacticus_Display_Message(trim(message))
          write (message,'(a,e12.6)') '        r = ',                                                blackHole%radialPosition()
          call Galacticus_Display_Message(trim(message))
          write (message,'(a,e12.6)') ' dV(r)/dr = ',Galactic_Structure_Rotation_Curve_Gradient(node,blackHole%radialPosition())
          call Galacticus_Display_Message(trim(message))
          write (message,'(a,e12.6)') '   a_{df} = ',dynamicalFrictionAcceleration
          call Galacticus_Display_Message(trim(message))
          call Galacticus_Display_Unindent('done')
          call Galacticus_Error_Report('rotation curve gradient is zero'//{introspection:location})
       end if
       rateScatteringDynamicalFriction=+2.0d0                         &
            &                          *dynamicalFrictionAcceleration &
            &                          *blackHole%radialPosition()    &
            &                          /rotationCurveGradient
       ! Take the most negative rate from dynamical friction or scattering of stars.
       rateScattering=min(rateScatteringDynamicalFriction,rateScatteringStars)
    else
       ! For sufficiently small radii assume that scattering of individual stars always dominates over dynamical friction.
       rateScattering=rateScatteringStars
    end if
    ! Sum the two contributions to the radial growth rate.
    standardGrowthRate=rateScattering+rateGravitationalWaves
    return
  end function standardGrowthRate
