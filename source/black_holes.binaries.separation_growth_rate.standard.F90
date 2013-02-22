!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a black hole binary separation growth rate which follows a modified version of
!% \cite{volonteri_assembly_2003}, including terms for dynamical friction, hardening due to scattering of stars and emission of
!% gravitational waves.

module Black_Hole_Binary_Separations_Standard
  !% Implements a black hole binary initial separation growth rate which follows a modified version of
  !% \cite{volonteri_assembly_2003}, including terms for dynamical friction, hardening due to scattering of stars and emission of
  !% gravitational waves.
  implicit none
  private
  public :: Black_Hole_Binary_Separation_Growth_Rate_Standard_Init
  
  ! Option indicating whether the density of stars in the loss cone should evolve or remain fixed.
  logical :: stellarDensityChangeBinaryMotion

  ! Option indicating whether or not to compute velocity dispersions.
  logical :: blackHoleBinariesComputeVelocityDispersion

contains

  !# <blackHoleBinarySeparationGrowthRateMethod>
  !#  <unitName>Black_Hole_Binary_Separation_Growth_Rate_Standard_Init</unitName>
  !# </blackHoleBinarySeparationGrowthRateMethod>
  subroutine Black_Hole_Binary_Separation_Growth_Rate_Standard_Init(blackHoleBinarySeparationGrowthRateMethod&
       &,Black_Hole_Binary_Separation_Growth_Rate_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: blackHoleBinarySeparationGrowthRateMethod
    procedure(double precision), pointer, intent(inout) :: Black_Hole_Binary_Separation_Growth_Rate_Get

    if (blackHoleBinarySeparationGrowthRateMethod == 'standard') then
       Black_Hole_Binary_Separation_Growth_Rate_Get => Black_Hole_Binary_Separation_Growth_Rate_Standard
       !@ <inputParameter>
       !@   <name>stellarDensityChangeBinaryMotion</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The change in density due to the black hole's motion.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarDensityChangeBinaryMotion',stellarDensityChangeBinaryMotion,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>blackHoleBinariesComputeVelocityDispersion</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not the velocity dispersion of dark matter and stars should be computed using Jeans equation
       !@     in black hole binary hardening calculations. If {\tt false}, then the velocity dispersions are assumed to equal
       !@     the characteristic velocity of dark matter and spheroid.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('blackHoleBinariesComputeVelocityDispersion',blackHoleBinariesComputeVelocityDispersion,defaultValue=.false.)
    end if
    return
  end subroutine Black_Hole_Binary_Separation_Growth_Rate_Standard_Init

  double precision function Black_Hole_Binary_Separation_Growth_Rate_Standard(thisBlackHoleComponent)
    !% Returns an initial separation growth rate for a binary black holes that follows a modified version of
    !% \cite{volonteri_assembly_2003}.
    use Galacticus_Nodes
    use Bessel_Functions
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use Galactic_Structure_Densities
    use Galactic_Structure_Options 
    use Galactic_Structure_Rotation_Curves
    use Galactic_Structure_Rotation_Curve_Gradients
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Velocity_Dispersions
    use Dark_Matter_Halo_Scales
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (nodeComponentBlackHole), intent(inout), pointer :: thisBlackHoleComponent
    type            (treeNode              ),                pointer :: thisNode
    class           (nodeComponentBlackHole),                pointer :: centralBlackHoleComponent
    class           (nodeComponentSpheroid ),                pointer :: thisSpheroidComponent
    double precision                        , parameter              :: hardeningRateDimensionless    =15.0d0
    double precision                        , parameter              :: outerRadiusMultiplier         =10.0d0
    double precision                        , parameter              :: dynamicalFrictionMinimumRadius=0.1d0
    double precision                                                 :: rateScattering, rateGravitationalWaves,&
         & dynamicalFrictionAcceleration, densitySpheroid, densityDarkMatter, dynamicalFrictionXSpheroid,&
         & dynamicalFrictionXDarkMatter, coulombLogarithmSpheroid , coulombLogarithmDarkMatter, densityStellar,&
         & stellarDensityFractionRemaining ,rateScatteringDynamicalFriction, rateScatteringStars, velocityDispersionSpheroid,&
         & velocityDispersionDarkMatter, radiusHardBinary,rotationCurveGradient
    character(len=24)                                                :: message
    
    ! Get the host node.
    thisNode                  => thisBlackHoleComponent%host()
    ! Get the primary (central) black hole of this node.
    centralBlackHoleComponent => thisNode%blackHole(instance=1)
    ! Return a zero separation growth rate if the black hole has non-positive radial position or either black hole (active or
    ! central) has negative mass.
    if     (                                                     &
         &   thisBlackHoleComponent   %radialPosition() <= 0.0d0 &
         &  .or.                                                 &
         &   thisBlackHoleComponent   %mass          () <= 0.0d0 &
         &  .or.                                                 &
         &   centralBlackHoleComponent%mass          () <= 0.0d0 &
         & ) then
       Black_Hole_Binary_Separation_Growth_Rate_Standard=0.0d0
       return
    end if
    ! Get the spheroid component.
    thisSpheroidComponent => thisNode%spheroid()
    ! Compute the velocity dispersion of stars and dark matter.
    if (blackHoleBinariesComputeVelocityDispersion) then
       velocityDispersionSpheroid  =Galactic_Structure_Velocity_Dispersion(                                                                       &
            &                                                                                                    thisNode,                        &
            &                                                              thisBlackHoleComponent%radialPosition(        ),                       &
            &                                                              thisSpheroidComponent %radius        (        )*outerRadiusMultiplier, &
            &                                                              componentTypeSpheroid                                                  &
            &                                                             )
       velocityDispersionDarkMatter=Galactic_Structure_Velocity_Dispersion(                                                                       &
            &                                                                                                    thisNode,                        &
            &                                                              thisBlackHoleComponent%radialPosition(        ),                       &
            &                                                              Dark_Matter_Halo_Virial_Radius       (thisNode)*outerRadiusMultiplier, &
            &                                                              componentTypeDarkHalo                                                  &
            &                                                             )
    else
       velocityDispersionSpheroid  =thisSpheroidComponent%velocity  (        )
       velocityDispersionDarkMatter=Dark_Matter_Halo_Virial_Velocity(thisNode)
    end if
    ! Compute the separation growth rate due to emission of gravitational waves.
    rateGravitationalWaves=-( 256.0d0                                &
         &             *gravitationalConstantGalacticus**3           &
         &             *  thisBlackHoleComponent   %mass       ()    &
         &             *  centralBlackHoleComponent%mass       ()    &
         &             *(                                            &
         &                thisBlackHoleComponent   %mass       ()    &
         &               +centralBlackHoleComponent%mass       ()    &
         &              )                                            &
         &            )                                              &
         &           /(                                              &
         &              5.0d0                                        &
         &             *(speedLight*milli)**5                        &
         &             *  thisBlackHoleComponent%radialPosition()**3 &
         &            )                                              &
         &           /Mpc_per_km_per_s_To_Gyr
    ! Compute the hard binary radius, where shrinking of the binary switches from dynamical friction to hardening due to strong
    ! scattering of individual stars.
    if (velocityDispersionSpheroid > 0.0d0) then
       radiusHardBinary= (                                    &
            &              gravitationalConstantGalacticus    &
            &             *(                                  &
            &                thisBlackHoleComponent   %mass() &
            &               +centralBlackHoleComponent%mass() &
            &              )                                  &
            &            )                                    &
            &           /(                                    &
            &              4.0d0                              &
            &             *velocityDispersionSpheroid**2      &
            &            )
    else
       radiusHardBinary=0.0d0
    end if
    ! First, check if the change in stellar density due to the binary's inward motion is to be computed.
    stellarDensityFractionRemaining=1.0d0
    ! If it does change, we first compute the fraction of that change according to Volonteri et al. (2003) else we set it as the
    ! normal density.
    if (stellarDensityChangeBinaryMotion .and. velocityDispersionSpheroid > 0.0d0)                                                        &
         &  stellarDensityFractionRemaining=(                                                  thisBlackHoleComponent   %radialPosition() &
         &                                   *                                                 velocityDispersionSpheroid**2              &
         &                                   *(4.0d0/3.0d0)/(gravitationalConstantGalacticus*( centralBlackHoleComponent%mass          () &
         &                                                                                    +thisBlackHoleComponent   %mass          () &
         &                                                                                   )                                            &
         &                                   *log                                            (                                            &
         &                                                   gravitationalConstantGalacticus*  thisBlackHoleComponent   %mass          () &
         &                                                                                    /4.0d0                                      &
         &                                                                                    /velocityDispersionSpheroid**2              &
         &                                                                                    /thisBlackHoleComponent   %radialPosition() &
         &                                                                                   )                                            &
         &                                                  )                                                                             &
         &                                  )**2
    ! Limit the density fraction to unity.
    stellarDensityFractionRemaining=min(stellarDensityFractionRemaining,1.0d0)
    ! Compute the stellar density, accounting for any loss.
    densityStellar= Galactic_Structure_Density(thisNode                                             , &
         &                                     [thisBlackHoleComponent%radialPosition(),0.0d0,0.0d0], &
         &                                     coordinateSystem=coordinateSystemCylindrical         , &
         &                                     componentType   =componentTypeSpheroid               , &
         &                                     massType        =massTypeStellar                       &
         &                                    )                                                       &
         &         *stellarDensityFractionRemaining
    ! Compute the hardening rate due to strong scattering of individual stars.
    if (velocityDispersionSpheroid > 0.0d0) then
       rateScatteringStars=-gravitationalConstantGalacticus            &
            &              *densityStellar                             &
            &              *thisBlackHoleComponent%radialPosition()**2 &
            &              *hardeningRateDimensionless                 &
            &              /velocityDispersionSpheroid                 &
            &              /Mpc_per_km_per_s_To_Gyr  
    else
       rateScatteringStars=0.0d0
    end if
    ! Check if the binary has sufficiently large separation that we should compute the rate of hardening due to dynamical friction.
    if (thisBlackHoleComponent%radialPosition() > dynamicalFrictionMinimumRadius*radiusHardBinary) then
       ! Compute the total density, including dark matter.
       densitySpheroid  =Galactic_Structure_Density(thisNode                                             , &
            &                                       [thisBlackHoleComponent%radialPosition(),0.0d0,0.0d0], &
            &                                       coordinateSystem=coordinateSystemCylindrical         , &
            &                                       massType        =massTypeGalactic                      &
            &                                      )                                                              
       densityDarkMatter=Galactic_Structure_Density(thisNode                                             , &
            &                                       [thisBlackHoleComponent%radialPosition(),0.0d0,0.0d0], &
            &                                       coordinateSystem=coordinateSystemCylindrical         , &
            &                                       massType        =massTypeDark                          &
            &                                      )
       ! Compute the Coulomb logarithms for dynamical friction.
       coulombLogarithmSpheroid  = (                                              &
            &                          thisBlackHoleComponent   %radialPosition() &
            &                       *  velocityDispersionSpheroid**2              &
            &                      )                                              &
            &                     /(                                              &
            &                        gravitationalConstantGalacticus              &
            &                       *(                                            &
            &                          thisBlackHoleComponent   %mass          () &
            &                         +centralBlackHoleComponent%mass          () &
            &                        )                                            &
            &                      )
       coulombLogarithmDarkMatter= (                                              &
            &                          thisBlackHoleComponent   %radialPosition() &
            &                       *  velocityDispersionDarkMatter**2            &
            &                      )                                              &
            &                     /(                                              &
            &                        gravitationalConstantGalacticus              &
            &                       *(                                            &
            &                          thisBlackHoleComponent   %mass          () &
            &                         +centralBlackHoleComponent%mass          () &
            &                        )                                            &
            &                      )
       ! Compute the rotation curve of the galaxy and the additional contribution from the active black hole. Add them in
       ! quadrature to get an estimate of the actual orbital speed of the black hole binary.
       ! Precompute the "X" term appearing in the dynamical friction formula.
       if (velocityDispersionSpheroid > 0.0d0) then
          dynamicalFrictionXSpheroid  = Galactic_Structure_Rotation_Curve(                                         &
               &                                                          thisNode                               , &
               &                                                          thisBlackHoleComponent%radialPosition()  &
               &                                                         )                                         &
               &                       /sqrt(2.0d0)                                                                &
               &                       /velocityDispersionSpheroid
       else
          dynamicalFrictionXSpheroid=0.0d0
       end if
       dynamicalFrictionXDarkMatter= Galactic_Structure_Rotation_Curve(                                         &
            &                                                          thisNode                               , &
            &                                                          thisBlackHoleComponent%radialPosition()  &
            &                                                         )                                         &
            &                 /sqrt(2.0d0)                                                                      &
            &                 /velocityDispersionDarkMatter
       ! Compute the acceleration due to dynamical friction.
       dynamicalFrictionAcceleration=-4.0d0                                                                                  & 
            &                        *Pi                                                                                     &
            &                        *gravitationalConstantGalacticus**2                                                     &
            &                        *thisBlackHoleComponent%mass()                                                          &
            &                        *0.5d0                                                                                  &
            &                        /Galactic_Structure_Rotation_Curve(thisNode,thisBlackHoleComponent%radialPosition())**2 &
            &                        /Mpc_per_km_per_s_To_Gyr                                                                &
            &                        *(                                                                                      &
            &                           densitySpheroid                                                                      &
            &                          *log(1.0d0+coulombLogarithmSpheroid**2)                                               &
            &                          *(                                                                                    &
            &                             erf(dynamicalFrictionXSpheroid)                                                    &
            &                            -(                                                                                  &
            &                               2.0d0                                                                            &
            &                              *dynamicalFrictionXSpheroid                                                       &
            &                              /dsqrt(Pi)                                                                        &
            &                              *exp(-dynamicalFrictionXSpheroid**2)                                              &
            &                             )                                                                                  &
            &                           )                                                                                    &
            &                          +densityDarkMatter                                                                    &
            &                          *log(1.0d0+coulombLogarithmDarkMatter**2)                                             &
            &                          *(                                                                                    &
            &                             erf(dynamicalFrictionXDarkMatter)                                                  &
            &                            -(                                                                                  &
            &                               2.0d0                                                                            &
            &                              *dynamicalFrictionXDarkMatter                                                     &
            &                              /dsqrt(Pi)                                                                        &
            &                              *exp(-dynamicalFrictionXDarkMatter**2)                                            &
            &                             )                                                                                  &
            &                           )                                                                                    &
            &                         )
       ! Compute the radial inflow velocity due to dynamical friction.
       rotationCurveGradient          =(                                                                                              &
            &                            Galactic_Structure_Rotation_Curve         (thisNode,thisBlackHoleComponent%radialPosition()) &
            &                           +                                                    thisBlackHoleComponent%radialPosition()  &
            &                           *Galactic_Structure_Rotation_Curve_Gradient(thisNode,thisBlackHoleComponent%radialPosition()) &
            &                          )
       if (rotationCurveGradient == 0.0d0) then
          call Galacticus_Display_Indent('dynamical friction calculation report')
          write (message,'(a,i12)  ') 'nodeIndex = ',thisNode%index()
          write (message,'(a,e12.6)') '     V(r) = ',Galactic_Structure_Rotation_Curve         (thisNode,thisBlackHoleComponent%radialPosition())
          call Galacticus_Display_Message(trim(message))
          write (message,'(a,e12.6)') '        r = ',                                                    thisBlackHoleComponent%radialPosition()
          call Galacticus_Display_Message(trim(message))
          write (message,'(a,e12.6)') ' dV(r)/dr = ',Galactic_Structure_Rotation_Curve_Gradient(thisNode,thisBlackHoleComponent%radialPosition())
          call Galacticus_Display_Message(trim(message))
          write (message,'(a,e12.6)') '   a_{df} = ',dynamicalFrictionAcceleration
          call Galacticus_Display_Message(trim(message))
          call Galacticus_Display_Unindent('done')
          call Galacticus_Error_Report('Black_Hole_Binary_Separation_Growth_Rate_Standard','rotation curve gradient is zero')
       end if
       rateScatteringDynamicalFriction= 2.0d0                                   &
            &                          *dynamicalFrictionAcceleration           &
            &                          *thisBlackHoleComponent%radialPosition() &
            &                          /rotationCurveGradient
       ! Take the most negative rate from dynamical friction or scattering of stars.
       rateScattering=min(rateScatteringDynamicalFriction,rateScatteringStars)
    else
       ! For sufficiently small radii assume that scattering of individual stars always dominates over dynamical friction.
       rateScattering=rateScatteringStars
    end if
    ! Sum the two contributions to the radial growth rate.
    Black_Hole_Binary_Separation_Growth_Rate_Standard=rateScattering+rateGravitationalWaves 
    return
  end function Black_Hole_Binary_Separation_Growth_Rate_Standard
     
end module Black_Hole_Binary_Separations_Standard
