!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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

  double precision function Black_Hole_Binary_Separation_Growth_Rate_Standard(thisNode)
    !% Returns an initial separation growth rate for a binary black holes that follows a modified version of
    !% \cite{volonteri_assembly_2003}.
    use Tree_Nodes
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
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, parameter              :: hardeningRateDimensionless    =15.0d0
    double precision, parameter              :: outerRadiusMultiplier         =10.0d0
    double precision, parameter              :: dynamicalFrictionMinimumRadius=0.1d0
    double precision                         :: rateScattering, rateGravitationalWaves, dynamicalFrictionAcceleration,&
         & densitySpheroid, densityDarkMatter, dynamicalFrictionXSpheroid, dynamicalFrictionXDarkMatter, coulombLogarithmSpheroid&
         &, coulombLogarithmDarkMatter, rotationCurveTotal, densityStellar, stellarDensityFractionRemaining &
         &,rateScatteringDynamicalFriction, rateScatteringStars, velocityDispersionSpheroid, velocityDispersionDarkMatter,&
         & radiusHardBinary

    ! Return a zero separation growth rate if the black hole has non-positive radial position or either black hole (active or
    ! central) has negative mass.
    if     (                                                                    &
         &   Tree_Node_Black_Hole_Radial_Position(thisNode           ) <= 0.0d0 &
         &  .or.                                                                &
         &   Tree_Node_Black_Hole_Mass           (thisNode           ) <= 0.0d0 &
         &  .or.                                                                &
         &   Tree_Node_Black_Hole_Mass           (thisNode,instance=1) <= 0.0d0 &
         & ) then
       Black_Hole_Binary_Separation_Growth_Rate_Standard=0.0d0
       return
    end if

    ! Compute the velocity dispersion of stars and dark matter.
    if (blackHoleBinariesComputeVelocityDispersion) then
       velocityDispersionSpheroid  =Galactic_Structure_Velocity_Dispersion(                                                                      &
            &                                                                                                   thisNode,                        &
            &                                                              Tree_Node_Black_Hole_Radial_Position(thisNode),                       &
            &                                                              Tree_Node_Spheroid_Radius           (thisNode)*outerRadiusMultiplier, &
            &                                                              componentTypeSpheroid                                                 &
            &                                                             )
       velocityDispersionDarkMatter=Galactic_Structure_Velocity_Dispersion(                                                                      &
            &                                                                                                   thisNode,                        &
            &                                                              Tree_Node_Black_Hole_Radial_Position(thisNode),                       &
            &                                                              Dark_Matter_Halo_Virial_Radius      (thisNode)*outerRadiusMultiplier, &
            &                                                              componentTypeDarkHalo                                                 &
            &                                                             )
    else
       velocityDispersionSpheroid  =Tree_Node_Spheroid_Velocity     (thisNode)
       velocityDispersionDarkMatter=Dark_Matter_Halo_Virial_Velocity(thisNode)
    end if
    ! Compute the separation growth rate due to emission of gravitational waves.
    rateGravitationalWaves=-( 256.0d0                                                  &
         &             *gravitationalConstantGalacticus**3                             &
         &             *  Tree_Node_Black_Hole_Mass           (thisNode           )    &
         &             *  Tree_Node_Black_Hole_Mass           (thisNode,instance=1)    &
         &             *(                                                              &
         &                Tree_Node_Black_Hole_Mass           (thisNode           )    &
         &               +Tree_Node_Black_Hole_Mass           (thisNode,instance=1)    &
         &              )                                                              &
         &            )                                                                &
         &           /(                                                                &
         &              5.0d0                                                          &
         &             *(speedLight*milli)**5                                          &
         &             *  Tree_Node_Black_Hole_Radial_Position(thisNode           )**3 &
         &            )                                                                &
         &           /Mpc_per_km_per_s_To_Gyr
    ! Compute the hard binary radius, where shrinking of the binary switches from dynamical friction to hardening due to strong
    ! scattering of individual stars.
    radiusHardBinary= (                                                         &
         &              gravitationalConstantGalacticus                         &
         &             *(                                                       &
         &                Tree_Node_Black_Hole_Mass       (thisNode           ) &
         &               +Tree_Node_Black_Hole_Mass       (thisNode,instance=1) &
         &              )                                                       &
         &            )                                                         &
         &           /(                                                         &
         &              4.0d0                                                   &
         &             *velocityDispersionSpheroid**2                           &
         &            )
    
    ! First, check if the change in stellar density due to the binary's inward motion is to be computed.
    stellarDensityFractionRemaining=1.0d0
    ! If it does change, we first compute the fraction of that change according to Volonteri et al. (2003) else we set it as the
    ! normal density.
    if (stellarDensityChangeBinaryMotion)                                                                                                                   &
         &  stellarDensityFractionRemaining=(                                                  Tree_Node_Black_Hole_Radial_Position(thisNode           )    &
         &                                   *                                                 velocityDispersionSpheroid**2                                &
         &                                   *(4.0d0/3.0d0)/(gravitationalConstantGalacticus*( Tree_Node_Black_Hole_Mass           (thisNode,instance=1)    &
         &                                                                                    +Tree_Node_Black_Hole_Mass           (thisNode           )    &
         &                                                                                   )                                                              &
         &                                   *log                                            (                                                              &
         &                                                   gravitationalConstantGalacticus*  Tree_Node_Black_Hole_Mass           (thisNode           )    &
         &                                                                                    /4.0d0                                                        &
         &                                                                                    /velocityDispersionSpheroid**2                                &
         &                                                                                    /Tree_Node_Black_Hole_Radial_Position(thisNode           )    &
         &                                                                                   )                                                              &
         &                                                  )                                                                                               &
         &                                  )**2
    ! Limit the density fraction to unity.
    stellarDensityFractionRemaining=min(stellarDensityFractionRemaining,1.0d0)

    ! Compute the stellar density, accounting for any loss.
    densityStellar= Galactic_Structure_Density(                                       thisNode               &
         &                                     ,[Tree_Node_Black_Hole_Radial_Position(thisNode),0.0d0,0.0d0] &
         &                                     ,coordinateSystem=coordinateSystemSpherical                   &
         &                                     ,massType        =massTypeStellar                             &
         &                                     ,componentType   =componentTypeSpheroid                       &
         &                                    )                                                              &
         &         *stellarDensityFractionRemaining
 
    ! Compute the hardening rate due to strong scattering of individual stars.
    rateScatteringStars=-gravitationalConstantGalacticus                   &
         &              *densityStellar                                    &
         &              *Tree_Node_Black_Hole_Radial_Position(thisNode)**2 &
         &              *hardeningRateDimensionless                        &
         &              /velocityDispersionSpheroid                        &
         &              /Mpc_per_km_per_s_To_Gyr  

    ! Check if the binary has sufficiently large separation that we should compute the rate of hardening due to dynamical friction.
    if (Tree_Node_Black_Hole_Radial_Position(thisNode) > dynamicalFrictionMinimumRadius*radiusHardBinary) then
       ! Compute the total density, including dark matter.
       densitySpheroid  =Galactic_Structure_Density(                                       thisNode               &
            &                                       ,[Tree_Node_Black_Hole_Radial_Position(thisNode),0.0d0,0.0d0] &
            &                                       ,coordinateSystem=coordinateSystemSpherical                   &
            &                                       ,massType        =massTypeGalactic                            &
            &                                      )                                                              
       densityDarkMatter=Galactic_Structure_Density(                                       thisNode               &
            &                                       ,[Tree_Node_Black_Hole_Radial_Position(thisNode),0.0d0,0.0d0] &
            &                                       ,coordinateSystem=coordinateSystemSpherical                   &
            &                                       ,massType        =massTypeDark                                &
            &                                      )
       ! Compute the Coulomb logarithms for dynamical friction.
       coulombLogarithmSpheroid  = (                                                                &
            &                          Tree_Node_Black_Hole_Radial_Position(thisNode           )    &
            &                       *  velocityDispersionSpheroid**2                                &
            &                      )                                                                &
            &                     /(                                                                &
            &                        gravitationalConstantGalacticus                                &
            &                       *(                                                              &
            &                          Tree_Node_Black_Hole_Mass           (thisNode           )    &
            &                         +Tree_Node_Black_Hole_Mass           (thisNode,instance=1)    &
            &                        )                                                              &
            &                      )
       coulombLogarithmDarkMatter= (                                                                &
            &                          Tree_Node_Black_Hole_Radial_Position(thisNode           )    &
            &                       *  velocityDispersionDarkMatter**2                              &
            &                      )                                                                &
            &                     /(                                                                &
            &                        gravitationalConstantGalacticus                                &
            &                       *(                                                              &
            &                          Tree_Node_Black_Hole_Mass           (thisNode           )    &
            &                         +Tree_Node_Black_Hole_Mass           (thisNode,instance=1)    &
            &                        )                                                              &
            &                      )
       ! Compute the rotation curve of the galaxy and the additional contribution from the active black hole. Add them in
       ! quadrature to get an estimate of the actual orbital speed of the black hole binary.
       ! Precompute the "X" term appearing in the dynamical friction formula.
       dynamicalFrictionXSpheroid  = Galactic_Structure_Rotation_Curve(                                                &
            &                                                                                                thisNode  &
            &                                                          ,Tree_Node_Black_Hole_Radial_Position(thisNode) &
            &                                                         )                                                &
            &                       /dsqrt(2.0d0)                                                                      &
            &                       /velocityDispersionSpheroid
       dynamicalFrictionXDarkMatter= Galactic_Structure_Rotation_Curve(                                                &
            &                                                                                                thisNode  &
            &                                                          ,Tree_Node_Black_Hole_Radial_Position(thisNode) &
            &                                                         )                                                &
            &                 /dsqrt(2.0d0)                                                                            &
            &                 /velocityDispersionDarkMatter
       ! Compute the acceleration due to dynamical friction.
       dynamicalFrictionAcceleration=-4.0d0                                                                                         & 
            &                        *Pi                                                                                            &
            &                        *gravitationalConstantGalacticus**2                                                            &
            &                        *Tree_Node_Black_Hole_Mass(thisNode)                                                           &
            &                        *0.5d0                                                                                         &
            &                        /Galactic_Structure_Rotation_Curve(thisNode,Tree_Node_Black_Hole_Radial_Position(thisNode))**2 &
            &                        /Mpc_per_km_per_s_To_Gyr                                                                       &
            &                        *(                                                                                             &
            &                           densitySpheroid                                                                             &
            &                          *log(1.0d0+coulombLogarithmSpheroid**2)                                                      &
            &                          *(                                                                                           &
            &                             erf(dynamicalFrictionXSpheroid)                                                           &
            &                            -(                                                                                         &
            &                               2.0d0                                                                                   &
            &                              *dynamicalFrictionXSpheroid                                                              &
            &                              /dsqrt(Pi)                                                                               &
            &                              *exp(-dynamicalFrictionXSpheroid**2)                                                     &
            &                             )                                                                                         &
            &                           )                                                                                           &
            &                          +densityDarkMatter                                                                           &
            &                          *log(1.0d0+coulombLogarithmDarkMatter**2)                                                    &
            &                          *(                                                                                           &
            &                             erf(dynamicalFrictionXDarkMatter)                                                         &
            &                            -(                                                                                         &
            &                               2.0d0                                                                                   &
            &                              *dynamicalFrictionXDarkMatter                                                            &
            &                              /dsqrt(Pi)                                                                               &
            &                              *exp(-dynamicalFrictionXDarkMatter**2)                                                   &
            &                             )                                                                                         &
            &                           )                                                                                           &
            &                         )
       ! Compute the radial inflow velocity due to dynamical friction.
       rateScatteringDynamicalFriction= 2.0d0                                                                                                 &
            &                          *dynamicalFrictionAcceleration                                                                         &
            &                          *Tree_Node_Black_Hole_Radial_Position(thisNode)                                                        &
            &                          /(                                                                                                     &
            &                             Galactic_Structure_Rotation_Curve         (thisNode,Tree_Node_Black_Hole_Radial_Position(thisNode)) &
            &                            +Tree_Node_Black_Hole_Radial_Position      (thisNode                                               ) &
            &                            *Galactic_Structure_Rotation_Curve_Gradient(thisNode,Tree_Node_Black_Hole_Radial_Position(thisNode)) &
            &                           )
       
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
