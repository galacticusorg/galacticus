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

!!{
Contains a program which tests heated dark matter profile implementations.
!!}

program Test_Dark_Matter_Profiles_Heated
  !!{
  Tests heated dark matter profile implementations. An isothermal dark matter halo is used, since analytic solutions are available
  for this case. Specifically, the initial radius in the unheated profile is given by $r_\mathrm{i}(r) = (r_\mathrm{h}^4/4 r^2 +
  r_\mathrm{h}^2)^{1/2}-r_\mathrm{h}^2/2 r$ where $r$ is the radius in the heated profile, and $r_\mathrm{h}=(\mathrm{G}
  M_\mathrm{v} / 2 Q r_\mathrm{v})^{1/2}$ is a characteristic heating radius. Here, $M_\mathrm{v}$, and $r_\mathrm{v}$ are the
  virial mass and radius of the halo respectively, and $Q r_\mathrm{i}^2$ is the specific heat input to the density profile, with
  $Q$ assumed to be a constant (as expected for tidal heating). Assuming no shell crossing, the enclosed mass in the final profile
  is simply $M(r) = M_\mathrm{v} r_\mathrm{i}(r)/r_\mathrm{v}$, from which the density of the final profile is found as $\rho(r) =
  (4 \pi r^2)^{-1} \mathrm{d} M(r) / \mathrm{d} r$. The velocity dispersion can also be found analytically in this
  case. Integrating the Jeans equation in the variables of the initial profile gives a result:
  \begin{eqnarray}
    \rho(r) \sigma^2(r) &=& \int_{r_\mathrm{i}(r)}^{r_\mathrm{h}} \mathrm{d}r^\prime_\mathrm{i} \frac{\mathrm{G} M(r^\prime_\mathrm{i})}{r^{\prime 2}_\mathrm{i}} \rho_\mathrm{i}(r^\prime_\mathrm{i}) \left(\frac{r^\prime_\mathrm{i}}{r}\right)^4 \nonumber \\
    &=& \frac{10}{6} \frac{\mathrm{G} M(r_\mathrm{h}) \rho_\mathrm{i}(r_\mathrm{h}}{r_\mathrm{h}} - \frac{(-3+18 y^4 - 6 y^6 + y^8 - 24 y^2 \log(y))}{6} \frac{\mathrm{G} M(r_\mathrm{i}) \rho_\mathrm{i}(r_\mathrm{i}}{r_\mathrm{i}},
  \end{eqnarray}
  where $y=r_\mathrm{i}/r_\mathrm{h}$.  
  !!}
  use :: Coordinates                     , only : coordinateSpherical
  use :: Cosmology_Parameters            , only : cosmologyParametersSimple
  use :: Cosmology_Functions             , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Particles           , only : darkMatterParticleSelfInteractingDarkMatter                      , darkMatterParticleCDM
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Virial_Density_Contrast         , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOHeated                                       , darkMatterProfileDMOHeatedMonotonic, darkMatterProfileDMOIsothermal          , darkMatterProfileHeatingTidal      , &
       &                                          darkMatterProfileDMOClass
  use :: Display                         , only : displayVerbositySet                                              , verbosityLevelStandard
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Galacticus_Nodes                , only : nodeClassHierarchyFinalize                                       , nodeClassHierarchyInitialize       , nodeComponentBasic                      , nodeComponentSatellite             , &
       &                                          treeNode
  use :: ISO_Varying_String              , only : varying_string
  use :: Input_Parameters                , only : inputParameters
  use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
  use :: Numerical_Constants_Math        , only : Pi
  use :: Mass_Distributions              , only : massDistributionSphericalHeated                                  , massDistributionClass              , massDistributionSpherical               , massDistributionHeatingClass       , &
       &                                          nonAnalyticSolversNumerical                                      , kinematicsDistributionHeated       , massDistributionSphericalHeatedMonotonic, kinematicsDistributionCollisionless, &
       &                                          kinematicsDistributionClass                                                                                                                                                          , &
       &                                          nonAnalyticSolversFallThroughDMO => nonAnalyticSolversFallThrough                                                                                                                    , &
       &                                          nonAnalyticSolversNumericalDMO   => nonAnalyticSolversNumerical
  use :: Unit_Tests                      , only : Assert                                                           , Unit_Tests_Begin_Group             , Unit_Tests_End_Group                    , Unit_Tests_Finish
  implicit none
  double precision                                                                , parameter    :: time                                        =13.8d+00
  double precision                                                                , parameter    :: massVirial                                  = 1.0d+10
  double precision                                                                , parameter    :: heatingSpecific                             = 1.0d+07
  double precision                                                                , parameter    :: coefficientSecondOrder                      = 0.0d+00
  double precision                                                                , parameter    :: correlationVelocityRadius                   =-0.3d+00
  double precision                                                                , parameter    :: toleranceRelativePotentialDifference        = 5.0d-03
  double precision                                                                , parameter    :: toleranceRelativeVelocityDispersion         = 1.0d-06
  double precision                                                                , parameter    :: toleranceRelativeVelocityDispersionMaximum  = 1.0d-03
  double precision                                                                , parameter    :: toleranceRelativePotential                  = 1.0d-03
  double precision                                                                , parameter    :: fractionRadiusFinalSmall                    = 1.0d-03
  logical                                                                         , parameter    :: velocityDispersionApproximate               =.false.
  logical                                                                         , parameter    :: tolerateVelocityMaximumFailure              =.false.
  logical                                                                         , parameter    :: tolerateEnclosedMassIntegrationFailure      =.false.
  logical                                                                         , parameter    :: tolerateVelocityDispersionFailure           =.false.
  logical                                                                         , parameter    :: toleratePotentialIntegrationFailure         =.false.
  class           (nodeComponentBasic                                            ), pointer      :: basic
  class           (nodeComponentSatellite                                        ), pointer      :: satellite
  double precision                                                                , dimension(3) :: radiusVirialFractional                      =[0.1d0,0.5d0,1.0d0]
  type            (cosmologyParametersSimple                                     )               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                )               :: cosmologyFunctions_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            )               :: darkMatterHaloScale_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)               :: virialDensityContrast_
  type            (treeNode                                                      )               :: node
  type            (varying_string                                                )               :: parameterFile
  type            (inputParameters                                               )               :: parameters
  type            (darkMatterProfileDMOIsothermal                                )               :: darkMatterProfileDMOIsothermal_
  type            (darkMatterProfileDMOHeated                                    ), target       :: darkMatterProfileDMOHeated_
  type            (darkMatterProfileDMOHeatedMonotonic                           ), target       :: darkMatterProfileDMOHeatedMonotonic_
  class           (darkMatterProfileDMOClass                                     ), pointer      :: darkMatterProfileDMO_
  type            (darkMatterProfileHeatingTidal                                 )               :: darkMatterProfileHeatingTidal_
  type            (darkMatterParticleCDM                                         )               :: darkMatterParticleCDM_
  type            (darkMatterParticleSelfInteractingDarkMatter                   )               :: darkMatterParticleSelfInteractingDarkMatter_
  class           (massDistributionSpherical                                     ), pointer      :: massDistributionSphericalHeated_
  class           (kinematicsDistributionClass                                   ), pointer      :: kinematicsDistributionHeated_
  class           (massDistributionHeatingClass                                  ), pointer      :: massDistributionHeatingTidal_
  class           (massDistributionClass                                         ), pointer      :: massDistributionIsothermal_
  type            (coordinateSpherical                                           )               :: coordinates                                                     , coordinatesHeated         , &
       &                                                                                            coordinatesInitial
  double precision                                                                               :: radiusVirial                                                    , radiusHeated              , &
       &                                                                                            massEnclosedAnalytic                                            , densityAnalytic           , &
       &                                                                                            radiusInitial                                                   , radiusInitialAnalytic     , &
       &                                                                                            massEnclosed                                                    , velocityDispersion        , &
       &                                                                                            potential                                                       , potentialAnalytic         , &
       &                                                                                            potentialZeroPoint                                              , potentialAnalyticZeroPoint, &
       &                                                                                            radius                                                          , toleranceRelative         , &
       &                                                                                            velocityDispersionAnalytic                                      , radiusScaled              , &
       &                                                                                            toleranceRelativeVelocityDispersionAssert                       , density
  integer                                                                                        :: i                                                               , profileType
  character       (len=5                                                         )               :: radiusLabel
  character       (len=128                                                       )               :: profileName

  ! Initialize.
  call displayVerbositySet(verbosityLevelStandard)
  parameterFile='testSuite/parameters/darkMatterProfileHeated.xml'
  parameters=inputParameters(parameterFile)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  call nodeClassHierarchyInitialize(parameters)
  ! Create the dark matter profiles.
  !![
  <referenceConstruct object="cosmologyParameters_"  >
   <constructor>
    cosmologyParametersSimple                                     (                                               &amp;
     &amp;                                                         OmegaMatter           = 0.2815d0             , &amp;
     &amp;                                                         OmegaBaryon           = 0.0465d0             , &amp;
     &amp;                                                         OmegaDarkEnergy       = 0.7185d0             , &amp;
     &amp;                                                         temperatureCMB        = 2.7800d0             , &amp;
     &amp;                                                         HubbleConstant        =69.3000d0               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"   >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                               &amp;
     &amp;                                                         cosmologyParameters_  =cosmologyParameters_    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_">
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                               &amp;
     &amp;                                                         tableStore            =.true.                , &amp;
     &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_     &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"                        >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                               &amp;
     &amp;                                                         cosmologyParameters_  =cosmologyParameters_  , &amp;
     &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_   , &amp;
     &amp;                                                         virialDensityContrast_=virialDensityContrast_  &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>  
  <referenceConstruct object="darkMatterParticleCDM_"                      >
   <constructor>
    darkMatterParticleCDM                             (                                                                                       &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticleSelfInteractingDarkMatter_">
   <constructor>
    darkMatterParticleSelfInteractingDarkMatter       (                                                                                       &amp;
     &amp;                                             crossSectionSelfInteraction         =1.0d0                                           , &amp;
     &amp;                                             darkMatterParticle_                 =darkMatterParticleCDM_                            &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  !!]
  darkMatterProfileHeatingTidal_      =darkMatterProfileHeatingTidal      (coefficientSecondOrder       ,coefficientSecondOrder       ,coefficientSecondOrder                                                                                   ,correlationVelocityRadius                                                          )
  darkMatterProfileDMOIsothermal_     =darkMatterProfileDMOIsothermal     (                                                                                                                                                                                                      darkMatterHaloScale_                               )
  darkMatterProfileDMOHeated_         =darkMatterProfileDMOHeated         (nonAnalyticSolversNumericalDMO,velocityDispersionApproximate,tolerateEnclosedMassIntegrationFailure,tolerateVelocityDispersionFailure,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,fractionRadiusFinalSmall,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,toleranceRelativePotential                           ,darkMatterProfileDMOIsothermal_,darkMatterProfileHeatingTidal_)
  darkMatterProfileDMOHeatedMonotonic_=darkMatterProfileDMOHeatedMonotonic(nonAnalyticSolversNumericalDMO                              ,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,darkMatterProfileDMOIsothermal_,darkMatterHaloScale_,darkMatterProfileHeatingTidal_)
  ! Set up the node.
  basic     => node%basic    (autoCreate=.true.)
  satellite => node%satellite(autoCreate=.true.)
  call basic    %massSet                  (massVirial     )
  call basic    %timeSet                  (time           )
  call basic    %timeLastIsolatedSet      (time           )
  call satellite%tidalHeatingNormalizedSet(heatingSpecific)
  ! Get the associated mass distribution.
  massDistributionIsothermal_   => darkMatterProfileDMOIsothermal_%get(node                                                                                                          )
  massDistributionHeatingTidal_ => darkMatterProfileHeatingTidal_ %get(node)
  ! Compute the characteristic radius for heating.
  radiusVirial=darkMatterHaloScale_%radiusVirial(node)
  radiusHeated=sqrt(gravitationalConstant_internal*massVirial/2.0d0/heatingSpecific/radiusVirial)
  ! Iterate over heated profile classes.
  do profileType=1,2
     select case (profileType)
     case (1)
        darkMatterProfileDMO_                     => darkMatterProfileDMOHeated_
        profileName                               =  "no shell crossing"
        toleranceRelative                         =  1.0d-6
        toleranceRelativeVelocityDispersionAssert =  1.0d-2
        select type (massDistributionIsothermal_)
        class is (massDistributionSpherical)
           allocate(massDistributionSphericalHeated :: massDistributionSphericalHeated_)
           allocate(kinematicsDistributionHeated    :: kinematicsDistributionHeated_   )
           select type (massDistributionSphericalHeated_)
           type is (massDistributionSphericalHeated         )
              !![
	      <referenceConstruct object="massDistributionSphericalHeated_" constructor="massDistributionSphericalHeated(nonAnalyticSolversNumerical,tolerateVelocityMaximumFailure,tolerateEnclosedMassIntegrationFailure,toleratePotentialIntegrationFailure,fractionRadiusFinalSmall,toleranceRelativePotential,massDistributionIsothermal_,massDistributionHeatingTidal_)"/>
	      !!]
           end select
           select type (kinematicsDistributionHeated_)
           type is (kinematicsDistributionHeated            )
              !![
              <referenceConstruct object="kinematicsDistributionHeated_"    constructor="kinematicsDistributionHeated   (nonAnalyticSolversNumerical,velocityDispersionApproximate,toleranceRelativeVelocityDispersion=1.0d-3,toleranceRelativeVelocityDispersionMaximum=1.0d-3)"/>
	      !!]
           end select
        end select
     case (2)
        darkMatterProfileDMO_                     => darkMatterProfileDMOHeatedMonotonic_
        profileName                               =  "monotonic perturbation"
        toleranceRelative                         =  7.0d-2
        toleranceRelativeVelocityDispersionAssert =  1.0d-1
        select type (massDistributionIsothermal_)
        class is (massDistributionSpherical)
           allocate(massDistributionSphericalHeatedMonotonic :: massDistributionSphericalHeated_)
           allocate(kinematicsDistributionCollisionless      :: kinematicsDistributionHeated_   )
           select type (massDistributionSphericalHeated_)
           type is (massDistributionSphericalHeatedMonotonic)
              !![
              <referenceConstruct object="massDistributionSphericalHeated_" constructor="massDistributionSphericalHeatedMonotonic(radiusVirial,nonAnalyticSolversNumerical,massDistributionIsothermal_,massDistributionHeatingTidal_)"/>
	      !!]
           end select
           select type (kinematicsDistributionHeated_)
           type is (kinematicsDistributionCollisionless)
              !![
              <referenceConstruct object="kinematicsDistributionHeated_"   constructor="kinematicsDistributionCollisionless     (                                                                                                  )"/>
	      !!]
           end select
        end select
     end select
     call Unit_Tests_Begin_Group("Heated dark matter profiles ("//trim(profileName)//")")
     ! Compute initial radius, enclosed mass, and density for a variety of radii and compare to the analytic solutions.
     potentialZeroPoint        =0.0d0
     potentialAnalyticZeroPoint=0.0d0
     do i=1,size(radiusVirialFractional)
        write (radiusLabel,'(f4.2)') radiusVirialFractional(i)
        radius               =+radiusVirial              &
             &                *radiusVirialFractional(i)
        coordinates          = [radius,0.0d0,0.0d0]
        select type (massDistributionSphericalHeated_)
        class is (massDistributionSphericalHeated)
           radiusInitial     =+massDistributionSphericalHeated_%radiusInitial           (radius                                                                       )
        class default
           radiusInitial     =-huge(0.0d0)
        end select
        massEnclosed         =+massDistributionSphericalHeated_%massEnclosedBySphere    (radius                                                                       )
        density              =+massDistributionSphericalHeated_%density                 (coordinates                                                                  )
        potential            =+massDistributionSphericalHeated_%potential               (coordinates                                                                  )
        velocityDispersion   =+kinematicsDistributionHeated_   %velocityDispersion1D    (coordinates,massDistributionSphericalHeated_,massDistributionSphericalHeated_)
        radiusInitialAnalytic=+sqrt(                            &
             &                      +  radiusHeated         **4 &
             &                      /4.0d0                      &
             &                      /  radius               **2 &
             &                      +  radiusHeated         **2 &
             &                     )                            &
             &                -        radiusHeated         **2 &
             &                /2.0d0                            &
             &                /        radius
        radiusScaled         =+radiusInitialAnalytic            &
             &                /radiusHeated
        coordinatesHeated    = [radiusHeated         ,0.0d0,0.0d0]
        coordinatesInitial   = [radiusInitialAnalytic,0.0d0,0.0d0]
        massEnclosedAnalytic =+massVirial                       &
             &                *        radiusInitialAnalytic    &
             &                /        radiusVirial
        densityAnalytic      =+massVirial                       &
             &                /4.0d0                            &
             &                /Pi                               &
             &                /        radiusVirial             &
             &                /        radius               **2 &
             &                *(                                &
             &                  -      radiusHeated         **4 &
             &                  /4.0d0                          &
             &                  /      radius               **3 &
             &                  /sqrt(                          &
             &                        +radiusHeated         **4 &
             &                        /4.0d0                    &
             &                        /radius               **2 &
             &                        +radiusHeated         **2 &
             &                       )                          &
             &                  +      radiusHeated         **2 &
             &                  /2.0d0                          &
             &                  /      radius               **2 &
             &                 )
        ! Evaluate the expression for the potential. Note that the atanh() term gives a complex result. However, the imaginary
        ! component of this is independent of radius and so cancels when finding potential differences. As such, we simply discard
        ! it here.
        potentialAnalytic=real(                                                                                                                  &
             &                 +gravitationalConstant_internal                                                                                   &
             &                 *massVirial                                                                                                       &
             &                 /4.0d0                                                                                                            &
             &                 /radiusVirial                                                                                                     &
             &                 *(                                                                                                                &
             &                   + radiusHeated**2                                                                                               &
             &                   -                               sqrt(4.0d0*radiusHeated**2+     radiusHeated**4/radius**2)*radius            )  &
             &                   *                               sqrt(      radiusHeated**2+4.0d0               *radius**2                    )  &
             &                   *(radiusHeated    +             sqrt(      radiusHeated**2+4.0d0               *radius**2                    )) &
             &                 *(                                                                                                                &
             &                   + radiusHeated                                                                                                  &
             &                   +(radiusHeated    +             sqrt(       radiusHeated**2+4.0d0              *radius**2                    )) &
             &                   *                  atanh(dcmplx(sqrt(       radiusHeated**2+4.0d0              *radius**2)/radiusHeated,0.0d0)) &
             &                  )                                                                                                                &
             &                 /   radiusHeated                                                                                                  &
             &                 /                                                                                 radius**2                       &
             &                 /(                                                                                                                &
             &                   +                                                          4.0d0               *radius**2                       &
             &                   + radiusHeated                                                                                                  &
             &                   *(radiusHeated                 +sqrt(      radiusHeated**2+4.0d0               *radius**2                    )) &
             &                  )                                                                                                                &
             &                )
        if (i == 1) then
           ! Compute potential zero points at our first radius - we compare only potential differences.
           potentialZeroPoint        =potential
           potentialAnalyticZeroPoint=potentialAnalytic
        end if
        potential        =+potential                  &
             &            -potentialZeroPoint
        potentialAnalytic=+potentialAnalytic          &
             &            -potentialAnalyticZeroPoint
        velocityDispersionAnalytic=+sqrt(                                                                                                                                                         &
             &                           +gravitationalConstant_internal                                                                                                                          &
             &                           *(                                                                                                                                                       &
             &                             +massDistributionIsothermal_%massEnclosedBySphere(radiusHeated         )*massDistributionIsothermal_%density(coordinatesHeated )/radiusHeated          &
             &                             *(+10.0d0                                                                                                      )                                       &
             &                             -massDistributionIsothermal_%massEnclosedBySphere(radiusInitialAnalytic)*massDistributionIsothermal_%density(coordinatesInitial)/radiusInitialAnalytic &
             &                             *(- 3.0d0+18.0d0*radiusScaled**4-6.0d0*radiusScaled**6+radiusScaled**8-24.0d0*radiusScaled**2*log(radiusScaled))                                       &
             &                            )                                                                                                                                                       &
             &                           /6.0d0                                                                                                                                                   &
             &                           /densityAnalytic                                                                                                                                         &
             &                          )
    
        if (profileType == 1) &
             & call Assert('r='//trim(radiusLabel)//'rᵥ; initial radius'     ,radiusInitial     ,radiusInitialAnalytic     ,relTol=toleranceRelative                        )
        call        Assert('r='//trim(radiusLabel)//'rᵥ; enclosed mass'      ,massEnclosed      ,massEnclosedAnalytic      ,relTol=toleranceRelative                        )
        call        Assert('r='//trim(radiusLabel)//'rᵥ; density'            ,density           ,densityAnalytic           ,relTol=toleranceRelative                        )
        call        Assert('r='//trim(radiusLabel)//'rᵥ; potential'          ,potential         ,potentialAnalytic         ,relTol=toleranceRelativePotentialDifference     )
        call        Assert('r='//trim(radiusLabel)//'rᵥ; velocity dispersion',velocityDispersion,velocityDispersionAnalytic,relTol=toleranceRelativeVelocityDispersionAssert)
     end do
     call Unit_Tests_End_Group()
     !![
     <objectDestructor name="massDistributionSphericalHeated_"/>
     <objectDestructor name="kinematicsDistributionHeated_"   />
     !!]
  end do
  call nodeClassHierarchyFinalize()
  call Unit_Tests_Finish         ()
end program Test_Dark_Matter_Profiles_Heated
