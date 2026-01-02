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

  !+ Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !!{
  Implementation of a satellite dynamical friction class which applies the core-stalling model of \cite{kaur_stalling_2018} to another dynamical
  friction class.
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Numerical_Interpolation , only : interpolator

  !![
  <satelliteDynamicalFriction name="satelliteDynamicalFrictionKaur2018">
   <description>
    A satellite dynamical friction class which applies the core-stalling model of \cite{kaur_stalling_2018} to another dynamical
    friction class. Specifically, the acceleration due to dynamical friction is given by:
    \begin{equation}
     \mathbf{a}_\mathrm{df} = \mathbf{a}^\prime_\mathrm{df} [S_\mathrm{trail}(r/r_*)+S_\mathrm{lead}(r/r_*)],
    \end{equation}
    where $\mathbf{a}^\prime_\mathrm{df}$ is the acceleration provided by the child dynamical friction class, and
    $S_\mathrm{trail}(x)$ and $S_\mathrm{lead}(x)$ are the torque suppression factors defined by
    \cite[][eqn.~79]{kaur_stalling_2018}. The functional forms for $S(x)$ are taken from Figure~8 of
    \cite{kaur_stalling_2018}---extrapolation in $\log(S)$ vs. $x$ is used to extend $S(x)$ to lower $x$ than shown in that figure,
    while for values of $x$ higher than that shown in that figure $S(x)$ is held constant at the maximum value shown.
  
    The stalling radius, $r_*$, is computed using equation~(55) of \cite{kaur_stalling_2018}. For dark matter halo profiles with a
    central cusp (for which that equation has no solution), the acceleration provided by the child dynamical friction class is
    returned unmodified.
   </description>
  </satelliteDynamicalFriction>
  !!]
  type, extends(satelliteDynamicalFrictionClass) :: satelliteDynamicalFrictionKaur2018
     !!{
     Implementation of a satellite dynamical friction class which applies the core-stalling model of \cite{kaur_stalling_2018}
     to another dynamical friction class.
     !!}
     private
     class           (satelliteDynamicalFrictionClass), pointer :: satelliteDynamicalFriction_         => null()
     class           (darkMatterProfileDMOClass      ), pointer :: darkMatterProfileDMO_               => null()
     type            (interpolator                   )          :: factorSuppressionLeadingLogarithmic          , factorSuppressionTrailingLogarithmic
     double precision                                           :: radiusDimensionlessMaximum
   contains
     final     ::                 kaur2018Destructor
     procedure :: acceleration => kaur2018Acceleration
  end type satelliteDynamicalFrictionKaur2018

  interface satelliteDynamicalFrictionKaur2018
     !!{
     Constructors for the kaur2018 satellite dynamical friction class.
     !!}
     module procedure kaur2018ConstructorParameters
     module procedure kaur2018ConstructorInternal
  end interface satelliteDynamicalFrictionKaur2018

  ! Submodule-scope variables used in root finding.
  class(satelliteDynamicalFrictionKaur2018), pointer :: self_
  !$omp threadprivate(self_)
  
contains

  function kaur2018ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteDynamicalFrictionKaur2018} satellite dynamical friction class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (satelliteDynamicalFrictionKaur2018)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_
    class(satelliteDynamicalFrictionClass   ), pointer       :: satelliteDynamicalFriction_

    !![
    <objectBuilder class="satelliteDynamicalFriction" name="satelliteDynamicalFriction_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"       name="darkMatterProfileDMO_"       source="parameters"/>
    !!]
    self=satelliteDynamicalFrictionKaur2018(satelliteDynamicalFriction_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteDynamicalFriction_"/>
    <objectDestructor name="darkMatterProfileDMO_"      />
    !!]
    return
  end function kaur2018ConstructorParameters

  function kaur2018ConstructorInternal(satelliteDynamicalFriction_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{satelliteDynamicalFrictionKaur2018} satellite dynamical friction class.
    !!}
    use :: Table_Labels, only : extrapolationTypeExtrapolate
    implicit none
    type (satelliteDynamicalFrictionKaur2018)                        :: self
    class(darkMatterProfileDMOClass         ), intent(in   ), target :: darkMatterProfileDMO_
    class(satelliteDynamicalFrictionClass   ), intent(in   ), target :: satelliteDynamicalFriction_
    !![
    <constructorAssign variables="*satelliteDynamicalFriction_, *darkMatterProfileDMO_"/>
    !!]
    
    ! Build interpolators for the suppression factor as a function of radius. These are extracted directly from the arXiv source
    ! of Figure 8 of Kaur & Sridhar (2018). Specifically we tabulate log₁₀(S) vs. r/r*.
    self%factorSuppressionLeadingLogarithmic =interpolator(                                                                                                                                            &
         &                                                 x                =[                                                                                                                         &
         &                                                                    +1.010644d0,+1.019908d0,+1.029171d0,+1.038434d0,+1.047697d0,+1.056961d0,+1.066224d0,+1.075488d0,+1.084751d0,+1.094014d0, &
         &                                                                    +1.103278d0,+1.112541d0,+1.121805d0,+1.131068d0,+1.140331d0,+1.149595d0,+1.158858d0,+1.168121d0,+1.177384d0,+1.186648d0, &
         &                                                                    +1.195912d0,+1.205175d0,+1.214438d0,+1.223702d0,+1.232965d0,+1.242228d0,+1.251491d0,+1.260755d0,+1.270018d0,+1.279281d0, &
         &                                                                    +1.288545d0,+1.297808d0,+1.307072d0,+1.316335d0,+1.325598d0,+1.334862d0,+1.344125d0,+1.353389d0,+1.362651d0,+1.371915d0  &
         &                                                                   ]                                                                                                                       , &
         &                                                 y                =[                                                                                                                         &
         &                                                                    -8.343948d0,-7.094500d0,-5.412923d0,-4.763084d0,-4.335815d0,-4.016717d0,-3.763314d0,-3.553291d0,-3.373846d0,-3.217850d0, &
         &                                                                    -3.081369d0,-2.960790d0,-2.853081d0,-2.756044d0,-2.667989d0,-2.587549d0,-2.513683d0,-2.445558d0,-2.382780d0,-2.324817d0, &
         &                                                                    -2.271252d0,-2.221622d0,-2.175557d0,-2.132687d0,-2.092756d0,-2.055441d0,-2.020441d0,-1.987525d0,-1.956529d0,-1.927339d0, &
         &                                                                    -1.899909d0,-1.874191d0,-1.850140d0,-1.827571d0,-1.806344d0,-1.786321d0,-1.767362d0,-1.749330d0,-1.732131d0,-1.715672d0  &
         &                                                                   ]                                                                                                                       , &
         &                                                 extrapolationType=extrapolationTypeExtrapolate                                                                                              &
         &                                                )

    self%factorSuppressionTrailingLogarithmic=interpolator(                                                                                                                                            &
         &                                                 x                =[                                                                                                                         &
         &                                                                    +0.522254d0,+0.685957d0,+0.692888d0,+0.699814d0,+0.706744d0,+0.713673d0,+0.720601d0,+0.727531d0,+0.734460d0,+0.741388d0, &
         &                                                                    +0.748317d0,+0.755247d0,+0.762175d0,+0.769104d0,+0.776032d0,+0.782962d0,+0.789890d0,+0.796819d0,+0.803747d0,+0.810677d0, &
         &                                                                    +0.817606d0,+0.824536d0,+0.831464d0,+0.838393d0,+0.845321d0,+0.852251d0,+0.859179d0,+0.866108d0,+0.873036d0,+0.879966d0, &
         &                                                                    +0.886895d0,+0.893823d0,+0.900752d0,+0.907682d0,+0.914610d0,+0.921539d0,+0.928467d0,+0.935397d0,+0.942325d0,+0.949254d0, &
         &                                                                    +0.956184d0,+0.963112d0,+0.970041d0,+0.976969d0,+0.983899d0,+0.990828d0,+0.997756d0,+1.004684d0,+1.011615d0,+1.018543d0, &
         &                                                                    +1.025471d0,+1.032401d0,+1.039330d0,+1.046258d0,+1.053187d0,+1.060117d0,+1.067045d0,+1.073974d0,+1.080904d0,+1.087832d0, &
         &                                                                    +1.094760d0,+1.101689d0,+1.108619d0,+1.115547d0,+1.122476d0,+1.129404d0,+1.136333d0,+1.143263d0,+1.150191d0,+1.157121d0, &
         &                                                                    +1.164048d0,+1.170978d0,+1.177906d0,+1.184836d0,+1.191765d0,+1.198693d0,+1.205622d0,+1.212552d0,+1.219480d0,+1.226409d0, &
         &                                                                    +1.233337d0,+1.240267d0,+1.247195d0,+1.254124d0,+1.261052d0,+1.267982d0,+1.274911d0,+1.281841d0,+1.288769d0,+1.295698d0, &
         &                                                                    +1.302626d0,+1.309556d0,+1.316485d0,+1.323413d0,+1.330341d0,+1.337271d0,+1.344200d0,+1.351128d0,+1.358057d0,+1.364985d0, &
         &                                                                    +1.371915d0                                                                                                              &
         &                                                                   ]                                                                                                                       , & 
         &                                                 y                =[                                                                                                                         &
         &                                                                    -5.749819d0,-4.263502d0,-4.142231d0,-4.068304d0,-3.982802d0,-3.882663d0,-3.780927d0,-3.698930d0,-3.663635d0,-3.677507d0, &
         &                                                                    -3.675094d0,-3.595371d0,-3.494156d0,-3.432558d0,-3.412177d0,-3.405475d0,-3.387281d0,-3.333843d0,-3.240666d0,-3.141325d0, &
         &                                                                    -3.060787d0,-3.006291d0,-2.979936d0,-2.972836d0,-2.950579d0,-2.886082d0,-2.783877d0,-2.662766d0,-2.558182d0,-2.496220d0, &
         &                                                                    -2.474866d0,-2.468547d0,-2.457627d0,-2.439241d0,-2.415005d0,-2.385856d0,-2.352331d0,-2.315803d0,-2.278373d0,-2.241376d0, &
         &                                                                    -2.204866d0,-2.169171d0,-2.139362d0,-2.121636d0,-2.118841d0,-2.131359d0,-2.145143d0,-2.101984d0,-1.967504d0,-1.782990d0, &
         &                                                                    -1.614223d0,-1.477538d0,-1.371583d0,-1.292781d0,-1.237069d0,-1.195038d0,-1.158892d0,-1.123267d0,-1.084326d0,-1.039829d0, &
         &                                                                    -0.990767d0,-0.939222d0,-0.886826d0,-0.834655d0,-0.783527d0,-0.734048d0,-0.686721d0,-0.641895d0,-0.599690d0,-0.560124d0, &
         &                                                                    -0.523058d0,-0.488197d0,-0.455263d0,-0.424047d0,-0.394342d0,-0.366009d0,-0.339047d0,-0.313405d0,-0.289064d0,-0.265974d0, &
         &                                                                    -0.244099d0,-0.223387d0,-0.203734d0,-0.185106d0,-0.167450d0,-0.150696d0,-0.134828d0,-0.119793d0,-0.105644d0,-0.092311d0, &
         &                                                                    -0.079793d0,-0.068074d0,-0.057085d0,-0.046772d0,-0.037067d0,-0.027953d0,-0.019376d0,-0.011286d0,-0.003630d0,+0.003627d0, &
         &                                                                    +0.010537d0                                                                                                              &
         &                                                                   ]                                                                                                                       , &
         &                                                 extrapolationType=extrapolationTypeExtrapolate                                                                                              &
         &                                                )
    self%radiusDimensionlessMaximum            =1.371915d0
    return
  end function kaur2018ConstructorInternal

  subroutine kaur2018Destructor(self)
    !!{
    Destructor for the simple cooling radius class.
    !!}
    implicit none
    type(satelliteDynamicalFrictionKaur2018), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"      />
    <objectDestructor name="self%satelliteDynamicalFriction_"/>
    !!]
    return
  end subroutine kaur2018Destructor

  function kaur2018Acceleration(self,node)
    !!{
    Return an acceleration for satellites due to dynamical friction using the core-stalling model of \cite{kaur_stalling_2018}.
    !!}
    use :: Coordinates               , only : coordinateSpherical      , assignment(=)
    use :: Galacticus_Nodes          , only : nodeComponentSatellite
    use :: Galactic_Structure_Options, only : coordinateSystemCartesian
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Root_Finder               , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    double precision                                    , dimension(3)          :: kaur2018Acceleration
    class           (satelliteDynamicalFrictionKaur2018), intent(inout), target :: self
    type            (treeNode                          ), intent(inout)         :: node
    class           (nodeComponentSatellite            ), pointer               :: satellite
    type            (treeNode                          ), pointer               :: nodeHost
    class           (massDistributionClass             ), pointer               :: massDistribution_
    double precision                                    , dimension(3)          :: position
    double precision                                    , parameter             :: toleranceAbsolute   =0.0d+0, toleranceRelative                   =1.0d-3
    double precision                                                            :: massSatellite              , densityHostCentral                         , &
         &                                                                         factorSuppression          , logSlopeDensityProfileDarkMatterHost       , &
         &                                                                         radiusOrbital              , radiusStalling                             , &
         &                                                                         radiusDimensionless
    type            (rootFinder                        )                        :: finder
    type            (coordinateSpherical               )                        :: coordinates
    
    ! Compute the base acceleration.    
    kaur2018Acceleration=+self%satelliteDynamicalFriction_%acceleration(node)
    ! Get host node and satellite component.
    nodeHost      => node     %mergesWith()
    satellite     => node     %satellite ()
    massSatellite =  satellite%boundMass ()
    position      =  satellite%position  ()
    radiusOrbital =  sqrt(sum(position**2))
    ! Check for physically-plausible satellite.
    if (massSatellite <= 0.0d0) return
    ! Check if the density profile has a finite density at the center. We do this by considering the logarithmic slope of the
    ! density profile. For cusped density profiles, we assume no stalling.
    coordinates                          =  [0.0d0,0.0d0,0.0d0]
    massDistribution_                    => self             %darkMatterProfileDMO_%get                  (nodeHost                      )
    logSlopeDensityProfileDarkMatterHost =  massDistribution_                      %densityGradientRadial(coordinates,logarithmic=.true.)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    if (logSlopeDensityProfileDarkMatterHost < 0.0d0) return
    ! Find the stalling radius.
    self_              => self
    massDistribution_  => nodeHost         %massDistribution(           )
    densityHostCentral =  massDistribution_%density         (coordinates)
    finder             =  rootFinder(                                                             &
         &                           rootFunction                 =radiusStallingRoot           , &
         &                           toleranceAbsolute            =toleranceAbsolute            , &
         &                           toleranceRelative            =toleranceRelative            , &
         &                           rangeExpandDownward          =0.5d0                        , &
         &                           rangeExpandUpward            =2.0d0                        , &
         &                           rangeExpandType              =rangeExpandMultiplicative    , &
         &                           rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
         &                           rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &
         &                          )
    radiusStalling=finder%find(rootGuess=radiusOrbital)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Compute the suppression factor.
    radiusDimensionless=+radiusOrbital  &
         &              /radiusStalling
    if (radiusDimensionless >= self%radiusDimensionlessMaximum) then
       factorSuppression=+ 1.0d0
    else
       factorSuppression=+10.0d0**self%factorSuppressionLeadingLogarithmic %interpolate(radiusDimensionless) &
            &            +10.0d0**self%factorSuppressionTrailingLogarithmic%interpolate(radiusDimensionless)
    end if
    ! Evaluate the dynamical friction acceleration.
    kaur2018Acceleration=+kaur2018Acceleration &
         &               *factorSuppression
    return
    
  contains
    
    double precision function radiusStallingRoot(radiusStalling)
      !!{
      Root function used in finding the stalling radius.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: radiusStalling
      
      radiusStallingRoot=+4.0d0                                                  &
           &             *Pi                                                     &
           &             /3.0d0                                                  &
           &             *densityHostCentral                                     &
           &             *radiusStalling    **3                                  &
           &             -massDistribution_%massEnclosedBySphere(radiusStalling) &
           &             -massSatellite
      return
    end function radiusStallingRoot
    
  end function kaur2018Acceleration
