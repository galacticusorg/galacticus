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
  Implements calculations of satellite merging times using the \cite{jiang_fitting_2008} method.
  !!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <satelliteMergingTimescales name="satelliteMergingTimescalesJiang2008">
   <description>
    A satellite merging timescales class which computes merging timescales using the dynamical friction calibration of
    \cite{jiang_fitting_2008}. \cite{jiang_fitting_2008} find that their fitting formula does not perfectly capture the results
    of N-body simulations, instead showing a scatter in the ratio $T_\mathrm{fit}/T_\mathrm{sim}$ where $T_\mathrm{fit}$ is the
    merging time from their fitting formula and $T_\mathrm{sim}$ is that measured from their N-body simulation. Furthermore,
    they show that the distribution of $T_\mathrm{fit}/T_\mathrm{sim}$ is well described by a log-normal with dispersion (in
    $\log[T_\mathrm{fit}/T_\mathrm{sim}]$) of $\sigma=0.4$. Random perturbations can be applied to the merger times returned by
    this implementation by setting $\sigma=${\normalfont \ttfamily [scatter]}$>0$ which will cause the merger time to be drawn
    from a log-normal distribution of width $\sigma$ with median equal to $T_\mathrm{fit}$.
   </description>
  </satelliteMergingTimescales>
  !!]
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesJiang2008
     !!{
     A class implementing the \cite{jiang_fitting_2008} method for satellite merging timescales.
     !!}
     private
     class          (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     class          (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     double precision                                    :: timescaleMultiplier
     ! Scatter (in log(T_merge)) to add to the merger times.
     double precision                                    :: scatter
   contains
     final     ::                     jiang2008Destructor
     procedure :: timeUntilMerging => jiang2008TimeUntilMerging
  end type satelliteMergingTimescalesJiang2008

  interface satelliteMergingTimescalesJiang2008
     !!{
     Constructors for the \cite{jiang_fitting_2008} merging timescale class.
     !!}
     module procedure jiang2008ConstructorParameters
     module procedure jiang2008ConstructorInternal
  end interface satelliteMergingTimescalesJiang2008

contains

  function jiang2008ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{jiang_fitting_2008} merging timescale class which builds the object from a parameter set.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : defaultBasicComponent
    use :: Input_Parameters, only : inputParameter       , inputParameters
    implicit none
    type            (satelliteMergingTimescalesJiang2008)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass          ), pointer       :: darkMatterProfileDMO_
    double precision                                                     :: scatter              , timescaleMultiplier

    if (.not.defaultBasicComponent%massIsGettable()) call Error_Report('this method requires that the "mass" property of the basic component be gettable'//{introspection:location})
    !![
    <inputParameter>
      <name>timescaleMultiplier</name>
      <defaultValue>0.75d0</defaultValue>
      <description>A multiplier for the merging timescale in dynamical friction timescale calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scatter</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Specifies whether or not to add random scatter to the dynamical friction timescales in the {\normalfont \ttfamily Jiang2008} satellite merging time implementation.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=satelliteMergingTimescalesJiang2008(timescaleMultiplier,scatter,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function jiang2008ConstructorParameters

  function jiang2008ConstructorInternal(timescaleMultiplier,scatter,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Constructor for the \cite{jiang_fitting_2008} merging timescale class.
    !!}
    implicit none
    type            (satelliteMergingTimescalesJiang2008)                        :: self
    double precision                                     , intent(in   )         :: timescaleMultiplier  , scatter
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass          ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="timescaleMultiplier, scatter, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function jiang2008ConstructorInternal

  subroutine jiang2008Destructor(self)
    !!{
    Destructor for the \refClass{satelliteMergingTimescalesJiang2008} satellite merging timescale class.
    !!}
    implicit none
    type(satelliteMergingTimescalesJiang2008), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine jiang2008Destructor

  double precision function jiang2008TimeUntilMerging(self,node,orbit)
    !!{
    Return the timescale for merging satellites using the \cite{jiang_fitting_2008} method.
    !!}
    use :: Error             , only : Error_Report
    use :: Galacticus_Nodes  , only : nodeComponentBasic                              , treeNode
    use :: Mass_Distributions, only : massDistributionClass
    use :: Satellite_Orbits  , only : Satellite_Orbit_Equivalent_Circular_Orbit_Radius, errorCodeNoEquivalentOrbit, errorCodeOrbitUnbound, errorCodeSuccess
    implicit none
    class           (satelliteMergingTimescalesJiang2008), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    type            (keplerOrbit                        ), intent(inout) :: orbit
    type            (treeNode                           ), pointer       :: nodeHost
    class           (nodeComponentBasic                 ), pointer       :: basicHost                            , basic
    class           (massDistributionClass              ), pointer       :: massDistribution_
    logical                                              , parameter     :: acceptUnboundOrbits          =.false.

    double precision                                     , parameter     :: C                            =0.43d0 , a            =0.94d0, &  !   Fitting parameters from Jiang's paper.
         &                                                                  b                            =0.60d0 , d            =0.60d0
    integer                                                              :: errorCode
    double precision                                                     :: equivalentCircularOrbitRadius        , massRatio           , &
         &                                                                  orbitalCircularity                   , radialScale         , &
         &                                                                  velocityScale

    ! Find the host node.
    if (node%isSatellite()) then
       nodeHost => node%parent
    else
       nodeHost => node%parent%firstChild
    end if
    ! Get the equivalent circular orbit.
    equivalentCircularOrbitRadius=Satellite_Orbit_Equivalent_Circular_Orbit_Radius(nodeHost,orbit,self%darkMatterHaloScale_,errorCode)
    ! Check error codes.
    select case (errorCode)
    case (errorCodeOrbitUnbound     )
       jiang2008TimeUntilMerging=satelliteMergeTimeInfinite
       return
    case (errorCodeNoEquivalentOrbit)
       ! Circularity is not defined. Assume instantaneous merging.
       jiang2008TimeUntilMerging=0.0d0
       return
    case (errorCodeSuccess          )
    case default
       call Error_Report('unrecognized error code'//{introspection:location})
    end select
    ! Get velocity scale.
    velocityScale=self%darkMatterHaloScale_%velocityVirial(nodeHost)
    radialScale  =self%darkMatterHaloScale_%radiusVirial  (nodeHost)
    ! Compute orbital circularity.
    massDistribution_  =>  self             %darkMatterProfileDMO_%get            (nodeHost                     )
    orbitalCircularity =  +orbit                                  %angularMomentum(                             ) &
         &                /massDistribution_                      %rotationCurve  (equivalentCircularOrbitRadius) &
         &                /equivalentCircularOrbitRadius
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Compute mass ratio (mass in host [not including satellite if the node is already a satellite] divided by mass in satellite).
    basic     =>  node     %basic()
    basicHost =>  nodeHost %basic()
    if (node%isSatellite()) then
       massRatio=+basicHost%mass () &
            &    /basic    %mass () &
            &    -1.0d0
    else
       massRatio=+basicHost%mass () &
            &    /basic    %mass ()
    end if
    ! Check for a non-zero mass ratio.
    if (massRatio <= 0.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       jiang2008TimeUntilMerging=0.0d0
    else
       ! Compute dynamical friction timescale.
       jiang2008TimeUntilMerging=+self%timescaleMultiplier                               &
            &                    *self%darkMatterHaloScale_%timescaleDynamical(nodeHost) &
            &                    *sqrt(equivalentCircularOrbitRadius/radialScale)        &
            &                    *((a*(orbitalCircularity**b)+d)/2.0d0/C)                &
            &                    *          massRatio                                    &
            &                    /log(1.0d0+massRatio)
       ! Add scatter if necessary.
       if (self%scatter > 0.0d0)                                                                          &
            & jiang2008TimeUntilMerging=+jiang2008TimeUntilMerging                                        &
            &                           *exp(                                                             &
            &                                +self%scatter                                                &
            &                                *node%hostTree%randomNumberGenerator_%standardNormalSample() &
            &                               )
    end if
    return
  end function jiang2008TimeUntilMerging
