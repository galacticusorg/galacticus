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

!!{
Implements a galactic high-pass filter for total star formation rate.
!!}

  use :: Star_Formation_Rates_Disks                , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids            , only : starFormationRateSpheroidsClass
  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass

  !![
  <galacticFilter name="galacticFilterStarFormationRate">
   <description>
   A galactic high-pass filter for star formation rate. Galaxies with a combined disk,
   spheroid, plus \gls{nsc} star formation rate greater than or equal to a mass-dependent threshold. The threshold is given by
   \begin{equation}
   \log_{10} \left( { \dot{\phi}_\mathrm{t} \over M_\odot\,\hbox{Gyr}^{-1}} \right) = \alpha_0 + \alpha_1  \left( \log_{10} M_\star - \log_{10} M_0 \right),
   \end{equation}
   where $M_0=${\normalfont \ttfamily [starFormationRateThresholdLogM0]}, $\alpha_0=${\normalfont
   \ttfamily [starFormationRateThresholdLogSFR0]}, and $\alpha_1=${\normalfont \ttfamily
   [starFormationRateThresholdLogSFR1]}.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterStarFormationRate
     !!{
     A galactic high-pass filter class for star formation rate.
     !!}
     private
     class           (starFormationRateDisksClass              ), pointer :: starFormationRateDisks_               => null()
     class           (starFormationRateSpheroidsClass          ), pointer :: starFormationRateSpheroids_           => null()
     class           (starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     double precision                                                     :: logM0                                          , logSFR0, &
          &                                                                  logSFR1
   contains
     final     ::           starFormationRateDestructor
     procedure :: passes => starFormationRatePasses
  end type galacticFilterStarFormationRate

  interface galacticFilterStarFormationRate
     !!{
     Constructors for the \refClass{galacticFilterStarFormationRate} galactic filter class.
     !!}
     module procedure starFormationRateConstructorParameters
     module procedure starFormationRateConstructorInternal
  end interface galacticFilterStarFormationRate

contains

  function starFormationRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterStarFormationRate} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterStarFormationRate          )                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (starFormationRateDisksClass              ), pointer       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass          ), pointer       :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass), pointer       :: starFormationRateNuclearStarClusters_
    double precision                                                           :: logM0                                , logSFR0, &
         &                                                                        logSFR1

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>logM0</name>
      <source>parameters</source>
      <defaultValue>10.0d0</defaultValue>
      <description>The parameter $\log_{10} M_0$ (with $M_0$ in units of $M_\odot$) appearing in the star formation rate threshold expression for the star formation rate galactic filter class.</description>
    </inputParameter>
    <inputParameter>
      <name>logSFR0</name>
      <source>parameters</source>
      <defaultValue>9.0d0</defaultValue>
      <description>The parameter $\alpha_0$ appearing in the star formation rate threshold expression for the star formation rate galactic filter class.</description>
    </inputParameter>
    <inputParameter>
      <name>logSFR1</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $\alpha_1$ appearing in the star formation rate threshold expression for the star formation rate galactic filter class.</description>
    </inputParameter>
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=galacticFilterStarFormationRate(logM0,logSFR0,logSFR1,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function starFormationRateConstructorParameters

  function starFormationRateConstructorInternal(logM0,logSFR0,logSFR1,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterStarFormationRate} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterStarFormationRate          )                        :: self
    double precision                                           , intent(in   )         :: logM0                                , logSFR0, &
         &                                                                                logSFR1
    class           (starFormationRateDisksClass              ), intent(in   ), target :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass          ), intent(in   ), target :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass), intent(in   ), target :: starFormationRateNuclearStarClusters_

    !![
    <constructorAssign variables="logM0, logSFR0, logSFR1, *starFormationRateDisks_, *starFormationRateSpheroids_, *starFormationRateNuclearStarClusters_"/>
    !!]

    return
  end function starFormationRateConstructorInternal

  subroutine starFormationRateDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterStarFormationRate} galactic filter class.
    !!}
    implicit none
    type(galacticFilterStarFormationRate), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateDisks_"              />
    <objectDestructor name="self%starFormationRateSpheroids_"          />
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine starFormationRateDestructor

  logical function starFormationRatePasses(self,node)
    !!{
    Implement an star formation rate galactic filter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC, treeNode
    implicit none
    class           (galacticFilterStarFormationRate), intent(inout)         :: self
    type            (treeNode                       ), intent(inout), target :: node
    class           (nodeComponentDisk              ), pointer               :: disk
    class           (nodeComponentSpheroid          ), pointer               :: spheroid
    class           (nodeComponentNSC               ), pointer               :: nuclearStarCluster
    double precision                                                         :: starFormationRate         , stellarMass, &
         &                                                                      starFormationRateThreshold

    disk               => node              %disk                                      (    )
    spheroid           => node              %spheroid                                  (    ) 
    nuclearStarCluster => node              %NSC                                       (    )
    stellarMass        = +disk              %massStellar                               (    ) &
         &               +spheroid          %massStellar                               (    ) &
         &               +nuclearStarCluster%massStellar                               (    )
    starFormationRate  = +self              %starFormationRateDisks_              %rate(node) &
         &               +self              %starFormationRateSpheroids_          %rate(node) &
         &               +self              %starFormationRateNuclearStarClusters_%rate(node)
    if (stellarMass > 0.0d0) then
       starFormationRateThreshold= +10.0d0**(                      &
            &                                +  self%logSFR0       &
            &                                +  self%logSFR1       &
            &                                *(                    &
            &                                  +log10(stellarMass) &
            &                                  -self%logM0         &
            &                                )                     &
            &                               )
       starFormationRatePasses=(starFormationRate >= starFormationRateThreshold)
    else
       starFormationRatePasses=.false.
    end if
    return
  end function starFormationRatePasses
