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
Implements a galactic (high- or low-pass) filter for total star formation rate with a non-parametric dependence on stellar mass.
!!}

  use :: Numerical_Interpolation                   , only : interpolator
  use :: Star_Formation_Rates_Disks                , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids            , only : starFormationRateSpheroidsClass
  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass

  !![
  <enumeration>
   <name>filterType</name>
   <description>Used to specify the type of filter for non-parametric star formation rate filtering.</description>
   <encodeFunction>yes</encodeFunction>
   <visibility>public</visibility>
   <entry label="lowPass"  description="Low-pass filter" />
   <entry label="highPass" description="High-pass filter"/>
  </enumeration>
  !!]

  !![
  <galacticFilter name="galacticFilterStarFormationRateNonParametric">
   <description>
    A galactic (high- or low-pass) filter for star formation rate. Galaxies with a combined disk, spheroid, plus \gls{nsc} star
    formation rate are passed if they are above or below (for {\normalfont \ttfamily [filterType])$=${\normalfont \ttfamily
    highPass} or {\normalfont \ttfamily lowPass} respectively) a mass-dependent threshold. The threshold is linearly interpolated
    in log({\normalfont \ttfamily [rateStarFormation]}) vs. log({\normalfont \ttfamily [massStellar]}).
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterStarFormationRateNonParametric
     !!{
     A galactic high-pass filter class for star formation rate.
     !!}
     private
     class           (starFormationRateDisksClass              ), pointer                   :: starFormationRateDisks_               => null()
     class           (starFormationRateSpheroidsClass          ), pointer                   :: starFormationRateSpheroids_           => null()
     class           (starFormationRateNuclearStarClustersClass), pointer                   :: starFormationRateNuclearStarClusters_ => null()
     double precision                                           , allocatable, dimension(:) :: massStellar                                    , rateStarFormation
     type            (enumerationFilterTypeType                )                            :: filterType
     type            (interpolator                             )                            :: interpolator_
   contains
     final     ::           starFormationRateNonParametricDestructor
     procedure :: passes => starFormationRateNonParametricPasses
  end type galacticFilterStarFormationRateNonParametric

  interface galacticFilterStarFormationRateNonParametric
     !!{
     Constructors for the \refClass{galacticFilterStarFormationRateNonParametric} galactic filter class.
     !!}
     module procedure starFormationRateNonParametricConstructorParameters
     module procedure starFormationRateNonParametricConstructorInternal
  end interface galacticFilterStarFormationRateNonParametric

contains

  function starFormationRateNonParametricConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterStarFormationRateNonParametric} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterStarFormationRateNonParametric)                              :: self
    type            (inputParameters                             ), intent(inout)               :: parameters
    class           (starFormationRateDisksClass                 ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass             ), pointer                     :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass   ), pointer                     :: starFormationRateNuclearStarClusters_
    double precision                                              , allocatable  , dimension(:) :: massStellar                          , rateStarFormation
    type            (varying_string                              )                              :: filterType

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>filterType</name>
      <source>parameters</source>
      <description>The type of filter: {\normalfont \ttfamily lowPass} or {\normalfont \ttfamily highPass}.</description>
    </inputParameter>
    <inputParameter>
      <name>massStellar</name>
      <source>parameters</source>
      <description>The list of stellar masses at which the star formation rate threshold is specified.</description>
    </inputParameter>
    <inputParameter>
      <name>rateStarFormation</name>
      <source>parameters</source>
      <description>The list of star formation rate thresholds for each stellar mass.</description>
    </inputParameter>
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=galacticFilterStarFormationRateNonParametric(enumerationFilterTypeEncode(char(filterType),includesPrefix=.false.),massStellar,rateStarFormation,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function starFormationRateNonParametricConstructorParameters

  function starFormationRateNonParametricConstructorInternal(filterType,massStellar,rateStarFormation,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterStarFormationRateNonParametric} galactic filter class.
    !!}
    use :: Table_Labels           , only : extrapolationTypeExtrapolate
    use :: Numerical_Interpolation, only : gsl_interp_linear
    implicit none
    type            (galacticFilterStarFormationRateNonParametric)                              :: self
    type            (enumerationFilterTypeType                   ), intent(in   )               :: filterType
    double precision                                              , intent(in   ), dimension(:) :: massStellar                          , rateStarFormation
    class           (starFormationRateDisksClass                 ), intent(in   ), target       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass             ), intent(in   ), target       :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass   ), intent(in   ), target       :: starFormationRateNuclearStarClusters_

    !![
    <constructorAssign variables="filterType, massStellar, rateStarFormation, *starFormationRateDisks_, *starFormationRateSpheroids_, *starFormationRateNuclearStarClusters_"/>
    !!]

    self%interpolator_=interpolator(log(massStellar),log(rateStarFormation),interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeExtrapolate)
    return
  end function starFormationRateNonParametricConstructorInternal

  subroutine starFormationRateNonParametricDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterStarFormationRateNonParametric} galactic filter class.
    !!}
    implicit none
    type(galacticFilterStarFormationRateNonParametric), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateDisks_"              />
    <objectDestructor name="self%starFormationRateSpheroids_"          />
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine starFormationRateNonParametricDestructor

  logical function starFormationRateNonParametricPasses(self,node) result(passes)
    !!{
    Implement an star formation rate galactic filter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC, treeNode
    implicit none
    class           (galacticFilterStarFormationRateNonParametric), intent(inout)         :: self
    type            (treeNode                                    ), intent(inout), target :: node
    class           (nodeComponentDisk                           ), pointer               :: disk
    class           (nodeComponentSpheroid                       ), pointer               :: spheroid
    class           (nodeComponentNSC                            ), pointer               :: nuclearStarCluster
    double precision                                                                      :: starFormationRate         , massStellar, &
         &                                                                                   starFormationRateThreshold

    disk               => node              %disk                                      (    )
    spheroid           => node              %spheroid                                  (    ) 
    nuclearStarCluster => node              %NSC                                       (    )
    massStellar        = +disk              %massStellar                               (    ) &
         &               +spheroid          %massStellar                               (    ) &
         &               +nuclearStarCluster%massStellar                               (    )
    starFormationRate  = +self              %starFormationRateDisks_              %rate(node) &
         &               +self              %starFormationRateSpheroids_          %rate(node) &
         &               +self              %starFormationRateNuclearStarClusters_%rate(node)
    if (massStellar > 0.0d0) then
       starFormationRateThreshold=exp(self%interpolator_%interpolate(log(massStellar)))
       passes                    = (starFormationRate >= starFormationRateThreshold) &
            &                     .eqv.&
            &                      (self%filterType   == filterTypeHighPass        )
    else
       passes                    =                                                     self%filterType == filterTypeLowPass
    end if
    return
  end function starFormationRateNonParametricPasses
