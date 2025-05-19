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
  Implements a merger tree operator which accumulates conditional mass functions for trees.
  !!}

  use    :: Cosmology_Functions              , only : cosmologyFunctionsClass
  !$ use :: OMP_Lib                          , only : omp_lock_kind
  use    :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass

  !![
  <mergerTreeOperator name="mergerTreeOperatorConditionalMF">
   <description>
    Provides a merger tree operator which accumulates conditional mass functions for trees. In
    addition to the cumulative mass function, 1$^\mathrm{st}$ through $n^\mathrm{th}$ most-massive
    progenitor mass functions, formation rate functions, and unevolved subhalo mass functions \citep{jiang_generating_2014}
    split by hierarchy depth are computed and output. Mass functions are accumulated in
    logarithmically-spaced bins of parent halo mass, logarithmically-spaced bins of mass ratio
    (the ratio of progenitor to parent halo mass), and at pairs of parent/progenitor
    redshifts. The following parameters control the operator:
    \begin{description}
    \item[{\normalfont \ttfamily countMassParent}] The number of bins in parent halo mass to use;
    \item[{\normalfont \ttfamily massParentMinimum}] The minimum parent halo mass to consider;
    \item[{\normalfont \ttfamily massParentMaximum}] The maximum parent halo mass to consider;
    \item[{\normalfont \ttfamily massRatioCount}] The number of bins in mass ratio to use;
    \item[{\normalfont \ttfamily massRatioMinimum}] The minimum mass ratio to consider;
    \item[{\normalfont \ttfamily massRatioMaximum}] The maximum mass ratio to consider;
    \item[{\normalfont \ttfamily redshiftsParent}] A list of redshifts at which to identify parent halos;
    \item[{\normalfont \ttfamily redshiftsProgenitor}] A corresponding list of redshifts at which to identify progenitor halos;
    \item[{\normalfont \ttfamily depthProgenitorPrimary}] The number of $i^\mathrm{th}$ most-massive progenitor mass functions to compute (starting from the 1$^\mathrm{st}$ most-massive);
    \item[{\normalfont \ttfamily depthHierarchySubhalo}] The maximum depth in the subhalo hierarchy for which to compute the unevolved subhalo mass function;
    \item[{\normalfont \ttfamily fractionTimeFormationRate}] The fraction of the current time over which to estimate the formation rate of halos when computing merger tree statistics;
    \item[{\normalfont \ttfamily nameGroupOutput}] The name of the \gls{hdf5} group to which mass functions will be written.
    \end{description}
    If the operator finds the named \gls{hdf5} group already in existence, it will accumulate its
    mass functions to those already written to the group, weighting by the inverse of the variance
    in each bin. The structure of the \gls{hdf5} group is as follows:
    \begin{verbatim}
    {
      DATASET "conditionalMassFunction" {
      COMMENT "Conditional mass functions []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Nratio, Nparent, Nz ) }
      }
      DATASET "conditionalMassFunctionError" {
      COMMENT "Conditional mass function errors []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Nratio, Nparent, Nz ) }
      }
      DATASET "conditionalMassFunctionCovariance" {
      COMMENT "Conditional mass function covariances []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Nratio, Nparent, Nz, Nratio, Nparent, Nz ) }
      }
      DATASET "formationRateFunction" {
      COMMENT "Formation rate functions []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( 2, Nratio, Nparent, Nz ) }
      }
      DATASET "formationRateFunctionError" {
      COMMENT "Formation rate function errors []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( 2, Nratio, Nparent, Nz ) }
      }
      DATASET "massParent" {
      COMMENT "Mass of parent node [Msolar]"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Nparent ) }
         ATTRIBUTE "unitsInSI" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SCALAR
            DATA {
            (0): 1.98892e+30
            }
          }
      }
      DATASET "massRatio" {
      COMMENT "Mass of ratio node [Msolar]"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Nratio ) }
         ATTRIBUTE "unitsInSI" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SCALAR
            DATA {
            (0): 1.98892e+30
            }
         }
      }
      DATASET "primaryProgenitorMassFunction" {
      COMMENT "Primary progenitor mass functions []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Ndepth, Nratio, Nparent, Nz ) }
      }
      DATASET "primaryProgenitorMassFunctionError" {
      COMMENT "Primary progenitor mass function errors []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Ndepth, Nratio, Nparent, Nz ) }
      }
      DATASET "redshiftParent" {
      COMMENT "Redshift of parent node []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Nz ) }
      }
      DATASET "redshiftProgenitor" {
      COMMENT "Redshift of progenitor node []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Nz ) }
      }
      DATASET "subhaloMassFunction" {
      COMMENT "Unevolved subhalo mass functions []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Nhierarchy, Nratio, Nparent ) }
      }
      DATASET "subhaloMassFunctionError" {
      COMMENT "Unevolved subhalo mass function errors []"
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( Nhierarchy, Nratio, Nparent ) }
      }
    }
    \end{verbatim}
    Where {\normalfont \ttfamily Nratio} is the number of bins in mass ratio, {\normalfont \ttfamily
    Nparent} is the number of bins in mass ratio, {\normalfont \ttfamily Nz} is the number of
    parent/progenitor redshift pairs, {\normalfont \ttfamily Ndepth} is the maximum depth in the
    ranking of most-massive progenitor mass functions, {\normalfont \ttfamily Nhierarchy} is the
    maximum depth in the subhalo hierarchy for subhalo mass functions. The first dimension of the
    {\normalfont \ttfamily formationRateFunction} dataset stores two different versions of the
    formation rate function. The first uses the mass of the forming halo at the time of formation,
    the second uses the mass of the node immediately prior to it becoming a subhalo.
  
    Mass functions are output as $\mathrm{d}N/\mathrm{d}\log_{10}m$ where $N$ is the number of halos per
    parent halo, and $m$ is the mass ratio.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorConditionalMF
     !!{
     A merger tree operator class which accumulates conditional mass functions for trees.
     !!}
     private
     class           (nbodyHaloMassErrorClass), pointer                             :: haloMassError_                       => null()
     class           (cosmologyFunctionsClass), pointer                             :: cosmologyFunctions_                  => null()
     double precision                         , allocatable, dimension(:          ) :: timeParents                                   , timeProgenitors                      , &
          &                                                                            redshiftsParent                               , redshiftsProgenitor                  , &
          &                                                                            massParents                                   , massRatios                           , &
          &                                                                            normalizationSubhaloMassFunction              , normalizationSubhaloMassFunctionError
     double precision                         , allocatable, dimension(:,:        ) :: normalization                                 , normalizationError
     double precision                         , allocatable, dimension(:,:,:,:    ) :: normalizationCovariance
     double precision                         , allocatable, dimension(:,:,:      ) :: conditionalMassFunction                       , conditionalMassFunctionError
     double precision                         , allocatable, dimension(:,:,:,:,:,:) :: conditionalMassFunctionCovariance          
     double precision                         , allocatable, dimension(:,:,:      ) :: subhaloMassFunction                           , subhaloMassFunctionError
     double precision                         , allocatable, dimension(:,:,:,:    ) :: primaryProgenitorMassFunction                 , primaryProgenitorMassFunctionError
     double precision                         , allocatable, dimension(:,:,:,:    ) :: formationRateFunction                         , formationRateFunctionError
     integer                                                                        :: countMassParent                               , depthProgenitorPrimary               , &
          &                                                                            massRatioCount                                , timeCount                            , &
          &                                                                            depthHierarchySubhalo                         , nodeHierarchyLevelMaximumID
     double precision                                                               :: massParentLogarithmicMinimum                  , massRatioLogarithmicMinimum          , &
          &                                                                            massParentLogarithmicBinWidthInverse          , massRatioLogarithmicBinWidthInverse  , &
          &                                                                            fractionTimeFormationRate                     , massParentMinimum                    , &
          &                                                                            massParentMaximum                             , massRatioMinimum                     , &
          &                                                                            massRatioMaximum
     logical                                                                        :: alwaysIsolatedHalosOnly                       , primaryProgenitorStatisticsValid     , &
          &                                                                            extendedStatistics                            , computeCovariances
     type            (varying_string         )                                      :: nameGroupOutput
     !$ integer      (omp_lock_kind          )                                      :: accumulateLock
   contains
     !![
     <methods>
       <method description="Compute weights for a halo in each bin of the mass function." method="binWeights" />
       <method description="Compute weights for a halo in each bin of a 2D mass function." method="binWeights2D" />
     </methods>
     !!]
     final     ::                        conditionalMFDestructor
     procedure :: operatePreEvolution => conditionalMFOperatePreEvolution
     procedure :: finalize            => conditionalMFFinalize
     procedure :: binWeights          => conditionalMFBinWeights
     procedure :: binWeights2D        => conditionalMFBinWeights2D
  end type mergerTreeOperatorConditionalMF

  interface mergerTreeOperatorConditionalMF
     !!{
     Constructors for the conditional mass function merger tree operator class.
     !!}
     module procedure conditionalMFConstructorParameters
     module procedure conditionalMFConstructorInternal
  end interface mergerTreeOperatorConditionalMF

contains

  function conditionalMFConstructorParameters(parameters) result(self)
    !!{
    Constructor for the conditional mass function merger tree operator class which takes a parameter set as input.
    !!}
    !$ use :: OMP_Lib          , only : OMP_Init_Lock
    implicit none
    type            (mergerTreeOperatorConditionalMF)                              :: self
    type            (inputParameters                ), intent(inout), target       :: parameters
    double precision                                 , allocatable  , dimension(:) :: redshiftsProgenitor      , redshiftsParent
    class           (cosmologyFunctionsClass        ), pointer                     :: cosmologyFunctions_
    class           (nbodyHaloMassErrorClass        ), pointer                     :: haloMassError_
    integer                                                                        :: countMassParent          , massRatioCount       , &
         &                                                                            depthProgenitorPrimary   , depthHierarchySubhalo
    double precision                                                               :: massParentMinimum        , massParentMaximum    , &
         &                                                                            massRatioMinimum         , massRatioMaximum     , &
         &                                                                            fractionTimeFormationRate
    logical                                                                        :: alwaysIsolatedHalosOnly  , extendedStatistics   , &
         &                                                                            computeCovariances
    type            (varying_string                 )                              :: nameGroupOutput

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="nbodyHaloMassError" name="haloMassError_"      source="parameters"/>
    <inputParameter>
      <name>countMassParent</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>The number of bins in parent mass when constructing conditional halo mass functions.</description>
    </inputParameter>
    <inputParameter>
      <name>massParentMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d10</defaultValue>
      <description>The minimum parent halo mass to bin when constructing conditional halo mass functions.</description>
    </inputParameter>
    <inputParameter>
      <name>massParentMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d15</defaultValue>
      <description>The maximum parent halo mass to bin when constructing conditional halo mass functions.</description>
    </inputParameter>
    <inputParameter>
      <name>massRatioCount</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>The number of bins in mass ratio when constructing conditional halo mass functions.</description>
    </inputParameter>
    <inputParameter>
      <name>massRatioMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d-4</defaultValue>
      <description>The minimum mass ratio to bin when constructing conditional halo mass functions.</description>
    </inputParameter>
    <inputParameter>
      <name>massRatioMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d1</defaultValue>
      <description>The maximum mass ratio to bin when constructing conditional halo mass functions.</description>
    </inputParameter>
    !!]
    allocate(redshiftsParent    (max(1,parameters%count('redshiftsParent'    ,zeroIfNotPresent=.true.))))
    !![
    <inputParameter>
      <name>redshiftsParent</name>
      <source>parameters</source>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The set of parent halo redshifts to use when constructing conditional halo mass functions.</description>
    </inputParameter>
    !!]
    allocate(redshiftsProgenitor(max(1,parameters%count('redshiftsProgenitor',zeroIfNotPresent=.true.))))
    !![
    <inputParameter>
      <name>redshiftsProgenitor</name>
      <source>parameters</source>
      <defaultValue>[1.0d0]</defaultValue>
      <description>The set of progenitor halo redshifts to use when constructing conditional halo mass functions.</description>
    </inputParameter>
    <inputParameter>
      <name>depthProgenitorPrimary</name>
      <source>parameters</source>
      <defaultValue>2</defaultValue>
      <description>The depth in progenitor ranking for which to store ranked progenitor mass functions. For example, a value of 2 means store mass functions for the most massive, and second most massive progenitor.</description>
    </inputParameter>
    <inputParameter>
      <name>depthHierarchySubhalo</name>
      <source>parameters</source>
      <defaultValue>2</defaultValue>
      <description>The depth in the subhalo hierarchy for which to store unevolved subhalo mass functions.</description>
    </inputParameter>
    <inputParameter>
      <name>fractionTimeFormationRate</name>
      <source>parameters</source>
      <defaultValue>0.01d0</defaultValue>
      <description>The fraction of the current time over which to estimate the formation rate of halos when computing merger tree statistics.</description>
    </inputParameter>
    <inputParameter>
      <name>alwaysIsolatedHalosOnly</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>Include only halos which have always been isolated when computing merger tree statistics?</description>
    </inputParameter>
    <inputParameter>
      <name>extendedStatistics</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>Compute extended statistics (formation rate function, $N^\mathrm{th}$ most-massive progenitor mass functions, etc.)?</description>
    </inputParameter>
    <inputParameter>
      <name>computeCovariances</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>Compute covariances for accumulated statistics?</description>
    </inputParameter>
    <inputParameter>
      <name>nameGroupOutput</name>
      <source>parameters</source>
      <defaultValue>var_str('conditionalMassFunction')</defaultValue>
      <description>The name of the HDF5 group to which the conditional mass function should be output.</description>
    </inputParameter>
    !!]
    ! Construct the instance.
    self=mergerTreeOperatorConditionalMF(                           &
         &                               countMassParent          , &
         &                               massParentMinimum        , &
         &                               massParentMaximum        , &
         &                               massRatioCount           , &
         &                               massRatioMinimum         , &
         &                               massRatioMaximum         , &
         &                               redshiftsParent          , &
         &                               redshiftsProgenitor      , &
         &                               depthProgenitorPrimary   , &
         &                               fractionTimeFormationRate, &
         &                               depthHierarchySubhalo    , &
         &                               alwaysIsolatedHalosOnly  , &
         &                               extendedStatistics       , &
         &                               computeCovariances       , &
         &                               nameGroupOutput          , &
         &                               cosmologyFunctions_      , &
         &                               haloMassError_             &
         &                              )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="haloMassError_"     />
    !!]
    return
  end function conditionalMFConstructorParameters

  function conditionalMFConstructorInternal(countMassParent,massParentMinimum,massParentMaximum,massRatioCount,massRatioMinimum,massRatioMaximum,redshiftsParent,redshiftsProgenitor,depthProgenitorPrimary,fractionTimeFormationRate,depthHierarchySubhalo,alwaysIsolatedHalosOnly,extendedStatistics,computeCovariances,nameGroupOutput,cosmologyFunctions_,haloMassError_) result(self)
    !!{
    Internal constructor for the conditional mass function merger tree operator class.
    !!}
    use :: Error            , only : Error_Report
    use :: Numerical_Ranges , only : Make_Range   , rangeTypeLogarithmic
    implicit none
    type            (mergerTreeOperatorConditionalMF)                              :: self
    double precision                                 , intent(in   ), dimension(:) :: redshiftsParent                 , redshiftsProgenitor
    integer                                          , intent(in   )               :: countMassParent                 , massRatioCount       , &
         &                                                                            depthProgenitorPrimary          , depthHierarchySubhalo
    double precision                                 , intent(in   )               :: massParentMinimum               , massParentMaximum    , &
         &                                                                            massRatioMinimum                , massRatioMaximum     , &
         &                                                                            fractionTimeFormationRate
    logical                                          , intent(in   )               :: alwaysIsolatedHalosOnly         , extendedStatistics   , &
         &                                                                            computeCovariances
    type            (varying_string                 ), intent(in   )               :: nameGroupOutput
    class           (cosmologyFunctionsClass        ), intent(in   ), target       :: cosmologyFunctions_
    class           (nbodyHaloMassErrorClass        ), intent(in   ), target       :: haloMassError_
    integer                                                                        :: i
    !![
    <constructorAssign variables="massParentMinimum, massParentMaximum, massRatioMinimum, massRatioMaximum, countMassParent, massRatioCount, depthProgenitorPrimary, depthHierarchySubhalo, fractionTimeFormationRate, alwaysIsolatedHalosOnly, extendedStatistics, computeCovariances, nameGroupOutput, *haloMassError_, *cosmologyFunctions_"/>
    !!]

    ! Store array sizes.
    self%timeCount                       =size(redshiftsParent)
    self%primaryProgenitorStatisticsValid=self%haloMassError_%errorZeroAlways()
    if (size(redshiftsProgenitor) /= self%timeCount) &
         & call Error_Report('mismatch in sizes of parent and progenitor redshift arrays'//{introspection:location})
    ! Allocate arrays.
    allocate(self%redshiftsParent    (self%timeCount        ))
    allocate(self%redshiftsProgenitor(self%timeCount        ))
    allocate(self%timeProgenitors    (self%timeCount        ))
    allocate(self%timeParents        (self%timeCount        ))
    allocate(self%massParents        (self%countMassParent+1))
    allocate(self%massRatios         (self%massRatioCount +1))
    allocate(                                          &
         &    self%normalization                       &
         &   (                                         &
         &    self%timeCount                         , &
         &    self%countMassParent                     &
         &   )                                         &
         &  )
    allocate(                                          &
         &    self%normalizationError                  &
         &   (                                         &
         &    self%timeCount                         , &
         &    self%countMassParent                     &
         &   )                                         &
         &  )
    allocate(                                          &
         &    self%conditionalMassFunction             &
         &   (                                         &
         &    self%timeCount                         , &
         &    self%countMassParent                   , &
         &    self%massRatioCount                      &
         &   )                                         &
         &  )
    allocate(                                          &
         &    self%conditionalMassFunctionError        &
         &   (                                         &
         &    self%timeCount                         , &
         &    self%countMassParent                   , &
         &    self%massRatioCount                      &
         &   )                                         &
         &  )
    if (self%computeCovariances) then
       allocate(                                          &
            &    self%conditionalMassFunctionCovariance   &
            &   (                                         &
            &    self%timeCount                         , &
            &    self%countMassParent                   , &
            &    self%massRatioCount                    , &
            &    self%timeCount                         , &
            &    self%countMassParent                   , &
            &    self%massRatioCount                      &
            &   )                                         &
            &  )
       allocate(                                          &
            &    self%normalizationCovariance             &
            &   (                                         &
            &    self%timeCount                         , &
            &    self%countMassParent                   , &
            &    self%timeCount                         , &
            &    self%countMassParent                     &
            &   )                                         &
            &  )
    end if
    if (self%extendedStatistics) then
       allocate(                                             &
            &    self%normalizationSubhaloMassFunction       &
            &   (                                            &
            &    self%countMassParent                        &
            &   )                                            &
            &  )
       allocate(                                             &
            &    self%normalizationSubhaloMassFunctionError  &
            &   (                                            &
            &    self%countMassParent                        &
            &   )                                            &
            &  )
       allocate(                                             &
            &    self%primaryProgenitorMassFunction          &
            &   (                                            &
            &    self%timeCount                            , &
            &    self%countMassParent                      , &
            &    self%massRatioCount                       , &
            &    self%depthProgenitorPrimary                 &
            &   )                                            &
            &  )
       allocate(                                             &
            &    self%primaryProgenitorMassFunctionError     &
            &   (                                            &
            &    self%timeCount                            , &
            &    self%countMassParent                      , &
            &    self%massRatioCount                       , &
            &    self%depthProgenitorPrimary                 &
            &   )                                            &
            &  )
       allocate(                                             &
            &    self%formationRateFunction                  &
            &   (                                            &
            &    self%timeCount                            , &
            &    self%countMassParent                      , &
            &    self%massRatioCount                       , &
            &    2                                           &
            &   )                                            &
            &  )
       allocate(                                             &
            &    self%formationRateFunctionError             &
            &   (                                            &
            &    self%timeCount                            , &
            &    self%countMassParent                      , &
            &    self%massRatioCount                       , &
            &    2                                           &
            &   )                                            &
            &  )
       allocate(                                             &
            &    self%subhaloMassFunction                    &
            &   (                                            &
            &    self%countMassParent                      , &
            &    self%massRatioCount                       , &
            &    self%depthHierarchySubhalo                  &
            &   )                                            &
            &  )
       allocate(                                             &
            &    self%subhaloMassFunctionError               &
            &   (                                            &
            &    self%countMassParent                      , &
            &    self%massRatioCount                       , &
            &    self%depthHierarchySubhalo                  &
            &   )                                            &
            &  )
    end if
    ! Construct bins for parent node mass.
    self%massParentLogarithmicMinimum        =  log( massParentMinimum)
    self%massParentLogarithmicBinWidthInverse= dble( countMassParent  ) &
         &                                    / log(                    &
         &                                          +massParentMaximum  &
         &                                          /massParentMinimum  &
         &                                         )
    self%massParents=Make_Range(                                        &
         &                      massParentMinimum                     , &
         &                      massParentMaximum                     , &
         &                      countMassParent                       , &
         &                      rangeType        =rangeTypeLogarithmic  &
         &                     )
    ! Construct bins for mass ratio.
    self%massRatioLogarithmicMinimum        =  log( massRatioMinimum)
    self%massRatioLogarithmicBinWidthInverse= dble( massRatioCount  ) &
         &                                   / log(                   &
         &                                         +massRatioMaximum  &
         &                                         /massRatioMinimum  &
         &                                        )
    self%massRatios=Make_Range(                                       &
         &                     massRatioMinimum                     , &
         &                     massRatioMaximum                     , &
         &                     massRatioCount                       , &
         &                     rangeType       =rangeTypeLogarithmic, &
         &                     rangeBinned     =.true.                &
         &                    )
    ! Construct arrays of times for progenitors.
    self%redshiftsProgenitor=redshiftsProgenitor
    do i=1,self%timeCount
       self%timeProgenitors(i)=                                      &
            & self%cosmologyFunctions_%cosmicTime(                   &
            &  self%cosmologyFunctions_%expansionFactorFromRedshift( &
            &   redshiftsProgenitor(i)                               &
            &  )                                                     &
            & )
    end do
    ! Construct arrays of times for parents.
    self%redshiftsParent=redshiftsParent
    do i=1,self%timeCount
       self%timeParents(i)=                                          &
            & self%cosmologyFunctions_%cosmicTime(                   &
            &  self%cosmologyFunctions_%expansionFactorFromRedshift( &
            &   redshiftsParent(i)                                   &
            &  )                                                     &
            & )
    end do
    ! Initialize mass function arrays.
    self   %normalization                        =0.0d0
    self   %normalizationError                   =0.0d0
    self   %conditionalMassFunction              =0.0d0
    self   %conditionalMassFunctionError         =0.0d0
    if (self%computeCovariances) then
       self%normalizationCovariance              =0.0d0
       self%conditionalMassFunctionCovariance    =0.0d0
    end if
    if (self%extendedStatistics) then
       self%normalizationSubhaloMassFunction     =0.0d0
       self%normalizationSubhaloMassFunctionError=0.0d0
       self%primaryProgenitorMassFunction        =0.0d0
       self%primaryProgenitorMassFunctionError   =0.0d0
       self%formationRateFunction                =0.0d0
       self%formationRateFunctionError           =0.0d0
       self%subhaloMassFunction                  =0.0d0
       self%subhaloMassFunctionError             =0.0d0
    end if
    ! Initialize OpenMP lock.
    !$ call OMP_Init_Lock(self%accumulateLock)
    ! Get required meta-property.
    !![
    <addMetaProperty component="basic" name="nodeHierarchyLevelMaximum" type="integer" id="self%nodeHierarchyLevelMaximumID"/>
    !!]
    return
  end function conditionalMFConstructorInternal

  subroutine conditionalMFDestructor(self)
    !!{
    Destructor for the merger tree operator function class.
    !!}
    !$ use :: OMP_Lib, only : OMP_Destroy_Lock
    implicit none
    type(mergerTreeOperatorConditionalMF), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%haloMassError_"     />
    !!]
    ! Destroy OpenMP lock.
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine conditionalMFDestructor

  subroutine conditionalMFOperatePreEvolution(self,tree)
    !!{
    Compute conditional mass function on {\normalfont \ttfamily tree}.
    !!}
    use    :: Error               , only : Error_Report
    use    :: Galacticus_Nodes    , only : mergerTree                   , nodeComponentBasic, treeNode
    use    :: Merger_Tree_Walkers , only : mergerTreeWalkerIsolatedNodes
    use    :: Numerical_Comparison, only : Values_Agree
    !$ use :: OMP_Lib             , only : OMP_Set_Lock                 , OMP_Unset_Lock
    implicit none
    class           (mergerTreeOperatorConditionalMF), intent(inout)                               , target :: self
    type            (mergerTree                     ), intent(inout)                               , target :: tree
    type            (treeNode                       ), pointer                                              :: node                   , nodeChild            , &
         &                                                                                                     nodeParent             , nodeParentChild      , &
         &                                                                                                     descendantNode         , nodeSibling
    type            (mergerTree                     ), pointer                                              :: treeCurrent
    class           (nodeComponentBasic             ), pointer                                              :: basic                  , basicChild           , &
         &                                                                                                     basicParent            , descendantBasic      , &
         &                                                                                                     basicParentChild       , basicSibling
    type            (mergerTreeWalkerIsolatedNodes  )                                                       :: treeWalker
    integer                                                                                                 :: i                      , binMassParent        , &
         &                                                                                                     binMassRatio           , iPrimary             , &
         &                                                                                                     jPrimary               , depthHierarchy       , &
         &                                                                                                     j                      , j2                   , &
         &                                                                                                     k2
    double precision                                                                                        :: branchBegin            , branchEnd            , &
         &                                                                                                     parentBranchBegin      , parentBranchEnd      , &
         &                                                                                                     massProgenitor         , massParent           , &
         &                                                                                                     branchMassInitial      , branchMassFinal      , &
         &                                                                                                     parentBranchMassInitial, parentBranchMassFinal, &
         &                                                                                                     massRatioLogarithmic   , timeUnevolved        , &
         &                                                                                                     massRatio              , massUnevolved
    logical                                                                                                 :: includeBranch
    double precision                                  , dimension(                                                                                             &
         &                                                        self%timeCount             ,                                                                 &
         &                                                        self%countMassParent       ,                                                                 &
         &                                                        self%depthProgenitorPrimary                                                                  &
         &                                                       )                                          :: primaryProgenitorMass
    double precision                                  , dimension(self%countMassParent                    ) :: weights1D
    double precision                                  , dimension(self%countMassParent,self%massRatioCount) :: weights2D

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Initialize primary progenitor masses to zero.
       primaryProgenitorMass=0.0d0
       ! Get root node of the tree.
       node => treeCurrent%nodeBase
       ! Accumulate normalization for subhalo mass function.
       if (self%extendedStatistics) then
          basic => node%basic()
          weights1D=self%binWeights(                                           &
               &                    basic%mass()                             , &
               &                    basic%time()                             , &
               &                    self%massParentLogarithmicMinimum        , &
               &                    self%massParentLogarithmicBinWidthInverse, &
               &                    self%countMassParent                       &
               &                   )
          !$ call OMP_Set_Lock(self%accumulateLock)
          self%normalizationSubhaloMassFunction     (:)=+self       %normalizationSubhaloMassFunction     (:) &
               &                                        +weights1D                                            &
               &                                        *treeCurrent%volumeWeight
          self%normalizationSubhaloMassFunctionError(:)=+self       %normalizationSubhaloMAssFunctionError(:) &
               &                                        +weights1D               **2                          &
               &                                        *treeCurrent%volumeWeight**2
          !$ call OMP_Unset_Lock(self%accumulateLock)
       end if
       ! Walk the tree, accumulating statistics.
       treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
       do while (treeWalker%next(node))
          ! Get the child node, and process if child exists.
          nodeChild => node%firstChild
          do while (associated(nodeChild))
             ! Check if child should be included.
             if (self%alwaysIsolatedHalosOnly) then
                basicChild    =>  nodeChild%basic                      (                                )
                includeBranch =  (basic    %integerRank0MetaPropertyGet(self%nodeHierarchyLevelMaximumID) == 0)
             else
                includeBranch =  .true.
             end if
             if (includeBranch) then
                ! Get the basic components.
                basic      => node     %basic()
                basicChild => nodeChild%basic()
                ! Determine range of times spanned by this branch.
                branchBegin=basicChild%time()
                branchEnd  =basic     %time()
                ! Iterate over times.
                do i=1,self%timeCount
                   ! Does the branch span a parent node time?
                   if     (                                                               &
                        &     branchBegin <= self%timeParents(i)                          &
                        &  .and.                                                          &
                        &   (                                                             &
                        &     branchEnd   >  self%timeParents(i)                          &
                        &    .or.                                                         &
                        &     (                                                           &
                        &       .not.associated(node%parent)                              &
                        &      .and.                                                      &
                        &       Values_Agree(branchEnd,self%timeParents(i),relTol=1.0d-6) &
                        &     )                                                           &
                        &   )                                                             &
                        & ) then
                      ! Get the masses on the branch.
                      branchMassInitial=basicChild%mass()
                      if (nodeChild%isPrimaryProgenitor()) then
                         branchMassFinal=basic%mass()
                         ! Remove the mass in any non-primary progenitors - we don't want to include
                         ! their mass in the estimated mass growth rate of this node.
                         nodeSibling => node%firstChild%sibling
                         do while (associated(nodeSibling))
                            basicSibling    => nodeSibling%basic()
                            branchMassFinal =  branchMassFinal-basicSibling%mass()
                            nodeSibling     => nodeSibling%sibling
                         end do
                         ! Do not let the parent mass decrease along the branch.
                         branchMassFinal=max(branchMassFinal,branchMassInitial)
                      else
                         branchMassFinal=branchMassInitial
                      end if
                      ! Interpolate to get the mass at the required time.
                      if (branchEnd == branchBegin) then
                         massParent=branchMassFinal
                      else
                         massParent=                     +branchMassInitial  &
                              &     +(branchMassFinal    -branchMassInitial) &
                              &     *(self%timeParents(i)-branchBegin      ) &
                              &     /(branchEnd          -branchBegin      )
                      end if
                      ! Find the bin to which this node accumulates.
                      weights1D     =self%binWeights(                                              &
                           &                         massParent                                  , &
                           &                         self%timeParents                         (i), &
                           &                         self%massParentLogarithmicMinimum           , &
                           &                         self%massParentLogarithmicBinWidthInverse   , &
                           &                         self%countMassParent                          &
                           &                        )
                      !$ call OMP_Set_Lock(self%accumulateLock)
                      self%normalization     (i,:)=+self       %normalization     (i,:) &
                           &                       +weights1D                           &
                           &                       *treeCurrent%volumeWeight
                      self%normalizationError(i,:)=+self       %normalizationError(i,:) &
                                    &              +weights1D               **2         &
                                    &              *treeCurrent%volumeWeight**2
                      if (self%computeCovariances) then
                         forall(j=1:self%countMassParent)
                            self%normalizationCovariance(i,j,i,:)=+self       %normalizationCovariance(i,j,i,:) &
                                 &                                +weights1D                          (  j    ) &
                                 &                                *weights1D                                    &
                                 &                                *treeCurrent%volumeWeight**2
                         end forall
                      end if
                      !$ call OMP_Unset_Lock(self%accumulateLock)
                   end if
                   ! Check if the branch spans the progenitor time.
                   if     (                                        &
                        &   branchBegin <= self%timeProgenitors(i) &
                        &  .and.                                   &
                        &   branchEnd   >  self%timeProgenitors(i) &
                        & ) then
                      ! Get the masses on the branch.
                      branchMassInitial=basicChild%mass()
                      if (nodeChild%isPrimaryProgenitor()) then
                         branchMassFinal=basic%mass()
                      else
                         branchMassFinal=branchMassInitial
                      end if
                      ! Interpolate to get the mass at the required time.
                      if (branchEnd == branchBegin) then
                         massProgenitor=branchMassFinal
                      else
                         massProgenitor=                         +branchMassInitial  &
                              &         +(branchMassFinal        -branchMassInitial) &
                              &         *(self%timeProgenitors(i)-branchBegin      ) &
                              &         /(branchEnd              -branchBegin      )
                      end if
                      ! Walk up the tree to find parents.
                      nodeParent => node
                      parentWalk : do while (associated(nodeParent))
                         ! Get the parent's child.
                         nodeParentChild   => nodeParent      %firstChild
                         ! Check if child should be included.
                         if (self%alwaysIsolatedHalosOnly) then
                            basicChild    =>  nodeParentChild%basic                 (                                )
                            includeBranch =  (basicChild     %integerRank0MetaPropertyGet(self%nodeHierarchyLevelMaximumID) == 0)
                         else
                            includeBranch =  .true.
                         end if
                         if (includeBranch) then
                            ! Get the basic components.
                            basicParent       => nodeParent      %basic     ()
                            basicParentChild  => nodeParentChild %basic     ()
                            ! Determine range of times spanned by this branch.
                            parentBranchBegin =  basicParentChild%time      ()
                            parentBranchEnd   =  basicParent     %time      ()
                            ! Does the branch span a parent node time?
                            if     (                                                                     &
                                 &   parentBranchBegin <= self%timeParents(i)                            &
                                 &  .and.                                                                &
                                 &   (                                                                   &
                                 &     parentBranchEnd >  self%timeParents(i)                            &
                                 &    .or.                                                               &
                                 &     (                                                                 &
                                 &       .not.associated(nodeParent%parent)                              &
                                 &      .and.                                                            &
                                 &       Values_Agree(parentBranchEnd,self%timeParents(i),relTol=1.0d-6) &
                                 &     )                                                                 &
                                 &   )                                                                   &
                                 & ) then
                               ! Get the masses on the parent branch.
                               parentBranchMassInitial=basicParentChild%mass()
                               parentBranchMassFinal  =     basicParent%mass()
                               ! Find the parent mass at the required time.
                               if (parentBranchEnd == parentBranchBegin) then
                                  massParent=parentBranchMassFinal
                               else
                                  massParent=                       +parentBranchMassInitial  &
                                       &     +(parentBranchMassFinal-parentBranchMassInitial) &
                                       &     *(self%timeParents(i)  -parentBranchBegin      ) &
                                       &     /(parentBranchEnd      -parentBranchBegin      )
                               end if
                               ! Accumulate to mass function array.
                               weights2D     =self%binWeights2D(                                              &
                                    &                           massParent                                  , &
                                    &                           self%timeParents                         (i), &
                                    &                           massProgenitor                              , &
                                    &                           self%timeProgenitors                     (i), &
                                    &                           self%massParentLogarithmicMinimum           , &
                                    &                           self%massParentLogarithmicBinWidthInverse   , &
                                    &                           self%countMassParent                        , &
                                    &                           self%massRatioLogarithmicMinimum            , &
                                    &                           self%massRatioLogarithmicBinWidthInverse    , &
                                    &                           self%massRatioCount                         , &
                                    &                           1                                             &
                                    &                          )
                               !$ call OMP_Set_Lock(self%accumulateLock)
                               self%conditionalMassFunction     (i,:,:)=+self          %conditionalMassFunction     (i,:,:) &
                                    &                                   +weights2D                                          &
                                    &                                   *treeCurrent   %volumeWeight
                               self%conditionalMassFunctionError(i,:,:)=+self          %conditionalMassFunctionError(i,:,:) &
                                    &                                   +weights2D                  **2                     &
                                    &                                   *treeCurrent   %volumeWeight**2
                               if (self%computeCovariances) then
                                  ! Accumulate only the upper triangle of the covariance matrix (technically, the upper triangle
                                  ! in the first mass ratio dimension) for speed. The lower triangle is constructed prior to
                                  ! output by simply copying the upper triangle.
                                  forall(j2=1:self%countMassParent)
                                     forall(k2=1:self%massRatioCount)
                                        self               %conditionalMassFunctionCovariance(i,:,k2:self%massRatioCount,i,j2,k2)= &
                                             & +self       %conditionalMassFunctionCovariance(i,:,k2:self%massRatioCount,i,j2,k2)  &
                                             & +weights2D                                    (  :,k2:self%massRatioCount        )  &
                                             & *weights2D                                    (                             j2,k2)  &
                                             & *treeCurrent%volumeWeight**2
                                     end forall
                                  end forall
                               end if
                               !$ call OMP_Unset_Lock(self%accumulateLock)
                               ! Check for formation.
                               if (self%extendedStatistics) then
                                  if (branchBegin > self%timeProgenitors(i)*(1.0d0-self%fractionTimeFormationRate) .and. .not.associated(nodeChild%firstChild)) then
                                     ! This is a newly formed halo, accumulate to formation rate arrays.
                                     weights2D     =self%binWeights2D(                                              &
                                          &                           massParent                                  , &
                                          &                           self%timeParents                         (i), &
                                          &                           branchMassInitial                           , &
                                          &                           basicChild%time                          ( ), &
                                          &                           self%massParentLogarithmicMinimum           , &
                                          &                           self%massParentLogarithmicBinWidthInverse   , &
                                          &                           self%countMassParent                        , &
                                          &                           self%massRatioLogarithmicMinimum            , &
                                          &                           self%massRatioLogarithmicBinWidthInverse    , &
                                          &                           self%massRatioCount                         , &
                                          &                           1                                             &
                                          &                          )
                                     !$ call OMP_Set_Lock(self%accumulateLock)
                                     self%formationRateFunction     (i,:,:,1)=+self          %formationRateFunction     (i,:,:,1) &
                                          &                                   +weights2D                                          &
                                          &                                   *treeCurrent   %volumeWeight
                                     self%formationRateFunctionError(i,:,:,1)=+self          %formationRateFunctionError(i,:,:,1) &
                                          &                                   +weights2D                  **2                     &
                                          &                                   *treeCurrent   %volumeWeight**2
                                     !$ call OMP_Unset_Lock(self%accumulateLock)
                                     ! Find the mass of this node just prior to it becoming a subhalo.
                                     descendantNode => node
                                     do while (associated(descendantNode%parent).and.associated(descendantNode%parent%firstChild,descendantNode))
                                        descendantNode => descendantNode%parent
                                     end do
                                     descendantBasic   => descendantNode%basic()
                                     weights2D     =self%binWeights2D(                                              &
                                          &                           massParent                                  , &
                                          &                           self%timeParents                         (i), &
                                          &                           descendantBasic%mass                     ( ), &
                                          &                           descendantBasic%time                     ( ), &
                                          &                           self%massParentLogarithmicMinimum           , &
                                          &                           self%massParentLogarithmicBinWidthInverse   , &
                                          &                           self%countMassParent                        , &
                                          &                           self%massRatioLogarithmicMinimum            , &
                                          &                           self%massRatioLogarithmicBinWidthInverse    , &
                                          &                           self%massRatioCount                         , &
                                          &                           1                                             &
                                          &                          )
                                     !$ call OMP_Set_Lock(self%accumulateLock)
                                     self%formationRateFunction     (i,:,:,2)=+self          %formationRateFunction     (i,:,:,2) &
                                          &                                   +weights2D                                          &
                                          &                                   *treeCurrent   %volumeWeight
                                     self%formationRateFunctionError(i,:,:,2)=+self          %formationRateFunctionError(i,:,:,2) &
                                          &                                   +weights2D                  **2                     &
                                          &                                   *treeCurrent   %volumeWeight**2
                                     !$ call OMP_Unset_Lock(self%accumulateLock)
                                  end if
                                  ! Accumulate to the primary progenitor mass array if necessary.
                                  if (self%primaryProgenitorStatisticsValid) then
                                     binMassParent=int(                                                      &
                                          &            +(+log(massParent)-self%massParentLogarithmicMinimum) &
                                          &            *self%massParentLogarithmicBinWidthInverse            &
                                          &           )                                                      &
                                          &        +1
                                     if (binMassParent >= 1 .and. binMassParent <= self%countMassParent) then
                                        massRatio=massProgenitor/massParent
                                        iPrimary =1
                                        do while (massRatio < primaryProgenitorMass(i,binMassParent,iPrimary))
                                           iPrimary=iPrimary+1
                                           if (iPrimary > self%depthProgenitorPrimary) exit
                                        end do
                                        if (iPrimary <= self%depthProgenitorPrimary) then
                                           if (iPrimary < self%depthProgenitorPrimary) then
                                              do jPrimary=self%depthProgenitorPrimary,iPrimary+1,-1
                                                 primaryProgenitorMass        (i,binMassParent,jPrimary  ) &
                                                      & =primaryProgenitorMass(i,binMassParent,jPrimary-1)
                                              end do
                                           end if
                                           primaryProgenitorMass(i,binMassParent,iPrimary)=massRatio
                                        end if
                                     end if
                                  end if
                               end if
                            end if
                         end if
                         nodeParent => nodeParent%parent
                      end do parentWalk
                   end if
                   ! Record the mass of the branch at the parent time.
                end do
                ! Accumulate unevolved subhalo mass functions.
                if (self%extendedStatistics) then
                   if (.not.associated(nodeChild%firstChild)) then
                      ! This is a branch tip. Follow until it to the final time, storing its mass just prior to becoming a subhalo, and
                      ! the hierarchy depth.
                      depthHierarchy =  0
                      descendantNode => nodeChild
                      do while (associated(descendantNode).and.depthHierarchy <= self%depthHierarchySubhalo)
                         if (associated(descendantNode%parent).and..not.descendantNode%isPrimaryProgenitor()) then
                            depthHierarchy=depthHierarchy+1
                            if (depthHierarchy == 1) then
                               descendantBasic => descendantNode %basic()
                               massUnevolved   =  descendantBasic%mass ()
                               timeUnevolved   =  descendantBasic%time ()
                            end if
                         end if
                         descendantNode => descendantNode%parent
                      end do
                      if (depthHierarchy > 0 .and. depthHierarchy <= self%depthHierarchySubhalo) then
                         basicParent => treeCurrent%nodeBase%basic()
                         weights2D     =self%binWeights2D(                                           &
                              &                           basicParent%mass()                       , &
                              &                           basicParent%time()                       , &
                              &                           massUnevolved                            , &
                              &                           timeUnevolved                            , &
                              &                           self%massParentLogarithmicMinimum        , &
                              &                           self%massParentLogarithmicBinWidthInverse, &
                              &                           self%countMassParent                     , &
                              &                           self%massRatioLogarithmicMinimum         , &
                              &                           self%massRatioLogarithmicBinWidthInverse , &
                              &                           self%massRatioCount                      , &
                              &                           1                                          &
                              &                          )
                         !$ call OMP_Set_Lock(self%accumulateLock)
                         self%subhaloMassFunction     (:,:,depthHierarchy)=+self          %subhaloMassFunction     (:,:,depthHierarchy) &
                              &                                            +weights2D                                                   &
                              &                                            *treeCurrent   %volumeWeight
                         self%subhaloMassFunctionError(:,:,depthHierarchy)=+self          %subhaloMassFunctionError(:,:,depthHierarchy) &
                              &                                            +weights2D                  **2                              &
                              &                                            *treeCurrent   %volumeWeight**2
                         !$ call OMP_Unset_Lock(self%accumulateLock)
                      end if
                   end if
                end if
             end if
             ! Move to the next child.
             nodeChild => nodeChild%sibling
          end do
       end do
       ! Store the computed primary progenitor mass functions.
       if (self%extendedStatistics.and.self%primaryProgenitorStatisticsValid) then
          do i=1,self%timeCount
             do binMassParent=1,self%countMassParent
                do iPrimary=1,self%depthProgenitorPrimary
                   if (primaryProgenitorMass(i,binMassParent,iPrimary) > 0.0d0) then
                      massRatioLogarithmic=log(primaryProgenitorMass(i,binMassParent,iPrimary))
                      binMassRatio         =int(                                                           &
                           &                     (massRatioLogarithmic-self%massRatioLogarithmicMinimum)   &
                           &                    *                                                          &
                           &                     self%massRatioLogarithmicBinWidthInverse                  &
                           &                   )                                                           &
                           &                +1
                      if (binMassRatio  >= 1 .and. binMassRatio  <= self%massRatioCount) then
                         !$ call OMP_Set_Lock(self%accumulateLock)
                         self          %primaryProgenitorMassFunction     (i,binMassParent,binMassRatio,iPrimary)= &
                              & +  self%primaryProgenitorMassFunction     (i,binMassParent,binMassRatio,iPrimary)  &
                              & +       primaryProgenitorMass             (i,binMassParent,             iPrimary)  &
                              & *treeCurrent%volumeWeight
                         self          %primaryProgenitorMassFunctionError(i,binMassParent,binMassRatio,iPrimary)= &
                              &    self%primaryProgenitorMassFunctionError(i,binMassParent,binMassRatio,iPrimary)  &
                              & +(                                                                                 &
                              &   +     primaryProgenitorMass             (i,binMassParent,             iPrimary)  &
                              &   *treeCurrent%volumeWeight                                                        &
                              &  )**2
                         !$ call OMP_Unset_Lock(self%accumulateLock)
                      end if
                   end if
                end do
             end do
          end do
       end if
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine conditionalMFOperatePreEvolution

  function conditionalMFBinWeights(self,mass,time,massLogarithmicMinimumBins,massLogarithmicWidthInverseBins,countBins)
    !!{
    Computes the weight that a given halo contributes to an array of bins.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (mergerTreeOperatorConditionalMF), intent(inout)        :: self
    double precision                                 , intent(in   )        :: mass                      , time                           , &
         &                                                                     massLogarithmicMinimumBins, massLogarithmicWidthInverseBins
    integer                                          , intent(in   )        :: countBins
    double precision                                 , dimension(countBins) :: conditionalMFBinWeights
    type            (treeNode                       ), pointer              :: node
    class           (nodeComponentBasic             ), pointer              :: basic
    double precision                                                        :: massError
    integer                                                                 :: i

    ! Construct a node and find the mass error.
    node  => treeNode      (                 )
    basic => node    %basic(autoCreate=.true.)
    call basic%massSet(mass)
    call basic%timeSet(time)
    massError=basic%mass()*self%haloMassError_%errorFractional(node)
    call node%destroy()
    deallocate(node)
    ! Handle zero errors.
    if (massError <= 0.0d0) then
       conditionalMFBinWeights=0.0d0
       i                      =int(                                         &
            &                      +(+log(mass)-massLogarithmicMinimumBins) &
            &                      *massLogarithmicWidthInverseBins         &
            &                     )                                         &
            &                  +1
       if (i >= 1 .and. i <= countBins) conditionalMFBinWeights(i)=1.0d0
    else
       ! Find the contribution to each bin.
       forall(i=1:countBins)
          conditionalMFBinWeights(i)=+(                                             &
               &                       +erf(                                        &
               &                            +(                                      &
               &                              +exp(                                 &
               &                                   +massLogarithmicMinimumBins      &
               &                                   +dble(i  )                       &
               &                                   /massLogarithmicWidthInverseBins &
               &                                  )                                 &
               &                              -mass                                 &
               &                             )                                      &
               &                            /sqrt(2.0d0)                            &
               &                            /massError                              &
               &                           )                                        &
               &                       -erf(                                        &
               &                            +(                                      &
               &                              +exp(                                 &
               &                                   +massLogarithmicMinimumBins      &
               &                                   +dble(i-1)                       &
               &                                   /massLogarithmicWidthInverseBins &
               &                                  )                                 &
               &                              -mass                                 &
               &                             )                                      &
               &                            /sqrt(2.0d0)                            &
               &                            /massError                              &
               &                           )                                        &
               &                      )                                             &
               &                     /2.0d0
       end forall
    end if
    return
  end function conditionalMFBinWeights

  function conditionalMFBinWeights2D(self,mass1,time1,mass2,time2,massLogarithmicMinimumBins1,massLogarithmicWidthInverseBins1,countBins1,massRatioLogarithmicMinimumBins2,massRatioLogarithmicWidthInverseBins2,countBins2,moment)
    !!{
    Computes the weight that a given halo contributes to a 2D array of bins.
    !!}
    use :: Error                , only : Error_Report
    use :: Galacticus_Nodes     , only : nodeComponentBasic     , treeNode
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (mergerTreeOperatorConditionalMF), intent(inout)                    :: self
    double precision                                 , intent(in   )                    :: mass1                                  , time1                                , &
         &                                                                                 massLogarithmicMinimumBins1            , massLogarithmicWidthInverseBins1     , &
         &                                                                                 mass2                                  , time2                                , &
         &                                                                                 massRatioLogarithmicMinimumBins2       , massRatioLogarithmicWidthInverseBins2
    integer                                          , intent(in   )                    :: countBins1                             , countBins2                           , &
         &                                                                                 moment
    double precision                                 , dimension(countBins1,countBins2) :: conditionalMFBinWeights2D
    type            (treeNode                       ), pointer                          :: node1                                  , node2
    class           (nodeComponentBasic             ), pointer                          :: basic1                                 , basic2
    double precision                                 , parameter                        :: integrationExtent               =10.0d0
    double precision                                                                    :: massError1                             , massError2                           , &
         &                                                                                 mass1LowerLimit                        , mass1UpperLimit                      , &
         &                                                                                 mass2LowerLimit                        , mass2UpperLimit                      , &
         &                                                                                 correlation                            , massError2Reduced
    integer                                                                             :: i                                      , j
    type            (integrator                     )                                   :: integrator_

    ! Validate moment.
    if (moment < 0 .or. moment > 2) call Error_Report('moment must be 0, 1, or 2'//{introspection:location})
    ! Construct nodes and find the mass errors and their correlation. Given the correlation coefficient, C₁₂, between the mass
    ! errors, σ₁ and σ₂, then once M₁ is fixed, M₂ is shifted by C₁₂ (σ₂/σ₁) (M₁-<M₁>), and has remaining variance (1-C₁₂²)σ₂²
    node1  => treeNode       (                 )
    node2  => treeNode       (                 )
    basic1 => node1    %basic(autoCreate=.true.)
    basic2 => node2    %basic(autoCreate=.true.)
    call basic1%massSet(mass1)
    call basic1%timeSet(time1)
    call basic2%massSet(mass2)
    call basic2%timeSet(time2)
    massError1       =basic1%mass()*self%haloMassError_%errorFractional(node1      )
    massError2       =basic2%mass()*self%haloMassError_%errorFractional(      node2)
    correlation      =              self%haloMassError_%correlation    (node1,node2)
    massError2Reduced=+massError2                 &
         &            *sqrt(1.0d0-correlation**2)
    call node1%destroy()
    call node2%destroy()
    deallocate(node1)
    deallocate(node2)
    ! Handle zero errors.
    if (massError1 <= 0.0d0 .or. massError2 <= 0.0d0) then
       ! We currently do not handle cases where only one error is zero.
       if (massError1 > 0.0d0 .or. massError2 > 0.0d0) call Error_Report('both mass errors must be zero or both must be non-zero'//{introspection:location})
       ! Find the bin contributed to.
       i      =int(                                                      &
            &      +(+log(      mass1)-massLogarithmicMinimumBins1     ) &
            &      *massLogarithmicWidthInverseBins1                     &
            &     )                                                      &
            &  +1
       j      =int(                                                      &
            &      +(+log(mass2/mass1)-massRatioLogarithmicMinimumBins2) &
            &      *massRatioLogarithmicWidthInverseBins2                &
            &     )                                                      &
            &  +1
       ! Accumulate to that bin.
       conditionalMFBinWeights2D(:,:)=0.0d0
       if     (                                    &
            &   i >= 1 .and. i <= countBins1       &
            &  .and.                               &
            &   j >= 1 .and. j <= countBins2       &
            & )                                    &
            & conditionalMFBinWeights2D(i,j)=(mass2/mass1)**moment
    else
       ! Find the contribution to each bin.
       integrator_=integrator(conditionalMFBinWeights2DIntegrand,toleranceAbsolute=1.0d-10,toleranceRelative=1.0d-03)
       do i=1,countBins1
          mass1LowerLimit=max(                                        &
               &              +exp(                                   &
               &                   +massLogarithmicMinimumBins1       &
               &                   +dble(i-1)                         &
               &                   /massLogarithmicWidthInverseBins1  &
               &                  )                                 , &
               &              +mass1                                  &
               &              -integrationExtent                      &
               &              *massError1                             &
               &             )
          mass1UpperLimit=min(                                        &
               &              +exp(                                   &
               &                   +massLogarithmicMinimumBins1       &
               &                   +dble(i  )                         &
               &                   /massLogarithmicWidthInverseBins1  &
               &                  )                                 , &
               &              +mass1                                  &
               &              +integrationExtent                      &
               &              *massError1                             &
               &             )
          if (mass1UpperLimit > mass1LowerLimit) then
             do j=1,countBins2
                mass2LowerLimit=+exp(                                       &
                     &               +massRatioLogarithmicMinimumBins2      &
                     &               +dble(j-1)                             &
                     &               /massRatioLogarithmicWidthInverseBins2 &
                     &              )                                       &
                     &          *mass1LowerLimit                            &
                     &          -integrationExtent                          &
                     &          *massError2                                 &
                     &          +correlation                                &
                     &          *(                                          &
                     &            +massError2                               &
                     &            /massError1                               &
                     &           )                                          &
                     &          *(                                          &
                     &            +mass1LowerLimit                          &
                     &            -mass1                                    &
                     &           )
                mass2UpperLimit=+exp(                                       &
                     &               +massRatioLogarithmicMinimumBins2      &
                     &               +dble(j  )                             &
                     &               /massRatioLogarithmicWidthInverseBins2 &
                     &              )                                       &
                     &          *mass1UpperLimit                            &
                     &          +integrationExtent                          &
                     &          *massError2                                 &
                     &          +correlation                                &
                     &          *(                                          &
                     &            +massError2                               &
                     &            /massError1                               &
                     &           )                                          &
                     &          *(                                          &
                     &            +mass1UpperLimit                          &
                     &            -mass1                                    &
                     &           )
                if     (                         &
                     &   mass2LowerLimit > mass2 &
                     &  .or.                     &
                     &   mass2UpperLimit < mass2 &
                     & ) then
                   conditionalMFBinWeights2D(i,j)=0.0d0
                else
                   conditionalMFBinWeights2D(i,j)=max(                                                        &
                        &                             integrator_%integrate(mass1LowerLimit,mass1UpperLimit), &
                        &                             0.0d0                                                   &
                        &                            )
                end if
             end do
          else
             conditionalMFBinWeights2D(i,:)=0.0d0
          end if
       end do
    end if
    return

  contains

    double precision function conditionalMFBinWeights2DIntegrand(mass1Primed)
      !!{
      Integrand used in finding the weight given to a bin in the space of parent mass vs. progenitor mass ratio.
      !!}
      use :: Error                   , only : Error_Report
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: mass1Primed
      double precision                :: mass2LowerLimit, mass2UpperLimit, &
           &                             mass2Shifted

      mass2LowerLimit                   =+exp(                                       &
           &                                  +massRatioLogarithmicMinimumBins2      &
           &                                  +dble(j-1)                             &
           &                                  /massRatioLogarithmicWidthInverseBins2 &
           &                                 )                                       &
           &                             *mass1Primed
      mass2UpperLimit                   =+exp(                                       &
           &                                  +massRatioLogarithmicMinimumBins2      &
           &                                  +dble(j  )                             &
           &                                  /massRatioLogarithmicWidthInverseBins2 &
           &                                 )                                       &
           &                             *mass1Primed
      ! Find the shifted mass2.
      mass2Shifted=+mass2         &
           &       +correlation   &
           &       *(             &
           &         +massError2  &
           &         /massError1  &
           &        )             &
           &       *(             &
           &         +mass1Primed &
           &         -mass1       &
           &        )
      ! Evaluate integrand for the relevant moment.
      select case (moment)
      case (0)
         conditionalMFBinWeights2DIntegrand=+(                        &
              &                               +erf(                   &
              &                                    +(                 &
              &                                      +mass2UpperLimit &
              &                                      -mass2Shifted    &
              &                                     )                 &
              &                                    /sqrt(2.0d0)       &
              &                                    /massError2Reduced &
              &                                   )                   &
              &                               -erf(                   &
              &                                    +(                 &
              &                                      +mass2LowerLimit &
              &                                      -mass2Shifted    &
              &                                     )                 &
              &                                    /sqrt(2.0d0)       &
              &                                    /massError2Reduced &
              &                                   )                   &
              &                              )                        &
              &                             /2.0d0
      case (1)
         conditionalMFBinWeights2DIntegrand=+(                          &
              &                               +mass2Shifted             &
              &                               *erf(                     &
              &                                    +(                   &
              &                                      +mass2UpperLimit   &
              &                                      -mass2Shifted      &
              &                                     )                   &
              &                                    /sqrt(2.0d0)         &
              &                                    /massError2Reduced   &
              &                                   )                     &
              &                               /2.0d0                    &
              &                               -massError2Reduced        &
              &                               *exp(                     &
              &                                    -(                   &
              &                                      +(                 &
              &                                        +mass2UpperLimit &
              &                                        -mass2Shifted    &
              &                                       )                 &
              &                                      /massError2Reduced &
              &                                     )**2                &
              &                                    /2.0d0               &
              &                                   )                     &
              &                               /sqrt(                    &
              &                                     +2.0d0              &
              &                                     *Pi                 &
              &                                    )                    &
              &                               -mass2Shifted             &
              &                               *erf(                     &
              &                                    +(                   &
              &                                      +mass2LowerLimit   &
              &                                      -mass2Shifted      &
              &                                     )                   &
              &                                    /sqrt(2.0d0)         &
              &                                    /massError2Reduced   &
              &                                   )                     &
              &                               /2.0d0                    &
              &                               +massError2Reduced        &
              &                               *exp(                     &
              &                                    -(                   &
              &                                      +(                 &
              &                                        +mass2LowerLimit &
              &                                        -mass2Shifted    &
              &                                       )                 &
              &                                      /massError2Reduced &
              &                                     )**2                &
              &                                    /2.0d0               &
              &                                   )                     &
              &                               /sqrt(                    &
              &                                     +2.0d0              &
              &                                     *Pi                 &
              &                                    )                    &
              &                              )
      case (2)
         conditionalMFBinWeights2DIntegrand=+(                          &
              &                               -massError2Reduced        &
              &                               *(                        &
              &                                 +mass2UpperLimit        &
              &                                 +mass2Shifted           &
              &                                )                        &
              &                               *exp(                     &
              &                                    -(                   &
              &                                      +(                 &
              &                                        +mass2UpperLimit &
              &                                        -mass2Shifted    &
              &                                       )                 &
              &                                      /massError2Reduced &
              &                                     )**2                &
              &                                    /2.0d0               &
              &                                   )                     &
              &                               /sqrt(                    &
              &                                     +2.0d0              &
              &                                     *Pi                 &
              &                                    )                    &
              &                               +(                        &
              &                                 +massError2Reduced**2   &
              &                                 +mass2Shifted     **2   &
              &                                )                        &
              &                               *erf(                     &
              &                                    +(                   &
              &                                      +mass2UpperLimit   &
              &                                      -mass2Shifted      &
              &                                     )                   &
              &                                    /sqrt(2.0d0)         &
              &                                    /massError2Reduced   &
              &                                   )                     &
              &                               /2.0d0                    &
              &                               +massError2Reduced        &
              &                               *(                        &
              &                                 +mass2LowerLimit        &
              &                                 +mass2Shifted           &
              &                                )                        &
              &                               *exp(                     &
              &                                    -(                   &
              &                                      +(                 &
              &                                        +mass2LowerLimit &
              &                                        -mass2Shifted    &
              &                                       )                 &
              &                                      /massError2Reduced &
              &                                     )**2                &
              &                                    /2.0d0               &
              &                                   )                     &
              &                               /sqrt(                    &
              &                                     +2.0d0              &
              &                                     *Pi                 &
              &                                    )                    &
              &                               -(                        &
              &                                 +massError2Reduced**2   &
              &                                 +mass2Shifted     **2   &
              &                                )                        &
              &                               *erf(                     &
              &                                    +(                   &
              &                                      +mass2LowerLimit   &
              &                                      -mass2Shifted      &
              &                                     )                   &
              &                                    /sqrt(2.0d0)         &
              &                                    /massError2Reduced   &
              &                                   )                     &
              &                               /2.0d0                    &
              &                              )
      case default
         conditionalMFBinWeights2DIntegrand=0.0d0
         call Error_Report('moment not supported'//{introspection:location})
      end select
      conditionalMFBinWeights2DIntegrand=+conditionalMFBinWeights2DIntegrand &
           &                             *exp(                               &
           &                                  -0.5d0                         &
           &                                  *(                             &
           &                                    +(                           &
           &                                      +mass1Primed               &
           &                                      -mass1                     &
           &                                     )                           &
           &                                    /massError1                  &
           &                                   )**2                          &
           &                                 )                               &
           &                             /sqrt(                              &
           &                                   +2.0d0                        &
           &                                   *Pi                           &
           &                                  )                              &
           &                             /massError1                         &
           &                             /mass1Primed**moment
      return
    end function conditionalMFBinWeights2DIntegrand

  end function conditionalMFBinWeights2D

  subroutine conditionalMFFinalize(self)
    !!{
    Outputs conditional mass function.
    !!}
    use    :: Output_HDF5                     , only : outputFile
    use    :: HDF5_Access                     , only : hdf5Access
    use    :: IO_HDF5                         , only : hdf5Object
    use    :: Numerical_Constants_Astronomical, only : massSolar
    !$ use :: OMP_Lib                         , only : OMP_Get_Num_Threads
    implicit none
    class           (mergerTreeOperatorConditionalMF), intent(inout)                         :: self
    type            (hdf5Object                     )                                        :: conditionalMassFunctionGroup       , massDataset
    double precision                                 , allocatable  , dimension(:          ) :: normalizationSubhaloMassFunction   , normalizationSubhaloMassFunctionError
    double precision                                 , allocatable  , dimension(:,:        ) :: normalization                      , normalizationError
    double precision                                 , allocatable  , dimension(:,:,:      ) :: conditionalMassFunction            , conditionalMassFunctionError
    double precision                                 , allocatable  , dimension(:,:,:      ) :: subhaloMassFunction                , subhaloMassFunctionError
    double precision                                 , allocatable  , dimension(:,:,:,:    ) :: primaryProgenitorMassFunction      , primaryProgenitorMassFunctionError
    double precision                                 , allocatable  , dimension(:,:,:,:    ) :: formationRateFunction              , formationRateFunctionError           , &
         &                                                                                      normalizationCovariance
    double precision                                 , allocatable  , dimension(:,:,:,:,:,:) :: conditionalMassFunctionCovariance
    integer                                                                                  :: i                                  , j                                    , &
         &                                                                                      iPrimary                           , k                                    , &
         &                                                                                      j1                                 , k1                                   , &
         &                                                                                      accumulationCount                  , j2

    ! Compute normalizations.
    self   %normalization                   =self%normalization                        /self%massRatioLogarithmicBinWidthInverse   /log(10.0d0)
    self   %normalizationError              =self%normalizationError                   /self%massRatioLogarithmicBinWidthInverse**2/log(10.0d0)**2
    if (self%computeCovariances) then
       self%normalizationCovariance         =self%normalizationCovariance              /self%massRatioLogarithmicBinWidthInverse**2/log(10.0d0)**2
    end if
    if (self%extendedStatistics) then
       self%normalizationSubhaloMassFunction=self%normalizationSubhaloMassFunction     /self%massRatioLogarithmicBinWidthInverse   /log(10.0d0)
       self%normalizationSubhaloMassFunction=self%normalizationSubhaloMassFunctionError/self%massRatioLogarithmicBinWidthInverse**2/log(10.0d0)**2
    end if
    ! Populate lower triangles of covariance matrices. We consider only the diagonal blocks of the output time dimensions, since
    ! throughout this module we assume no correlation between times. Also note that the indices on the parent mass dimensions must
    ! be switched in this assignment to ensure that we copy the correct transposed section of the covariance matrix.
    if (self%computeCovariances) then
       forall(i=1:self%timeCount)
          forall(j1=1:self%countMassParent)
             forall(j2=1:self%countMassParent)
                forall(k1=2:self%massRatioCount)
                   self%conditionalMassFunctionCovariance(i,j1,1:k1-1,i,j2,k1)=self%conditionalMassFunctionCovariance(i,j2,k1,i,j1,1:k1-1)
                end forall
             end forall
          end forall
       end forall
    end if
    ! Output the data.
    !$ call hdf5Access%set()
    ! Check if our output group already exists.
    if (outputFile%hasGroup(char(self%nameGroupOutput))) then
       ! Our group does exist. Read existing mass functions, add them to our own, then write back to file.
       conditionalMassFunctionGroup=outputFile%openGroup(char(self%nameGroupOutput),'Conditional mass functions of merger trees.',objectsOverwritable=.true.,overwriteOverride=.true.)
       allocate(normalization               ,mold=self%normalization               )
       allocate(normalizationError          ,mold=self%normalizationError          )
       allocate(conditionalMassFunctionError,mold=self%conditionalMassFunctionError)
       call conditionalMassFunctionGroup%readAttribute('accumulationCount'           ,accumulationCount           )
       call conditionalMassFunctionGroup%readDataset  ('normalization'               ,normalization               )
       call conditionalMassFunctionGroup%readDataset  ('normalizationError'          ,normalizationError          )
       call conditionalMassFunctionGroup%readDataset  ('conditionalMassFunction'     ,conditionalMassFunction     )
       call conditionalMassFunctionGroup%readDataset  ('conditionalMassFunctionError',conditionalMassFunctionError)
       if ( self%computeCovariances) then
          allocate(normalizationCovariance          ,mold=self%normalizationCovariance          )
          allocate(conditionalMassFunctionCovariance,mold=self%conditionalMassFunctionCovariance)
          call conditionalMassFunctionGroup%readDataset('normalizationCovariance'          ,normalizationCovariance          )
          call conditionalMassFunctionGroup%readDataset('conditionalMassFunctionCovariance',conditionalMassFunctionCovariance)
       end if
       if (self%extendedStatistics) then
          allocate(normalizationSubhaloMassFunction     ,mold=self%normalizationSubhaloMassFunction     )
          allocate(normalizationSubhaloMassFunctionError,mold=self%normalizationSubhaloMassFunctionError)
          allocate(primaryProgenitorMassFunction        ,mold=self%primaryProgenitorMassFunction        )
          allocate(primaryProgenitorMassFunctionError   ,mold=self%primaryProgenitorMassFunctionError   )
          allocate(formationRateFunction                ,mold=self%formationRateFunction                )
          allocate(formationRateFunctionError           ,mold=self%formationRateFunctionError           )
          allocate(subhaloMassFunction                  ,mold=self%subhaloMassFunction                  )
          allocate(subhaloMassFunctionError             ,mold=self%subhaloMassFunctionError             )
          call conditionalMassFunctionGroup%readDataset('normalizationSubhaloMassFunction'     ,normalizationSubhaloMassFunction     )
          call conditionalMassFunctionGroup%readDataset('normalizationSubhaloMassFunctionError',normalizationSubhaloMassFunctionError)
          call conditionalMassFunctionGroup%readDataset('primaryProgenitorMassFunction'        ,primaryProgenitorMassFunction        )
          call conditionalMassFunctionGroup%readDataset('primaryProgenitorMassFunctionError'   ,primaryProgenitorMassFunctionError   )
          call conditionalMassFunctionGroup%readDataset('formationRateFunction'                ,formationRateFunction                )
          call conditionalMassFunctionGroup%readDataset('formationRateFunctionError'           ,formationRateFunctionError           )
          call conditionalMassFunctionGroup%readDataset('subhaloMassFunction'                  ,subhaloMassFunction                  )
          call conditionalMassFunctionGroup%readDataset('subhaloMassFunctionError'             ,subhaloMassFunctionError             )
       end if
       ! Accumulate the conditional mass functions.
       self       %conditionalMassFunction          =self%conditionalMassFunction           +conditionalMassFunction
       self       %conditionalMassFunctionError     =self%conditionalMassFunctionError      +conditionalMassFunctionError
       if (self%computeCovariances) then
          self    %conditionalMassFunctionCovariance=self%conditionalMassFunctionCovariance +conditionalMassFunctionCovariance
       end if
       if (self%extendedStatistics) then
          self    %formationRateFunction            =self%formationRateFunction             +formationRateFunction
          self    %formationRateFunctionError       =self%formationRateFunctionError        +formationRateFunctionError
          self    %subhaloMassFunction              =self%subhaloMassFunction               +subhaloMassFunction
          self    %subhaloMassFunctionError         =self%subhaloMassFunctionError          +subhaloMassFunctionError
          if (self%primaryProgenitorStatisticsValid) then
             self%primaryProgenitorMassFunction     =self%primaryProgenitorMassFunction     +primaryProgenitorMassFunction
             self%primaryProgenitorMassFunctionError=self%primaryProgenitorMassFunctionError+primaryProgenitorMassFunctionError
          else
             self%primaryProgenitorMassFunction     =-1.0d0
             self%primaryProgenitorMassFunctionError=-1.0d0
          end if
       end if
       ! Accumulate normalizations.
       self   %normalization                        =self%normalization                        +normalization
       self   %normalizationError                   =self%normalizationError                   +normalizationError
       if (self%computeCovariances) then
          self%normalizationCovariance              =self%normalizationCovariance              +normalizationCovariance
       end if
       if (self%extendedStatistics) then
          self%normalizationSubhaloMassFunction     =self%normalizationSubhaloMassFunction     +normalizationSubhaloMassFunction
          self%normalizationSubhaloMassFunctionError=self%normalizationSubhaloMassFunctionError+normalizationSubhaloMassFunctionError
       end if
       if (self%computeCovariances) then
          deallocate(normalizationCovariance           )
          deallocate(conditionalMassFunctionCovariance )
       end if
       if (self%extendedStatistics) then
          deallocate(normalizationSubhaloMassFunction     )
          deallocate(normalizationSubhaloMassFunctionError)
          deallocate(primaryProgenitorMassFunction        )
          deallocate(primaryProgenitorMassFunctionError   )
          deallocate(formationRateFunction                )
          deallocate(formationRateFunctionError           )
          deallocate(subhaloMassFunction                  )
          deallocate(subhaloMassFunctionError             )
       end if
    else
       ! Our group does not already exist. Simply write the data.
       conditionalMassFunctionGroup=outputFile%openGroup(char(self%nameGroupOutput),'Conditional mass functions of merger trees.',objectsOverwritable=.true.,overwriteOverride=.true.)
       call conditionalMassFunctionGroup%writeDataset  (self%massParents                       ,"massParent"                        ,"Mass of parent node [Msolar]"              ,datasetReturned=massDataset)
       call massDataset                 %writeAttribute(massSolar                              ,"unitsInSI"                                                                                                  )
       call conditionalMassFunctionGroup%writeDataset  (self%massRatios                        ,"massRatio"                         ,"Mass of ratio node [Msolar]"               ,datasetReturned=massDataset)
       call massDataset                 %writeAttribute(massSolar                              ,"unitsInSI"                                                                                                  )
       call conditionalMassFunctionGroup%writeDataset  (self%redshiftsParent                   ,"redshiftParent"                    ,"Redshift of parent node []"                                            )
       call conditionalMassFunctionGroup%writeDataset  (self%redshiftsProgenitor               ,"redshiftProgenitor"                ,"Redshift of progenitor node []"                                        )
       accumulationCount=0
    end if
    ! Increment the count of the number of mass functions which have been accumulated.
    accumulationCount=accumulationCount+1
    ! Normalize the conditional mass functions if this is the final thread to perform accumulation.
    !$ if (accumulationCount == OMP_Get_Num_Threads()) then
       if (self%computeCovariances) then
          ! If computing covariances, normalize the covariance first because it requires access to the unnormalized conditional mass
          ! function.
          do i=1,self%timeCount
             do j=1,self%countMassParent
                if (self%normalization(i,j) <= 0.0d0) cycle
                do k=1,self%massRatioCount
                   do j1=1,self%countMassParent
                      if (self%normalization(i,j1) <= 0.0d0) cycle
                      do k1=1,self%massRatioCount
                         self%conditionalMassFunctionCovariance(i,j,k,i,j1,k1)=+self%conditionalMassFunctionCovariance(i,j,k,i,j1,k1)    &
                              &                                                /self%normalization                    (i,j          )    &
                              &                                                /self%normalization                    (      i,j1   )    &
                              &                                                +self%normalizationCovariance          (i,j,  i,j1   )    &
                              &                                                *self%conditionalMassFunction          (i,j,k        )    &
                              &                                                *self%conditionalMassFunction          (      i,j1,k1)    &
                              &                                                /self%normalization                    (i,j          )**2 &
                              &                                                /self%normalization                    (      i,j1   )**2
                      end do
                   end do
                end do
             end do
          end do
       end if
       do i=1,self%timeCount
          do j=1,self%countMassParent
             if (self%normalization(i,j) <= 0.0d0) cycle
             where (self%conditionalMassFunction(i,j,:) > 0.0d0)
                self%conditionalMassFunctionError(i,j,:)=+sqrt(                                             &
                     &                                         +self%conditionalMassFunctionError(i,j,:)    &
                     &                                         +self%normalizationError          (i,j  )    &
                     &                                         *self%conditionalMassFunction     (i,j,:)**2 &
                     &                                         /self%normalization               (i,j  )**2 &
                     &                                        )                                             &
                     &                                   /      self%normalization               (i,j  )
             end where
             self%conditionalMassFunction(i,j,:  )=+self%conditionalMassFunction(i,j,:  ) &
                  &                                /self%normalization          (i,j    )
             if (self%extendedStatistics) then
                where (self%formationRateFunction(i,j,:,:) > 0.0d0)
                   self%formationRateFunctionError(i,j,:,:)=+      self%formationRateFunction     (i,j,:,:)    &
                        &                                   /      self%normalization             (i,j    )    &
                        &                                   *sqrt(                                             &
                        &                                         +self%formationRateFunctionError(i,j,:,:)    &
                        &                                         /self%formationRateFunction     (i,j,:,:)**2 &
                        &                                         +self%normalizationError        (i,j    )    &
                        &                                         /self%normalization             (i,j    )**2 &
                        &                                        )
                end where
                self%formationRateFunction(i,j,:,:)=self%formationRateFunction(i,j,:,:)/self%normalization(i,j)
                if (self%primaryProgenitorStatisticsValid) then
                   do iPrimary=1,self%depthProgenitorPrimary
                      where (self%primaryProgenitorMassFunction(i,j,:,:) > 0.0d0)
                         self%primaryProgenitorMassFunctionError(i,j,:,:)=+      self%primaryProgenitorMassFunction     (i,j,:,:)    &
                              &                                           /      self%normalization                     (i,j    )    &
                              &                                           *sqrt(                                                     &
                              &                                                 +self%primaryProgenitorMassFunctionError(i,j,:,:)    &
                              &                                                 /self%primaryProgenitorMassFunction     (i,j,:,:)**2 &
                              &                                                 +self%normalizationError                (i,j    )    &
                              &                                                 /self%normalization                     (i,j    )**2 &
                              &                                                )
                      end where
                      self%primaryProgenitorMassFunction(i,j,:,iPrimary)=self%primaryProgenitorMassFunction(i,j,:,iPrimary)/self%normalization(i,j)
                   end do
                end if
             end if
             ! Convert from squared error to error for normalization.
             self%normalizationError(i,j)=sqrt(self%normalizationError(i,j))
          end do
       end do
       ! Normalize subhalo mass functions.
       if (self%extendedStatistics) then
          do j=1,self%countMassParent
             if (self%normalizationSubhaloMassFunction(j) > 0.0d0) then
                where (self%subhaloMassFunction(j,:,:) > 0.0d0)
                   self%subhaloMassFunctionError(j,:,:)=+sqrt(                                                      &
                        &                                     +self%subhaloMassFunctionError             (j,:,:)    &
                        &                                     +self%normalizationSubhaloMassFunctionError(j    )    &
                        &                                     *self%subhaloMassFunction                  (j,:,:)**2 &
                        &                                     /self%normalizationSubhaloMassFunction     (j    )**2 &
                        &                                    )                                                      &
                        &                               /      self%normalizationSubhaloMassFunction     (j    )
                end where
                self%subhaloMassFunction     (j,:,:)=     self%subhaloMassFunction     (j,:,:) /self%normalizationSubhaloMassFunction(j)
                ! Convert from squared error to error for normalization.
                self%normalizationSubhaloMassFunctionError(j)=sqrt(self%normalizationSubhaloMassFunctionError(j))
            end if
          end do
       end if
    !$ end if
    call    conditionalMassFunctionGroup%writeAttribute(accumulationCount                       ,"accumulationCount"                                                                                   )
    call    conditionalMassFunctionGroup%writeDataset  (self%normalization                      ,"normalization"                          ,"Normalization for conditional mass functions []"           )
    call    conditionalMassFunctionGroup%writeDataset  (self%normalizationError                 ,"normalizationError"                     ,"Normalization error for conditional mass functions []"     )
    call    conditionalMassFunctionGroup%writeDataset  (self%conditionalMassFunction            ,"conditionalMassFunction"                ,"Conditional mass functions []"                             )
    call    conditionalMassFunctionGroup%writeDataset  (self%conditionalMassFunctionError       ,"conditionalMassFunctionError"           ,"Conditional mass function errors []"                       )
    if (self%computeCovariances) then
       call conditionalMassFunctionGroup%writeDataset  (self%normalizationCovariance            ,"normalizationCovariance"                ,"Normalization covariance for conditional mass functions []")
       call conditionalMassFunctionGroup%writeDataset  (self%conditionalMassFunctionCovariance  ,"conditionalMassFunctionCovariance"      ,"Conditional mass function covariance []"                   )
    end if
    if (self%extendedStatistics) then
       call conditionalMassFunctionGroup%writeDataset  (self%normalizationSubhaloMassFunction     ,"normalizationSubhaloMassFunction"     ,"Normalization for subhalo mass functions []"               )
       call conditionalMassFunctionGroup%writeDataset  (self%normalizationSubhaloMassFunctionError,"normalizationSubhaloMassFunctionError","Normalization for subhalo mass functions []"               )
       call conditionalMassFunctionGroup%writeDataset  (self%primaryProgenitorMassFunction        ,"primaryProgenitorMassFunction"        ,"Primary progenitor mass functions []"                      )
       call conditionalMassFunctionGroup%writeDataset  (self%primaryProgenitorMassFunctionError   ,"primaryProgenitorMassFunctionError"   ,"Primary progenitor mass function errors []"                )
       call conditionalMassFunctionGroup%writeDataset  (self%formationRateFunction                ,"formationRateFunction"                ,"Formation rate functions []"                               )
       call conditionalMassFunctionGroup%writeDataset  (self%formationRateFunctionError           ,"formationRateFunctionError"           ,"Formation rate function errors []"                         )
       call conditionalMassFunctionGroup%writeDataset  (self%subhaloMassFunction                  ,"subhaloMassFunction"                  ,"Unevolved subhalo mass functions []"                       )
       call conditionalMassFunctionGroup%writeDataset  (self%subhaloMassFunctionError             ,"subhaloMassFunctionError"             ,"Unevolved subhalo mass function errors []"                 )
    end if
    call    outputFile                  %flush         (                                                                                                                                               )
    !$ call hdf5Access%unset()
    return
  end subroutine conditionalMFFinalize
