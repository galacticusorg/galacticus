  !! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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
  
  !% Contains a module which implements a merger tree operator which accumulates conditional mass functions for trees.
  
  use Statistics_NBody_Halo_Mass_Errors

  !# <mergerTreeOperator name="mergerTreeOperatorConditionalMF" defaultThreadPrivate="yes">
  !#  <description>
  !#   Provides a merger tree operator which accumulates conditional mass functions for trees. In
  !#   addition to the cumulative mass function, 1$^{\rm st}$ through $n^{\rm th}$ most-massive
  !#   progenitor mass functions, formation rate functions, and unevolved subhalo mass functions \citep{jiang_generating_2014}
  !#   split by hierarchy depth are computed and output. Mass functions are accumulated in
  !#   logarithmically-spaced bins of parent halo mass, logarithmically-spaced bins of mass ratio
  !#   (the ratio of progenitor to parent halo mass), and at pairs of parent/progenitor
  !#   redshifts. The following parameters control the operator:
  !#   \begin{description}
  !#   \item[{\normalfont \ttfamily parentMassCount}] The number of bins in parent halo mass to use;
  !#   \item[{\normalfont \ttfamily parentMassMinimum}] The minimum parent halo mass to consider;
  !#   \item[{\normalfont \ttfamily parentMassMaximum}] The maximum parent halo mass to consider;
  !#   \item[{\normalfont \ttfamily massRatioCount}] The number of bins in mass ratio to use;
  !#   \item[{\normalfont \ttfamily massRatioMinimum}] The minimum mass ratio to consider;
  !#   \item[{\normalfont \ttfamily massRatioMaximum}] The maximum mass ratio to consider;
  !#   \item[{\normalfont \ttfamily parentRedshifts}] A list of redshifts at which to identify parent halos;
  !#   \item[{\normalfont \ttfamily progenitorRedshifts}] A corresponding list of redshifts at which to identify progenitor halos;
  !#   \item[{\normalfont \ttfamily primaryProgenitorDepth}] The number of $i^{\rm th}$ most-massive progenitor mass functions to compute (starting from the 1$^{\rm st}$ most-massive);
  !#   \item[{\normalfont \ttfamily subhaloHierarchyDepth}] The maximum depth in the subhalo hierarchy for which to compute the unevolved subhalo mass function;
  !#   \item[{\normalfont \ttfamily formationRateTimeFraction}] The fraction of the current time over which to estimate the formation rate of halos when computing merger tree statistics;
  !#   \item[{\normalfont \ttfamily outputGroupName}] The name of the \gls{hdf5} group to which mass functions will be written.
  !#   \end{description}  
  !#   If the operator finds the named \gls{hdf5} group already in existance, it will accumulate its
  !#   mass functions to those already written to the group, weighting by the inverse of the variance
  !#   in each bin. The structure of the \gls{hdf5} group is as follows:
  !#   \begin{verbatim}
  !#   {
  !#     DATASET "conditionalMassFunction" {
  !#     COMMENT "Conditional mass functions []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Nratio, Nparent, Nz ) }
  !#     }
  !#     DATASET "conditionalMassFunctionError" {
  !#     COMMENT "Conditional mass function errors []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Nratio, Nparent, Nz ) }
  !#     }
  !#     DATASET "formationRateFunction" {
  !#     COMMENT "Formation rate functions []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( 2, Nratio, Nparent, Nz ) }
  !#     }
  !#     DATASET "formationRateFunctionError" {
  !#     COMMENT "Formation rate function errors []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( 2, Nratio, Nparent, Nz ) }
  !#     }
  !#     DATASET "massParent" {
  !#     COMMENT "Mass of parent node [Msolar]"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Nparent ) }
  !#        ATTRIBUTE "unitsInSI" {
  !#           DATATYPE  H5T_IEEE_F64LE
  !#           DATASPACE  SCALAR
  !#           DATA {
  !#           (0): 1.98892e+30
  !#           }
  !#         }
  !#     }
  !#     DATASET "massRatio" {
  !#     COMMENT "Mass of ratio node [Msolar]"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Nratio ) }
  !#        ATTRIBUTE "unitsInSI" {
  !#           DATATYPE  H5T_IEEE_F64LE
  !#           DATASPACE  SCALAR
  !#           DATA {
  !#           (0): 1.98892e+30
  !#           }
  !#        }
  !#     }
  !#     DATASET "primaryProgenitorMassFunction" {
  !#     COMMENT "Primary progenitor mass functions []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Ndepth, Nratio, Nparent, Nz ) }
  !#     }
  !#     DATASET "primaryProgenitorMassFunctionError" {
  !#     COMMENT "Primary progenitor mass function errors []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Ndepth, Nratio, Nparent, Nz ) }
  !#     }
  !#     DATASET "redshiftParent" {
  !#     COMMENT "Redshift of parent node []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Nz ) }
  !#     }
  !#     DATASET "redshiftProgenitor" {
  !#     COMMENT "Redshift of progenitor node []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Nz ) }
  !#     }
  !#     DATASET "subhaloMassFunction" {
  !#     COMMENT "Unevolved subhalo mass functions []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Nhierarchy, Nratio, Nparent ) }
  !#     }
  !#     DATASET "subhaloMassFunctionError" {
  !#     COMMENT "Unevolved subhalo mass function errors []"
  !#        DATATYPE  H5T_IEEE_F64LE
  !#        DATASPACE  SIMPLE { ( Nhierarchy, Nratio, Nparent ) }
  !#     }
  !#   }
  !#   \end{verbatim}
  !#   Where {\normalfont \ttfamily Nratio} is the number of bins in mass ratio, {\normalfont \ttfamily
  !#   Nparent} is the number of bins in mass ratio, {\normalfont \ttfamily Nz} is the number of
  !#   parent/progenitor redshift pairs, {\normalfont \ttfamily Ndepth} is the maximum depth in the
  !#   ranking of most-massive progenitor mass functions, {\normalfont \ttfamily Nhierarchy} is the
  !#   maximum depth in the subhalo hierarchy for subhalo mass functions. The first dimension of the
  !#   {\normalfont \ttfamily formationRateFunction} dataset stores two different versions of the
  !#   formation rate function. The first uses the mass of the forming halo at the time of formation,
  !#   the second uses the mass of the node immediately prior to it becoming a subhalo.
  !#  
  !#   Mass functions are output as ${\rm d}N/{\rm d}\log_{10}m$ where $N$ is the number of halos per
  !#   parent halo, and $m$ is the mass ratio.  
  !#  </description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorConditionalMF
     !% A merger tree operator class which accumulates conditional mass functions for trees.
     private
     class           (nbodyHaloMassErrorClass), pointer                         :: haloMassError_
     double precision                         , allocatable, dimension(:      ) :: timeParents                         , timeProgenitors                    , &
          &                                                                        parentRedshifts                     , progenitorRedshifts                , &
          &                                                                        massParents                         , massRatios                         , &
          &                                                                        normalizationSubhaloMassFunction
     double precision                         , allocatable, dimension(:,:    ) :: normalization
     double precision                         , allocatable, dimension(:,:,:  ) :: conditionalMassFunction             , conditionalMassFunctionError
     double precision                         , allocatable, dimension(:,:,:  ) :: subhaloMassFunction                 , subhaloMassFunctionError
     double precision                         , allocatable, dimension(:,:,:,:) :: primaryProgenitorMassFunction       , primaryProgenitorMassFunctionError
     double precision                         , allocatable, dimension(:,:,:,:) :: formationRateFunction               , formationRateFunctionError
     integer                                                                    :: parentMassCount                     , primaryProgenitorDepth             , &
          &                                                                        massRatioCount                      , timeCount                          , &
          &                                                                        subhaloHierarchyDepth
     double precision                                                           :: massParentLogarithmicMinimum        , massRatioLogarithmicMinimum        , &
          &                                                                        massParentLogarithmicBinWidthInverse, massRatioLogarithmicBinWidthInverse, &
          &                                                                        formationRateTimeFraction
     logical                                                                    :: alwaysIsolatedHalosOnly             , primaryProgenitorStatisticsValid   , &
          &                                                                        extendedStatistics
     type            (varying_string         )                                  :: outputGroupName
   contains
     !@ <objectMethods>
     !@   <object>mergerTreeOperatorConditionalMF</object>
     !@   <objectMethod>
     !@     <method>binWeights</method>
     !@     <type>\doubleone</type>
     !@     <arguments>\doublezero\ mass\argin,\doublezero\ time\argin,\doublezero\ massLogarithmicMinimumBins\argin,\doublezero\ massLogarithmicWidthInverseBins\argin,\intzero\ countBins\argin</arguments>
     !@     <description>Compute weights for a halo in each bin of the mass function.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>binWeights2D</method>
     !@     <type>\doubletwo</type>
     !@     <arguments>\doublezero\ mass1\argin,\doublezero\ time1\argin,\doublezero\ mass2\argin,\doublezero\ time2\argin,\doublezero\ massLogarithmicMinimumBins1\argin,\doublezero\ massLogarithmicWidthInverseBins1\argin,\intzero\ countBins2\argin,\doublezero\ massLogarithmicMinimumBins2\argin,\doublezero\ massLogarithmicWidthInverseBins2\argin,\intzero\ countBins2\argin</arguments>
     !@     <description>Compute weights for a halo in each bin of a 2D mass function.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                 conditionalMFDestructor
     procedure :: operate      => conditionalMFOperate
     procedure :: finalize     => conditionalMFFinalize
     procedure :: binWeights   => conditionalMFBinWeights
     procedure :: binWeights2D => conditionalMFBinWeights2D
  end type mergerTreeOperatorConditionalMF
  
  interface mergerTreeOperatorConditionalMF
     !% Constructors for the conditional mass function merger tree operator class.
     module procedure conditionalMFConstructorParameters
     module procedure conditionalMFConstructorInternal
  end interface mergerTreeOperatorConditionalMF

contains
  
  function conditionalMFConstructorParameters(parameters)
    !% Constructor for the conditional mass function merger tree operator class which takes a parameter set as input.
    use Memory_Management
    implicit none
    type            (mergerTreeOperatorConditionalMF)                              :: conditionalMFConstructorParameters
    type            (inputParameters                ), intent(inout), target       :: parameters
    double precision                                 , allocatable  , dimension(:) :: progenitorRedshifts               , parentRedshifts
    class           (nbodyHaloMassErrorClass        ), pointer                     :: haloMassError_
    integer                                                                        :: parentMassCount                   , massRatioCount       , &
         &                                                                            primaryProgenitorDepth            , subhaloHierarchyDepth
    double precision                                                               :: parentMassMinimum                 , parentMassMaximum    , &
         &                                                                            massRatioMinimum                  , massRatioMaximum     , &
         &                                                                            formationRateTimeFraction
    logical                                                                           alwaysIsolatedHalosOnly           , extendedStatistics
    type            (varying_string                 )                              :: outputGroupName
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)
    !# <objectBuilder class="nbodyHaloMassError" name="haloMassError_" source="parameters"/>
    !# <inputParameter>
    !#   <name>parentMassCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins in parent mass when constructing conditional halo mass functions.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>parentMassMinimum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d10</defaultValue>
    !#   <description>The minimum parent halo mass to bin when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>parentMassMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d15</defaultValue>
    !#   <description>The maximum parent halo mass to bin when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massRatioCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins in mass ratio when constructing conditional halo mass functions.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massRatioMinimum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d-4</defaultValue>
    !#   <description>The minimum mass ratio to bin when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massRatioMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d1</defaultValue>
    !#   <description>The maximum mass ratio to bin when constructing conditional halo mass functions.</description>
    !#   <type>logical</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    call allocateArray(parentRedshifts    ,[max(1,parameters%count('parentRedshifts'    ,zeroIfNotPresent=.true.))])
    !# <inputParameter>
    !#   <name>parentRedshifts</name>
    !#   <source>parameters</source>
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The set of parent halo redshifts to use when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1..</cardinality>
    !# </inputParameter>
    call allocateArray(progenitorRedshifts,[max(1,parameters%count('progenitorRedshifts',zeroIfNotPresent=.true.))])
    !# <inputParameter>
    !#   <name>progenitorRedshifts</name>
    !#   <source>parameters</source>
    !#   <defaultValue>[1.0d0]</defaultValue>
    !#   <description>The set of progenitor halo redshifts to use when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>primaryProgenitorDepth</name>
    !#   <source>parameters</source>
    !#   <defaultValue>2</defaultValue>
    !#   <description>The depth in progenitor ranking for which to store ranked progenitor mass functions. For example, a value of 2 means store mass functions for the most massive, and second most massive progenitor.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>subhaloHierarchyDepth</name>
    !#   <source>parameters</source>
    !#   <defaultValue>2</defaultValue>
    !#   <description>The depth in the subhalo hierarchy for which to store unevolved subhalo mass functions.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>formationRateTimeFraction</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <description>The fraction of the current time over which to estimate the formation rate of halos when computing merger tree statistics.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alwaysIsolatedHalosOnly</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Include only halos which have always been isolated when computing merger tree statistics?</description>
    !#   <type>boolean</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>extendedStatistics</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>Compute extended statistics (formation rate function, $N^{\rm th}$ most-massive progenitor mass functions, etc.)?</description>
    !#   <type>boolean</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputGroupName</name>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('conditionalMassFunction')</defaultValue>
    !#   <description>The name of the HDF5 group to which the conditional mass function should be output.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    ! Construct the instance.
    conditionalMFConstructorParameters=conditionalMFConstructorInternal(                           &
         &                                                              parentMassCount          , &
         &                                                              parentMassMinimum        , &
         &                                                              parentMassMaximum        , &
         &                                                              massRatioCount           , &
         &                                                              massRatioMinimum         , &
         &                                                              massRatioMaximum         , &
         &                                                              parentRedshifts          , &
         &                                                              progenitorRedshifts      , &
         &                                                              primaryProgenitorDepth   , &
         &                                                              formationRateTimeFraction, &
         &                                                              subhaloHierarchyDepth    , &
         &                                                              alwaysIsolatedHalosOnly  , &
         &                                                              extendedStatistics       , &
         &                                                              outputGroupName          , &
         &                                                              haloMassError_             &
         &                                                             )
    return
  end function conditionalMFConstructorParameters

  function conditionalMFConstructorInternal(parentMassCount,parentMassMinimum,parentMassMaximum,massRatioCount,massRatioMinimum,massRatioMaximum,parentRedshifts,progenitorRedshifts,primaryProgenitorDepth,formationRateTimeFraction,subhaloHierarchyDepth,alwaysIsolatedHalosOnly,extendedStatistics,outputGroupName,haloMassError_)
    !% Internal constructor for the conditional mass function merger tree operator class.
    use Cosmology_Functions
    use Numerical_Ranges
    use Memory_Management
    use Galacticus_Error
    implicit none
    type            (mergerTreeOperatorConditionalMF)                              :: conditionalMFConstructorInternal
    double precision                                 , intent(in   ), dimension(:) :: parentRedshifts                 , progenitorRedshifts
    integer                                          , intent(in   )               :: parentMassCount                 , massRatioCount       , &
         &                                                                            primaryProgenitorDepth          , subhaloHierarchyDepth
    double precision                                 , intent(in   )               :: parentMassMinimum               , parentMassMaximum    , &
         &                                                                            massRatioMinimum                , massRatioMaximum     , &
         &                                                                            formationRateTimeFraction
    logical                                          , intent(in   )               :: alwaysIsolatedHalosOnly         , extendedStatistics
    type            (varying_string                 ), intent(in   )               :: outputGroupName
    class           (nbodyHaloMassErrorClass        ), intent(in   ), target       :: haloMassError_
    class           (cosmologyFunctionsClass        ), pointer                     :: cosmologyFunctions_
    integer                                                                        :: i

    ! Store array sizes.
    conditionalMFConstructorInternal%timeCount                         =  size(parentRedshifts)
    conditionalMFConstructorInternal%parentMassCount                   =  parentMassCount
    conditionalMFConstructorInternal%massRatioCount                    =  massRatioCount
    conditionalMFConstructorInternal%primaryProgenitorDepth            =  primaryProgenitorDepth
    conditionalMFConstructorInternal%subhaloHierarchyDepth             =  subhaloHierarchyDepth
    conditionalMFConstructorInternal%formationRateTimeFraction         =  formationRateTimeFraction
    conditionalMFConstructorInternal%alwaysIsolatedHalosOnly           =  alwaysIsolatedHalosOnly
    conditionalMFConstructorInternal%extendedStatistics                =  extendedStatistics
    conditionalMFConstructorInternal%outputGroupName                   =  outputGroupName
    conditionalMFConstructorInternal%haloMassError_                    => haloMassError_
    conditionalMFConstructorInternal%primaryProgenitorStatisticsValid  =  conditionalMFConstructorInternal%haloMassError_%errorZeroAlways()
    if (size(progenitorRedshifts) /= conditionalMFConstructorInternal%timeCount) &
         & call Galacticus_Error_Report('conditionalMFConstructorInternal','mismatch in sizes of parent and progenitor redshift arrays')
    ! Check for required property attributes.
    if (alwaysIsolatedHalosOnly.and..not.defaultMergingStatisticsComponent%nodeHierarchyLevelMaximumIsGettable())                                              &
         & call Galacticus_Error_Report                                                                                                                        &
         &      (                                                                                                                                              &
         &       'conditionalMFConstructorInternal'                                                                                                         ,  &
         &       'statistics of always isolated halos require a merging statistics component that provides a gettable "nodeHierarchyLevelMaximum" property.'// &
         &       Galacticus_Component_List(                                                                                                                    &
         &                                 'mergingStatistics'                                                                                               , &
         &                                  defaultMergingStatisticsComponent%nodeHierarchyLevelMaximumAttributeMatch(requireGettable=.true.)                  &
         &                                )                                                                                                                    &
         &      )     
    ! Allocate arrays.
    call allocateArray(conditionalMFConstructorInternal%parentRedshifts    ,[conditionalMFConstructorInternal%timeCount        ])
    call allocateArray(conditionalMFConstructorInternal%progenitorRedshifts,[conditionalMFConstructorInternal%timeCount        ])
    call allocateArray(conditionalMFConstructorInternal%timeProgenitors    ,[conditionalMFConstructorInternal%timeCount        ])
    call allocateArray(conditionalMFConstructorInternal%timeParents        ,[conditionalMFConstructorInternal%timeCount        ])
    call allocateArray(conditionalMFConstructorInternal%massParents        ,[conditionalMFConstructorInternal%parentMassCount+1])
    call allocateArray(conditionalMFConstructorInternal%massRatios         ,[conditionalMFConstructorInternal%massRatioCount +1])
    call allocateArray(                                                                      &
         &            conditionalMFConstructorInternal%normalization                     , &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                     &
         &           ]                                                                     &
         &          )
    call allocateArray(                                                                      &
         &            conditionalMFConstructorInternal%conditionalMassFunction           , &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                   , &
         &            conditionalMFConstructorInternal%massRatioCount                      &
         &           ]                                                                     &
         &          )
    call allocateArray(                                                                      &
         &            conditionalMFConstructorInternal%conditionalMassFunctionError      , &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                   , &
         &            conditionalMFConstructorInternal%massRatioCount                      &
         &           ]                                                                     &
         &          )
    if (conditionalMFConstructorInternal%extendedStatistics) then
       call allocateArray(                                                                      &
            &            conditionalMFConstructorInternal%normalizationSubhaloMassFunction  , &
            &           [                                                                     &
            &            conditionalMFConstructorInternal%parentMassCount                     &
            &           ]                                                                     &
            &          )
       call allocateArray(                                                                      &
            &            conditionalMFConstructorInternal%primaryProgenitorMassFunction     , &
            &           [                                                                     &
            &            conditionalMFConstructorInternal%timeCount                         , &
            &            conditionalMFConstructorInternal%parentMassCount                   , &
            &            conditionalMFConstructorInternal%massRatioCount                    , &
            &            conditionalMFConstructorInternal%primaryProgenitorDepth              &
            &           ]                                                                     &
            &          )
       call allocateArray(                                                                      &
            &            conditionalMFConstructorInternal%primaryProgenitorMassFunctionError, &
            &           [                                                                     &
            &            conditionalMFConstructorInternal%timeCount                         , &
            &            conditionalMFConstructorInternal%parentMassCount                   , &
            &            conditionalMFConstructorInternal%massRatioCount                    , &
            &            conditionalMFConstructorInternal%primaryProgenitorDepth              &
            &           ]                                                                     &
            &          )
       call allocateArray(                                                                      &
            &            conditionalMFConstructorInternal%formationRateFunction             , &
            &           [                                                                     &
            &            conditionalMFConstructorInternal%timeCount                         , &
            &            conditionalMFConstructorInternal%parentMassCount                   , &
            &            conditionalMFConstructorInternal%massRatioCount                    , &
            &            2                                                                    &
            &           ]                                                                     &
            &          )
       call allocateArray(                                                                      &
            &            conditionalMFConstructorInternal%formationRateFunctionError        , &
            &           [                                                                     &
            &            conditionalMFConstructorInternal%timeCount                         , &
            &            conditionalMFConstructorInternal%parentMassCount                   , &
            &            conditionalMFConstructorInternal%massRatioCount                    , &
            &            2                                                                    &
            &           ]                                                                     &
            &          )
       call allocateArray(                                                                      &
            &            conditionalMFConstructorInternal%subhaloMassFunction               , &
            &           [                                                                     &
            &            conditionalMFConstructorInternal%parentMassCount                   , &
            &            conditionalMFConstructorInternal%massRatioCount                    , &
            &            conditionalMFConstructorInternal%subhaloHierarchyDepth               &
            &           ]                                                                     &
            &          )
       call allocateArray(                                                                      &
            &            conditionalMFConstructorInternal%subhaloMassFunctionError          , &
            &           [                                                                     &
            &            conditionalMFConstructorInternal%parentMassCount                   , &
            &            conditionalMFConstructorInternal%massRatioCount                    , &
            &            conditionalMFConstructorInternal%subhaloHierarchyDepth               &
            &           ]                                                                     &
            &          )
    end if
    ! Construct bins for parent node mass.
    conditionalMFConstructorInternal%massParentLogarithmicMinimum        =  log( parentMassMinimum)
    conditionalMFConstructorInternal%massParentLogarithmicBinWidthInverse= dble( parentMassCount  ) &
         &                                                                / log(                    &
         &                                                                      +parentMassMaximum  &
         &                                                                      /parentMassMinimum  &
         &                                                                     )
    conditionalMFConstructorInternal%massParents=Make_Range(                                        &
         &                                                  parentMassMinimum                     , &
         &                                                  parentMassMaximum                     , &
         &                                                  parentMassCount                       , &
         &                                                  rangeType        =rangeTypeLogarithmic  &
         &                                                 )
    ! Construct bins for mass ratio.
    conditionalMFConstructorInternal%massRatioLogarithmicMinimum        =  log( massRatioMinimum)
    conditionalMFConstructorInternal%massRatioLogarithmicBinWidthInverse= dble( massRatioCount  ) &
         &                                                               / log(                   &
         &                                                                     +massRatioMaximum  &
         &                                                                     /massRatioMinimum  &
         &                                                                    )
    conditionalMFConstructorInternal%massRatios=Make_Range(                                       &
         &                                                 massRatioMinimum                     , &
         &                                                 massRatioMaximum                     , &
         &                                                 massRatioCount                       , &
         &                                                 rangeType       =rangeTypeLogarithmic, &
         &                                                 rangeBinned     =.true.                &
         &                                                )
    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()
    ! Construct arrays of times for progenitors.
    conditionalMFConstructorInternal%progenitorRedshifts=progenitorRedshifts
    do i=1,conditionalMFConstructorInternal%timeCount
       conditionalMFConstructorInternal%timeProgenitors(i)=              &
            & cosmologyFunctions_%cosmicTime(                   &
            &  cosmologyFunctions_%expansionFactorFromRedshift( &
            &   progenitorRedshifts(i)                          &
            &  )                                                &
            & )
    end do
    ! Construct arrays of times for parents.
    conditionalMFConstructorInternal%parentRedshifts=parentRedshifts
    do i=1,conditionalMFConstructorInternal%timeCount
       conditionalMFConstructorInternal%timeParents(i)=         & 
            & cosmologyFunctions_%cosmicTime(                   &
            &  cosmologyFunctions_%expansionFactorFromRedshift( &
            &   parentRedshifts(i)                              &
            &  )                                                &
            & )
    end do
    ! Initialize mass function arrays.
    conditionalMFConstructorInternal   %normalization                     =0.0d0
    conditionalMFConstructorInternal   %conditionalMassFunction           =0.0d0
    conditionalMFConstructorInternal   %conditionalMassFunctionError      =0.0d0
    if (conditionalMFConstructorInternal%extendedStatistics) then
       conditionalMFConstructorInternal%normalizationSubhaloMassFunction  =0.0d0
       conditionalMFConstructorInternal%primaryProgenitorMassFunction     =0.0d0
       conditionalMFConstructorInternal%primaryProgenitorMassFunctionError=0.0d0
       conditionalMFConstructorInternal%formationRateFunction             =0.0d0
       conditionalMFConstructorInternal%formationRateFunctionError        =0.0d0
       conditionalMFConstructorInternal%subhaloMassFunction               =0.0d0
       conditionalMFConstructorInternal%subhaloMassFunctionError          =0.0d0
    end if
    return
  end function conditionalMFConstructorInternal

  subroutine conditionalMFDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorConditionalMF), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    !# <objectDestructor name="self%haloMassError_"/>
    return
  end subroutine conditionalMFDestructor

  subroutine conditionalMFOperate(self,tree)
    !% Compute conditional mass function on {\normalfont \ttfamily tree}.
    use Galacticus_Nodes
    use Input_Parameters
    use Memory_Management
    use Numerical_Comparison
    use Galacticus_Error
    implicit none
    class           (mergerTreeOperatorConditionalMF), intent(inout)                                  :: self
    type            (mergerTree                     ), intent(inout)                         , target :: tree
    type            (treeNode                       ), pointer                                        :: node                          , nodeChild                 , &
         &                                                                                               nodeParent                    , nodeParentChild           , &
         &                                                                                               descendentNode                , nodeSibling
    type            (mergerTree                     ), pointer                                        :: treeCurrent
    class           (nodeComponentBasic             ), pointer                                        :: basic                         , basicChild                , &
         &                                                                                               basicParent                   , descendentBasic           , &
         &                                                                                               basicParentChild              , basicSibling
    class           (nodeComponentMergingStatistics ), pointer                                        :: mergingStatistics
    integer                                                                                           :: i                             , binMassParent             , &
         &                                                                                               binMassRatio                  , iPrimary                  , &
         &                                                                                               jPrimary                      , depthHierarchy
    double precision                                                                                  :: branchBegin                   , branchEnd                 , &
         &                                                                                               parentBranchBegin             , parentBranchEnd           , &
         &                                                                                               massProgenitor                , massParent                , &
         &                                                                                               branchMassInitial             , branchMassFinal           , &
         &                                                                                               parentBranchMassInitial       , parentBranchMassFinal     , &
         &                                                                                               massRatioLogarithmic          , timeUnevolved             , &
         &                                                                                               massRatio                     , massUnevolved
    double precision                                  , dimension(                                                                                                   &
         &                                                        self%timeCount             ,                                                                       &
         &                                                        self%parentMassCount       ,                                                                       &
         &                                                        self%primaryProgenitorDepth                                                                        &
         &                                                       )                                    :: primaryProgenitorMass
    double precision                                  , dimension(self%parentMassCount) :: weights1D
    double precision                                  , dimension(self%parentMassCount,self%massRatioCount) :: weights2D, weights2DError
    
    ! Iterate over trees.
    treeCurrent => tree    
    do while (associated(treeCurrent))
       ! Initialize primary progenitor masses to zero.
       primaryProgenitorMass=0.0d0
       ! Get root node of the tree.       
       node => treeCurrent%baseNode
       ! Accumulate normalization for subhalo mass function.
       if (self%extendedStatistics) then
          basic => node%basic()
          weights1D=self%binWeights(                                           &
               &                    basic%mass()                             , &
               &                    basic%time()                             , &
               &                    self%massParentLogarithmicMinimum        , &
               &                    self%massParentLogarithmicBinWidthInverse, &
               &                    self%parentMassCount                       &
               &                   )
          !$omp critical(conditionalMassFunctionAccumulate)
          self%normalizationSubhaloMassFunction(:)=+self       %normalizationSubhaloMassFunction(:) &
               &                                   +weights1D                                       &
               &                                   *treeCurrent%volumeWeight
          !$omp end critical(conditionalMassFunctionAccumulate)
       end if
       ! Walk the tree, accumulating statistics.
       do while (associated(node))
          ! Get the child node, and process if child exists.
          nodeChild => node%firstChild
          do while (associated(nodeChild))
             ! Check if child should be included.
             mergingStatistics => nodeChild%mergingStatistics()
             if     (                                                         &
                  &   .not.self             %alwaysIsolatedHalosOnly          &
                  &  .or.                                                     &
                  &        mergingStatistics%nodeHierarchyLevelMaximum() == 0 &
                  & ) then
                ! Get the basic components.
                basic      => node     %basic()
                basicChild => nodeChild%basic()
                ! Determine range of times spanned by this branch.
                branchBegin=basicChild%time()
                branchEnd  =basic     %time()
                ! Does the branch span a progenitor node time?
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
                      weights1D=self%binWeights(                                              &
                           &                    massParent                                  , &
                           &                    self%timeParents                         (i), &
                           &                    self%massParentLogarithmicMinimum           , &
                           &                    self%massParentLogarithmicBinWidthInverse   , &
                           &                    self%parentMassCount                          &
                           &                   )
                      !$omp critical(conditionalMassFunctionAccumulate)
                      self%normalization(i,:)=+self       %normalization(i,:) &
                           &                  +weights1D                      &
                           &                  *treeCurrent%volumeWeight
                      !$omp end critical(conditionalMassFunctionAccumulate)
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
                                 &                           self%parentMassCount                        , &
                                 &                           self%massRatioLogarithmicMinimum            , &
                                 &                           self%massRatioLogarithmicBinWidthInverse    , &
                                 &                           self%massRatioCount                         , &
                                 &                           1                                             &
                                 &                          )
                            weights2DError=self%binWeights2D(                                              &
                                 &                           massParent                                  , &
                                 &                           self%timeParents                         (i), &
                                 &                           massProgenitor                              , &
                                 &                           self%timeProgenitors                     (i), &
                                 &                           self%massParentLogarithmicMinimum           , &
                                 &                           self%massParentLogarithmicBinWidthInverse   , &
                                 &                           self%parentMassCount                        , &
                                 &                           self%massRatioLogarithmicMinimum            , &
                                 &                           self%massRatioLogarithmicBinWidthInverse    , &
                                 &                           self%massRatioCount                         , &
                                 &                           2                                             &
                                 &                          )
                            !$omp critical(conditionalMassFunctionAccumulate)
                            self%conditionalMassFunction     (i,:,:)=+self          %conditionalMassFunction     (i,:,:) &
                                 &                                   +weights2D                                          &
                                 &                                   *treeCurrent   %volumeWeight                                  
                            self%conditionalMassFunctionError(i,:,:)=+self          %conditionalMassFunctionError(i,:,:) &
                                 &                                   +weights2DError                                     &
                                 &                                   *treeCurrent   %volumeWeight**2
                            !$omp end critical(conditionalMassFunctionAccumulate)
                            ! Check for formation.
                            if (self%extendedStatistics) then
                               if (branchBegin > self%timeProgenitors(i)*(1.0d0-self%formationRateTimeFraction) .and. .not.associated(nodeChild%firstChild)) then
                                  ! This is a newly formed halo, accumulate to formation rate arrays.
                                  weights2D     =self%binWeights2D(                                              &
                                       &                           massParent                                  , &
                                       &                           self%timeParents                         (i), &
                                       &                           branchMassInitial                           , &
                                       &                           basicChild%time                          ( ), &
                                       &                           self%massParentLogarithmicMinimum           , &
                                       &                           self%massParentLogarithmicBinWidthInverse   , &
                                       &                           self%parentMassCount                        , &
                                       &                           self%massRatioLogarithmicMinimum            , &
                                       &                           self%massRatioLogarithmicBinWidthInverse    , &
                                       &                           self%massRatioCount                         , &
                                       &                           1                                             &
                                       &                          )
                                  weights2DError=self%binWeights2D(                                              &
                                       &                           massParent                                  , &
                                       &                           self%timeParents                         (i), &
                                       &                           branchMassInitial                           , &
                                       &                           basicChild%time                          ( ), &
                                       &                           self%massParentLogarithmicMinimum           , &
                                       &                           self%massParentLogarithmicBinWidthInverse   , &
                                       &                           self%parentMassCount                        , &
                                       &                           self%massRatioLogarithmicMinimum            , &
                                       &                           self%massRatioLogarithmicBinWidthInverse    , &
                                       &                           self%massRatioCount                         , &
                                       &                           2                                             &
                                       &                          )
                                  !$omp critical(conditionalMassFunctionAccumulate)
                                  self%formationRateFunction     (i,:,:,1)=+self          %formationRateFunction     (i,:,:,1) &
                                       &                                   +weights2D                                          &
                                       &                                   *treeCurrent   %volumeWeight               
                                  self%formationRateFunctionError(i,:,:,1)=+self          %formationRateFunctionError(i,:,:,1) &
                                       &                                   +weights2DError                                     &
                                       &                                   *treeCurrent   %volumeWeight**2                                  
                                  !$omp end critical(conditionalMassFunctionAccumulate)
                                  ! Find the mass of this node just prior to it becoming a subhalo.
                                  descendentNode => node
                                  do while (associated(descendentNode%parent).and.associated(descendentNode%parent%firstChild,descendentNode))
                                     descendentNode => descendentNode%parent
                                  end do
                                  descendentBasic   => descendentNode%basic()
                                  weights2D     =self%binWeights2D(                                              &
                                       &                           massParent                                  , &
                                       &                           self%timeParents                         (i), &
                                       &                           descendentBasic%mass                     ( ), &
                                       &                           descendentBasic%time                     ( ), &
                                       &                           self%massParentLogarithmicMinimum           , &
                                       &                           self%massParentLogarithmicBinWidthInverse   , &
                                       &                           self%parentMassCount                        , &
                                       &                           self%massRatioLogarithmicMinimum            , &
                                       &                           self%massRatioLogarithmicBinWidthInverse    , &
                                       &                           self%massRatioCount                         , &
                                       &                           1                                             &
                                       &                          ) 
                                  weights2DError=self%binWeights2D(                                              &
                                       &                           massParent                                  , &
                                       &                           self%timeParents                         (i), &
                                       &                           descendentBasic%mass                     ( ), &
                                       &                           descendentBasic%time                     ( ), &
                                       &                           self%massParentLogarithmicMinimum           , &
                                       &                           self%massParentLogarithmicBinWidthInverse   , &
                                       &                           self%parentMassCount                        , &
                                       &                           self%massRatioLogarithmicMinimum            , &
                                       &                           self%massRatioLogarithmicBinWidthInverse    , &
                                       &                           self%massRatioCount                         , &
                                       &                           2                                             &
                                       &                          )
                                  !$omp critical(conditionalMassFunctionAccumulate)
                                  self%formationRateFunction     (i,:,:,2)=+self          %formationRateFunction     (i,:,:,2) &
                                       &                                   +weights2D                                          &
                                       &                                   *treeCurrent   %volumeWeight               
                                  self%formationRateFunctionError(i,:,:,2)=+self          %formationRateFunctionError(i,:,:,2) &
                                       &                                   +weights2DError                                     &
                                       &                                   *treeCurrent   %volumeWeight**2                                  
                                  !$omp end critical(conditionalMassFunctionAccumulate)
                               end if
                               ! Accumulate to the primary progenitor mass array if necessary.
                               if (self%primaryProgenitorStatisticsValid) then
                                  binMassParent=int(                                                      &
                                       &            +(+log(massParent)-self%massParentLogarithmicMinimum) &
                                       &            *self%massParentLogarithmicBinWidthInverse            &
                                       &           )                                                      &
                                       &        +1
                                  if (binMassParent >= 1 .and. binMassParent <= self%parentMassCount) then
                                     massRatio=massProgenitor/massParent
                                     iPrimary =1
                                     do while (massRatio < primaryProgenitorMass(i,binMassParent,iPrimary))
                                        iPrimary=iPrimary+1
                                        if (iPrimary > self%primaryProgenitorDepth) exit
                                     end do
                                     if (iPrimary <= self%primaryProgenitorDepth) then
                                        if (iPrimary < self%primaryProgenitorDepth) then
                                           do jPrimary=self%primaryProgenitorDepth,iPrimary+1,-1
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
                         nodeParent => nodeParent%parent
                      end do parentWalk
                   end if
                   ! Record the mass of the branch at the parent time.
                end do
                ! Accumulate unevoled subhalo mass functions.
                if (self%extendedStatistics) then
                   if (.not.associated(nodeChild%firstChild)) then
                      ! This is a branch tip. Follow until it to the final time, storing its mass just prior to becoming a subhalo, and
                      ! the hierarchy depth.                
                      depthHierarchy =  0
                      descendentNode => nodeChild
                      do while (associated(descendentNode).and.depthHierarchy <= self%subhaloHierarchyDepth)
                         if (associated(descendentNode%parent).and..not.descendentNode%isPrimaryProgenitor()) then
                            depthHierarchy=depthHierarchy+1
                            if (depthHierarchy == 1) then
                               descendentBasic => descendentNode %basic()
                               massUnevolved   =  descendentBasic%mass ()
                               timeUnevolved   =  descendentBasic%time ()
                            end if
                         end if
                         descendentNode => descendentNode%parent
                      end do
                      if (depthHierarchy > 0 .and. depthHierarchy <= self%subhaloHierarchyDepth) then
                         basicParent => treeCurrent%baseNode%basic()
                         weights2D     =self%binWeights2D(                                           &
                              &                           basicParent%mass()                       , &
                              &                           basicParent%time()                       , &
                              &                           massUnevolved                            , &
                              &                           timeUnevolved                            , &
                              &                           self%massParentLogarithmicMinimum        , &
                              &                           self%massParentLogarithmicBinWidthInverse, &
                              &                           self%parentMassCount                     , &
                              &                           self%massRatioLogarithmicMinimum         , &
                              &                           self%massRatioLogarithmicBinWidthInverse , &
                              &                           self%massRatioCount                      , &
                              &                           1                                          &
                              &                          )
                         weights2DError=self%binWeights2D(                                           &
                              &                           basicParent%mass()                       , &
                              &                           basicParent%time()                       , &
                              &                           massUnevolved                            , &
                              &                           timeUnevolved                            , &
                              &                           self%massParentLogarithmicMinimum        , &
                              &                           self%massParentLogarithmicBinWidthInverse, &
                              &                           self%parentMassCount                     , &
                              &                           self%massRatioLogarithmicMinimum         , &
                              &                           self%massRatioLogarithmicBinWidthInverse , &
                              &                           self%massRatioCount                      , &
                              &                           2                                          &
                              &                          )
                         !$omp critical(conditionalMassFunctionAccumulate)
                         self%subhaloMassFunction     (:,:,depthHierarchy)=+self          %subhaloMassFunction     (:,:,depthHierarchy) &
                              &                                            +weights2D                                                   &
                              &                                            *treeCurrent   %volumeWeight               
                         self%subhaloMassFunctionError(:,:,depthHierarchy)=+self          %subhaloMassFunctionError(:,:,depthHierarchy) &
                              &                                            +weights2DError                                              &
                              &                                            *treeCurrent   %volumeWeight**2                                  
                         !$omp end critical(conditionalMassFunctionAccumulate)
                      end if
                   end if
                end if
             end if
             ! Move to the next child.
             nodeChild => nodeChild%sibling
          end do
          ! Move to the next node.
          node => node%walkTree()
       end do
       ! Store the computed primary progenitor mass functions.
       if (self%extendedStatistics.and.self%primaryProgenitorStatisticsValid) then
          do i=1,self%timeCount
             do binMassParent=1,self%parentMassCount
                do iPrimary=1,self%primaryProgenitorDepth
                   if (primaryProgenitorMass(i,binMassParent,iPrimary) > 0.0d0) then
                      massRatioLogarithmic=log(primaryProgenitorMass(i,binMassParent,iPrimary))
                      binMassRatio         =int(                                                           &
                           &                     (massRatioLogarithmic-self%massRatioLogarithmicMinimum)   &
                           &                    *                                                          &
                           &                     self%massRatioLogarithmicBinWidthInverse                  &
                           &                   )                                                           &
                           &                +1
                      if (binMassRatio  >= 1 .and. binMassRatio  <= self%massRatioCount) then
                         !$omp critical(conditionalMassFunctionAccumulate)
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
                         !$omp end critical(conditionalMassFunctionAccumulate)
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
  end subroutine conditionalMFOperate

  function conditionalMFBinWeights(self,mass,time,massLogarithmicMinimumBins,massLogarithmicWidthInverseBins,countBins)
    !% Computes the weight that a given halo contributes to an array of bins.
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
    !% Computes the weight that a given halo contributes to a 2D array of bins.
    use FGSL
    use Numerical_Integration
    use Galacticus_Error
    implicit none
    class           (mergerTreeOperatorConditionalMF), intent(inout)                    :: self
    double precision                                 , intent(in   )                    :: mass1                           , time1                                , &
         &                                                                                 massLogarithmicMinimumBins1     , massLogarithmicWidthInverseBins1     , &
         &                                                                                 mass2                           , time2                                , &
         &                                                                                 massRatioLogarithmicMinimumBins2, massRatioLogarithmicWidthInverseBins2
    integer                                          , intent(in   )                    :: countBins1                      , countBins2                           , &
         &                                                                                 moment
    double precision                                 , dimension(countBins1,countBins2) :: conditionalMFBinWeights2D
    type            (treeNode                       ), pointer                          :: node1                           , node2
    class           (nodeComponentBasic             ), pointer                          :: basic1                          , basic2
    double precision                                 , parameter                        :: integrationExtent               =10.0d0
    double precision                                                                    :: massError1                      , massError2                           , &
         &                                                                                 mass1LowerLimit                 , mass1UpperLimit                      , &
         &                                                                                 mass2LowerLimit                 , mass2UpperLimit
    integer                                                                             :: i                               , j
    type            (fgsl_function                  )                                   :: integrandFunction
    type            (fgsl_integration_workspace     )                                   :: integrationWorkspace

    ! Validate moment.
    if (moment < 0 .or. moment > 2) call Galacticus_Error_Report('conditionalMFBinWeights2D','moment must be 0, 1, or 2')
    ! Construct nodes and find the mass errors.
    node1  => treeNode       (                 )
    node2  => treeNode       (                 )
    basic1 => node1    %basic(autoCreate=.true.)
    basic2 => node2    %basic(autoCreate=.true.)
    call basic1%massSet(mass1)
    call basic1%timeSet(time1)
    call basic2%massSet(mass2)
    call basic2%timeSet(time2)
    massError1=basic1%mass()*self%haloMassError_%errorFractional(node1)
    massError2=basic2%mass()*self%haloMassError_%errorFractional(node2)
    call node1%destroy()
    call node2%destroy()
    deallocate(node1)
    deallocate(node2)
    ! Handle zero errors.
    if (massError1 <= 0.0d0 .or. massError2 <= 0.0d0) then
       ! We currently do not handle cases where only one error is zero.
       if (massError1 > 0.0d0 .or. massError2 > 0.0d0) call Galacticus_Error_Report('conditionalMFBinWeights2D','both mass errors must be zero or both must be non-zero')
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
                     &          *massError2
                mass2UpperLimit=+exp(                                       &
                     &               +massRatioLogarithmicMinimumBins2      &
                     &               +dble(j  )                             &
                     &               /massRatioLogarithmicWidthInverseBins2 &
                     &              )                                       &
                     &          *mass1UpperLimit                            &
                     &          +integrationExtent                          &
                     &          *massError2
                if     (                         &
                     &   mass2LowerLimit > mass2 &
                     &  .or.                     &
                     &   mass2UpperLimit < mass2 &
                     & ) then                   
                   conditionalMFBinWeights2D(i,j)=0.0d0 
                else
                   conditionalMFBinWeights2D(i,j)=max(                                                      &
                        &                             Integrate(                                            &
                        &                                       mass1LowerLimit                           , &
                        &                                       mass1UpperLimit                           , &
                        &                                       conditionalMFBinWeights2DIntegrand        , &
                        &                                       integrandFunction                         , &
                        &                                       integrationWorkspace                      , &
                        &                                       toleranceAbsolute                 =1.0d-10, &
                        &                                       toleranceRelative                 =1.0d-03  &
                        &                                      )                                          , &
                        &                             0.0d0                                                 &
                        &                            )
                   call Integrate_Done(integrandFunction,integrationWorkspace)
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
      !% Integrand used in finding the weight given to a bin in the space of parent mass vs. progenitor mass ratio.
      use Numerical_Constants_Math
      use Galacticus_Error
      implicit none
      double precision, intent(in   ) :: mass1Primed
      double precision                :: mass2LowerLimit, mass2UpperLimit
      
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
      ! Evaluate integrand for the relevant moment.
      select case (moment)
      case (0)
         conditionalMFBinWeights2DIntegrand=+(                        &
              &                               +erf(                   &
              &                                    +(                 &
              &                                      +mass2UpperLimit &
              &                                      -mass2           &
              &                                     )                 &
              &                                    /sqrt(2.0d0)       &
              &                                    /massError2        &
              &                                   )                   &
              &                               -erf(                   &
              &                                    +(                 &
              &                                      +mass2LowerLimit &
              &                                      -mass2           &
              &                                     )                 &
              &                                    /sqrt(2.0d0)       &
              &                                    /massError2        &
              &                                   )                   &
              &                              )                        &
              &                             /2.0d0
      case (1)
         conditionalMFBinWeights2DIntegrand=+(                          &
              &                               +mass2                    &
              &                               *erf(                     &
              &                                    +(                   &
              &                                      +mass2UpperLimit   &
              &                                      -mass2             &
              &                                     )                   &
              &                                    /sqrt(2.0d0)         &
              &                                    /massError2          &
              &                                   )                     &
              &                               /2.0d0                    &
              &                               -massError2               &
              &                               *exp(                     &
              &                                    -(                   &
              &                                      +(                 &
              &                                        +mass2UpperLimit &
              &                                        -mass2           &
              &                                       )                 &
              &                                      /massError2        &
              &                                     )**2                &
              &                                    /2.0d0               &
              &                                   )                     &
              &                               /sqrt(                    &
              &                                     +2.0d0              &
              &                                     *Pi                 &
              &                                    )                    &
              &                               -mass2                    &
              &                               *erf(                     &
              &                                    +(                   &
              &                                      +mass2LowerLimit   &
              &                                      -mass2             &
              &                                     )                   &
              &                                    /sqrt(2.0d0)         &
              &                                    /massError2          &
              &                                   )                     &
              &                               /2.0d0                    &
              &                               +massError2               &
              &                               *exp(                     &
              &                                    -(                   &
              &                                      +(                 &
              &                                        +mass2LowerLimit &
              &                                        -mass2           &
              &                                       )                 &
              &                                      /massError2        &
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
              &                               -massError2               &
              &                               *(                        &
              &                                 +mass2UpperLimit        &
              &                                 +mass2                  &
              &                                )                        &
              &                               *exp(                     &
              &                                    -(                   &
              &                                      +(                 &
              &                                        +mass2UpperLimit &
              &                                        -mass2           &
              &                                       )                 &
              &                                      /massError2        &
              &                                     )**2                &
              &                                    /2.0d0               &
              &                                   )                     &
              &                               /sqrt(                    &
              &                                     +2.0d0              &
              &                                     *Pi                 &
              &                                    )                    &
              &                               +(                        &
              &                                 +massError2**2          &
              &                                 +mass2     **2          &
              &                                )                        &
              &                               *erf(                     &
              &                                    +(                   &
              &                                      +mass2UpperLimit   &
              &                                      -mass2             &
              &                                     )                   &
              &                                    /sqrt(2.0d0)         &
              &                                    /massError2          &
              &                                   )                     &
              &                               /2.0d0                    &
              &                               +massError2               &
              &                               *(                        &
              &                                 +mass2LowerLimit        &
              &                                 +mass2                  &
              &                                )                        &
              &                               *exp(                     &
              &                                    -(                   &
              &                                      +(                 &
              &                                        +mass2LowerLimit &
              &                                        -mass2           &
              &                                       )                 &
              &                                      /massError2        &
              &                                     )**2                &
              &                                    /2.0d0               &
              &                                   )                     &
              &                               /sqrt(                    &
              &                                     +2.0d0              &
              &                                     *Pi                 &
              &                                    )                    &
              &                               -(                        &
              &                                 +massError2**2          &
              &                                 +mass2     **2          &
              &                                )                        &
              &                               *erf(                     &
              &                                    +(                   &
              &                                      +mass2LowerLimit   &
              &                                      -mass2             &
              &                                     )                   &
              &                                    /sqrt(2.0d0)         &
              &                                    /massError2          &
              &                                   )                     &
              &                               /2.0d0                    &
              &                              )
      case default
         conditionalMFBinWeights2DIntegrand=0.0d0
         call Galacticus_Error_Report('conditionalMFBinWeights2DIntegrand','moment not supported')
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
    !% Outputs conditional mass function.
    use ISO_Varying_String
    use IO_HDF5
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    use Memory_Management
    implicit none
    class           (mergerTreeOperatorConditionalMF), intent(inout)                     :: self
    type            (hdf5Object                     )                                    :: conditionalMassFunctionGroup       , massDataset
    double precision                                 , allocatable  , dimension(:      ) :: normalizationSubhaloMassFunction
    double precision                                 , allocatable  , dimension(:,:    ) :: normalization
    double precision                                 , allocatable  , dimension(:,:,:  ) :: conditionalMassFunction            , conditionalMassFunctionError      , &
         &                                                                                  conditionalMassFunctionWeight
    double precision                                 , allocatable  , dimension(:,:,:  ) :: subhaloMassFunction                , subhaloMassFunctionError          , &
         &                                                                                  subhaloMassFunctionWeight
    double precision                                 , allocatable  , dimension(:,:,:,:) :: primaryProgenitorMassFunction      , primaryProgenitorMassFunctionError, &
         &                                                                                  primaryProgenitorMassFunctionWeight
    double precision                                 , allocatable  , dimension(:,:,:,:) :: formationRateFunction              , formationRateFunctionError        , &
         &                                                                                  formationRateFunctionWeight
    integer                                                                              :: i                                  , j                                 , &
         &                                                                                  iPrimary

    ! Compute normalizations.
    self                             %normalization                   =self%normalization                   /self%massRatioLogarithmicBinWidthInverse/log(10.0d0)
    if (self%extendedStatistics) self%normalizationSubhaloMassFunction=self%normalizationSubhaloMassFunction/self%massRatioLogarithmicBinWidthInverse/log(10.0d0)
    ! Output the data.
    !$omp critical(HDF5_Access)
    ! Check if our output group already exists.
    if (galacticusOutputFile%hasGroup(char(self%outputGroupName))) then
       ! Our group does exist. Read existing mass functions, add them to our own, then write back to file.
       conditionalMassFunctionGroup=galacticusOutputFile%openGroup(char(self%outputGroupName),'Conditional mass functions of merger trees.',objectsOverwritable=.true.,overwriteOverride=.true.)
       call allocateArray(normalization               ,shape(self%normalization               ))
       call allocateArray(conditionalMassFunction     ,shape(self%conditionalMassFunction     ))
       call allocateArray(conditionalMassFunctionError,shape(self%conditionalMassFunctionError))
       call conditionalMassFunctionGroup%readDataset('normalization'               ,normalization               )
       call conditionalMassFunctionGroup%readDataset('conditionalMassFunction'     ,conditionalMassFunction     )
       call conditionalMassFunctionGroup%readDataset('conditionalMassFunctionError',conditionalMassFunctionError)
       if (self%extendedStatistics) then
          call allocateArray(normalizationSubhaloMassFunction  ,shape(self%normalizationSubhaloMassFunction  ))
          call allocateArray(primaryProgenitorMassFunction     ,shape(self%primaryProgenitorMassFunction     ))
          call allocateArray(primaryProgenitorMassFunctionError,shape(self%primaryProgenitorMassFunctionError))
          call allocateArray(formationRateFunction             ,shape(self%formationRateFunction             ))
          call allocateArray(formationRateFunctionError        ,shape(self%formationRateFunctionError        ))
          call allocateArray(subhaloMassFunction               ,shape(self%subhaloMassFunction               ))
          call allocateArray(subhaloMassFunctionError          ,shape(self%subhaloMassFunctionError          ))
          call conditionalMassFunctionGroup%readDataset('normalizationSubhaloMassFunction'  ,normalizationSubhaloMassFunction  )
          call conditionalMassFunctionGroup%readDataset('primaryProgenitorMassFunction'     ,primaryProgenitorMassFunction     )
          call conditionalMassFunctionGroup%readDataset('primaryProgenitorMassFunctionError',primaryProgenitorMassFunctionError)
          call conditionalMassFunctionGroup%readDataset('formationRateFunction'             ,formationRateFunction             )
          call conditionalMassFunctionGroup%readDataset('formationRateFunctionError'        ,formationRateFunctionError        )
          call conditionalMassFunctionGroup%readDataset('subhaloMassFunction'               ,subhaloMassFunction               )
          call conditionalMassFunctionGroup%readDataset('subhaloMassFunctionError'          ,subhaloMassFunctionError          )
       end if
       ! Accumulate the conditional mass functions.
       do i=1,self%timeCount
          do j=1,self%parentMassCount
             self%conditionalMassFunction     (i,j,:  )=self%conditionalMassFunction     (i,j,:  )+ conditionalMassFunction     (i,j,:  )*normalization(i,j) 
             self%conditionalMassFunctionError(i,j,:  )=self%conditionalMassFunctionError(i,j,:  )+(conditionalMassFunctionError(i,j,:  )*normalization(i,j))**2 
             if (self%extendedStatistics) then
                self%formationRateFunction       (i,j,:,:)=self%formationRateFunction       (i,j,:,:)+ formationRateFunction       (i,j,:,:)*normalization(i,j) 
                self%formationRateFunctionError  (i,j,:,:)=self%formationRateFunctionError  (i,j,:,:)+(formationRateFunctionError  (i,j,:,:)*normalization(i,j))**2 
                if (self%primaryProgenitorStatisticsValid) then
                   do iPrimary=1,self%primaryProgenitorDepth
                      self%primaryProgenitorMassFunction     (i,j,:,iPrimary)=self%primaryProgenitorMassFunction     (i,j,:,iPrimary)+ primaryProgenitorMassFunction     (i,j,:,iPrimary)*normalization(i,j)
                      self%primaryProgenitorMassFunctionError(i,j,:,iPrimary)=self%primaryProgenitorMassFunctionError(i,j,:,iPrimary)+(primaryProgenitorMassFunctionError(i,j,:,iPrimary)*normalization(i,j))**2 
                   end do
                else
                   ! Primary progenitor statistics are not valid (because of non-zero errors on halo masses), so set to an unphysical value.
                   self%primaryProgenitorMassFunction     (:,:,:,:)=-1.0d0
                   self%primaryProgenitorMassFunctionError(:,:,:,:)=-1.0d0
                end if
             end if
          end do
       end do
       ! Accumulate subhalo mass functions. 
       if (self%extendedStatistics) then
          do j=1,self%parentMassCount
             self%subhaloMassFunction     (j,:,:)=self%subhaloMassFunction     (j,:,:)+ subhaloMassFunction     (j,:,:)*normalizationSubhaloMassFunction(j)
             self%subhaloMassFunctionError(j,:,:)=self%subhaloMassFunctionError(j,:,:)+(subhaloMassFunctionError(j,:,:)*normalizationSubhaloMassFunction(j))**2 
          end do
       end if
       ! Accumulate normalizations.
       self                             %normalization                   =self%normalization                   +normalization
       if (self%extendedStatistics) self%normalizationSubhaloMassFunction=self%normalizationSubhaloMassFunction+normalizationSubhaloMassFunction
       call    deallocateArray(normalization                     )
       call    deallocateArray(conditionalMassFunction           )
       call    deallocateArray(conditionalMassFunctionError      )
       if (self%extendedStatistics) then
          call deallocateArray(normalizationSubhaloMassFunction  )
          call deallocateArray(primaryProgenitorMassFunction     )
          call deallocateArray(primaryProgenitorMassFunctionError)
          call deallocateArray(formationRateFunction             )
          call deallocateArray(formationRateFunctionError        )
          call deallocateArray(subhaloMassFunction               )
          call deallocateArray(subhaloMassFunctionError          )
       end if
    else
       ! Our group does not already exist. Simply write the data.
       conditionalMassFunctionGroup=galacticusOutputFile%openGroup(char(self%outputGroupName),'Conditional mass functions of merger trees.',objectsOverwritable=.true.,overwriteOverride=.true.)
       call conditionalMassFunctionGroup%writeDataset  (self%massParents                       ,"massParent"                        ,"Mass of parent node [Msolar]"              ,datasetReturned=massDataset)
       call massDataset                 %writeAttribute(massSolar                              ,"unitsInSI"                                                                                                  )
       call massDataset                 %close         (                                                                                                                                                     )
       call conditionalMassFunctionGroup%writeDataset  (self%massRatios                        ,"massRatio"                         ,"Mass of ratio node [Msolar]"               ,datasetReturned=massDataset)
       call massDataset                 %writeAttribute(massSolar                              ,"unitsInSI"                                                                                                  )
       call massDataset                 %close         (                                                                                                                                                     )
       call conditionalMassFunctionGroup%writeDataset  (self%parentRedshifts                   ,"redshiftParent"                    ,"Redshift of parent node []"                                            )
       call conditionalMassFunctionGroup%writeDataset  (self%progenitorRedshifts               ,"redshiftProgenitor"                ,"Redshift of progenitor node []"                                        )
    end if
    ! Create weight datasets.
    call    allocateArray(conditionalMassFunctionWeight      ,shape(self%conditionalMassFunction      ))
    conditionalMassFunctionWeight         =0.0d0
    if (self%extendedStatistics) then
       call allocateArray(primaryProgenitorMassFunctionWeight,shape(self%primaryProgenitorMassFunction))
       call allocateArray(formationRateFunctionWeight        ,shape(self%formationRateFunction        ))
       call allocateArray(subhaloMassFunctionWeight          ,shape(self%subhaloMassFunction          ))
       primaryProgenitorMassFunctionWeight=0.0d0
       formationRateFunctionWeight        =0.0d0
       subhaloMassFunctionWeight          =0.0d0
    end if
    ! Normalize the conditional mass functions.
    do i=1,self%timeCount
       do j=1,self%parentMassCount
          if (self%normalization(i,j) > 0.0d0) then
             self%conditionalMassFunction      (i,j,:  )=     self%conditionalMassFunction     (i,j,:  ) /self%normalization(i,j)
             self%conditionalMassFunctionError (i,j,:  )=sqrt(self%conditionalMassFunctionError(i,j,:  ))/self%normalization(i,j)
             conditionalMassFunctionWeight     (i,j,:  )=                                                 self%normalization(i,j)
             if (self%extendedStatistics) then
                self%formationRateFunction        (i,j,:,:)=     self%formationRateFunction       (i,j,:,:) /self%normalization(i,j)
                self%formationRateFunctionError   (i,j,:,:)=sqrt(self%formationRateFunctionError  (i,j,:,:))/self%normalization(i,j)
                formationRateFunctionWeight       (i,j,:,:)=                                                 self%normalization(i,j)
                if (self%primaryProgenitorStatisticsValid) then
                   do iPrimary=1,self%primaryProgenitorDepth
                      self%primaryProgenitorMassFunction      (i,j,:,iPrimary)=     self%primaryProgenitorMassFunction     (i,j,:,iPrimary) /self%normalization(i,j)
                      self%primaryProgenitorMassFunctionError (i,j,:,iPrimary)=sqrt(self%primaryProgenitorMassFunctionError(i,j,:,iPrimary))/self%normalization(i,j)
                      primaryProgenitorMassFunctionWeight     (i,j,:,iPrimary)=                                                              self%normalization(i,j)
                   end do
                end if
             end if
          end if
       end do
    end do
    ! Normalize subhalo mass functions.
    if (self%extendedStatistics) then
       do j=1,self%parentMassCount
          if (self%normalizationSubhaloMassFunction(j) > 0.0d0) then
             self%subhaloMassFunction      (j,:,:)=     self%subhaloMassFunction     (j,:,:)/self%normalizationSubhaloMassFunction(j)
             self%subhaloMassFunctionError (j,:,:)=sqrt(self%subhaloMassFunctionError(j,:,:)/self%normalizationSubhaloMassFunction(j))
             subhaloMassFunctionWeight     (j,:,:)=                                          self%normalizationSubhaloMassFunction(j)
          end if
       end do
    end if
    call    conditionalMassFunctionGroup%writeDataset  (self%normalization                      ,"normalization"                      ,"Normalization for conditional mass functions []")
    call    conditionalMassFunctionGroup%writeDataset  (self%conditionalMassFunction            ,"conditionalMassFunction"            ,"Conditional mass functions []"                  )
    call    conditionalMassFunctionGroup%writeDataset  (self%conditionalMassFunctionError       ,"conditionalMassFunctionError"       ,"Conditional mass function errors []"            )
    call    conditionalMassFunctionGroup%writeDataset  (     conditionalMassFunctionWeight      ,"conditionalMassFunctionWeight"      ,"Conditional mass function weights []"           )
    if (self%extendedStatistics) then
       call conditionalMassFunctionGroup%writeDataset  (self%normalizationSubhaloMassFunction   ,"normalizationSubhaloMassFunction"   ,"Normalization for subhalo mass functions []"    )
       call conditionalMassFunctionGroup%writeDataset  (self%primaryProgenitorMassFunction      ,"primaryProgenitorMassFunction"      ,"Primary progenitor mass functions []"           )
       call conditionalMassFunctionGroup%writeDataset  (self%primaryProgenitorMassFunctionError ,"primaryProgenitorMassFunctionError" ,"Primary progenitor mass function errors []"     )
       call conditionalMassFunctionGroup%writeDataset  (     primaryProgenitorMassFunctionWeight,"primaryProgenitorMassFunctionWeight","Primary progenitor mass function weights []"    )
       call conditionalMassFunctionGroup%writeDataset  (self%formationRateFunction              ,"formationRateFunction"              ,"Formation rate functions []"                    )
       call conditionalMassFunctionGroup%writeDataset  (self%formationRateFunctionError         ,"formationRateFunctionError"         ,"Formation rate function errors []"              )
       call conditionalMassFunctionGroup%writeDataset  (     formationRateFunctionWeight        ,"formationRateFunctionWeight"        ,"Formation rate function weights []"             )
       call conditionalMassFunctionGroup%writeDataset  (self%subhaloMassFunction                ,"subhaloMassFunction"                ,"Unevolved subhalo mass functions []"            )
       call conditionalMassFunctionGroup%writeDataset  (self%subhaloMassFunctionError           ,"subhaloMassFunctionError"           ,"Unevolved subhalo mass function errors []"      )
       call conditionalMassFunctionGroup%writeDataset  (     subhaloMassFunctionWeight          ,"subhaloMassFunctionWeight"          ,"Unevolved subhalo mass function weights []"     )
    end if
    call    conditionalMassFunctionGroup%close         (                                                                                                                                )    
    call galacticusOutputFile        %flush         (                                                                                                                                                     )
    !$omp end critical(HDF5_Access)
    ! Deallocate weight arrays.
    call deallocateArray(conditionalMassFunctionWeight      )
    if (self%extendedStatistics) then
       call deallocateArray(primaryProgenitorMassFunctionWeight)
       call deallocateArray(formationRateFunctionWeight        )
       call deallocateArray(subhaloMassFunctionWeight          )
    end if
    return
  end subroutine conditionalMFFinalize
