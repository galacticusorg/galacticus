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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program which tests the attenuation scatter and random error output analysis distribution operators.
!!}

program Test_Output_Analyses_Attenuation_Scatter
  !!{RST
  Tests the :galacticus-class:`outputAnalysisDistributionOperatorAttenuationScatter` distribution operator, and the
  application of a random error distribution operator to a distribution, each against an independent calculation.
  !!}
  use            :: Cosmology_Functions                  , only : cosmologyFunctionsClass
  use            :: Display                              , only : displayVerbositySet                                 , verbosityLevelStandard
  use            :: Dust_Attenuations                    , only : dustAttenuationClass                                , dustAttenuationStellarMassRedshift
  use            :: Error                                , only : Error_Handler_Register
  use            :: Events_Hooks                         , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities           , only : Functions_Global_Set
  use            :: Galacticus_Nodes                     , only : mergerTree                                          , nodeClassHierarchyInitialize                     , nodeComponentBasic, nodeComponentDisk, &
       &                                                          treeNode
  use            :: Input_Parameters                     , only : inputParameters
  use, intrinsic :: ISO_C_Binding                        , only : c_size_t
  use            :: Node_Components                      , only : Node_Components_Initialize                          , Node_Components_Thread_Initialize
  use            :: Output_Analyses_Options              , only : outputAnalysisPropertyTypeLog10
  use            :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorAttenuationScatter, outputAnalysisDistributionOperatorRandomErrorFixed
  use            :: Unit_Tests                           , only : Assert                                              , Unit_Tests_Begin_Group                           , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                                                               , parameter                   :: countBins           =120    , countQuantiles=200000, &
       &                                                                                                 countSamples        =2000
  double precision                                                      , parameter                   :: propertyMinimum     =40.0d0 , widthBin      =0.05d0, &
       &                                                                                                 rootVarianceScatter =0.25d0 , rootVarianceError=0.1d0
  type            (inputParameters                                     )                              :: parameters
  class           (dustAttenuationClass                                ), pointer                     :: dustAttenuation_
  class           (cosmologyFunctionsClass                             ), pointer                     :: cosmologyFunctions_
  type            (outputAnalysisDistributionOperatorAttenuationScatter)                              :: scatter           , scatterNone
  type            (outputAnalysisDistributionOperatorRandomErrorFixed  )                              :: randomError
  type            (mergerTree                                          ), target                      :: tree
  type            (treeNode                                            ), pointer                     :: node
  double precision                                                      , dimension(countBins)        :: propertyValueMinimum, propertyValueMaximum, &
       &                                                                                                 distribution        , distributionExpected, &
       &                                                                                                 distributionSource
  double precision                                                                                    :: attenuationMean     , luminosityAttenuated, &
       &                                                                                                 luminosityUnattenuated, epsilon            , &
       &                                                                                                 time                , value
  integer                                                                                             :: i                   , j                   , &
       &                                                                                                 k
  integer         (c_size_t                                            ), parameter                   :: outputIndex         =1_c_size_t

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  call Unit_Tests_Begin_Group("Output analysis distribution operators: attenuation scatter and random error")
  parameters=inputParameters('testSuite/parameters/dustAttenuationStellarMassRedshift.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="dustAttenuation"    name="dustAttenuation_"    source="parameters"/>
  <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
  !!]
  ! Bins in log10 luminosity.
  do i=1,countBins
     propertyValueMinimum(i)=propertyMinimum+dble(i-1)*widthBin
     propertyValueMaximum(i)=propertyMinimum+dble(i  )*widthBin
  end do
  ! A galaxy of 10^9 Msun at z=1.5. With the parameters used, its mean attenuation is positive but small compared to the scatter,
  ! so that there is a substantial probability of zero attenuation.
  time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(1.5d0))
  node  => treeNode(hostTree=tree)
  block
    class(nodeComponentBasic), pointer :: basic
    class(nodeComponentDisk ), pointer :: disk
    basic => node%basic(autoCreate=.true.)
    disk  => node%disk (autoCreate=.true.)
    call basic%massSet       (1.0d12)
    call basic%timeSet       (time  )
    call disk %massStellarSet(1.0d9 )
    call disk %radiusSet     (3.0d-3)
  end block
  select type (dustAttenuation_)
  class is (dustAttenuationStellarMassRedshift)
     attenuationMean=dustAttenuation_%attenuationMean(node)
  class default
     attenuationMean=0.0d0
     call Assert('attenuation is of the expected class',.false.,.true.)
  end select
  call Assert('mean attenuation is positive, and comparable to the scatter',attenuationMean > 0.0d0 .and. attenuationMean < 2.0d0*rootVarianceScatter,.true.)
  luminosityUnattenuated=42.0d0+0.3d0*widthBin
  luminosityAttenuated  =luminosityUnattenuated-0.4d0*max(attenuationMean,0.0d0)

  ! Attenuation scatter, against a histogram of the attenuated luminosity evaluated at evenly spaced quantiles of the scatter.
  call Unit_Tests_Begin_Group("Attenuation scatter")
  scatter    =outputAnalysisDistributionOperatorAttenuationScatter(rootVarianceScatter,dustAttenuation_)
  scatterNone=outputAnalysisDistributionOperatorAttenuationScatter(0.0d0              ,dustAttenuation_)
  distribution=scatter%operateScalar(luminosityAttenuated,outputAnalysisPropertyTypeLog10,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
  distributionExpected=0.0d0
  do k=1,countQuantiles
     epsilon=rootVarianceScatter*normalInverse((dble(k)-0.5d0)/dble(countQuantiles))
     value  =luminosityUnattenuated-0.4d0*max(attenuationMean+epsilon,0.0d0)
     j      =floor((value-propertyMinimum)/widthBin)+1
     if (j >= 1 .and. j <= countBins) distributionExpected(j)=distributionExpected(j)+1.0d0/dble(countQuantiles)
  end do
  call Assert('probability is conserved'                     ,sum(distribution),1.0d0,absTol=1.0d-9)
  call Assert('distribution matches a sampled histogram'     ,distribution,distributionExpected,absTol=2.0d-4)
  j=floor((luminosityUnattenuated-propertyMinimum)/widthBin)+1
  call Assert('probability of zero attenuation is in the unattenuated bin',distribution(j) >= 0.5d0*erfc(attenuationMean/rootVarianceScatter/sqrt(2.0d0)),.true.)
  call Assert('no probability above the unattenuated bin'    ,maxval(distribution(j+1:)),0.0d0)
  distribution=scatterNone%operateScalar(luminosityAttenuated,outputAnalysisPropertyTypeLog10,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
  j=floor((luminosityAttenuated-propertyMinimum)/widthBin)+1
  call Assert('without scatter, all weight is in the attenuated bin',[sum(distribution),distribution(j)],[1.0d0,1.0d0])
  call Unit_Tests_End_Group()

  ! Random error applied to a distribution, against operateScalar averaged over points spread uniformly across each source bin.
  call Unit_Tests_Begin_Group("Random error applied to a distribution")
  randomError=outputAnalysisDistributionOperatorRandomErrorFixed(rootVarianceError)
  distributionSource=0.0d0
  distributionSource(40)=0.7d0
  distributionSource(41)=0.3d0
  distribution=randomError%operateDistribution(distributionSource,outputAnalysisPropertyTypeLog10,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
  distributionExpected=0.0d0
  do j=40,41
     do k=1,countSamples
        value               =propertyValueMinimum(j)+widthBin*(dble(k)-0.5d0)/dble(countSamples)
        distributionExpected=+distributionExpected                                                                                                         &
             &               +distributionSource(j)                                                                                                         &
             &               *randomError%operateScalar(value,outputAnalysisPropertyTypeLog10,propertyValueMinimum,propertyValueMaximum,outputIndex,node) &
             &               /dble(countSamples)
     end do
  end do
  call Assert('weight is conserved'                          ,sum(distribution),1.0d0,absTol=1.0d-9)
  call Assert('distribution matches averaged scalar operation',distribution,distributionExpected,absTol=1.0d-7)
  call Unit_Tests_End_Group()

  !![
  <objectDestructor name="dustAttenuation_"   />
  <objectDestructor name="cosmologyFunctions_"/>
  !!]
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

contains

  double precision function normalInverse(p)
    !!{RST
    Invert the cumulative distribution of the standard normal distribution, by bisection.
    !!}
    implicit none
    double precision, intent(in   ) :: p
    double precision                :: xLow, xHigh, x
    integer                         :: iteration

    xLow =-10.0d0
    xHigh=+10.0d0
    do iteration=1,100
       x=0.5d0*(xLow+xHigh)
       if (0.5d0*erfc(-x/sqrt(2.0d0)) < p) then
          xLow =x
       else
          xHigh=x
       end if
    end do
    normalInverse=0.5d0*(xLow+xHigh)
    return
  end function normalInverse

end program Test_Output_Analyses_Attenuation_Scatter
