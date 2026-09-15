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

program Test_Output_Analysis_Capabilities
  !!{RST
  Tests the capability queries by which an output analysis discovers that a property operator or
  distribution normalizer can not be used in a given role. An ``unoperator`` is applied to bin centers
  once the model has run, where no node and no output index exist, so an operator requiring either can
  not serve in that role; a cross-correlation accumulates no binned distribution, so a normalizer which
  operates on the distribution can not be used with it. Both queries must also see through the
  ``sequence`` classes, or wrapping an operator in a sequence would evade the check.
  !!}
  use :: Display                                 , only : displayVerbositySet                          , verbosityLevelStandard
  use :: Output_Analysis_Molecular_Ratios        , only : outputAnalysisMolecularRatioObreschkow2009   , outputAnalysisMolecularRatioClass
  use :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerIdentity , outputAnalysisDistributionNormalizerSequence, &
       &                                                  outputAnalysisDistributionNormalizerUnitarity, normalizerList                              , &
       &                                                  outputAnalysisDistributionNormalizerClass
  use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10      , outputAnalysisPropertyOperatorHIMass        , &
       &                                                  outputAnalysisPropertyOperatorSequence       , propertyOperatorList                        , &
       &                                                  outputAnalysisPropertyOperatorClass
  use :: Unit_Tests                              , only : Assert                                       , Unit_Tests_Begin_Group                      , &
       &                                                  Unit_Tests_End_Group                         , Unit_Tests_Finish
  implicit none
  class(outputAnalysisMolecularRatioClass            ), pointer :: molecularRatio_
  class(outputAnalysisPropertyOperatorClass          ), pointer :: operatorAntiLog10_       , operatorHIMass_                 , &
       &                                                           operatorAntiLog10Second_
  class(outputAnalysisDistributionNormalizerClass    ), pointer :: normalizerIdentity_      , normalizerUnitarity_            , &
       &                                                           normalizerIdentitySecond_
  type (outputAnalysisPropertyOperatorSequence       )          :: operatorSequencePlain_   , operatorSequenceWithHIMass_
  type (outputAnalysisDistributionNormalizerSequence )          :: normalizerSequencePlain_ , normalizerSequenceWithUnitarity_
  type (propertyOperatorList                         ), pointer :: operatorsPlain           , operatorsWithHIMass
  type (normalizerList                               ), pointer :: normalizersPlain         , normalizersWithUnitarity

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group('Output analysis capabilities')
  ! Build the operators. The molecular ratio parameters are arbitrary - only the capability is under test.
  allocate(outputAnalysisMolecularRatioObreschkow2009 :: molecularRatio_          )
  allocate(outputAnalysisPropertyOperatorAntiLog10    :: operatorAntiLog10_       )
  allocate(outputAnalysisPropertyOperatorAntiLog10    :: operatorAntiLog10Second_ )
  allocate(outputAnalysisPropertyOperatorHIMass       :: operatorHIMass_          )
  select type (molecularRatio_)
  type is (outputAnalysisMolecularRatioObreschkow2009)
     molecularRatio_=outputAnalysisMolecularRatioObreschkow2009(0.4d0,0.5d0,0.2d0,1.6d0,0.8d0,0.8d0,0.8d0,0.0d0)
  end select
  select type (operatorHIMass_)
  type is (outputAnalysisPropertyOperatorHIMass)
     operatorHIMass_=outputAnalysisPropertyOperatorHIMass(molecularRatio_)
  end select
  ! One sequence of operators which need nothing, and one which also contains an operator needing a node. Each
  ! sequence owns its list, which its destructor deallocates, so the two lists are allocated separately.
  allocate(operatorsPlain     )
  allocate(operatorsWithHIMass)
  allocate(operatorsWithHIMass%next)
  operatorsPlain          %operator_ => operatorAntiLog10_
  operatorsWithHIMass     %operator_ => operatorAntiLog10Second_
  operatorsWithHIMass%next%operator_ => operatorHIMass_
  operatorSequencePlain_     =outputAnalysisPropertyOperatorSequence(operatorsPlain     )
  operatorSequenceWithHIMass_=outputAnalysisPropertyOperatorSequence(operatorsWithHIMass)

  call Unit_Tests_Begin_Group('Property operators')
  call Assert('antiLog10 needs no node'                          ,operatorAntiLog10_         %isNodeDependent  (),.false.)
  call Assert('antiLog10 needs no output index'                  ,operatorAntiLog10_         %isOutputDependent(),.false.)
  call Assert('HI mass needs a node'                             ,operatorHIMass_            %isNodeDependent  (),.true. )
  call Assert('HI mass needs no output index'                    ,operatorHIMass_            %isOutputDependent(),.false.)
  call Assert('a sequence of independent operators needs no node',operatorSequencePlain_     %isNodeDependent  (),.false.)
  call Assert('a sequence containing HI mass needs a node'       ,operatorSequenceWithHIMass_%isNodeDependent  (),.true. )
  call Unit_Tests_End_Group()

  ! Build the normalizers, and a sequence containing one which operates on the distribution itself.
  allocate(outputAnalysisDistributionNormalizerIdentity  :: normalizerIdentity_      )
  allocate(outputAnalysisDistributionNormalizerIdentity  :: normalizerIdentitySecond_)
  allocate(outputAnalysisDistributionNormalizerUnitarity :: normalizerUnitarity_     )
  allocate(normalizersPlain        )
  allocate(normalizersWithUnitarity)
  allocate(normalizersWithUnitarity%next)
  normalizersPlain             %normalizer_ => normalizerIdentity_
  normalizersWithUnitarity     %normalizer_ => normalizerIdentitySecond_
  normalizersWithUnitarity%next%normalizer_ => normalizerUnitarity_
  normalizerSequencePlain_        =outputAnalysisDistributionNormalizerSequence(normalizersPlain        )
  normalizerSequenceWithUnitarity_=outputAnalysisDistributionNormalizerSequence(normalizersWithUnitarity)

  call Unit_Tests_Begin_Group('Distribution normalizers')
  call Assert('identity needs no distribution'                        ,normalizerIdentity_             %requiresDistribution(),.false.)
  call Assert('unitarity needs the distribution'                      ,normalizerUnitarity_            %requiresDistribution(),.true. )
  call Assert('a sequence of independent normalizers needs none'      ,normalizerSequencePlain_        %requiresDistribution(),.false.)
  call Assert('a sequence containing unitarity needs the distribution',normalizerSequenceWithUnitarity_%requiresDistribution(),.true. )
  call Unit_Tests_End_Group()

  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Output_Analysis_Capabilities
