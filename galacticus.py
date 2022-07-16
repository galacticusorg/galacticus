from ctypes import *
# Load the shared library into ctypes.
libname = "./galacticus/lib/libgalacticus.so"
c_lib = CDLL(libname)
c_lib.libGalacticusInitL()
c_lib.cosmologyFunctionsMatterDarkEnergyL.restype  = c_void_p
c_lib.cosmologyFunctionsMatterDarkEnergyL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.cosmologyFunctionsMatterLambdaL.restype  = c_void_p
c_lib.cosmologyFunctionsMatterLambdaL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsStaticUniverseL.restype  = c_void_p
c_lib.cosmologyFunctionsStaticUniverseL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsComovingVolumeElementRedshiftL.restype  = c_double
c_lib.cosmologyFunctionsComovingVolumeElementRedshiftL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsComovingVolumeElementTimeL.restype  = c_double
c_lib.cosmologyFunctionsComovingVolumeElementTimeL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsCosmicTimeL.restype  = c_double
c_lib.cosmologyFunctionsCosmicTimeL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsDensityScalingEarlyTimeL.argtypes = [ c_void_p, c_int, c_double, POINTER(c_double), POINTER(c_double) ]
c_lib.cosmologyFunctionsDistanceAngularL.restype  = c_double
c_lib.cosmologyFunctionsDistanceAngularL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsDistanceComovingL.restype  = c_double
c_lib.cosmologyFunctionsDistanceComovingL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsDistanceComovingConvertL.restype  = c_double
c_lib.cosmologyFunctionsDistanceComovingConvertL.argtypes = [ c_void_p, c_int, c_int ]
c_lib.cosmologyFunctionsDistanceLuminosityL.restype  = c_double
c_lib.cosmologyFunctionsDistanceLuminosityL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsDistanceParticleHorizonComovingL.restype  = c_double
c_lib.cosmologyFunctionsDistanceParticleHorizonComovingL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsDominationEpochMatterL.restype  = c_double
c_lib.cosmologyFunctionsDominationEpochMatterL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsEpochTimeL.restype  = c_double
c_lib.cosmologyFunctionsEpochTimeL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsEpochValidateL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsEqualityEpochMatterCurvatureL.restype  = c_double
c_lib.cosmologyFunctionsEqualityEpochMatterCurvatureL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsEqualityEpochMatterDarkEnergyL.restype  = c_double
c_lib.cosmologyFunctionsEqualityEpochMatterDarkEnergyL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsEqualityEpochMatterRadiationL.restype  = c_double
c_lib.cosmologyFunctionsEqualityEpochMatterRadiationL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsEquationOfStateDarkEnergyL.restype  = c_double
c_lib.cosmologyFunctionsEquationOfStateDarkEnergyL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsExpansionFactorL.restype  = c_double
c_lib.cosmologyFunctionsExpansionFactorL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsExpansionFactorFromRedshiftL.restype  = c_double
c_lib.cosmologyFunctionsExpansionFactorFromRedshiftL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsExpansionRateL.restype  = c_double
c_lib.cosmologyFunctionsExpansionRateL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsExponentDarkEnergyL.restype  = c_double
c_lib.cosmologyFunctionsExponentDarkEnergyL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsHubbleParameterEpochalL.restype  = c_double
c_lib.cosmologyFunctionsHubbleParameterEpochalL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL.restype  = c_double
c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsMatterDensityEpochalL.restype  = c_double
c_lib.cosmologyFunctionsMatterDensityEpochalL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL.restype  = c_double
c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsOmegaMatterEpochalL.restype  = c_double
c_lib.cosmologyFunctionsOmegaMatterEpochalL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL.restype  = c_double
c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsRedshiftFromExpansionFactorL.restype  = c_double
c_lib.cosmologyFunctionsRedshiftFromExpansionFactorL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsTemperatureCMBEpochalL.restype  = c_double
c_lib.cosmologyFunctionsTemperatureCMBEpochalL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsTimeAtDistanceComovingL.restype  = c_double
c_lib.cosmologyFunctionsTimeAtDistanceComovingL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.cosmologyFunctionsTimeBigCrunchL.restype  = c_double
c_lib.cosmologyFunctionsTimeBigCrunchL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyFunctionsDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyParametersSimpleL.restype  = c_void_p
c_lib.cosmologyParametersSimpleL.argtypes = [ c_double, c_double, c_double, c_double, c_double ]
c_lib.cosmologyParametersHubbleConstantL.restype  = c_double
c_lib.cosmologyParametersHubbleConstantL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyParametersOmegaBaryonL.restype  = c_double
c_lib.cosmologyParametersOmegaBaryonL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyParametersOmegaCurvatureL.restype  = c_double
c_lib.cosmologyParametersOmegaCurvatureL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyParametersOmegaDarkEnergyL.restype  = c_double
c_lib.cosmologyParametersOmegaDarkEnergyL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyParametersOmegaMatterL.restype  = c_double
c_lib.cosmologyParametersOmegaMatterL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyParametersOmegaRadiationL.restype  = c_double
c_lib.cosmologyParametersOmegaRadiationL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyParametersDensityCriticalL.restype  = c_double
c_lib.cosmologyParametersDensityCriticalL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyParametersTemperatureCMBL.restype  = c_double
c_lib.cosmologyParametersTemperatureCMBL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologyParametersDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.darkMatterParticleCDML.restype  = c_void_p
c_lib.darkMatterParticleCDML.argtypes = [  ]
c_lib.darkMatterParticleDecayingDarkMatterL.restype  = c_void_p
c_lib.darkMatterParticleDecayingDarkMatterL.argtypes = [ c_double, c_double, c_void_p, c_int ]
c_lib.darkMatterParticleFuzzyDarkMatterL.restype  = c_void_p
c_lib.darkMatterParticleFuzzyDarkMatterL.argtypes = [ c_double, c_double ]
c_lib.darkMatterParticleSelfInteractingDarkMatterL.restype  = c_void_p
c_lib.darkMatterParticleSelfInteractingDarkMatterL.argtypes = [ c_double, c_void_p, c_int ]
c_lib.darkMatterParticleWDMThermalL.restype  = c_void_p
c_lib.darkMatterParticleWDMThermalL.argtypes = [ c_double, c_double, c_void_p, c_int ]
c_lib.darkMatterParticleMassL.restype  = c_double
c_lib.darkMatterParticleMassL.argtypes = [ c_void_p, c_int ]
c_lib.darkMatterParticleDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.darkMatterHaloScaleVirialDensityContrastDefinitionL.restype  = c_void_p
c_lib.darkMatterHaloScaleVirialDensityContrastDefinitionL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterHaloScaleDensityMeanL.restype  = c_double
c_lib.darkMatterHaloScaleDensityMeanL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterHaloScaleDensityMeanGrowthRateL.restype  = c_double
c_lib.darkMatterHaloScaleDensityMeanGrowthRateL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterHaloScaleRadiusVirialL.restype  = c_double
c_lib.darkMatterHaloScaleRadiusVirialL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterHaloScaleRadiusVirialGradientLogarithmicMassL.restype  = c_double
c_lib.darkMatterHaloScaleRadiusVirialGradientLogarithmicMassL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterHaloScaleRadiusVirialGrowthRateL.restype  = c_double
c_lib.darkMatterHaloScaleRadiusVirialGrowthRateL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterHaloScaleTemperatureVirialL.restype  = c_double
c_lib.darkMatterHaloScaleTemperatureVirialL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterHaloScaleTimescaleDynamicalL.restype  = c_double
c_lib.darkMatterHaloScaleTimescaleDynamicalL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterHaloScaleVelocityVirialL.restype  = c_double
c_lib.darkMatterHaloScaleVelocityVirialL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterHaloScaleVelocityVirialGrowthRateL.restype  = c_double
c_lib.darkMatterHaloScaleVelocityVirialGrowthRateL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterHaloScaleDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.darkMatterProfileDMOBurkertL.restype  = c_void_p
c_lib.darkMatterProfileDMOBurkertL.argtypes = [ c_void_p, c_int ]
c_lib.darkMatterProfileDMOEinastoL.restype  = c_void_p
c_lib.darkMatterProfileDMOEinastoL.argtypes = [ c_void_p, c_int ]
c_lib.darkMatterProfileDMONFWL.restype  = c_void_p
c_lib.darkMatterProfileDMONFWL.argtypes = [ c_bool, c_void_p, c_int ]
c_lib.darkMatterProfileDMOPenarrubia2010L.restype  = c_void_p
c_lib.darkMatterProfileDMOPenarrubia2010L.argtypes = [ c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_void_p, c_int ]
c_lib.darkMatterProfileDMOSIDMCoreNFWL.restype  = c_void_p
c_lib.darkMatterProfileDMOSIDMCoreNFWL.argtypes = [ c_double, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMOSIDMIsothermalL.restype  = c_void_p
c_lib.darkMatterProfileDMOSIDMIsothermalL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMOZhao1996L.restype  = c_void_p
c_lib.darkMatterProfileDMOZhao1996L.argtypes = [ c_double, c_double, c_double, c_void_p, c_int ]
c_lib.darkMatterProfileDMOAcceleratorL.restype  = c_void_p
c_lib.darkMatterProfileDMOAcceleratorL.argtypes = [ c_double, c_double, c_void_p, c_int ]
c_lib.darkMatterProfileDMOAccretionFlowL.restype  = c_void_p
c_lib.darkMatterProfileDMOAccretionFlowL.argtypes = [ c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMODecayingL.restype  = c_void_p
c_lib.darkMatterProfileDMODecayingL.argtypes = [ c_void_p, c_int, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMOFiniteResolutionL.restype  = c_void_p
c_lib.darkMatterProfileDMOFiniteResolutionL.argtypes = [ c_double, c_double, c_bool, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMOFiniteResolutionNFWL.restype  = c_void_p
c_lib.darkMatterProfileDMOFiniteResolutionNFWL.argtypes = [ c_double, c_double, c_bool, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMOHeatedL.restype  = c_void_p
c_lib.darkMatterProfileDMOHeatedL.argtypes = [ c_int, c_bool, c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMOHeatedMonotonicL.restype  = c_void_p
c_lib.darkMatterProfileDMOHeatedMonotonicL.argtypes = [ c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_double, c_double, c_int ]
c_lib.darkMatterProfileDMOIsothermalL.restype  = c_void_p
c_lib.darkMatterProfileDMOIsothermalL.argtypes = [ c_void_p, c_int ]
c_lib.darkMatterProfileDMOMultipleL.restype  = c_void_p
c_lib.darkMatterProfileDMOMultipleL.argtypes = [ c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMOTruncatedL.restype  = c_void_p
c_lib.darkMatterProfileDMOTruncatedL.argtypes = [ c_double, c_double, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMOTruncatedExponentialL.restype  = c_void_p
c_lib.darkMatterProfileDMOTruncatedExponentialL.argtypes = [ c_double, c_double, c_double, c_double, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.darkMatterProfileDMOCircularVelocityL.restype  = c_double
c_lib.darkMatterProfileDMOCircularVelocityL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMOCircularVelocityMaximumL.restype  = c_double
c_lib.darkMatterProfileDMOCircularVelocityMaximumL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterProfileDMODensityL.restype  = c_double
c_lib.darkMatterProfileDMODensityL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMODensityLogSlopeL.restype  = c_double
c_lib.darkMatterProfileDMODensityLogSlopeL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMOEnclosedMassL.restype  = c_double
c_lib.darkMatterProfileDMOEnclosedMassL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMOEnergyL.restype  = c_double
c_lib.darkMatterProfileDMOEnergyL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterProfileDMOFreeFallRadiusIncreaseRateL.restype  = c_double
c_lib.darkMatterProfileDMOFreeFallRadiusIncreaseRateL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMOFreefallRadiusL.restype  = c_double
c_lib.darkMatterProfileDMOFreefallRadiusL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMOKSpaceL.restype  = c_double
c_lib.darkMatterProfileDMOKSpaceL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMOPotentialL.restype  = c_double
c_lib.darkMatterProfileDMOPotentialL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMORadialMomentL.restype  = c_double
c_lib.darkMatterProfileDMORadialMomentL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMORadialVelocityDispersionL.restype  = c_double
c_lib.darkMatterProfileDMORadialVelocityDispersionL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMORadiusCircularVelocityMaximumL.restype  = c_double
c_lib.darkMatterProfileDMORadiusCircularVelocityMaximumL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterProfileDMORadiusEnclosingDensityL.restype  = c_double
c_lib.darkMatterProfileDMORadiusEnclosingDensityL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMORadiusEnclosingMassL.restype  = c_double
c_lib.darkMatterProfileDMORadiusEnclosingMassL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMORadiusFromSpecificAngularMomentumL.restype  = c_double
c_lib.darkMatterProfileDMORadiusFromSpecificAngularMomentumL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.darkMatterProfileDMORotationNormalizationL.restype  = c_double
c_lib.darkMatterProfileDMORotationNormalizationL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.darkMatterProfileDMODestructorL.argtypes = [ c_void_p, c_int ]
c_lib.darkMatterProfileHeatingDDML.restype  = c_void_p
c_lib.darkMatterProfileHeatingDDML.argtypes = [ c_double, c_double ]
c_lib.darkMatterProfileHeatingDDMv2L.restype  = c_void_p
c_lib.darkMatterProfileHeatingDDMv2L.argtypes = [ c_void_p, c_int, c_bool, c_bool, c_double ]
c_lib.darkMatterProfileHeatingImpulsiveOutflowL.restype  = c_void_p
c_lib.darkMatterProfileHeatingImpulsiveOutflowL.argtypes = [ c_double, c_void_p ]
c_lib.darkMatterProfileHeatingNullL.restype  = c_void_p
c_lib.darkMatterProfileHeatingNullL.argtypes = [  ]
c_lib.darkMatterProfileHeatingSummationL.restype  = c_void_p
c_lib.darkMatterProfileHeatingSummationL.argtypes = [ c_void_p ]
c_lib.darkMatterProfileHeatingTidalL.restype  = c_void_p
c_lib.darkMatterProfileHeatingTidalL.argtypes = [ c_double, c_double, c_double, c_double ]
c_lib.darkMatterProfileHeatingTwoBodyRelaxationL.restype  = c_void_p
c_lib.darkMatterProfileHeatingTwoBodyRelaxationL.argtypes = [ c_double, c_double, c_double, c_double ]
c_lib.darkMatterProfileHeatingSpecificEnergyL.restype  = c_double
c_lib.darkMatterProfileHeatingSpecificEnergyL.argtypes = [ c_void_p, c_int, c_void_p, c_double, c_void_p, c_int ]
c_lib.darkMatterProfileHeatingSpecificEnergyGradientL.restype  = c_double
c_lib.darkMatterProfileHeatingSpecificEnergyGradientL.argtypes = [ c_void_p, c_int, c_void_p, c_double, c_void_p, c_int ]
c_lib.darkMatterProfileHeatingSpecificEnergyIsEverywhereZeroL.restype  = c_bool
c_lib.darkMatterProfileHeatingSpecificEnergyIsEverywhereZeroL.argtypes = [ c_void_p, c_int, c_void_p, c_void_p, c_int ]
c_lib.darkMatterProfileHeatingDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.intergalacticMediumFilteringMassGnedin2000L.restype  = c_void_p
c_lib.intergalacticMediumFilteringMassGnedin2000L.argtypes = [ c_bool, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.intergalacticMediumFilteringMassFractionBaryonsL.restype  = c_double
c_lib.intergalacticMediumFilteringMassFractionBaryonsL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.intergalacticMediumFilteringMassFractionBaryonsGradientMassL.restype  = c_double
c_lib.intergalacticMediumFilteringMassFractionBaryonsGradientMassL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.intergalacticMediumFilteringMassFractionBaryonsRateOfChangeL.restype  = c_double
c_lib.intergalacticMediumFilteringMassFractionBaryonsRateOfChangeL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.intergalacticMediumFilteringMassMassFilteringL.restype  = c_double
c_lib.intergalacticMediumFilteringMassMassFilteringL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumFilteringMassMassFilteringRateOfChangeL.restype  = c_double
c_lib.intergalacticMediumFilteringMassMassFilteringRateOfChangeL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumFilteringMassDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.intergalacticMediumStateRecFastL.restype  = c_void_p
c_lib.intergalacticMediumStateRecFastL.argtypes = [ c_void_p, c_int, c_void_p, c_int ]
c_lib.intergalacticMediumStateFileL.restype  = c_void_p
c_lib.intergalacticMediumStateFileL.argtypes = [ c_char_p, c_void_p, c_int, c_void_p, c_int ]
c_lib.intergalacticMediumStateInstantReionizationL.restype  = c_void_p
c_lib.intergalacticMediumStateInstantReionizationL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_double, c_double ]
c_lib.intergalacticMediumStateInternalL.restype  = c_void_p
c_lib.intergalacticMediumStateInternalL.argtypes = [ c_void_p, c_int, c_void_p, c_int ]
c_lib.intergalacticMediumStateSimpleL.restype  = c_void_p
c_lib.intergalacticMediumStateSimpleL.argtypes = [ c_double, c_double, c_double, c_void_p, c_int, c_void_p, c_int ]
c_lib.intergalacticMediumStateDoublyIonizedHeliumFractionL.restype  = c_double
c_lib.intergalacticMediumStateDoublyIonizedHeliumFractionL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateElectronFractionL.restype  = c_double
c_lib.intergalacticMediumStateElectronFractionL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateElectronScatteringOpticalDepthL.restype  = c_double
c_lib.intergalacticMediumStateElectronScatteringOpticalDepthL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateElectronScatteringTimeL.restype  = c_double
c_lib.intergalacticMediumStateElectronScatteringTimeL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateMassJeansL.restype  = c_double
c_lib.intergalacticMediumStateMassJeansL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateNeutralHeliumFractionL.restype  = c_double
c_lib.intergalacticMediumStateNeutralHeliumFractionL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateNeutralHydrogenFractionL.restype  = c_double
c_lib.intergalacticMediumStateNeutralHydrogenFractionL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateSinglyIonizedHeliumFractionL.restype  = c_double
c_lib.intergalacticMediumStateSinglyIonizedHeliumFractionL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateSinglyIonizedHydrogenFractionL.restype  = c_double
c_lib.intergalacticMediumStateSinglyIonizedHydrogenFractionL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateTemperatureL.restype  = c_double
c_lib.intergalacticMediumStateTemperatureL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.intergalacticMediumStateDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.nbodyHaloMassErrorSOHaloFinderL.restype  = c_void_p
c_lib.nbodyHaloMassErrorSOHaloFinderL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_double ]
c_lib.nbodyHaloMassErrorTrenti2010L.restype  = c_void_p
c_lib.nbodyHaloMassErrorTrenti2010L.argtypes = [ c_double ]
c_lib.nbodyHaloMassErrorFriendsOfFriendsL.restype  = c_void_p
c_lib.nbodyHaloMassErrorFriendsOfFriendsL.argtypes = [ c_double ]
c_lib.nbodyHaloMassErrorNullL.restype  = c_void_p
c_lib.nbodyHaloMassErrorNullL.argtypes = [  ]
c_lib.nbodyHaloMassErrorPowerLawL.restype  = c_void_p
c_lib.nbodyHaloMassErrorPowerLawL.argtypes = [ c_double, c_double, c_double ]
c_lib.nbodyHaloMassErrorCorrelationL.restype  = c_double
c_lib.nbodyHaloMassErrorCorrelationL.argtypes = [ c_void_p, c_int, c_void_p, c_void_p ]
c_lib.nbodyHaloMassErrorErrorFractionalL.restype  = c_double
c_lib.nbodyHaloMassErrorErrorFractionalL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.nbodyHaloMassErrorErrorZeroAlwaysL.restype  = c_bool
c_lib.nbodyHaloMassErrorErrorZeroAlwaysL.argtypes = [ c_void_p, c_int ]
c_lib.nbodyHaloMassErrorDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.criticalOverdensityKitayamaSuto1996L.restype  = c_void_p
c_lib.criticalOverdensityKitayamaSuto1996L.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.criticalOverdensityMarsh2016FDML.restype  = c_void_p
c_lib.criticalOverdensityMarsh2016FDML.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_bool ]
c_lib.criticalOverdensityEnvironmentalL.restype  = c_void_p
c_lib.criticalOverdensityEnvironmentalL.argtypes = [ c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.criticalOverdensityFixedL.restype  = c_void_p
c_lib.criticalOverdensityFixedL.argtypes = [ c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.criticalOverdensityPeakBackgroundSplitL.restype  = c_void_p
c_lib.criticalOverdensityPeakBackgroundSplitL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.criticalOverdensityRenormalizeL.restype  = c_void_p
c_lib.criticalOverdensityRenormalizeL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgyL.restype  = c_void_p
c_lib.criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgyL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_bool ]
c_lib.criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstntL.restype  = c_void_p
c_lib.criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstntL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_bool ]
c_lib.criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgyL.restype  = c_void_p
c_lib.criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgyL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_bool ]
c_lib.criticalOverdensityBarkana2001WDML.restype  = c_void_p
c_lib.criticalOverdensityBarkana2001WDML.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_bool ]
c_lib.criticalOverdensityCollapsingMassL.restype  = c_double
c_lib.criticalOverdensityCollapsingMassL.argtypes = [ c_void_p, c_int ]
c_lib.criticalOverdensityGradientMassL.restype  = c_double
c_lib.criticalOverdensityGradientMassL.argtypes = [ c_void_p, c_int ]
c_lib.criticalOverdensityGradientTimeL.restype  = c_double
c_lib.criticalOverdensityGradientTimeL.argtypes = [ c_void_p, c_int ]
c_lib.criticalOverdensityIsMassDependentL.restype  = c_bool
c_lib.criticalOverdensityIsMassDependentL.argtypes = [ c_void_p, c_int ]
c_lib.criticalOverdensityIsNodeDependentL.restype  = c_bool
c_lib.criticalOverdensityIsNodeDependentL.argtypes = [ c_void_p, c_int ]
c_lib.criticalOverdensityTimeOfCollapseL.restype  = c_double
c_lib.criticalOverdensityTimeOfCollapseL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.criticalOverdensityValueL.restype  = c_double
c_lib.criticalOverdensityValueL.argtypes = [ c_void_p, c_int ]
c_lib.criticalOverdensityDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.haloEnvironmentFixedL.restype  = c_void_p
c_lib.haloEnvironmentFixedL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_double ]
c_lib.haloEnvironmentLogNormalL.restype  = c_void_p
c_lib.haloEnvironmentLogNormalL.argtypes = [ c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloEnvironmentNormalL.restype  = c_void_p
c_lib.haloEnvironmentNormalL.argtypes = [ c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloEnvironmentUniformL.restype  = c_void_p
c_lib.haloEnvironmentUniformL.argtypes = [  ]
c_lib.haloEnvironmentCdfL.restype  = c_double
c_lib.haloEnvironmentCdfL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.haloEnvironmentEnvironmentMassL.restype  = c_double
c_lib.haloEnvironmentEnvironmentMassL.argtypes = [ c_void_p, c_int ]
c_lib.haloEnvironmentEnvironmentRadiusL.restype  = c_double
c_lib.haloEnvironmentEnvironmentRadiusL.argtypes = [ c_void_p, c_int ]
c_lib.haloEnvironmentOverdensityIsSettableL.restype  = c_bool
c_lib.haloEnvironmentOverdensityIsSettableL.argtypes = [ c_void_p, c_int ]
c_lib.haloEnvironmentOverdensityLinearL.restype  = c_double
c_lib.haloEnvironmentOverdensityLinearL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.haloEnvironmentOverdensityLinearGradientTimeL.restype  = c_double
c_lib.haloEnvironmentOverdensityLinearGradientTimeL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.haloEnvironmentOverdensityLinearMaximumL.restype  = c_double
c_lib.haloEnvironmentOverdensityLinearMaximumL.argtypes = [ c_void_p, c_int ]
c_lib.haloEnvironmentOverdensityLinearMinimumL.restype  = c_double
c_lib.haloEnvironmentOverdensityLinearMinimumL.argtypes = [ c_void_p, c_int ]
c_lib.haloEnvironmentOverdensityLinearSetL.argtypes = [ c_void_p, c_int, c_void_p, c_double ]
c_lib.haloEnvironmentOverdensityNonLinearL.restype  = c_double
c_lib.haloEnvironmentOverdensityNonLinearL.argtypes = [ c_void_p, c_int, c_void_p ]
c_lib.haloEnvironmentPdfL.restype  = c_double
c_lib.haloEnvironmentPdfL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.haloEnvironmentVolumeFractionOccupiedL.restype  = c_double
c_lib.haloEnvironmentVolumeFractionOccupiedL.argtypes = [ c_void_p, c_int ]
c_lib.haloEnvironmentDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologicalMassVarianceFilteredPowerL.restype  = c_void_p
c_lib.cosmologicalMassVarianceFilteredPowerL.argtypes = [  ]
c_lib.cosmologicalMassVariancePeakBackgroundSplitL.restype  = c_void_p
c_lib.cosmologicalMassVariancePeakBackgroundSplitL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.cosmologicalMassVarianceGrowthIsMassDependentL.restype  = c_bool
c_lib.cosmologicalMassVarianceGrowthIsMassDependentL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologicalMassVarianceMassL.restype  = c_double
c_lib.cosmologicalMassVarianceMassL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.cosmologicalMassVariancePowerNormalizationL.restype  = c_double
c_lib.cosmologicalMassVariancePowerNormalizationL.argtypes = [ c_void_p, c_int ]
c_lib.cosmologicalMassVarianceRootVarianceL.restype  = c_double
c_lib.cosmologicalMassVarianceRootVarianceL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.cosmologicalMassVarianceRootVarianceAndLogarithmicGradientL.argtypes = [ c_void_p, c_int, c_double, c_double, POINTER(c_double), POINTER(c_double) ]
c_lib.cosmologicalMassVarianceRootVarianceLogarithmicGradientL.restype  = c_double
c_lib.cosmologicalMassVarianceRootVarianceLogarithmicGradientL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.cosmologicalMassVarianceRootVarianceLogarithmicGradientTimeL.restype  = c_double
c_lib.cosmologicalMassVarianceRootVarianceLogarithmicGradientTimeL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.cosmologicalMassVarianceSigma8L.restype  = c_double
c_lib.cosmologicalMassVarianceSigma8L.argtypes = [ c_void_p, c_int ]
c_lib.cosmologicalMassVarianceDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.excursionSetBarrierCriticalOverdensityL.restype  = c_void_p
c_lib.excursionSetBarrierCriticalOverdensityL.argtypes = [ c_void_p, c_int, c_void_p, c_int ]
c_lib.excursionSetBarrierLinearL.restype  = c_void_p
c_lib.excursionSetBarrierLinearL.argtypes = [ c_double, c_double ]
c_lib.excursionSetBarrierQuadraticL.restype  = c_void_p
c_lib.excursionSetBarrierQuadraticL.argtypes = [ c_double, c_double, c_double ]
c_lib.excursionSetBarrierRemapShethMoTormenL.restype  = c_void_p
c_lib.excursionSetBarrierRemapShethMoTormenL.argtypes = [ c_double, c_double, c_double, c_int, c_void_p, c_int ]
c_lib.excursionSetBarrierRemapScaleL.restype  = c_void_p
c_lib.excursionSetBarrierRemapScaleL.argtypes = [ c_double, c_int, c_void_p, c_int ]
c_lib.excursionSetBarrierBarrierL.restype  = c_double
c_lib.excursionSetBarrierBarrierL.argtypes = [ c_void_p, c_int, c_double, c_double, c_void_p, c_bool ]
c_lib.excursionSetBarrierBarrierGradientL.restype  = c_double
c_lib.excursionSetBarrierBarrierGradientL.argtypes = [ c_void_p, c_int, c_double, c_double, c_void_p, c_bool ]
c_lib.excursionSetBarrierDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.excursionSetFirstCrossingFarahiL.restype  = c_void_p
c_lib.excursionSetFirstCrossingFarahiL.argtypes = [ c_double, c_char_p, c_int, c_int, c_int, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.excursionSetFirstCrossingFarahiMidpointL.restype  = c_void_p
c_lib.excursionSetFirstCrossingFarahiMidpointL.argtypes = [ c_double, c_char_p, c_int, c_int, c_int, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.excursionSetFirstCrossingZhangHuiL.restype  = c_void_p
c_lib.excursionSetFirstCrossingZhangHuiL.argtypes = [ c_void_p, c_int ]
c_lib.excursionSetFirstCrossingZhangHuiHighOrderL.restype  = c_void_p
c_lib.excursionSetFirstCrossingZhangHuiHighOrderL.argtypes = [ c_void_p, c_int ]
c_lib.excursionSetFirstCrossingLinearBarrierL.restype  = c_void_p
c_lib.excursionSetFirstCrossingLinearBarrierL.argtypes = [ c_void_p, c_int, c_void_p, c_int ]
c_lib.excursionSetFirstCrossingCoordinatedMPIL.argtypes = [ c_void_p, c_int, c_bool ]
c_lib.excursionSetFirstCrossingProbabilityL.restype  = c_double
c_lib.excursionSetFirstCrossingProbabilityL.argtypes = [ c_void_p, c_int, c_double, c_double, c_void_p ]
c_lib.excursionSetFirstCrossingRateL.restype  = c_double
c_lib.excursionSetFirstCrossingRateL.argtypes = [ c_void_p, c_int, c_double, c_double, c_double, c_void_p ]
c_lib.excursionSetFirstCrossingRateNonCrossingL.restype  = c_double
c_lib.excursionSetFirstCrossingRateNonCrossingL.argtypes = [ c_void_p, c_int, c_double, c_double, c_void_p ]
c_lib.excursionSetFirstCrossingDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.haloMassFunctionBhattacharya2011L.restype  = c_void_p
c_lib.haloMassFunctionBhattacharya2011L.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_double, c_double, c_double, c_double ]
c_lib.haloMassFunctionDespali2015L.restype  = c_void_p
c_lib.haloMassFunctionDespali2015L.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionOndaroMallea2021L.restype  = c_void_p
c_lib.haloMassFunctionOndaroMallea2021L.argtypes = [ POINTER(c_double), POINTER(c_double), c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionPressSchechterL.restype  = c_void_p
c_lib.haloMassFunctionPressSchechterL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionReed2007L.restype  = c_void_p
c_lib.haloMassFunctionReed2007L.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionRodriguezPuebla2016L.restype  = c_void_p
c_lib.haloMassFunctionRodriguezPuebla2016L.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionShethTormenL.restype  = c_void_p
c_lib.haloMassFunctionShethTormenL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_double, c_double, c_double ]
c_lib.haloMassFunctionTinker2008L.restype  = c_void_p
c_lib.haloMassFunctionTinker2008L.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionTinker2008GenericL.restype  = c_void_p
c_lib.haloMassFunctionTinker2008GenericL.argtypes = [ c_double, c_double, c_double, c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionEnvironmentAveragedL.restype  = c_void_p
c_lib.haloMassFunctionEnvironmentAveragedL.argtypes = [ c_bool, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionEnvironmentalL.restype  = c_void_p
c_lib.haloMassFunctionEnvironmentalL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionErrorConvolvedL.restype  = c_void_p
c_lib.haloMassFunctionErrorConvolvedL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_double ]
c_lib.haloMassFunctionFofBiasL.restype  = c_void_p
c_lib.haloMassFunctionFofBiasL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_double, c_double, c_bool, c_double ]
c_lib.haloMassFunctionSimpleSystematicL.restype  = c_void_p
c_lib.haloMassFunctionSimpleSystematicL.argtypes = [ c_double, c_double, c_void_p, c_int, c_void_p, c_int ]
c_lib.haloMassFunctionDifferentialL.restype  = c_double
c_lib.haloMassFunctionDifferentialL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.haloMassFunctionIntegratedL.restype  = c_double
c_lib.haloMassFunctionIntegratedL.argtypes = [ c_void_p, c_int, c_double, c_double, c_double ]
c_lib.haloMassFunctionMassFractionL.restype  = c_double
c_lib.haloMassFunctionMassFractionL.argtypes = [ c_void_p, c_int, c_double, c_double, c_double ]
c_lib.haloMassFunctionDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.linearGrowthBaryonsDarkMatterL.restype  = c_void_p
c_lib.linearGrowthBaryonsDarkMatterL.argtypes = [ c_double, c_double, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.linearGrowthCollisionlessMatterL.restype  = c_void_p
c_lib.linearGrowthCollisionlessMatterL.argtypes = [ c_void_p, c_int, c_void_p, c_int ]
c_lib.linearGrowthNonClusteringBaryonsDarkMatterL.restype  = c_void_p
c_lib.linearGrowthNonClusteringBaryonsDarkMatterL.argtypes = [ c_void_p, c_int, c_void_p, c_int ]
c_lib.linearGrowthIsWavenumberDependentL.restype  = c_bool
c_lib.linearGrowthIsWavenumberDependentL.argtypes = [ c_void_p, c_int ]
c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL.restype  = c_double
c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL.argtypes = [ c_void_p, c_int ]
c_lib.linearGrowthLogarithmicDerivativeWavenumberL.restype  = c_double
c_lib.linearGrowthLogarithmicDerivativeWavenumberL.argtypes = [ c_void_p, c_int ]
c_lib.linearGrowthValueL.restype  = c_double
c_lib.linearGrowthValueL.argtypes = [ c_void_p, c_int ]
c_lib.linearGrowthDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.powerSpectrumPrimordialCosmologicalCubeL.restype  = c_void_p
c_lib.powerSpectrumPrimordialCosmologicalCubeL.argtypes = [ c_double, c_double, c_void_p, c_int ]
c_lib.powerSpectrumPrimordialPiecewisePowerLawL.restype  = c_void_p
c_lib.powerSpectrumPrimordialPiecewisePowerLawL.argtypes = [ POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double) ]
c_lib.powerSpectrumPrimordialPowerLawL.restype  = c_void_p
c_lib.powerSpectrumPrimordialPowerLawL.argtypes = [ c_double, c_double, c_double, c_double, c_bool ]
c_lib.powerSpectrumPrimordialLogarithmicDerivativeL.restype  = c_double
c_lib.powerSpectrumPrimordialLogarithmicDerivativeL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.powerSpectrumPrimordialPowerL.restype  = c_double
c_lib.powerSpectrumPrimordialPowerL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.powerSpectrumPrimordialDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.powerSpectrumPrimordialTransferredFileL.restype  = c_void_p
c_lib.powerSpectrumPrimordialTransferredFileL.argtypes = [ c_char_p, c_void_p, c_int, c_void_p, c_int ]
c_lib.powerSpectrumPrimordialTransferredSimpleL.restype  = c_void_p
c_lib.powerSpectrumPrimordialTransferredSimpleL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.powerSpectrumPrimordialTransferredGrowthIsWavenumberDependentL.restype  = c_bool
c_lib.powerSpectrumPrimordialTransferredGrowthIsWavenumberDependentL.argtypes = [ c_void_p, c_int ]
c_lib.powerSpectrumPrimordialTransferredLogarithmicDerivativeL.restype  = c_double
c_lib.powerSpectrumPrimordialTransferredLogarithmicDerivativeL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.powerSpectrumPrimordialTransferredPowerL.restype  = c_double
c_lib.powerSpectrumPrimordialTransferredPowerL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.powerSpectrumPrimordialTransferredDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.powerSpectrumWindowFunctionETHOSL.restype  = c_void_p
c_lib.powerSpectrumWindowFunctionETHOSL.argtypes = [ c_double, c_double, c_void_p, c_int ]
c_lib.powerSpectrumWindowFunctionLagrangianChan2017L.restype  = c_void_p
c_lib.powerSpectrumWindowFunctionLagrangianChan2017L.argtypes = [ c_void_p, c_int, c_double ]
c_lib.powerSpectrumWindowFunctionSharpKSpaceL.restype  = c_void_p
c_lib.powerSpectrumWindowFunctionSharpKSpaceL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.powerSpectrumWindowFunctionSmoothKSpaceL.restype  = c_void_p
c_lib.powerSpectrumWindowFunctionSmoothKSpaceL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.powerSpectrumWindowFunctionTopHatL.restype  = c_void_p
c_lib.powerSpectrumWindowFunctionTopHatL.argtypes = [ c_void_p, c_int ]
c_lib.powerSpectrumWindowFunctionTopHatGeneralizedL.restype  = c_void_p
c_lib.powerSpectrumWindowFunctionTopHatGeneralizedL.argtypes = [ c_double, c_void_p, c_int ]
c_lib.powerSpectrumWindowFunctionTopHatSharpKHybridL.restype  = c_void_p
c_lib.powerSpectrumWindowFunctionTopHatSharpKHybridL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.powerSpectrumWindowFunctionTopHatSmoothedL.restype  = c_void_p
c_lib.powerSpectrumWindowFunctionTopHatSmoothedL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.powerSpectrumWindowFunctionAmplitudeIsMassIndependentL.restype  = c_bool
c_lib.powerSpectrumWindowFunctionAmplitudeIsMassIndependentL.argtypes = [ c_void_p, c_int ]
c_lib.powerSpectrumWindowFunctionValueL.restype  = c_double
c_lib.powerSpectrumWindowFunctionValueL.argtypes = [ c_void_p, c_int, c_double, c_double ]
c_lib.powerSpectrumWindowFunctionWavenumberMaximumL.restype  = c_double
c_lib.powerSpectrumWindowFunctionWavenumberMaximumL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.powerSpectrumWindowFunctionDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.transferFunctionAxionCAMBL.restype  = c_void_p
c_lib.transferFunctionAxionCAMBL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_double, c_int ]
c_lib.transferFunctionBBKSL.restype  = c_void_p
c_lib.transferFunctionBBKSL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionBBKSWDML.restype  = c_void_p
c_lib.transferFunctionBBKSWDML.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionBode2001L.restype  = c_void_p
c_lib.transferFunctionBode2001L.argtypes = [ c_void_p, c_int, c_double, c_double, c_double, c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionCAMBL.restype  = c_void_p
c_lib.transferFunctionCAMBL.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_double, c_int ]
c_lib.transferFunctionCLASSCDML.restype  = c_void_p
c_lib.transferFunctionCLASSCDML.argtypes = [ c_void_p, c_int, c_void_p, c_int, c_void_p, c_int, c_double, c_int ]
c_lib.transferFunctionETHOSDML.restype  = c_void_p
c_lib.transferFunctionETHOSDML.argtypes = [ c_void_p, c_int, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionEisensteinHu1999L.restype  = c_void_p
c_lib.transferFunctionEisensteinHu1999L.argtypes = [ c_double, c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionHu2000FDML.restype  = c_void_p
c_lib.transferFunctionHu2000FDML.argtypes = [ c_void_p, c_int, c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionMurgia2017L.restype  = c_void_p
c_lib.transferFunctionMurgia2017L.argtypes = [ c_void_p, c_int, c_double, c_double, c_double, c_double, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionAcceleratorL.restype  = c_void_p
c_lib.transferFunctionAcceleratorL.argtypes = [ c_void_p, c_int, c_int ]
c_lib.transferFunctionFileL.restype  = c_void_p
c_lib.transferFunctionFileL.argtypes = [ c_char_p, c_double, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionFileFuzzyDarkMatterL.restype  = c_void_p
c_lib.transferFunctionFileFuzzyDarkMatterL.argtypes = [ c_char_p, c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionFuzzyDML.restype  = c_void_p
c_lib.transferFunctionFuzzyDML.argtypes = [ c_void_p, c_int, c_double, c_double, c_double, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.transferFunctionIdentityL.restype  = c_void_p
c_lib.transferFunctionIdentityL.argtypes = [ c_double ]
c_lib.transferFunctionEpochTimeL.restype  = c_double
c_lib.transferFunctionEpochTimeL.argtypes = [ c_void_p, c_int ]
c_lib.transferFunctionFractionModeMassL.restype  = c_double
c_lib.transferFunctionFractionModeMassL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.transferFunctionHalfModeMassL.restype  = c_double
c_lib.transferFunctionHalfModeMassL.argtypes = [ c_void_p, c_int ]
c_lib.transferFunctionLogarithmicDerivativeL.restype  = c_double
c_lib.transferFunctionLogarithmicDerivativeL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.transferFunctionQuarterModeMassL.restype  = c_double
c_lib.transferFunctionQuarterModeMassL.argtypes = [ c_void_p, c_int ]
c_lib.transferFunctionValueL.restype  = c_double
c_lib.transferFunctionValueL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.transferFunctionDestructorL.argtypes = [ c_void_p, c_int ]
c_lib.virialDensityContrastBryanNorman1998L.restype  = c_void_p
c_lib.virialDensityContrastBryanNorman1998L.argtypes = [ c_void_p, c_int, c_void_p, c_int ]
c_lib.virialDensityContrastKitayamaSuto1996L.restype  = c_void_p
c_lib.virialDensityContrastKitayamaSuto1996L.argtypes = [ c_void_p, c_int ]
c_lib.virialDensityContrastFixedL.restype  = c_void_p
c_lib.virialDensityContrastFixedL.argtypes = [ c_double, c_int, c_double, c_void_p, c_int, c_void_p, c_int ]
c_lib.virialDensityContrastFriendsOfFriendsL.restype  = c_void_p
c_lib.virialDensityContrastFriendsOfFriendsL.argtypes = [ c_double, c_double ]
c_lib.virialDensityContrastPercolationL.restype  = c_void_p
c_lib.virialDensityContrastPercolationL.argtypes = [ c_double, c_void_p, c_int, c_void_p ]
c_lib.virialDensityContrastSphericalCollapseBrynsDrkMttrDrkEnrgyL.restype  = c_void_p
c_lib.virialDensityContrastSphericalCollapseBrynsDrkMttrDrkEnrgyL.argtypes = [ c_bool, c_int, c_void_p, c_int, c_void_p, c_int, c_void_p, c_int ]
c_lib.virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstntL.restype  = c_void_p
c_lib.virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstntL.argtypes = [ c_bool, c_void_p, c_int ]
c_lib.virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgyL.restype  = c_void_p
c_lib.virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgyL.argtypes = [ c_bool, c_int, c_void_p, c_int ]
c_lib.virialDensityContrastDensityContrastL.restype  = c_double
c_lib.virialDensityContrastDensityContrastL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.virialDensityContrastDensityContrastRateOfChangeL.restype  = c_double
c_lib.virialDensityContrastDensityContrastRateOfChangeL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.virialDensityContrastIsMassDependentL.restype  = c_bool
c_lib.virialDensityContrastIsMassDependentL.argtypes = [ c_void_p, c_int ]
c_lib.virialDensityContrastTurnAroundOverVirialRadiiL.restype  = c_double
c_lib.virialDensityContrastTurnAroundOverVirialRadiiL.argtypes = [ c_void_p, c_int, c_double ]
c_lib.virialDensityContrastDestructorL.argtypes = [ c_void_p, c_int ]

class cosmologyFunctions:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.cosmologyFunctionsDestructorL(self._glcObj,self._classID)

    def timeBigCrunch(self):
        return c_lib.cosmologyFunctionsTimeBigCrunchL(self._glcObj,self._classID)

    def timeAtDistanceComoving(self,comovingDistance):
        return c_lib.cosmologyFunctionsTimeAtDistanceComovingL(self._glcObj,self._classID,comovingDistance)

    def temperatureCMBEpochal(self,time=None,expansionFactor=None,collapsingPhase=None):
        if not collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsTemperatureCMBEpochalL(self._glcObj,self._classID,None,None,None)
        elif not collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsTemperatureCMBEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,None)
        elif not collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsTemperatureCMBEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None)
        elif not collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsTemperatureCMBEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsTemperatureCMBEpochalL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsTemperatureCMBEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsTemperatureCMBEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsTemperatureCMBEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
    

    def redshiftFromExpansionFactor(self,expansionFactor):
        return c_lib.cosmologyFunctionsRedshiftFromExpansionFactorL(self._glcObj,self._classID,expansionFactor)

    def omegaMatterRateOfChange(self,time=None,expansionFactor=None,collapsingPhase=None):
        if not collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL(self._glcObj,self._classID,None,None,None)
        elif not collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL(self._glcObj,self._classID,byref(c_double(time)),None,None)
        elif not collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None)
        elif not collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaMatterRateOfChangeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
    

    def omegaMatterEpochal(self,time=None,expansionFactor=None,collapsingPhase=None):
        if not collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaMatterEpochalL(self._glcObj,self._classID,None,None,None)
        elif not collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaMatterEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,None)
        elif not collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaMatterEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None)
        elif not collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaMatterEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaMatterEpochalL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaMatterEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaMatterEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaMatterEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
    

    def omegaDarkEnergyEpochal(self,time=None,expansionFactor=None,collapsingPhase=None):
        if not collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL(self._glcObj,self._classID,None,None,None)
        elif not collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,None)
        elif not collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None)
        elif not collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsOmegaDarkEnergyEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
    

    def matterDensityEpochal(self,time=None,expansionFactor=None,collapsingPhase=None):
        if not collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsMatterDensityEpochalL(self._glcObj,self._classID,None,None,None)
        elif not collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsMatterDensityEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,None)
        elif not collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsMatterDensityEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None)
        elif not collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsMatterDensityEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsMatterDensityEpochalL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsMatterDensityEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsMatterDensityEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsMatterDensityEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
    

    def hubbleParameterRateOfChange(self,time=None,expansionFactor=None,collapsingPhase=None):
        if not collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL(self._glcObj,self._classID,None,None,None)
        elif not collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL(self._glcObj,self._classID,byref(c_double(time)),None,None)
        elif not collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None)
        elif not collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsHubbleParameterRateOfChangeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
    

    def hubbleParameterEpochal(self,time=None,expansionFactor=None,collapsingPhase=None):
        if not collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsHubbleParameterEpochalL(self._glcObj,self._classID,None,None,None)
        elif not collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsHubbleParameterEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,None)
        elif not collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsHubbleParameterEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None)
        elif not collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsHubbleParameterEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsHubbleParameterEpochalL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsHubbleParameterEpochalL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsHubbleParameterEpochalL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsHubbleParameterEpochalL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
    

    def exponentDarkEnergy(self,time=None,expansionFactor=None):
        if not expansionFactor and not time:
            return c_lib.cosmologyFunctionsExponentDarkEnergyL(self._glcObj,self._classID,None,None)
        elif not expansionFactor and time:
            return c_lib.cosmologyFunctionsExponentDarkEnergyL(self._glcObj,self._classID,byref(c_double(time)),None)
        elif expansionFactor and not time:
            return c_lib.cosmologyFunctionsExponentDarkEnergyL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)))
        elif expansionFactor and time:
            return c_lib.cosmologyFunctionsExponentDarkEnergyL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)))
    

    def expansionRate(self,expansionFactor):
        return c_lib.cosmologyFunctionsExpansionRateL(self._glcObj,self._classID,expansionFactor)

    def expansionFactorFromRedshift(self,redshift):
        return c_lib.cosmologyFunctionsExpansionFactorFromRedshiftL(self._glcObj,self._classID,redshift)

    def expansionFactor(self,time):
        return c_lib.cosmologyFunctionsExpansionFactorL(self._glcObj,self._classID,time)

    def equationOfStateDarkEnergy(self,time=None,expansionFactor=None):
        if not expansionFactor and not time:
            return c_lib.cosmologyFunctionsEquationOfStateDarkEnergyL(self._glcObj,self._classID,None,None)
        elif not expansionFactor and time:
            return c_lib.cosmologyFunctionsEquationOfStateDarkEnergyL(self._glcObj,self._classID,byref(c_double(time)),None)
        elif expansionFactor and not time:
            return c_lib.cosmologyFunctionsEquationOfStateDarkEnergyL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)))
        elif expansionFactor and time:
            return c_lib.cosmologyFunctionsEquationOfStateDarkEnergyL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)))
    

    def equalityEpochMatterRadiation(self,requestType=None):
        if not requestType:
            return c_lib.cosmologyFunctionsEqualityEpochMatterRadiationL(self._glcObj,self._classID,None)
        elif requestType:
            return c_lib.cosmologyFunctionsEqualityEpochMatterRadiationL(self._glcObj,self._classID,byref(c_int(requestType)))
    

    def equalityEpochMatterDarkEnergy(self,requestType=None):
        if not requestType:
            return c_lib.cosmologyFunctionsEqualityEpochMatterDarkEnergyL(self._glcObj,self._classID,None)
        elif requestType:
            return c_lib.cosmologyFunctionsEqualityEpochMatterDarkEnergyL(self._glcObj,self._classID,byref(c_int(requestType)))
    

    def equalityEpochMatterCurvature(self,requestType=None):
        if not requestType:
            return c_lib.cosmologyFunctionsEqualityEpochMatterCurvatureL(self._glcObj,self._classID,None)
        elif requestType:
            return c_lib.cosmologyFunctionsEqualityEpochMatterCurvatureL(self._glcObj,self._classID,byref(c_int(requestType)))
    

    def epochValidate(self,timeIn=None,expansionFactorIn=None,collapsingIn=None,timeOut=None,expansionFactorOut=None,collapsingOut=None):
        if not collapsingIn and not collapsingOut and not expansionFactorIn and not expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,None,None,None,None)
        elif not collapsingIn and not collapsingOut and not expansionFactorIn and not expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,None,byref(c_double(timeOut)),None,None)
        elif not collapsingIn and not collapsingOut and not expansionFactorIn and not expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,None,None,None,None)
        elif not collapsingIn and not collapsingOut and not expansionFactorIn and not expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,None,byref(c_double(timeOut)),None,None)
        elif not collapsingIn and not collapsingOut and not expansionFactorIn and expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,None,None,byref(c_double(expansionFactorOut)),None)
        elif not collapsingIn and not collapsingOut and not expansionFactorIn and expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,None,byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),None)
        elif not collapsingIn and not collapsingOut and not expansionFactorIn and expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,None,None,byref(c_double(expansionFactorOut)),None)
        elif not collapsingIn and not collapsingOut and not expansionFactorIn and expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,None,byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),None)
        elif not collapsingIn and not collapsingOut and expansionFactorIn and not expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),None,None,None,None)
        elif not collapsingIn and not collapsingOut and expansionFactorIn and not expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),None,byref(c_double(timeOut)),None,None)
        elif not collapsingIn and not collapsingOut and expansionFactorIn and not expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),None,None,None,None)
        elif not collapsingIn and not collapsingOut and expansionFactorIn and not expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),None,byref(c_double(timeOut)),None,None)
        elif not collapsingIn and not collapsingOut and expansionFactorIn and expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),None,None,byref(c_double(expansionFactorOut)),None)
        elif not collapsingIn and not collapsingOut and expansionFactorIn and expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),None,byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),None)
        elif not collapsingIn and not collapsingOut and expansionFactorIn and expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),None,None,byref(c_double(expansionFactorOut)),None)
        elif not collapsingIn and not collapsingOut and expansionFactorIn and expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),None,byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),None)
        elif not collapsingIn and collapsingOut and not expansionFactorIn and not expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,None,None,None,byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and not expansionFactorIn and not expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,None,byref(c_double(timeOut)),None,byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and not expansionFactorIn and not expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,None,None,None,byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and not expansionFactorIn and not expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,None,byref(c_double(timeOut)),None,byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and not expansionFactorIn and expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,None,None,byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and not expansionFactorIn and expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,None,byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and not expansionFactorIn and expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,None,None,byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and not expansionFactorIn and expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,None,byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and expansionFactorIn and not expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),None,None,None,byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and expansionFactorIn and not expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),None,byref(c_double(timeOut)),None,byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and expansionFactorIn and not expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),None,None,None,byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and expansionFactorIn and not expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),None,byref(c_double(timeOut)),None,byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and expansionFactorIn and expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),None,None,byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and expansionFactorIn and expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),None,byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and expansionFactorIn and expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),None,None,byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif not collapsingIn and collapsingOut and expansionFactorIn and expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),None,byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif collapsingIn and not collapsingOut and not expansionFactorIn and not expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingIn)),None,None,None)
        elif collapsingIn and not collapsingOut and not expansionFactorIn and not expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingIn)),byref(c_double(timeOut)),None,None)
        elif collapsingIn and not collapsingOut and not expansionFactorIn and not expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,byref(c_bool(collapsingIn)),None,None,None)
        elif collapsingIn and not collapsingOut and not expansionFactorIn and not expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,byref(c_bool(collapsingIn)),byref(c_double(timeOut)),None,None)
        elif collapsingIn and not collapsingOut and not expansionFactorIn and expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingIn)),None,byref(c_double(expansionFactorOut)),None)
        elif collapsingIn and not collapsingOut and not expansionFactorIn and expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingIn)),byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),None)
        elif collapsingIn and not collapsingOut and not expansionFactorIn and expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,byref(c_bool(collapsingIn)),None,byref(c_double(expansionFactorOut)),None)
        elif collapsingIn and not collapsingOut and not expansionFactorIn and expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,byref(c_bool(collapsingIn)),byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),None)
        elif collapsingIn and not collapsingOut and expansionFactorIn and not expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),None,None,None)
        elif collapsingIn and not collapsingOut and expansionFactorIn and not expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),byref(c_double(timeOut)),None,None)
        elif collapsingIn and not collapsingOut and expansionFactorIn and not expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),None,None,None)
        elif collapsingIn and not collapsingOut and expansionFactorIn and not expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),byref(c_double(timeOut)),None,None)
        elif collapsingIn and not collapsingOut and expansionFactorIn and expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),None,byref(c_double(expansionFactorOut)),None)
        elif collapsingIn and not collapsingOut and expansionFactorIn and expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),None)
        elif collapsingIn and not collapsingOut and expansionFactorIn and expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),None,byref(c_double(expansionFactorOut)),None)
        elif collapsingIn and not collapsingOut and expansionFactorIn and expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),None)
        elif collapsingIn and collapsingOut and not expansionFactorIn and not expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingIn)),None,None,byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and not expansionFactorIn and not expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingIn)),byref(c_double(timeOut)),None,byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and not expansionFactorIn and not expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,byref(c_bool(collapsingIn)),None,None,byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and not expansionFactorIn and not expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,byref(c_bool(collapsingIn)),byref(c_double(timeOut)),None,byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and not expansionFactorIn and expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingIn)),None,byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and not expansionFactorIn and expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingIn)),byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and not expansionFactorIn and expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,byref(c_bool(collapsingIn)),None,byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and not expansionFactorIn and expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),None,byref(c_bool(collapsingIn)),byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and expansionFactorIn and not expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),None,None,byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and expansionFactorIn and not expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),byref(c_double(timeOut)),None,byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and expansionFactorIn and not expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),None,None,byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and expansionFactorIn and not expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),byref(c_double(timeOut)),None,byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and expansionFactorIn and expansionFactorOut and not timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),None,byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and expansionFactorIn and expansionFactorOut and not timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,None,byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and expansionFactorIn and expansionFactorOut and timeIn and not timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),None,byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
        elif collapsingIn and collapsingOut and expansionFactorIn and expansionFactorOut and timeIn and timeOut:
            return c_lib.cosmologyFunctionsEpochValidateL(self._glcObj,self._classID,byref(c_double(timeIn)),byref(c_double(expansionFactorIn)),byref(c_bool(collapsingIn)),byref(c_double(timeOut)),byref(c_double(expansionFactorOut)),byref(c_bool(collapsingOut)))
    

    def epochTime(self,time=None,expansionFactor=None,collapsingPhase=None):
        if not collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsEpochTimeL(self._glcObj,self._classID,None,None,None)
        elif not collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsEpochTimeL(self._glcObj,self._classID,byref(c_double(time)),None,None)
        elif not collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsEpochTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None)
        elif not collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsEpochTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsingPhase and not expansionFactor and not time:
            return c_lib.cosmologyFunctionsEpochTimeL(self._glcObj,self._classID,None,None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and not expansionFactor and time:
            return c_lib.cosmologyFunctionsEpochTimeL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and not time:
            return c_lib.cosmologyFunctionsEpochTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
        elif collapsingPhase and expansionFactor and time:
            return c_lib.cosmologyFunctionsEpochTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsingPhase)))
    

    def dominationEpochMatter(self,dominateFactor):
        return c_lib.cosmologyFunctionsDominationEpochMatterL(self._glcObj,self._classID,dominateFactor)

    def distanceParticleHorizonComoving(self,time):
        return c_lib.cosmologyFunctionsDistanceParticleHorizonComovingL(self._glcObj,self._classID,time)

    def distanceLuminosity(self,time):
        return c_lib.cosmologyFunctionsDistanceLuminosityL(self._glcObj,self._classID,time)

    def distanceComovingConvert(self,output,distanceLuminosity=None,distanceModulus=None,distanceModulusKCorrected=None,redshift=None):
        if not distanceLuminosity and not distanceModulus and not distanceModulusKCorrected and not redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,None,None,None,None)
        elif not distanceLuminosity and not distanceModulus and not distanceModulusKCorrected and redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,None,None,None,byref(c_double(redshift)))
        elif not distanceLuminosity and not distanceModulus and distanceModulusKCorrected and not redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,None,None,byref(c_double(distanceModulusKCorrected)),None)
        elif not distanceLuminosity and not distanceModulus and distanceModulusKCorrected and redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,None,None,byref(c_double(distanceModulusKCorrected)),byref(c_double(redshift)))
        elif not distanceLuminosity and distanceModulus and not distanceModulusKCorrected and not redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,None,byref(c_double(distanceModulus)),None,None)
        elif not distanceLuminosity and distanceModulus and not distanceModulusKCorrected and redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,None,byref(c_double(distanceModulus)),None,byref(c_double(redshift)))
        elif not distanceLuminosity and distanceModulus and distanceModulusKCorrected and not redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,None,byref(c_double(distanceModulus)),byref(c_double(distanceModulusKCorrected)),None)
        elif not distanceLuminosity and distanceModulus and distanceModulusKCorrected and redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,None,byref(c_double(distanceModulus)),byref(c_double(distanceModulusKCorrected)),byref(c_double(redshift)))
        elif distanceLuminosity and not distanceModulus and not distanceModulusKCorrected and not redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,byref(c_double(distanceLuminosity)),None,None,None)
        elif distanceLuminosity and not distanceModulus and not distanceModulusKCorrected and redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,byref(c_double(distanceLuminosity)),None,None,byref(c_double(redshift)))
        elif distanceLuminosity and not distanceModulus and distanceModulusKCorrected and not redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,byref(c_double(distanceLuminosity)),None,byref(c_double(distanceModulusKCorrected)),None)
        elif distanceLuminosity and not distanceModulus and distanceModulusKCorrected and redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,byref(c_double(distanceLuminosity)),None,byref(c_double(distanceModulusKCorrected)),byref(c_double(redshift)))
        elif distanceLuminosity and distanceModulus and not distanceModulusKCorrected and not redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,byref(c_double(distanceLuminosity)),byref(c_double(distanceModulus)),None,None)
        elif distanceLuminosity and distanceModulus and not distanceModulusKCorrected and redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,byref(c_double(distanceLuminosity)),byref(c_double(distanceModulus)),None,byref(c_double(redshift)))
        elif distanceLuminosity and distanceModulus and distanceModulusKCorrected and not redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,byref(c_double(distanceLuminosity)),byref(c_double(distanceModulus)),byref(c_double(distanceModulusKCorrected)),None)
        elif distanceLuminosity and distanceModulus and distanceModulusKCorrected and redshift:
            return c_lib.cosmologyFunctionsDistanceComovingConvertL(self._glcObj,self._classID,output,byref(c_double(distanceLuminosity)),byref(c_double(distanceModulus)),byref(c_double(distanceModulusKCorrected)),byref(c_double(redshift)))
    

    def distanceComoving(self,time):
        return c_lib.cosmologyFunctionsDistanceComovingL(self._glcObj,self._classID,time)

    def distanceAngular(self,time):
        return c_lib.cosmologyFunctionsDistanceAngularL(self._glcObj,self._classID,time)

    def densityScalingEarlyTime(self,dominateFactor,densityPower,expansionFactorDominant,OmegaDominant=None):
        if not OmegaDominant:
            return c_lib.cosmologyFunctionsDensityScalingEarlyTimeL(self._glcObj,self._classID,dominateFactor,densityPower,expansionFactorDominant,None)
        elif OmegaDominant:
            return c_lib.cosmologyFunctionsDensityScalingEarlyTimeL(self._glcObj,self._classID,dominateFactor,densityPower,expansionFactorDominant,byref(c_double(OmegaDominant)))
    

    def cosmicTime(self,expansionFactor,collapsingPhase=None):
        if not collapsingPhase:
            return c_lib.cosmologyFunctionsCosmicTimeL(self._glcObj,self._classID,expansionFactor,None)
        elif collapsingPhase:
            return c_lib.cosmologyFunctionsCosmicTimeL(self._glcObj,self._classID,expansionFactor,byref(c_bool(collapsingPhase)))
    

    def comovingVolumeElementTime(self,time):
        return c_lib.cosmologyFunctionsComovingVolumeElementTimeL(self._glcObj,self._classID,time)

    def comovingVolumeElementRedshift(self,time):
        return c_lib.cosmologyFunctionsComovingVolumeElementRedshiftL(self._glcObj,self._classID,time)

class darkMatterProfileDMO:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.darkMatterProfileDMODestructorL(self._glcObj,self._classID)

    def rotationNormalization(self,node):
        return c_lib.darkMatterProfileDMORotationNormalizationL(self._glcObj,self._classID,node)

    def radiusFromSpecificAngularMomentum(self,node,specificAngularMomentum):
        return c_lib.darkMatterProfileDMORadiusFromSpecificAngularMomentumL(self._glcObj,self._classID,node,specificAngularMomentum)

    def radiusEnclosingMass(self,node,mass):
        return c_lib.darkMatterProfileDMORadiusEnclosingMassL(self._glcObj,self._classID,node,mass)

    def radiusEnclosingDensity(self,node,density):
        return c_lib.darkMatterProfileDMORadiusEnclosingDensityL(self._glcObj,self._classID,node,density)

    def radiusCircularVelocityMaximum(self,node):
        return c_lib.darkMatterProfileDMORadiusCircularVelocityMaximumL(self._glcObj,self._classID,node)

    def radialVelocityDispersion(self,node,radius):
        return c_lib.darkMatterProfileDMORadialVelocityDispersionL(self._glcObj,self._classID,node,radius)

    def radialMoment(self,node,moment,radiusMinimum=None,radiusMaximum=None):
        if not radiusMaximum and not radiusMinimum:
            return c_lib.darkMatterProfileDMORadialMomentL(self._glcObj,self._classID,node,moment,None,None)
        elif not radiusMaximum and radiusMinimum:
            return c_lib.darkMatterProfileDMORadialMomentL(self._glcObj,self._classID,node,moment,byref(c_double(radiusMinimum)),None)
        elif radiusMaximum and not radiusMinimum:
            return c_lib.darkMatterProfileDMORadialMomentL(self._glcObj,self._classID,node,moment,None,byref(c_double(radiusMaximum)))
        elif radiusMaximum and radiusMinimum:
            return c_lib.darkMatterProfileDMORadialMomentL(self._glcObj,self._classID,node,moment,byref(c_double(radiusMinimum)),byref(c_double(radiusMaximum)))
    

    def potential(self,node,radius,status=None):
        if not status:
            return c_lib.darkMatterProfileDMOPotentialL(self._glcObj,self._classID,node,radius,None)
        elif status:
            return c_lib.darkMatterProfileDMOPotentialL(self._glcObj,self._classID,node,radius,byref(c_int(status)))
    

    def kSpace(self,node,wavenumber):
        return c_lib.darkMatterProfileDMOKSpaceL(self._glcObj,self._classID,node,wavenumber)

    def freefallRadius(self,node,time):
        return c_lib.darkMatterProfileDMOFreefallRadiusL(self._glcObj,self._classID,node,time)

    def freeFallRadiusIncreaseRate(self,node,time):
        return c_lib.darkMatterProfileDMOFreeFallRadiusIncreaseRateL(self._glcObj,self._classID,node,time)

    def energy(self,node):
        return c_lib.darkMatterProfileDMOEnergyL(self._glcObj,self._classID,node)

    def enclosedMass(self,node,radius):
        return c_lib.darkMatterProfileDMOEnclosedMassL(self._glcObj,self._classID,node,radius)

    def densityLogSlope(self,node,radius):
        return c_lib.darkMatterProfileDMODensityLogSlopeL(self._glcObj,self._classID,node,radius)

    def density(self,node,radius):
        return c_lib.darkMatterProfileDMODensityL(self._glcObj,self._classID,node,radius)

    def circularVelocityMaximum(self,node):
        return c_lib.darkMatterProfileDMOCircularVelocityMaximumL(self._glcObj,self._classID,node)

    def circularVelocity(self,node,radius):
        return c_lib.darkMatterProfileDMOCircularVelocityL(self._glcObj,self._classID,node,radius)

class darkMatterProfileHeating:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.darkMatterProfileHeatingDestructorL(self._glcObj,self._classID)

    def specificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_):
        return c_lib.darkMatterProfileHeatingSpecificEnergyIsEverywhereZeroL(self._glcObj,self._classID,node,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID)

    def specificEnergyGradient(self,node,radius,darkMatterProfileDMO_):
        return c_lib.darkMatterProfileHeatingSpecificEnergyGradientL(self._glcObj,self._classID,node,radius,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID)

    def specificEnergy(self,node,radius,darkMatterProfileDMO_):
        return c_lib.darkMatterProfileHeatingSpecificEnergyL(self._glcObj,self._classID,node,radius,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID)

class excursionSetFirstCrossing:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.excursionSetFirstCrossingDestructorL(self._glcObj,self._classID)

    def rateNonCrossing(self,variance,time,node):
        return c_lib.excursionSetFirstCrossingRateNonCrossingL(self._glcObj,self._classID,variance,time,node)

    def rate(self,variance,varianceProgenitor,time,node):
        return c_lib.excursionSetFirstCrossingRateL(self._glcObj,self._classID,variance,varianceProgenitor,time,node)

    def probability(self,variance,time,node):
        return c_lib.excursionSetFirstCrossingProbabilityL(self._glcObj,self._classID,variance,time,node)

    def coordinatedMPI(self,state):
        return c_lib.excursionSetFirstCrossingCoordinatedMPIL(self._glcObj,self._classID,state)

class transferFunction:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.transferFunctionDestructorL(self._glcObj,self._classID)

    def value(self,wavenumber):
        return c_lib.transferFunctionValueL(self._glcObj,self._classID,wavenumber)

    def quarterModeMass(self,status=None):
        if not status:
            return c_lib.transferFunctionQuarterModeMassL(self._glcObj,self._classID,None)
        elif status:
            return c_lib.transferFunctionQuarterModeMassL(self._glcObj,self._classID,byref(c_int(status)))
    

    def logarithmicDerivative(self,wavenumber):
        return c_lib.transferFunctionLogarithmicDerivativeL(self._glcObj,self._classID,wavenumber)

    def halfModeMass(self,status=None):
        if not status:
            return c_lib.transferFunctionHalfModeMassL(self._glcObj,self._classID,None)
        elif status:
            return c_lib.transferFunctionHalfModeMassL(self._glcObj,self._classID,byref(c_int(status)))
    

    def fractionModeMass(self,fraction,status=None):
        if not status:
            return c_lib.transferFunctionFractionModeMassL(self._glcObj,self._classID,fraction,None)
        elif status:
            return c_lib.transferFunctionFractionModeMassL(self._glcObj,self._classID,fraction,byref(c_int(status)))
    

    def epochTime(self):
        return c_lib.transferFunctionEpochTimeL(self._glcObj,self._classID)

class criticalOverdensity:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.criticalOverdensityDestructorL(self._glcObj,self._classID)

    def value(self,time=None,expansionFactor=None,collapsing=None,mass=None,node=None):
        if not collapsing and not expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,None,None,None,None)
        elif not collapsing and not expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,None)
        elif not collapsing and not expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,None,None,None,byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,None,None,byref(c_double(mass)),None)
        elif not collapsing and not expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_double(mass)),None)
        elif not collapsing and not expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,None,None,byref(c_double(mass)),byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_double(mass)),byref(c_void_p(node)))
        elif not collapsing and expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,byref(c_void_p(node)))
        elif not collapsing and expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,byref(c_void_p(node)))
        elif not collapsing and expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_double(mass)),None)
        elif not collapsing and expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_double(mass)),None)
        elif not collapsing and expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_double(mass)),byref(c_void_p(node)))
        elif not collapsing and expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and not expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and not expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and not expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and not expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and not expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and not expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
    

    def timeOfCollapse(self,criticalOverdensity,mass=None,node=None):
        if not mass and not node:
            return c_lib.criticalOverdensityTimeOfCollapseL(self._glcObj,self._classID,criticalOverdensity,None,None)
        elif not mass and node:
            return c_lib.criticalOverdensityTimeOfCollapseL(self._glcObj,self._classID,criticalOverdensity,None,byref(c_void_p(node)))
        elif mass and not node:
            return c_lib.criticalOverdensityTimeOfCollapseL(self._glcObj,self._classID,criticalOverdensity,byref(c_double(mass)),None)
        elif mass and node:
            return c_lib.criticalOverdensityTimeOfCollapseL(self._glcObj,self._classID,criticalOverdensity,byref(c_double(mass)),byref(c_void_p(node)))
    

    def isNodeDependent(self):
        return c_lib.criticalOverdensityIsNodeDependentL(self._glcObj,self._classID)

    def isMassDependent(self):
        return c_lib.criticalOverdensityIsMassDependentL(self._glcObj,self._classID)

    def gradientTime(self,time=None,expansionFactor=None,collapsing=None,mass=None,node=None):
        if not collapsing and not expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,None,None,None,None)
        elif not collapsing and not expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,None)
        elif not collapsing and not expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,None,None,None,byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,None,None,byref(c_double(mass)),None)
        elif not collapsing and not expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_double(mass)),None)
        elif not collapsing and not expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,None,None,byref(c_double(mass)),byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_double(mass)),byref(c_void_p(node)))
        elif not collapsing and expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,byref(c_void_p(node)))
        elif not collapsing and expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,byref(c_void_p(node)))
        elif not collapsing and expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_double(mass)),None)
        elif not collapsing and expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_double(mass)),None)
        elif not collapsing and expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_double(mass)),byref(c_void_p(node)))
        elif not collapsing and expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and not expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and not expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and not expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and not expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and not expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and not expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityGradientTimeL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
    

    def gradientMass(self,time=None,expansionFactor=None,collapsing=None,mass=None,node=None):
        if not collapsing and not expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,None,None,None,None)
        elif not collapsing and not expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,None)
        elif not collapsing and not expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,None,None,None,byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,None,None,byref(c_double(mass)),None)
        elif not collapsing and not expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_double(mass)),None)
        elif not collapsing and not expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,None,None,byref(c_double(mass)),byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_double(mass)),byref(c_void_p(node)))
        elif not collapsing and expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,byref(c_void_p(node)))
        elif not collapsing and expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,byref(c_void_p(node)))
        elif not collapsing and expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_double(mass)),None)
        elif not collapsing and expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_double(mass)),None)
        elif not collapsing and expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_double(mass)),byref(c_void_p(node)))
        elif not collapsing and expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and not expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and not expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and not expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and not expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and not expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and not expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and expansionFactor and not mass and not node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and expansionFactor and not mass and not node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and expansionFactor and not mass and node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and expansionFactor and not mass and node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_void_p(node)))
        elif collapsing and expansionFactor and mass and not node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and expansionFactor and mass and not node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),None)
        elif collapsing and expansionFactor and mass and node and not time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
        elif collapsing and expansionFactor and mass and node and time:
            return c_lib.criticalOverdensityGradientMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_double(mass)),byref(c_void_p(node)))
    

    def collapsingMass(self,time=None,expansionFactor=None,collapsing=None,node=None):
        if not collapsing and not expansionFactor and not node and not time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,None,None,None,None)
        elif not collapsing and not expansionFactor and not node and time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,byref(c_double(time)),None,None,None)
        elif not collapsing and not expansionFactor and node and not time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,None,None,None,byref(c_void_p(node)))
        elif not collapsing and not expansionFactor and node and time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_void_p(node)))
        elif not collapsing and expansionFactor and not node and not time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None)
        elif not collapsing and expansionFactor and not node and time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None)
        elif not collapsing and expansionFactor and node and not time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_void_p(node)))
        elif not collapsing and expansionFactor and node and time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_void_p(node)))
        elif collapsing and not expansionFactor and not node and not time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None)
        elif collapsing and not expansionFactor and not node and time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None)
        elif collapsing and not expansionFactor and node and not time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_void_p(node)))
        elif collapsing and not expansionFactor and node and time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_void_p(node)))
        elif collapsing and expansionFactor and not node and not time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None)
        elif collapsing and expansionFactor and not node and time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None)
        elif collapsing and expansionFactor and node and not time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_void_p(node)))
        elif collapsing and expansionFactor and node and time:
            return c_lib.criticalOverdensityCollapsingMassL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_void_p(node)))
    

class powerSpectrumWindowFunction:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.powerSpectrumWindowFunctionDestructorL(self._glcObj,self._classID)

    def wavenumberMaximum(self,smoothingMass):
        return c_lib.powerSpectrumWindowFunctionWavenumberMaximumL(self._glcObj,self._classID,smoothingMass)

    def value(self,wavenumber,smoothingMass):
        return c_lib.powerSpectrumWindowFunctionValueL(self._glcObj,self._classID,wavenumber,smoothingMass)

    def amplitudeIsMassIndependent(self):
        return c_lib.powerSpectrumWindowFunctionAmplitudeIsMassIndependentL(self._glcObj,self._classID)

class nbodyHaloMassError:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.nbodyHaloMassErrorDestructorL(self._glcObj,self._classID)

    def errorZeroAlways(self):
        return c_lib.nbodyHaloMassErrorErrorZeroAlwaysL(self._glcObj,self._classID)

    def errorFractional(self,node):
        return c_lib.nbodyHaloMassErrorErrorFractionalL(self._glcObj,self._classID,node)

    def correlation(self,node1,node2):
        return c_lib.nbodyHaloMassErrorCorrelationL(self._glcObj,self._classID,node1,node2)

class excursionSetBarrier:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.excursionSetBarrierDestructorL(self._glcObj,self._classID)

    def barrierGradient(self,variance,time,node,rateCompute):
        return c_lib.excursionSetBarrierBarrierGradientL(self._glcObj,self._classID,variance,time,node,rateCompute)

    def barrier(self,variance,time,node,rateCompute):
        return c_lib.excursionSetBarrierBarrierL(self._glcObj,self._classID,variance,time,node,rateCompute)

class haloMassFunction:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.haloMassFunctionDestructorL(self._glcObj,self._classID)

    def massFraction(self,time,massLow,massHigh,node=None):
        if not node:
            return c_lib.haloMassFunctionMassFractionL(self._glcObj,self._classID,time,massLow,massHigh,None)
        elif node:
            return c_lib.haloMassFunctionMassFractionL(self._glcObj,self._classID,time,massLow,massHigh,byref(c_void_p(node)))
    

    def integrated(self,time,massLow,massHigh,node=None):
        if not node:
            return c_lib.haloMassFunctionIntegratedL(self._glcObj,self._classID,time,massLow,massHigh,None)
        elif node:
            return c_lib.haloMassFunctionIntegratedL(self._glcObj,self._classID,time,massLow,massHigh,byref(c_void_p(node)))
    

    def differential(self,time,mass,node=None):
        if not node:
            return c_lib.haloMassFunctionDifferentialL(self._glcObj,self._classID,time,mass,None)
        elif node:
            return c_lib.haloMassFunctionDifferentialL(self._glcObj,self._classID,time,mass,byref(c_void_p(node)))
    

class intergalacticMediumFilteringMass:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.intergalacticMediumFilteringMassDestructorL(self._glcObj,self._classID)

    def massFilteringRateOfChange(self,time):
        return c_lib.intergalacticMediumFilteringMassMassFilteringRateOfChangeL(self._glcObj,self._classID,time)

    def massFiltering(self,time):
        return c_lib.intergalacticMediumFilteringMassMassFilteringL(self._glcObj,self._classID,time)

    def fractionBaryonsRateOfChange(self,mass,time):
        return c_lib.intergalacticMediumFilteringMassFractionBaryonsRateOfChangeL(self._glcObj,self._classID,mass,time)

    def fractionBaryonsGradientMass(self,mass,time):
        return c_lib.intergalacticMediumFilteringMassFractionBaryonsGradientMassL(self._glcObj,self._classID,mass,time)

    def fractionBaryons(self,mass,time):
        return c_lib.intergalacticMediumFilteringMassFractionBaryonsL(self._glcObj,self._classID,mass,time)

class darkMatterParticle:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.darkMatterParticleDestructorL(self._glcObj,self._classID)

    def mass(self):
        return c_lib.darkMatterParticleMassL(self._glcObj,self._classID)

class virialDensityContrast:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.virialDensityContrastDestructorL(self._glcObj,self._classID)

    def turnAroundOverVirialRadii(self,mass,time=None,expansionFactor=None,collapsing=None):
        if not collapsing and not expansionFactor and not time:
            return c_lib.virialDensityContrastTurnAroundOverVirialRadiiL(self._glcObj,self._classID,mass,None,None,None)
        elif not collapsing and not expansionFactor and time:
            return c_lib.virialDensityContrastTurnAroundOverVirialRadiiL(self._glcObj,self._classID,mass,byref(c_double(time)),None,None)
        elif not collapsing and expansionFactor and not time:
            return c_lib.virialDensityContrastTurnAroundOverVirialRadiiL(self._glcObj,self._classID,mass,None,byref(c_double(expansionFactor)),None)
        elif not collapsing and expansionFactor and time:
            return c_lib.virialDensityContrastTurnAroundOverVirialRadiiL(self._glcObj,self._classID,mass,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsing and not expansionFactor and not time:
            return c_lib.virialDensityContrastTurnAroundOverVirialRadiiL(self._glcObj,self._classID,mass,None,None,byref(c_bool(collapsing)))
        elif collapsing and not expansionFactor and time:
            return c_lib.virialDensityContrastTurnAroundOverVirialRadiiL(self._glcObj,self._classID,mass,byref(c_double(time)),None,byref(c_bool(collapsing)))
        elif collapsing and expansionFactor and not time:
            return c_lib.virialDensityContrastTurnAroundOverVirialRadiiL(self._glcObj,self._classID,mass,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)))
        elif collapsing and expansionFactor and time:
            return c_lib.virialDensityContrastTurnAroundOverVirialRadiiL(self._glcObj,self._classID,mass,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)))
    

    def isMassDependent(self):
        return c_lib.virialDensityContrastIsMassDependentL(self._glcObj,self._classID)

    def densityContrastRateOfChange(self,mass,time=None,expansionFactor=None,collapsing=None):
        if not collapsing and not expansionFactor and not time:
            return c_lib.virialDensityContrastDensityContrastRateOfChangeL(self._glcObj,self._classID,mass,None,None,None)
        elif not collapsing and not expansionFactor and time:
            return c_lib.virialDensityContrastDensityContrastRateOfChangeL(self._glcObj,self._classID,mass,byref(c_double(time)),None,None)
        elif not collapsing and expansionFactor and not time:
            return c_lib.virialDensityContrastDensityContrastRateOfChangeL(self._glcObj,self._classID,mass,None,byref(c_double(expansionFactor)),None)
        elif not collapsing and expansionFactor and time:
            return c_lib.virialDensityContrastDensityContrastRateOfChangeL(self._glcObj,self._classID,mass,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsing and not expansionFactor and not time:
            return c_lib.virialDensityContrastDensityContrastRateOfChangeL(self._glcObj,self._classID,mass,None,None,byref(c_bool(collapsing)))
        elif collapsing and not expansionFactor and time:
            return c_lib.virialDensityContrastDensityContrastRateOfChangeL(self._glcObj,self._classID,mass,byref(c_double(time)),None,byref(c_bool(collapsing)))
        elif collapsing and expansionFactor and not time:
            return c_lib.virialDensityContrastDensityContrastRateOfChangeL(self._glcObj,self._classID,mass,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)))
        elif collapsing and expansionFactor and time:
            return c_lib.virialDensityContrastDensityContrastRateOfChangeL(self._glcObj,self._classID,mass,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)))
    

    def densityContrast(self,mass,time=None,expansionFactor=None,collapsing=None):
        if not collapsing and not expansionFactor and not time:
            return c_lib.virialDensityContrastDensityContrastL(self._glcObj,self._classID,mass,None,None,None)
        elif not collapsing and not expansionFactor and time:
            return c_lib.virialDensityContrastDensityContrastL(self._glcObj,self._classID,mass,byref(c_double(time)),None,None)
        elif not collapsing and expansionFactor and not time:
            return c_lib.virialDensityContrastDensityContrastL(self._glcObj,self._classID,mass,None,byref(c_double(expansionFactor)),None)
        elif not collapsing and expansionFactor and time:
            return c_lib.virialDensityContrastDensityContrastL(self._glcObj,self._classID,mass,byref(c_double(time)),byref(c_double(expansionFactor)),None)
        elif collapsing and not expansionFactor and not time:
            return c_lib.virialDensityContrastDensityContrastL(self._glcObj,self._classID,mass,None,None,byref(c_bool(collapsing)))
        elif collapsing and not expansionFactor and time:
            return c_lib.virialDensityContrastDensityContrastL(self._glcObj,self._classID,mass,byref(c_double(time)),None,byref(c_bool(collapsing)))
        elif collapsing and expansionFactor and not time:
            return c_lib.virialDensityContrastDensityContrastL(self._glcObj,self._classID,mass,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)))
        elif collapsing and expansionFactor and time:
            return c_lib.virialDensityContrastDensityContrastL(self._glcObj,self._classID,mass,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)))
    

class haloEnvironment:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.haloEnvironmentDestructorL(self._glcObj,self._classID)

    def volumeFractionOccupied(self):
        return c_lib.haloEnvironmentVolumeFractionOccupiedL(self._glcObj,self._classID)

    def pdf(self,overdensity):
        return c_lib.haloEnvironmentPdfL(self._glcObj,self._classID,overdensity)

    def overdensityNonLinear(self,node):
        return c_lib.haloEnvironmentOverdensityNonLinearL(self._glcObj,self._classID,node)

    def overdensityLinearSet(self,node,overdensity):
        return c_lib.haloEnvironmentOverdensityLinearSetL(self._glcObj,self._classID,node,overdensity)

    def overdensityLinearMinimum(self):
        return c_lib.haloEnvironmentOverdensityLinearMinimumL(self._glcObj,self._classID)

    def overdensityLinearMaximum(self):
        return c_lib.haloEnvironmentOverdensityLinearMaximumL(self._glcObj,self._classID)

    def overdensityLinearGradientTime(self,node):
        return c_lib.haloEnvironmentOverdensityLinearGradientTimeL(self._glcObj,self._classID,node)

    def overdensityLinear(self,node,presentDay=None):
        if not presentDay:
            return c_lib.haloEnvironmentOverdensityLinearL(self._glcObj,self._classID,node,None)
        elif presentDay:
            return c_lib.haloEnvironmentOverdensityLinearL(self._glcObj,self._classID,node,byref(c_bool(presentDay)))
    

    def overdensityIsSettable(self):
        return c_lib.haloEnvironmentOverdensityIsSettableL(self._glcObj,self._classID)

    def environmentRadius(self):
        return c_lib.haloEnvironmentEnvironmentRadiusL(self._glcObj,self._classID)

    def environmentMass(self):
        return c_lib.haloEnvironmentEnvironmentMassL(self._glcObj,self._classID)

    def cdf(self,overdensity):
        return c_lib.haloEnvironmentCdfL(self._glcObj,self._classID,overdensity)

class linearGrowth:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.linearGrowthDestructorL(self._glcObj,self._classID)

    def value(self,time=None,expansionFactor=None,collapsing=None,normalize=None,component=None,wavenumber=None):
        if not collapsing and not component and not expansionFactor and not normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,None,None,None,None)
        elif not collapsing and not component and not expansionFactor and not normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,None,None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and not expansionFactor and not normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,None,None)
        elif not collapsing and not component and not expansionFactor and not normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and not expansionFactor and normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,None,byref(c_int(normalize)),None,None)
        elif not collapsing and not component and not expansionFactor and normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,None,byref(c_int(normalize)),None,byref(c_double(wavenumber)))
        elif not collapsing and not component and not expansionFactor and normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_int(normalize)),None,None)
        elif not collapsing and not component and not expansionFactor and normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_int(normalize)),None,byref(c_double(wavenumber)))
        elif not collapsing and not component and expansionFactor and not normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,None,None)
        elif not collapsing and not component and expansionFactor and not normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and expansionFactor and not normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,None,None)
        elif not collapsing and not component and expansionFactor and not normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and expansionFactor and normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_int(normalize)),None,None)
        elif not collapsing and not component and expansionFactor and normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_int(normalize)),None,byref(c_double(wavenumber)))
        elif not collapsing and not component and expansionFactor and normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_int(normalize)),None,None)
        elif not collapsing and not component and expansionFactor and normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_int(normalize)),None,byref(c_double(wavenumber)))
        elif not collapsing and component and not expansionFactor and not normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,None,None,byref(c_int(component)),None)
        elif not collapsing and component and not expansionFactor and not normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,None,None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and not expansionFactor and not normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,byref(c_int(component)),None)
        elif not collapsing and component and not expansionFactor and not normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and not expansionFactor and normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,None,byref(c_int(normalize)),byref(c_int(component)),None)
        elif not collapsing and component and not expansionFactor and normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,None,byref(c_int(normalize)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and not expansionFactor and normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_int(normalize)),byref(c_int(component)),None)
        elif not collapsing and component and not expansionFactor and normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_int(normalize)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and expansionFactor and not normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,byref(c_int(component)),None)
        elif not collapsing and component and expansionFactor and not normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and expansionFactor and not normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,byref(c_int(component)),None)
        elif not collapsing and component and expansionFactor and not normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and expansionFactor and normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_int(normalize)),byref(c_int(component)),None)
        elif not collapsing and component and expansionFactor and normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_int(normalize)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and expansionFactor and normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_int(normalize)),byref(c_int(component)),None)
        elif not collapsing and component and expansionFactor and normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_int(normalize)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and not component and not expansionFactor and not normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,None,None)
        elif collapsing and not component and not expansionFactor and not normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,None,byref(c_double(wavenumber)))
        elif collapsing and not component and not expansionFactor and not normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,None,None)
        elif collapsing and not component and not expansionFactor and not normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,None,byref(c_double(wavenumber)))
        elif collapsing and not component and not expansionFactor and normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_int(normalize)),None,None)
        elif collapsing and not component and not expansionFactor and normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_int(normalize)),None,byref(c_double(wavenumber)))
        elif collapsing and not component and not expansionFactor and normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_int(normalize)),None,None)
        elif collapsing and not component and not expansionFactor and normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_int(normalize)),None,byref(c_double(wavenumber)))
        elif collapsing and not component and expansionFactor and not normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None,None)
        elif collapsing and not component and expansionFactor and not normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None,byref(c_double(wavenumber)))
        elif collapsing and not component and expansionFactor and not normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None,None)
        elif collapsing and not component and expansionFactor and not normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None,byref(c_double(wavenumber)))
        elif collapsing and not component and expansionFactor and normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(normalize)),None,None)
        elif collapsing and not component and expansionFactor and normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(normalize)),None,byref(c_double(wavenumber)))
        elif collapsing and not component and expansionFactor and normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(normalize)),None,None)
        elif collapsing and not component and expansionFactor and normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(normalize)),None,byref(c_double(wavenumber)))
        elif collapsing and component and not expansionFactor and not normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,byref(c_int(component)),None)
        elif collapsing and component and not expansionFactor and not normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and not expansionFactor and not normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,byref(c_int(component)),None)
        elif collapsing and component and not expansionFactor and not normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and not expansionFactor and normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_int(normalize)),byref(c_int(component)),None)
        elif collapsing and component and not expansionFactor and normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_int(normalize)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and not expansionFactor and normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_int(normalize)),byref(c_int(component)),None)
        elif collapsing and component and not expansionFactor and normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_int(normalize)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and expansionFactor and not normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_int(component)),None)
        elif collapsing and component and expansionFactor and not normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and expansionFactor and not normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_int(component)),None)
        elif collapsing and component and expansionFactor and not normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and expansionFactor and normalize and not time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(normalize)),byref(c_int(component)),None)
        elif collapsing and component and expansionFactor and normalize and not time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(normalize)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and expansionFactor and normalize and time and not wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(normalize)),byref(c_int(component)),None)
        elif collapsing and component and expansionFactor and normalize and time and wavenumber:
            return c_lib.linearGrowthValueL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(normalize)),byref(c_int(component)),byref(c_double(wavenumber)))
    

    def logarithmicDerivativeWavenumber(self,time=None,expansionFactor=None,collapsing=None,component=None,wavenumber=None):
        if not collapsing and not component and not expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,None,None,None,None)
        elif not collapsing and not component and not expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,None,None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and not expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,None)
        elif not collapsing and not component and not expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and not component and expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and not component and expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,byref(c_double(wavenumber)))
        elif not collapsing and component and not expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,None,None,byref(c_int(component)),None)
        elif not collapsing and component and not expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,None,None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and not expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_int(component)),None)
        elif not collapsing and component and not expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_int(component)),None)
        elif not collapsing and component and expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_int(component)),None)
        elif not collapsing and component and expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and not component and not expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not component and not expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,byref(c_double(wavenumber)))
        elif collapsing and not component and not expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not component and not expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,byref(c_double(wavenumber)))
        elif collapsing and not component and expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and not component and expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_double(wavenumber)))
        elif collapsing and not component and expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and not component and expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_double(wavenumber)))
        elif collapsing and component and not expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_int(component)),None)
        elif collapsing and component and not expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and not expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_int(component)),None)
        elif collapsing and component and not expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(component)),None)
        elif collapsing and component and expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(component)),None)
        elif collapsing and component and expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeWavenumberL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(component)),byref(c_double(wavenumber)))
    

    def logarithmicDerivativeExpansionFactor(self,time=None,expansionFactor=None,collapsing=None,component=None,wavenumber=None):
        if not collapsing and not component and not expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,None,None,None,None)
        elif not collapsing and not component and not expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,None,None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and not expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,None)
        elif not collapsing and not component and not expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),None,None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and not component and expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,None,byref(c_double(wavenumber)))
        elif not collapsing and not component and expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,None)
        elif not collapsing and not component and expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,None,byref(c_double(wavenumber)))
        elif not collapsing and component and not expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,None,None,byref(c_int(component)),None)
        elif not collapsing and component and not expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,None,None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and not expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_int(component)),None)
        elif not collapsing and component and not expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),None,None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_int(component)),None)
        elif not collapsing and component and expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif not collapsing and component and expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_int(component)),None)
        elif not collapsing and component and expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),None,byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and not component and not expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not component and not expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),None,byref(c_double(wavenumber)))
        elif collapsing and not component and not expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,None)
        elif collapsing and not component and not expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),None,byref(c_double(wavenumber)))
        elif collapsing and not component and expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and not component and expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_double(wavenumber)))
        elif collapsing and not component and expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,None)
        elif collapsing and not component and expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),None,byref(c_double(wavenumber)))
        elif collapsing and component and not expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_int(component)),None)
        elif collapsing and component and not expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,None,byref(c_bool(collapsing)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and not expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_int(component)),None)
        elif collapsing and component and not expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),None,byref(c_bool(collapsing)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and expansionFactor and not time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(component)),None)
        elif collapsing and component and expansionFactor and not time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,None,byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(component)),byref(c_double(wavenumber)))
        elif collapsing and component and expansionFactor and time and not wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(component)),None)
        elif collapsing and component and expansionFactor and time and wavenumber:
            return c_lib.linearGrowthLogarithmicDerivativeExpansionFactorL(self._glcObj,self._classID,byref(c_double(time)),byref(c_double(expansionFactor)),byref(c_bool(collapsing)),byref(c_int(component)),byref(c_double(wavenumber)))
    

    def isWavenumberDependent(self,component=None):
        if not component:
            return c_lib.linearGrowthIsWavenumberDependentL(self._glcObj,self._classID,None)
        elif component:
            return c_lib.linearGrowthIsWavenumberDependentL(self._glcObj,self._classID,byref(c_int(component)))
    

class darkMatterHaloScale:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.darkMatterHaloScaleDestructorL(self._glcObj,self._classID)

    def velocityVirialGrowthRate(self,node):
        return c_lib.darkMatterHaloScaleVelocityVirialGrowthRateL(self._glcObj,self._classID,node)

    def velocityVirial(self,node):
        return c_lib.darkMatterHaloScaleVelocityVirialL(self._glcObj,self._classID,node)

    def timescaleDynamical(self,node):
        return c_lib.darkMatterHaloScaleTimescaleDynamicalL(self._glcObj,self._classID,node)

    def temperatureVirial(self,node):
        return c_lib.darkMatterHaloScaleTemperatureVirialL(self._glcObj,self._classID,node)

    def radiusVirialGrowthRate(self,node):
        return c_lib.darkMatterHaloScaleRadiusVirialGrowthRateL(self._glcObj,self._classID,node)

    def radiusVirialGradientLogarithmicMass(self,node):
        return c_lib.darkMatterHaloScaleRadiusVirialGradientLogarithmicMassL(self._glcObj,self._classID,node)

    def radiusVirial(self,node):
        return c_lib.darkMatterHaloScaleRadiusVirialL(self._glcObj,self._classID,node)

    def densityMeanGrowthRate(self,node):
        return c_lib.darkMatterHaloScaleDensityMeanGrowthRateL(self._glcObj,self._classID,node)

    def densityMean(self,node):
        return c_lib.darkMatterHaloScaleDensityMeanL(self._glcObj,self._classID,node)

class cosmologicalMassVariance:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.cosmologicalMassVarianceDestructorL(self._glcObj,self._classID)

    def sigma8(self):
        return c_lib.cosmologicalMassVarianceSigma8L(self._glcObj,self._classID)

    def rootVarianceLogarithmicGradientTime(self,mass,time):
        return c_lib.cosmologicalMassVarianceRootVarianceLogarithmicGradientTimeL(self._glcObj,self._classID,mass,time)

    def rootVarianceLogarithmicGradient(self,mass,time):
        return c_lib.cosmologicalMassVarianceRootVarianceLogarithmicGradientL(self._glcObj,self._classID,mass,time)

    def rootVarianceAndLogarithmicGradient(self,mass,time,rootVariance,rootVarianceLogarithmicGradient):
        return c_lib.cosmologicalMassVarianceRootVarianceAndLogarithmicGradientL(self._glcObj,self._classID,mass,time,rootVariance,rootVarianceLogarithmicGradient)

    def rootVariance(self,mass,time):
        return c_lib.cosmologicalMassVarianceRootVarianceL(self._glcObj,self._classID,mass,time)

    def powerNormalization(self):
        return c_lib.cosmologicalMassVariancePowerNormalizationL(self._glcObj,self._classID)

    def mass(self,rootVariance,time):
        return c_lib.cosmologicalMassVarianceMassL(self._glcObj,self._classID,rootVariance,time)

    def growthIsMassDependent(self):
        return c_lib.cosmologicalMassVarianceGrowthIsMassDependentL(self._glcObj,self._classID)

class intergalacticMediumState:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.intergalacticMediumStateDestructorL(self._glcObj,self._classID)

    def temperature(self,time):
        return c_lib.intergalacticMediumStateTemperatureL(self._glcObj,self._classID,time)

    def singlyIonizedHydrogenFraction(self,time):
        return c_lib.intergalacticMediumStateSinglyIonizedHydrogenFractionL(self._glcObj,self._classID,time)

    def singlyIonizedHeliumFraction(self,time):
        return c_lib.intergalacticMediumStateSinglyIonizedHeliumFractionL(self._glcObj,self._classID,time)

    def neutralHydrogenFraction(self,time):
        return c_lib.intergalacticMediumStateNeutralHydrogenFractionL(self._glcObj,self._classID,time)

    def neutralHeliumFraction(self,time):
        return c_lib.intergalacticMediumStateNeutralHeliumFractionL(self._glcObj,self._classID,time)

    def massJeans(self,time):
        return c_lib.intergalacticMediumStateMassJeansL(self._glcObj,self._classID,time)

    def electronScatteringTime(self,opticalDepth,assumeFullyIonized=None):
        if not assumeFullyIonized:
            return c_lib.intergalacticMediumStateElectronScatteringTimeL(self._glcObj,self._classID,opticalDepth,None)
        elif assumeFullyIonized:
            return c_lib.intergalacticMediumStateElectronScatteringTimeL(self._glcObj,self._classID,opticalDepth,byref(c_bool(assumeFullyIonized)))
    

    def electronScatteringOpticalDepth(self,time,assumeFullyIonized=None):
        if not assumeFullyIonized:
            return c_lib.intergalacticMediumStateElectronScatteringOpticalDepthL(self._glcObj,self._classID,time,None)
        elif assumeFullyIonized:
            return c_lib.intergalacticMediumStateElectronScatteringOpticalDepthL(self._glcObj,self._classID,time,byref(c_bool(assumeFullyIonized)))
    

    def electronFraction(self,time):
        return c_lib.intergalacticMediumStateElectronFractionL(self._glcObj,self._classID,time)

    def doublyIonizedHeliumFraction(self,time):
        return c_lib.intergalacticMediumStateDoublyIonizedHeliumFractionL(self._glcObj,self._classID,time)

class powerSpectrumPrimordialTransferred:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.powerSpectrumPrimordialTransferredDestructorL(self._glcObj,self._classID)

    def power(self,wavenumber,time):
        return c_lib.powerSpectrumPrimordialTransferredPowerL(self._glcObj,self._classID,wavenumber,time)

    def logarithmicDerivative(self,wavenumber,time):
        return c_lib.powerSpectrumPrimordialTransferredLogarithmicDerivativeL(self._glcObj,self._classID,wavenumber,time)

    def growthIsWavenumberDependent(self):
        return c_lib.powerSpectrumPrimordialTransferredGrowthIsWavenumberDependentL(self._glcObj,self._classID)

class powerSpectrumPrimordial:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.powerSpectrumPrimordialDestructorL(self._glcObj,self._classID)

    def power(self,wavenumber):
        return c_lib.powerSpectrumPrimordialPowerL(self._glcObj,self._classID,wavenumber)

    def logarithmicDerivative(self,wavenumber):
        return c_lib.powerSpectrumPrimordialLogarithmicDerivativeL(self._glcObj,self._classID,wavenumber)

class cosmologyParameters:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1

    # Destructor
    def __del__(self):
        c_lib.cosmologyParametersDestructorL(self._glcObj,self._classID)

    def temperatureCMB(self):
        return c_lib.cosmologyParametersTemperatureCMBL(self._glcObj,self._classID)

    def densityCritical(self):
        return c_lib.cosmologyParametersDensityCriticalL(self._glcObj,self._classID)

    def OmegaRadiation(self):
        return c_lib.cosmologyParametersOmegaRadiationL(self._glcObj,self._classID)

    def OmegaMatter(self):
        return c_lib.cosmologyParametersOmegaMatterL(self._glcObj,self._classID)

    def OmegaDarkEnergy(self):
        return c_lib.cosmologyParametersOmegaDarkEnergyL(self._glcObj,self._classID)

    def OmegaCurvature(self):
        return c_lib.cosmologyParametersOmegaCurvatureL(self._glcObj,self._classID)

    def OmegaBaryon(self):
        return c_lib.cosmologyParametersOmegaBaryonL(self._glcObj,self._classID)

    def HubbleConstant(self,units=None):
        if not units:
            return c_lib.cosmologyParametersHubbleConstantL(self._glcObj,self._classID,None)
        elif units:
            return c_lib.cosmologyParametersHubbleConstantL(self._glcObj,self._classID,byref(c_int(units)))
    

class transferFunctionAxionCAMB(transferFunction):

    # Constructor
    def __init__(self,darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_,redshift,axionCambCountPerDecade):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.transferFunctionAxionCAMBL(darkMatterParticle_._glcObj,darkMatterParticle_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,redshift,axionCambCountPerDecade)

class transferFunctionAccelerator(transferFunction):

    # Constructor
    def __init__(self,transferFunction_,tablePointsPerDecade):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 11
    
        self._glcObj = c_lib.transferFunctionAcceleratorL(transferFunction_._glcObj,transferFunction_._classID,tablePointsPerDecade)

class virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy(virialDensityContrast):

    # Constructor
    def __init__(self,tableStore,energyFixedAt,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 8
    
        self._glcObj = c_lib.virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgyL(tableStore,energyFixedAt,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class powerSpectrumWindowFunctionTopHatSmoothed(powerSpectrumWindowFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,sigma):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 8
    
        self._glcObj = c_lib.powerSpectrumWindowFunctionTopHatSmoothedL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,sigma)

class powerSpectrumWindowFunctionTopHatSharpKHybrid(powerSpectrumWindowFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,normalization,radiiRatio):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 7
    
        self._glcObj = c_lib.powerSpectrumWindowFunctionTopHatSharpKHybridL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,normalization,radiiRatio)

class powerSpectrumWindowFunctionTopHatGeneralized(powerSpectrumWindowFunction):

    # Constructor
    def __init__(self,mu,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 6
    
        self._glcObj = c_lib.powerSpectrumWindowFunctionTopHatGeneralizedL(mu,cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class powerSpectrumWindowFunctionTopHat(powerSpectrumWindowFunction):

    # Constructor
    def __init__(self,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.powerSpectrumWindowFunctionTopHatL(cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class powerSpectrumWindowFunctionSmoothKSpace(powerSpectrumWindowFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,beta,normalization):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.powerSpectrumWindowFunctionSmoothKSpaceL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,beta,normalization)

class powerSpectrumWindowFunctionSharpKSpace(powerSpectrumWindowFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,normalization):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.powerSpectrumWindowFunctionSharpKSpaceL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,normalization)

class powerSpectrumWindowFunctionLagrangianChan2017(powerSpectrumWindowFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,f):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.powerSpectrumWindowFunctionLagrangianChan2017L(cosmologyParameters_._glcObj,cosmologyParameters_._classID,f)

class powerSpectrumWindowFunctionETHOS(powerSpectrumWindowFunction):

    # Constructor
    def __init__(self,cW,beta,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.powerSpectrumWindowFunctionETHOSL(cW,beta,cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class virialDensityContrastBryanNorman1998(virialDensityContrast):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.virialDensityContrastBryanNorman1998L(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class powerSpectrumPrimordialTransferredSimple(powerSpectrumPrimordialTransferred):

    # Constructor
    def __init__(self,powerSpectrumPrimordial_,transferFunction_,linearGrowth_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.powerSpectrumPrimordialTransferredSimpleL(powerSpectrumPrimordial_._glcObj,powerSpectrumPrimordial_._classID,transferFunction_._glcObj,transferFunction_._classID,linearGrowth_._glcObj,linearGrowth_._classID)

class powerSpectrumPrimordialTransferredFile(powerSpectrumPrimordialTransferred):

    # Constructor
    def __init__(self,fileName,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.powerSpectrumPrimordialTransferredFileL(fileName,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class transferFunctionBode2001(transferFunction):

    # Constructor
    def __init__(self,transferFunctionCDM,epsilon,eta,nu,time,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.transferFunctionBode2001L(transferFunctionCDM._glcObj,transferFunctionCDM._classID,epsilon,eta,nu,time,cosmologyParameters_._glcObj,cosmologyParameters_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class powerSpectrumPrimordialPowerLaw(powerSpectrumPrimordial):

    # Constructor
    def __init__(self,index_,running,runningRunning,wavenumberReference,runningSmallScalesOnly):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.powerSpectrumPrimordialPowerLawL(index_,running,runningRunning,wavenumberReference,runningSmallScalesOnly)

class powerSpectrumPrimordialPiecewisePowerLaw(powerSpectrumPrimordial):

    # Constructor
    def __init__(self,index_,running,runningRunning,wavenumberReference):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.powerSpectrumPrimordialPiecewisePowerLawL(index_,running,runningRunning,wavenumberReference)

class powerSpectrumPrimordialCosmologicalCube(powerSpectrumPrimordial):

    # Constructor
    def __init__(self,lengthCube,wavenumberMinimumFactor,powerSpectrumPrimordial_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.powerSpectrumPrimordialCosmologicalCubeL(lengthCube,wavenumberMinimumFactor,powerSpectrumPrimordial_._glcObj,powerSpectrumPrimordial_._classID)

class transferFunctionBBKSWDM(transferFunction):

    # Constructor
    def __init__(self,transferFunctionCDM,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.transferFunctionBBKSWDML(transferFunctionCDM._glcObj,transferFunctionCDM._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class nbodyHaloMassErrorTrenti2010(nbodyHaloMassError):

    # Constructor
    def __init__(self,massParticle,correlationNormalization=None,correlationMassExponent=None,correlationRedshiftExponent=None,cosmologyFunctions_=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
        if cosmologyFunctions_:
            cosmologyFunctions__glcObj =cosmologyFunctions_._glcObj
            cosmologyFunctions__classID=cosmologyFunctions_._classID
        else:
            cosmologyFunctions__glcObj =None
            cosmologyFunctions__classID=None
    
        if not correlationMassExponent and not correlationNormalization and not correlationRedshiftExponent and not cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,None,None,None,None,None)
        elif not correlationMassExponent and not correlationNormalization and not correlationRedshiftExponent and cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,None,None,None,byref(c_void_p(cosmologyFunctions__glcObj)),byref(c_int(cosmologyFunctions__classID)))
        elif not correlationMassExponent and not correlationNormalization and correlationRedshiftExponent and not cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,None,None,byref(c_double(correlationRedshiftExponent)),None,None)
        elif not correlationMassExponent and not correlationNormalization and correlationRedshiftExponent and cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,None,None,byref(c_double(correlationRedshiftExponent)),byref(c_void_p(cosmologyFunctions__glcObj)),byref(c_int(cosmologyFunctions__classID)))
        elif not correlationMassExponent and correlationNormalization and not correlationRedshiftExponent and not cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,byref(c_double(correlationNormalization)),None,None,None,None)
        elif not correlationMassExponent and correlationNormalization and not correlationRedshiftExponent and cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,byref(c_double(correlationNormalization)),None,None,byref(c_void_p(cosmologyFunctions__glcObj)),byref(c_int(cosmologyFunctions__classID)))
        elif not correlationMassExponent and correlationNormalization and correlationRedshiftExponent and not cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,byref(c_double(correlationNormalization)),None,byref(c_double(correlationRedshiftExponent)),None,None)
        elif not correlationMassExponent and correlationNormalization and correlationRedshiftExponent and cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,byref(c_double(correlationNormalization)),None,byref(c_double(correlationRedshiftExponent)),byref(c_void_p(cosmologyFunctions__glcObj)),byref(c_int(cosmologyFunctions__classID)))
        elif correlationMassExponent and not correlationNormalization and not correlationRedshiftExponent and not cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,None,byref(c_double(correlationMassExponent)),None,None,None)
        elif correlationMassExponent and not correlationNormalization and not correlationRedshiftExponent and cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,None,byref(c_double(correlationMassExponent)),None,byref(c_void_p(cosmologyFunctions__glcObj)),byref(c_int(cosmologyFunctions__classID)))
        elif correlationMassExponent and not correlationNormalization and correlationRedshiftExponent and not cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,None,byref(c_double(correlationMassExponent)),byref(c_double(correlationRedshiftExponent)),None,None)
        elif correlationMassExponent and not correlationNormalization and correlationRedshiftExponent and cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,None,byref(c_double(correlationMassExponent)),byref(c_double(correlationRedshiftExponent)),byref(c_void_p(cosmologyFunctions__glcObj)),byref(c_int(cosmologyFunctions__classID)))
        elif correlationMassExponent and correlationNormalization and not correlationRedshiftExponent and not cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,byref(c_double(correlationNormalization)),byref(c_double(correlationMassExponent)),None,None,None)
        elif correlationMassExponent and correlationNormalization and not correlationRedshiftExponent and cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,byref(c_double(correlationNormalization)),byref(c_double(correlationMassExponent)),None,byref(c_void_p(cosmologyFunctions__glcObj)),byref(c_int(cosmologyFunctions__classID)))
        elif correlationMassExponent and correlationNormalization and correlationRedshiftExponent and not cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,byref(c_double(correlationNormalization)),byref(c_double(correlationMassExponent)),byref(c_double(correlationRedshiftExponent)),None,None)
        elif correlationMassExponent and correlationNormalization and correlationRedshiftExponent and cosmologyFunctions_:
            self._glcObj = c_lib.nbodyHaloMassErrorTrenti2010L(massParticle,byref(c_double(correlationNormalization)),byref(c_double(correlationMassExponent)),byref(c_double(correlationRedshiftExponent)),byref(c_void_p(cosmologyFunctions__glcObj)),byref(c_int(cosmologyFunctions__classID)))
    

class nbodyHaloMassErrorSOHaloFinder(nbodyHaloMassError):

    # Constructor
    def __init__(self,darkMatterHaloScale_,darkMatterProfileDMO_,massParticle):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.nbodyHaloMassErrorSOHaloFinderL(darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,massParticle)

class nbodyHaloMassErrorPowerLaw(nbodyHaloMassError):

    # Constructor
    def __init__(self,normalization,exponent,fractionalErrorHighMass):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.nbodyHaloMassErrorPowerLawL(normalization,exponent,fractionalErrorHighMass)

class nbodyHaloMassErrorNull(nbodyHaloMassError):

    # Constructor
    def __init__(self):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.nbodyHaloMassErrorNullL()

class nbodyHaloMassErrorFriendsOfFriends(nbodyHaloMassError):

    # Constructor
    def __init__(self,massParticle):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.nbodyHaloMassErrorFriendsOfFriendsL(massParticle)

class transferFunctionFileFuzzyDarkMatter(transferFunction):

    # Constructor
    def __init__(self,fileName,redshift,cosmologyParameters_,cosmologyFunctions_,darkMatterParticle_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 13
    
        self._glcObj = c_lib.transferFunctionFileFuzzyDarkMatterL(fileName,redshift,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID)

class linearGrowthNonClusteringBaryonsDarkMatter(linearGrowth):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.linearGrowthNonClusteringBaryonsDarkMatterL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class linearGrowthCollisionlessMatter(linearGrowth):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.linearGrowthCollisionlessMatterL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class linearGrowthBaryonsDarkMatter(linearGrowth):

    # Constructor
    def __init__(self,redshiftInitial,redshiftInitialDelta,cambCountPerDecade,cosmologyParameters_,cosmologyFunctions_,intergalacticMediumState_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.linearGrowthBaryonsDarkMatterL(redshiftInitial,redshiftInitialDelta,cambCountPerDecade,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,intergalacticMediumState_._glcObj,intergalacticMediumState_._classID)

class transferFunctionEisensteinHu1999(transferFunction):

    # Constructor
    def __init__(self,neutrinoNumberEffective,neutrinoMassSummed,darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 8
    
        self._glcObj = c_lib.transferFunctionEisensteinHu1999L(neutrinoNumberEffective,neutrinoMassSummed,darkMatterParticle_._glcObj,darkMatterParticle_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class intergalacticMediumStateSimple(intergalacticMediumState):

    # Constructor
    def __init__(self,reionizationRedshift,reionizationTemperature,preReionizationTemperature,cosmologyFunctions_,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.intergalacticMediumStateSimpleL(reionizationRedshift,reionizationTemperature,preReionizationTemperature,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class intergalacticMediumStateRecFast(intergalacticMediumState):

    # Constructor
    def __init__(self,cosmologyFunctions_,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.intergalacticMediumStateRecFastL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class intergalacticMediumStateInternal(intergalacticMediumState):

    # Constructor
    def __init__(self,cosmologyFunctions_,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.intergalacticMediumStateInternalL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class intergalacticMediumStateInstantReionization(intergalacticMediumState):

    # Constructor
    def __init__(self,cosmologyFunctions_,cosmologyParameters_,preReionizationState,reionizationTemperature,presentDayTemperature,reionizationRedshift=None,electronScatteringOpticalDepth=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        if not electronScatteringOpticalDepth and not reionizationRedshift:
            self._glcObj = c_lib.intergalacticMediumStateInstantReionizationL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,preReionizationState._glcObj,preReionizationState._classID,reionizationTemperature,presentDayTemperature,None,None)
        elif not electronScatteringOpticalDepth and reionizationRedshift:
            self._glcObj = c_lib.intergalacticMediumStateInstantReionizationL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,preReionizationState._glcObj,preReionizationState._classID,reionizationTemperature,presentDayTemperature,byref(c_double(reionizationRedshift)),None)
        elif electronScatteringOpticalDepth and not reionizationRedshift:
            self._glcObj = c_lib.intergalacticMediumStateInstantReionizationL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,preReionizationState._glcObj,preReionizationState._classID,reionizationTemperature,presentDayTemperature,None,byref(c_double(electronScatteringOpticalDepth)))
        elif electronScatteringOpticalDepth and reionizationRedshift:
            self._glcObj = c_lib.intergalacticMediumStateInstantReionizationL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,preReionizationState._glcObj,preReionizationState._classID,reionizationTemperature,presentDayTemperature,byref(c_double(reionizationRedshift)),byref(c_double(electronScatteringOpticalDepth)))
    

class intergalacticMediumStateFile(intergalacticMediumState):

    # Constructor
    def __init__(self,fileName,cosmologyFunctions_,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.intergalacticMediumStateFileL(fileName,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class transferFunctionCAMB(transferFunction):

    # Constructor
    def __init__(self,darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_,redshift,cambCountPerDecade):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.transferFunctionCAMBL(darkMatterParticle_._glcObj,darkMatterParticle_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,redshift,cambCountPerDecade)

class intergalacticMediumFilteringMassGnedin2000(intergalacticMediumFilteringMass):

    # Constructor
    def __init__(self,timeTooEarlyIsFatal,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,intergalacticMediumState_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.intergalacticMediumFilteringMassGnedin2000L(timeTooEarlyIsFatal,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,linearGrowth_._glcObj,linearGrowth_._classID,intergalacticMediumState_._glcObj,intergalacticMediumState_._classID)

class transferFunctionHu2000FDM(transferFunction):

    # Constructor
    def __init__(self,transferFunctionCDM,time,cosmologyParameters_,cosmologyFunctions_,darkMatterParticle_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 9
    
        self._glcObj = c_lib.transferFunctionHu2000FDML(transferFunctionCDM._glcObj,transferFunctionCDM._classID,time,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID)

class virialDensityContrastFriendsOfFriends(virialDensityContrast):

    # Constructor
    def __init__(self,linkingLength,densityRatio):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.virialDensityContrastFriendsOfFriendsL(linkingLength,densityRatio)

class haloMassFunctionTinker2008Generic(haloMassFunction):

    # Constructor
    def __init__(self,normalization,a,b,c,cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 10
    
        self._glcObj = c_lib.haloMassFunctionTinker2008GenericL(normalization,a,b,c,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class haloMassFunctionTinker2008(haloMassFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,cosmologyFunctions_,virialDensityContrast_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 8
    
        self._glcObj = c_lib.haloMassFunctionTinker2008L(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,virialDensityContrast_._glcObj,virialDensityContrast_._classID)

class haloMassFunctionSimpleSystematic(haloMassFunction):

    # Constructor
    def __init__(self,alpha,beta,cosmologyParameters_,referenceMassFunction):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 15
    
        self._glcObj = c_lib.haloMassFunctionSimpleSystematicL(alpha,beta,cosmologyParameters_._glcObj,cosmologyParameters_._classID,referenceMassFunction._glcObj,referenceMassFunction._classID)

class haloMassFunctionShethTormen(haloMassFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,a,p,normalization):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 7
    
        self._glcObj = c_lib.haloMassFunctionShethTormenL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID,a,p,normalization)

class haloMassFunctionRodriguezPuebla2016(haloMassFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 6
    
        self._glcObj = c_lib.haloMassFunctionRodriguezPuebla2016L(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class haloMassFunctionReed2007(haloMassFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.haloMassFunctionReed2007L(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID)

class haloMassFunctionPressSchechter(haloMassFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologicalMassVariance_,excursionSetFirstCrossing_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.haloMassFunctionPressSchechterL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,excursionSetFirstCrossing_._glcObj,excursionSetFirstCrossing_._classID)

class haloMassFunctionOndaroMallea2021(haloMassFunction):

    # Constructor
    def __init__(self,coefficientsN,coefficientsA,cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,haloMassFunction_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.haloMassFunctionOndaroMallea2021L(coefficientsN,coefficientsA,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID,haloMassFunction_._glcObj,haloMassFunction_._classID)

class haloMassFunctionFofBias(haloMassFunction):

    # Constructor
    def __init__(self,massFunctionIntrinsic,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,massParticle,linkingLength,linkingLengthIsComoving,massInfiniteToMassSharpEdge):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 14
    
        self._glcObj = c_lib.haloMassFunctionFofBiasL(massFunctionIntrinsic._glcObj,massFunctionIntrinsic._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,massParticle,linkingLength,linkingLengthIsComoving,massInfiniteToMassSharpEdge)

class haloMassFunctionErrorConvolved(haloMassFunction):

    # Constructor
    def __init__(self,massFunctionIntrinsic,cosmologyParameters_,nbodyHaloMassError_,errorFractionalMaximum):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 13
    
        self._glcObj = c_lib.haloMassFunctionErrorConvolvedL(massFunctionIntrinsic._glcObj,massFunctionIntrinsic._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,nbodyHaloMassError_._glcObj,nbodyHaloMassError_._classID,errorFractionalMaximum)

class haloMassFunctionEnvironmental(haloMassFunction):

    # Constructor
    def __init__(self,haloMassFunctionConditioned_,haloMassFunctionUnconditioned_,haloEnvironment_,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 12
    
        self._glcObj = c_lib.haloMassFunctionEnvironmentalL(haloMassFunctionConditioned_._glcObj,haloMassFunctionConditioned_._classID,haloMassFunctionUnconditioned_._glcObj,haloMassFunctionUnconditioned_._classID,haloEnvironment_._glcObj,haloEnvironment_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class haloMassFunctionEnvironmentAveraged(haloMassFunction):

    # Constructor
    def __init__(self,includeUnoccupiedVolume,haloMassFunctionConditioned_,haloMassFunctionUnconditioned_,haloEnvironment_,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 11
    
        self._glcObj = c_lib.haloMassFunctionEnvironmentAveragedL(includeUnoccupiedVolume,haloMassFunctionConditioned_._glcObj,haloMassFunctionConditioned_._classID,haloMassFunctionUnconditioned_._glcObj,haloMassFunctionUnconditioned_._classID,haloEnvironment_._glcObj,haloEnvironment_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class haloMassFunctionDespali2015(haloMassFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_,virialDensityContrast_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.haloMassFunctionDespali2015L(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID,virialDensityContrast_._glcObj,virialDensityContrast_._classID)

class haloMassFunctionBhattacharya2011(haloMassFunction):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,a,p,q,normalization):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.haloMassFunctionBhattacharya2011L(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID,a,p,q,normalization)

class transferFunctionIdentity(transferFunction):

    # Constructor
    def __init__(self,time):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 15
    
        self._glcObj = c_lib.transferFunctionIdentityL(time)

class haloEnvironmentUniform(haloEnvironment):

    # Constructor
    def __init__(self):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.haloEnvironmentUniformL()

class haloEnvironmentNormal(haloEnvironment):

    # Constructor
    def __init__(self,time,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,criticalOverdensity_,radiusEnvironment=None,massEnvironment=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        if not massEnvironment and not radiusEnvironment:
            self._glcObj = c_lib.haloEnvironmentNormalL(time,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID,None,None)
        elif not massEnvironment and radiusEnvironment:
            self._glcObj = c_lib.haloEnvironmentNormalL(time,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID,byref(c_double(radiusEnvironment)),None)
        elif massEnvironment and not radiusEnvironment:
            self._glcObj = c_lib.haloEnvironmentNormalL(time,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID,None,byref(c_double(massEnvironment)))
        elif massEnvironment and radiusEnvironment:
            self._glcObj = c_lib.haloEnvironmentNormalL(time,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID,byref(c_double(radiusEnvironment)),byref(c_double(massEnvironment)))
    

class haloEnvironmentLogNormal(haloEnvironment):

    # Constructor
    def __init__(self,radiusEnvironment,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,criticalOverdensity_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.haloEnvironmentLogNormalL(radiusEnvironment,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID)

class haloEnvironmentFixed(haloEnvironment):

    # Constructor
    def __init__(self,cosmologyFunctions_,linearGrowth_,overdensity,radiusEnvironment=None,massEnvironment=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        if not massEnvironment and not radiusEnvironment:
            self._glcObj = c_lib.haloEnvironmentFixedL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,linearGrowth_._glcObj,linearGrowth_._classID,overdensity,None,None)
        elif not massEnvironment and radiusEnvironment:
            self._glcObj = c_lib.haloEnvironmentFixedL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,linearGrowth_._glcObj,linearGrowth_._classID,overdensity,byref(c_double(radiusEnvironment)),None)
        elif massEnvironment and not radiusEnvironment:
            self._glcObj = c_lib.haloEnvironmentFixedL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,linearGrowth_._glcObj,linearGrowth_._classID,overdensity,None,byref(c_double(massEnvironment)))
        elif massEnvironment and radiusEnvironment:
            self._glcObj = c_lib.haloEnvironmentFixedL(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,linearGrowth_._glcObj,linearGrowth_._classID,overdensity,byref(c_double(radiusEnvironment)),byref(c_double(massEnvironment)))
    

class transferFunctionFile(transferFunction):

    # Constructor
    def __init__(self,fileName,redshift,cosmologyParameters_,cosmologyFunctions_,transferFunctionReference=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 12
        if transferFunctionReference:
            transferFunctionReference_glcObj =transferFunctionReference._glcObj
            transferFunctionReference_classID=transferFunctionReference._classID
        else:
            transferFunctionReference_glcObj =None
            transferFunctionReference_classID=None
    
        if not transferFunctionReference:
            self._glcObj = c_lib.transferFunctionFileL(fileName,redshift,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,None,None)
        elif transferFunctionReference:
            self._glcObj = c_lib.transferFunctionFileL(fileName,redshift,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,byref(c_void_p(transferFunctionReference_glcObj)),byref(c_int(transferFunctionReference_classID)))
    

class excursionSetFirstCrossingZhangHuiHighOrder(excursionSetFirstCrossing):

    # Constructor
    def __init__(self,excursionSetBarrier_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.excursionSetFirstCrossingZhangHuiHighOrderL(excursionSetBarrier_._glcObj,excursionSetBarrier_._classID)

class excursionSetFirstCrossingZhangHui(excursionSetFirstCrossing):

    # Constructor
    def __init__(self,excursionSetBarrier_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.excursionSetFirstCrossingZhangHuiL(excursionSetBarrier_._glcObj,excursionSetBarrier_._classID)

class excursionSetFirstCrossingLinearBarrier(excursionSetFirstCrossing):

    # Constructor
    def __init__(self,excursionSetBarrier_,cosmologicalMassVariance_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.excursionSetFirstCrossingLinearBarrierL(excursionSetBarrier_._glcObj,excursionSetBarrier_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID)

class excursionSetFirstCrossingFarahiMidpoint(excursionSetFirstCrossing):

    # Constructor
    def __init__(self,timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.excursionSetFirstCrossingFarahiMidpointL(timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,excursionSetBarrier_._glcObj,excursionSetBarrier_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID)

class excursionSetFirstCrossingFarahi(excursionSetFirstCrossing):

    # Constructor
    def __init__(self,timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.excursionSetFirstCrossingFarahiL(timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,excursionSetBarrier_._glcObj,excursionSetBarrier_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID)

class virialDensityContrastKitayamaSuto1996(virialDensityContrast):

    # Constructor
    def __init__(self,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.virialDensityContrastKitayamaSuto1996L(cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class excursionSetBarrierRemapShethMoTormen(excursionSetBarrier):

    # Constructor
    def __init__(self,a,b,c,applyTo,excursionSetBarrier_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.excursionSetBarrierRemapShethMoTormenL(a,b,c,applyTo,excursionSetBarrier_._glcObj,excursionSetBarrier_._classID)

class excursionSetBarrierRemapScale(excursionSetBarrier):

    # Constructor
    def __init__(self,factor,applyTo,excursionSetBarrier_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.excursionSetBarrierRemapScaleL(factor,applyTo,excursionSetBarrier_._glcObj,excursionSetBarrier_._classID)

class excursionSetBarrierQuadratic(excursionSetBarrier):

    # Constructor
    def __init__(self,coefficientConstant,coefficientLinear,coefficientQuadratic):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.excursionSetBarrierQuadraticL(coefficientConstant,coefficientLinear,coefficientQuadratic)

class excursionSetBarrierLinear(excursionSetBarrier):

    # Constructor
    def __init__(self,coefficientConstant,coefficientLinear):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.excursionSetBarrierLinearL(coefficientConstant,coefficientLinear)

class excursionSetBarrierCriticalOverdensity(excursionSetBarrier):

    # Constructor
    def __init__(self,criticalOverdensity_,cosmologicalMassVariance_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.excursionSetBarrierCriticalOverdensityL(criticalOverdensity_._glcObj,criticalOverdensity_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID)

class transferFunctionMurgia2017(transferFunction):

    # Constructor
    def __init__(self,transferFunctionCDM,alpha,beta,gamma,time,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 10
    
        self._glcObj = c_lib.transferFunctionMurgia2017L(transferFunctionCDM._glcObj,transferFunctionCDM._classID,alpha,beta,gamma,time,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class darkMatterProfileHeatingTwoBodyRelaxation(darkMatterProfileHeating):

    # Constructor
    def __init__(self,massParticle,lengthSoftening,timeStart,efficiency):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 7
    
        self._glcObj = c_lib.darkMatterProfileHeatingTwoBodyRelaxationL(massParticle,lengthSoftening,timeStart,efficiency)

class darkMatterProfileHeatingTidal(darkMatterProfileHeating):

    # Constructor
    def __init__(self,coefficientSecondOrder0,coefficientSecondOrder1,coefficientSecondOrder2,correlationVelocityRadius):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 6
    
        self._glcObj = c_lib.darkMatterProfileHeatingTidalL(coefficientSecondOrder0,coefficientSecondOrder1,coefficientSecondOrder2,correlationVelocityRadius)

class darkMatterProfileHeatingSummation(darkMatterProfileHeating):

    # Constructor
    def __init__(self,heatSources):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.darkMatterProfileHeatingSummationL(heatSources)

class darkMatterProfileHeatingNull(darkMatterProfileHeating):

    # Constructor
    def __init__(self):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.darkMatterProfileHeatingNullL()

class darkMatterProfileHeatingImpulsiveOutflow(darkMatterProfileHeating):

    # Constructor
    def __init__(self,impulsiveEnergyFactor,galacticStructure_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.darkMatterProfileHeatingImpulsiveOutflowL(impulsiveEnergyFactor,galacticStructure_)

class darkMatterProfileHeatingDDMv2(darkMatterProfileHeating):

    # Constructor
    def __init__(self,darkMatterParticle_,heating,massLoss,gamma):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.darkMatterProfileHeatingDDMv2L(darkMatterParticle_._glcObj,darkMatterParticle_._classID,heating,massLoss,gamma)

class darkMatterProfileHeatingDDM(darkMatterProfileHeating):

    # Constructor
    def __init__(self,lifetime,massSplitting):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.darkMatterProfileHeatingDDML(lifetime,massSplitting)

class virialDensityContrastPercolation(virialDensityContrast):

    # Constructor
    def __init__(self,linkingLength,cosmologyFunctions_,percolationObjects_,recursiveSelf=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
        if recursiveSelf:
            recursiveSelf_glcObj =recursiveSelf._glcObj
            recursiveSelf_classID=recursiveSelf._classID
        else:
            recursiveSelf_glcObj =None
            recursiveSelf_classID=None
    
        if not recursiveSelf:
            self._glcObj = c_lib.virialDensityContrastPercolationL(linkingLength,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,percolationObjects_,None,None)
        elif recursiveSelf:
            self._glcObj = c_lib.virialDensityContrastPercolationL(linkingLength,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,percolationObjects_,byref(c_void_p(recursiveSelf_glcObj)),byref(c_int(recursiveSelf_classID)))
    

class darkMatterProfileDMOZhao1996(darkMatterProfileDMO):

    # Constructor
    def __init__(self,alpha,beta,gamma,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 8
    
        self._glcObj = c_lib.darkMatterProfileDMOZhao1996L(alpha,beta,gamma,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMOTruncatedExponential(darkMatterProfileDMO):

    # Constructor
    def __init__(self,radiusFractionalDecay,alpha,beta,gamma,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 19
    
        self._glcObj = c_lib.darkMatterProfileDMOTruncatedExponentialL(radiusFractionalDecay,alpha,beta,gamma,nonAnalyticSolver,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMOTruncated(darkMatterProfileDMO):

    # Constructor
    def __init__(self,radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 18
    
        self._glcObj = c_lib.darkMatterProfileDMOTruncatedL(radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,nonAnalyticSolver,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMOSIDMIsothermal(darkMatterProfileDMO):

    # Constructor
    def __init__(self,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterParticle_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 7
    
        self._glcObj = c_lib.darkMatterProfileDMOSIDMIsothermalL(darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID)

class darkMatterProfileDMOSIDMCoreNFW(darkMatterProfileDMO):

    # Constructor
    def __init__(self,factorRadiusCore,darkMatterHaloScale_,darkMatterParticle_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 6
    
        self._glcObj = c_lib.darkMatterProfileDMOSIDMCoreNFWL(factorRadiusCore,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID)

class darkMatterProfileDMOPenarrubia2010(darkMatterProfileDMO):

    # Constructor
    def __init__(self,alpha,beta,gamma,betaStripped,muRadius,etaRadius,muVelocity,etaVelocity,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.darkMatterProfileDMOPenarrubia2010L(alpha,beta,gamma,betaStripped,muRadius,etaRadius,muVelocity,etaVelocity,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMONFW(darkMatterProfileDMO):

    # Constructor
    def __init__(self,velocityDispersionUseSeriesExpansion,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.darkMatterProfileDMONFWL(velocityDispersionUseSeriesExpansion,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMOMultiple(darkMatterProfileDMO):

    # Constructor
    def __init__(self,darkMatterProfileDMOHost_,darkMatterProfileDMOSatellite_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 17
    
        self._glcObj = c_lib.darkMatterProfileDMOMultipleL(darkMatterProfileDMOHost_._glcObj,darkMatterProfileDMOHost_._classID,darkMatterProfileDMOSatellite_._glcObj,darkMatterProfileDMOSatellite_._classID)

class darkMatterProfileDMOIsothermal(darkMatterProfileDMO):

    # Constructor
    def __init__(self,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 16
    
        self._glcObj = c_lib.darkMatterProfileDMOIsothermalL(darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMOHeatedMonotonic(darkMatterProfileDMO):

    # Constructor
    def __init__(self,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_,radiusFractionMinimum,radiusFractionMaximum,countPerDecadeRadius):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 15
    
        self._glcObj = c_lib.darkMatterProfileDMOHeatedMonotonicL(nonAnalyticSolver,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID,darkMatterProfileHeating_._glcObj,darkMatterProfileHeating_._classID,radiusFractionMinimum,radiusFractionMaximum,countPerDecadeRadius)

class darkMatterProfileDMOHeated(darkMatterProfileDMO):

    # Constructor
    def __init__(self,nonAnalyticSolver,velocityDispersionApproximate,toleranceRelativeVelocityDispersion,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 14
    
        self._glcObj = c_lib.darkMatterProfileDMOHeatedL(nonAnalyticSolver,velocityDispersionApproximate,toleranceRelativeVelocityDispersion,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID,darkMatterProfileHeating_._glcObj,darkMatterProfileHeating_._classID)

class darkMatterProfileDMOFiniteResolutionNFW(darkMatterProfileDMO):

    # Constructor
    def __init__(self,lengthResolution,massResolution,resolutionIsComoving,nonAnalyticSolver,darkMatterHaloScale_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 13
    
        self._glcObj = c_lib.darkMatterProfileDMOFiniteResolutionNFWL(lengthResolution,massResolution,resolutionIsComoving,nonAnalyticSolver,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class darkMatterProfileDMOFiniteResolution(darkMatterProfileDMO):

    # Constructor
    def __init__(self,lengthResolution,massResolution,resolutionIsComoving,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 12
    
        self._glcObj = c_lib.darkMatterProfileDMOFiniteResolutionL(lengthResolution,massResolution,resolutionIsComoving,nonAnalyticSolver,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class darkMatterProfileDMOEinasto(darkMatterProfileDMO):

    # Constructor
    def __init__(self,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.darkMatterProfileDMOEinastoL(darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMODecaying(darkMatterProfileDMO):

    # Constructor
    def __init__(self,darkMatterParticle_,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 11
    
        self._glcObj = c_lib.darkMatterProfileDMODecayingL(darkMatterParticle_._glcObj,darkMatterParticle_._classID,nonAnalyticSolver,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMOBurkert(darkMatterProfileDMO):

    # Constructor
    def __init__(self,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.darkMatterProfileDMOBurkertL(darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMOAccretionFlow(darkMatterProfileDMO):

    # Constructor
    def __init__(self,toleranceRelativePotential,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,accretionFlows_,darkMatterProfileDMO_,darkMatterHaloScale_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 10
    
        self._glcObj = c_lib.darkMatterProfileDMOAccretionFlowL(toleranceRelativePotential,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,criticalOverdensity_._glcObj,criticalOverdensity_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,accretionFlows_,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID,darkMatterHaloScale_._glcObj,darkMatterHaloScale_._classID)

class darkMatterProfileDMOAccelerator(darkMatterProfileDMO):

    # Constructor
    def __init__(self,toleranceRelative,factorRadiusMaximum,darkMatterProfileDMO_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 9
    
        self._glcObj = c_lib.darkMatterProfileDMOAcceleratorL(toleranceRelative,factorRadiusMaximum,darkMatterProfileDMO_._glcObj,darkMatterProfileDMO_._classID)

class virialDensityContrastSphericalCollapseBrynsDrkMttrDrkEnrgy(virialDensityContrast):

    # Constructor
    def __init__(self,tableStore,energyFixedAt,cosmologyParameters_,cosmologyFunctions_,intergalacticMediumFilteringMass_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 6
    
        self._glcObj = c_lib.virialDensityContrastSphericalCollapseBrynsDrkMttrDrkEnrgyL(tableStore,energyFixedAt,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,intergalacticMediumFilteringMass_._glcObj,intergalacticMediumFilteringMass_._classID)

class darkMatterParticleWDMThermal(darkMatterParticle):

    # Constructor
    def __init__(self,mass,degreesOfFreedomEffective,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.darkMatterParticleWDMThermalL(mass,degreesOfFreedomEffective,cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class darkMatterParticleSelfInteractingDarkMatter(darkMatterParticle):

    # Constructor
    def __init__(self,crossSectionSelfInteraction,darkMatterParticle_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.darkMatterParticleSelfInteractingDarkMatterL(crossSectionSelfInteraction,darkMatterParticle_._glcObj,darkMatterParticle_._classID)

class darkMatterParticleFuzzyDarkMatter(darkMatterParticle):

    # Constructor
    def __init__(self,mass,densityFraction):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.darkMatterParticleFuzzyDarkMatterL(mass,densityFraction)

class darkMatterParticleDecayingDarkMatter(darkMatterParticle):

    # Constructor
    def __init__(self,lifetime,massSplitting,darkMatterParticle_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.darkMatterParticleDecayingDarkMatterL(lifetime,massSplitting,darkMatterParticle_._glcObj,darkMatterParticle_._classID)

class darkMatterParticleCDM(darkMatterParticle):

    # Constructor
    def __init__(self):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.darkMatterParticleCDML()

class transferFunctionFuzzyDM(transferFunction):

    # Constructor
    def __init__(self,transferFunctionCDM,beta,gamma,time,cosmologyParameters_,cosmologyFunctions_,darkMatterParticle_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 14
    
        self._glcObj = c_lib.transferFunctionFuzzyDML(transferFunctionCDM._glcObj,transferFunctionCDM._classID,beta,gamma,time,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID)

class darkMatterHaloScaleVirialDensityContrastDefinition(darkMatterHaloScale):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,recursiveConstruct=None,recursiveSelf=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
        if recursiveSelf:
            recursiveSelf_glcObj =recursiveSelf._glcObj
            recursiveSelf_classID=recursiveSelf._classID
        else:
            recursiveSelf_glcObj =None
            recursiveSelf_classID=None
    
        if not recursiveConstruct and not recursiveSelf:
            self._glcObj = c_lib.darkMatterHaloScaleVirialDensityContrastDefinitionL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,virialDensityContrast_._glcObj,virialDensityContrast_._classID,None,None,None)
        elif not recursiveConstruct and recursiveSelf:
            self._glcObj = c_lib.darkMatterHaloScaleVirialDensityContrastDefinitionL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,virialDensityContrast_._glcObj,virialDensityContrast_._classID,None,byref(c_void_p(recursiveSelf_glcObj)),byref(c_int(recursiveSelf_classID)))
        elif recursiveConstruct and not recursiveSelf:
            self._glcObj = c_lib.darkMatterHaloScaleVirialDensityContrastDefinitionL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,virialDensityContrast_._glcObj,virialDensityContrast_._classID,byref(c_bool(recursiveConstruct)),None,None)
        elif recursiveConstruct and recursiveSelf:
            self._glcObj = c_lib.darkMatterHaloScaleVirialDensityContrastDefinitionL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,virialDensityContrast_._glcObj,virialDensityContrast_._classID,byref(c_bool(recursiveConstruct)),byref(c_void_p(recursiveSelf_glcObj)),byref(c_int(recursiveSelf_classID)))
    

class transferFunctionETHOSDM(transferFunction):

    # Constructor
    def __init__(self,transferFunctionCDM,alpha,beta,gamma,sigma,tau,kPeak,hPeak,h2,time,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 7
    
        self._glcObj = c_lib.transferFunctionETHOSDML(transferFunctionCDM._glcObj,transferFunctionCDM._classID,alpha,beta,gamma,sigma,tau,kPeak,hPeak,h2,time,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy(criticalOverdensity):

    # Constructor
    def __init__(self,linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,tableStore,normalization=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 9
    
        if not normalization:
            self._glcObj = c_lib.criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgyL(linearGrowth_._glcObj,linearGrowth_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,tableStore,None)
        elif normalization:
            self._glcObj = c_lib.criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgyL(linearGrowth_._glcObj,linearGrowth_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,tableStore,byref(c_double(normalization)))
    

class criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(criticalOverdensity):

    # Constructor
    def __init__(self,linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,tableStore,normalization=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 8
    
        if not normalization:
            self._glcObj = c_lib.criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstntL(linearGrowth_._glcObj,linearGrowth_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,tableStore,None)
        elif normalization:
            self._glcObj = c_lib.criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstntL(linearGrowth_._glcObj,linearGrowth_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,tableStore,byref(c_double(normalization)))
    

class criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy(criticalOverdensity):

    # Constructor
    def __init__(self,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,intergalacticMediumFilteringMass_,tableStore,normalization=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 7
    
        if not normalization:
            self._glcObj = c_lib.criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgyL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,intergalacticMediumFilteringMass_._glcObj,intergalacticMediumFilteringMass_._classID,tableStore,None)
        elif normalization:
            self._glcObj = c_lib.criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgyL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,intergalacticMediumFilteringMass_._glcObj,intergalacticMediumFilteringMass_._classID,tableStore,byref(c_double(normalization)))
    

class criticalOverdensityRenormalize(criticalOverdensity):

    # Constructor
    def __init__(self,criticalOverdensity_,cosmologyFunctions_,cosmologicalMassVariance_,cosmologicalMassVarianceReference_,linearGrowth_,massMatch=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 6
    
        if not massMatch:
            self._glcObj = c_lib.criticalOverdensityRenormalizeL(criticalOverdensity_._glcObj,criticalOverdensity_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,cosmologicalMassVarianceReference_._glcObj,cosmologicalMassVarianceReference_._classID,linearGrowth_._glcObj,linearGrowth_._classID,None)
        elif massMatch:
            self._glcObj = c_lib.criticalOverdensityRenormalizeL(criticalOverdensity_._glcObj,criticalOverdensity_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,cosmologicalMassVarianceReference_._glcObj,cosmologicalMassVarianceReference_._classID,linearGrowth_._glcObj,linearGrowth_._classID,byref(c_double(massMatch)))
    

class criticalOverdensityPeakBackgroundSplit(criticalOverdensity):

    # Constructor
    def __init__(self,criticalOverdensity_,haloEnvironment_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 5
    
        self._glcObj = c_lib.criticalOverdensityPeakBackgroundSplitL(criticalOverdensity_._glcObj,criticalOverdensity_._classID,haloEnvironment_._glcObj,haloEnvironment_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,linearGrowth_._glcObj,linearGrowth_._classID)

class criticalOverdensityMarsh2016FDM(criticalOverdensity):

    # Constructor
    def __init__(self,criticalOverdensityCDM,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,linearGrowth_,useFittingFunction):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.criticalOverdensityMarsh2016FDML(criticalOverdensityCDM._glcObj,criticalOverdensityCDM._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,linearGrowth_._glcObj,linearGrowth_._classID,useFittingFunction)

class criticalOverdensityKitayamaSuto1996(criticalOverdensity):

    # Constructor
    def __init__(self,linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.criticalOverdensityKitayamaSuto1996L(linearGrowth_._glcObj,linearGrowth_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID)

class criticalOverdensityFixed(criticalOverdensity):

    # Constructor
    def __init__(self,criticalOverdensity_,linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 4
    
        self._glcObj = c_lib.criticalOverdensityFixedL(criticalOverdensity_,linearGrowth_._glcObj,linearGrowth_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID)

class criticalOverdensityEnvironmental(criticalOverdensity):

    # Constructor
    def __init__(self,a,criticalOverdensity_,haloEnvironment_,cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.criticalOverdensityEnvironmentalL(a,criticalOverdensity_._glcObj,criticalOverdensity_._classID,haloEnvironment_._glcObj,haloEnvironment_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,linearGrowth_._glcObj,linearGrowth_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID)

class criticalOverdensityBarkana2001WDM(criticalOverdensity):

    # Constructor
    def __init__(self,criticalOverdensityCDM,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,linearGrowth_,useFittingFunction):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 10
    
        self._glcObj = c_lib.criticalOverdensityBarkana2001WDML(criticalOverdensityCDM._glcObj,criticalOverdensityCDM._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,darkMatterParticle_._glcObj,darkMatterParticle_._classID,linearGrowth_._glcObj,linearGrowth_._classID,useFittingFunction)

class virialDensityContrastFixed(virialDensityContrast):

    # Constructor
    def __init__(self,densityContrastValue,densityType,turnAroundOverVirialRadius,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.virialDensityContrastFixedL(densityContrastValue,densityType,turnAroundOverVirialRadius,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class cosmologyParametersSimple(cosmologyParameters):

    # Constructor
    def __init__(self,OmegaMatter,OmegaBaryon,OmegaDarkEnergy,temperatureCMB,HubbleConstant):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.cosmologyParametersSimpleL(OmegaMatter,OmegaBaryon,OmegaDarkEnergy,temperatureCMB,HubbleConstant)

class transferFunctionBBKS(transferFunction):

    # Constructor
    def __init__(self,darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.transferFunctionBBKSL(darkMatterParticle_._glcObj,darkMatterParticle_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class cosmologyFunctionsStaticUniverse(cosmologyFunctions):

    # Constructor
    def __init__(self,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 3
    
        self._glcObj = c_lib.cosmologyFunctionsStaticUniverseL(cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class cosmologyFunctionsMatterLambda(cosmologyFunctions):

    # Constructor
    def __init__(self,cosmologyParameters_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.cosmologyFunctionsMatterLambdaL(cosmologyParameters_._glcObj,cosmologyParameters_._classID)

class cosmologyFunctionsMatterDarkEnergy(cosmologyFunctions):

    # Constructor
    def __init__(self,cosmologyParameters_,darkEnergyEquationOfStateW0,darkEnergyEquationOfStateW1):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
    
        self._glcObj = c_lib.cosmologyFunctionsMatterDarkEnergyL(cosmologyParameters_._glcObj,cosmologyParameters_._classID,darkEnergyEquationOfStateW0,darkEnergyEquationOfStateW1)

class virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(virialDensityContrast):

    # Constructor
    def __init__(self,tableStore,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 7
    
        self._glcObj = c_lib.virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstntL(tableStore,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class cosmologicalMassVariancePeakBackgroundSplit(cosmologicalMassVariance):

    # Constructor
    def __init__(self,haloEnvironment_,cosmologicalMassVariance_,cosmologyParameters_,cosmologyFunctions_):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 2
    
        self._glcObj = c_lib.cosmologicalMassVariancePeakBackgroundSplitL(haloEnvironment_._glcObj,haloEnvironment_._classID,cosmologicalMassVariance_._glcObj,cosmologicalMassVariance_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID)

class cosmologicalMassVarianceFilteredPower(cosmologicalMassVariance):

    # Constructor
    def __init__(self,sigma8=None,cosmologicalMassVarianceReference=None,powerSpectrumPrimordialTransferredReference=None,wavenumberReference=None,tolerance=None,toleranceTopHat=None,nonMonotonicIsFatal=None,monotonicInterpolation=None,truncateAtParticleHorizon=None,cosmologyParameters_=None,cosmologyFunctions_=None,linearGrowth_=None,transferFunction_=None,powerSpectrumPrimordialTransferred_=None,powerSpectrumWindowFunction_=None,powerSpectrumWindowFunctionTopHat_=None):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 1
        if cosmologicalMassVarianceReference:
            cosmologicalMassVarianceReference_glcObj =cosmologicalMassVarianceReference._glcObj
            cosmologicalMassVarianceReference_classID=cosmologicalMassVarianceReference._classID
        else:
            cosmologicalMassVarianceReference_glcObj =None
            cosmologicalMassVarianceReference_classID=None
    
        if powerSpectrumPrimordialTransferredReference:
            powerSpectrumPrimordialTransferredReference_glcObj =powerSpectrumPrimordialTransferredReference._glcObj
            powerSpectrumPrimordialTransferredReference_classID=powerSpectrumPrimordialTransferredReference._classID
        else:
            powerSpectrumPrimordialTransferredReference_glcObj =None
            powerSpectrumPrimordialTransferredReference_classID=None
    
        if transferFunction_:
            transferFunction__glcObj =transferFunction_._glcObj
            transferFunction__classID=transferFunction_._classID
        else:
            transferFunction__glcObj =None
            transferFunction__classID=None
    
        if powerSpectrumWindowFunctionTopHat_:
            powerSpectrumWindowFunctionTopHat__glcObj =powerSpectrumWindowFunctionTopHat_._glcObj
            powerSpectrumWindowFunctionTopHat__classID=powerSpectrumWindowFunctionTopHat_._classID
        else:
            powerSpectrumWindowFunctionTopHat__glcObj =None
            powerSpectrumWindowFunctionTopHat__classID=None
    
        if not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif not cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),None,None,byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and not powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),None,None,byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and not powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),None,None)
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and not sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(None,byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and not transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),None,None,c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and not wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),None,c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
        elif cosmologicalMassVarianceReference and powerSpectrumPrimordialTransferredReference and powerSpectrumWindowFunctionTopHat_ and sigma8 and transferFunction_ and wavenumberReference:
            self._glcObj = c_lib.cosmologicalMassVarianceFilteredPowerL(byref(c_double(sigma8)),byref(c_void_p(cosmologicalMassVarianceReference_glcObj)),byref(c_int(cosmologicalMassVarianceReference_classID)),byref(c_void_p(powerSpectrumPrimordialTransferredReference_glcObj)),byref(c_int(powerSpectrumPrimordialTransferredReference_classID)),byref(c_double(wavenumberReference)),c_double(tolerance),c_double(toleranceTopHat),c_bool(nonMonotonicIsFatal),c_bool(monotonicInterpolation),c_bool(truncateAtParticleHorizon),c_void_p(cosmologyParameters_._glcObj),c_int(cosmologyParameters_._classID),c_void_p(cosmologyFunctions_._glcObj),c_int(cosmologyFunctions_._classID),c_void_p(linearGrowth_._glcObj),c_int(linearGrowth_._classID),byref(c_void_p(transferFunction__glcObj)),byref(c_int(transferFunction__classID)),c_void_p(powerSpectrumPrimordialTransferred_._glcObj),c_int(powerSpectrumPrimordialTransferred_._classID),c_void_p(powerSpectrumWindowFunction_._glcObj),c_int(powerSpectrumWindowFunction_._classID),byref(c_void_p(powerSpectrumWindowFunctionTopHat__glcObj)),byref(c_int(powerSpectrumWindowFunctionTopHat__classID)))
    

class transferFunctionCLASSCDM(transferFunction):

    # Constructor
    def __init__(self,darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_,redshift,countPerDecade):
        # Assign class ID so relevant pointers can be constructed on the Fortran side.
        self._classID = 6
    
        self._glcObj = c_lib.transferFunctionCLASSCDML(darkMatterParticle_._glcObj,darkMatterParticle_._classID,cosmologyParameters_._glcObj,cosmologyParameters_._classID,cosmologyFunctions_._glcObj,cosmologyFunctions_._classID,redshift,countPerDecade)

