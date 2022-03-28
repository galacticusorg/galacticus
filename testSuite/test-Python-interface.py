# Simple test of Python interface to libgalacticus
import galacticus
import ctypes

cosmologyParameters = galacticus.cosmologyParametersSimple(0.3,0.045,0.7,2.78,70.0)
OmegaMatter = cosmologyParameters.OmegaMatter()
print "Omega_Matter = "+str(OmegaMatter)
HubbleConstant = cosmologyParameters.HubbleConstant()
print "Hubble_Constant = "+str(HubbleConstant)

cosmologyFunctions = galacticus.cosmologyFunctionsMatterLambda(cosmologyParameters)
time = cosmologyFunctions.cosmicTime(1.0)
print "time at a=1 is "+str(time)+" Gyr"
expansionFactor = cosmologyFunctions.expansionFactor(time)
print "expansion factor at time "+str(time)+" Gyr is "+str(expansionFactor)


OmegaMatterEpochal = cosmologyFunctions.omegaMatterEpochal(time=6.0)
print "Omega_Matter at t=6.0 is "+str(OmegaMatterEpochal)


OmegaMatterEpochal = cosmologyFunctions.omegaMatterEpochal(expansionFactor=0.5)
print "Omega_Matter at a=0.5 is "+str(OmegaMatterEpochal)

powerSpectrumPrimordial = galacticus.powerSpectrumPrimordialPowerLaw(0.965,0.0,0.0,1.0,False)
power = powerSpectrumPrimordial.power(2.0)
print "Primordial power spectrum at k=2/Mpc is "+str(power)

powerSpectrumWindowFunction = galacticus.powerSpectrumWindowFunctionTopHat(cosmologyParameters)
windowFunction = powerSpectrumWindowFunction.value(2.0,1.0e12)
print "Window function for 1.0e12 mass halo at k=2/Mpc is "+str(windowFunction)

darkMatterParticle = galacticus.darkMatterParticleCDM()

transferFunction = galacticus.transferFunctionCAMB(darkMatterParticle,cosmologyParameters,cosmologyFunctions,0.0,0)
transferValue = transferFunction.value(2.0)
print "CAMB transfer function at k=2/Mpc is "+str(transferValue)

transferFunctionFile = galacticus.transferFunctionFile(ctypes.c_char_p("/home/abenson/Galacticus/datasets/dynamic/largeScaleStructure/transfer_function_CAMB_7Z0TsZiN4r667I9DosWt8..hdf5"),0.0,cosmologyParameters,cosmologyFunctions,transferFunction)
transferValue = transferFunction.value(2.0)
print "File transfer function at k=2/Mpc is "+str(transferValue)

linearGrowth = galacticus.linearGrowthCollisionlessMatter(cosmologyParameters,cosmologyFunctions)
growthFunction = linearGrowth.value(time=6.0)
print "Growth function at t=6.0 is "+str(growthFunction)

powerSpectrumTransferred = galacticus.powerSpectrumPrimordialTransferredSimple(powerSpectrumPrimordial,transferFunction,linearGrowth)
powerTransferred = powerSpectrumTransferred.power(wavenumber=2.0,time=6.0)
print "Transferred power at k=2.0 and t=6.0 is "+str(powerTransferred)

cosmologicalMassVariance = galacticus.cosmologicalMassVarianceFilteredPower(sigma8=0.8,tolerance=1.0e-4,toleranceTopHat=1.0e-4,nonMonotonicIsFatal=True,monotonicInterpolation=False,truncateAtParticleHorizon=False,cosmologyParameters_=cosmologyParameters,cosmologyFunctions_=cosmologyFunctions,linearGrowth_=linearGrowth,powerSpectrumPrimordialTransferred_=powerSpectrumTransferred,powerSpectrumWindowFunction_=powerSpectrumWindowFunction)
rootVariance = cosmologicalMassVariance.rootVariance(mass=1.0e12,time=13.8)
print "Root variance at present day for M=1.0e12 is "+str(rootVariance)

criticalOverdensity = galacticus.criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(linearGrowth,cosmologyFunctions,cosmologicalMassVariance,darkMatterParticle,True)
deltaCrit = criticalOverdensity.value(time=13.8)
print "Critical overdensity at present day is "+str(deltaCrit)

virialDensityContrast = galacticus.virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(True,cosmologyFunctions)
contrast = virialDensityContrast.densityContrast(mass=1.0e12,time=13.8)
print "Virial density contrast now is "+str(contrast)

haloMassFunction = galacticus.haloMassFunctionShethTormen(cosmologyParameters,cosmologicalMassVariance,criticalOverdensity,0.707,0.3,0.322183)
haloMass = 1.0e10
hmf = haloMassFunction.differential(13.8,haloMass)*haloMass
print "Halo mass function is "+str(hmf)
