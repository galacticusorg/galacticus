# Simple test of Python interface to libgalacticus
import galacticus
import ctypes
import numpy as np

# Simple tests of the Python API.
# Andrew Benson (23-April-2026)

# Constructs various objects and asserts that their methods return result that match expectations.
# Writes PASS/FAIL for each test.

# Cosmological parameters.
print("--- cosmologyParametersSimple ---")
cosmologyParameters = galacticus.cosmologyParametersSimple(0.3,0.045,0.7,2.78,70.0)
## Ω_m
OmegaMatter = cosmologyParameters.OmegaMatter()
status      = "PASS" if np.isclose(OmegaMatter, 0.3, rtol=1.0e-6) else "FAIL"
print(f'   {status}: Ω_m = {OmegaMatter:.2f}')
## H₀
HubbleConstant = cosmologyParameters.HubbleConstant()
status      = "PASS" if np.isclose(HubbleConstant, 70.0, rtol=1.0e-6) else "FAIL"
print(f'   {status}: H₀ = {HubbleConstant:.2f}')

# Cosmological functions.
print("--- cosmologyFunctionsMatterLambda ---")
cosmologyFunctions = galacticus.cosmologyFunctionsMatterLambda(cosmologyParameters)
time   = cosmologyFunctions.cosmicTime(1.0)
status = "PASS" if np.isclose(time, 13.466717044688263, rtol=1.0e-6) else "FAIL"
print(f'   {status}: t(a=1) = {time:.4f} Gyr')
expansionFactor = cosmologyFunctions.expansionFactor(time)
status = "PASS" if np.isclose(expansionFactor, 1.0, rtol=1.0e-6) else "FAIL"
print(f'   {status}: a(t={time:.2f}Gyr) = {expansionFactor:.2f}')
OmegaMatterEpochal = cosmologyFunctions.omegaMatterEpochal(time=6.0)
status = "PASS" if np.isclose(OmegaMatterEpochal, 0.7576520892903698, rtol=1.0e-6) else "FAIL"
print(f'   {status}: Ω_m(t=6.0Gyr) = {OmegaMatterEpochal:.2f}')
OmegaMatterEpochal = cosmologyFunctions.omegaMatterEpochal(expansionFactor=0.5)
status = "PASS" if np.isclose(OmegaMatterEpochal, 0.774193548387097, rtol=1.0e-6) else "FAIL"
print(f'   {status}: Ω_m(a=0.5) = {OmegaMatterEpochal:.2f}')

# Primordial power spectrum.
print("--- powerSpectrumPrimordialPowerLaw ---")
powerSpectrumPrimordial = galacticus.powerSpectrumPrimordialPowerLaw(0.965,0.0,0.0,1.0,False)
power = powerSpectrumPrimordial.power(2.0)
status = "PASS" if np.isclose(power, 1.9520635215524493, rtol=1.0e-6) else "FAIL"
print(f'   {status}: P₀(k=2 Mpc⁻¹) = {power:.2f} Mpc³')

# Window functions.
print("--- powerSpectrumWindowFunctionTopHat ---")
powerSpectrumWindowFunction = galacticus.powerSpectrumWindowFunctionTopHat(cosmologyParameters)
windowFunction = powerSpectrumWindowFunction.value(2.0,1.0e12,time)
status = "PASS" if np.isclose(windowFunction, 0.1780975268538394, rtol=1.0e-6) else "FAIL"
print(f'   {status}: W(k=2 Mpc⁻¹|M=10¹²M☉) = {windowFunction:.2f}')

# Dark matter particle.
darkMatterParticle = galacticus.darkMatterParticleCDM()

# Transfer functions.
print("--- transferFunctionCAMB ---")
transferFunction = galacticus.transferFunctionCAMB(darkMatterParticle,cosmologyParameters,cosmologyFunctions,0,0.0,0)
transferValue = transferFunction.value(2.0)
status = "PASS" if np.isclose(transferValue, 14637.776794245852, rtol=1.0e-6) else "FAIL"
print(f'   {status}: T(k=2 Mpc⁻¹) = {transferValue:.2f}')

# Linear growth.
print("--- linearGrowthCollisionlessMatter ---")
linearGrowth = galacticus.linearGrowthCollisionlessMatter(cosmologyParameters,cosmologyFunctions)
growthFunction = linearGrowth.value(time=6.0)
status = "PASS" if np.isclose(growthFunction, 0.6282911249247264, rtol=1.0e-6) else "FAIL"
print(f'   {status}: D(t=6 Gyr) = {growthFunction:.2f}')

# Power spectrum transferred
print("--- powerSpectrumPrimordialTransferredSimple ---")
powerSpectrumTransferred = galacticus.powerSpectrumPrimordialTransferredSimple(powerSpectrumPrimordial,transferFunction,linearGrowth)
powerTransferred = powerSpectrumTransferred.power(wavenumber=2.0,time=6.0)
status = "PASS" if np.isclose(powerTransferred, 165107209.2923229, rtol=1.0e-6) else "FAIL"
print(f'   {status}: P(k=2 Mpc⁻¹,t=6 Gyr) = {powerTransferred:.2f}')

# Cosmological mass variance.
print("--- cosmologicalMassVarianceFilteredPower ---")
cosmologicalMassVariance = galacticus.cosmologicalMassVarianceFilteredPower(sigma8=0.8,tolerance=1.0e-4,toleranceTopHat=1.0e-4,nonMonotonicIsFatal=True,integrationFailureIsFatal=True,monotonicInterpolation=False,rootVarianceLogarithmicGradientTolerance=1.0e-4,truncateAtParticleHorizon=False,storeTabulations=True,cosmologyParameters_=cosmologyParameters,cosmologyFunctions_=cosmologyFunctions,linearGrowth_=linearGrowth,powerSpectrumPrimordialTransferred_=powerSpectrumTransferred,powerSpectrumWindowFunction_=powerSpectrumWindowFunction)
rootVariance = cosmologicalMassVariance.rootVariance(mass=1.0e12,time=13.8)
status = "PASS" if np.isclose(rootVariance, 2.1968385715548044, rtol=1.0e-6) else "FAIL"
print(f'   {status}: σ(M=10¹²M☉) = {rootVariance:.2f}')

# Critical overdensity.
print("--- criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt ---")
criticalOverdensity = galacticus.criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(linearGrowth,cosmologyFunctions,cosmologicalMassVariance,darkMatterParticle,True)
deltaCrit = criticalOverdensity.value(time=13.8)
status = "PASS" if np.isclose(deltaCrit, 1.6750993420127407, rtol=1.0e-6) else "FAIL"
print(f'   {status}: δ_c(a=1) = {deltaCrit:.2f}')

# Virial density contrast.
print("--- virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt ---")
virialDensityContrast = galacticus.virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(True,cosmologyFunctions)
contrast = virialDensityContrast.densityContrast(mass=1.0e12,time=13.8)
status = "PASS" if np.isclose(contrast, 350.506868239381, rtol=1.0e-6) else "FAIL"
print(f'   {status}: Δ_vir(a=1) = {contrast:.2f}')

# Halo mass function.
print("--- haloMassFunctionShethTormen ---")
haloMassFunction = galacticus.haloMassFunctionShethTormen(cosmologyParameters,cosmologicalMassVariance,criticalOverdensity,0.707,0.3,0.322183)
haloMass = 1.0e10
hmf = haloMassFunction.differential(13.8,haloMass)*haloMass
status = "PASS" if np.isclose(hmf, 0.1040973175348119, rtol=1.0e-6) else "FAIL"
print(f'   {status}: dn/dlnM(M=10¹⁰M☉,a=1) = {hmf:.4f}')
