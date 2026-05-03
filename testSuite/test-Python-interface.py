# Simple test of Python interface to libgalacticus
import galacticus
import ctypes
import sys
import numpy as np
from contextlib import contextmanager

# Simple tests of the Python API.
# Andrew Benson (23-April-2026)

# Constructs various objects and asserts that their methods return results
# that match expectations.  Writes PASS/FAIL for each test, and exits with a
# non-zero status if any test failed (so CI catches regressions even if the
# script is interrupted before the failure summary is printed).

_failures = 0


def check(label, actual, expected, rtol=1.0e-6, fmt=".2f", unit=""):
    """Compare actual vs expected; print PASS/FAIL and tally failures."""
    global _failures
    suffix = f" {unit}" if unit else ""
    if np.isclose(actual, expected, rtol=rtol):
        print(f'   PASS: {label} = {actual:{fmt}}{suffix}')
    else:
        _failures += 1
        print(f'   FAIL: {label} = {actual:{fmt}}{suffix} '
              f'(expected {expected:{fmt}}, Δ={actual-expected:.3g})')


def check_eq(label, actual, expected):
    """Exact-equality variant for non-numeric values (strings, etc.)."""
    global _failures
    if actual == expected:
        print(f'   PASS: {label} = {actual!r}')
    else:
        _failures += 1
        print(f'   FAIL: {label} = {actual!r} (expected {expected!r})')


@contextmanager
def safe_section(name):
    """Wrap a test section so an unexpected exception is reported as FAIL
    rather than aborting the rest of the suite."""
    global _failures
    print(f"--- {name} ---")
    try:
        yield
    except Exception as exc:
        _failures += 1
        print(f'   FAIL: section raised {type(exc).__name__}: {exc}')


# Cosmological parameters.
with safe_section("cosmologyParametersSimple"):
    cosmologyParameters = galacticus.cosmologyParametersSimple(0.3,0.045,0.7,2.78,70.0)
    check("Ω_m", cosmologyParameters.OmegaMatter()   ,  0.3)
    check("H₀" , cosmologyParameters.HubbleConstant(), 70.0)

# Cosmological functions.
with safe_section("cosmologyFunctionsMatterLambda"):
    cosmologyFunctions = galacticus.cosmologyFunctionsMatterLambda(cosmologyParameters)
    time = cosmologyFunctions.cosmicTime(1.0)
    check("t(a=1)"             , time                                                       , 13.466717044688263, fmt=".4f", unit="Gyr")
    check(f"a(t={time:.2f}Gyr)", cosmologyFunctions.expansionFactor(time)                   ,  1.0)
    check("Ω_m(t=6.0Gyr)"      , cosmologyFunctions.omegaMatterEpochal(time=6.0)            ,  0.7576520892903698)
    check("Ω_m(a=0.5)"         , cosmologyFunctions.omegaMatterEpochal(expansionFactor=0.5) ,  0.774193548387097)

# Primordial power spectrum.
with safe_section("powerSpectrumPrimordialPowerLaw"):
    powerSpectrumPrimordial = galacticus.powerSpectrumPrimordialPowerLaw(0.965,0.0,0.0,1.0,False)
    check("P₀(k=2 Mpc⁻¹)", powerSpectrumPrimordial.power(2.0), 1.9520635215524493, unit="Mpc³")

# Window functions.
with safe_section("powerSpectrumWindowFunctionTopHat"):
    powerSpectrumWindowFunction = galacticus.powerSpectrumWindowFunctionTopHat(cosmologyParameters)
    check("W(k=2 Mpc⁻¹|M=10¹²M☉)", powerSpectrumWindowFunction.value(2.0,1.0e12,time), 0.1780975268538394)

# Dark matter particle.
with safe_section("darkMatterParticleCDM"):
    darkMatterParticle = galacticus.darkMatterParticleCDM()

# Transfer functions.
with safe_section("transferFunctionCAMB"):
    transferFunction = galacticus.transferFunctionCAMB(darkMatterParticle,cosmologyParameters,cosmologyFunctions,0,0.0,0)
    check("T(k=2 Mpc⁻¹)", transferFunction.value(2.0), 14637.776794245852)

# Linear growth.
with safe_section("linearGrowthCollisionlessMatter"):
    linearGrowth = galacticus.linearGrowthCollisionlessMatter(cosmologyParameters,cosmologyFunctions)
    check("D(t=6 Gyr)", linearGrowth.value(time=6.0), 0.6282911249247264)

# Power spectrum transferred.  Tolerance is loosened here because the result
# accumulates the integration tolerances of every upstream object.
with safe_section("powerSpectrumPrimordialTransferredSimple"):
    powerSpectrumTransferred = galacticus.powerSpectrumPrimordialTransferredSimple(powerSpectrumPrimordial,transferFunction,linearGrowth)
    check("P(k=2 Mpc⁻¹,t=6 Gyr)", powerSpectrumTransferred.power(wavenumber=2.0,time=6.0), 165107209.2923229, rtol=3.0e-6)

# Cosmological mass variance.
with safe_section("cosmologicalMassVarianceFilteredPower"):
    cosmologicalMassVariance = galacticus.cosmologicalMassVarianceFilteredPower(sigma8=0.8,tolerance=1.0e-4,toleranceTopHat=1.0e-4,nonMonotonicIsFatal=True,integrationFailureIsFatal=True,monotonicInterpolation=False,rootVarianceLogarithmicGradientTolerance=1.0e-4,truncateAtParticleHorizon=False,storeTabulations=True,cosmologyParameters_=cosmologyParameters,cosmologyFunctions_=cosmologyFunctions,linearGrowth_=linearGrowth,powerSpectrumPrimordialTransferred_=powerSpectrumTransferred,powerSpectrumWindowFunction_=powerSpectrumWindowFunction)
    check("σ(M=10¹²M☉)", cosmologicalMassVariance.rootVariance(mass=1.0e12,time=13.8), 2.1968385715548044)

# Critical overdensity.
with safe_section("criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt"):
    criticalOverdensity = galacticus.criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(linearGrowth,cosmologyFunctions,cosmologicalMassVariance,darkMatterParticle,True)
    check("δ_c(a=1)", criticalOverdensity.value(time=13.8), 1.6750993420127407)

# Virial density contrast.
with safe_section("virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt"):
    virialDensityContrast = galacticus.virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(True,cosmologyFunctions)
    check("Δ_vir(a=1)", virialDensityContrast.densityContrast(mass=1.0e12,time=13.8), 350.506868239381)

# Halo mass function.
with safe_section("haloMassFunctionShethTormen"):
    haloMassFunction = galacticus.haloMassFunctionShethTormen(cosmologyParameters,cosmologicalMassVariance,criticalOverdensity,0.707,0.3,0.322183)
    haloMass = 1.0e10
    check("dn/dlnM(M=10¹⁰M☉,a=1)", haloMassFunction.differential(13.8,haloMass)*haloMass, 0.1040973175348119, fmt=".4f")

# Random number generator.  Exercises method returns of `integer` (poisson),
# `integer(c_long)` (range/sample/seed), and the `integer(c_long)` argument
# path (sample(n=...)) plus the optional-arg branching for c_long.
with safe_section("randomNumberGeneratorGSL"):
    rng = galacticus.randomNumberGeneratorGSL(seed_=219,ompThreadOffset=False,mpiRankOffset=False)
    # TODO: replace dummy expectations below with golden values from a real run.
    check_eq("seed"        , rng.seed()                  ,         219)  # integer(c_long) return
    check_eq("rangeMinimum", rng.rangeMinimum()          ,           0)  # integer(c_long) return
    check_eq("rangeMaximum", rng.rangeMaximum()          ,  4294967295)  # integer(c_long) return
    check_eq("sample()"    , rng.sample()                ,           57953729)  # integer(c_long) return, no args
    check_eq("sample(n=10)", rng.sample(n=10)            ,           1)  # integer(c_long) optional arg
    check_eq("poisson(5.0)", rng.poissonSample(mean=5.0) ,           3)  # plain `integer` return
    check   ("uniform"     , rng.uniformSample()         ,         0.271657)  # double precision baseline
    
# Output times — exercises `integer(c_size_t)` return + arg.
with safe_section("outputTimesUniformSpacingInRedshift"):
    outputTimes = galacticus.outputTimesUniformSpacingInRedshift(0.0,5.0,10,cosmologyFunctions)
    # TODO: replace dummy expectations below with golden values from a real run.
    check_eq("count()"           , outputTimes.count()             ,  10)  # integer(c_size_t) return
    check   ("time(indexOutput=1)", outputTimes.time(indexOutput=1),  1.15473815)  # integer(c_size_t) arg, double return
    check   ("redshift(idx=10)"   , outputTimes.redshift(indexOutput=10),  0.0)

# Dark matter profile concentration — confirms partial-method exposure works
# (densityContrastDefinition/darkMatterProfileDMODefinition return class(...)
# and were auto-skipped, but the constructor and other methods should still
# be reachable).  Methods of this class take type(treeNode) which Python
# can't construct yet, so we just verify the constructor succeeds
darkMatterHaloScale = galacticus.darkMatterHaloScaleVirialDensityContrastDefinition(cosmologyParameters,cosmologyFunctions,virialDensityContrast)
darkMatterProfileDMO = galacticus.darkMatterProfileDMOIsothermal(darkMatterHaloScale)
with safe_section("darkMatterProfileConcentrationFixed"):
    concentration = galacticus.darkMatterProfileConcentrationFixed(10.0,virialDensityContrast,darkMatterProfileDMO)
    # No assertion — the constructor not raising and the subsequent destructor
    # running cleanly is the test.

# Initial mass function — exercises `type(varying_string)` return-type path
# (Fortran-side static c_char buffer + Python-side .decode("utf-8")).
with safe_section("initialMassFunctionSalpeter1955"):
    imf = galacticus.initialMassFunctionSalpeter1955()
    # TODO: replace dummy expectations below with golden values from a real run.
    check_eq("label"      , imf.label()                            ,    "Salpeter1955")  # varying_string return
    check   ("massMinimum", imf.massMinimum()                      ,    0.1)
    check   ("massMaximum", imf.massMaximum()                      ,  125.0)
    check   ("phi(M=1)"   , imf.phi(massInitial=1.0)               ,    0.170384)
    check   ("N(0.1..125)", imf.numberCumulative(massLower=0.1,massUpper=125.0), 2.82531)

# Final summary and exit code.
print(f"--- {_failures} failure(s) ---")
sys.exit(1 if _failures else 0)
