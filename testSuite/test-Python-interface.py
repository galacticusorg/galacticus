# Simple test of Python interface to libgalacticus
import galacticus
import ctypes
import inspect
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

# Output times list — exercises the 1D deferred-shape numeric array path
# (the constructor's `times` arg is `dimension(:)`).  Python passes a list
# / numpy array; the wrapper converts via np.ascontiguousarray and
# generates the (data_pointer, c_size_t(count)) pair the bind(c) function
# expects.  Without dimension support this whole class would have been
# auto-skipped at constructor-arg validation time.
with safe_section("outputTimesList"):
    outputTimesL = galacticus.outputTimesList([0.5, 1.0, 2.0, 5.0, 13.0], cosmologyFunctions)
    # TODO: replace dummy expectations below with golden values from a real run.
    check_eq("count() (5 entries)"     , outputTimesL.count()                   , 5)
    check   ("time(indexOutput=1)"     , outputTimesL.time(indexOutput=1)       , 0.5)
    check   ("time(indexOutput=5)"     , outputTimesL.time(indexOutput=5)       , 13.0)

# Dark matter profile concentration — exercises the class(FooClass) return
# path: densityContrastDefinition() returns class(virialDensityContrastClass),
# which the wrapper resolves into the matching Python subclass via
# virialDensityContrast._from_classID.  The returned object is owned by
# Galacticus (not by Python), so its destructor is a no-op (_owned=False)
# and we should still be able to call methods on it.
darkMatterHaloScale = galacticus.darkMatterHaloScaleVirialDensityContrastDefinition(cosmologyParameters,cosmologyFunctions,virialDensityContrast)
darkMatterProfileDMO = galacticus.darkMatterProfileDMOIsothermal(darkMatterHaloScale)
with safe_section("darkMatterProfileConcentrationFixed"):
    concentration = galacticus.darkMatterProfileConcentrationFixed(10.0,virialDensityContrast,darkMatterProfileDMO)
    vdcReturned   = concentration.densityContrastDefinition()
    # The returned wrapper should be the same concrete subclass we passed in.
    check_eq("class(...) return type"   , type(vdcReturned).__name__           , type(virialDensityContrast).__name__)
    check_eq("class(...) return _owned" , getattr(vdcReturned, '_owned', True) , False)
    # And the wrapped object should be usable for further method calls; the
    # value should match the original instance (same Galacticus object).
    check   ("Δ_vir(via returned)"      , vdcReturned.densityContrast(mass=1.0e12,time=13.8), 350.506868239381)

# Node property extractor — exercises `type(enumerationXxxType)` return path
# (the inner method gives a derived enum type, the wrapper lifts the %ID
# component out as c_int).
with safe_section("nodePropertyExtractorNodeMajorMergerTime"):
    npe = galacticus.nodePropertyExtractorNodeMajorMergerTime()
    # TODO: replace dummy expectations below with golden values from a real run.
    check_eq("type()"    , npe.type()    , 0)  # type(enumerationOutputAnalysisPropertyType    Type) → c_int
    check_eq("quantity()", npe.quantity(), 0)  # type(enumerationOutputAnalysisPropertyQuantityType) → c_int

# Abstract-intermediate class-argument path — `galacticFilterHighPass`'s
# constructor takes `class(nodePropertyExtractorScalar)`, which is an
# *abstract* derived type extending the registered functionClass base
# `nodePropertyExtractorClass`.  The wrapper now accepts these through
# the base's `nodePropertyExtractorGetPtr` and narrows the polymorphic
# pointer to the intermediate via a `select type` block at runtime.
# Without the abstract-intermediate path, the whole galacticFilter
# *Pass family (high / low / interval) would have been rejected at
# constructor-arg validation time.  We pass a concrete
# `nodePropertyExtractorNodeMajorMergerTime` (which extends
# `nodePropertyExtractorScalar`) — the select-type match must succeed
# at runtime for construction to complete without Error_Report.
with safe_section("galacticFilterHighPass (abstract-intermediate arg)"):
    gFilter = galacticus.galacticFilterHighPass(
        threshold              = 1.0e12,
        nodePropertyExtractor_ = npe,
    )
    check_eq("constructed type",
             type(gFilter).__name__, 'galacticFilterHighPass')
    # Confirm the sibling impls in the same family also expose
    # constructors that accept the intermediate — they share the
    # same pipeline path, so a single end-to-end construction is
    # enough; the others just need to be importable.
    check_eq("galacticFilterLowPass exposed",
             hasattr(galacticus, 'galacticFilterLowPass'), True)
    check_eq("galacticFilterIntervalPass exposed",
             hasattr(galacticus, 'galacticFilterIntervalPass'), True)
    # And the property-extractor-tree intermediates (DescendantNode,
    # HostNode) take `class(nodePropertyExtractorScalar)` too —
    # confirm they survived the constructor-arg validation.
    check_eq("nodePropertyExtractorDescendantNode exposed",
             hasattr(galacticus, 'nodePropertyExtractorDescendantNode'), True)
    check_eq("nodePropertyExtractorHostNode exposed",
             hasattr(galacticus, 'nodePropertyExtractorHostNode'), True)

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

# Radiative transfer photon packet — exercises both fixed-size dimensional
# shapes in one round-trip: positionSet/directionSet take a `dimension(3)`
# *input* and position()/direction() return a `dimension(3)` numpy array.
# The constructor takes only scalar doubles, so this is fully self-contained.
with safe_section("radiativeTransferPhotonPacketSimple"):
    photon = galacticus.radiativeTransferPhotonPacketSimple(
        wavelength=500.0, wavelengthMinimum=400.0, wavelengthMaximum=700.0,
        luminosity=1.0e30,
    )
    photon.positionSet ([1.0, 2.0, 3.0])
    photon.directionSet([0.0, 0.0, 1.0])
    pos = photon.position ()
    dir = photon.direction()
    # Plain round-trip — values stored go straight back out, so exact
    # equality is reliable here.
    check_eq("position[0]" , float(pos[0]), 1.0)
    check_eq("position[1]" , float(pos[1]), 2.0)
    check_eq("position[2]" , float(pos[2]), 3.0)
    check_eq("direction[2]", float(dir[2]), 1.0)
    # Returned object should be a numpy array of length 3.
    check_eq("position type" , type(pos).__name__, 'ndarray')
    check_eq("position size" , pos.size           , 3)
    # Wrong-size input should raise ValueError (size guard in the wrapper).
    try:
        photon.positionSet([1.0, 2.0])
    except ValueError as exc:
        check_eq("ValueError on size mismatch", "expects 3" in str(exc), True)
    else:
        check_eq("ValueError on size mismatch", "no exception raised", "ValueError")

# Non-central χ² (degree 3) distribution — exercises the Python-keyword
# escape on a constructor argument.  The Fortran constructor takes
# `lambda` as a parameter; without renaming the Python signature to
# `lambda_` the wrapper module fails to import (SyntaxError on
# `def __init__(self, lambda, ...)`).
with safe_section("distributionFunction1DNonCentralChiDegree3"):
    chi3 = galacticus.distributionFunction1DNonCentralChiDegree3(
        lambda_=2.5, randomNumberGenerator_=rng,
    )
    # density() at the mode-ish region should be finite and positive;
    # the value itself isn't load-bearing — surviving the call is.
    d = chi3.density(x=2.5)
    check_eq("density(x=2.5) is finite", np.isfinite(d) and d > 0.0, True)

# Supernovae Type Ia power-law DTD (differential) — exercises the
# Python-keyword escape on a method NAME.  The Fortran class declares
# `<method name="yield">`; `yield` is a reserved word in Python, so the
# wrapper renames it to `yield_` (PEP 8).  Just confirm the rename:
# calling it would need a stellar age / metallicity dataset we don't
# initialise here.
with safe_section("supernovaeTypeIaPowerLawDTDDifferential"):
    sn1a = galacticus.supernovaeTypeIaPowerLawDTDDifferential(
        timeMinimum=0.04, exponent=-1.0, normalization=2.0e-3,
    )
    check_eq("yield_ method exposed" , hasattr(sn1a, 'yield_'), True)
    check_eq("'yield' not exposed"   , hasattr(sn1a, 'yield' ), False)

# Spherical computational-domain volume integrator — exercises the
# procedure-pointer-arg skip.  The class's `integrate(integrand)`
# method takes `procedure(...)` which the pipeline can't translate, so
# the wrapper drops just that method while keeping the rest of the
# class.  `volume()` is a plain double-precision return so it survives.
with safe_section("computationalDomainVolumeIntegratorSpherical"):
    cdom = galacticus.computationalDomainVolumeIntegratorSpherical([1.0, 5.0])
    check   ("volume()"               , cdom.volume(), (4.0/3.0)*np.pi*(5.0**3 - 1.0**3))
    check_eq("integrate() not exposed", hasattr(cdom, 'integrate'), False)

# Empirical UniverseMachine node operator — exercises the
# continuation-character strip in the Internal-constructor arg-name
# capture: the source declares the constructor across many
# continuation lines (24 scalar args + 3 functionClass deps), and
# without stripping `&` from each captured token the args leaked into
# the emitted `<referenceConstruct>` directive and broke XML parsing.
# Constructing the object end-to-end confirms the multi-line opener
# was parsed cleanly.
with safe_section("nodeOperatorEmpiricalGalaxyUniverseMachine"):
    um = galacticus.nodeOperatorEmpiricalGalaxyUniverseMachine(
        massStellarFinal=-1.0, fractionMassSpheroid=0.0, fractionMassDisk=1.0,
        epsilon_0=-1.435, epsilon_a= 1.831, epsilon_lna= 1.368, epsilon_z=-0.217,
        M_0      =12.035, M_a      = 4.556, M_lna      = 4.417, M_z      =-0.731,
        alpha_0  = 1.963, alpha_a  =-2.316, alpha_lna  =-1.732, alpha_z  = 0.178,
        beta_0   = 0.482, beta_a   =-0.841, beta_z     =-0.471,
        gamma_0  =-1.034, gamma_a  =-3.100, gamma_z    =-1.055,
        delta_0  = 0.411,
        redshiftMaximum=15.0, massHaloMinimum=1.0e10,
        cosmologyParameters_  =cosmologyParameters,
        cosmologyFunctions_   =cosmologyFunctions,
        virialDensityContrast_=virialDensityContrast,
    )
    check_eq("constructed type",
             type(um).__name__, 'nodeOperatorEmpiricalGalaxyUniverseMachine')

# Position-interpolated node operator — exercises the kind_int8 (=
# selected_int_kind(18), 64-bit) integer array path.  Without the
# kind_int8 → c_long mapping, `nodeIndicesReport` would be emitted as
# `integer(c_int)` (32-bit) and mismatch the inner constructor.  We
# include a value > 2^31 to make the 64-bit-ness load-bearing.
with safe_section("nodeOperatorPositionInterpolated"):
    indices = np.array([1, 42, 1_000_000_000_000], dtype=np.int64)
    nopi = galacticus.nodeOperatorPositionInterpolated(
        lengthBox=100.0, wrapPeriodic=True,
        nodeIndicesReport=indices,
        cosmologyFunctions_=cosmologyFunctions,
    )
    check_eq("constructed type",
             type(nopi).__name__, 'nodeOperatorPositionInterpolated')

# `mergerTreeImporter` — exercises the kind-aliased integer method
# *return* types.  `treeCount` / `nodeCount` / `subhaloTraceCount`
# declare `integer(kind=c_size_t)` and `treeIndex` declares
# `integer(kind=kind_int8)`.  Before the alias normalisation in
# `_normalize_method_return_type`, the generator's return-type switch
# only had branches for the unprefixed `integer(c_size_t)` /
# `integer(c_long)` forms, so these methods fell through and were
# dropped with an "unsupported method return type" caution.
#
# Calling these end-to-end needs an HDF5 merger-tree file (they each
# call `galacticusForestIndicesRead` on first use, which opens the
# file), so the meaningful check is the wrapper-exposure level — the
# same pattern used for `radiativeTransferMatter` above.  We also pin
# the negative case (`subhaloTraceCount` stays skipped because its
# `class(nodeData)` arg is in the deferred non-fc-hierarchy bucket) so
# accidentally widening the predicate without the related work is
# caught.
with safe_section("mergerTreeImporter (kind-aliased integer returns)"):
    check_eq("mergerTreeImporterGalacticus exposed",
             hasattr(galacticus, 'mergerTreeImporterGalacticus'), True)
    for method in ('treeCount', 'treeIndex', 'nodeCount'):
        check_eq(f"mergerTreeImporterGalacticus.{method} exposed",
                 hasattr(galacticus.mergerTreeImporterGalacticus, method),
                 True)
    # `subhaloTraceCount` survived the return-type fix but still has the
    # `class(nodeData)` arg blocker (non-fc hierarchy — out of scope).
    check_eq("subhaloTraceCount stays skipped (class(nodeData) arg)",
             hasattr(galacticus.mergerTreeImporterGalacticus,
                     'subhaloTraceCount'),
             False)

# Fixed-length character-array constructor argument — exercises the
# `character(len=N), dimension(:)` pipeline path.  The
# `radiativeTransferMatterAtomic` constructor takes
# `character(len=2), dimension(:) :: elements`; without the new
# code-gen path it (and its parent `radiativeTransferMatter`, plus
# `computationalDomain` which transitively depends on it) would have
# been rejected at constructor-arg validation time.  The constructor
# itself needs nine atomic-physics functionClass dependencies that
# we don't build here — confirming the wrapper symbol exists is the
# meaningful end-to-end check that the emission succeeded.
#
# `computationalDomain` has three impls but two (Cartesian3D,
# cylindrical) take `dimension(3,2)` / `dimension(2,2)` boundaries
# arrays which the pipeline still doesn't support; only the
# spherical impl (with `dimension(2)`) registers, and that's the one
# we check for.
with safe_section("radiativeTransferMatter (character len=N array path)"):
    check_eq("radiativeTransferMatterAtomic exposed",
             hasattr(galacticus, 'radiativeTransferMatterAtomic'), True)
    check_eq("computationalDomainSpherical exposed",
             hasattr(galacticus, 'computationalDomainSpherical'),  True)

# Multivariate normal — exercises the 2D deferred-shape numeric array
# constructor argument path: the inner constructor takes
# `double precision, dimension(:,:) :: covariance`, which the wrapper
# now flattens at the bind(c) boundary, ships as
# (data_pointer, shape[0], shape[1]), and reshapes inside (Fortran
# column-major; the Python side hands ctypes an `np.asfortranarray`
# buffer so the layouts match without a transpose).  Without 2D
# support, this whole class would have failed constructor-arg
# validation.
with safe_section("distributionFunctionMultivariateNormal"):
    mvn = galacticus.distributionFunctionMultivariateNormal(
        mean=[0.0, 0.0],
        covariance=[[1.0, 0.0], [0.0, 1.0]],   # 2D identity; passed as nested list
        errorAbsolute     =1.0e-6,
        errorRelative     =1.0e-6,
        countTrialsMaximum=100,
        randomNumberGenerator_=rng,
    )
    check_eq("constructed type",
             type(mvn).__name__, 'distributionFunctionMultivariateNormal')

# `outputAnalysisTargetDataStandard` — a tiny functionClass that bundles the
# seven coupled-but-individually-optional fields the 1D function output-analyses
# (mean / scatter / volume) used to expose as separate constructor args.
# Bundling them collapses the wrapper-pipeline's optional-argument branching
# from 2^12+ to 2^6 on each affected outer constructor (see source/output.analyses.target_data.F90).
# This test exercises three things:
#   1. Round-trip from Python: construct an instance via the wrapper, call a
#      method on it, get a value back.  `hasTarget()` returns true iff both
#      target arrays are allocated, so it gives us a clean boolean witness.
#   2. The bundle is genuinely partial-fill safe: defaults work, label-only
#      construction works, full construction works.
#   3. The refactored outer-class constructors actually carry `targetData_`
#      in their Python signature (we don't try to instantiate one — that
#      would need a half-dozen functionClass dependencies — just confirm
#      the signature).
with safe_section("outputAnalysisTargetDataStandard"):
    # Default-construct: every field omitted.  hasTarget() is False because
    # neither the value nor the covariance array was allocated.
    td_empty = galacticus.outputAnalysisTargetDataStandard()
    check_eq("hasTarget() default-constructed", td_empty.hasTarget(), False)

    # Labels + log-scale flags only — still no target dataset, so
    # hasTarget() remains False.
    td_partial = galacticus.outputAnalysisTargetDataStandard(
        xAxisLabel = 'log10(M_halo / M_sun)',
        yAxisLabel = 'M_HI / M_sun',
        xAxisIsLog = True,
        yAxisIsLog = True,
    )
    check_eq("hasTarget() label-only construction", td_partial.hasTarget(), False)

    # Fully populated, including a 3-element value array and 3x3 covariance
    # (which exercises the 1D + 2D array boundary plumbing in the wrapper).
    td_full = galacticus.outputAnalysisTargetDataStandard(
        xAxisLabel       = 'log10(M_halo / M_sun)',
        yAxisLabel       = 'log10(M_HI / M_sun)',
        targetLabel      = 'Test target',
        xAxisIsLog       = True,
        yAxisIsLog       = True,
        valueTarget      = [9.0, 9.5, 10.0],
        covarianceTarget = [[0.1, 0.0, 0.0],
                            [0.0, 0.2, 0.0],
                            [0.0, 0.0, 0.3]],
    )
    check_eq("hasTarget() fully populated", td_full.hasTarget(), True)

    # And confirm the refactored outer constructor exposes the bundled object
    # under the expected keyword name.  This is what proves the
    # 2^N branching reduction is wired up: each refactored constructor now
    # has ONE optional argument for the bundle instead of seven.
    sig = inspect.signature(galacticus.outputAnalysisVolumeFunction1D.__init__)
    check_eq("outputAnalysisVolumeFunction1D exposes targetData_",
             'targetData_' in sig.parameters, True)
    sig = inspect.signature(galacticus.outputAnalysisMeanFunction1D.__init__)
    check_eq("outputAnalysisMeanFunction1D exposes targetData_",
             'targetData_' in sig.parameters, True)
    sig = inspect.signature(galacticus.outputAnalysisScatterFunction1D.__init__)
    check_eq("outputAnalysisScatterFunction1D exposes targetData_",
             'targetData_' in sig.parameters, True)

# Stellar population spectra postprocessor builder (lookup impl) —
# exercises BOTH new array-arg paths in a single constructor:
#   * `type(varying_string), dimension(:) :: names` — Python list of
#     str is encoded to ASCII at runtime, padded to the longest
#     element, and shipped as a contiguous S{N} byte buffer plus
#     count + per-element-length companions; the Fortran wrapper
#     repacks into `type(varying_string), dimension(:)`.
#   * `type(stellarPopulationSpectraPostprocessorList), dimension(:)
#     :: postprocessors` — Galacticus's idiom for an array of
#     polymorphic pointers; the Python wrapper builds parallel
#     `(c_void_p * n)` and `(c_int * n)` ctypes arrays from each
#     instance's `_glcObj` / `_classID`, and the Fortran wrapper
#     rebuilds the wrapper-list via `stellarPopulationSpectra`
#     `PostprocessorGetPtr` element-by-element.
# `build(descriptor=...)` then exercises the class(...) return path on
# top of all that, so the round-trip lands back in the right Python
# subclass via _from_classID.
with safe_section("stellarPopulationSpectraPostprocessorBuilderLookup"):
    ppIdentity = galacticus.stellarPopulationSpectraPostprocessorIdentity ()
    ppInoue    = galacticus.stellarPopulationSpectraPostprocessorInoue2014()
    builder    = galacticus.stellarPopulationSpectraPostprocessorBuilderLookup(
        names          = ['default', 'igm'         ],
        postprocessors = [ppIdentity, ppInoue      ],
    )
    ppDefault = builder.build(descriptor='default')
    ppIGM     = builder.build(descriptor='igm'    )
    check_eq("build('default') resolves to identity subclass",
             type(ppDefault).__name__,
             'stellarPopulationSpectraPostprocessorIdentity')
    check_eq("build('igm') resolves to Inoue2014 subclass",
             type(ppIGM    ).__name__,
             'stellarPopulationSpectraPostprocessorInoue2014')

# Abstract-intermediate class-argument path (massDistribution family) —
# `massDistributionSphericalScaler`'s constructor takes
# `class(massDistributionSpherical)`, the abstract intermediate that
# extends the registered `massDistributionClass`.  We construct a
# concrete spherical (`massDistributionZero` — dimensionless, no other
# deps) and pass it through; the wrapper must route via
# `massDistributionGetPtr` and successfully narrow to
# `massDistributionSpherical` for construction to complete.
#
# Also confirms the related impls in the same family that hit the
# intermediate-arg path are exposed.  `massDistributionCylindricalScaler`
# takes `class(massDistributionCylindrical)` (a separate intermediate
# under the same base) so it stresses the second branch of the
# hierarchy walk.
with safe_section("massDistributionSphericalScaler (abstract-intermediate arg)"):
    massDistributionZero    = galacticus.massDistributionZero(dimensionless=True)
    massDistributionScaler  = galacticus.massDistributionSphericalScaler(
        factorScalingLength = 2.0,
        factorScalingMass   = 3.0,
        massDistribution_   = massDistributionZero,
    )
    check_eq("constructed type",
             type(massDistributionScaler).__name__,
             'massDistributionSphericalScaler')
    # Other impls in the spherical / cylindrical intermediate families
    # — every one of these required the abstract-intermediate path to
    # survive constructor-arg validation.
    for impl in ('massDistributionCylindricalScaler',
                 'massDistributionSphericalDecaying',
                 'massDistributionSphericalHeated',
                 'massDistributionSphericalTruncated',
                 'massDistributionSphericalFiniteResolution'):
        check_eq(f"{impl} exposed", hasattr(galacticus, impl), True)

# Explicit-lower-bound fixed-size 1D array — exercises the
# `dimension(L:U)` shape (e.g. `dimension(0:2)`).  The wrapper accepts
# any contiguous sequence whose length matches U-L+1 and ships it to
# the inner constructor, whose dummy has the explicit lower bound;
# Fortran copies elementwise on entry, so the lower bound is irrelevant
# at the wire level — only the element count matters.  Without this
# support, `haloMassFunctionOndaroMallea2021` would have been rejected
# at constructor-arg validation time.
#
# coefficientsN is `dimension(0:2)` (size 3); coefficientsA is
# `dimension(0:1)` (size 2).  We pass the matching counts to confirm
# the size validator picks them up correctly — and check that passing
# a wrong-length array raises ValueError with the corrected size in
# the diagnostic.
with safe_section("haloMassFunctionOndaroMallea2021 (dimension(0:2))"):
    haloMassFunctionOM = galacticus.haloMassFunctionOndaroMallea2021(
        coefficientsN             = [-0.04, 0.03, -0.001],
        coefficientsA             = [ 1.0 , 0.2          ],
        cosmologyParameters_      = cosmologyParameters,
        cosmologicalMassVariance_ = cosmologicalMassVariance,
        linearGrowth_             = linearGrowth,
        haloMassFunction_         = haloMassFunction,
    )
    check_eq("constructed type",
             type(haloMassFunctionOM).__name__,
             'haloMassFunctionOndaroMallea2021')
    # Wrong-length coefficientsA should be rejected by the wrapper
    # with a message that includes the correct expected size (2).
    try:
        galacticus.haloMassFunctionOndaroMallea2021(
            coefficientsN             = [-0.04, 0.03, -0.001],
            coefficientsA             = [ 1.0 , 0.2 ,  0.05 ],   # wrong size
            cosmologyParameters_      = cosmologyParameters,
            cosmologicalMassVariance_ = cosmologicalMassVariance,
            linearGrowth_             = linearGrowth,
            haloMassFunction_         = haloMassFunction,
        )
    except ValueError as exc:
        check_eq("ValueError on wrong-length dimension(0:1)",
                 "expects 2" in str(exc), True)
    else:
        check_eq("ValueError on wrong-length dimension(0:1)",
                 "no exception raised", "ValueError")

# `logical, dimension(:)` constructor arg path — `outputMask` on the
# stellar-luminosity property extractor.  bind(c) ships the values as a
# `logical(c_bool), dimension(*)` buffer plus a c_size_t count; the
# wrapper allocates a default-kind `logical, dimension(:)` local and
# populates it via an elemental `logical()` cast before the inner call.
#
# Constructing one of these impls end-to-end needs filter / stellar
# state we don't initialise here, so the meaningful check is that the
# wrapper symbol exists and exposes `outputMask` in its signature
# (parallels the `radiativeTransferMatter` smoke test above).
with safe_section("nodePropertyExtractor* (logical(:) outputMask)"):
    for impl in ('nodePropertyExtractorLuminosityStellar',
                 'nodePropertyExtractorLmnstyStllrCF2000',
                 'nodePropertyExtractorLmnstyEmssnLineAGN',
                 'nodePropertyExtractorLmnstyEmssnLinePanuzzo2003'):
        check_eq(f"{impl} exposed", hasattr(galacticus, impl), True)
        sig = inspect.signature(getattr(galacticus, impl).__init__)
        check_eq(f"{impl}: outputMask in signature",
                 'outputMask' in sig.parameters, True)

# Multi-rank fixed-shape numeric array constructor arg —
# `double precision, dimension(N, M[, K…])`.  The wrapper converts the
# user's nested-list / numpy input to a Fortran-order (column-major)
# buffer, validates `.shape` matches the expected tuple, and passes
# the contiguous data pointer to the bind(c) function whose dummy is
# declared with the full fixed shape.  No count companions are
# inserted — every axis size is baked into the bind(c) declaration.
#
# `computationalDomainVolumeIntegratorCartesian3D` takes a single
# `boundaries: dimension(3,2)` arg encoding {x,y,z}{min,max}.  We
# construct with a unit-thick cuboid (1*2*3=6) and round-trip `volume()`
# to confirm both that the column-major byte layout matches Fortran's
# convention and that the shape validator allows the correct shape
# (a wrong shape should raise ValueError).
with safe_section("computationalDomainVolumeIntegratorCartesian3D (dimension(3,2))"):
    cdomInteg = galacticus.computationalDomainVolumeIntegratorCartesian3D(
        boundaries = [[0.0, 1.0],    # xMin, xMax
                      [0.0, 2.0],    # yMin, yMax
                      [0.0, 3.0]],   # zMin, zMax
    )
    check("volume()", cdomInteg.volume(), 1.0 * 2.0 * 3.0)
    # Wrong-shape input is rejected with a message naming the expected
    # tuple (proves the validator uses array_shape, not just .size).
    try:
        galacticus.computationalDomainVolumeIntegratorCartesian3D(
            boundaries = [[0.0, 1.0, 5.0],    # wrong axis count
                          [0.0, 2.0, 5.0]],
        )
    except ValueError as exc:
        check_eq("ValueError on wrong shape",
                 "expects shape (3, 2)" in str(exc), True)
    else:
        check_eq("ValueError on wrong shape", "no exception raised", "ValueError")

# Ambiguous-Internal-constructor disambiguation —
# `<constructor internal="..."/>` in libraryClasses.xml picks one
# constructor when the impl exposes more than one
# Internal-suffixed module procedure.  Without the hint, the generator
# can't choose and skips the impl entirely.
# `darkMatterProfileConcentrationDuttonMaccio2014` has two:
# `InternalType` (keyed by a fit-type enum) and `InternalDefined`
# (raw coefficients).  We expose the higher-level `InternalType` form.
with safe_section("darkMatterProfileConcentrationDuttonMaccio2014 (internal=…)"):
    check_eq("darkMatterProfileConcentrationDuttonMaccio2014 exposed",
             hasattr(galacticus, 'darkMatterProfileConcentrationDuttonMaccio2014'),
             True)
    sig = inspect.signature(
        galacticus.darkMatterProfileConcentrationDuttonMaccio2014.__init__)
    # The chosen InternalType form takes `fitType` (enum) +
    # cosmologyParameters_ / cosmologyFunctions_; the rejected
    # InternalDefined form would have surfaced `a1`/`a2`/...
    check_eq("fitType in signature (InternalType chosen)",
             'fitType' in sig.parameters, True)
    check_eq("a1 not in signature (InternalDefined rejected)",
             'a1' in sig.parameters, False)

# Null-filled constructor arg path — `<argument name="..." value="null"/>`
# overrides in libraryClasses.xml tell the wrapper to drop callback-
# injection escape-hatch args (a procedure-pointer plus paired class(*)
# state slots) that the parameter-driven path of the impl already passes
# as null.  Without the override, the procedure-pointer arg blocks the
# whole impl at constructor-arg validation time.
#
# Both impls share three null-filled args
# (initializationFunction / initializationSelf / initializationArgument);
# constructing one of them end-to-end needs a stack of other
# functionClass deps we don't initialise here, so the meaningful end-to-
# end check is that the wrapper symbols exist (parallels the
# `radiativeTransferMatter` smoke test above).
with safe_section("value='null' constructor-arg overrides"):
    for impl in ('massDistributionSphericalAdiabaticGnedin2004',
                 'massDistributionSphericalSIDMIsothermalBaryons'):
        check_eq(f"{impl} exposed", hasattr(galacticus, impl), True)
        # And the three null-filled args really are absent from the
        # Python signature — that's the contract the override promises.
        sig = inspect.signature(getattr(galacticus, impl).__init__)
        for arg_name in ('initializationFunction',
                         'initializationSelf',
                         'initializationArgument'):
            check_eq(f"{impl}: {arg_name} not in signature",
                     arg_name in sig.parameters, False)

# Final summary and exit code.
print(f"--- {_failures} failure(s) ---")
sys.exit(1 if _failures else 0)
