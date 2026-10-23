#!/usr/bin/env python3
"""Check broadband luminosities computed from SEDs for lightcone output.

The model outputs a lightcone with fixed-age star formation histories, for which a single SED template serves every
lightcone crossing. Observed-frame luminosities are computed from a rest-frame SED by integrating it through the filter
response blueshifted into the rest frame of the galaxy, with any redshift-dependent postprocessing (here, IGM attenuation)
folded into the filter rather than the SED. Checked here are that:

 * the rest- and observed-frame luminosities agree with an independent integration of the output rest-frame SED through
   the filter response curves;
 * the observed-frame luminosities with IGM attenuation agree with the same independent integration, with the attenuation of
   Inoue et al. (2014) evaluated independently and applied at the observed wavelength;
 * the IGM attenuation matters for the galaxies in the model, so that the above agreement is a test of how it is applied;
 * the template path, rather than the per-galaxy fallback, was the one used;
 * an SED which does not cover the blueshifted filter is reported as an error.

Andrew Benson, with assistance from Claude.
"""
import os
import re
import subprocess
import sys

import h5py
import numpy as np
import lxml.etree as ET
import scipy.constants as constants

# Tolerance on the agreement with the independent integration. Galacticus performs the filter integration adaptively to a
# relative tolerance of 4×10⁻³.
TOLERANCE_INDEPENDENT = 5.0e-3
# Minimum largest effect of IGM attenuation which the model must show for the comparison with the independent integration to
# test its application.
EFFECT_IGM_MINIMUM    = 4.0*TOLERANCE_INDEPENDENT

# This script runs with the working directory set to testSuite/, while the model is run from the repository root, so the two
# need different paths to the same file.
outputDirectory = "outputs/test-luminosities-from-SED-lightcone"
outputPath      = os.path.join(outputDirectory,"galacticus.hdf5"        )
coveragePath    = os.path.join(outputDirectory,"galacticusCoverage.hdf5")
os.makedirs(outputDirectory,exist_ok=True)
for path in ( outputPath, coveragePath ):
    if os.path.exists(path):
        os.remove(path)

def runModel(arguments,logName):
    """Run a model, returning its exit status and log."""
    status = subprocess.run("cd ..; ./Galacticus.exe "+arguments,shell=True,capture_output=True,stdin=subprocess.DEVNULL)
    log    = status.stdout.decode()+status.stderr.decode()
    with open(os.path.join(outputDirectory,logName),"w") as logFile:
        logFile.write(log)
    return status.returncode, log

print("Running model...")
status, log = runModel("testSuite/parameters/test-luminosities-from-SED-lightcone.xml","galacticus.log")
if status != 0 or re.search(r"fatal|aborted|Galacticus experienced an error in the GSL library",log,re.IGNORECASE) or not os.path.exists(outputPath):
    print("FAILED: model run:")
    print(log[-4000:])
    sys.exit(0)
print("SUCCESS: model run")

def unitsInSI(dataset):
    """Return the SI units of a dataset."""
    return dataset.attrs["units"]["unitsInSI"]

def filterResponse(name):
    """Read a filter response curve, or construct it for the special top-hat filters."""
    if name.startswith("adaptiveResolutionTopHat_"):
        # Constructed as in `Filter_Response_Load()`.
        wavelengthCentral, width = map(float,name.split("_")[1:3])
        wavelengths = wavelengthCentral+np.array([-width/2.0-width/100.0,-width/2.0,+width/2.0,+width/2.0+width/100.0])
        return wavelengths, np.array([0.0,1.0,1.0,0.0])
    document = ET.parse(os.path.join(os.environ["GALACTICUS_DATA_PATH"],"static","filters",name+".xml"))
    response = np.array([list(map(float,datum.text.split())) for datum in document.findall("response/datum")])
    return response[:,0], response[:,1]

def integrate(wavelengths,function):
    """Integrate a function, given as a callable, over a range of wavelengths with measure dλ/λ, subdividing every interval
    between the given wavelengths finely so that the result is limited only by the interpolation of the integrand."""
    countSubdivision = 64
    fine = np.concatenate([np.linspace(wavelengths[i],wavelengths[i+1],countSubdivision,endpoint=False) for i in range(len(wavelengths)-1)]+[wavelengths[-1:]])
    # `numpy.trapz` was renamed `numpy.trapezoid` in NumPy 2.0.
    trapezoid = np.trapezoid if hasattr(np,"trapezoid") else np.trapz
    return trapezoid(function(fine)/fine,fine)

# The IGM attenuation model of Inoue et al. (2014; MNRAS; 442; 1805). The coefficients of their Table 2 are read from the
# Galacticus source, while their fitting functions for the optical depth are implemented again here.
def inoue2014Coefficients(name,count):
    """Read a table of coefficients for the Inoue et al. (2014) model from the Galacticus source."""
    with open(os.path.join(os.environ["GALACTICUS_EXEC_PATH"],"source","stellar_populations","spectra","postprocess","Inoue2014.F90")) as sourceFile:
        source = sourceFile.read()
    table = re.search(name+r"=reshape\(.*?\[(.*?)\]",source,re.DOTALL).group(1)
    return np.array([float(value.replace("d","e")) for value in re.findall(r"[0-9.]+d[-+][0-9]+",table)]).reshape(-1,count)
coefficientsLAF = inoue2014Coefficients("aLAF",3)
coefficientsDLA = inoue2014Coefficients("aDLA",2)
# The Lyman series limit wavelength of hydrogen, including the reduced mass correction [Å].
wavelengthLymanLimit = 1.0e10/constants.Rydberg/(1.0-constants.m_e/(constants.m_p+constants.m_e))
wavelengthsLyman     = wavelengthLymanLimit/(1.0-1.0/np.arange(2,41)**2)

def inoue2014Attenuation(wavelength,redshift):
    """Return the IGM attenuation of Inoue et al. (2014) at the given rest-frame wavelengths [Å] for a source at the given redshift."""
    wavelengthObserved = wavelength*(1.0+redshift)
    opticalDepth       = np.zeros_like(wavelength)
    # Lyman series absorption by the Lyman-α forest and damped Lyman-α systems (equations 21 and 22).
    for wavelengthLyman, (a1LAF, a2LAF, a3LAF), (a1DLA, a2DLA) in zip(wavelengthsLyman,coefficientsLAF,coefficientsDLA):
        x       = wavelengthObserved/wavelengthLyman
        absorbs = (x >= 1.0) & (x <= 1.0+redshift)
        opticalDepth += np.where(absorbs,np.where(x < 2.2,a1LAF*x**1.2,np.where(x < 5.7,a2LAF*x**3.7,a3LAF*x**5.5)),0.0)
        opticalDepth += np.where(absorbs,np.where(x < 3.0,a1DLA*x**2.0,                              a2DLA*x**3.0 ),0.0)
    # Lyman continuum absorption (equations 25 and 28). Only the low-redshift forms are needed here.
    if redshift >= 1.2:
        print(f"FAILED: Lyman continuum absorption is implemented in this test only for redshifts below 1.2, but a galaxy has z={redshift:.4f}")
        sys.exit(0)
    x        = wavelengthObserved/wavelengthLymanLimit
    absorbs  = x < 1.0+redshift
    opticalDepth += np.where(absorbs,0.3248*(x**1.2-(1.0+redshift)**(-0.9)*x**2.1),0.0)
    opticalDepth += np.where(absorbs,0.2113*(1.0+redshift)**2.0-0.07661*(1.0+redshift)**2.3/x**0.3-0.1347*x**2.0,0.0)
    return np.exp(-opticalDepth)

# A top-hat filter spanning the Lyman limit in the rest frame of the galaxies, where IGM attenuation is substantial.
topHat = "adaptiveResolutionTopHat_1300_200"

with h5py.File(outputPath,"r") as model:
    try:
        nodes = model["Lightcone/Output1/nodeData"]
    except KeyError:
        print("FAILED: no lightcone output found")
        sys.exit(0)
    # Read the properties of the galaxies.
    try:
        redshift = nodes["redshift"][:]
        # Routes for computing luminosities: whether IGM attenuation is applied, the frame, and the filters.
        routes   = {
            "diskStellarSED:inoue2014:{filter}:observed": ( True , "observed", [ "Galex_FUV", "Galex_NUV",           topHat ] ),
            "diskStellarSED:sequence:{filter}:observed" : ( True , "observed", [ "Galex_FUV", "Galex_NUV",           topHat ] ),
            "diskStellarSED:identity:{filter}:observed" : ( False, "observed", [ "Galex_FUV", "Galex_NUV", "SDSS_r", topHat ] ),
            "diskStellarSED:identity:{filter}:rest"     : ( False, "rest"    , [ "Galex_FUV", "Galex_NUV", "SDSS_r", topHat ] )
        }
        luminosities = {
            route: {
                filter: nodes[route.format(filter=filter)][:]
                for filter in routes[route][2]
            }
            for route in routes
        }
        unitsLuminosity = unitsInSI(nodes["diskStellarSED:identity:SDSS_r:rest"])
        sed             = nodes["diskStellarSED:identity"]
        wavelengthsSED  = nodes["diskStellarSED:identityColumnValues"][:]
        unitsSED        = unitsInSI(sed)
        sed             = sed[:]
    except KeyError as e:
        print(f"FAILED: expected dataset missing from the output: {e}")
        sys.exit(0)
    # Find for which SEDs with redshift-dependent postprocessing a template was constructed. The template file name is recorded
    # in the parameters group of the SED property extractor which built it, along with the name of its postprocessor.
    templateFiles = {}
    def findTemplates(name,object):
        postprocessor = object.attrs.get("stellarPopulationSpectraPostprocessor") if isinstance(object,h5py.Group) else None
        if postprocessor in ( b"inoue2014", b"sequence" ):
            templateFiles.setdefault(postprocessor.decode(),[]).extend([object.attrs[attribute] for attribute in object.attrs if attribute.startswith("meta:sedMatrixFileName")])
    model["Parameters"].visititems(findTemplates)

countGalaxies = len(redshift)
selected      = luminosities["diskStellarSED:identity:{filter}:rest"]["SDSS_r"] > 0.0
countLuminous = int(selected.sum())
if countLuminous == 0:
    print("FAILED: no galaxy in the lightcone has a non-zero luminosity, so nothing is being tested")
    sys.exit(0)
print(f"SUCCESS: {countLuminous} of {countGalaxies} galaxies in the lightcone are luminous, spanning redshifts {redshift[selected].min():.3f} to {redshift[selected].max():.3f}")

# Templates must have been used for the SEDs with redshift-dependent postprocessing. Without the separation of that
# postprocessing into the filter there is no such template, and each galaxy is computed individually.
for postprocessor in ( "inoue2014", "sequence" ):
    if len(templateFiles.get(postprocessor,[])) == 0:
        print(f"FAILED: no SED template was constructed for the '{postprocessor}' postprocessor, so the per-galaxy fallback was used")
    else:
        print(f"SUCCESS: an SED template was used for the '{postprocessor}' postprocessor")

def compare(label,values,references,tolerance):
    """Compare two sets of luminosities, reporting the outcome."""
    difference = np.abs(values[selected]-references[selected])/references[selected]
    if np.any(~np.isfinite(difference)) or np.any(difference > tolerance):
        indexWorst = np.nanargmax(np.where(np.isfinite(difference),difference,np.inf))
        print(f"FAILED: {label}: maximum relative difference is {difference[indexWorst]:.3e} (tolerance {tolerance:.1e}) at z={redshift[selected][indexWorst]:.4f}")
    else:
        print(f"SUCCESS: {label}: maximum relative difference is {difference.max():.3e} (tolerance {tolerance:.1e})")

# Integrate the output rest-frame SED through the filters independently. The luminosities are in units of the AB zero point,
# the SED in units of Solar luminosities per Hz, and the ratio of those gives the AB zero point in the units of the SED.
zeroPointAB = unitsLuminosity/unitsSED
filters     = sorted(set(sum([routes[route][2] for route in routes],[])))
responses   = {filter: filterResponse(filter) for filter in filters}
for route, ( attenuated, frame, filtersRoute ) in routes.items():
    for filter in filtersRoute:
        wavelengthsFilter, responseFilter = responses[filter]
        normalization = zeroPointAB*integrate(wavelengthsFilter,lambda wavelength: np.interp(wavelength,wavelengthsFilter,responseFilter))
        independent   = np.zeros(countGalaxies)
        for i in np.nonzero(selected)[0]:
            # In the observed frame the filter is blueshifted into the rest frame of the galaxy.
            expansionFactor = 1.0/(1.0+redshift[i]) if frame == "observed" else 1.0
            wavelengths     = np.union1d(wavelengthsFilter*expansionFactor,wavelengthsSED)
            wavelengths     = wavelengths[(wavelengths >= wavelengthsFilter[0]*expansionFactor) & (wavelengths <= wavelengthsFilter[-1]*expansionFactor)]
            attenuation     = (lambda wavelength: inoue2014Attenuation(wavelength,redshift[i])) if attenuated else (lambda wavelength: 1.0)
            independent[i]  = integrate(
                wavelengths,
                lambda wavelength: np.interp(wavelength/expansionFactor,wavelengthsFilter,responseFilter)*np.interp(wavelength,wavelengthsSED,sed[i])*attenuation(wavelength)
            )/normalization
        compare(f"{route.format(filter=filter)} against an independent integration of the SED",luminosities[route][filter],independent,TOLERANCE_INDEPENDENT)

# The IGM attenuation must have an effect large enough that the above comparison tests its application.
for filter in routes["diskStellarSED:inoue2014:{filter}:observed"][2]:
    effectIGM = np.max(1.0-luminosities["diskStellarSED:inoue2014:{filter}:observed"][filter][selected]/luminosities["diskStellarSED:identity:{filter}:observed"][filter][selected])
    print(f"   IGM attenuation reduces the {filter} luminosity by up to {effectIGM:.3e}")
    if filter == topHat:
        if effectIGM < EFFECT_IGM_MINIMUM:
            print(f"FAILED: IGM attenuation reduces the top-hat luminosity by at most {effectIGM:.3e}, too little for the comparison with the independent integration to test its application (minimum {EFFECT_IGM_MINIMUM:.1e})")
        else:
            print(f"SUCCESS: IGM attenuation reduces the top-hat luminosity by up to {effectIGM:.3f}")

# An SED which does not cover the blueshifted filter must be reported as an error.
print("Running model with insufficient SED wavelength coverage...")
status, log = runModel("testSuite/parameters/test-luminosities-from-SED-lightcone.xml testSuite/parameters/test-luminosities-from-SED-lightcone-coverage.xml","galacticusCoverage.log")
if status == 0:
    print("FAILED: a model whose SED does not cover the blueshifted filter ran without error")
elif not re.search(r"does not cover the wavelength range",log):
    print("FAILED: a model whose SED does not cover the blueshifted filter failed, but not with the expected error:")
    print(log[-4000:])
else:
    print("SUCCESS: an SED which does not cover the blueshifted filter is reported as an error")

sys.exit(0)
