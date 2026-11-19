#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize  import brentq, minimize_scalar
from scipy.special   import gamma, gammainc, gammaincc

# Check ram pressure stripping of the gas of a satellite's disk by the `koppen2018` ramPressureStripping class against an
# independent calculation.
#
# A single satellite, carrying a gaseous disk, orbits on an eccentric orbit in a static host which carries a hot halo (see
# testSuite/parameters/ramPressureKoppen2018.xml). The gas disk is a test mass, and the satellite's own hot halo is negligible, so
# that the gravitational potential of the satellite is that of its NFW halo alone. The `outputTrajectory` node operator records,
# after every accepted ODE step, the maximum, P, and time integral, I, of the ram pressure experienced by the satellite, the scale
# length, R_d, and gas mass of its disk, and the degree of truncation of the gas, x_gas.
#
# At each record the stripping radius, r_s, is found here from the criterion of Köppen et al. (2018), using the NFW potential, and
# the surface density of the untruncated exponential disk (which is that of the truncated disk within its truncation radius, and so
# at r_s whenever stripping occurs). This gives the degree of truncation, x_s, of a disk truncated
# at r_s, and the period of vertical oscillations, T_vert, at r_s. The equation
#
#     dx_gas/dt = max(0, x_s - x_gas)/T_vert
#
# is then integrated here, using SciPy's ODE solver with x_s and T_vert interpolated between records, and compared with x_gas
# recorded by the model. Since no other process changes the gas mass, the gas mass must also equal M_0 exp(-x_gas), where M_0 is the
# initial gas mass.
#
# The model is run with no decay of P and I ([efficiencyDecay]=0), and with decay ([efficiencyDecay]=1). With decay, r_s computed
# from the untruncated disk may exceed the truncation radius of the gas, in which case the gas is already truncated within r_s, so
# x_s < x_gas and no stripping occurs, as for the truncated disk in the model.
#
# Gas arriving in the disk is assumed to have the untruncated distribution, so that the untruncated mass of the gas, M_gas exp(x_gas),
# grows by the mass of gas which arrives. This is checked for gas which cools onto the disk, and for gas which arrives in a merger
# (from a sub-satellite which merges with the satellite after its gas has been partly stripped).
#
# Andrew Benson, Claude

# Event identifier for records made after each accepted ODE step, as enumerated in
# source/nodes/operators/output/trajectory_events.F90.
eventDifferentialEvolutionPostStep = 1

# Physical constants and units, matching Galacticus' definitions (from the GSL, and the sidereal year).
gravitationalConstantSI = 6.673e-11
massSolar               = 1.98892e30
megaParsec              = 3.08567758135e22
gigaYear                = 3.15581497635456e16
G                       = gravitationalConstantSI*massSolar/1.0e6/megaParsec # (km/s)² Mpc/M☉.
MpcPerKmPerSToGyr       = megaParsec/1.0e3/gigaYear

# Properties of the satellite (see testSuite/parameters/ramPressureKoppen2018Tree.xml) and of the model.
massVirial                    = 1.0e10
radiusScale                   = 0.0057
massGasInitial                = 1.0e6
densityContrast               = 329.6207717164546
HubbleConstant                = 67.36
OmegaMatter                   = 0.3153
OmegaBaryon                   = 0.0493
radiusDecayOverRadiusTruncate = 0.1
truncationLarge               = 25.0

# Tolerances, set from the measured agreement with a factor of ten or so in hand. The degree of truncation agrees with the reference
# to 3.6e-4 (limited by the interpolation of x_s and T_vert between records, and by the self-gravity of the gas disk, which is
# neglected here), and the gas mass with M_0 exp(-x_gas) to 8e-11 (limited by the ODE solver).
toleranceTruncation = 4.0e-3
toleranceMass       = 1.0e-9
# The enclosed mass of the satellite includes its gas disk (of mass 10⁻⁴ that of the halo, but a larger fraction of the halo mass
# enclosed within small radii), which is neglected in the comparison with the NFW profile. They agree to 4.8e-4.
toleranceProfile    = 5.0e-3
# Tolerance for the growth of the untruncated mass of the gas by cooling (limited by the ODE solver; the measured agreement is 2e-11).
toleranceMassCooling = 1.0e-9
# Tolerance for the growth of the untruncated mass of the gas by a merger.
toleranceMassMerger  = 1.0e-9

# Mass of gas (M☉) in the disk of a sub-satellite which merges with the satellite, and the time (Gyr) at which it merges.
massGasMerger = 3.0e5
timeMerger    = 4.0

# Virial radius and NFW normalization of the satellite. The dark matter of the halo is a fraction 1-Ω_b/Ω_M of its total mass.
densityCritical = 3.0*HubbleConstant**2/8.0/np.pi/G
radiusVirial    = (3.0*massVirial/4.0/np.pi/densityContrast/OmegaMatter/densityCritical)**(1.0/3.0)
concentration   = radiusVirial/radiusScale
massNormalization = massVirial*(1.0-OmegaBaryon/OmegaMatter)/(np.log(1.0+concentration)-concentration/(1.0+concentration))

def massEnclosed(radius):
    x = radius/radiusScale
    return massNormalization*(np.log(1.0+x)-x/(1.0+x))

def potential(radius):
    return -G*massNormalization*np.log(1.0+radius/radiusScale)/radius

def accelerationVerticalMaximum(radius, radiusDisk):
    """Maximum over height of the vertical acceleration, in (km/s)²/Mpc, found over the same range of heights as Galacticus."""
    heightMinimum = 1.0e-3*radiusDisk
    heightMaximum = max(10.0*radiusDisk, 3.0*radius)
    def accelerationNegative(heightLogarithmic):
        height = np.exp(heightLogarithmic)
        r      = np.sqrt(radius**2+height**2)
        return -G*massEnclosed(r)*height/r**3
    # Locate the maximum on a coarse grid, then refine.
    grid  = np.linspace(np.log(heightMinimum), np.log(heightMaximum), 101)
    i     = np.argmin([accelerationNegative(h) for h in grid])
    lower = grid[max(i-1, 0)]
    upper = grid[min(i+1, len(grid)-1)]
    result = minimize_scalar(accelerationNegative, bounds=(lower, upper), method="bounded", options={"xatol": 1.0e-8})
    return -result.fun

def fractionRetained(radiusTruncate):
    """Fraction of mass retained by an exponential disk truncated at the given radius (in units of its scale length), computed from
    the radial moment, I₁, of the truncated surface density."""
    b = 1.0/radiusDecayOverRadiusTruncate
    a = 1.0+b/radiusTruncate
    inner = gammainc(2.0, radiusTruncate)
    outer = radiusTruncate**(-b)*np.exp(b)*a**(-(b+2.0))*gammaincc(b+2.0, a*radiusTruncate)*gamma(b+2.0)
    return inner+outer

def criterionMargin(radius, radiusDisk, forceMaximum, momentumIntegrated):
    """The margin by which the Köppen et al. (2018) criterion is satisfied (positive where gas is stripped)."""
    surfaceDensity = massGasInitial/2.0/np.pi/radiusDisk**2*np.exp(-radius/radiusDisk)
    forceRestoring = surfaceDensity*accelerationVerticalMaximum(radius, radiusDisk)
    velocitySquared = max(-2.0*potential(radius)-G*massEnclosed(radius)/radius, 0.0)
    momentumEscape  = surfaceDensity*np.sqrt(velocitySquared)
    if momentumIntegrated > momentumEscape:
        return (1.0-momentumEscape/momentumIntegrated)-forceRestoring/forceMaximum*(1.0-0.5*momentumEscape/momentumIntegrated)
    else:
        return (1.0-momentumEscape/momentumIntegrated)-forceRestoring/forceMaximum*0.5

def strippingTarget(radiusDisk, forceMaximum, forceIntegrated):
    """Return the degree of truncation, x_s, toward which the gas relaxes, and the timescale (Gyr) for it to do so."""
    if forceMaximum <= 0.0 or forceIntegrated <= 0.0:
        return 0.0, np.inf
    # Convert the time-integrated ram pressure to units of momentum per unit area.
    momentumIntegrated = forceIntegrated/MpcPerKmPerSToGyr
    radiusMinimum = 1.0e-2*radiusDisk
    radiusMaximum = 3.0e+1*radiusDisk
    marginMinimum = criterionMargin(radiusMinimum, radiusDisk, forceMaximum, momentumIntegrated)
    marginMaximum = criterionMargin(radiusMaximum, radiusDisk, forceMaximum, momentumIntegrated)
    if marginMaximum < 0.0:
        return 0.0, np.inf
    elif marginMinimum >= 0.0:
        radiusStripping = radiusMinimum
        truncation      = truncationLarge
    else:
        radiusStripping = brentq(criterionMargin, radiusMinimum, radiusMaximum, args=(radiusDisk, forceMaximum, momentumIntegrated), xtol=1.0e-12, rtol=1.0e-10)
        truncation      = min(-np.log(max(fractionRetained(radiusStripping/radiusDisk), np.exp(-truncationLarge))), truncationLarge)
    # Period of vertical oscillations. For a spherical potential ∂²Φ/∂z² at z=0 is GM(<R)/R³.
    periodVertical = 2.0*np.pi/np.sqrt(G*massEnclosed(radiusStripping)/radiusStripping**3)*MpcPerKmPerSToGyr
    return truncation, periodVertical

def check(efficiencyDecay):
    """Run the model with the given decay efficiency of the ram pressure history, and check the stripping of the gas."""
    label = f"efficiencyDecay={efficiencyDecay}"
    with open("parameters/ramPressureKoppen2018.xml") as file:
        parameters = file.read()
    parameters = parameters.replace('<efficiencyDecay value="0.0"/>', f'<efficiencyDecay value="{efficiencyDecay}"/>')
    parameters = parameters.replace("testSuite/outputs/ramPressureKoppen2018.hdf5", f"testSuite/outputs/ramPressureKoppen2018_{efficiencyDecay}.hdf5")
    with open(f"outputs/ramPressureKoppen2018_{efficiencyDecay}.xml", "w") as file:
        file.write(parameters)
    status = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/outputs/ramPressureKoppen2018_{efficiencyDecay}.xml", shell=True)
    if status.returncode != 0:
        print(f"FAILED: model run ({label})")
        return
    print(f"SUCCESS: model run ({label})")

    # Read the trajectory of the satellite.
    try:
        model     = h5py.File(f"outputs/ramPressureKoppen2018_{efficiencyDecay}.hdf5", "r")
        nodeData  = model["nodeTrajectories/Trajectory1/nodeData"]
        selection = nodeData["trajectoryEvent"][:] == eventDifferentialEvolutionPostStep
        data      = {name: nodeData[name][:][selection] for name in ("time", "diskRadius", "diskMassGas", "diskTruncationGas", "satelliteRamPressureMaximum", "satelliteRamPressureTimeIntegrated", "massProfile")}
    except (OSError, KeyError) as error:
        print(f"FAILED: unable to read the model trajectory ({label}: {error})")
        return
    order = np.argsort(data["time"])
    for name in data:
        data[name] = data[name][order]

    # Check that the enclosed mass profile of the satellite is that of the NFW halo assumed here.
    fractionsProfile = [0.05, 0.2, 1.0]
    errorProfile = max(np.max(np.abs(data["massProfile"][:, i]/massEnclosed(fractionsProfile[i]*radiusVirial)-1.0)) for i in range(3))
    if errorProfile > toleranceProfile:
        print(f"FAILED: enclosed mass of the satellite differs from that of the assumed NFW halo by up to {errorProfile:.3e} ({label})")
    else:
        print(f"SUCCESS: enclosed mass of the satellite agrees with that of the assumed NFW halo (to {errorProfile:.3e}; {label})")

    # Check that the gas mass is consistent with the degree of truncation.
    errorMass = np.max(np.abs(data["diskMassGas"]/(massGasInitial*np.exp(-data["diskTruncationGas"]))-1.0))
    if errorMass > toleranceMass:
        print(f"FAILED: gas mass differs from M_0 exp(-x_gas) by up to {errorMass:.3e} ({label})")
    else:
        print(f"SUCCESS: gas mass agrees with M_0 exp(-x_gas) (to {errorMass:.3e}; {label})")

    # Compute the target degree of truncation and timescale at each record.
    targets = np.array([strippingTarget(data["diskRadius"][i], data["satelliteRamPressureMaximum"][i], data["satelliteRamPressureTimeIntegrated"][i]) for i in range(len(data["time"]))])
    truncationTarget = targets[:, 0]
    rateTarget       = np.where(np.isfinite(targets[:, 1]), 1.0/targets[:, 1], 0.0)

    # Integrate the evolution of the degree of truncation.
    def rate(time, truncation):
        target = np.interp(time, data["time"], truncationTarget)
        return [max(target-truncation[0], 0.0)*np.interp(time, data["time"], rateTarget)]
    solution = solve_ivp(rate, (data["time"][0], data["time"][-1]), [data["diskTruncationGas"][0]], t_eval=data["time"], rtol=1.0e-8, atol=1.0e-10, max_step=np.min(np.diff(data["time"])[np.diff(data["time"]) > 0]))
    truncationReference = solution.y[0]

    errorTruncation = np.max(np.abs(data["diskTruncationGas"]-truncationReference)/np.maximum(truncationReference, 1.0))
    print(f"Final degree of truncation of the gas ({label}): model {data['diskTruncationGas'][-1]:.6f}, reference {truncationReference[-1]:.6f} (target {truncationTarget[-1]:.6f})")
    if data["diskTruncationGas"][-1] <= 0.0:
        print(f"FAILED: the gas of the disk was not stripped ({label})")
    elif errorTruncation > toleranceTruncation:
        print(f"FAILED: degree of truncation of the gas differs from the reference by up to {errorTruncation:.3e} ({label})")
    else:
        print(f"SUCCESS: degree of truncation of the gas agrees with the reference (to {errorTruncation:.3e}; {label})")

def checkCooling():
    """Run the model with gas cooling onto the disk from the satellite's own hot halo, and check that the untruncated mass of the
    gas, M_gas exp(x_gas), grows by the mass of gas which has cooled, as required by the assumption that arriving gas has the
    untruncated distribution."""
    label = "with cooling"
    with open("parameters/ramPressureKoppen2018Tree.xml") as file:
        tree = file.read()
    tree = tree.replace("<mass>1.0e3</mass>", "<mass>1.0e6</mass>").replace("<angularMomentum>1.0e3</angularMomentum>", "<angularMomentum>1.0e6</angularMomentum>")
    with open("outputs/ramPressureKoppen2018TreeCooling.xml", "w") as file:
        file.write(tree)
    with open("parameters/ramPressureKoppen2018.xml") as file:
        parameters = file.read()
    parameters = parameters.replace("testSuite/parameters/ramPressureKoppen2018Tree.xml", "testSuite/outputs/ramPressureKoppen2018TreeCooling.xml")
    parameters = parameters.replace("testSuite/outputs/ramPressureKoppen2018.hdf5", "testSuite/outputs/ramPressureKoppen2018Cooling.hdf5")
    parameters = parameters.replace(
        '    <nodeOperator value="ramPressureHistory"               >',
        '    <nodeOperator value="massCooled"/>\n'
        '    <nodeOperator value="CGMCoolingHeating"><component value="disk"/></nodeOperator>\n'
        '    <nodeOperator value="ramPressureHistory"               >'
    )
    # The cooling rate is declared at the top level, as it is used by both the cooling and the cooled mass node operators.
    parameters = parameters.replace("</parameters>", '<coolingRate value="simple"><timeScale value="3.0"/></coolingRate>\n</parameters>')
    parameters = parameters.replace('<nodePropertyExtractor value="diskTruncationGas" />', '<nodePropertyExtractor value="diskTruncationGas" />\n<nodePropertyExtractor value="massCooled"/>')
    with open("outputs/ramPressureKoppen2018Cooling.xml", "w") as file:
        file.write(parameters)
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/outputs/ramPressureKoppen2018Cooling.xml", shell=True)
    if status.returncode != 0:
        print(f"FAILED: model run ({label})")
        return
    print(f"SUCCESS: model run ({label})")
    try:
        model     = h5py.File("outputs/ramPressureKoppen2018Cooling.hdf5", "r")
        nodeData  = model["nodeTrajectories/Trajectory1/nodeData"]
        selection = nodeData["trajectoryEvent"][:] == eventDifferentialEvolutionPostStep
        data      = {name: nodeData[name][:][selection] for name in ("time", "diskMassGas", "diskTruncationGas", "massCooledCGM")}
    except (OSError, KeyError) as error:
        print(f"FAILED: unable to read the model trajectory ({label}: {error})")
        return
    massUntruncated = data["diskMassGas"]*np.exp(np.maximum(data["diskTruncationGas"], 0.0))
    errorMass       = np.max(np.abs(massUntruncated/(massGasInitial+data["massCooledCGM"])-1.0))
    print(f"Cooled mass {np.max(data['massCooledCGM']):.4e} M☉; maximum degree of truncation of the gas {np.max(data['diskTruncationGas']):.4f}, final {data['diskTruncationGas'][np.argmax(data['time'])]:.4f}")
    if np.max(data["massCooledCGM"]) < 0.1*massGasInitial or np.max(data["diskTruncationGas"]) <= 0.0:
        print(f"FAILED: too little gas cooled, or the gas was not stripped ({label})")
    elif errorMass > toleranceMassCooling:
        print(f"FAILED: untruncated gas mass, M_gas exp(x_gas), differs from M_0 plus the cooled mass by up to {errorMass:.3e} ({label})")
    else:
        print(f"SUCCESS: untruncated gas mass, M_gas exp(x_gas), equals M_0 plus the cooled mass (to {errorMass:.3e}; {label})")

def checkMerger():
    """Run the model with a sub-satellite, carrying a gaseous disk, which merges with the satellite after the gas of the
    satellite's disk has been partly stripped, and with the gas of the merging galaxy sent to the disk of the satellite. Check that
    the untruncated mass of the gas of the satellite's disk, M_gas exp(x_gas), grows by the mass of gas which arrives in the merger,
    as required by the assumption that arriving gas has the untruncated distribution."""
    label = "with a merger"
    with open("parameters/ramPressureKoppen2018Tree.xml") as file:
        tree = file.read()
    # Make the sub-satellite a satellite of the satellite.
    nodeSatellite = "<index>3</index>\n    <parent>1</parent>\n    <firstChild>-1</firstChild>\n    <sibling>-1</sibling>\n"
    if nodeSatellite not in tree:
        print(f"FAILED: unable to find the satellite in the merger tree ({label})")
        return
    tree = tree.replace(nodeSatellite, nodeSatellite+"    <firstSatellite>4</firstSatellite>\n")
    # The sub-satellite, on an orbit close to the center of the satellite. It merges at a preset time, set by the
    # `satellite` timestep, after the gas of the satellite has begun to be stripped. Its gas disk, of mass ΔM, is compact, so that it
    # is not stripped by the (negligible) hot halo of the satellite.
    subSatellite = f"""
  <node>
    <index>4</index>
    <parent>3</parent>
    <firstChild>-1</firstChild>
    <sibling>-1</sibling>
    <basic>
      <time>1.0</time>
      <timeLastIsolated>1.0</timeLastIsolated>
      <mass>1.0e8</mass>
    </basic>
    <darkMatterProfile>
      <scale>0.001</scale>
    </darkMatterProfile>
    <hotHalo>
      <isInitialized>true</isInitialized>
      <mass>1.0e1</mass>
      <outerRadius>0.01</outerRadius>
      <angularMomentum>1.0e1</angularMomentum>
    </hotHalo>
    <disk>
      <massStellar>0.0</massStellar>
      <massGas>{massGasMerger}</massGas>
      <abundancesStellar>
	<metals>0.0</metals>
      </abundancesStellar>
      <abundancesGas>
	<metals>0.0</metals>
      </abundancesGas>
      <angularMomentum>{massGasMerger*1.0e-3}</angularMomentum>
    </disk>
    <satellite>
      <boundMass>1.0e8</boundMass>
      <tidalHeatingNormalized>0.0</tidalHeatingNormalized>
      <position>0.005</position>
      <position>0.0</position>
      <position>0.0</position>
      <velocity>0.0</velocity>
      <velocity>60.0</velocity>
      <velocity>0.0</velocity>
      <tidalTensorPathIntegrated>
        <x00>0.0</x00>
        <x01>0.0</x01>
        <x02>0.0</x02>
        <x11>0.0</x11>
        <x12>0.0</x12>
        <x22>0.0</x22>
      </tidalTensorPathIntegrated>
      <timeOfMerging>{timeMerger}</timeOfMerging>
      <destructionTime>-1.0</destructionTime>
      <virialOrbit>
        <radius>0.005</radius>
        <velocityRadial>0.0</velocityRadial>
        <velocityTangential>60.0</velocityTangential>
        <massHost>1.0e10</massHost>
        <massSatellite>1.0e8</massSatellite>
      </virialOrbit>
    </satellite>
  </node>
</tree>"""
    tree = tree.replace("</tree>", subSatellite)
    with open("outputs/ramPressureKoppen2018TreeMerger.xml", "w") as file:
        file.write(tree)
    with open("parameters/ramPressureKoppen2018.xml") as file:
        parameters = file.read()
    parameters = parameters.replace("testSuite/parameters/ramPressureKoppen2018Tree.xml", "testSuite/outputs/ramPressureKoppen2018TreeMerger.xml")
    parameters = parameters.replace("testSuite/outputs/ramPressureKoppen2018.hdf5", "testSuite/outputs/ramPressureKoppen2018Merger.hdf5")
    # Send all gas and stars of the merging galaxy to the disk of the satellite, with every merger treated as minor.
    parameters = parameters.replace(
        '<mergerRemnantSize          value="null"/>',
        '<mergerRemnantSize          value="null"/>\n'
        '  <mergerMassMovements value="simple">\n'
        '    <destinationGasMinorMerger   value="disk"/>\n'
        '    <destinationStarsMinorMerger value="disk"/>\n'
        '    <massRatioMajorMerger        value="1.0e3"/>\n'
        '  </mergerMassMovements>'
    )
    with open("outputs/ramPressureKoppen2018Merger.xml", "w") as file:
        file.write(parameters)
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/outputs/ramPressureKoppen2018Merger.xml", shell=True)
    if status.returncode != 0:
        print(f"FAILED: model run ({label})")
        return
    print(f"SUCCESS: model run ({label})")
    try:
        model     = h5py.File("outputs/ramPressureKoppen2018Merger.hdf5", "r")
        nodeData  = model["nodeTrajectories/Trajectory1/nodeData"]
        selection = nodeData["trajectoryEvent"][:] == eventDifferentialEvolutionPostStep
        data      = {name: nodeData[name][:][selection] for name in ("time", "diskMassGas", "diskTruncationGas")}
    except (OSError, KeyError) as error:
        print(f"FAILED: unable to read the model trajectory ({label}: {error})")
        return
    order = np.argsort(data["time"])
    for name in data:
        data[name] = data[name][order]
    before = data["time"] < timeMerger
    after  = data["time"] > timeMerger
    if not np.any(before) or not np.any(after):
        print(f"FAILED: no records before and after the merger ({label})")
        return
    truncationBefore = data["diskTruncationGas"][before][-1]
    truncationAfter  = data["diskTruncationGas"][after ][ 0]
    massBefore       = data["diskMassGas"      ][before][-1]
    massAfter        = data["diskMassGas"      ][after ][ 0]
    print(f"Across the merger ({label}): gas mass {massBefore:.6e} → {massAfter:.6e} M☉, degree of truncation of the gas {truncationBefore:.6f} → {truncationAfter:.6f}")
    massUntruncated = data["diskMassGas"]*np.exp(np.maximum(data["diskTruncationGas"], 0.0))
    errorBefore = np.max(np.abs(massUntruncated[before]/ massGasInitial               -1.0))
    errorAfter  = np.max(np.abs(massUntruncated[after ]/(massGasInitial+massGasMerger)-1.0))
    if truncationBefore < 0.1 or truncationAfter >= truncationBefore or abs(massAfter-massBefore-massGasMerger) > 1.0e-6*massGasMerger:
        print(f"FAILED: the gas was not stripped before the merger, or the merger did not add gas to the disk and reduce its truncation ({label})")
    elif max(errorBefore, errorAfter) > toleranceMassMerger:
        print(f"FAILED: untruncated gas mass, M_gas exp(x_gas), differs from M_0 (before the merger; by up to {errorBefore:.3e}) or M_0+ΔM (after; by up to {errorAfter:.3e}) ({label})")
    else:
        print(f"SUCCESS: untruncated gas mass, M_gas exp(x_gas), equals M_0 before the merger and M_0+ΔM after it (to {max(errorBefore, errorAfter):.3e}; {label})")

# Without decay, P and I never decrease, so the gas is stripped progressively. With decay, the gas is stripped near each
# pericentric passage, and not otherwise.
for efficiencyDecay in (0.0, 1.0):
    check(efficiencyDecay)
checkCooling()
checkMerger()
