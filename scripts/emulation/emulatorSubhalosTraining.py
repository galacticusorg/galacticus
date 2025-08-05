# Train a normalizing flow emulator on Galacticus subhalo populations.
import argparse
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import regularizers
import numpy as np
import matplotlib.pyplot as plt
import tensorflow_probability as tfp
import h5py
from scipy import stats
from emulator_utils import norm_transform, norm_transform_inv, Coupling, RealNVP

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='emulatorSubhalosTraining.py',description='Train an emulator on Galacticus subhalo populations.')
parser.add_argument('filename',nargs='+')
parser.add_argument('--outputDirectory', action='store',default='.',help='the directory to which emulator weights (and plots) should be output')
args = parser.parse_args()

# Read all training data from file.
## Initialize arrays that will accumulate training data.
countSubhalosMean_             = np.array([])
weights                        = np.array([])
massInfallNormalized           = np.array([])
concentrationNormalized        = np.array([])
massBoundNormalized            = np.array([])
redshiftLastIsolatedNormalized = np.array([])
truncationRadiusNormalized     = np.array([])
projectedRadiusNormalized      = np.array([])
# Readi in Galacticus data from each file.
for fileName in args.filename:
    model                      = h5py.File(fileName, 'r')
    mergerTreeConstructorGroup = model['Parameters/mergerTreeConstructor/mergerTreeConstructor']
    mergerTreeBuildMassesGroup = model['Parameters/mergerTreeBuildMasses'                      ]
    massResolutionGroup        = model['Parameters/mergerTreeMassResolution'                   ]
    outputGroup                = model['Outputs/Output1/nodeData'                              ]
    # Get masses and counts of trees.
    massTree              = mergerTreeBuildMassesGroup.attrs['massTree'      ][0]
    countTree             = mergerTreeBuildMassesGroup.attrs['treeCount'     ][0]
    redshiftTree          = mergerTreeConstructorGroup.attrs['redshiftBase'  ]
    massResolution        = massResolutionGroup       .attrs['massResolution']
    # Read all required properties.
    weight                = outputGroup['nodeSubsamplingWeight'       ][:]
    treeIndex             = outputGroup['mergerTreeIndex'             ][:]
    isCentral             = outputGroup['nodeIsIsolated'              ][:]
    massInfall            = outputGroup['massHaloEnclosedCurrent'     ][:]
    massBound             = outputGroup['satelliteBoundMass'          ][:]
    concentration         = outputGroup['concentration'               ][:]
    truncationRadius      = outputGroup['radiusTidalTruncationNFW'    ][:]
    scaleRadius           = outputGroup['darkMatterProfileScale'      ][:]
    redshiftLastIsolated  = outputGroup['redshiftLastIsolated'        ][:]
    positionOrbitalX      = outputGroup['positionOrbitalX'            ][:]
    positionOrbitalY      = outputGroup['positionOrbitalY'            ][:]
    positionOrbitalZ      = outputGroup['positionOrbitalZ'            ][:]
    radiusVirial          = outputGroup['darkMatterOnlyRadiusVirial'  ][:]
    velocityVirial        = outputGroup['darkMatterOnlyVelocityVirial'][:]
    # Compute derived quantities.
    radiusProjected       = np.sqrt(positionOrbitalX**2+positionOrbitalY**2                    )
    radiusOrbital         = np.sqrt(positionOrbitalX**2+positionOrbitalY**2+positionOrbitalZ**2)
    # Identify centrals and subhalos.
    subhalos              = (isCentral == 0) & (massInfall > 2.0*massResolution)
    centrals              = (isCentral == 1)
    # Count subhalos per tree, and find the mean.
    countSubhalos         = np.zeros(countTree)
    for j in range(countTree):
        selectTree       = (isCentral == 0) & (massInfall > 2.0*massResolution) & (treeIndex == j+1)
        countSubhalos[j] = np.sum(weight[selectTree])
    # Fix any overmassive subhalos to their infall mass. This can happen (rarely) because of the choice of halo mass definition.
    overMassive            = (massBound > massInfall)
    massBound[overMassive] = massInfall[overMassive]
    # Compute normalized quantities.
    radiusVirialHost_               = radiusVirial[centrals][0]
    massInfallNormalized_           = np.log10(massInfall          [subhalos]/massTree                   )
    massBoundNormalized_            = np.log10(massBound           [subhalos]/massInfall       [subhalos])
    concentrationNormalized_        =          concentration       [subhalos]
    redshiftLastIsolatedNormalized_ =          redshiftLastIsolated[subhalos]
    radiusOrbitalNormalized_        = np.log10(radiusOrbital       [subhalos]/radiusVirialHost_          )
    truncationRadiusNormalized_     = np.log10(truncationRadius    [subhalos]/radiusVirialHost_          )
    projectedRadiusNormalized_      = np.log10(radiusProjected     [subhalos]/radiusVirialHost_          )
    weight_                         = weight[subhalos]
    # Accumulate normalized quantities.
    countSubhalosMean_              = np.append(countSubhalosMean_            , countSubhalos                  )
    weights                         = np.append(weights                       , weight_                        )
    massInfallNormalized            = np.append(massInfallNormalized          , massInfallNormalized_          )
    concentrationNormalized         = np.append(concentrationNormalized       , concentrationNormalized_       )
    massBoundNormalized             = np.append(massBoundNormalized           , massBoundNormalized_           )
    redshiftLastIsolatedNormalized  = np.append(redshiftLastIsolatedNormalized, redshiftLastIsolatedNormalized_)
    truncationRadiusNormalized      = np.append(truncationRadiusNormalized    , truncationRadiusNormalized_    )
    projectedRadiusNormalized       = np.append(projectedRadiusNormalized     , projectedRadiusNormalized_     )

# Find the number number of subhalos per tree.
countSubhalosMean = np.mean(countSubhalosMean_)
    
# Pack all subhalo data into a 6D array.
data = np.array(
    list(
        zip(
            massInfallNormalized          ,
            concentrationNormalized       ,
            massBoundNormalized           ,
            redshiftLastIsolatedNormalized,
            truncationRadiusNormalized    ,
            projectedRadiusNormalized
        )
    )
)

# Shift the data into a unit cube, append the weights, and shuffle.
data_min, data_max, dataNormalized = norm_transform(data,-1,1)
dataNormalizedWeighted             = np.hstack((dataNormalized, np.expand_dims(weights,1)))
np.random.shuffle(dataNormalizedWeighted)

# Store normalization data to a file. This will be used to un-normalize the emulator output during inference.
minArray = np.nanmin(data, axis = 0)
maxArray = np.nanmax(data, axis = 0)
with h5py.File(args.outputDirectory+'/normalizationData.hdf5', 'w') as normalizationFile:
    normalizationFile.create_dataset("minArray",data=minArray)
    normalizationFile.create_dataset("maxArray",data=maxArray)
    normalizationFile.attrs['massTree'         ] = massTree
    normalizationFile.attrs['redshiftTree'     ] = redshiftTree
    normalizationFile.attrs['radiusVirialHost' ] = radiusVirialHost_
    normalizationFile.attrs['massResolution'   ] = massResolution
    normalizationFile.attrs['countSubhalosMean'] = countSubhalosMean

# Train the emulator.
model = RealNVP(num_coupling_layers=12)
model.compile(optimizer=keras.optimizers.Adam(learning_rate=0.0001))
history = model.fit(
    dataNormalizedWeighted, batch_size=256, epochs=50, verbose=2, validation_split=0.2
)
model.save_weights(args.outputDirectory+'/emulator.weights.h5')

# Build a new emulator from the training weights.
emulator = RealNVP(num_coupling_layers=12)
emulator.build(input_shape= (6))
emulator.load_weights(args.outputDirectory+'/emulator.weights.h5')

# Specify the number of subhalos to sample.
countSample = 500

# Generate a sample of subhalos from the emulator, an unnormalize.
samples = emulator.distribution.sample(countSample)
x, _    = emulator.predict(samples)
xt      = norm_transform_inv(x, minArray, maxArray, -1, 1)

# Extract to named arrays for convenience.
massInfallNormalized           = xt[:,0]
concentrationNormalized        = xt[:,1]
massBoundNormalized            = xt[:,2]
redshiftLastIsolatedNormalized = xt[:,3]
truncationRadiusNormalized     = xt[:,4]
projectedRadiusNormalized      = xt[:,5]

# Create a filter to impose physical constraints on emulator halos.
filter =                                                                                          \
    (massInfallNormalized           >   np.log10(2.0*massResolution/massTree)                 ) & \
    (massBoundNormalized            <   np.log10(1.0e9/(massTree*10.0**massInfallNormalized ))) & \
    (massBoundNormalized            <=  0.0                                                   ) & \
    (massBoundNormalized            >  -massInfallNormalized+np.log10(massResolution/massTree)) & \
    (redshiftLastIsolatedNormalized >=  redshiftTree                                          )

# Generate a weighted subsample of the original data with same number of subhalos as emulator subhalo population after filter is applied
weight_   = weight[subhalos]
selection = np.arange(0, weight_.size, 1, dtype=int)
subsample = np.random.choice(selection, size=countSample, replace=True, p=weight_/np.sum(weight_))

# Plot the loss function vs. training epoch.
plt.figure()
plt.plot(history.history["loss"    ])
plt.plot(history.history["val_loss"])
plt.legend(["train", "validation"], loc="upper right")
plt.ylabel("loss function", fontsize = 'large')
plt.xlabel("epoch"        , fontsize = 'large')
plt.savefig(args.outputDirectory+'/lossFunction.pdf')
plt.clf()

# Plot CDF of concentrations.
concentrations = np.linspace(3, 15, 100)
cdfGalacticus  = []
cdfEmulator    = []
for concentration in concentrations:
    cdfGalacticus.append(np.sum(data                   [subsample, 1] < concentration))
    cdfEmulator  .append(np.sum(concentrationNormalized[filter      ] < concentration))
cdfGalacticus                  = np.array(cdfGalacticus)/len(data                   [subsample, 1])
cdfEmulator                    = np.array(cdfEmulator  )/len(concentrationNormalized[filter      ])
concentrationDifferenceMaximum = concentrations[np.argmax(np.abs(cdfGalacticus-cdfEmulator))]
plt.figure()
plt.plot(concentrations, cdfGalacticus, 'k-', label = 'Galacticus')
plt.plot(concentrations, cdfEmulator  , 'r-', label = 'Emulator'  )
plt.axvline(x = concentrationDifferenceMaximum, linestyle = '--', color = 'grey')
plt.ylim(0, 1)
plt.xlabel('Concentration', fontsize = 'large')
plt.ylabel('CDF'          , fontsize = 'large')
plt.legend()
plt.savefig(args.outputDirectory+'/concentrationCDF.pdf')
plt.clf()

# Compute two-sample K-S tests for each property.
print('Two-sample K-S tests:')
print('        infall mass: ', stats.ks_2samp(data[subsample, 0], xt[filter, 1]))
print('      concentration: ', stats.ks_2samp(data[subsample, 1], xt[filter, 1]))
print('         bound mass: ', stats.ks_2samp(data[subsample, 2], xt[filter, 2]))
print('    infall redshift: ', stats.ks_2samp(data[subsample, 3], xt[filter, 3]))
print('  truncation radius: ', stats.ks_2samp(data[subsample, 4], xt[filter, 4]))
print('   projected radius: ', stats.ks_2samp(data[subsample, 5], xt[filter, 5]))

# Make 2D density plots.
from scipy.stats import gaussian_kde
## Concentration.
concentration_density_galacticus     = np.vstack([data[:, 0][subsample], data[:,1][subsample]])
concentration_density_emulated       = np.vstack([xt[filter, 0], xt[filter, 1]])
z1_galacticus                        = gaussian_kde(concentration_density_galacticus)(concentration_density_galacticus)
z1_emulated                          = gaussian_kde(concentration_density_emulated)(concentration_density_emulated)
## Bound mass.
mass_bound_density_galacticus        = np.vstack([data[:, 0][subsample], data[:, 2][subsample]])
mass_bound_density_emulated          = np.vstack([xt[filter, 0], xt[filter, 2]])
z2_galacticus                        = gaussian_kde(mass_bound_density_galacticus)(mass_bound_density_galacticus)
z2_emulated                          = gaussian_kde(mass_bound_density_emulated)(mass_bound_density_emulated)
## Infall redshift.
redshift_infall_density_galacticus   = np.vstack([data[:, 0][subsample], data[:, 3][subsample]])
redshift_infall_density_emulated     = np.vstack([xt[filter, 0], xt[filter, 3]])
z3_galacticus                        = gaussian_kde(redshift_infall_density_galacticus)(redshift_infall_density_galacticus)
z3_emulated                          = gaussian_kde(redshift_infall_density_emulated)(redshift_infall_density_emulated)
## Projected radius.
orbital_radius_density_galacticus    = np.vstack([data[:, 0][subsample], data[:, 4][subsample]])
orbital_radius_density_emulated      = np.vstack([xt[filter, 0], xt[filter, 4]])
z4_galacticus                        = gaussian_kde(orbital_radius_density_galacticus)(orbital_radius_density_galacticus)
z4_emulated                          = gaussian_kde(orbital_radius_density_emulated)(orbital_radius_density_emulated)
## Truncation radius.
truncation_radius_density_galacticus = np.vstack([data[:, 0][subsample], data[:, 5][subsample]])
truncation_radius_density_emulated   = np.vstack([xt[filter, 0], xt[filter, 5]])
z5_galacticus                        = gaussian_kde(truncation_radius_density_galacticus)(truncation_radius_density_galacticus)
z5_emulated                          = gaussian_kde(truncation_radius_density_emulated)(truncation_radius_density_emulated)
## Create the plot.
f, axes = plt.subplots(5, 2)
f.set_size_inches(15, 18)
axes[0, 0].scatter(data[:, 0][subsample], data[:, 1][subsample], c = z1_galacticus, s=9)
axes[0, 0].set(title="Galacticus", xlabel="mass infall", ylabel="concentration")
axes[0, 1].scatter(xt[filter, 0], xt[filter, 1], c = z1_emulated, s=9)
axes[0, 1].set(title="Generated", xlabel="mass infall", ylabel="concentration")
axes[1, 0].scatter(data[:, 0][subsample], data[:, 2][subsample], c = z2_galacticus, s=9)
axes[1, 0].set(title="Galacticus", xlabel="mass infall", ylabel="mass bound")
axes[1, 1].scatter(xt[filter, 0], xt[filter, 2], c = z2_emulated, s=9)
axes[1, 1].set(title="Generated", xlabel="mass infall", ylabel="mass bound")
axes[2, 0].scatter(data[:, 0][subsample], data[:, 3][subsample], c = z3_galacticus, s=9)
axes[2, 0].set(title="Galacticus", xlabel="mass infall", ylabel="redshift infall")
axes[2, 1].scatter(xt[filter, 0], xt[filter, 3], c = z3_emulated, s=9)
axes[2, 1].set(title="Generated", xlabel="mass infall", ylabel="redshift infall")
axes[3, 0].scatter(data[:, 0][subsample], data[:, 4][subsample], c = z4_galacticus, s=9)
axes[3, 0].set(title="Galacticus", xlabel="mass infall", ylabel="orbital radius")
axes[3, 1].scatter(xt[filter, 0], xt[filter, 4], c = z4_emulated, s=9)
axes[3, 1].set(title="Generated", xlabel= "mass infall", ylabel="orbital radius")
axes[4, 0].scatter(data[:, 0][subsample], data[:, 5][subsample], c = z5_galacticus, s=9)
axes[4, 0].set(title="Galacticus", xlabel="mass infall", ylabel="truncation radius")
axes[4, 0].set_ylim([-4.0, 0])
axes[4, 1].scatter(xt[filter, 0], xt[filter, 5], c = z5_emulated, s=9)
axes[4, 1].set(title="Generated", xlabel="mass infall", ylabel="truncation radius")
axes[4, 1].set_ylim([-4.0, 0])
plt.savefig(args.outputDirectory+'/density2D.pdf')
plt.clf()

# Make 1D PDF plots.
f, axes = plt.subplots(6)
f.set_size_inches(15, 20)
axes[0].hist(data[:     , 0][subsample], bins = 70, range = (-5, 0), label = 'Galacticus', fill = True , edgecolor = 'blue'  )
axes[0].hist(xt  [filter, 0]           , bins = 70, range = (-5, 0), label = 'Emulated'  , fill = False, edgecolor = 'orange')
axes[0].set(title = 'Infall Mass'      )
axes[0].legend()
axes[1].hist(data[:     , 1][subsample], bins = 70, range = (0, 30), label = 'Galacticus', fill = True , edgecolor = 'blue'  )
axes[1].hist(xt  [filter, 1]           , bins = 70, range = (0, 30), label = 'Emulated'  , fill = False, edgecolor = 'orange')
axes[1].set(title = 'Concentration'    )
axes[1].legend()
axes[2].hist(data[:     , 2][subsample], bins = 70, range = (-4, 0), label = 'Galacticus', fill = True , edgecolor = 'blue'  )
axes[2].hist(xt  [filter, 2]           , bins = 70, range = (-4, 0), label = 'Emulated'  , fill = False, edgecolor = 'orange')
axes[2].set(title = 'Bound Mass'       )
axes[2].legend()
axes[3].hist(data[:     , 3][subsample], bins = 70, range = (0, 10), label = 'Galacticus', fill = True , edgecolor = 'blue'  )
axes[3].hist(xt  [filter, 3]           , bins = 70, range = (0, 10), label = 'Emulated'  , fill = False, edgecolor = 'orange')
axes[3].set(title = 'Infall Redshift'  )
axes[3].legend()
axes[4].hist(data[:     , 4][subsample], bins = 70, range = (-2, 2), label = 'Galacticus', fill = True , edgecolor = 'blue'  )
axes[4].hist(xt  [filter, 4]           , bins = 70, range = (-2, 2), label = 'Emulated'  , fill = False, edgecolor = 'orange')
axes[4].set(title = 'Orbital radius'   )
axes[4].legend()
axes[5].hist(data[:     , 5][subsample], bins = 70, range = (-4, 0), label = 'Galacticus', fill = True , edgecolor = 'blue'  )
axes[5].hist(xt  [filter, 5]           , bins = 70, range = (-4, 0), label = 'Emulated'  , fill = False, edgecolor = 'orange')
axes[5].set(title = 'Truncation radius')
axes[5].legend()
plt.savefig(args.outputDirectory+'/density1D.pdf')
plt.clf()

# Make a plot showing the negative binomial distribution for subhalo numbers.
widthIntrinsic     = 0.18
probabilitySuccess = 1.0/(1.0+countSubhalosMean*widthIntrinsic**2)
rateStopping       = 1.0/widthIntrinsic**2
count              = np.arange(stats.nbinom.ppf(0.01, rateStopping, probabilitySuccess),stats.nbinom.ppf(0.99, rateStopping, probabilitySuccess))
plt.plot(count, stats.nbinom.pmf(count, rateStopping, probabilitySuccess), 'ko', label='Negative binomial PMF')
plt.savefig(args.outputDirectory+'/countSubhaloPMF.pdf')
plt.clf()

print('SUCCESS: subhalo emulator training')
