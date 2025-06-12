# Use a subhalo emulator to perform inference on a quad-lensing system.
import argparse
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import regularizers
import numpy as np
import tensorflow_probability as tfp
import h5py
from samana.forward_model import forward_model
import os
import sys
import random
from scipy import stats
from emulator_utils import norm_transform_inv, Coupling, RealNVP

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='emulatorSubhaloQuadLensingInference.py',description='Use an emulator to perform inference on a quad-lensing system.')
parser.add_argument('--outputDirectory', action='store',default='.',help='the directory to which emulator weights (and plots) were output and to which results should be output')
args = parser.parse_args()

# Reading in the data and converting it from a string to a list of lists and floats.
with h5py.File(args.outputDirectory+'/normalizationData.hdf5', 'r') as normalizationFile:
    minArray          = normalizationFile['minArray'][:]
    maxArray          = normalizationFile['maxArray'][:]
    massTree          = normalizationFile.attrs['massTree'         ]
    redshiftTree      = normalizationFile.attrs['redshiftTree'     ]
    radiusVirialHost  = normalizationFile.attrs['radiusVirialHost' ]
    massResolution    = normalizationFile.attrs['massResolution'   ]
    countSubhalosMean = normalizationFile.attrs['countSubhalosMean']

# Build an emulator and read in weights.
emulator = RealNVP(num_coupling_layers=12)
emulator.build(input_shape= (6))
emulator.load_weights(args.outputDirectory+'/emulator.weights.h5')

# Creating function that takes emulator weights (and other required data) to produce a population of un-normalized emulated subhalos.
def emulator_data(emulator = emulator, minArray = minArray, maxArray = maxArray, massTree = massTree, massResolution = massResolution, radiusVirialHost = radiusVirialHost, countSubhalosMean = countSubhalosMean, count_iterations = 1):

    # Define properties of the negative binomial distribution used to sample the number of subhalos.
    widthIntrinsic     = 0.18
    probabilitySuccess = 1.0/(1.0+countSubhalosMean*widthIntrinsic**2)
    rateStopping       = 1.0/widthIntrinsic**2
    countSubhalos      = np.arange(stats.nbinom.ppf(0, rateStopping, probabilitySuccess),stats.nbinom.ppf(0.9999999999999999, rateStopping, probabilitySuccess))
    countSubhalosPMF   = np.nan_to_num(stats.nbinom.pmf(countSubhalos,rateStopping,probabilitySuccess))
    # Sample a number of subhalos.
    countSubhalo       = int(np.random.choice(countSubhalos, p = countSubhalosPMF))
    # Build a sample of subhalos.
    iteration          =  0
    iterationsMaximum  = 10
    data               = np.empty((0,6))
    while len(data) < countSubhalo and iteration < iterationsMaximum:
        samples    = emulator.distribution.sample(countSubhalo)
        x, _       = emulator.predict(samples, batch_size=65336)
        xt         = norm_transform_inv(x, minArray, maxArray, -1, 1)
        filter     = (xt[:,0] > np.log10(2.0*massResolution/massTree)) & (xt[:,2] <= 0.0) & (xt[:,2] > -xt[:,0]+np.log10(massResolution/massTree)) & (xt[:,3] >= redshiftTree) & (xt[:,2] < np.log10(1e9/(massTree*10**xt[:,0])))
        data       = np.vstack((data,xt[filter]))
    if len(data) < countSubhalo:
        sys.exit('failed to generate sufficient subhalos')
    # Truncate the number of subhalos to the required number.
    data = data[:countSubhalo]
    # Extract unnormalized subhalo properties.
    massInfall       = massTree        *(10.0**data[:,0])
    concentration    =                         data[:,1]
    massBound        = massInfall      *(10.0**data[:,2])
    redshift         =                         data[:,3]
    radiusTruncation = radiusVirialHost*(10.0**data[:,4])
    radiusProjected  = radiusVirialHost*(10.0**data[:,5])

    # Find (x,y) coordinates.
    phi   = 2.0*np.pi*np.random.uniform(size=len(radiusProjected))
    x     = radiusProjected*np.cos(phi)
    y     = radiusProjected*np.sin(phi)

    return massInfall, concentration, massBound, redshift, radiusTruncation, x, y

# Construct input parameters for the forward_model() function of samana.
job_index                   = 1
n_keep                      = 2
summary_statistic_tolerance = 1e5

from samana.Data.Mocks.baseline_smooth_mock import BaselineSmoothMockModel
from samana.Data.Mocks.baseline_smooth_mock import BaselineSmoothMock
data_class        = BaselineSmoothMock()
model             = BaselineSmoothMockModel
preset_model_name = 'DMEmulator'

# Sample options.
kwargs_sample_realization = {}
kwargs_sample_realization['LOS_normalization'        ] = ['FIXED',  0.0         ]
kwargs_sample_realization['log_m_host'               ] = ['FIXED', 13.3         ]
kwargs_sample_realization['cone_opening_angle_arcsec'] = ['FIXED',  8.0         ]
kwargs_sample_realization['log_mlow'                 ] = ['FIXED',  6.0         ]
kwargs_sample_realization['log_mhigh'                ] = ['FIXED',  9.0         ]
kwargs_sample_realization['sigma_sub'                ] = ['FIXED',  0.0         ]
kwargs_sample_realization['log_mc'                   ] = ['FIXED',  4.0         ]
kwargs_sample_realization['emulator_data_function'   ] = ['FIXED', emulator_data]

# Source options.
kwargs_sample_source      = {'source_size_pc': ['FIXED', 5]}
kwargs_sample_macro_fixed = {
    'a4_a':         ['FIXED'   , 0.0              ], 
    'a3_a':         ['FIXED'   , 0.0              ],
    'delta_phi_m3': ['GAUSSIAN', -np.pi/6, np.pi/6]
}
kwargs_model_class = {'shapelets_order': 10}

# Run the forward model analysis.
forward_model(args.outputDirectory+"/", job_index, n_keep, data_class, model, preset_model_name, 
              kwargs_sample_realization, kwargs_sample_source, kwargs_sample_macro_fixed,
              tolerance=summary_statistic_tolerance, log_mlow_mass_sheets = 6.0, kwargs_model_class = kwargs_model_class, verbose=False, test_mode=False)

fluxes = np.genfromtxt(args.outputDirectory+'/job_'+str(job_index)+'/fluxes.txt')
print('Flux ratios:')
for i in range(len(fluxes)):
    f2f1 = fluxes[i][1]/fluxes[i][0]
    f3f1 = fluxes[i][2]/fluxes[i][0]
    f4f1 = fluxes[i][3]/fluxes[i][0]
    print(f'   {f2f1}, {f3f1}, {f4f1}')


print('SUCCESS')
