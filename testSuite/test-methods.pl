#!/usr/bin/env perl

# Run a set of short Galacticus models spanning a full range of method options to ensure
# that they at least run to completion.
# Andrew Benson (04-Sep-2010)

# Remove automatically generated files to force them to be regenerated (and, therefore, to test that the generating code works
# correctly).
#
# FSPS stellar population synthesis code and associated file.
system("rm -f data/stellarPopulations/SSP_Spectra_Conroy-et-al_v2.1_imfSalpeter.hdf5");
system("rm -rf aux/FSPS_v2.3");

# Simply run the models.
system("cd ..; scripts/aux/Run_Galacticus.pl testSuite/test-methods.xml");

exit;
