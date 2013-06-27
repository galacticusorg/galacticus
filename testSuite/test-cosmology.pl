#!/usr/bin/env perl
use strict;
use warnings;

# Run a set of short Galacticus models to explore different cosmological models.
# Andrew Benson (05-Sep-2010)

# Simply run the models.
system("cd ..; scripts/aux/Run_Galacticus.pl testSuite/test-cosmology.xml");

exit;
