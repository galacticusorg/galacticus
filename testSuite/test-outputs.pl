#!/usr/bin/env perl
use strict;
use warnings;

# Run a set of short Galacticus models exploring a range of output options from Galacticus to ensure
# that they at least run to completion.
# Andrew Benson (26-Jan-2011)

# Simply run the models.
system("cd ..; scripts/aux/Run_Galacticus.pl testSuite/test-outputs.xml");

exit;
