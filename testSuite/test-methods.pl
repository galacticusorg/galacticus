#!/usr/bin/env perl

# Run a set of short Galacticus models spanning a full range of method options to ensure
# that they at least run to completion.
# Andrew Benson (04-Sep-2010)

# Simply run the models.
system("cd ..; scripts/aux/Run_Galacticus.pl testSuite/test-methods.xml");

exit;
