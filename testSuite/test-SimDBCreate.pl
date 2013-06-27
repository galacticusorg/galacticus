#!/usr/bin/env perl
use strict;
use warnings;

# Run a script to create a SimDB document descrbining Galacticus.
# Andrew Benson (25-January-2013)

# Run the creation script.
system("cd ..; scripts/aux/createSimDBDocument.pl");

exit;
