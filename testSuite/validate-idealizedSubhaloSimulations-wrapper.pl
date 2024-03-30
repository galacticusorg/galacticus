#!/usr/bin/env perl
use strict;
use warnings;

# A simple wrapper script for validate-idealizedSubhaloSimulations.pl
# that runs each model in turn. This is needed to work around an issue
# in the PDL::IO::HDF5 library.
# Andrew Benson (22-March-2024)

for(my $i=0;$i<8;++$i) {
    system("./validate-idealizedSubhaloSimulations.pl ".$i);
}

exit;
