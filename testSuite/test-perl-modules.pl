#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";

# Test Perl modules.
# Andrew Benson (7-July-2013)

# Check individual compilation.
system("cd ..; find . -type f -name \"*.pl\" | xargs -n 1 perl -c 1>perl.tmp 2>&1; sed -r s/failed/FAILED/g perl.tmp; rm perl.tmp");
system("cd ..; find perl -type f -name \"*.pm\" | xargs -n 1 perl -c 1>perl.tmp 2>&1; sed -r s/failed/FAILED/g perl.tmp; rm perl.tmp");

exit;
