#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/../perl';

# Test Galacticus version environment variable usage.
# Andrew Benson (1-November-2013)

# Search files for incorrect version environment.
system("cd ..; find . -name \"*.pm\" -or -name \"*.pl\" -or -name \"*.F90\" -or -name \"*.Inc\" | grep -v test-version.pl | xargs grep GALACTICUS_ROOT | grep -v GALACTICUS_ROOT_V094");
print "FAILED: incorrect version of GALACTICUS_ROOT environment variable used\n"
    if ( $? == 0 );

exit;
