#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
 $ENV{"GALACTICUS_ROOT_V093"} = "/";
}
unshift(@INC,$galacticusPath."perl"); 

# Test Galacticus version environment variable usage.
# Andrew Benson (1-November-2013)

# Search files for incorrect version environment.
system("cd ..; find . -name \"*.pm\" -or -name \"*.pl\" -or -name \"*.F90\" -or -name \"*.Inc\" | grep -v test-version.pl | xargs grep GALACTICUS_ROOT | grep -v GALACTICUS_ROOT_V093");
print "FAILED: incorrect version of GALACTICUS_ROOT environment variable used\n"
    if ( $? == 0 );

exit;
