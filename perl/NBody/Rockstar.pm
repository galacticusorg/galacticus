package NBody::Rockstar;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';

# Module providing functions for manipulating Rockstar parameter files.

sub readParameters {
    # Read Rockstar parameter file and return a hash of parameter values.
    my $rockstarParametersFileName = shift();
    # Extract Rockstar parameters.
    my %rockstarParameters;
    die("NBody::Rockstar::readParameters(): cannot find Rockstar parameter file")
	unless ( -e $rockstarParametersFileName );
    open(my $rockstarParametersFile,$rockstarParametersFileName);
    while ( my $line = <$rockstarParametersFile> ) {
	if ( $line =~ m/^([a-zA-Z0-9\._]+)\s*=\s*([a-zA-Z0-9\._\+\-\/]+)/ ) {	    
	    $rockstarParameters{$1} =  $2;
	    $rockstarParameters{$1} =~ s/^"//;
	    $rockstarParameters{$1} =~ s/"$//;
	}
    }
    close($rockstarParametersFile);
    # Return the hash.
    return %rockstarParameters;
}

1;
