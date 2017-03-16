package NBody::Gadget;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';

# Module providing functions for manipulating Gadget parameter files.

sub readParameters {
    # Read Gadget parameter file and return a hash of parameter values.
    my $gadgetParametersFileName = shift();
    # Extract Gadget parameters.
    my %gadgetParameters;
    die("NBody::Gadget::readParameters(): cannot find Gadget parameter file")
	unless ( -e $gadgetParametersFileName );
    open(my $gadgetParametersFile,$gadgetParametersFileName);
    while ( my $line = <$gadgetParametersFile> ) {
	if ( $line =~ m/^([a-zA-Z0-9\._]+)\s+([a-zA-Z0-9\._\+\-\/]+)/ ) {
	    $gadgetParameters{$1} = $2;
	}
    }
    close($gadgetParametersFile);
    # Return the hash.
    return %gadgetParameters;
}

sub writeParameters {
    # Write a Gadget parameter file given a hash of parameters.
    my %gadgetParameters         = %{shift()};
    my $gadgetParametersFileName =   shift() ;
    open(my $gadgetParametersFile,">".$gadgetParametersFileName);
    print $gadgetParametersFile $_." ".$gadgetParameters{$_}."\n"
	foreach ( sort(keys(%gadgetParameters)) );
    close($gadgetParametersFile);
}

1;
