# Provides a module which returns the path to the Galacticus install location.

package Galacticus::Path;
use strict;
use warnings;
use Exporter qw(import);
use Cwd;
our $VERSION = 1.00;
our @EXPORT  = qw(galacticusPath galacticusPathExport);

sub galacticusPath {
    # Return the path to the Galacticus install location.
    my $galacticusPath;
    # First check the environment variable.
    if ( exists($ENV{'GALACTICUS_ROOT_V094'}) ) {
	# A path is specified via the environment variable, so simply use it.
	$galacticusPath  = $ENV{'GALACTICUS_ROOT_V094'};
	$galacticusPath .= "/"
	    unless ( $galacticusPath =~ m/\/$/ );
    } else {
	# No environment variable exists to specify the path, so we're required to attempt to determine it. First guess is that
	# it's the current directory.
	$galacticusPath = cwd();
	while ( $galacticusPath ne "" ) {
	    # Check for existance of main source file.
	    if ( -e $galacticusPath."/source/Galacticus.F90" ) {
		# Found, so this is the valid path.
		last;
	    } else {
		# Not found, try one level up in the directory hierarchy.
		$galacticusPath =~ s/\/[^\/]+$//;
	    }
	}
	if ( $galacticusPath eq "" ) {
	    # No install path found, so undefine our result.
	    undef($galacticusPath);
	} else {
	    # Add a terminating separator for convenience.
	    $galacticusPath .= "/";
	}
    }
    # Check we found the path and return it.
    die("Galacticus::Path::galacticusPath(): unable to determine Galacticus path")
	unless ( defined($galacticusPath) );
    return $galacticusPath;
}

sub galacticusPathExport {
    # Export the Galacticus path to the environment.
    $ENV{'GALACTICUS_ROOT_V094'} = &galacticusPath()
	unless ( exists($ENV{'GALACTICUS_ROOT_V094'}) );
    return 0;
}

1;
