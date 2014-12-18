# Contains a Perl module which handles downloading and compiling Cloudy.

package Cloudy;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 

sub Initialize {
    # Download and compile Cloudy so that it is ready for use.
    # Specify Cloudy version.
    my $cloudyVersion = "c13.03";
    # Specify Cloudy path.
    my $cloudyPath    = $galacticusPath."aux/".$cloudyVersion;
    # Download the code.
    unless ( -e $cloudyPath.".tar.gz" ) {
	print "Cloudy::Initialize: downloading Cloudy code.\n";
	system("wget \"http://data.nublado.org/cloudy_releases/c13/".$cloudyVersion.".tar.gz\" -O ".$cloudyPath.".tar.gz");
	die("Cloudy::Initialize: FATAL - failed to download Cloudy code.") 
	    unless ( -e $cloudyPath.".tar.gz" );
    }
    # Unpack the code.
    unless ( -e $cloudyPath ) {
	print "Cloudy::Initialize: unpacking Cloudy code.\n";
	system("tar -x -v -z -C aux -f ".$cloudyPath.".tar.gz");
	die("Cloudy::Initialize: FATAL - failed to unpack Cloudy code.")
	    unless ( -e $cloudyPath );
    }    
    # Build the code.
    unless ( -e $cloudyPath."/source/cloudy.exe" ) {
	print "Cloudy::Initialize: compiling Cloudy code.\n";
	system("cd ".$cloudyPath."/source; chmod u=wrx configure.sh capabilities.pl; make");
	die("Cloudy::Initialize: FATAL - failed to build Cloudy code.")
	    unless ( -e $cloudyPath."/source/cloudy.exe" );
    }
    # Return path and version.
    return $cloudyPath."/", $cloudyVersion;
}

1;
