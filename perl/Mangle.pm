# Contains a Perl module which downloads and compiles the "mangle" code.

package Mangle;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 

sub Build {
    # Ensure that we have a compiled mangle.
    unless ( -e $galacticusPath."aux/mangle/bin/ransack" ) {
	# Ensure that we have the mangle source.
	unless ( -e $galacticusPath."aux/mangle" ) {	   
	    # Clone mangle
	    system("cd ".$galacticusPath."aux; git clone https://github.com/mollyswanson/mangle.git");
	    die("Mangle::Build: failed to clone mangle")
		unless ( $? == 0 && -e $galacticusPath."aux/mangle" );
	}
	# Build mangle.
	system("cd ".$galacticusPath."aux/mangle/src; ./configure; make cleanest; make");
	die("Mangle::Build: failed to build mangle")
	    unless ( $? == 0 && -e $galacticusPath."aux/mangle/bin/ransack" );
    }
}

1;

