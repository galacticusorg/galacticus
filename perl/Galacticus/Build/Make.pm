# Contains a Perl module which provides utilities for Galacticus' Make system.

package Make;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use Switch;

sub Module_Name {
    # Return the module name associated with a file.
    my $fileName = shift;
    (my %options) = @_
	if ( scalar(@_) > 1 );
    my $moduleName;
    if ( exists($options{'default'}) ) {
	switch ( $options{'default'} ) {
	    case ( 'self' ) {
		$moduleName = $fileName;
	    }
	}
    }
    my $moduleFileName = "./work/build/".$fileName.".m";
    if ( -e $moduleFileName ) {
	open(my $moduleHandle,$moduleFileName);
	($moduleName =<$moduleHandle>) =~ s/\.\/work\/build\/(.*)\.mod\n$/$1/
	    unless ( eof($moduleHandle) );
	close($moduleHandle);
    }
    return $moduleName;
}

1;
