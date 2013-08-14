# Contains a Perl module which implements checking of whether a file has changed with respect to another version.

package File_Changes;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use File::Copy;
use File::Compare;

sub Update {
    # Checks if a file1 is different from file2 and, if it is, replaces file2 with file1. Otherwise, file1 is deleted.

    # Get file names.
    my $oldFile = $_[0];
    my $newFile = $_[1];

    # Check if the files are the same.
    if (compare($oldFile,$newFile) == 0) {
	# They are: remove the new file.
	unlink($newFile);
    } else {
	# They are not: replace the old file.
	move($newFile,$oldFile);
    }

}

1;
