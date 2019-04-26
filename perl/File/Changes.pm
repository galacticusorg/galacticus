# Contains a Perl module which implements checking of whether a file has changed with respect to another version.

package File::Changes;
use strict;
use warnings;
use File::Copy;
use File::Compare;

sub Update {
    # Checks if a file1 is different from file2 and, if it is, replaces file1 with file2. Otherwise, file2 is deleted.

    # Get file names.
    my $oldFile = shift();
    my $newFile = shift();
    (my %options) = @_
	if ( scalar(@_) > 1 );

    # Check if the old file exists.
    if ( ! -e $oldFile ) {
	# It does not, so move the new file to the old file.
	move($newFile,$oldFile);	
    # Check if the files are the same.
    } elsif (compare($oldFile,$newFile) == 0) {
	# They are: remove the new file.
	unlink($newFile);
    } else {
	# They are not: replace the old file.
	move($newFile,$oldFile);
    }

    # If instructed, touch an update file to indicate that work was done.
    system("touch ".$oldFile.".up")
	if ( exists($options{'proveUpdate'}) && $options{'proveUpdate'} eq "yes" );
}

1;
