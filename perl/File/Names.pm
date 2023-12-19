# Contains a Perl module which implements obtaining file names from file handles.

package File::Names;
use strict;
use warnings;

sub Get_Name {
    # Gets the file name associated with a given handle.
    my $fileHandle = shift();
    my $fileNumber = fileno($fileHandle);
    return readlink("/proc/".$$."/fd/".$fileNumber);
}

1;
