# Contains a Perl module which implements matching lines in a file.

package File::Match;
use strict;
use warnings;

sub Get_Matching_Lines {
    # Return a structure containing all lines from a file that match a regex,
    my $fileName = shift;
    my $regEx    = shift;
    # Open the file, and read each line.
    my @matches;
    open(my $fileHandle,$fileName);
    while ( my $line = <$fileHandle> ) {
	if ( my @submatches = $line =~ $regEx ) {
	    push(
		@matches,
		{
		    line       => $line,
		    submatches => \@submatches
		}
		);
	}
    }
    close($fileName);
    return @matches;
}

1;
