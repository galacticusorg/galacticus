# Parse command line options.

package Options;
use strict;
use warnings;

sub Parse_Options {
    # Parses command line options of the form:
    #   --optionName value
    # and places them into a hash (which can have default values preloaded into it). Supports double
    # quotes in the "value" of options to allow for values that contain spaces,
    my @arguments = @{$_[0]};
    my $options   =   $_[1];
    my $iArg = 0;
    while ( $iArg < scalar(@arguments)-1 ) {
	++$iArg;
	if ( $arguments[$iArg] =~ m/^\-\-(.*)/ ) {
	    my $argument = $1;
	    my $value = $arguments[$iArg+1];
	    ++$iArg;
	    if ( $arguments[$iArg] =~ m/^\"/ ) {
		until ( $value =~ m/\"$/ ) {
		    ++$iArg;
		    $value .= $arguments[$iArg];
		}
		$value =~ s/^\"(.*)\"$/$1/;
	    }
	    $options->{$argument} = $value;
	}
    }
}

1;
