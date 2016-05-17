#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
my $haveColor = eval
{
    require Term::ANSIColor;
    Term::ANSIColor->import();
};

# Remap line numbers in preprocessed files to their original sources.
# Andrew Benson (11-April-2016)

# Get the name of the reference file.
die("Usage: postProcess.pl <preprocessedSource>")
    unless ( scalar(@ARGV) == 1 );
my $preprocessedSourceName = $ARGV[0];

# Determine if interactive.
$haveColor = -t STDOUT ? $haveColor : 0;

# Initalize a map.
my @map;

# Initialize a hash of (possibly) unused functions.
my %unusedFunctions;

# Open and read the file.
my $lineNumber = 0;
push(
    @map,
    {
	source       => $preprocessedSourceName,
	line         => 1,
	lineOriginal => 1
    }
    );
open(my $file,$preprocessedSourceName);
while ( my $line = <$file> ) {
    ++$lineNumber;
    if ( $line =~ m/^\!\-\-\>\s+(\d+)\s+\"([a-zA-Z0-9_\.\/\(\):]+)\"\s*$/ ) {
	push(
	    @map,
	    {
		source       => $2,
		line         => $1,
		lineOriginal => $lineNumber+1 # We add 1 here because the line marker actually refers to the next line in the file.
	    }
	    );	 
    }
    # Capture unused function attributes.
    if ( $line =~ m/^\s*!\$GLC\s+function\s+attributes\s+unused\s*::\s*([a-zA-Z0-9_,\s]+)\s*$/ ) {
	(my $functions = $1) =~ s/\s*$//;
	$unusedFunctions{lc($_)} = 1
	    foreach ( split(/\s*,\s*/,$functions) );
    }
}
close($file);

# Do the remapping.
my $buffer;
my $status = 0;
my $functionName;
while ( my $line = <STDIN> ) {
    if ( $line =~ m/^([a-zA-Z0-9_\.\/]+\.p\.F90):(\d+):(\d+):\s*$/ ) {
	my $fileName     = $1;
	my $lineOriginal = $2;
	my $flag         = $3;
	my $source;
	foreach ( @map ) {
	    last
		if ( $_->{'lineOriginal'} > $lineOriginal );
	    $source     = $_->{'source'};
	    my $lineNumber;
	    if ( $source =~ m/\(\)$/ ) {
		$lineNumber = "meta";
	    } else {
		$lineNumber = $lineOriginal-$_->{'lineOriginal'}+$_->{'line'};
	    }
	    $line = $source.":".$lineNumber.":".$flag."\n";
	}
	print $buffer
	    if ( $buffer );
	undef($buffer);
    }
    my $dropBuffer = 0;
    # <workaround type="gfortran" PR="58175" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58175"/>
    $dropBuffer = 1
	if ( $line =~ m/Only array FINAL procedures declared for derived type/ );
    # Handle unused function attributes.
    if ( $line =~ m/^\s*subroutine\s+([a-z0-9_]+)/i ) {
	$functionName = lc($1);
    }
    if ( $line =~ m/\[\-Wunused\-function\]/ && defined($functionName) ) {
	$dropBuffer = 1
	    if ( exists($unusedFunctions{lc($functionName)}) );
	undef($functionName);
    }
    # Determine when to print the buffered output.
    my $printBuffer = 0;
    $printBuffer = 1
	if ( $line =~ m/^(Error|Warning):/ );
    $status = 1
	if ( $line =~ m/^Error:/ );
    if ( $haveColor ) {
    	if ( $line =~ m/^Warning:\s/ ) {
    	    my $warning = colored(['bright_magenta bold'],"Warning: ");
    	    $line =~ s/^Warning:\s/$warning/;
	    my $bold  = color('bold' );
	    my $reset = color('reset');
	    $line =~ s/\'([^\']+)\'/\'$bold$1$reset\'/g
    	}
    	if ( $line =~ m/^\s*\^\s*$/ ) {
    	    my $arrow = colored(['bright_green bold'],"^");
    	    $line =~ s/\^/$arrow/;
    	}
    	if ( $line =~ m/^\s*1\s*$/ ) {
    	    my $number = colored(['bright_green bold'],"1");
    	    $line =~ s/1/$number/;
    	}
    }
    $buffer .= $line;
    if ( $dropBuffer ) {
	undef($buffer);
    } elsif ( $printBuffer ) {
	print $buffer;
	undef($buffer);
    }
}
print $buffer
    if ( $buffer );

exit $status;
