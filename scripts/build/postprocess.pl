#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Fortran::Utils;
use File::Slurp qw(slurp);
use utf8;
use open ":std", ":encoding(UTF-8)";
my $haveColor = eval
{
use Term::ANSIColor;
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

# Initialize a hash of (possibly) unused functions.
my %unusedFunctions;

# Initialize a hash of variables which are asserted to be initialized.
my %initializedVariables;

# Initialize a hash of variables for which "pointer may outlive target" is ignored.
my %ignoreOutlives;

# Initialize a hash of variables for which "Unused PRIVATE module variable" is ignored.
my %ignoreUnused;

# Initialize a structure for unused variables.
my $unusedVariables;

# Initialize a structure for interoperable variables.
my $interoperableVariables;

# Parse the map of line numbers.
my @map =
    (
     {
	 source       => $preprocessedSourceName,
	 line         => 1,
	 lineOriginal => 1
     }
    );
open(my $mapFile,$preprocessedSourceName.".lmap");
while ( my $line = <$mapFile> ) {
    if ( $line =~ m/^\!\-\-\>\s+(\d+)\s+(\d+)\s+\"([a-zA-Z0-9_\-\.\/\(\):]+)\"\s*$/ ) {
	push(
	    @map,
	    {
		source       => $3,
		line         => $1,
		lineOriginal => $2
	    }
	    );	 
    }
}
close($mapFile);

# Open and read the file.
my $lineNumber = 0;
my $unitName;
open(my $file,$preprocessedSourceName);
while ( my $line = <$file> ) {
    ++$lineNumber;
    # Detect functions/subroutines/submodule procedure.
    foreach my $type ( 'subroutine', 'function', 'moduleProcedure' ) {
	if ( my @matches = ( $line =~ $Fortran::Utils::unitOpeners{$type}->{'regEx'} ) ) {		
	    $unitName = lc($matches[$Fortran::Utils::unitOpeners{$type}->{'unitName'}]);
	}
    }
    # Capture unused function attributes.
    if ( $line =~ m/^\s*!\$GLC\s+function\s+attributes\s+unused\s*::\s*([a-zA-Z0-9_,\s]+)\s*$/ ) {
	(my $functions = $1) =~ s/\s*$//;
	$unusedFunctions{lc($_)} = 1
	    foreach ( split(/\s*,\s*/,$functions) );
    }
    # Capture bogus uninitialized variable attributes.
    if ( $line =~ m/^\s*!\$GLC\s+attributes\s+initialized\s*::\s*([a-zA-Z0-9_,\s]+)\s*$/ ) {
	(my $variables = $1) =~ s/\s*$//;
	$initializedVariables{lc($_)} = 1
	    foreach ( split(/\s*,\s*/,$variables) );
    }
    # Capture pointer outlive ignores.
    if ( $line =~ m/^\s*!\$GLC\s+ignore\s+outlive\s*::\s*([a-zA-Z0-9_,\s]+)\s*$/ ) {
	(my $variables = $1) =~ s/\s*$//;
	$ignoreOutlives{lc($_)} = 1
	    foreach ( split(/\s*,\s*/,$variables) );
    }
    # Capture unused module ignores.
    if ( $line =~ m/^\s*!\$GLC\s+ignore\s+unused\s*::\s*([a-zA-Z0-9_,\s]+)\s*$/ ) {
	(my $variables = $1) =~ s/\s*$//;
	$ignoreUnused{lc($_)} = 1
	    foreach ( split(/\s*,\s*/,$variables) );
    }
    # Capture unused variable attributes.
    if ( $line =~ m/^\s*!\$GLC\s+attributes\s+unused\s*::\s*([a-zA-Z0-9_,\s]+)\s*$/ ) {
	(my $variables = $1) =~ s/\s*$//;
	$unusedVariables->{$unitName}->{lc($_)} = 1
	    foreach ( split(/\s*,\s*/,$variables) );
    }
    # Capture C-interoperable attributes.
    if ( $line =~ m/^\s*!\$GLC\s+attributes\s+interoperable\s*::\s*([a-zA-Z0-9_,\s]+)\s*$/ ) {
	(my $variables = $1) =~ s/\s*$//;
	$interoperableVariables->{lc($unitName)}->{lc($_)} = 1
	    foreach ( split(/\s*,\s*/,$variables) );
    }
}
close($file);

# Do the remapping.
my $buffer;
my $status = 0;
my $functionName;
my $procedureName;
my $pointerName;
my %bogusUninitialized;
my $lastDropped = 0;
while ( my $line = <STDIN> ) {
    if ( $line =~ m/^([a-zA-Z0-9_\.\/]+\.p\.F90):(\d+):([\d\-]+):\s*$/ ) {
	my $fileName     = $1;
	my $lineOriginal = $2;
	my $flag         = $3;
	my $source;
	foreach ( @map ) {
	    last
		if ( $_->{'lineOriginal'} > $lineOriginal );
	    $source     = $_->{'source'};
	    my $lineDescriptor;
	    if ( $source =~ m/\(\)$/ ) {
		$lineDescriptor = "auto-generated code (no line number)";
	    } else {
		my $lineNumber = $lineOriginal-$_->{'lineOriginal'}+$_->{'line'};
		$lineDescriptor = "line ".$lineNumber;
	    }
    	    $line = $source."; ".$lineDescriptor." [preprocessed line ".$lineOriginal."]; code ".$flag."\n";
	}
	print $buffer
	    if ( $buffer );
	undef($buffer);
    }
    my $dropBuffer = 0;
    # <workaround type="gfortran" PR="58175" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58175"/>
    $dropBuffer = 1
	if ( $line =~ m/Only array FINAL procedures declared for derived type/ );
    # <workaround type="gfortran" PR="86117" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86117"/>
    if ( $line =~ m/Warning:\s+['‘]MEM\[\(struct\s+([a-zA-Z0-9_]+)[a-z0-9_\s\*]+\)[a-z0-9_&]+\s+\+\s+\d+B\]['’] (is|may be) used uninitialized( in this function)?/ ) {
	$dropBuffer = 1;
	# In these cases we must add the problem structure to a list that we will ignore other "uninitialized" complaints about.
	$bogusUninitialized{$1} = 1;
    }
    # <workaround type="gfortran" PR="86117" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86117"/>
    if ( $line =~ m/Warning:\s+['‘](MEM\s<real[a-zA-Z0-9_\.<>\(\)=\[\]:\s\*]+>\s[a-zA-Z0-9_\.<>\(\)=\[\]:\s\*]+)['’] (is|may be) used uninitialized( in this function)?/ ) {
	$dropBuffer = 1;
	# In these cases we must add the problem structure to a list that we will ignore other "uninitialized" complaints about.
	$bogusUninitialized{$1} = 1;
    }
    if ( $line =~ m/Warning:\s+['‘]([a-zA-Z0-9_\.<>\(\)=\[\]:\s\*]+)['’]\s+may be used uninitialized( in this function)?/ ) {
	my @elements = split(/\./,$1);
	foreach my $symbolName ( keys(%bogusUninitialized) ) {
	    $dropBuffer = 1
		if ( grep {$_ eq $symbolName} @elements );
	}
    }
    # Handle unused function attributes.
    if ( $line =~ m/^\s*\d+\s*\|\s*subroutine\s+([a-z0-9_]+)/i ) {
	$functionName = lc($1);
    }
    if ( $line =~ m/\[\-Wunused\-function\]/ && defined($functionName) ) {
	$dropBuffer = 1
	    if ( exists($unusedFunctions{lc($functionName)}) );
	undef($functionName);
    }
    # Handle unused variable attributes.
    if ( $line =~ m/^\s*\d+\s*\|\s*(.+)/ ) {
	my $opener = $1;
	foreach my $type ( 'subroutine', 'function' ) {
	    if ( my @matches = ( $opener =~ $Fortran::Utils::unitOpeners{$type}->{'regEx'} ) ) {		
		$procedureName = lc($matches[$Fortran::Utils::unitOpeners{$type}->{'unitName'}]);
	    }
	}
    }
    if ( $line =~ m/^<stdin>:\d+:\d+:/ ) {
	$procedureName = "stdin";
    }
    if ( $line =~ m/^Warning: Unused dummy argument ['‘]([a-zA-Z0-9_]+)['’] at \(\d+\) \[\-Wunused\-dummy-argument\]/ && defined($procedureName) ) {
	my $variableName = $1;
	$dropBuffer = 1
	    if ( $procedureName eq "stdin" || exists($unusedVariables->{$procedureName}->{$variableName}) );
    }
    if ( $line =~ m/^Warning: Dummy argument ['‘]([a-zA-Z0-9_]+)['’] at \(\d+\) was declared INTENT\(OUT\) but was not set \[\-Wunused\-dummy-argument\]/ && defined($procedureName) ) {
	my $variableName = $1;
	$dropBuffer = 1
	    if ( $procedureName eq "stdin" || exists($unusedVariables->{$procedureName}->{$variableName}) );
    }
    if ( $line =~ m/^Warning: Derived-type dummy argument ['‘]([a-zA-Z0-9_]+)['’] at \(\d+\) was declared INTENT\(OUT\) but was not set and does not have a default initializer \[\-Wunused\-dummy\-argument\]/ && defined($procedureName) ) {
	my $variableName = $1;
	$dropBuffer = 1
	    if ( $procedureName eq "stdin" || exists($unusedVariables->{$procedureName}->{$variableName}) );
    }
    # Handle explicit C-interoperable dummy arguments.
    if ( $line =~ m/^Warning: Variable ['‘]([a-zA-Z0-9_]+)['’] at \(\d+\) is a dummy argument of the BIND\(C\) procedure ['‘]([a-z_]+)['’] but may not be C interoperable \[\-Wc\-binding\-type\]/ && defined($procedureName) ) {
	my $variableName = $1;
	my $inProcedure  = $2;
	$dropBuffer = 1
	    if (  exists($interoperableVariables->{lc($inProcedure)}->{lc($variableName)}) );
    }
    # Handle uninitialized variable attributes.
    if ( $line =~ /Warning: ['‘](\(\*)??([a-zA-Z0-9_]+)[a-zA-Z0-9_\.\[\]\s\{\}\):]*['’] may be used uninitialized( in this function)?? \[\-Wmaybe\-uninitialized\]/ ) {
	$dropBuffer = 1
	    if ( exists($initializedVariables{lc($2)}) );
    }
    if ( $line =~ /Warning: ['‘]([a-zA-Z0-9_]+)[a-zA-Z0-9_\.\[\]]*['’] is used uninitialized( in this function)?? \[\-Wuninitialized\]/ ) {
	$dropBuffer = 1
	    if ( exists($initializedVariables{lc($1)}) );
    }
    if ( $line =~ /note: ['‘]([a-zA-Z0-9_]+)[a-zA-Z0-9_\.\[\]]*['’]( was)?? declared here/ ) {
	$dropBuffer = 1
	    if ( exists($initializedVariables{lc($1)}) );
    }
    # Handle ignore "Unused PRIVATE module variable" warnings.
    if ( $line =~ m/^\s*Warning: Unused PRIVATE module variable ['‘]([a-zA-Z0-9_]+)['’] declared at \(1\) \[\-Wunused\-value\]/ ) {
	my $variableName = lc($1);
	$dropBuffer = 1
	    if ( exists($ignoreUnused{$variableName}) );
    }
    # Handle ignore "pointer may outlive target" warnings.
    if ( $line =~ m/^\s*\d+\s*\|\s*([a-z0-9_]+)\s*=>\s*[a-z0-9_]+/i ) {
	$pointerName = lc($1);
    }
    if ( $line =~ m/\[\-Wtarget\-lifetime\]/ ) {
	$dropBuffer = 1
	    if ( exists($ignoreOutlives{lc($pointerName)}) );
	undef($pointerName);
    }
    # Handle "note:"s.
    if ( $line =~m/^note:/ ) {
	$dropBuffer = 1
	    if ( $lastDropped );
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
	    $line =~ s/(['‘])([^['’]]+)(['’])/$1$bold$2$reset$3\'/g
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
	$lastDropped = 1;
    } elsif ( $printBuffer ) {
	print $buffer;
	undef($buffer);
	$lastDropped = 0;
    }
}
print $buffer
    if ( $buffer );

exit $status;
