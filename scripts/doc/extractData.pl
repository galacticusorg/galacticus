#!/usr/bin/env perl
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;
use LaTeX::Encode;
use Fcntl qw(SEEK_SET);
use UNIVERSAL;
use Fortran::Utils;
use Storable qw(dclone);
use List::ExtraUtils;
use Scalar::Util qw(reftype);

# Scan Fortran90 source code and extract various useful data from "!@" lines.
# Andrew Benson 12-Mar-2010

# Get source directory
die 'Usage: extractData.pl <sourceDir> <outputRoot>'
    unless ( scalar(@ARGV) == 2 );
my $sourceDir  = $ARGV[0];
my $outputRoot = $ARGV[1];
# Create XML object.
my $xml = new XML::Simple;
# Parse all class definition files.
my $classes;
opendir(my $directoryClasses,"./work/build/");
while ( my $fileName = readdir($directoryClasses) ) {
    next
	unless ( $fileName =~ m/\.classes\.xml$/ );
    my $classesFile = $xml->XMLin("./work/build/".$fileName);
    foreach my $className ( keys(%{$classesFile}) ) {
	$classes->{$className} = $classesFile->{$className};
    }
}
# Process child/parent relationships, and check for missing method descriptions.
foreach my $className ( sort(keys(%{$classes})) ) {
    $classes->{$className}->{'isFunctionClass'} = 0;
    my $parentName = $classes->{$className}->{'extends'};
    if ( $parentName ) {
	# Add this class to the parent's list of children.
	push(@{$classes->{$parentName}->{'children'}},$className);
	# Detect functionClass classes.
	my $ancestorName = $parentName;
	while ( $ancestorName ) {
	    $classes->{$className}->{'isFunctionClass'} = 1
		if ( $ancestorName eq "functionClass" );
	    $ancestorName = exists($classes->{$ancestorName}->{'extends'}) ? $classes->{$ancestorName}->{'extends'} : undef();
	}
	# Remove any "missing" methods which are defined in an ancestor class.
	unless ( reftype($classes->{$className}->{'missingMethods'}) ) {
	    my @missingMethods = split(" ",$classes->{$className}->{'missingMethods'});
	    if ( scalar(@missingMethods) > 0 ) {
		my @remainingMissingMethods;
		foreach my $missingMethod ( @missingMethods ) {
		    my $methodFound  = 0;
		    my $ancestorName = $parentName;
		    while ( $ancestorName ) {
			if ( grep {lc($_->{'method'}) eq lc($missingMethod)} &List::ExtraUtils::as_array($classes->{$ancestorName}->{'descriptions'}) ) {
			    $methodFound = 1;
			    last;
			} elsif ( defined($classes->{$ancestorName}->{'genericUses'}) && grep {lc($_) eq lc($missingMethod)} split(" ",$classes->{$ancestorName}->{'genericUses'}) ) {
			    $methodFound = 1;
			    last;
			}
			$ancestorName = exists($classes->{$ancestorName}->{'extends'}) ? $classes->{$ancestorName}->{'extends'} : undef();
		    }
		    push(@remainingMissingMethods,$missingMethod)
			unless ( $methodFound );
		}
		$classes->{$className}->{'missingMethods'} = join(" ",@remainingMissingMethods);
	    }
	}
    }
    # Check for any unresolved missing methods.
    unless ( reftype($classes->{$className}->{'missingMethods'}) ) {
	my @missingMethods = split(" ",$classes->{$className}->{'missingMethods'});
	if ( scalar(@missingMethods) > 0 ) {
	    print "Warning: missing method descriptions in class '".$className."':\n";
	    foreach my $missingMethod ( @missingMethods ) {
		print "\t".$missingMethod."\n";
	    }
	}
    }
}

# Open output files.
open(my $methodsFile    ,">".$outputRoot."Methods.tex");
# Write method descriptions.
my @functionClassExcludes = ( "allowedParameters", "autoHook", "deepCopy", "deepCopyReset", "descriptor", "hashedDescriptor", "objectType", "stateRestore", "stateStore" );
foreach my $className ( sort(keys(%{$classes})) ) {    
    print $methodsFile "\\subsection{\\large {\\normalfont \\ttfamily ".latex_encode($classes->{$className}->{'name'})."}}\\label{class:".$classes->{$className}->{'name'}."}\\hyperdef{class}{".$classes->{$className}->{'name'}."}{}\n\n";
    print $methodsFile "\\emph{Physics model:} \\refPhysics{".$classes->{$className}->{'name'}."}\n\n"
	if ( $classes->{$className}->{'isFunctionClass'} );
    print $methodsFile "\\noindent\\emph{Parent class:} \\refClass{".$classes->{$className}->{'extends'}."}\n\n"
	if ( exists($classes->{$className}->{'extends'}) );
    if ( exists($classes->{$className}->{'children'}) ) {
	my @sortedChildren = sort(@{$classes->{$className}->{'children'}});
	print $methodsFile "\\noindent\\emph{Child classes:}\n\n\\begin{tabular}{ll}\n";
	for(my $i=0;$i<scalar(@sortedChildren);$i+=2) {
	    print $methodsFile    "\\refClass{".$sortedChildren[$i  ]."}";
	    print $methodsFile " & \\refClass{".$sortedChildren[$i+1]."}"
		if ( $i+1 < scalar(@sortedChildren) );
	    print $methodsFile "\\\\\n";
	}
	print $methodsFile "\\end{tabular}\n\n";
    }
    if ( exists($classes->{$className}->{'descriptions'}) || $classes->{$className}->{'isFunctionClass'} ) {
	print $methodsFile "\\begin{description}\n";
	print $methodsFile "\\item[] All standard \\refClass{functionClass} methods (see \\S\\ref{sec:functionClassAll})\n"
	    if ( $classes->{$className}->{'isFunctionClass'} );

	if ( exists($classes->{$className}->{'descriptions'}) ) {
	    my @methods;
	    foreach my $method ( &List::ExtraUtils::as_array($classes->{$className}->{'descriptions'}) ) {
		push(@methods,$method)
		    unless ( $classes->{$className}->{'isFunctionClass'} && grep {$_ eq $method->{'method'}} @functionClassExcludes );
	    }
	    foreach my $method ( @methods ) {
		print $methodsFile "\\item[]{\\normalfont \\ttfamily ";
		if ( exists($method->{'type'}) ) {
		    (my $methodLabel = $method->{'method'}) =~ s/([^\\])_/$1\\_/g;
		    my $description = $method->{'description'};
		    chomp($description);
		    print $methodsFile $methodLabel."} ".$description."\n";
		    print $methodsFile "\\begin{itemize}\n";
		    print $methodsFile "\\item Return type: ".&declarationBuilder($method->{'type'}, variables => 0)."\n";
		    my @argumentList = &List::ExtraUtils::as_array($method->{'argumentList'});
		    my @arguments    = &List::ExtraUtils::as_array($method->{'arguments'   });
		    for(my $i=0;$i<scalar(@argumentList);++$i) {
			print $methodsFile "\\item Interface: ";
			if ( reftype($argumentList[$i]) ) {
			    print $methodsFile "{\\normalfont \\ttfamily ()}\\\\\n";
			} else {
			    my $argumentListEscaped = latex_encode($argumentList[$i]);
			    print $methodsFile "{\\normalfont \\ttfamily (".$argumentListEscaped.")}\\\\\n";
			    foreach my $argument ( &List::ExtraUtils::as_array($arguments[$i]->{'argument'}) ) {
				print $methodsFile &declarationBuilder($argument)."\\\\\n";
			    }
			}
		    }
		    print $methodsFile "\\end{itemize}\n";
		} else {
		    print "Warning: missing function type for method '".$method->{'method'}."' of class '".$className."'\n"; 
		}
	    }
	}
	print $methodsFile "\\end{description}\n";
    }
}
close($methodsFile);

exit;

sub declarationBuilder {
    my $type = shift();
    my %options = @_;
    my $variables = delete $options{'variables'} // 1;
    my $declaration;
    $declaration .= "{\\normalfont \\ttfamily ";
    if ( reftype($type) ) {
	$declaration .=
	    $type->{'intrinsic'}.
	    (exists($type->{'type'      }) && ! reftype($type->{'type'      })               ? "("   .               latex_encode                                 ($type->{'type'      } ).")" : "").
	    (exists($type->{'attributes'}) &&   defined($type->{'attributes'})               ? ", "  .join(", ",map {latex_encode($_)} &List::ExtraUtils::as_array($type->{'attributes'}))     : "").
	    (exists($type->{'variables' }) &&   defined($type->{'variables' }) && $variables ? " :: ".join(", ",map {latex_encode($_)} &List::ExtraUtils::as_array($type->{'variables' }))     : "");
    } elsif ( $type eq "subroutine" ) {
	$declaration .= "void";
    } else {
	die("unknown type");
    }
    $declaration .= "}";
    return $declaration;
}
