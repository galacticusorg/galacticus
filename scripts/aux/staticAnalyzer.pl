#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Galacticus::Build::SourceTree;
use Fortran::Utils;

# Perform static analysis of Fortran files.
# Andrew Benson (28-February-2023)

# Get the file to process.
die('Usage: staticAnalyzer.pl <fileName>')
    unless ( scalar(@ARGV) == 1 );
my $fileName = $ARGV[0];

# Parse the file.
my $tree = &Galacticus::Build::SourceTree::ParseFile($fileName);
# Walk the tree.
my $node         = $tree;
my $depth        = 0;
my $status       = 0;
my @types;
my @constructors;
my @destructors;
while ( $node ) {
    # Class/type pointers in derived types should be null initialized.
    if ( $node->{'type'} eq "declaration" ) {
	foreach my $declaration ( @{$node->{'declarations'}} ) {
	    if ( $node->{'parent'}->{'type'} eq "type" && ( $declaration->{'intrinsic'} eq "type" || $declaration->{'intrinsic'} eq "class" ) && grep {$_ eq "pointer"} @{$declaration->{'attributes'}} ) {
		for(my $i=0;$i<scalar(@{$declaration->{'variables'}});++$i) {
		    next
			if ( $declaration->{'variables'}->[$i] =~ m/=>null\(\)$/ );
		    (my $typeName = $node->{'parent'}->{'opener'}) =~ s/.*::\s*([a-zA-Z0-9_]+).*/$1/;
		    chomp($typeName);
		    print "Pointer variable '".$declaration->{'variableNames'}->[$i]."' in type '".$typeName."' in file '".$fileName."' is not null initialized\n";
		    $status = 1;
		}
	    }
	}
    }
    # Look for duplicated assignments in `constructorAssign` directives.
    if ( $node->{'type'} eq "constructorAssign" ) {
	my @variables = split(/\s*,\s*/,$node->{'directive'}->{'variables'});
	my %countAssignments;
	foreach my $variable ( @variables ) {
	    ++$countAssignments{$variable};
	}
	foreach my $variable ( keys(%countAssignments) ) {
	    if ( $countAssignments{$variable} > 1 ) {
		print "Duplicated assignment of '".$variable."' in `constructorAssign` directive in file '".$fileName."'\n";
		$status = 1;
	    }
	}
    }
    # Look for empty constructors.
    ## First identify any type definitions.
    if ( $node->{'type'} eq "type" ) {
	push(@types,lc($node->{'name'}));
    }
    ## Next identify constructor functions for a known type.
    if ( @types && $node->{'type'} eq "interface" && defined($node->{'name'}) && grep {$_ eq lc($node->{'name'})} @types ) {
	my $nodeChild = $node->{'firstChild'};
	while ( $nodeChild ) {
	    if ( $nodeChild->{'type'} eq "moduleProcedure" ) {
		push(@constructors,map {lc($_)} @{$nodeChild->{'names'}});
	    }
	    $nodeChild = $nodeChild->{'sibling'};
	}
    }
    ## Last, check any constructor functions.
    if ( $node->{'type'} eq "function" && grep {$_ eq lc($node->{'name'})} @constructors ) {
	# Check for arguments.
	if ( $node->{'opener'} =~ m/function\s+[a-zA-Z0-9_]+\(\)/ ) {
	    # No arguments - check for empty function.
	    if ( &functionIsEmpty($node) ) {
		print "Empty constructor function '".$node->{'name'}."' in file '".$fileName."'\n";
		$status = 1;
	    }	    
	}
    }
    # Look for empty finalizers.
    ## First identify finalizer functions.
    if ( $node->{'type'} eq "declaration" ) {
	foreach my $declaration ( @{$node->{'declarations'}} ) {
	    if ( $declaration->{'intrinsic'} eq "final" ) {
		push(@destructors,map {lc($_)} @{$declaration->{'variables'}});
	    }
	}
    }
    ## Then check finalizer functions.
    if ( $node->{'type'} eq "subroutine" && grep {$_ eq lc($node->{'name'})} @destructors ) {
	if ( &functionIsEmpty($node) ) {
	    print "Empty finalizer function '".$node->{'name'}."' in file '".$fileName."'\n";
	    $status = 1;
	}
    }
    # Walk to the next node in the tree.
    $node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
}

exit $status;

sub functionIsEmpty {
    # Return true if the provided function is empty (other than attribute ignores and a return statement.
    my $node    = shift();
    my $isEmpty = 1;
    my $depth   = 0;
    $node       = $node->{'firstChild'};
    while ( $node && $depth >= 0 ) {
	# Ignore module uses and declaration sections.
	if ( $node->{'type'} eq "moduleUse" || $node->{'type'} eq "declaration" ) {
	    $node = $node->{'sibling'};
	} else {
	    if ( $node->{'type'} eq 'code' ) {
		open(my $content,"<",\$node->{'content'});
		foreach my $line ( <$content> ) {
		    next
			if ( $line =~ m/^\s*$/ );
		    next
			if ( $line =~ m/^\s*![^\$]/ );
		    next
			if ( $line =~ m/^\s*use(\s*,\s*intrinsic)??\s*::/ );
		    next
			if ( $line =~ m/^\s*implicit\s+none\s*$/ );
		    next
			if ( $line =~ m/^\s*return\s*$/ );
		    next
			if ( $line =~ $Fortran::Utils::variableDeclarationRegEx );
		    # Non-ignorable content.
		    $isEmpty = 0;
		}
		close($content);
	    } else {
		# Some other structure encountered (e.g. a directive) - function is not empty.
		$isEmpty = 0;
		last;
	    }
	    # Walk to the next node in the tree.
	    $node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
	}
    }
    return $isEmpty;
}
