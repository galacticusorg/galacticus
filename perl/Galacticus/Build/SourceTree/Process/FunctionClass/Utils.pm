# Contains a Perl module which provides utility functions for processing of functionClass directives.

package Galacticus::Build::SourceTree::Process::FunctionClass::Utils;
use strict;
use warnings;
use utf8;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use List::Uniq ':all';
use Exporter 'import';

our @EXPORT_OK = qw(
    LaTeX_Breakable
    trimlc
    striplc
    lctrim
    stripVariableName
    declarationRank
);

sub Class_Dependencies {
    # Find which type a functionClass implementation extends.
    my $classNode     = shift();
    my $directiveName = shift();
    my $classDepth    = 0;
    my $class;
    my @dependencies;
    while ( $classNode ) {
	# Collect class directives.
	if ( $classNode->{'type'} eq $directiveName ) {
	    $class->{'node'} = $classNode;
	    $class->{$_    } = $classNode->{'directive'}->{$_}
	        foreach ( sort(keys(%{$classNode->{'directive'}})) );
	}
	if ( $classNode->{'type'} eq "type" ) {
	    # Parse class openers to find dependencies.
	    if (
		$classNode->{'opener'} =~ m/^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*$directiveName([a-z0-9_]+)\s*$/i
		&&
		defined($2)
		) {
		$class->{'extends'} = $2;
		$class->{'type'   } = $directiveName.$4;
		push(@dependencies,$class->{'extends'});
		# Also determine if any other members of this class are used in this type definition, and add suitable dependencies.
		my $childNode = $classNode->{'firstChild'};
		while ( $childNode ) {
		    if ( $childNode->{'type'} eq "declaration" ) {
			foreach my $declaration ( @{$childNode->{'declarations'}} ) {
			    push(@dependencies,$declaration->{'type'})
				if ( ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" ) && $declaration->{'type'} =~ m/^$directiveName/ && $declaration->{'type'} !~ m/Class\s*$/ );
			}
		    }
		    $childNode = $childNode->{'sibling'};
		}
		last;
	    }
	}
	$classNode = &Galacticus::Build::SourceTree::Walk_Tree($classNode,\$classDepth);
    }
    @dependencies = uniq(sort(@dependencies));
    return ($class, @dependencies);
}

sub LaTeX_Breakable {
    my $text = shift;
    $text =~ s/([a-z])([A-Z])/$1\\-$2/g;
    return $text;
}

sub trimlc {
    (my $result = lc(shift())) =~ s/^\s+|\s+$//g;
    return $result;
}

sub striplc {
    (my $result = lc(shift())) =~ s/\s//g;
    return $result;
}

sub lctrim {
    # Trim trailing whitespace and return lowercased.
    my $string = shift();
    $string =~ s/\s*$//;
    return lc($string);
}

sub stripVariableName {
    # Strip away anything (e.g. array indices, assignment operators) after the variable name.
    (my $name = shift()) =~ s/^([a-zA-Z0-9_]+).*/$1/;
    return $name;
}

sub declarationRank {
    # Return the rank (number of array dimensions) of a variable declaration.
    my $declaration = shift();
    return 0
	unless ( grep {$_ =~ m/^dimension\s*\(/} @{$declaration->{'attributes'}} );
    my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,:\s]+)\)/} @{$declaration->{'attributes'}});
    return ($dimensionDeclarator =~ tr/,//)+1;
}

1;
