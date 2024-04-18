# Contains a Perl module which provides utility functions for processing of functionClass directives.

package Galacticus::Build::SourceTree::Process::FunctionClass::Utils;
use strict;
use warnings;
use utf8;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use List::Uniq ':all';

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

1;
