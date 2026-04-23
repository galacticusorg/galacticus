# Contains a Perl module which provides descriptor parameter utilities for functionClass directives.

package Galacticus::Build::SourceTree::Process::FunctionClass::Descriptor;
use strict;
use warnings;
use utf8;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use List::ExtraUtils;
use Galacticus::Build::SourceTree::Process::FunctionClass::Utils qw(lctrim trimlc);
use Exporter 'import';

our @EXPORT_OK = qw(potentialDescriptorParameters);

# Alias the shared state storables database from the parent package.
our $stateStorables;
*stateStorables = \$Galacticus::Build::SourceTree::Process::FunctionClass::stateStorables;

sub potentialDescriptorParameters {
    # Process variable declarations for potential parameters to include in descriptors.
    my  $declarations     = shift();
    my  $nonAbstractClass = shift();
    my  $class            = shift();
    my  $potentialNames   = shift();
    foreach my $declaration ( &List::ExtraUtils::as_array($declarations) ) {
	# Identify object pointers.
	push(@{$potentialNames->{'objects'}},map {$_ =~ s/\s*([a-zA-Z0-9_]+).*/$1/; $_} @{$declaration->{'variables'}})
	    if
	    (
	     $declaration->{'intrinsic'} eq "class"
	     &&
	     (
	      (grep {&lctrim($declaration->{'type'}) eq lc($_)} keys                       (%{$stateStorables->{'functionClasses'       }}))
	      ||
	      (grep {&lctrim($declaration->{'type'}) eq lc($_)} &List::ExtraUtils::as_array(  $stateStorables->{'functionClassInstances'} ))
	     )
	     &&
	     grep {$_ eq "pointer"} @{$declaration->{'attributes'}}
	    );
	# Identify stateful types.
	push(@{$potentialNames->{'statefulTypes'}},$declaration)
	    if
	    (
	     $declaration->{'intrinsic'} eq "type"
	     &&
	     $declaration->{'type'     } =~ m/^stateful(Integer|Double|Logical)\s*$/i
	    );
	# Identify enumerations.
	push(@{$potentialNames->{'enumerations'}},$declaration)
	    if
	    (
	     $declaration->{'intrinsic'} eq "type"
	     &&
	     $declaration->{'type'     } =~ m/^enumeration[a-z0-9_]+type\s*$/i
	    );
	# Identify regular parameters.
	push(@{$potentialNames->{'parameters'}},$declaration)
	    if
	    (
	     (grep {$_ eq $declaration->{'intrinsic'}} ( "integer", "logical", "double precision", "character" ))
	     ||
	     (
	      $declaration->{'intrinsic'}  eq "type"
	      &&
	      trimlc($declaration->{'type'     }) eq "varying_string"
	     )
	    );
	$nonAbstractClass->{'hasCustomDescriptor'} = 1
	    if
	    (
	     $declaration->{'intrinsic'} eq "procedure"
	     &&
	     $declaration->{'variables'}->[0] =~ m/^descriptor=>/
	    );
    }
    # Identify linked lists parameters.
    if ( defined($class) ) {
	if ( exists($class->{'linkedList'}) ) {
	    foreach my $object ( split(" ",$class->{'linkedList'}->{'object'}) ) {
		push(@{$potentialNames->{'linkedListObjects'}},$object)
		    unless ( grep {$_ eq $object} @{$potentialNames->{'linkedListObjects'}} );
		$potentialNames->{'linkedLists'}->{$object} = $class->{'linkedList'};
	    }
	}
    }
}

1;
