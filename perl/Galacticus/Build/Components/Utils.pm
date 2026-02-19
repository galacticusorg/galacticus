# Contains a Perl module which implements hooks for the component build system.

package Galacticus::Build::Components::Utils;
use strict;
use warnings;
use utf8;
use Data::Dumper;
use List::Util qw(max);
use Exporter qw(import);
our $VERSION = 1.00;
our @EXPORT_OK = qw($verbosityLevel @booleanLabel %intrinsicTypes %intrinsicNulls %outputTypeMap $classNameLengthMax $implementationNameLengthMax $fullyQualifiedNameLengthMax $propertyNameLengthMax $implementationPropertyNameLengthMax $linkedDataNameLengthMax applyDefaults isIntrinsic isOutputIntrinsic offsetName padClass padImplementation padFullyQualified padPropertyName padImplementationPropertyName padLinkedData);

# Define a hash into which modules can insert their hooks.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     utils => 
     {
	 gather =>
	     [
	      \&Label_Lengths
	     ]
     }
    );

# Global verbosity level.
our $verbosityLevel = 1;

# Boolean labels.
our @booleanLabel = ( "false", "true" );

# Intrinsic types.
our %intrinsicTypes =
    (
     "integer"     => "integer"                ,
     "longInteger" => "integer(kind=kind_int8)",
     "logical"     => "logical"                ,
     "double"      => "double precision"       ,
     "void"        => "void"
    );

# Null values for intrinsics.
our %intrinsicNulls =
    (
     "double"      => "0.0d0",
     "integer"     => "0",
     "longInteger" => "0_kind_int8",
     "logical"     => ".false."
    );

# Mapping for internal to output types.
our %outputTypeMap =
    (
     "double"      => "double" ,
     "integer"     => "integer",
     "longInteger" => "integer"
    );

# Maximum lengths of labels (used for formatting).
our $classNameLengthMax                     ;
our $implementationNameLengthMax            ; 
our $fullyQualifiedNameLengthMax            ;
our $propertyNameLengthMax                  ;
our $implementationPropertyNameLengthMax    ;
our $linkedDataNameLengthMax             = 0;

sub Label_Lengths {
    # Determine the lengths of various types of label for use in formatting.
    my $build = shift();
    # Find maximum lengths.
    $classNameLengthMax          = max map {                      length($_->{'class'})} &List::ExtraUtils::hashList($build->{'components'});
    $implementationNameLengthMax = max map {length($_->{'name' })                      } &List::ExtraUtils::hashList($build->{'components'});
    $fullyQualifiedNameLengthMax = max map {length($_->{'name' })+length($_->{'class'})} &List::ExtraUtils::hashList($build->{'components'});
    # Get property label lengths.
    $implementationPropertyNameLengthMax = 0;
    my $propertyNameLengthMax            = 0;
    foreach my $component ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	$propertyNameLengthMax               = max map {length($_->{'name' })} &List::ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' );
	$implementationPropertyNameLengthMax = 
	    max 
	    (
	     $implementationPropertyNameLengthMax,
	     map 
	     {
		 length($component->{'name' })
		  +
		 length($component->{'class'})
		  +
		 length($_        ->{'name' })
	     }
	     &List::ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' )
	    );
	$propertyNameLengthMax = 
	    max 
	    (
	     $propertyNameLengthMax,
	     map 
	     {
		 length($_        ->{'name' })
	     }
	     &List::ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' )
	    )
	    if ( exists($component->{'properties'}->{'property'}) );
    }
    # Add variables for use at run-time.
        push
	(
	 @{$build->{'variables'}},
	 {
	     intrinsic  => "integer"                                                                                                       ,
	     variables  => [ "propertyNameLengthMax=".$propertyNameLengthMax ]
	 }
	);
    # Report.
    if ( $verbosityLevel >= 1 ) {
	print "         --> Maximum label lengths:\n";
	print "            -->           Class: ".$classNameLengthMax         ."\n";
	print "            -->  Implementation: ".$implementationNameLengthMax."\n";
	print "            --> Fully-qualified: ".$fullyQualifiedNameLengthMax."\n";
	print "            -->        Property: ".$propertyNameLengthMax      ."\n";
    }
}

sub isIntrinsic {
    # Return true if the given type matches an intrinsic type.
    my $type = shift();
    return (grep {$_ eq $type} keys(%intrinsicTypes)) == 1 ? 1 : 0;
}

sub isOutputIntrinsic {
    # Return true if the given type matches an outputtable intrinsic type.
    my $type = shift();
    return (grep {$_ eq $type} ("double","integer","longInteger")) == 1 ? 1 : 0;
}

sub offsetName {
    # Return the name of the variable used to store the offset of a property into the ODE solver arrays.
    my %shortName =
	(
	 all      => "all",
	 active   => "atv",
	 inactive => "itv"
	);
    if ( scalar(@_) == 3 ) {
	my $status        = shift();
	my $componentName = shift();
	my $propertyName  = shift();
	die("unrecognized status")
	    unless ( $status eq "all" || $status eq "active" || $status eq "inactive" );
	return "offset".ucfirst($shortName{$status}).ucfirst($componentName).ucfirst($propertyName);
    } elsif ( scalar(@_) == 4 ) {
	my $status   = shift();
	my $class    = shift();
	my $member   = shift();
	my $property = shift();
	die("unrecognized status")
	    unless ( $status eq "all" || $status eq "active" || $status eq "inactive" );
	return "offset".ucfirst($shortName{$status}).ucfirst($class->{'name'}).ucfirst($member->{'name'}).ucfirst($property->{'name'});
    } else {
	die("offsetName(): incorrect number of arguments");
    }
}

sub padClass {
    # Pad a class name to give nicely aligned formatting in the output code.
    return &padGeneric($classNameLengthMax,@_);
}

sub padImplementation {
    # Pad an implementation name to give nicely aligned formatting in the output code.
    return &padGeneric($implementationNameLengthMax,@_);
}

sub padFullyQualified {
    # Pad a fully-qualified name to give nicely aligned formatting in the output code.
    return &padGeneric($fullyQualifiedNameLengthMax,@_);
}

sub padPropertyName {
    # Pad a property name to give nicely aligned formatting in the output code.
    return &padGeneric($propertyNameLengthMax,@_);
}

sub padImplementationPropertyName {
    # Pad a implementation + property name to give nicely aligned formatting in the output code.
    return &padGeneric($implementationPropertyNameLengthMax,@_);
}

sub padLinkedData {
    # Pad a linked data name to give nicely aligned formatting in the output code.
    return &padGeneric($linkedDataNameLengthMax,@_);
}

sub padGeneric {
    # Pad a string to give nicely aligned formatting in the output code.
    my $length     = shift;
    my $text       = shift;
    my @extraPad   = @{$_[0]};
    my $padLength  = $length+$extraPad[0];
    $padLength     = $extraPad[1]
	if ($extraPad[1] > $padLength);
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

sub applyDefaults {
    # Applies a set of default values to a nested data structure using recursion.
    my $object  = shift();
    my $name    = shift();
    my $default = shift();
    if ( ref($default) ) {
	# The default we've been passed is actually a list of deeper default settings. Iterate over them if the currently named
	# element exists ion the current data structure and apply them by calling ourself recursively.
	if ( exists($object->{$name}) ) {
	    foreach my $subObject ( &List::ExtraUtils::as_array($object->{$name}) ) {
		&applyDefaults($subObject,$_,$default->{$_})
		    foreach ( keys(%{$default}) );
	    }
	# Alternatively if the special name "ALL" is used, then iterate over all members of the object.    
	} elsif ( $name eq "ALL" ) {
	    foreach my $subObject ( &List::ExtraUtils::hashList($object) ) {
		&applyDefaults($subObject,$_,$default->{$_})
		    foreach ( keys(%{$default}) );
	    }	    
	}
    } else {
	# We've been passed an actual default. Iterate over all named objects in our data structure and apply the default if
	# necessary.
	foreach ( &List::ExtraUtils::as_array($object) ) {
	    if ( $default =~ m/^boolean/ ) {
		# In the case of a boolean default, we also translate any preset value into the associated boolean.
		if ( exists($_->{$name}) ) {
		    $_->{$name} = $_->{$name} eq        "true" ? 1 : 0;
		} else {
		    $_->{$name} = $default    eq "booleanTrue" ? 1 : 0;
		}
	    } else {
		# A regular (non-boolean) default.
		$_->{$name} = $default
		    unless ( exists($_->{$name}) );
	    }
	}
    }
}

sub argumentList {
    # Return a list of all arguments extracted from a list of variable descriptors.
    my @descriptors = @_;
    my @arguments;
    foreach my $descriptor ( @descriptors ) {
	push(
	    @arguments,
	    @{$descriptor->{'variables'}}
	    )
	    if (
		# Regular variables with intent.
		(
		 exists($descriptor->{'attributes'})
		 &&
		 grep {$_ =~ m/^intent\s*\(\s*(in|inout|out)\s*\)/} @{$descriptor->{'attributes'}}		      
		)
		||
		# Procedure pointers
		(
		 $descriptor->{'intrinsic'} eq "procedure"
		 &&
		 (
		  (
		   exists($descriptor->{'attributes'})
		   &&
		   grep {$_ eq "pointer"} @{$descriptor->{'attributes'}}
		  )
		  ||
		  (
		   exists($descriptor->{'isArgument'})
		   &&
		   $descriptor->{'isArgument'}
		  )
		 )		 
		)
		||
		# External variables.
		$descriptor->{'intrinsic'} eq "external"
	    );
    }
    return @arguments;
}

1;
