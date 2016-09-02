# Contains a Perl module which handles deferred attributes of properties.

package Galacticus::Build::Components::Properties::Deferred;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::NullFunctions qw(createNullFunction);

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     propertiesDeferred => 
     {
	 propertyIteratedFunctions =>
	     [
	      \&Properties_Deferred_Pointers
	     ]
     }
    );

# Record of previously created pointers.
my %createdPointers;

sub Properties_Deferred_Pointers {
    # Generate procedure pointers for deferred attributes of properties.
    my $build    = shift();
    my $class    = shift();
    my $member   = shift();
    my $property = shift();
    # Skip properties which are not deferred.
    return
	if ( $property->{'attributes' }->{'isDeferred'} eq "" );
    # Determine the type of object to which this function will be bound.
    my $selfType = $property->{'attributes'}->{'bindsTo'} eq "top" ? "generic"        : $class->{'name'}                           ;
    # Determine where to attach.
    my $attachTo = $property->{'attributes'}->{'bindsTo'} eq "top" ? $class->{'name'} : $class->{'name'}.ucfirst($member->{'name'});
    # Adjectives for attributes.
    my %attributeAdjective =
	(
	 get  => "isGettable" ,
	 set  => "isSettable" ,
	 rate => "isEvolvable"
	);
    # Iterate over attributes.
    foreach my $attribute ( split(/:/,$property->{'attributes' }->{'isDeferred'}) ) {	
	# Determine function name.
	my $functionLabel = lcfirst($attachTo).ucfirst($property->{'name'}).ucfirst($attribute);
	# Determine if this attribute is deferred and has not yet had a procedure pointer created.
	next
	    unless ( exists($attributeAdjective{$attribute}) && $property->{'attributes' }->{$attributeAdjective{$attribute}} && ! exists($createdPointers{$functionLabel}) );
	# Construct the template function.
	my $template =
	    $attribute eq "get" 
	    ?
	    $class->{'name'}.ucfirst($member->{'name'}).ucfirst($property->{'name'}).ucfirst($attribute)
	    :
	    &createNullFunction($build,{selfType => $selfType, attribute => $attribute, property => $property, intent => "inout"});
	# Generate the procedure pointer and a boolean to indicate if is has been attached.
	push(
	    @{$build->{'variables'}},
	    {
		intrinsic  => "procedure",
		type       => $template,
		attributes => [ "pointer" ],
		variables  => [ $functionLabel."Deferred" ]
	    },
	    {
		intrinsic  => "logical",
		variables  => [ $functionLabel."IsAttachedValue=.false." ]
	    },
	    );
	# Record that this procedure pointer has been created.
	$createdPointers{$functionLabel} = 1;
    }
}

1;
