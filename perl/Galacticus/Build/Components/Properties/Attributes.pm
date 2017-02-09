# Contains a Perl module which handles setting of component properties during build.

package Galacticus::Build::Components::Properties::Attributes;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Text::Template 'fill_in_string';
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::Properties::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     propertiesAttributes => 
     {
	 classIteratedFunctions =>
	     [
	      \&Attributes_Match
	     ]
     }
    );

sub Attributes_Match {
    # Generate functions which return lists of class implementations that match a given set of attribute requirements.
    my $build    = shift();
    $code::class = shift();
    # Initialize a structure of properties.
    my $properties;
    # Iterate over class member implementations.
    foreach my $member ( @{$code::class->{'members'}} ) {
	# Iterate over component and parents.
	my $parentMember = $member;
	while ( defined($parentMember) ) {
	    # Iterate over the properties of this implementation.
	    foreach my $property ( &List::ExtraUtils::hashList($parentMember->{'properties'}->{'property'}) ) {
		# Record attributes.
		$properties->{$property->{'name'}}->{'members'}->{$member->{'name'}}->{$_} = $property->{'attributes'}->{$Galacticus::Build::Components::Properties::Utils::attributeAdjective{$_}}
		   foreach ( "set", "get", "rate" ); 
	    }
	    $parentMember = exists($parentMember->{'extends'}) ? $parentMember->{'extends'}->{'implementation'} : undef();
	}
    }
    # Iterate over properties, creating a function for each.
    foreach $code::property ( &List::ExtraUtils::hashList($properties, keyAs => "name") ) {
	# Create the function.
	my $function =
	{
	    type        => "type(varying_string), allocatable, dimension(:) => matches",
	    name        => $code::class->{'name'}.ucfirst($code::property->{'name'})."AttributeMatch",
	    description => "Return a text list of component implementations in the {\\normalfont \\ttfamily ".$code::class->{'name'}."} class that have the desired attributes for the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property",
	    modules =>
		[
		 "ISO_Varying_String"
		],
	    variables =>
		[
		 {
		     intrinsic  => "logical",
		     attributes => [ "intent(in   )", "optional" ],
		     variables  => [ "requireSettable", "requireGettable", "requireEvolvable" ]
		 },
		 {
		     intrinsic  => "logical",
		     variables  => [ "requireSettableActual", "requireGettableActual", "requireEvolvableActual" ]
		 },
		 {
		     intrinsic  => "type",
		     type       => "varying_string",
		     attributes => [ "allocatable", "dimension(:)" ],

		     variables  => [ "temporaryList" ]
		 }
		]
	};
	# Create the function code.
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
requireSettableActual =.false.
requireGettableActual =.false.
requireEvolvableActual=.false.
if (present(requireSettable )) requireSettableActual =requireSettable
if (present(requireGettable )) requireGettableActual =requireGettable
if (present(requireEvolvable)) requireEvolvableActual=requireEvolvable
CODE
	# Iterate over component implementations.
	foreach $code::member ( &List::ExtraUtils::hashList($code::property->{'members'}, keyAs => 'name' ) ) {
	    @code::logic = ();
	    foreach ( "get", "set", "rate" ) {
		push(@code::logic,".not.require".($_ eq "rate" ? "Evolvable" : ucfirst($_)."table")."Actual" )
		    unless ( $code::member->{$_} );
	    }
	    if ( scalar(@code::logic) > 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if ({join(".and.",@logic)}) then
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(matches)) then
   call Move_Alloc(matches,temporaryList)
   allocate(matches(size(temporaryList)+1))
   matches(1:size(temporaryList))=temporaryList
   deallocate(temporaryList)
else
   allocate(matches(1))
end if
matches(size(matches))='{$member->{'name'}}'
CODE
	    if ( scalar(@code::logic) > 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
	     }
	}	    
	# Insert a type-binding for this function into the node component class type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		descriptor  => $function,
		pass        => "nopass",
		name        => $code::property->{'name'}."AttributeMatch"
	    }
	    );
    }
}

1;
