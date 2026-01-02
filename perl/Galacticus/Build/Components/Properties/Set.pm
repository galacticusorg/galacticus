# Contains a Perl module which handles setting of component properties during build.

package Galacticus::Build::Components::Properties::Set;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Text::Template 'fill_in_string';

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     propertiesSet => 
     {
	 functions =>
	     [
	      \&Build_Class_Setters
	     ],
	 propertyIteratedFunctions =>
	     [
	      \&Bind_Set_Functions ,
	      \&Build_Set_Functions
	     ]
     }
    );

sub Build_Class_Setters {
    # Insert settability functions for each property into the associated component class.
    my $build = shift();
    # Iterate over components.
    foreach my $component ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	foreach my $property ( &List::ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    # Insert an "isSettable" function into the base class.
	    my $functionName = $property->{'name'}."IsSettable";
	    unless ( grep {$_->{'name'} eq $functionName} @{$build->{'types'}->{'nodeComponent'.ucfirst($component->{'class'})}->{'boundFunctions'}} ) {
	    	push(
	    	    @{$build->{'types'}->{'nodeComponent'.ucfirst($component->{'class'})}->{'boundFunctions'}},
	    	    {
	    		type        => "procedure"    ,
	    		pass        => "nopass"       ,
	    		name        => $functionName  , 
	    		function    => "Boolean_False",
	    		returnType  => "\\logicalzero",
	    		arguments   => ""             ,
	    		description => "Specify whether the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$component->{'class'}."} component is settable."
	    	    }
	    	    );
	    }
	}
    }
}

sub Bind_Set_Functions {
    # Bind compile-time specified custom set functions which bind at the component level to the component implementation.
    my $build    = shift();
    my $class    = shift();
    my $member   = shift();
    my $property = shift();
    push(
	@{$build->{'types'}->{'nodeComponent'.ucfirst($class->{'name'}).ucfirst($member->{'name'})}->{'boundFunctions'}},
	{type => "procedure", name => $property->{'name'}."Set", function => $property->{'setFunction'}->{'content'}}
	)
	if (
	                                   $property->{'attributes' }->{'isGettable'}  &&
	    !                              $property->{'setFunction'}->{'build'     }  &&
	    ! grep {$_ eq "set"} split(":",$property->{'attributes' }->{'isDeferred'})
	)
}

sub Build_Set_Functions {
    # Build set functions for non-deferred properties which bind at the component level.
    my $build       = shift();
    my $class       = shift();
    my $member      = shift();
    $code::property = shift();
    # Skip this property if it is not settable, or if a compile-time custom get function was specified.
    return
	if
	(
	    $code::property->{'attributes' }->{'isVirtual' }
	 ||
	  ! $code::property->{'attributes' }->{'isSettable'}
	 ||
	  ! $code::property->{'setFunction'}->{'build'     }   
	);
    # Determine a suffix for the function and method names. If the set attribute of this property is deferred, a suffix of "Value"
    # is used - this can then be accessed by the deferred functions if necessary.
    my $suffix   = (grep {$_ eq "set"} split(":",$code::property->{'attributes' }->{'isDeferred'})) ? "Value" : "";
    # Determine the type of the property.
    (my $propertyTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'data'});
    @{$propertyTypeDescriptor->{'variables' }} = ( "setValue"      );
    @{$propertyTypeDescriptor->{'attributes'}} = ( "intent(in   )" );
    push(@{$propertyTypeDescriptor->{'attributes'}},"dimension(".join(",",(":") x $code::property->{'data'}->{'rank'}).")")
	if ( $code::property->{'data'}->{'rank'} > 0 );
    # Build the function.
    my $implementationTypeName = "nodeComponent".ucfirst($class->{'name'}).ucfirst($member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $class->{'name'}.ucfirst($member->{'name'}).ucfirst($code::property->{'name'})."Set".$suffix,
	description => "Set the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of an {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     $propertyTypeDescriptor
	    ]
    };
    # Build the function.
    if ( $code::property->{'data'}->{'rank'} == 0 ) {
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
self%{$property->{'name'}}Data=setValue
CODE
    } elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
if (.not.allocated(self%{$property->{'name'}}Data)) then
      !![
      <allocate variable="self%{$property->{'name'}}Data" size="setValue" rank="{$property->{'data'}->{'rank'}}"/>
      !!]
else
   if (size(self%{$property->{'name'}}Data) /= size(setValue)) then
      deallocate(self%{$property->{'name'}}Data)
      !![
      <allocate variable="self%{$property->{'name'}}Data" size="setValue" rank="{$property->{'data'}->{'rank'}}"/>
      !!]
   end if
end if
self%{$property->{'name'}}Data=setValue
CODE
    }
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::property->{'name'}."Set".$suffix
	}
	);  
}

1;
