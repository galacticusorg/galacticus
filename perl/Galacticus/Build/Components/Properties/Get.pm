# Contains a Perl module which handles getting of component properties.

package Galacticus::Build::Components::Properties::Get;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::Properties::Utils;
use Text::Template 'fill_in_string';

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     propertiesGet => 
     {
	 propertyIteratedFunctions =>
	     [
	      \&Bind_Get_Functions ,
	      \&Build_Get_Functions
	     ]
     }
    );

sub Bind_Get_Functions {
    # Bind compile-time specified custom get functions which bind at the component level to the component implementation.
    my $build    = shift();
    my $class    = shift();
    my $member   = shift();
    my $property = shift();
    push(
	@{$build->{'types'}->{'nodeComponent'.ucfirst($class->{'name'}).ucfirst($member->{'name'})}->{'boundFunctions'}},
	{type => "procedure", name => $property->{'name'}, function => $property->{'getFunction'}->{'content'}}
	)
	if (
	                                   $property->{'attributes' }->{'isGettable'}  &&
	    !                              $property->{'getFunction'}->{'build'     }  &&
	    ! grep {$_ eq "get"} split(":",$property->{'attributes' }->{'isDeferred'})
	)
}

sub Build_Get_Functions {
    # Build get functions for non-deferred properties which bind at the component level.
    my $build       = shift();
    my $class       = shift();
    my $member      = shift();
    $code::property = shift();
    # Skip this property if it is not gettable, or if a compile-time custom get function was specified.
    return
	if
	(
	    $code::property->{'attributes' }->{'isVirtual' }
	 ||
	  ! $code::property->{'attributes' }->{'isGettable'}
	 ||
	  ! $code::property->{'getFunction'}->{'build'     }   
	);
    # Determine a suffix for the function and method names. If the get attribute of this property is deferred, a suffix of "Value"
    # is used - this can then be accessed by the deferred functions if necessary.
    my $suffix   = (grep {$_ eq "get"} split(":",$code::property->{'attributes' }->{'isDeferred'})) ? "Value" : "";
    # Determine the return type of the function.
    (my $functionTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'data'});
    my $functionType = 
	                                                                    $functionTypeDescriptor->{'intrinsic' }            .
	(exists($functionTypeDescriptor->{'type'      }) ? "(" .            $functionTypeDescriptor->{'type'      }  .")" : "").
	(exists($functionTypeDescriptor->{'attributes'}) ? ", ".join(", ",@{$functionTypeDescriptor->{'attributes'}})     : "");
    # Build the function.
    my $implementationTypeName = "nodeComponent".ucfirst($class->{'name'}).ucfirst($member->{'name'});
    my $function =
    {
	type        => $functionType." => propertyValue",
	name        => $class->{'name'}.ucfirst($member->{'name'}).ucfirst($code::property->{'name'})."Get".$suffix,
	description => "Get the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of an {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	variables   =>
	[
	 {
	     intrinsic  => "class",
	     type       => $implementationTypeName,
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 }
	]
    };
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
propertyValue=self%{$property->{'name'}}Data
CODE
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::property->{'name'}.$suffix
	}
	);  
}

1;
