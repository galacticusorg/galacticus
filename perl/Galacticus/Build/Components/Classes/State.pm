# Contains a Perl module which provides state variables for component classes.

package Galacticus::Build::Components::Classes::State;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;
use Galacticus::Build::Components::Classes::MetaProperties;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classState =>
     {
	 classIteratedFunctions =>
	     [
	      \&Class_State  ,
	      \&Class_Size_Of
	     ]
     }
    );

sub Class_State {
    # Generate variables which store state for component classes.
    my $build    = shift();
    $code::class = shift();
    my $classTypeName = "nodeComponent".ucfirst($code::class->{'name'});
    push
	(
	 @{$build->{'variables'}},
	 # Object that will record which implementation of the class is to be used by default.
	 {
	     intrinsic  => "class"                                                  ,
	     type       => $classTypeName                                           ,
	     attributes => [ "public", "allocatable", "target" ]                    ,
	     variables  => [ "default".ucfirst($code::class->{'name'})."Component" ]
	 });
    push
	(
	 @{$build->{'variables'}},
	 # Object used to allocate component.
	 {
	     intrinsic  => "type"                                                   ,
	     type       => $classTypeName                                           ,
	     variables  => [ $code::class->{'name'}."Class" ]
	 }
	)
	if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );


    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	    push
		(
		 @{$build->{'variables'}},
		 # Arrays used to store meta-property indices.
		 {
		     intrinsic  => "type"                                                                                                                            ,
		     type       => "varying_string"                                                                                                                  ,
		     attributes => [ "allocatable", "dimension(:)" ]                                                                                                 ,
		     variables  => [ $code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyLabels"           ]
		 },
		 {
		     intrinsic  => "type"                                                                                                                            ,
		     type       => "varying_string"                                                                                                                  ,
		     attributes => [ "allocatable", "dimension(:)" ]                                                                                                 ,
		     variables  => [ $code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyNames"            ]
		 },
		 {
		     intrinsic  => "logical"                                                                                                                         ,
		     attributes => [ "allocatable", "dimension(:)" ]                                                                                                 ,
		     variables  => [ $code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyCreator"          ]
		 },
		 {
		     intrinsic  => "integer"                                                                                                                         ,
		     variables  => [ $code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyCount         =0" ]
		 }
		);
	    push
		(
		 @{$build->{'variables'}},
		 {
		     intrinsic  => "logical"                                                                                                                         ,
		     attributes => [ "allocatable", "dimension(:)" ]                                                                                                 ,
		     variables  => [ $code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyEvolvable"        ]
		 },
		 {
		     intrinsic  => "integer"                                                                                                                         ,
		     variables  => [ $code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyEvolvableCount=0" ]
		 }
		)
		if ( $metaPropertyType->{'label'} eq "float" && $metaPropertyType->{'rank'} == 0 );	    
	}
    }
}

sub Class_Size_Of {
    # Generate a function to return the size of the component implementation in bytes.
    my $build            = shift();
    $code::class         = shift();
    $code::classTypeName = "nodeComponent".ucfirst($code::class->{'name'});
    my $function =
    {
	type        => "integer(c_size_t)",
	name        => $code::classTypeName."SizeOf",
	description => "Return the size in bytes of a ".$code::classTypeName." component.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $code::classTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    ],
	content     => $code::classTypeName."SizeOf=sizeof(self)\n"
    };
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$code::classTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "sizeOf", 
	}
	);
}

1;
