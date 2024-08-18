# Contains a Perl module which provides handling of meta-properties for component classes.

package Galacticus::Build::Components::Classes::MetaProperties;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classMetaProperties =>
     {
	 classIteratedFunctions =>
	     [
	      \&Class_Meta_Property_Get,
	      \&Class_Meta_Property_Set
	     ]
     }
    );

# Define types of meta-property to add.
our @metaPropertyTypes =
    (
     {
	 label     => "float"           ,
	 intrinsic => "double precision",
	 rank      => 0
     },
     {
	 label     => "float"           ,
	 intrinsic => "double precision",
	 rank      => 1
     },
     {
	 label     => "longInteger"     ,
	 intrinsic => "integer"         ,
	 type      => "kind_int8"       ,
	 rank      => 0
     },
     {
	 label     => "longInteger"     ,
	 intrinsic => "integer"         ,
	 type      => "kind_int8"       ,
	 rank      => 1
     },
     {
	 label     => "integer"         ,
	 intrinsic => "integer"         ,
	 rank      => 0
     },
     {
	 label     => "integer"         ,
	 intrinsic => "integer"         ,
	 rank      => 1
     }
    );

sub Class_Meta_Property_Get {
    # Generate a function to return the value of an indexed meta-property.
    my $build             = shift();
    $code::class          = shift();
    $code::classTypeName  = "nodeComponent".ucfirst($code::class->{'name'});
    # Iterate over meta-property types.
    foreach my $metaPropertyType ( @metaPropertyTypes ) {
	my $content = "";
	if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	    $content              = "if (.not.".$code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyCreator(metaPropertyID)) call metaPropertyNoCreator('".$code::class->{'name'}."',char(".$code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyLabels(metaPropertyID)),'".$metaPropertyType->{'label'}."',".$metaPropertyType->{'rank'}.")\n";
	    $content             .= "allocate(value_(".join(",",map {"size(self%".$metaPropertyType->{'label'}."Rank".$metaPropertyType->{'rank'}."MetaProperties(metaPropertyID)%values,dim=".$_.")"} 1..$metaPropertyType->{'rank'})."))\n"
		if ($metaPropertyType->{'rank'} > 0);
	    $content             .= "value_=self%".$metaPropertyType->{'label'}."Rank".$metaPropertyType->{'rank'}."MetaProperties(metaPropertyID)".($metaPropertyType->{'rank'} > 0 ? "%values" : "")."\n";
	}
	my $function =
	{
	    type        => $metaPropertyType->{'intrinsic'}.(exists($metaPropertyType->{'type'}) ? "(".$metaPropertyType->{'type'}.")" : "").($metaPropertyType->{'rank'} > 0 ? ", allocatable, dimension(".join(",",":" x $metaPropertyType->{'rank'}).")" : "")." => value_",
	    name        => $code::classTypeName.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyGet",
	    description => "Return the value of a rank-".$metaPropertyType->{'rank'}." ".$metaPropertyType->{'label'}." meta-property of a ".$code::classTypeName." component given its ID.",
	    modules     =>
		[
		 "ISO_Varying_String"
		],
	    variables   =>
		[
		 {
		     intrinsic  => "class",
		     type       => $code::classTypeName,
		     attributes => [ "intent(in   )" ],
		     variables  => [ "self" ]
		 },
		 {
		     intrinsic  => "integer",
		     attributes => [ "intent(in   )" ],
		     variables  => [ "metaPropertyID" ]
		 }
		],
	    content     => $content
	};
	# Insert a type-binding for this function.
	push(
	    @{$build->{'types'}->{$code::classTypeName}->{'boundFunctions'}},
	    {
		type        => "procedure",
		descriptor  => $function,
		name        => $metaPropertyType->{'label'}."Rank".$metaPropertyType->{'rank'}."MetaPropertyGet", 
	    }
	    );
    }
}

sub Class_Meta_Property_Set {
    # Generate a function to set the value of an indexed meta-property.
    my $build             = shift();
    $code::class          = shift();
    $code::classTypeName  = "nodeComponent".ucfirst($code::class->{'name'});
    # Iterate over meta-property types.
    foreach my $metaPropertyType ( @metaPropertyTypes ) {
	my $content = "";
	if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	    $content              = "if (.not.".$code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyCreator(metaPropertyID)) call metaPropertyNoCreator('".$code::class->{'name'}."',char(".$code::class->{'name'}.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertyLabels(metaPropertyID)),'".$metaPropertyType->{'label'}."',".$metaPropertyType->{'rank'}.")\n";
	    $content             .= "self%".$metaPropertyType->{'label'}."Rank".$metaPropertyType->{'rank'}."MetaProperties(metaPropertyID)".($metaPropertyType->{'rank'} > 0 ? "%values" : "")."=metaPropertyValue\n";
	}
	my @attributes = ( "intent(in   )" );
	push(@attributes,"dimension(".join(",",":" x $metaPropertyType->{'rank'}).")")
	    if ( $metaPropertyType->{'rank'} > 0 );
	my $metaPropertyVariable =
	{
	     intrinsic  => $metaPropertyType->{'intrinsic'},
	     attributes => \@attributes,
	     variables  => [ "metaPropertyValue" ]
	};
	$metaPropertyVariable->{'type'} = $metaPropertyType->{'type'}
	     if ( exists($metaPropertyType->{'type'}) );
	my $function =
	{
	    type        => "void",
	    name        => $code::classTypeName.ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaPropertySet",
	    description => "Set the value of a rank-".$metaPropertyType->{'rank'}." ".$metaPropertyType->{'label'}." meta-property of a ".$code::classTypeName." component given its ID.",
	    modules     =>
		[
		 "ISO_Varying_String"
		],
	    variables   =>
		[
		 {
		     intrinsic  => "class",
		     type       => $code::classTypeName,
		     attributes => [ "intent(inout)" ],
		     variables  => [ "self" ]
		 },
		 {
		     intrinsic  => "integer",
		     attributes => [ "intent(in   )" ],
		     variables  => [ "metaPropertyID" ]
		 }
		],
	    content     => $content
	};
	push(@{$function->{'variables'}},$metaPropertyVariable);
	# Insert a type-binding for this function.
	push(
	    @{$build->{'types'}->{$code::classTypeName}->{'boundFunctions'}},
	    {
		type        => "procedure",
		descriptor  => $function,
		name        => $metaPropertyType->{'label'}."Rank".$metaPropertyType->{'rank'}."MetaPropertySet", 
	    }
	    );
    }
}

1;
