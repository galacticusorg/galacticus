# Contains a Perl module which handles building the base types of the class hierarchy during build.

package Galacticus::Build::Components::Components;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::Classes::MetaProperties;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     baseTypes => 
     {
	 types =>
	     [
	      \&Build_Node_Component_Class
	     ]
     }
    );

sub Build_Node_Component_Class {
    # Build the "nodeComponent" base class from which all other component classes and implementations inherit.
    my $build = shift();
    # Define type-bound functions.
    my @typeBoundFunctions = 
	(
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "type"                                                                                                 ,
	     function    => "Node_Component_Generic_Type"                                                                          ,
	     description => "Return the type of this object."                                                                      ,
	     returnType  => "\\textcolor{red}{\\textless type(varying\\_string)\\textgreater}"                                     ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "host"                                                                                                 ,
	     function    => "Node_Component_Host_Node"                                                                             ,
	     description => "Return a pointer to the host {\\normalfont \\ttfamily treeNode} object."                              ,
	     returnType  => "\\textcolor{red}{\\textless *type(treeNode)\\textgreater}"                                            ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "destroy"                                                                                              ,
	     function    => "Node_Component_Generic_Destroy"                                                                       ,
	     description => "Destroy the object."                                                                                  ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeCount"                                                                                       ,
	     function    => "Node_Component_Serialize_Count_Zero"                                                                  ,
	     description => "Return a count of the number of evolvable quantities to be evolved."                                  ,
	     returnType  => "\\intzero"                                                                                            ,
	     arguments   => "\\intzero\\ propertyType\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializationOffsets"                                                                                 ,
	     function    => "Node_Component_Serialization_Offsets"                                                                 ,
	     description => "Set offsets into serialization arrays."                                                               ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeValues"                                                                                      ,
	     function    => "Node_Component_Serialize_Null"                                                                        ,
	     description => "Serialize the evolvable quantities to an array."                                                      ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\doubleone\\ array\\argout, \\intzero\\ propertyType\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeNonNegative"                                                                                 ,
	     function    => "Node_Component_Serialize_NonNegative_Null"                                                            ,
	     description => "Serialize the non-negative status of evolvable quantities to an array."                               ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\logicalone\\ array\\argout"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "deserializeRaw"                                                                                       ,
	     function    => "Node_Component_Read_Raw_Null"                                                                         ,
	     description => "Read properties from raw file."                                                                       ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ fileHandle\\argin"
	 },	 
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "deserializeValues"                                                                                    ,
	     function    => "Node_Component_Deserialize_Null"                                                                      ,
	     description => "Deserialize the evolvable quantities from an array."                                                  ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\doubleone\\ array\\argin, \\intzero\\ propertyType\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "odeStepRatesInitialize"                                                                               ,
	     function    => "Node_Component_ODE_Step_Initialize_Null"                                                              ,
	     description => "Initialize rates for evolvable properties."                                                           ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "odeStepScalesInitialize"                                                                              ,
	     function    => "Node_Component_ODE_Step_Initialize_Null"                                                              ,
	     description => "Initialize scales for evolvable properties."                                                          ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""          
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeASCII"                                                                                       ,
	     function    => "Node_Component_Dump_Null"                                                                             ,
	     description => "Generate an ASCII dump of all properties."                                                            ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeXML"                                                                                         ,
	     function    => "Node_Component_Dump_XML_Null"                                                                         ,
	     description => "Generate an XML dump of all properties."                                                              ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeRaw"                                                                                         ,
	     function    => "Node_Component_Dump_Raw_Null"                                                                         ,
	     description => "Generate a binary dump of all properties."                                                            ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ fileHandle\\argin"                   
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "outputCount"                                                                                          ,
	     function    => "Node_Component_Output_Count_Null"                                                                     ,
	     description => "Compute a count of outputtable properties."                                                           ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ integerPropertyCount\\arginout, \\intzero\\ doublePropertyCount\\arginout, \\doublezero\\ time\\argin, \\intzero\\ instance\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "outputNames"                                                                                          ,
	     function    => "Node_Component_Output_Names_Null"                                                                     ,
	     description => "Generate names of outputtable properties."                                                            ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ integerProperty\\arginout, \\textcolor{red}{\\textless type(outputPropertyInteger)(:)\\textgreater} integerProperties\\arginout, \\intzero\\ doubleProperty\\arginout, \\textcolor{red}{\\textless type(otuputPropertyDouble)(:)\\textgreater} doubleProperties\\arginout, \\doublezero\\ time\\argin, \\intzero\\ instance\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "output"                                                                                               ,
	     function    => "Node_Component_Output_Null"                                                                           ,
	     description => "Generate values of outputtable properties."                                                           ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ integerProperty\\arginout, \\intzero\\ integerBufferCount\\arginout, \\textcolor{red}{\\textless type(outputPropertyInteger)(:)\\textgreater} integerProperties\\arginout, \\intzero doubleProperty\\arginout, \\intzero\\ doubleBufferCount\\arginout, \\textcolor{red}{\\textless type(outputPropertyDouble)(:)\\textgreater} doubleProperties\\arginout, \\doublezero\\ time\\argin, \\intzero\\ instance\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "massDistribution"                                                                                     ,
	     function    => "Node_Component_Mass_Distribution_Null"                                                                ,
	     description => "Return the mass distribution for this component."                                                     ,
	     returnType  => "\\textcolor{red}{\\textless class(massDistribution)\\textgreater}"                                    ,
	     arguments   => "\\textcolor{red}{\\textless type(enumerationComponentTypeType)\\textgreater} [componentType]\\argin, \\textcolor{red}{\\textless type(enumeratioMassTypeType)\\textgreater} [massType]\\argin, \\textcolor{red}{\\textless type(enumeratioWeightByType)\\textgreater} [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "massBaryonic"                                                                                         ,
	     function    => "Node_Component_Mass_Baryonic_Null"                                                                    ,
	     description => "Return the total baryonic mass for this component."                                                   ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "density"                                                                                              ,
	     function    => "Node_Component_Density_Null"                                                                          ,
	     description => "Compute the density."                                                                                 ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\textcolor{red}{\\textless double(3)\\textgreater} positionSpherical\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\enumWeightBy\\ [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "densitySphericalAverage"                                                                              ,
	     function    => "Node_Component_Density_Spherical_Average_Null"                                                        ,
	     description => "Compute the spherically-averaged density."                                                            ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\doublezero\\ radius\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\enumWeightBy\\ [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "surfaceDensity"                                                                                       ,
	     function    => "Node_Component_Surface_Density_Null"                                                                  ,
	     description => "Compute the surface density."                                                                         ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\textcolor{red}{\\textless double(3)\\textgreater} positionCylindrical\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\enumWeightBy\\ [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin"
	 }
	);
    # Add meta-property methods.
    foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	push(
	    @typeBoundFunctions,
	    {
		type        => "procedure"                                                                                                             ,
		name        => "add".ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaProperty"                           ,
		function    => "Node_Component_Generic_Add_".ucfirst($metaPropertyType->{'label'})."_Rank".$metaPropertyType->{'rank'}."_Meta_Property",
		description => "Add a rank-".$metaPropertyType->{'rank'}." ".$metaPropertyType->{'label'}." meta-property to this class."              ,
		returnType  => "\\intzero"                                                                                                             ,
		arguments   => "\\textcolor{red}{\\textless type(varying\_string)\\textgreater label, \\textcolor{red}{\\textless character(len=*)\\textgreater name".($metaPropertyType->{'label'} eq "flat" && $metaPropertyType->{'rank'} == 0 ? ", \\logicalzero\ [isEvolvable]" : "").", \\logicalzero\ [isCreator]"
	    }
	    );
    }
    # Specify the data content.
    my @dataContent =
	(
	 {
	     intrinsic  =>   "type",
	     type       =>   "treeNode"            ,
	     attributes => [ "pointer" , "public" ],
	     variables  => [ "hostNode => null()" ]
	 }
	);
    foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	my $typeDef;
	if ( $metaPropertyType->{'rank'} == 0 ) {
	    $typeDef->{'intrinsic'} = $metaPropertyType->{'intrinsic'};
	    $typeDef->{'type'     } = $metaPropertyType->{'type'     }
	        if ( exists($metaPropertyType->{'type'}) );
	} else {
	    $typeDef->{'intrinsic'} = "type";
	    $typeDef->{'type'     } = $metaPropertyType->{'label'}."Rank".$metaPropertyType->{'rank'}."MetaProperty";
	}
	$typeDef->{'attributes'} = [ "allocatable", "dimension(:)" ];
	$typeDef->{'variables' } = [ $metaPropertyType->{'label'}."Rank".$metaPropertyType->{'rank'}."MetaProperties" ];
	push(
	    @dataContent,
	    $typeDef
	    );
    }
    # Create the nodeComponent class.
    $build->{'types'}->{'nodeComponent'} = {
	name           => "nodeComponent"                           ,
	comment        => "A class for components in \\glspl{node}.",
	isPublic       => 1                                         ,
	boundFunctions => \@typeBoundFunctions,
	dataContent    => \@dataContent
    };
}

1;
