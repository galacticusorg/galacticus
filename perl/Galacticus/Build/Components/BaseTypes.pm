# Contains a Perl module which handles building the base types of the class hierarchy during build.

package BaseTypes;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;
require Galacticus::Build::Components::Utils;

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
	     arguments   => ""
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
	     arguments   => "\\doubleone\\ array\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "readRaw"                                                                                              ,
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
	     arguments   => "\\doubleone\\ array\\argout"
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
	     name        => "dump"                                                                                                 ,
	     function    => "Node_Component_Dump_Null"                                                                             ,
	     description => "Generate an ASCII dump of all properties."                                                            ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "dumpXML"                                                                                              ,
	     function    => "Node_Component_Dump_XML_Null"                                                                         ,
	     description => "Generate an XML dump of all properties."                                                              ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "dumpRaw"                                                                                              ,
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
	     arguments   => "\\intzero\\ integerProperty\\arginout, \\textcolor{red}{\\textless char[*](:)\\textgreater} integerPropertyNames\\arginout, \\textcolor{red}{\\textless char[*](:)\\textgreater} integerPropertyComments\\arginout, \\doubleone\\ integerPropertyUnitsSI\\arginout, \\intzero\\ doubleProperty\\arginout, \\textcolor{red}{\\textless char[*](:)\\textgreater} doublePropertyNames\\arginout, \\textcolor{red}{\\textless char[*](:)\\textgreater} doublePropertyComments\\arginout, \\doubleone\\ doublePropertyUnitsSI\\arginout, \\doublezero\\ time\\argin, \\intzero\\ instance\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "output"                                                                                               ,
	     function    => "Node_Component_Output_Null"                                                                           ,
	     description => "Generate values of outputtable properties."                                                           ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ integerProperty\\arginout, \\intzero\\ integerBufferCount\\arginout, \\inttwo\\ integerBuffer\\arginout, \\intzero doubleProperty\\arginout, \\intzero\\ doubleBufferCount\\arginout, \\doubletwo\\ doubleBuffer\\arginout, \\doublezero\\ time\\argin, \\intzero\\ instance\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "enclosedMass"                                                                                         ,
	     function    => "Node_Component_Enclosed_Mass_Null"                                                                    ,
	     description => "Compute the mass enclosed within a radius."                                                           ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\doublezero\\ radius\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\enumWeightBy\\ [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin, \\logicalzero\\ [haloLoaded]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "density"                                                                                              ,
	     function    => "Node_Component_Density_Null"                                                                          ,
	     description => "Compute the density."                                                                                 ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\textcolor{red}{\\textless double(3)\\textgreater} positionSpherical\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\enumWeightBy\\ [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin, \\logicalzero\\ [haloLoaded]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "surfaceDensity"                                                                                       ,
	     function    => "Node_Component_Surface_Density_Null"                                                                  ,
	     description => "Compute the surface density."                                                                         ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\textcolor{red}{\\textless double(3)\\textgreater} positionCylindrical\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\enumWeightBy\\ [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin, \\logicalzero\\ [haloLoaded]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "potential"                                                                                            ,
	     function    => "Node_Component_Potential_Null"                                                                        ,
	     description => "Compute the gravitational potential."                                                                 ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\doublezero\\ radius\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\logicalzero\\ [haloLoaded]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "rotationCurve"                                                                                        ,
	     function    => "Node_Component_Rotation_Curve_Null"                                                                   ,
	     description => "Compute the rotation curve."                                                                          ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\doublezero\\ radius\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\logicalzero\\ [haloLoaded]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "rotationCurveGradient"                                                                                ,
	     function    => "Node_Component_Rotation_Curve_Gradient_Null"                                                          ,
	     description => "Compute the rotation curve gradient."                                                                 ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\doublezero\\ radius\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\logicalzero\\ [haloLoaded]\\argin"
	 }
	);
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
