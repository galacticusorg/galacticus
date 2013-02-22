# Contains a Perl module which implements processing of directives associated with the generation of C bindings in the Galacticus
# build system.

package BindingsC;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;
use DateTime;
use Data::Dumper;
use Switch;
use Scalar::Util 'reftype';
require Galacticus::Build::Hooks;
require Galacticus::Build::Dependencies;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     cBinding   => {parse => \&CBinding_Parse_Directive  , generate => \&CBinding_Generate_Output  },
     cTemplate  => {parse => \&CTemplate_Parse_Directive , generate => \&CTemplate_Generate_Output },
     cInterface => {parse => \&CInterface_Parse_Directive, generate => \&CInterface_Generate_Output}
    );

sub CBinding_Parse_Directive {
    # Parse content for a "cBinding" directive.
    my $buildData = shift;

    # Assert that we have a currentDocument.
    die("Galacticus::Build::BindingsC::CBinding_Parse_Directive: no currentDocument present"  )
	unless ( exists($buildData->{'currentDocument'} ) );

    # Extract the unit name.
    my $unitName = $buildData->{'currentDocument'}->{'unitName'};

    # Construct the code for the binding.
    $buildData->{'cBinding'}->{'code'} .= "bind(c,name='".$unitName."') :: ".$unitName."\n"
	if ( $buildData->{'codeType'} eq "c" );

}

sub CBinding_Generate_Output {
    # Generate output for a "cBinding" directive.
    my $buildData = shift;

    # Assert that we have a file name and directive present.
    die("Galacticus::Build::BindingsC::CBinding_Generate_Output: no fileName present" )
	unless ( exists($buildData->{'fileName'}) );
    die("Galacticus::Build::BindingsC::CBinding_Generate_Output: no directive present")
	unless ( exists($buildData->{'directive'}) );

    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::BindingsC\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "!  Time: ".$now."\n";

    # Iterate over all labels, and add them to the content.
    $buildData->{'content'} .= $buildData->{'cBinding'}->{'code'};
}

sub CTemplate_Parse_Directive {
    # Parse component declarations to build C-binding templates.
    my $buildData = shift;

    # Get the class name.
    my $className = $buildData->{'currentDocument'}->{'class'};

    # Construct a record of component classes.
    $buildData->{'componentClasses'}->{$className}->{'name'} = $className;

    # Store all properties.
    foreach my $propertyName ( keys(%{$buildData->{'currentDocument'}->{'properties'}->{'property'}}) ) {
	my $property = $buildData->{'currentDocument'}->{'properties'}->{'property'}->{$propertyName};
	if ( exists($buildData->{'componentClasses'}->{$className}->{'properties'}->{'property'}->{$propertyName}) ) {
	    if (
		exists($property->{$propertyName}->{'attributes'}->{'isGettable'})          &&
		       $property->{$propertyName}->{'attributes'}->{'isGettable'} eq "true" 
		) {
		$buildData->{'componentClasses'}->{$className}->{'properties'}->{'property'}->{$propertyName}->{'attributes'}->{'isGettable'} = "true"
	    }
	} else {
	    $buildData->{'componentClasses'}->{$className}->{'properties'}->{'property'}->{$propertyName} = $property;
	}
    }
}

sub CTemplate_Generate_Output {
    # Generate output for a "cTemplate" directive.
    my $buildData = shift;

    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;

    # Add a header.
    $buildData->{'content'}  = "// Generated automatically by Galacticus::Build::BindingsC\n";
    $buildData->{'content'} .= "//  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "//  Time: ".$now."\n\n";

    # Initialize external declarations.
    my $externalDeclarations = "extern \"C\"\n";
    $externalDeclarations   .= "{\n";
    
    # Insert class definitions.
    my $classDefinitions;
    foreach my $componentClassName ( keys(%{$buildData->{'componentClasses'}}) ) {
	my $componentClass = $buildData->{'componentClasses'}->{$componentClassName};
	# Insert an external declaration for the get function.
	$externalDeclarations .= "  void* Node_Component_".ucfirst($componentClassName)."_Get(void *thisNode);\n";
	# Build the class.
	$classDefinitions .= "// Define a class for ".$componentClassName." components.\n";
	$classDefinitions .= "class nodeComponent".ucfirst($componentClassName)."\n";
	$classDefinitions .= "{\n";
	$classDefinitions .= "   void *selfNode;\n";
	$classDefinitions .= "public:\n";
	$classDefinitions .= " nodeComponent".ucfirst($componentClassName)."(void *thisNode);\n";
	# Create methods for scalar real properties.
	foreach my $propertyName ( keys(%{$componentClass->{'properties'}->{'property'}}) ) {
	    my $property = $componentClass->{'properties'}->{'property'}->{$propertyName};
	    if (
		$property->{'type'      }                 eq "real" &&
		$property->{'rank'      }                 ==      0 &&
		$property->{'attributes'}->{'isGettable'} eq "true"
		) {
		$classDefinitions .= "  double ".$propertyName."() {return Node_Component_".ucfirst($componentClassName)."_".ucfirst($propertyName)."(selfNode);}\n";
		# Insert an external declaration for the get function.
		$externalDeclarations .= "  double Node_Component_".ucfirst($componentClassName)."_".ucfirst($propertyName)."(void *selfNode);\n";
	    }
	}
	$classDefinitions .= "};\n\n";
	# Create a constructor function.
	$classDefinitions .= "// Define a constructor for ".$componentClassName." components.\n";
	$classDefinitions .= "nodeComponent".ucfirst($componentClassName)."::nodeComponent".ucfirst($componentClassName)." (void *thisNode) {\n";
	$classDefinitions .= "  selfNode = thisNode;\n";
	$classDefinitions .= "};\n\n";
    }

    # Close external declarations.
    $externalDeclarations .= "}\n";

    # Insert the content.
    $buildData->{'content'} .= $externalDeclarations;
    $buildData->{'content'} .= $classDefinitions;
}

sub CInterface_Parse_Directive {
    # Parse component declarations to build C-binding interfaces.
    my $buildData = shift;

    # Get the class name.
    my $className = $buildData->{'currentDocument'}->{'class'};

    # Construct a record of component classes.
    $buildData->{'componentClasses'}->{$className}->{'name'} = $className;

    # Store all properties.
    foreach my $propertyName ( keys(%{$buildData->{'currentDocument'}->{'properties'}->{'property'}}) ) {
	my $property = $buildData->{'currentDocument'}->{'properties'}->{'property'}->{$propertyName};
	if ( exists($buildData->{'componentClasses'}->{$className}->{'properties'}->{'property'}->{$propertyName}) ) {
	    if (
		exists($property->{$propertyName}->{'attributes'}->{'isGettable'})          &&
		       $property->{$propertyName}->{'attributes'}->{'isGettable'} eq "true" 
		) {
		$buildData->{'componentClasses'}->{$className}->{'properties'}->{'property'}->{$propertyName}->{'attributes'}->{'isGettable'} = "true"
	    }
	} else {
	    $buildData->{'componentClasses'}->{$className}->{'properties'}->{'property'}->{$propertyName} = $property;
	}
    }
}

sub CInterface_Generate_Output {
    # Generate output for a "cInterface" directive.
    my $buildData = shift;

    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::BindingsC\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "!  Time: ".$now."\n\n";

    # Iterate over all component classes.
    foreach my $componentClassName ( keys(%{$buildData->{'componentClasses'}}) ) {
	my $componentClass = $buildData->{'componentClasses'}->{$componentClassName};
	# Create interfaces for scalar real properties.
	foreach my $propertyName ( keys(%{$componentClass->{'properties'}->{'property'}}) ) {
	    my $property = $componentClass->{'properties'}->{'property'}->{$propertyName};
	    if (
		$property->{'type'      }                 eq "real" &&
		$property->{'rank'      }                 ==      0 &&
		$property->{'attributes'}->{'isGettable'} eq "true"
		) {
		$buildData->{'content'} .= "function cNode_Component_".ucfirst($componentClassName)."_".ucfirst($propertyName)."(cSelfNode) bind(c,name='Node_Component_".ucfirst($componentClassName)."_".ucfirst($propertyName)."')\n";
		$buildData->{'content'} .= "  use, intrinsic :: ISO_C_Binding\n";
		$buildData->{'content'} .= "  use Galacticus_Nodes\n";
		$buildData->{'content'} .= "  implicit none\n";
		$buildData->{'content'} .= "  real(c_double)          :: cNode_Component_".ucfirst($componentClassName)."_".ucfirst($propertyName)."\n";
		$buildData->{'content'} .= "  type(c_ptr),    value   :: cSelfNode\n";
		$buildData->{'content'} .= "  type(treeNode), pointer :: selfNode\n";
		$buildData->{'content'} .= "  class(nodeComponent".ucfirst($componentClassName)."), pointer :: self\n\n";
		$buildData->{'content'} .= "  call c_f_pointer(cSelfNode,selfNode)\n";
		$buildData->{'content'} .= "  self => selfNode%".$componentClassName."()\n";
		$buildData->{'content'} .= "  cNode_Component_".ucfirst($componentClassName)."_".ucfirst($propertyName)."=self%".$propertyName."()\n";
		$buildData->{'content'} .= "  return\n";
		$buildData->{'content'} .= "end function cNode_Component_".ucfirst($componentClassName)."_".ucfirst($propertyName)."\n\n";
	    }
	}
    }
}

1;
