# Contains a Perl module which handles the object hierarchy.

package Hierarchy;
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
use Data::Dumper;
use Text::Template 'fill_in_string';
require List::ExtraUtils;
require Fortran::Utils;
require Galacticus::Build::Components::Utils;
require Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     hierarchy => 
     {
	 functions =>
	     [
	      \&Hierarchy_Initialization
	     ]
     }
    );

sub Hierarchy_Initialization {
    # Generate a function to initialize the node class hierarchy.
    my $build = shift();
    # Generate the function code.
    my $functionCode;
    $functionCode .= "  subroutine Node_Class_Hierarchy_Initialize()\n";
    $functionCode .= "    !% Initialize the \\glc\\\ object system.\n";
    $functionCode .= "    use Input_Parameters\n";
    $functionCode .= "    use ISO_Varying_String\n";
    $functionCode .= "    use Memory_Management\n";
    $functionCode .= "    implicit none\n";


    # AJB BUILDING THE FUNCTION HERE
    my $function =
    {
	type        => "void",
	name        => "Node_Class_Hierarchy_Initialize",
	description => "Initialize the \\glc\\\ node/component class hierarchy.",
	modules     => 
	    [ 
	      "Input_Parameters"  ,
	      "ISO_Varying_String",
	      "Memory_Management" 
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "methodSelection", "message" ]
	     },
	     map
	     {
		 {
		     intrinsic  => "type",
		     type       => "nodeComponent".$_,
		     variables  => [ "default".$_."Component" ]
		 }
	     }
	     @{$build->{'componentIdList'}}
	    ],
    };
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
if (.not.moduleIsInitialized) then\n";
  !$omp critical (Node_Class_Hierarchy_Initialize)
  if (.not.moduleIsInitialized) then
    ! Parameters controlling output.
CODE
    foreach $code::className ( @{$build->{'componentClassList'}} ) {
	@code::members = sort(@{$build->{'componentClasses'}->{$code::className}->{'members'}});
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    ! Insert a function call to get the parameter controlling the choice of implementation for this class.
    !@ <inputParameter>
    !@   <name>treeNodeMethod".ucfirst($componentClass)."</name>
    !@   <defaultValue>".$defaultMethod."</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@    Specifies the implementation to be used for the ".$componentClass." component of nodes.
    !@   </description>
    !@   <type>string</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('treeNodeMethod".&Utils::padClass(ucfirst($componentClass)."'",[1,0]).",methodSelection,defaultValue='".&Utils::padImplementation($defaultMethod."'",[1,0]).")
    if (.not.allocated(default{&Utils::padClass(ucfirst($className)."Component",[9,0])})) then
      message='unrecognized method "'//methodSelection//'" for "{$className}" component'
      message=message//char(10)//'  available methods are:'
      {
       join(
            "\n      ",
            map 
             {"message=message//char(10)//'    ".$_."'"} 
             @members
           )
      }
      call Galacticus_Error_Report('Node_Class_Hierarchy_Initialize',message)
    end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    ! Record that the module is now initialized.
    moduleIsInitialized=.true.
  end if
  !$omp end critical (Node_Class_Hierarchy_Initialize)
end if\n";
CODE

    print Dumper($function);
    
    
    
    # Generate variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "type",
	     type       => "varying_string",
	     variables  => [ "methodSelection", "message" ]
	 },
	 map
	 {
	     {
		 intrinsic  => "type",
		 type       => "nodeComponent".$_,
		 variables  => [ "default".$_."Component" ]
	     }
	 }
	 @{$build->{'componentIdList'}}
	);
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Check for already initialized.
    $functionCode .= "   if (.not.moduleIsInitialized) then\n";
    $functionCode .= "      !\$omp critical (Galacticus_Nodes_Initialize)\n";
    $functionCode .= "      if (.not.moduleIsInitialized) then\n";
    # Record of output conditions seen.
    my %outputConditions;
    # Iterate over all component classes.
    $build->{'content'} .= "  ! Parameters controlling output.\n\n";
    foreach my $componentClass ( @{$build->{'componentClassList'}} ) {
	# Identify the default implementation.
    	my $defaultMethod;
    	foreach my $implementationName ( @{$build->{'componentClasses'}->{$componentClass}->{'members'}} ) {
    	    my $fullName = ucfirst($componentClass).ucfirst($implementationName);	    
    	    $defaultMethod = $implementationName
		if ( $build->{'components'}->{$fullName}->{'isDefault'} );
    	}
	die("Hierarchy_Initialization: no default method was found for '".$componentClass."' class")
	    unless ( defined($defaultMethod) );
	# Insert a function call to get the parameter controlling the choice of implementation for this class.
        $functionCode .= "       !@ <inputParameter>\n";
        $functionCode .= "       !@   <name>treeNodeMethod".ucfirst($componentClass)."</name>\n";
        $functionCode .= "       !@   <defaultValue>".$defaultMethod."</defaultValue>\n";
        $functionCode .= "       !@   <attachedTo>module</attachedTo>\n";
        $functionCode .= "       !@   <description>\n";
        $functionCode .= "       !@    Specifies the implementation to be used for the ".$componentClass." component of nodes.\n";
        $functionCode .= "       !@   </description>\n";
        $functionCode .= "       !@   <type>string</type>\n";
        $functionCode .= "       !@   <cardinality>1</cardinality>\n";
        $functionCode .= "       !@ </inputParameter>\n";
    	$functionCode .= "       call Get_Input_Parameter('treeNodeMethod".&Utils::padClass(ucfirst($componentClass)."'",[1,0]).",methodSelection,defaultValue='".&Utils::padImplementation($defaultMethod."'",[1,0]).")\n";
    	foreach my $implementationName ( @{$build->{'componentClasses'}->{$componentClass}->{'members'}} ) {
    	    my $fullName  = ucfirst($componentClass).ucfirst($implementationName);
	    my $component = $build->{'components'}->{$fullName};
    	    $functionCode .= "       if (methodSelection == '".&Utils::padImplementation($implementationName."'",[1,0]).") then\n";
	    $functionCode .= "          allocate(default".&Utils::padClass(ucfirst($componentClass)."Component",[9,0]).",source=default".&Utils::padFullyQualified($fullName."Component",[9,0]).")\n";
	    $functionCode .= "          nodeComponent".&Utils::padFullyQualified($fullName."IsActive",[8,0])."=.true.\n";
	    until ( $fullName eq "" ) {
		if ( exists($build->{'components'}->{$fullName}->{'extends'}) ) {
		    $fullName = ucfirst($build->{'components'}->{$fullName}->{'extends'}->{'class'}).ucfirst($build->{'components'}->{$fullName}->{'extends'}->{'name'});
		    $functionCode .= "          nodeComponent".&Utils::padFullyQualified($fullName."IsActive",[8,0])."=.true.\n";
		} else {
		    $fullName = "";
		}
	    }
	    $functionCode .= "    end if\n";
	    # Insert code to read and parameters controlling outputs.
	    foreach my $property ( &ExtraUtils::hashList($component->{'properties'}->{'property'}) ) {
		# Check for output and output condition.
		if (
		    exists($property->{'output'}               )                       &&
		    exists($property->{'output'}->{'condition'})                       &&
		    $property->{'output'}->{'condition'} =~ m/\[\[([^\]]+)\]\]/
		    )
		{		    
		    my $parameterName = $1;
		    unless ( exists($outputConditions{$parameterName}) ) {
			$functionCode .= "       !@ <inputParameter>\n";
			$functionCode .= "       !@   <name>".$parameterName."</name>\n";
			$functionCode .= "       !@   <defaultValue>false</defaultValue>\n";
			$functionCode .= "       !@   <attachedTo>module</attachedTo>\n";
			$functionCode .= "       !@   <description>\n";
			$functionCode .= "       !@    Specifies whether the {\\normalfont \\ttfamily ".$property->{'name'}."} method of the {\\normalfont \\ttfamily ".$implementationName."} implemention of the {\\normalfont \\ttfamily ".$componentClass."} component class should be output.\n";
			$functionCode .= "       !@   </description>\n";
			$functionCode .= "       !@   <type>string</type>\n";
			$functionCode .= "       !@   <cardinality>1</cardinality>\n";
			$functionCode .= "       !@ </inputParameter>\n";
			$functionCode .= "       call Get_Input_Parameter('".$parameterName."',".$parameterName.",defaultValue=.false.)\n";
			# Add a module-scope variable to store the output status of this property.
			push(
			    @{$build->{'variables'}},
			    {
				intrinsic => "logical",
				variables => [ $parameterName ]
			    }
			    );
			$outputConditions{$parameterName} = 1;
		    }
		}
	    }
    	}
    	$functionCode .= "       if (.not.allocated(default".&Utils::padClass(ucfirst($componentClass)."Component",[9,0]).")) then\n";
    	$functionCode .= "          message='unrecognized method \"'//methodSelection//'\" for \"".$componentClass."\" component'\n";
	$functionCode .= "          message=message//char(10)//'  available methods are:'\n";
    	foreach my $implementationName ( sort(@{$build->{'componentClasses'}->{$componentClass}->{'members'}}) ) {
	    $functionCode .= "          message=message//char(10)//'    ".$implementationName."'\n";
	}
    	$functionCode .= "          call Galacticus_Error_Report('Galacticus_Nodes_Initialize',message)\n";
    	$functionCode .= "       end if\n";
    }
    $build->{'content'} .= "\n";
    $functionCode .= "         ! Record that the module is now initialized.\n";
    $functionCode .= "         moduleIsInitialized=.true.\n";
    $functionCode .= "       end if\n";
    $functionCode .= "       !\$omp end critical (Galacticus_Nodes_Initialize)\n";
    $functionCode .= "    end if\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Class_Hierarchy_Initialize\n";	
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
}

1;
