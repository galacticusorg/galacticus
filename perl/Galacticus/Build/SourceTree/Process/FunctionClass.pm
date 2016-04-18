# Contains a Perl module which implements processing of functionClass directives.

package FunctionClass;
use strict;
use warnings;
use utf8;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use Data::Dumper;
use XML::Simple;
use Sort::Topological qw(toposort);
use LaTeX::Encode;
use Scalar::Util qw(reftype);
require List::ExtraUtils;
require Fortran::Utils;
require Galacticus::Build::SourceTree::Hooks;
require Galacticus::Build::SourceTree;
require Galacticus::Build::SourceTree::Parse::ModuleUses;

# Insert hooks for our functions.
$Hooks::processHooks{'functionClass'} = \&Process_FunctionClass;

sub Process_FunctionClass {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Initialize code directive locations.
    my $directiveLocations;
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "functionClass" ) {
	    # Assert that our parent is a module (for now).
	    die("Process_FunctionClass: parent node must be a module")
		unless ( $node->{'parent'}->{'type'} eq "module" );
	    # Extract the directive.
	    my $directive = $node->{'directive'};
	    # Get code directive locations if we do not have them.
	    $directiveLocations = $xml->XMLin($galacticusPath.$ENV{'BUILDPATH'}."/Code_Directive_Locations.xml")
		unless ( $directiveLocations );	    
	    # Find methods.
	    my %methods;
	    if ( exists($directive->{'method'}) ) {
		if ( exists($directive->{'method'}->{'name'}) && ! reftype($directive->{'method'}->{'name'}) ) {
		    $methods{$directive->{'method'}->{'name'}} = $directive->{'method'};
		} else {
		    %methods = %{$directive->{'method'}};
		}
	    }
	    # Find class locations.
	    my @classLocations = &ExtraUtils::as_array($directiveLocations->{$directive->{'name'}}->{'file'})
		if ( exists($directiveLocations->{$directive->{'name'}}) );
	    # Parse classes.
	    my %dependencies;
	    my %classes;
	    foreach my $classLocation ( @classLocations ) {
		my $classTree  = &SourceTree::ParseFile($classLocation);
		&SourceTree::ProcessTree($classTree, errorTolerant => 1);
		my $classNode  = $classTree;
		my $classDepth = 0;
		my $className;
		my $class;
		$class->{'file'} = $classLocation;
		while ( $classNode ) {
		    # Collect class directives.
		    if ( $classNode->{'type'} eq $directive->{'name'} ) {
			$class->{$_} = $classNode->{'directive'}->{$_}
			    foreach ( keys(%{$classNode->{'directive'}}) );
		    }
		    if ( $classNode->{'type'} eq "type" ) {
			# Parse class openers to find dependencies.
			if ( $classNode->{'opener'} =~ m/^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*$directive->{'name'}([a-z0-9_]+)\s*$/i ) {
			    my $extends = $2;
			    my $type    = $directive->{'name'}.$4;
			    $class->{'type'   } = $type;
			    $class->{'extends'} = $extends;
			    push(@{$dependencies{$extends}},$type);
			}
		    }
		    $classNode = &SourceTree::Walk_Tree($classNode,\$classDepth);
		}
		# Store tree.
		$class->{'tree'} = $classTree;
		# Set defaults.
		$class->{'abstract'            } = "no"
		    unless ( exists($class->{'abstract'}            ) );
		$class->{'defaultThreadPrivate'} = "default"
		    unless ( exists($class->{'defaultThreadPrivate'}) );
		# Store to set of all classes.
		$classes{$class->{'type'}} = $class;
	    }
	    my @unsortedClasses = keys(%classes);
	    my @sortedClasses   = toposort(sub { @{$dependencies{$_[0]} || []}; }, \@unsortedClasses );
	    my @classes = map($classes{$_},@sortedClasses);
	    # Create a set of non-abstract classes.
	    my @nonAbstractClasses;
	    foreach ( @classes ) {
		push(
		    @nonAbstractClasses,
		    $_
		    )
		    unless ( exists($_->{'abstract'}) && $_->{'abstract'} eq "yes" );
	    }
	    # If the function is stateful, add methods to store and retrieve state.
	    if ( exists($directive->{'stateful'}) && $directive->{'stateful'} eq "yes" ) {
		$methods{'stateStore'} =
		{
		    description => "Store the state of the object to file.",
		    type        => "void",
		    pass        => "yes",
		    modules     => "FGSL",
		    argument    => [ "integer, intent(in   ) :: stateFile", "type(fgsl_file), intent(in   ) :: fgslStateFile" ],
		    # <workaround type="gfortran" PR="41209" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=41209"/>
		    code        => join("",map {"if (sizeof(".$_.")<0.and.sizeof(".$_.")>0) then\nend if\n"} ('self', 'stateFile', 'fgslStateFile') )
		};
 		$methods{'stateRestore'} =
		{
		    description => "Restore the state of the object to file.",
		    type        => "void",
		    pass        => "yes",
		    modules     => "FGSL",
		    argument    => [ "integer, intent(in   ) :: stateFile", "type(fgsl_file), intent(in   ) :: fgslStateFile" ],
		    code        => join("",map {"if (sizeof(".$_.")<0.and.sizeof(".$_.")>0) then\nend if\n"} ('self', 'stateFile', 'fgslStateFile') )
		};
  		$methods{'stateSnapshot'} =
		{
		    description => "Snapshot the state of the object.",
		    type        => "void",
		    pass        => "yes",
		    argument    => [ ],
		    code        => join("",map {"if (sizeof(".$_.")<0.and.sizeof(".$_.")>0) then\nend if\n"} ('self') )
		};
	    }
	    # If the function requires calculation reset, add method to do so.
	    if ( exists($directive->{'calculationReset'}) && $directive->{'calculationReset'} eq "yes" ) {
		$methods{'calculationReset'} =
		{
		    description => "Reset the calculation state of the object.",
		    type        => "void",
		    pass        => "yes",
		    argument    => [ "type(treeNode), intent(inout), pointer :: thisNode" ],
		    code        => join("",map {"if (sizeof(".$_.")<0.and.sizeof(".$_.")>0) then\nend if\n"} ('self','thisNode') )
		};
	    }
	    # Add "isFinalizable" method.
	    $methods{'isFinalizable'} = 
	    {
		description => "Return true if this object can be finalized.",
		type        => "logical",
		pass        => "yes",
		code        => $directive->{'name'}."isFinalizable=.not.self%isDefault\n"
	    };
	    # Add "descriptor" method.
	    my $descriptorCode;
	    $descriptorCode .= "select type (self)\n";
	    foreach ( @nonAbstractClasses ) {
		(my $label = $_->{'name'}) =~ s/^$directive->{'name'}//;
		$label = lcfirst($label)
		    unless ( $label =~ m/^[A-Z]{2,}/ );
		$descriptorCode .= "type is (".$_->{'name'}.")\n";
		$descriptorCode .= " call descriptor%addParameter('".$directive->{'name'}."Method','".$label."')\n";	
	    }
	    $descriptorCode .= "end select\n";
	    $methods{'descriptor'} = 
	    {
		description => "Return an input parameter list descriptor which could be used to recreate this object.",
		type        => "void",
		pass        => "yes",
		modules     => "Input_Parameters2",
		argument    => [ "type(inputParameters), intent(inout) :: descriptor" ],
		code        => $descriptorCode
	    };
	    # Determine if any methods request that C-bindings be produced.
	    my %methodsCBound;
	    foreach ( keys(%methods) ) {
		$methodsCBound{$_} = $methods{$_}
		    if ( exists($methods{$_}->{'bindC'}) && $methods{$_}->{'bindC'} eq "true" );
	    }
	    # Initialize new nodes.
	    my $preContains = 
		[
		 {
		     type       => "code" ,
		     content    => ""     ,
		     firstChild => undef(),
		     sibling    => undef(),
		     parent     => undef(),
		     source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
		     line       => 1
		 }
		];
	    my $postContains =
		[
		 {
		     type       => "code" ,
		     content    => ""     ,
		     firstChild => undef(),
		     sibling    => undef(),
		     parent     => undef(),
		     source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
		     line       => 1
		 }
		];

	    # Add variable tracking module initialization status.
	    $preContains->[0]->{'content'} .= "   logical, private :: moduleInitialized=.false.\n\n";

	    # Generate the base class.
	    &SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}."Class","public");
	    &SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}        ,"public");
	    $preContains->[0]->{'content'} .= "   type :: ".$directive->{'name'}."Class\n";
	    $preContains->[0]->{'content'} .= "    private\n";
	    $preContains->[0]->{'content'} .= "    logical :: isDefault=.false.\n";
	    foreach ( &ExtraUtils::as_array($directive->{'data'}) ) {
		if ( reftype($_) ) {
		    $_->{'scope'} = "self"
			unless ( exists($_->{'scope'}) );
		    $preContains->[0]->{'content'} .= $_->{'content'}."\n"
			if (  $_->{'scope'} eq "self" );
		} else {
		    $preContains->[0]->{'content'} .= $_."\n";
		}
	    }
	    $preContains->[0]->{'content'} .= "    contains\n";
	    $preContains->[0]->{'content'} .= "    !@ <objectMethods>\n";
	    $preContains->[0]->{'content'} .= "    !@   <object>".$directive->{'name'}."Class</object>\n";
	    foreach my $methodName ( keys(%methods) ) {
		my $method = $methods{$methodName};
		my $argumentList = "";
		my @arguments;
		if ( exists($method->{'argument'}) ) {
		    if ( UNIVERSAL::isa($method->{'argument'},"ARRAY") ) {
			push(@arguments,@{$method->{'argument'}});
		    } else {
			push(@arguments,  $method->{'argument'} );
		    }
		}
		my $separator = "";
		foreach my $argument ( @arguments ) {
		    foreach my $intrinsic ( keys(%Fortran_Utils::intrinsicDeclarations) ) {
			my $declarator = $Fortran_Utils::intrinsicDeclarations{$intrinsic};
			if ( my @matches = $argument =~ m/$declarator->{'regEx'}/ ) {
			    my $intrinsicName =                          $declarator->{'intrinsic' }  ;
			    my $type          =                 $matches[$declarator->{'type'      }] ;
			    my $attributeList =                 $matches[$declarator->{'attributes'}] ;
			    $attributeList =~ s/^\s*,?\s*//;
			    $attributeList =~ s/\s*$//;
			    my @attributes = &Fortran_Utils::Extract_Variables($attributeList, keepQualifiers => 1, removeSpaces => 1);
			    my @variables     = split(/\s*,\s*/,$matches[$declarator->{'variables' }]);
			    foreach my $variable ( @variables ) {
				$argumentList .= $separator."\\textcolor{red}{\\textless ".$intrinsicName;
				$argumentList .= "(".$type.")"
				    if ( defined($type) );
				$argumentList .= "\\textgreater} ".$variable;
				foreach my $attribute ( @attributes ) {
				    $argumentList .= "\\argin"
					if ( $attribute eq "intent(in)" );
				    $argumentList .= "\\argout"
					if ( $attribute eq "intent(out)" );
				    $argumentList .= "\\arginout"
					if ( $attribute eq "intent(inout)" );
				}
				$separator     = ",";
			    }
			}
		    }
		}
		$preContains->[0]->{'content'} .= "    !@   <objectMethod>\n";
		$preContains->[0]->{'content'} .= "    !@     <method>".$methodName."</method>\n";
		$preContains->[0]->{'content'} .= "    !@     <type>".$method->{'type'}."</type>\n";
		$preContains->[0]->{'content'} .= "    !@     <arguments>".latex_encode($argumentList)."</arguments>\n";
		$preContains->[0]->{'content'} .= "    !@     <description>".$method->{'description'}."</description>\n";
		$preContains->[0]->{'content'} .= "    !@   </objectMethod>\n";
	    }
	    $preContains->[0]->{'content'} .= "    !@ </objectMethods>\n";
	    my $methodTable = Text::Table->new(
		{
		    is_sep => 1,
		    body   => "    procedure"
		},
		{
		    align  => "left"
		},
		{
		    is_sep => 1,
		    body   => " :: "
		},
		{
		    align  => "left"
		},
		{
		    is_sep => 1,
		    body   => " => ",
		},
		{
		    align  => "left"
		}
		);    
	    foreach ( keys(%methods) ) {
		my $method = $methods{$_};
		my $extension = "Null";
		$extension = ""
		    if ( exists($method->{'code'}) );
		$methodTable->add("",$_,$directive->{'name'}.ucfirst($_).$extension);
	    }
	    $preContains->[0]->{'content'} .= $methodTable->table();
	    $preContains->[0]->{'content'} .= "   end type ".$directive->{'name'}."Class\n\n";
	    # Insert any module-scope class content.
	    foreach ( &ExtraUtils::as_array($directive->{'data'}) ) {
		if ( reftype($_) ) {
		    if ( exists($_->{'scope'}) && $_->{'scope'} eq "module" ) {
			$preContains->[0]->{'content'} .= $_->{'content'}."\n";
			if ( exists($_->{'threadprivate'}) && $_->{'threadprivate'} eq "yes" && $_->{'content'} =~ m/::\s*(.*)$/ ) {
			    $preContains->[0]->{'content'} .= "   !\$omp threadprivate(".$1.")\n";
			}
		    }
		}
	    }

	    # Generate class constructors
	    $preContains->[0]->{'content'} .= "   interface ".$directive->{'name'}."\n";
	    $preContains->[0]->{'content'} .= "    module procedure ".$directive->{'name'}."ConstructorDefault\n";
	    $preContains->[0]->{'content'} .= "    module procedure ".$directive->{'name'}."ConstructorParameters\n";
	    $preContains->[0]->{'content'} .= "   end interface ".$directive->{'name'}."\n";
	    # Add method name parameter.
	    $preContains->[0]->{'content'} .= "   ! Method name parameter.\n";
	    $preContains->[0]->{'content'} .= "   type(varying_string) :: ".$directive->{'name'}."Method\n\n";
	    my $usesNode =
	    {
		type      => "moduleUse",
		moduleUse =>
		{
		    ISO_Varying_String =>
		    {
			intrinsic => 0,
			all       => 1
		    }
		}
	    };
	    &ModuleUses::AddUses($node->{'parent'},$usesNode);
	    if ( $tree->{'type'} eq "file" ) {
		(my $fileName = $tree->{'name'}) =~ s/\.F90$/.p/;
		open(my $parametersFile,">>".$ENV{'BUILDPATH'}."/".$fileName);
		print $parametersFile $directive->{'name'}."Method\n";
		close($parametersFile);
	    }
	    # Add default implementation.
	    $directive->{'defaultThreadPrivate'} = "no"
		unless ( exists($directive->{'defaultThreadPrivate'}) );
	    my $requireThreadPublicDefault  = 0;    
	    foreach my $className ( keys(%classes) ) {
		my $class = $classes{$className};
		$class->{'defaultThreadPrivate'} = $directive->{'defaultThreadPrivate'}
		    if ( $class->{'defaultThreadPrivate'} eq "default" );
		$requireThreadPublicDefault  = 1
		    if ( $class->{'defaultThreadPrivate'} eq "no"      );
	    }
	    $preContains->[0]->{'content'} .= "   ! Default ".$directive->{'name'}." object.\n";
	    $preContains->[0]->{'content'} .= "   class(".$directive->{'name'}."Class), private , pointer :: ".$directive->{'name'}."Default       => null()\n";
	    $preContains->[0]->{'content'} .= "   !\$omp threadprivate(".$directive->{'name'}."Default)\n";
	    $preContains->[0]->{'content'} .= "   class(".$directive->{'name'}."Class), private , pointer :: ".$directive->{'name'}."PublicDefault => null()\n"
		if ( $requireThreadPublicDefault  == 1 );
	    $preContains->[0]->{'content'} .= "\n";
	    # Create default constructor.
	    $postContains->[0]->{'content'} .= "   function ".$directive->{'name'}."ConstructorDefault()\n";
	    $postContains->[0]->{'content'} .= "      !% Return a pointer to the default {\\normalfont \\ttfamily ".$directive->{'name'}."} object.\n";
	    $postContains->[0]->{'content'} .= "      implicit none\n";
	    $postContains->[0]->{'content'} .= "      class(".$directive->{'name'}."Class), pointer :: ".$directive->{'name'}."ConstructorDefault\n\n";
	    $postContains->[0]->{'content'} .= "      if (.not.associated(".$directive->{'name'}."Default)) call ".$directive->{'name'}."Initialize()\n";
	    $postContains->[0]->{'content'} .= "      ".$directive->{'name'}."ConstructorDefault => ".$directive->{'name'}."Default\n";
	    $postContains->[0]->{'content'} .= "      return\n";
	    $postContains->[0]->{'content'} .= "   end function ".$directive->{'name'}."ConstructorDefault\n\n";
	    # Create XML constructor.
	    $postContains->[0]->{'content'} .= "   function ".$directive->{'name'}."ConstructorParameters(parameters,copyInstance)\n";
	    $postContains->[0]->{'content'} .= "      !% Return a pointer to a newly created {\\normalfont \\ttfamily ".$directive->{'name'}."} object as specified by the provided parameters.\n";
	    $postContains->[0]->{'content'} .= "      use Input_Parameters2\n";
	    $postContains->[0]->{'content'} .= "      use Galacticus_Error\n";
	    $postContains->[0]->{'content'} .= "      implicit none\n";
	    $postContains->[0]->{'content'} .= "      class  (".$directive->{'name'}."Class), pointer :: ".$directive->{'name'}."ConstructorParameters\n";
	    $postContains->[0]->{'content'} .= "      type   (inputParameters), intent(inout)           :: parameters\n";
	    $postContains->[0]->{'content'} .= "      integer                 , intent(in   ), optional :: copyInstance\n";
	    $postContains->[0]->{'content'} .= "      type   (inputParameters)                          :: subParameters\n";
	    $postContains->[0]->{'content'} .= "      type   (varying_string )                          :: message      , instanceName\n\n";
	    $postContains->[0]->{'content'} .= "      call parameters%value('".$directive->{'name'}."Method',instanceName,copyInstance=copyInstance)\n";
	    $postContains->[0]->{'content'} .= "      subParameters=parameters%subParameters('".$directive->{'name'}."Method',copyInstance=copyInstance)\n";
	    $postContains->[0]->{'content'} .= "      select case (char(instanceName))\n";
	    foreach my $class ( @nonAbstractClasses ) {
		(my $name = $class->{'name'}) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$postContains->[0]->{'content'} .= "     case ('".$name."')\n";
		$postContains->[0]->{'content'} .= "        allocate(".$class->{'name'}." :: ".$directive->{'name'}."ConstructorParameters)\n";
		$postContains->[0]->{'content'} .= "        select type (".$directive->{'name'}."ConstructorParameters)\n";
		$postContains->[0]->{'content'} .= "          type is (".$class->{'name'}.")\n";
		$postContains->[0]->{'content'} .= "            ".$directive->{'name'}."ConstructorParameters=".$class->{'name'}."(subParameters)\n";
		$postContains->[0]->{'content'} .= "         end select\n";
	    }
	    $postContains->[0]->{'content'} .= "      case default\n";
	    $postContains->[0]->{'content'} .= "         message='Unrecognized type \"'//trim(instanceName)//'\" Available options are:'\n";
	    my @classNames;
	    push(@classNames,$_->{'name'})
		foreach ( @nonAbstractClasses );
	    foreach ( sort(@classNames) ) {
		(my $name = $_) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$postContains->[0]->{'content'} .= "         message=message//char(10)//'   -> ".$name."'\n";
	    }
	    $postContains->[0]->{'content'} .= "         call Galacticus_Error_Report('".$directive->{'name'}."ConstructorParameters',message)\n";
	    $postContains->[0]->{'content'} .= "      end select\n";
	    $postContains->[0]->{'content'} .= "      return\n";
	    $postContains->[0]->{'content'} .= "   end function ".$directive->{'name'}."ConstructorParameters\n\n";
	    
	    # Insert class code.
	    foreach my $class ( @classes ) {
		&SourceTree::SetVisibility($node->{'parent'},$class->{'type'},"public")	
		    if ( grep {$_->{'type'} eq $class->{'type'}} @nonAbstractClasses );
		my $classTree = $class->{'tree'};
		my $classNode = $classTree->{'firstChild'};
		my $contained = 0;
		while ( $classNode ) {
		    # Check for composition in nonprivate instances.
		    if ( $class->{'defaultThreadPrivate'} eq "no" ) {
			my $subNode = $classNode;
			while ( $subNode ) {
			    if ( $subNode->{'type'} eq "objectBuilder" ) {
				print 
				    "WARN: instance '"                                                                                          .
				    $class    ->{'type'}                                                                                  .
				    "' of function class '"                                                                               .
				    $directive->{'name'}                                                                                  .
				    "' is not default thread private, but composites other objects which may be default thread private.\n";
			    }
			    $subNode = &SourceTree::Walk_Tree($subNode);
			}
		    }		    
		    if ( $classNode->{'type'} eq "contains" ) {
			$classNode = $classNode->{'firstChild'};
			$contained = 1;
		    }
		    if ( $contained ) {
			push(@{$postContains},$classNode);
		    } else {
			if ( $classNode->{'type'} eq "moduleUse" ) {
			    &ModuleUses::AddUses($node->{'parent'},$classNode);
			} else {			    
			    push(@{$preContains },$classNode);
			}
		    }
		    $classNode = $classNode->{'sibling'};
		}
	    }
	    
	    # Create initialization function.
	    $postContains->[0]->{'content'} .= "   subroutine ".$directive->{'name'}."Initialize()\n";
	    $postContains->[0]->{'content'} .= "      !% Initialize the default {\\normalfont \\ttfamily ".$directive->{'name'}."} object.\n";
	    $postContains->[0]->{'content'} .= "      use ISO_Varying_String\n";
	    $postContains->[0]->{'content'} .= "      use Input_Parameters2\n";
	    $postContains->[0]->{'content'} .= "      use Galacticus_Error\n";
	    $postContains->[0]->{'content'} .= "      use IO_HDF5\n";
	    $postContains->[0]->{'content'} .= "      implicit none\n";
	    $postContains->[0]->{'content'} .= "      type(inputParameters) :: subParameters\n";
	    $postContains->[0]->{'content'} .= "      type(varying_string ) :: message\n\n";
	    $postContains->[0]->{'content'} .= "      !\$omp critical (".$directive->{'name'}."Initialization)\n";
	    $postContains->[0]->{'content'} .= "      if (.not.moduleInitialized) then\n";
	    $postContains->[0]->{'content'} .= "         !@ <inputParameter>\n";
	    $postContains->[0]->{'content'} .= "         !@   <name>".$directive->{'name'}."Method</name>\n";
	    $postContains->[0]->{'content'} .= "         !@   <defaultValue>".$directive->{'default'}."</defaultValue>\n";
	    $postContains->[0]->{'content'} .= "         !@   <attachedTo>module</attachedTo>\n";
	    $postContains->[0]->{'content'} .= "         !@   <description>\n";
	    $postContains->[0]->{'content'} .= "         !@     The method to be used for {\\normalfont \\ttfamily ".$directive->{'name'}."}.\n";
	    $postContains->[0]->{'content'} .= "         !@   </description>\n";
	    $postContains->[0]->{'content'} .= "         !@   <type>string</type>\n";
	    $postContains->[0]->{'content'} .= "         !@   <cardinality>1</cardinality>\n";
	    $postContains->[0]->{'content'} .= "         !@ </inputParameter>\n";
	    $postContains->[0]->{'content'} .= "         call globalParameters%value('".$directive->{'name'}."Method',".$directive->{'name'}."Method,defaultValue=var_str('".$directive->{'default'}."'))\n";
	    $postContains->[0]->{'content'} .= "         moduleInitialized=.true.\n";
	    $postContains->[0]->{'content'} .= "      end if\n";
	    $postContains->[0]->{'content'} .= "      subParameters=globalParameters%subParameters('".$directive->{'name'}."Method',requirePresent=.false.)\n";
	    $postContains->[0]->{'content'} .= "      select case (char(".$directive->{'name'}."Method))\n";
	    foreach my $class ( @nonAbstractClasses ) {
		(my $name = $class->{'name'}) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$postContains->[0]->{'content'} .= "     case ('".$name."')\n";
		if ( $class->{'defaultThreadPrivate'} eq "yes" ) {
		    $postContains->[0]->{'content'} .= "        parametersObjectBuildIsPrivate=.true.\n";
		    $postContains->[0]->{'content'} .= "        allocate(".$class->{'name'}." :: ".$directive->{'name'}."Default)\n";
		    $postContains->[0]->{'content'} .= "        select type (".$directive->{'name'}."Default)\n";
		    $postContains->[0]->{'content'} .= "          type is (".$class->{'name'}.")\n";
		    $postContains->[0]->{'content'} .= "            ".$directive->{'name'}."Default=".$class->{'name'}."(subParameters)\n";
		    $postContains->[0]->{'content'} .= "         end select\n";
		} else {
		    $postContains->[0]->{'content'} .= "        parametersObjectBuildIsPrivate=.false.\n";
		    $postContains->[0]->{'content'} .= "        if (.not.associated(".$directive->{'name'}."PublicDefault)) then\n";
		    $postContains->[0]->{'content'} .= "           allocate(".$class->{'name'}." :: ".$directive->{'name'}."PublicDefault)\n";
		    $postContains->[0]->{'content'} .= "           select type (".$directive->{'name'}."PublicDefault)\n";
		    $postContains->[0]->{'content'} .= "           type is (".$class->{'name'}.")\n";
		    $postContains->[0]->{'content'} .= "             ".$directive->{'name'}."PublicDefault=".$class->{'name'}."(subParameters)\n";
		    $postContains->[0]->{'content'} .= "           end select\n";
		    $postContains->[0]->{'content'} .= "        end if\n";
		    $postContains->[0]->{'content'} .= "         ".$directive->{'name'}."Default => ".$directive->{'name'}."PublicDefault\n";
		}
	    }
	    $postContains->[0]->{'content'} .= "      case default\n";
	    $postContains->[0]->{'content'} .= "         message='Unrecognized option for [".$directive->{'name'}."Method](='//".$directive->{'name'}."Method//'). Available options are:'\n";
	    foreach ( sort(@classNames) ) {
		(my $name = $_) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$postContains->[0]->{'content'} .= "        message=message//char(10)//'   -> ".$name."'\n";
	    }
	    $postContains->[0]->{'content'} .= "         call Galacticus_Error_Report('".$directive->{'name'}."Initialize',message)\n";
	    $postContains->[0]->{'content'} .= "      end select\n";
	    $postContains->[0]->{'content'} .= "      ".$directive->{'name'}."Default%isDefault=.true.\n";
	    $postContains->[0]->{'content'} .= "      !\$omp end critical (".$directive->{'name'}."Initialization)\n";
	    $postContains->[0]->{'content'} .= "      return\n";
	    $postContains->[0]->{'content'} .= "   end subroutine ".$directive->{'name'}."Initialize\n\n";

	    # Create global state store/restore functions.
	    if ( exists($directive->{'stateful'}) && $directive->{'stateful'} eq "yes" ) {
		&SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}.$_,"public")
		    foreach ( "DoStateStore", "DoStateRetrieve", "DoStateSnapshot" );
		$postContains->[0]->{'content'} .= "  !# <galacticusStateStoreTask>\n";
		$postContains->[0]->{'content'} .= "  !#  <unitName>".$directive->{'name'}."DoStateStore</unitName>\n";
		$postContains->[0]->{'content'} .= "  !# </galacticusStateStoreTask>\n";
		$postContains->[0]->{'content'} .= "  subroutine ".$directive->{'name'}."DoStateStore(stateFile,fgslStateFile)\n";
		$postContains->[0]->{'content'} .= "    !% Store the state to file.\n";
		$postContains->[0]->{'content'} .= "    use FGSL\n";
		$postContains->[0]->{'content'} .= "    implicit none\n";
		$postContains->[0]->{'content'} .= "    integer           , intent(in   ) :: stateFile\n";
		$postContains->[0]->{'content'} .= "    type   (fgsl_file), intent(in   ) :: fgslStateFile\n";
		$postContains->[0]->{'content'} .= "    class  (".$directive->{'name'}."Class), pointer :: default\n\n";
		$postContains->[0]->{'content'} .= "    default => ".$directive->{'name'}."()\n";
		$postContains->[0]->{'content'} .= "    call default%stateStore(stateFile,fgslStateFile)\n";
		$postContains->[0]->{'content'} .= "    return\n";
		$postContains->[0]->{'content'} .= "  end subroutine ".$directive->{'name'}."DoStateStore\n\n";
		$postContains->[0]->{'content'} .= "  !# <galacticusStateRetrieveTask>\n";
		$postContains->[0]->{'content'} .= "  !#  <unitName>".$directive->{'name'}."DoStateRetrieve</unitName>\n";
		$postContains->[0]->{'content'} .= "  !# </galacticusStateRetrieveTask>\n";
		$postContains->[0]->{'content'} .= "  subroutine ".$directive->{'name'}."DoStateRetrieve(stateFile,fgslStateFile)\n";
		$postContains->[0]->{'content'} .= "    !% Retrieve the state from file.\n";
		$postContains->[0]->{'content'} .= "    use FGSL\n";
		$postContains->[0]->{'content'} .= "    implicit none\n";
		$postContains->[0]->{'content'} .= "    integer           , intent(in   ) :: stateFile\n";
		$postContains->[0]->{'content'} .= "    type   (fgsl_file), intent(in   ) :: fgslStateFile\n";
		$postContains->[0]->{'content'} .= "    class  (".$directive->{'name'}."Class), pointer :: default\n\n";
		$postContains->[0]->{'content'} .= "    default => ".$directive->{'name'}."()\n";
		$postContains->[0]->{'content'} .= "    call default%stateRestore(stateFile,fgslStateFile)\n";
		$postContains->[0]->{'content'} .= "    return\n";
		$postContains->[0]->{'content'} .= "  end subroutine ".$directive->{'name'}."DoStateRetrieve\n\n";
		$postContains->[0]->{'content'} .= "  !# <galacticusStateSnapshotTask>\n";
		$postContains->[0]->{'content'} .= "  !#  <unitName>".$directive->{'name'}."DoStateSnapshot</unitName>\n";
		$postContains->[0]->{'content'} .= "  !# </galacticusStateSnapshotTask>\n";
		$postContains->[0]->{'content'} .= "  subroutine ".$directive->{'name'}."DoStateSnapshot()\n";
		$postContains->[0]->{'content'} .= "    !% Snapshot the object.\n";
		$postContains->[0]->{'content'} .= "    implicit none\n";
		$postContains->[0]->{'content'} .= "    class  (".$directive->{'name'}."Class), pointer :: default\n\n";
		$postContains->[0]->{'content'} .= "    default => ".$directive->{'name'}."()\n";
		$postContains->[0]->{'content'} .= "    call default%stateSnapshot()\n";
		$postContains->[0]->{'content'} .= "    return\n";
		$postContains->[0]->{'content'} .= "  end subroutine ".$directive->{'name'}."DoStateSnapshot\n\n";
	    }
	    
	    # Create global calculation reset function.
	    if ( exists($directive->{'calculationReset'}) && $directive->{'calculationReset'} eq "yes" ) {
		&SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}."DoCalculationReset","public");
		$postContains->[0]->{'content'} .= "  !# <calculationResetTask>\n";
		$postContains->[0]->{'content'} .= "  !#  <unitName>".$directive->{'name'}."DoCalculationReset</unitName>\n";
		$postContains->[0]->{'content'} .= "  !# </calculationResetTask>\n";
		$postContains->[0]->{'content'} .= "  subroutine ".$directive->{'name'}."DoCalculationReset(thisNode)\n";
		$postContains->[0]->{'content'} .= "    !% Store the state to file.\n";
		$postContains->[0]->{'content'} .= "    implicit none\n";
		$postContains->[0]->{'content'} .= "    type (treeNode), pointer, intent(inout) :: thisNode\n";
		$postContains->[0]->{'content'} .= "    class(".$directive->{'name'}."Class), pointer :: default\n\n";
		$postContains->[0]->{'content'} .= "    default => ".$directive->{'name'}."()\n";
		$postContains->[0]->{'content'} .= "    call default%calculationReset(thisNode)\n";
		$postContains->[0]->{'content'} .= "    return\n";
		$postContains->[0]->{'content'} .= "  end subroutine ".$directive->{'name'}."DoCalculationReset\n\n";
	    }

	    # Create functions.
	    foreach my $methodName ( keys(%methods) ) {
		my $method = $methods{$methodName};
		# Insert arguments.
		my @arguments;
		if ( exists($method->{'argument'}) ) {
		    if ( UNIVERSAL::isa($method->{'argument'},"ARRAY") ) {
			push(@arguments,@{$method->{'argument'}});
		    } else {
			push(@arguments,  $method->{'argument'} );
		    }
		}
		my $pass = "yes";
		$pass = $method->{'pass'}
		    if ( exists($method->{'pass'}) );
		my $argumentList = "";
		my $argumentCode;
		if ( $pass eq "yes" ) {
		    $argumentCode .= "      class(".$directive->{'name'}."Class), intent(inout)";
		    $argumentCode .= ", target"
			if ( exists($method->{'selfTarget'}) && $method->{'selfTarget'} eq "yes" );
		    $argumentCode .= " :: self\n";
		}
		my $separator = "";
		foreach my $argument ( @arguments ) {
		    (my $variables = $argument) =~ s/^.*::\s*(.*?)\s*$/$1/;
		    $argumentList .= $separator.$variables;
		    $argumentCode .= "      ".$argument."\n";
		    $separator     = ",";
		}
		my $type;
		my $category;
		my $self;
		my $extension = "Null";
		$extension = ""
		    if ( exists($method->{'code'}) );
		if ( $method->{'type'} eq "void" ) {
		    $category = "subroutine";
		    $type     = "";
		    $self     = "";
		} elsif ( $method->{'type'} =~ m/^class/ ) {
		    $category = "function";
		    $type     = "";
		    $self     = "      ".$method->{'type'}.", pointer :: ".$directive->{'name'}.ucfirst($methodName).$extension."\n";
		} elsif ( $method->{'type'} =~ m/^type/ ) {
		    $category = "function";
		    $type     = "";
		    $self     = "      ".$method->{'type'}." :: ".$directive->{'name'}.ucfirst($methodName).$extension."\n";
		} else {
		    $category = "function";
		    $type     = $method->{'type'}." ";
		    $self     = "";
		}
		$postContains->[0]->{'content'} .= "   ".$type.$category." ".$directive->{'name'}.ucfirst($methodName).$extension."(self";
		$postContains->[0]->{'content'} .= ",".$argumentList
		    unless ( $argumentList eq "" );
		$postContains->[0]->{'content'} .= ")\n";
		$postContains->[0]->{'content'} .= "      !% ".$method->{'description'}."\n";
		if ( exists($method->{'code'}) ) {
		    if ( exists($method->{'modules'}) ) {
   		$postContains->[0]->{'content'} .= "      use ".$_."\n"
   		    foreach ( split(/\s+/,$method->{'modules'}) );
		    }
		} else {
		    $postContains->[0]->{'content'} .= "      use Galacticus_Error\n";
		}
		$postContains->[0]->{'content'} .= "      implicit none\n";
		$postContains->[0]->{'content'} .= $self;
		$postContains->[0]->{'content'} .= $argumentCode;
		if ( exists($method->{'code'}) ) {
		    my $code = "      ".$method->{'code'};
		    $code =~ s/\n/\n      /g;
		    $postContains->[0]->{'content'} .= $code."\n";
		} else {
		    $postContains->[0]->{'content'} .= "      call Galacticus_Error_Report('".$methodName."Null','this is a null method - initialize the ".$directive->{'name'}." object before use')\n";
		    if ( $category eq "function" ) {
			# Avoid warnings about unset function values.
			$postContains->[0]->{'content'} .= "      ".$directive->{'name'}.ucfirst($methodName).$extension."=";
			my $setValue;
			if ( $method->{'type'} =~ m/^class/ ) {
			    $setValue = "> null()";
			} elsif ( $method->{'type'} =~ m/^type\s*\(\s*(.*)\s*\)/ ) {
			    $setValue = $2."()";
			} elsif ( $method->{'type'} =~ m/^integer/ ) {
			    $setValue = "0";
			} elsif ( $method->{'type'} =~ m/^double\s+precision/ ) {
			    $setValue = "0.0d0";
			}
			die("Process_FunctionClass(): do not know how to set '".$method->{'type'}."'")
			    unless ( defined($setValue) );
			$postContains->[0]->{'content'} .= $setValue."\n";
		    }
		    $postContains->[0]->{'content'} .= "      return\n";
		    # <workaround type="gfortran" PR="41209" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=41209"/>
		    $postContains->[0]->{'content'} .= join("",map {"if (sizeof(".$_.")<0.and.sizeof(".$_.")>0) then\nend if\n"} split(/,/,$argumentList eq "" ? "self" : "self,".$argumentList));
		}
		$postContains->[0]->{'content'} .= "   end ".$category." ".$directive->{'name'}.ucfirst($methodName).$extension."\n\n";
	    }
 	    
	    # Generate C-bindings if required.
	    if ( %methodsCBound ) {
		# Insert a wrapper class to permit passing of polymorphic pointers between Fortran and C++.
		$preContains->[0]->{'content'} .= "   type :: ".$directive->{'name'}."Wrapper\n";
		$preContains->[0]->{'content'} .= "     class(".$directive->{'name'}."Class), pointer :: wrappedObject\n";
		$preContains->[0]->{'content'} .= "   end type ".$directive->{'name'}."Wrapper\n\n";
		# C-bound default constructor. Here, we use a wrapper object which contains a pointer to the default polymorphic Fortran
		# object. This wrapper is then passed back to the calling C++ function so that it can be stored in the appropriate C++
		# class.
		$postContains->[0]->{'content'} .= "   function ".$directive->{'name'}."_C() bind(c,name='".$directive->{'name'}."')\n";
		$postContains->[0]->{'content'} .= "     use, intrinsic :: ISO_C_Binding\n";
		$postContains->[0]->{'content'} .= "     implicit none\n";
		$postContains->[0]->{'content'} .= "     type(c_ptr) :: ".$directive->{'name'}."_C\n";
		$postContains->[0]->{'content'} .= "     type(".$directive->{'name'}."Wrapper), pointer :: wrapper\n";
		$postContains->[0]->{'content'} .= "      if (.not.associated(".$directive->{'name'}."Default)) call ".$directive->{'name'}."Initialize()\n";
		$postContains->[0]->{'content'} .= "       allocate(wrapper)\n";
		$postContains->[0]->{'content'} .= "       wrapper%wrappedObject => ".$directive->{'name'}."Default\n";
		$postContains->[0]->{'content'} .= "       ".$directive->{'name'}."_C=c_loc(wrapper)\n";
		$postContains->[0]->{'content'} .= "     return\n";
		$postContains->[0]->{'content'} .= "   end function ".$directive->{'name'}."_C\n\n";
		# C-bound destructor. We simply deallocate the wrapper object, letting the associated finalizor clean up the Fortran
		# object.
		$postContains->[0]->{'content'} .= "   subroutine ".$directive->{'name'}."Destructor_C(wrapperC) bind(c,name='".$directive->{'name'}."Destructor')\n";
		$postContains->[0]->{'content'} .= "     use, intrinsic :: ISO_C_Binding\n";
		$postContains->[0]->{'content'} .= "     implicit none\n";
		$postContains->[0]->{'content'} .= "     type(c_ptr), intent(in   ), value :: wrapperC\n";
		$postContains->[0]->{'content'} .= "     type(".$directive->{'name'}."Wrapper), pointer :: wrapper\n\n";
		$postContains->[0]->{'content'} .= "     call c_f_pointer(wrapperC,wrapper)\n";
		$postContains->[0]->{'content'} .= "     deallocate(wrapper)\n";
		$postContains->[0]->{'content'} .= "     return\n";
		$postContains->[0]->{'content'} .= "   end subroutine ".$directive->{'name'}."Destructor_C\n\n";
		# Generate method functions.
		foreach my $methodName ( keys(%methodsCBound) ) {
		    my @arguments;
		    if ( exists($methodsCBound{$methodName}->{'argument'}) ) {
			if ( UNIVERSAL::isa($methodsCBound{$methodName}->{'argument'},"ARRAY") ) {
			    push(@arguments,@{$methodsCBound{$methodName}->{'argument'}});
			} else {
			    push(@arguments,  $methodsCBound{$methodName}->{'argument'} );
			}
		    }
		    my $separator    = "";
		    my $argumentList = "";
		    foreach my $argument ( @arguments ) {
			foreach my $intrinsic ( keys(%Fortran_Utils::intrinsicDeclarations) ) {
			    my $declarator = $Fortran_Utils::intrinsicDeclarations{$intrinsic};
			    if ( my @matches = $argument =~ m/$declarator->{'regEx'}/ ) {
				my $intrinsicName =                          $declarator->{'intrinsic' }  ;
				my $type          =                 $matches[$declarator->{'type'      }] ;
				my $attributeList =                 $matches[$declarator->{'attributes'}] ;
				$attributeList =~ s/^\s*,?\s*//;
				$attributeList =~ s/\s*$//;
				my @attributes = &Fortran_Utils::Extract_Variables($attributeList, keepQualifiers => 1, removeSpaces => 1);
				foreach my $attribute ( @attributes ) {
				    die("Galacticus::Build::Functions::Functions_Generate_Output:  attribute not supported for C++-binding")
					unless ( $attribute eq "intent(in)" );
				}
				my @variables     = split(/\s*,\s*/,$matches[$declarator->{'variables' }]);
				die("Galacticus::Build::Functions::Functions_Generate_Output: non-standard kinds are not supported for C++-binding")
				    if ( defined($type) );
				$argumentList .= $separator.join(",",@variables);
				$separator     = ",";
			    }
			}
		    }
		    $postContains->[0]->{'content'} .= "  real(c_double) function ".$methodName."_C(";
		    $postContains->[0]->{'content'} .= "wrapperC".$separator
			if ( $methodsCBound{$methodName}->{'pass'} eq "yes" );
		    $postContains->[0]->{'content'} .= $argumentList.") bind(c,name='".$methodName."_C')\n";
		    $postContains->[0]->{'content'} .= "     use, intrinsic :: ISO_C_Binding\n";
		    $postContains->[0]->{'content'} .= "     implicit none\n";
		    $postContains->[0]->{'content'} .= "     type(c_ptr), intent(in   ), value :: wrapperC\n";
		    foreach my $argument( @arguments ) {
			(my $argumentInteroperable = $argument) =~ s/(\s*::)/, value$1/;
			$argumentInteroperable =~ s/^\s*double\s+precision/real(c_double)/;
			$postContains->[0]->{'content'} .= "     ".$argumentInteroperable."\n"
		    }
		    $postContains->[0]->{'content'} .= "     type(".$directive->{'name'}."Wrapper), pointer :: wrapper\n";
		    $postContains->[0]->{'content'} .= "     call c_f_pointer(wrapperC,wrapper)\n";
		    $postContains->[0]->{'content'} .= "     ".$methodName."_C=wrapper\%wrappedObject\%".$methodName."(".$argumentList.")\n";
		    $postContains->[0]->{'content'} .= "     return\n";
		    $postContains->[0]->{'content'} .= "   end function ".$methodName."_C\n\n";
		}
		# Iterate over methods and generate the necessary code.
		my $externCode;
		my $classCode;
		my $methodCode;
		foreach my $methodName ( keys(%methodsCBound) ) {
		    my $type;
		    if ( $methodsCBound{$methodName}->{'type'} eq "double precision" ) {
			$type = "double";
		    } else {
			die("Galacticus::Build::Functions::Functions_Generate_Output: type unsupported for C++-binding");
		    }
		    my $separator     = "";
		    my $fullSeparator = "";
		    my $argumentList  = "";
		    my $variableList  = "";
		    my $fullList      = "";
		    if ( $methodsCBound{$methodName}->{'pass'} eq "yes" ) {
			$argumentList .= $separator."void*";
			$variableList .= $separator."fortranSelf";
			$separator     = ",";
		    }
		    my @arguments;
		    if ( exists($methodsCBound{$methodName}->{'argument'}) ) {
			if ( UNIVERSAL::isa($methodsCBound{$methodName}->{'argument'},"ARRAY") ) {
			    push(@arguments,@{$methodsCBound{$methodName}->{'argument'}});
			} else {
			    push(@arguments,  $methodsCBound{$methodName}->{'argument'} );
			}
		    }
		    foreach my $argument ( @arguments ) {
			foreach my $intrinsic ( keys(%Fortran_Utils::intrinsicDeclarations) ) {
			    my $declarator = $Fortran_Utils::intrinsicDeclarations{$intrinsic};
			    if ( my @matches = $argument =~ m/$declarator->{'regEx'}/ ) {
				my $intrinsicName =                          $declarator->{'intrinsic' }  ;
				my $type          =                 $matches[$declarator->{'type'      }] ;
				my $attributeList =                 $matches[$declarator->{'attributes'}] ;
				$attributeList =~ s/^\s*,?\s*//;
				$attributeList =~ s/\s*$//;
				my @attributes = &Fortran_Utils::Extract_Variables($attributeList, keepQualifiers => 1, removeSpaces => 1);
				foreach my $attribute ( @attributes ) {
				    die("Galacticus::Build::Functions::Functions_Generate_Output:  attribute not supported for C++-binding")
					unless ( $attribute eq "intent(in)" );
				}
				my @variables     = split(/\s*,\s*/,$matches[$declarator->{'variables' }]);
				die("Galacticus::Build::Functions::Functions_Generate_Output: non-standard kinds are not supported for C++-binding")
				    if ( defined($type) );
				my $cType;
				if ( $intrinsicName eq "double precision" ) {
				    $cType = "double";
				} else {
				    die("Galacticus::Build::Functions::Functions_Generate_Output: type not supported for C++-binding");
				}
				$argumentList .=     $separator.join(",",map($cType       ,1..scalar(@variables)));
				$variableList .=     $separator.join(",",                            @variables  );
				$fullList     .= $fullSeparator.join(",",map($cType." ".$_,          @variables ));
				$separator     = ",";
				$fullSeparator = ",";
			    }
			}
		    }
		    # Build extern and class declarations.
		    $externCode .= " ".$type." ".$methodName."_C(".$argumentList.");\n";
		    my $classArgumentList = $argumentList;
		    $classArgumentList =~ s/^void\*,?//
			if ( $methodsCBound{$methodName}->{'pass'} eq "yes" );
		    $classCode  .= " ".$type." ".$methodName."(".$classArgumentList.");\n";
		    # Build the method.
		    $methodCode .= $type." ".$directive->{'name'}."Class::".$methodName." (".$fullList.") {\n";
		    $methodCode .= " return ".$methodName."_C(".$variableList.");\n";
		    $methodCode .= "}\n\n";
		}
		my $cBindings;
		$cBindings  = "// Generated automatically by Galacticus::Build::SourceTree::Process::FunctionClass\n";
		$cBindings .= "//  From: ".$tree->{'name'}."\n"
		    if ( $tree->{'type'} eq "file" );
		# Generate external linkage for creator, destructor, and method functions.
		$cBindings .= "extern \"C\"\n";
		$cBindings .= "{\n";
		$cBindings .= " void* ".$directive->{'name'}."();\n";
		$cBindings .= " void ".$directive->{'name'}."Destructor(void*);\n";
		$cBindings .= $externCode;
		$cBindings .= "}\n\n";
		# Create a class for this object.		
		$cBindings .= "class ".$directive->{'name'}."Class {\n";
		$cBindings .= "  void *fortranSelf;\n";
		$cBindings .= " public:\n";
		$cBindings .= " ".$directive->{'name'}."Class ();\n";
		$cBindings .= " ~".$directive->{'name'}."Class ();\n";
		$cBindings .= $classCode;
		$cBindings .= "};\n\n";	
		# Create a creator.
		$cBindings .= $directive->{'name'}."Class::".$directive->{'name'}."Class () {\n";
		$cBindings .= " fortranSelf=".$directive->{'name'}."();\n";
		$cBindings .= "};\n\n";
		# Create a destructor.
		$cBindings .= $directive->{'name'}."Class::~".$directive->{'name'}."Class () {\n";
		$cBindings .= " ".$directive->{'name'}."Destructor(fortranSelf);\n";
		$cBindings .= "};\n\n";
		# Create methods.
		$cBindings .= $methodCode;
		open(cHndl,">".$ENV{'BUILDPATH'}."/".$directive->{'name'}.".h");
		print cHndl $cBindings;
		close(cHndl);
	    }
	    
	    # Generate documentation.
	    my $documentation = "\\subsubsection{".$directive->{'descriptiveName'}."}\\label{sec:methods".ucfirst($directive->{'name'})."}\n\n";
	    $documentation   .= "Additional implementations for ".lc($directive->{'descriptiveName'})." are added using the {\\normalfont \\ttfamily ".$directive->{'name'}."} class.\n";
	    $documentation   .= "The implementation should be placed in a file containing the directive:\n";
	    $documentation   .= "\\begin{verbatim}\n";
	    $documentation   .= "!# <".$directive->{'name'}." name=\"".$directive->{'name'}."MyImplementation\">\n";
	    $documentation   .= "!# <description>A short description of the implementation.</description>\n";
	    $documentation   .= "!# </".$directive->{'name'}.">\n";
	    $documentation   .= "\\end{verbatim}\n";
	    $documentation   .= "where {\\normalfont \\ttfamily MyImplementation} is an appropriate name for the implemention. This file should be treated as a regular Fortran module, but without the initial {\\normalfont \\ttfamily module} and final {\\normalfont \\ttfamily end module} lines. That is, it may contain {\\normalfont \\ttfamily use} statements and variable declarations prior to the {\\normalfont \\ttfamily contains} line, and should contain all functions required by the implementation after that line. Function names should begin with {\\normalfont \\ttfamily ".&LaTeX_Breakable($directive->{'name'}."MyImplementation")."}. The file \\emph{must} define a type that extends the {\\normalfont \\ttfamily ".$directive->{'name'}."Class} class (or extends another type which is itself an extension of the {\\normalfont \\ttfamily ".$directive->{'name'}."Class} class), containing any data needed by the implementation along with type-bound functions required by the implementation. The following type-bound functions are required (unless inherited from the parent type):\n";
	    $documentation   .= "\\begin{description}\n";
	    # Create functions.
	    foreach my $methodName ( keys(%methods) ) {
		my $method = $methods{$methodName};
		$documentation   .= "\\item[{\\normalfont \\ttfamily ".$methodName."}] ".$method->{'description'};
		if ( exists($method->{'code'}) ) {
		    $documentation .= " A default implementation exists. If overridden the following interface must be used:\n";
		} else {
		    $documentation .= " Must have the following interface:\n";
		}
		$documentation   .= "\\begin{lstlisting}[language=Fortran,basicstyle=\\small\\ttfamily,escapechar=@,breaklines,prebreak=\\&,postbreak=\\&\\space\\space,columns=flexible,keepspaces=true,breakautoindent=true,breakindent=10pt]\n";
		# Insert arguments.
		my @arguments;
		if ( exists($method->{'argument'}) ) {
		    if ( UNIVERSAL::isa($method->{'argument'},"ARRAY") ) {
			push(@arguments,@{$method->{'argument'}});
		    } else {
			push(@arguments,  $method->{'argument'} );
		    }
		}
		unshift(@arguments,"class(".$directive->{'name'}."Class), intent(inout) :: self");
		my $argumentList = "";
		my $separator    = "";
		my @argumentDefinitions;
		foreach my $argument ( @arguments ) {
		    if ( $argument =~ $Fortran_Utils::variableDeclarationRegEx ) {
			my $intrinsic     = $1;
			my $type          = $2;
			my $attributeList = $3;
			my $variableList  = $4;
			my @variables  = &Fortran_Utils::Extract_Variables($variableList,keepQualifiers => 1,lowerCase => 0);
			my $declaration =
			{
			    intrinsic  => $intrinsic,
			    attributes => $attributeList,
			    variables  => \@variables
			}; 
			if ( defined($type) ) {
			    $type =~ s/\((.*)\)/$1/;
			    $declaration->{'type'} = $type;
			}
			if ( defined($attributeList) ) {
			    $attributeList =~ s/^\s*,\s*//;
			    my @attributes = &Fortran_Utils::Extract_Variables($attributeList,keepQualifiers => 1);
			    $declaration->{'attributes'} = \@attributes;
			}
			push(@argumentDefinitions,$declaration);
		    } else {
			print "Argument does not match expected pattern:\n\t".$argument."\n";
			die("Functions_Generate_Output: argument parse error");
		    }
		    (my $variables = $argument) =~ s/^.*::\s*(.*?)\s*$/$1/;
		    $argumentList .= $separator.$variables;
		    $separator     = ",";
		}
		my $type;
		my $category;
		if ( $method->{'type'} eq "void" ) {
		    $category = "subroutine";
		    $type     = "";
		} else {
		    $category = "function";
		    $type     = $method->{'type'}." ";
		}
		$documentation .= "   ".$type.$category." myImplementation".ucfirst($methodName)."(";
		$documentation .= $argumentList
		    unless ( $argumentList eq "" );
		$documentation .= ")\n";
		$documentation .= &Fortran_Utils::Format_Variable_Defintions(\@argumentDefinitions);
		$documentation .= "   end ".$type.$category." myImplementation".ucfirst($methodName)."\n";
		$documentation .= "\\end{lstlisting}\n\n";
	    }
	    $documentation   .= "\\end{description}\n\n";	    
	    $documentation   .= "Existing implementations are:\n";
	    $documentation   .= "\\begin{description}\n";
	    foreach my $className ( keys(%classes) ) {
		my $class = $classes{$className};
		$documentation   .= "\\item[{\\normalfont \\ttfamily ".$class->{'name'}."}] ".$class->{'description'};
		$documentation   .= " \\iflabelexists{phys:".$directive->{'name'}.":".$class->{'name'}."}{See \\S\\ref{phys:".$directive->{'name'}.":".$class->{'name'}."}.}{}\n";
	    }
	    $documentation   .= "\\end{description}\n\n";
	    system("mkdir -p doc/methods");
	    open(my $docHndl,">doc/methods/".$directive->{'name'}.".tex");
	    print $docHndl $documentation;
	    close($docHndl);
		
	    # Insert into tree.
	    &SourceTree::InsertPreContains ($node->{'parent'},$preContains );
	    &SourceTree::InsertPostContains($node->{'parent'},$postContains);

	}
	$node = &SourceTree::Walk_Tree($node,\$depth);
    }
}

sub LaTeX_Breakable {
    my $text = shift;
    $text =~ s/([a-z])([A-Z])/$1\\-$2/g;
    return $text;
}

1;
