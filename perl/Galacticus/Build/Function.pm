# Contains a Perl module which implements processing of "function" directives in the Galacticus build system.

package Function;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
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
use Sort::Topological qw(toposort);
use Fortran::Utils;
require Galacticus::Build::Hooks;
require Galacticus::Build::Dependencies;
require Fortran::Utils;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     function        => {parse => \&Functions_Parse_Directive        , generate => \&Functions_Generate_Output        },
     functionModules => {parse => \&Functions_Modules_Parse_Directive, generate => \&Functions_Modules_Generate_Output}
    );

sub Functions_Parse_Directive {
    # Parse content for a "function" directive.
    my $buildData = shift;

    # Assert that we have a prefix, currentDocument and directive.
    die("Galacticus::Build::Function::Functions_Parse_Directive: no currentDocument present" )
	unless ( exists($buildData->{'currentDocument'}           ) );
    die("Galacticus::Build::Function::Functions_Parse_Directive: no name present"            )
	unless ( exists($buildData->{'currentDocument'}->{'name'} ) );

    # Store the name and associated file name
    my $directive = $buildData->{'directive'};
    my $className = $buildData->{'currentDocument'}->{'name'};
    my $fileName  = $buildData->{'currentFileName'};
    push(@{$buildData->{$directive}->{'classes'}},{name => $className, file => $fileName, description => $buildData->{'currentDocument'}->{'description'}});

}

sub Functions_Modules_Parse_Directive {
    # Parse content for a "functionModules" directive.
    my $buildData = shift;

    # Assert that we have a prefix, currentDocument and directive.
    die("Galacticus::Build::Function::Functions_Parse_Directive: no currentDocument present" )
	unless ( exists($buildData->{'currentDocument'}           ) );
    die("Galacticus::Build::Function::Functions_Parse_Directive: no name present"            )
	unless ( exists($buildData->{'currentDocument'}->{'name'} ) );

    # Store the name and associated file name
    my $directive = $buildData->{'directive'};
    my $className = $buildData->{'currentDocument'}->{'name'};
    my $fileName  = $buildData->{'currentFileName'};
    push(@{$buildData->{$directive}->{'functionModules'}},{name => $className, file => $fileName});
}

sub Functions_Generate_Output {
    # Generate output for a "function" directive.
    my $buildData = shift;

    # Assert that we have a file name and directive present.
    die("Galacticus::Build::Function::Functions_Parse_Directive: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );
    die("Galacticus::Build::Function::Functions_Parse_Directive: no directive present")
	unless ( exists($buildData->{'directive'}) );

    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;

    # Specify unit opening regexs.
    my %unitOpeners = (
	# Find module openings, avoiding module procedures.
	module             => { unitName => 1, regEx => "^\\s*module\\s+(?!procedure\\s)([a-z0-9_]+)" },
	# Find program openings.
	program            => { unitName => 1, regEx => "^\\s*program\\s+([a-z0-9_]+)" },
	# Find subroutine openings, allowing for pure, elemental and recursive subroutines.
	subroutine         => { unitName => 2, regEx => "^\\s*(pure\\s+|elemental\\s+|recursive\\s+)*\\s*subroutine\\s+([a-z0-9_]+)"},
	# Find function openings, allowing for pure, elemental and recursive functions, and different function types.
	function           => { unitName => 5, regEx => "^\\s*(pure\\s+|elemental\\s+|recursive\\s+)*\\s*(real|integer|double\\s+precision|character|logical)*\\s*(\\((kind|len)=[\\w\\d]*\\))*\\s*function\\s+([a-z0-9_]+)"},
	# Find interfaces.
	interface          => { unitName => 2, regEx => "^\\s*(abstract\\s+)??interface\\s+([a-z0-9_\\(\\)\\/\\+\\-\\*\\.=]*)"},
	# Find types.
	type               => { unitName => 3, regEx => "^\\s*type\\s*(,\\s*abstract\\s*|,\\s*public\\s*|,\\s*private\\s*|,\\s*extends\\s*\\([a-zA-Z0-9_]+\\)\\s*)*(::)??\\s*([a-z0-9_]+)\\s*\$"}
	);

    # Specify unit closing regexs.
    my %unitClosers = (
	module             => { unitName => 1, regEx => "^\\s*end\\s+module\\s+([a-z0-9_]+)" },
	program            => { unitName => 1, regEx => "^\\s*end\\s+program\\s+([a-z0-9_]+)" },
	subroutine         => { unitName => 1, regEx => "^\\s*end\\s+subroutine\\s+([a-z0-9_]+)"},
	function           => { unitName => 1, regEx => "^\\s*end\\s+function\\s+([a-z0-9_]+)"},
	interface          => { unitName => 1, regEx => "^\\s*end\\s+interface"},
	type               => { unitName => 1, regEx => "^\\s*end\\s+type\\s+([a-z0-9_]+)"}
	);

    # Specify regexs for intrinsic variable declarations.
    my %intrinsicDeclarations = (
	integer   => { intrinsic => "integer", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)integer(?-i)\\s*(\\(\\s*kind\\s*=\\s*[a-zA-Z0-9_]+\\s*\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
	real      => { intrinsic => "real", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)real(?-i)\\s*(\\(\\s*kind\\s*=\\s*[a-zA-Z0-9_]+\\s*\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9\\._,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
	double    => { intrinsic => "double precision", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)double\\s+precision(?-i)\\s*(\\(\\s*kind\\s*=\\s*[a-zA-Z0-9_]+\\s*\\))*([\\sa-zA-Z0-9_,:=\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9\\._,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
	logical   => { intrinsic => "logical", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)logical(?-i)\\s*(\\(\\s*kind\\s*=\\s*[a-zA-Z0-9_]+\\s*\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
	character => { intrinsic => "character", type => 1, attributes => 4, variables => 5, regEx => "^\\s*(?i)character(?-i)\\s*(\\((\\s*(len|kind)\\s*=\\s*[a-zA-Z0-9_,\\+\\-\\*\\(\\)]+\\s*)+\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
	procedure => { intrinsic => "procedure", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)procedure(?-i)\\s*(\\(\\s*[a-zA-Z0-9_]*\\s*\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9_,:=>\\+\\-\\*\\/\\(\\)]+)\\s*\$" },
	);

    # Extract the directive.
    my $directive = $buildData->{'directive'};

    # Find methods.
    my @methods;
    if ( exists($buildData->{'method'}) ) {
	if ( UNIVERSAL::isa($buildData->{'method'},"ARRAY") ) {
	    push(@methods,@{$buildData->{'method'}});
	} elsif ( UNIVERSAL::isa($buildData->{'method'},"HASH") ) {
	    push(@methods,$buildData->{'method'});
	} else {
	    push(@methods,  $buildData->{'method'} );
	}
    }

    # Determine if any methods request that C-bindings be produced.
    my @methodsCBound;
    foreach ( @methods ) {
	push(@methodsCBound,$_)
	    if ( exists($_->{'bindC'}) && $_->{'bindC'} eq "true" );
    }

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::Function\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "!  Time: ".$now."\n\n";

    # Add public functions.
    $buildData->{'content'} .= "   public :: ".$directive.",".$directive."Class";
    $buildData->{'content'} .= ", ".$_->{'name'}
	foreach ( @{$buildData->{$directive}->{'classes'}});
    $buildData->{'content'} .= "\n\n";

    # Add variable tracking module initialization status.
    $buildData->{'content'} .= "   logical, private :: moduleInitialized=.false.\n\n";

    # Generate the function object.
    $buildData->{'content'} .= "   type :: ".$directive."Class\n";
    $buildData->{'content'} .= "    private\n";
    $buildData->{'content'} .= "    contains\n";
    $buildData->{'content'} .= "    !@ <objectMethods>\n";
    $buildData->{'content'} .= "    !@   <object>".$directive."Class</object>\n";
    foreach my $method ( @methods ) {
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
	    foreach my $intrinsic ( keys(%intrinsicDeclarations) ) {
		my $declarator = $intrinsicDeclarations{$intrinsic};
		if ( my @matches = $argument =~ m/$declarator->{'regEx'}/ ) {
		    my $intrinsicName =                          $declarator->{'intrinsic' }    ;
		    my $type          =                 $matches[$declarator->{'type'      }-1] ;
		    my $attributeList =                 $matches[$declarator->{'attributes'}-1] ;
		    $attributeList =~ s/^\s*,?\s*//;
		    $attributeList =~ s/\s*$//;
		    my @attributes = &Fortran_Utils::Extract_Variables($attributeList, keepQualifiers => 1, removeSpaces => 1);
		    my @variables     = split(/\s*,\s*/,$matches[$declarator->{'variables' }-1]);
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
	$buildData->{'content'} .= "    !@   <objectMethod>\n";
	$buildData->{'content'} .= "    !@     <method>".$method->{'name'}."</method>\n";
	$buildData->{'content'} .= "    !@     <type>".$method->{'type'}."</type>\n";
	$buildData->{'content'} .= "    !@     <arguments>".$argumentList."</arguments>\n";
	$buildData->{'content'} .= "    !@     <description>".$method->{'description'}."</description>\n";
	$buildData->{'content'} .= "    !@   </objectMethod>\n";
    }
    $buildData->{'content'} .= "    !@ </objectMethods>\n";
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
    foreach ( @methods ) {
	my $extension = "Null";
	$extension = ""
	    if ( exists($_->{'code'}) );
	$methodTable->add("",$_->{'name'},$_->{'name'}.$extension);
    }
    $buildData->{'content'} .= $methodTable->table();
    $buildData->{'content'} .= "   end type ".$directive."Class\n\n";

    # Insert interface to class constructors.
    $buildData->{'content'} .= "   interface ".$directive."\n";
    $buildData->{'content'} .= "    module procedure ".$directive."ConstructorDefault\n";
    $buildData->{'content'} .= "    module procedure ".$directive."ConstructorNamed\n";
    $buildData->{'content'} .= "   end interface\n";

    # Scan implementation code to determine dependencies.
    my %dependencies;
    my %classes;
    foreach my $class ( @{$buildData->{$directive}->{'classes'}} ) {
	open(my $classFile,$class->{'file'});
	until ( eof($classFile) ) {
	    &Fortran_Utils::Get_Fortran_Line($classFile,my $rawLine, my $processedLine, my $bufferedComments);
	    if ( $processedLine =~ m/^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*$directive([a-z0-9_]+)\s*$/i ) {
		my $extends = $2;
		my $type    = $directive.$4;
		$class  ->{'type'   } = $type;
		$class  ->{'extends'} = $extends;
		$classes  {$type    } = $class;
		push(@{$dependencies{$extends}},$type);
	    }
	}
    }
    my @unsortedClasses = keys(%classes);
    my @sortedClasses = toposort(sub { @{$dependencies{$_[0]} || []}; }, \@unsortedClasses );
    @{$buildData->{$directive}->{'classes'}} = map($classes{$_},@sortedClasses);

    # Insert pre-contains implementation code.
    foreach my $class ( @{$buildData->{$directive}->{'classes'}} ) {
	my $unitDepth = 0;
	open(my $classFile,$class->{'file'});
	until ( eof($classFile) ) {
	    &Fortran_Utils::Get_Fortran_Line($classFile,my $rawLine, my $processedLine, my $bufferedComments);
	    foreach my $unitType ( keys(%unitOpeners) ) {
		++$unitDepth
		    if ( $processedLine =~ m/$unitOpeners{$unitType}->{"regEx"}/i );
		--$unitDepth
		    if ( $processedLine =~ m/$unitClosers{$unitType}->{"regEx"}/i );
	    }
	    next
		if ( $processedLine =~ m/^\s*use\s+[a-zA-Z0-9_\s:\,]+/ );
	    last
		if ( $processedLine =~ m/^\s*contains\s*$/i && $unitDepth == 0 );
	    $buildData->{'content'} .= $rawLine;
	}
	close($classFile);
    }

    # Add method name parameter.
    $buildData->{'content'} .= "   ! Method name parameter.\n";
    $buildData->{'content'} .= "   type(varying_string) :: ".$directive."Method\n\n";

    # Add default implementation.
    $buildData->{'content'} .= "   ! Default ".$directive." object.\n";
    $buildData->{'content'} .= "   class(".$directive."Class), private , pointer :: ".$directive."Default => null()\n";
    $buildData->{'content'} .= "   !\$omp threadprivate(".$directive."Default)\n"
	if ( exists($buildData->{'defaultThreadPrivate'}) && $buildData->{'defaultThreadPrivate'} eq "yes" );
    $buildData->{'content'} .= "\n";

    # If we need to generate C-bindings, insert a wrapper class to permit passing of polymorphic pointers between Fortran and C++.
    if ( @methodsCBound ) {
	$buildData->{'content'} .= "   type :: ".$directive."Wrapper\n";
	$buildData->{'content'} .= "     class(".$directive."Class), pointer :: wrappedObject\n";
	$buildData->{'content'} .= "   end type ".$directive."Wrapper\n\n";
    }

    # Insert "contains" separator.
    $buildData->{'content'} .= "contains\n\n";

    # Create default constructor.
    $buildData->{'content'} .= "   function ".$directive."ConstructorDefault()\n";
    $buildData->{'content'} .= "      !% Return a pointer to the default {\\tt ".$directive."} object.\n";
    $buildData->{'content'} .= "      implicit none\n";
    $buildData->{'content'} .= "      class(".$directive."Class), pointer :: ".$directive."ConstructorDefault\n\n";
    $buildData->{'content'} .= "      if (.not.associated(".$directive."Default)) call ".$directive."Initialize()\n";
    $buildData->{'content'} .= "      ".$directive."ConstructorDefault => ".$directive."Default\n";
    $buildData->{'content'} .= "      return\n";
    $buildData->{'content'} .= "   end function ".$directive."ConstructorDefault\n\n";

    # Create named constructor.
    $buildData->{'content'} .= "   function ".$directive."ConstructorNamed(typeName)\n";
    $buildData->{'content'} .= "      !% Return a pointer to a newly created {\\tt ".$directive."} object of the specified type.\n";
    $buildData->{'content'} .= "      use ISO_Varying_String\n";
    $buildData->{'content'} .= "      use Galacticus_Error\n";
    $buildData->{'content'} .= "      implicit none\n";
    $buildData->{'content'} .= "      class(".$directive."Class), pointer :: ".$directive."ConstructorNamed\n\n";
    $buildData->{'content'} .= "      character(len=*), intent(in   ) :: typeName\n";
    $buildData->{'content'} .= "      type(varying_string) :: message\n\n";
    $buildData->{'content'} .= "      select case (trim(typeName))\n";
    foreach my $class ( @{$buildData->{$directive}->{'classes'}} ) {
	(my $name = $class->{'name'}) =~ s/^$directive//;
	$name = lcfirst($name);
	$buildData->{'content'} .= "     case ('".$name."')\n";
	$buildData->{'content'} .= "        allocate(".$class->{'name'}." :: ".$directive."ConstructorNamed)\n";
	$buildData->{'content'} .= "        select type (".$directive."ConstructorNamed)\n";
	$buildData->{'content'} .= "          type is (".$class->{'name'}.")\n";
	$buildData->{'content'} .= "            ".$directive."ConstructorNamed=".$class->{'name'}."()\n";
	$buildData->{'content'} .= "         end select\n";
    }
    $buildData->{'content'} .= "      case default\n";
    $buildData->{'content'} .= "         message='Unrecognized type \"'//trim(typeName)//'\" Available options are:'\n";
    my @classNames;
    push(@classNames,$_->{'name'})
	foreach ( @{$buildData->{$directive}->{'classes'}} );
    foreach ( sort(@classNames) ) {
	(my $name = $_) =~ s/^$directive//;
	$buildData->{'content'} .= "        message=message//char(10)//'   -> ".lcfirst($name)."'\n";
    }
    $buildData->{'content'} .= "         call Galacticus_Error_Report('".$directive."ConstructorNamed',message)\n";
    $buildData->{'content'} .= "      end select\n";
    $buildData->{'content'} .= "      return\n";
    $buildData->{'content'} .= "   end function ".$directive."ConstructorNamed\n\n";

    # Create initialization function.
    $buildData->{'content'} .= "   subroutine ".$directive."Initialize()\n";
    $buildData->{'content'} .= "      !% Initialize the default {\\tt ".$directive."} object.\n";
    $buildData->{'content'} .= "      use ISO_Varying_String\n";
    $buildData->{'content'} .= "      use Input_Parameters\n";
    $buildData->{'content'} .= "      use Galacticus_Error\n";
    $buildData->{'content'} .= "      implicit none\n";
    $buildData->{'content'} .= "      type(varying_string) :: message\n\n";
    $buildData->{'content'} .= "      if (.not.moduleInitialized) then\n";
    $buildData->{'content'} .= "         !\$omp critical (".$directive."Initialization)\n";
    $buildData->{'content'} .= "         if (.not.moduleInitialized) then\n";
    $buildData->{'content'} .= "            !@ <inputParameter>\n";
    $buildData->{'content'} .= "            !@   <name>".$directive."Method</name>\n";
    $buildData->{'content'} .= "            !@   <defaultValue>".$buildData->{'default'}."</defaultValue>\n";
    $buildData->{'content'} .= "            !@   <attachedTo>module</attachedTo>\n";
    $buildData->{'content'} .= "            !@   <description>\n";
    $buildData->{'content'} .= "            !@     The method to be used for {\\tt ".$directive."}.\n";
    $buildData->{'content'} .= "            !@   </description>\n";
    $buildData->{'content'} .= "            !@   <type>string</type>\n";
    $buildData->{'content'} .= "            !@   <cardinality>1</cardinality>\n";
    $buildData->{'content'} .= "            !@ </inputParameter>\n";
    $buildData->{'content'} .= "            call Get_Input_Parameter('".$directive."Method',".$directive."Method,defaultValue='".$buildData->{'default'}."')\n";
    $buildData->{'content'} .= "            moduleInitialized=.true.\n";
    $buildData->{'content'} .= "         end if\n";
    $buildData->{'content'} .= "         !\$omp end critical (".$directive."Initialization)\n";
    $buildData->{'content'} .= "      end if\n";
    $buildData->{'content'} .= "      select case (char(".$directive."Method))\n";
    foreach my $class ( @{$buildData->{$directive}->{'classes'}} ) {
	(my $name = $class->{'name'}) =~ s/^$directive//;
	$name = lcfirst($name);
	$buildData->{'content'} .= "     case ('".$name."')\n";
	$buildData->{'content'} .= "        allocate(".$class->{'name'}." :: ".$directive."Default)\n";
	$buildData->{'content'} .= "        select type (".$directive."Default)\n";
	$buildData->{'content'} .= "          type is (".$class->{'name'}.")\n";
	$buildData->{'content'} .= "            ".$directive."Default=".$class->{'name'}."()\n";
	$buildData->{'content'} .= "         end select\n";
    }
    $buildData->{'content'} .= "      case default\n";
    $buildData->{'content'} .= "         message='Unrecognized option for [".$directive."Method](='//".$directive."Method//'). Available options are:'\n";
    foreach ( sort(@classNames) ) {
	(my $name = $_) =~ s/^$directive//;
	$buildData->{'content'} .= "        message=message//char(10)//'   -> ".lcfirst($name)."'\n";
    }
    $buildData->{'content'} .= "         call Galacticus_Error_Report('".$directive."Initialize',message)\n";
    $buildData->{'content'} .= "      end select\n";
    $buildData->{'content'} .= "      return\n";
    $buildData->{'content'} .= "   end subroutine ".$directive."Initialize\n\n";

    # Create functions.
    foreach my $method ( @methods ) {
	# Insert arguments.
	my @arguments;
	if ( exists($method->{'argument'}) ) {
	    if ( UNIVERSAL::isa($method->{'argument'},"ARRAY") ) {
		push(@arguments,@{$method->{'argument'}});
	    } else {
		push(@arguments,  $method->{'argument'} );
	    }
	}
	my $argumentList = "";
	my $argumentCode = "      class(".$directive."Class), intent(inout) :: self\n";
	my $separator = "";
	foreach my $argument ( @arguments ) {
	    (my $variables = $argument) =~ s/^.*::\s*(.*?)\s*$/$1/;
	    $argumentList .= $separator.$variables;
	    $argumentCode .= "      ".$argument."\n";
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
	my $extension = "Null";
	$extension = ""
	    if ( exists($method->{'code'}) );
	$buildData->{'content'} .= "   ".$type.$category." ".$method->{'name'}.$extension."(self";
	$buildData->{'content'} .= ",".$argumentList
	    unless ( $argumentList eq "" );
	$buildData->{'content'} .= ")\n";
	$buildData->{'content'} .= "      !% ".$method->{'description'}."\n";
	if ( exists($method->{'code'}) ) {
	    if ( exists($method->{'modules'}) ) {
		$buildData->{'content'} .= "      use ".$_."\n"
		    foreach ( split(/\s+/,$method->{'modules'}) );
	    }
	} else {
	    $buildData->{'content'} .= "      use Galacticus_Error\n";
	}
	$buildData->{'content'} .= "      implicit none\n";
	$buildData->{'content'} .= $argumentCode;
	if ( exists($method->{'code'}) ) {
	    my $code = "      ".$method->{'code'};
	    $code =~ s/\n/\n      /g;
	    $buildData->{'content'} .= $code."\n";
	} else {
	    $buildData->{'content'} .= "      call Galacticus_Error_Report('".$method->{'name'}."Null','this is a null method - initialize the ".$directive." object before use')\n";
	}
	$buildData->{'content'} .= "      return\n";
	$buildData->{'content'} .= "   end ".$category." ".$method->{'name'}.$extension."\n\n";
    }
    # Generate C-bindings if required.
    if ( @methodsCBound ) {
	# C-bound default constructor. Here, we use a wrapper object which contains a pointer to the default polymorphic Fortran
	# object. This wrapper is then passed back to the calling C++ function so that it can be stored in the appropriate C++
	# class.
	$buildData->{'content'} .= "   function ".$directive."_C() bind(c,name='".$directive."')\n";
	$buildData->{'content'} .= "     implicit none\n";
	$buildData->{'content'} .= "     type(c_ptr) :: ".$directive."_C\n";
	$buildData->{'content'} .= "     type(".$directive."Wrapper), pointer :: wrapper\n";
	$buildData->{'content'} .= "      if (.not.associated(".$directive."Default)) call ".$directive."Initialize()\n";
	$buildData->{'content'} .= "       allocate(wrapper)\n";
	$buildData->{'content'} .= "       wrapper%wrappedObject => ".$directive."Default\n";
	$buildData->{'content'} .= "       ".$directive."_C=c_loc(wrapper)\n";
	$buildData->{'content'} .= "     return\n";
	$buildData->{'content'} .= "   end function ".$directive."_C\n\n";
	# C-bound destructor. We simply deallocate the wrapper object, letting the associated finalizor clean up the Fortran
	# object.
	$buildData->{'content'} .= "   subroutine ".$directive."Destructor_C(wrapperC) bind(c,name='".$directive."Destructor')\n";
	$buildData->{'content'} .= "     implicit none\n";
	$buildData->{'content'} .= "     type(c_ptr), intent(in   ), value :: wrapperC\n";
	$buildData->{'content'} .= "     type(".$directive."Wrapper), pointer :: wrapper\n\n";
	$buildData->{'content'} .= "     call c_f_pointer(wrapperC,wrapper)\n";
	$buildData->{'content'} .= "     deallocate(wrapper)\n";
	$buildData->{'content'} .= "     return\n";
	$buildData->{'content'} .= "   end subroutine ".$directive."Destructor_C\n\n";
	# Generate method functions.
	foreach ( @methodsCBound ) {
	    my @arguments;
	    if ( exists($_->{'argument'}) ) {
		if ( UNIVERSAL::isa($_->{'argument'},"ARRAY") ) {
		    push(@arguments,@{$_->{'argument'}});
		} else {
		    push(@arguments,  $_->{'argument'} );
		}
	    }
	    my $separator    = "";
	    my $argumentList = "";
	    foreach my $argument ( @arguments ) {
		foreach my $intrinsic ( keys(%intrinsicDeclarations) ) {
		    my $declarator = $intrinsicDeclarations{$intrinsic};
		    if ( my @matches = $argument =~ m/$declarator->{'regEx'}/ ) {
			my $intrinsicName =                          $declarator->{'intrinsic' }    ;
			my $type          =                 $matches[$declarator->{'type'      }-1] ;
			my $attributeList =                 $matches[$declarator->{'attributes'}-1] ;
			$attributeList =~ s/^\s*,?\s*//;
			$attributeList =~ s/\s*$//;
			my @attributes = &Fortran_Utils::Extract_Variables($attributeList, keepQualifiers => 1, removeSpaces => 1);
			foreach my $attribute ( @attributes ) {
			    die("Galacticus::Build::Functions::Functions_Generate_Output:  attribute not supported for C++-binding")
				unless ( $attribute eq "intent(in)" );
			}
			my @variables     = split(/\s*,\s*/,$matches[$declarator->{'variables' }-1]);
			die("Galacticus::Build::Functions::Functions_Generate_Output: non-standard kinds are not supported for C++-binding")
			    if ( defined($type) );
			$argumentList .= $separator.join(",",@variables);
			$separator     = ",";
		    }
		}
	    }
	    $buildData->{'content'} .= "  double precision function ".$_->{'name'}."_C(";
	    $buildData->{'content'} .= "wrapperC".$separator
		if ( $_->{'pass'} eq "yes" );
	    $buildData->{'content'} .= $argumentList.") bind(c,name='".$_->{'name'}."_C')\n";
	    $buildData->{'content'} .= "     implicit none\n";
	    $buildData->{'content'} .= "     type(c_ptr), intent(in   ), value :: wrapperC\n";
	    foreach my $argument( @arguments ) {
		(my $argumentInteroperable = $argument) =~ s/(\s*::)/, value$1/;
		$buildData->{'content'} .= "     ".$argumentInteroperable."\n"
	    }
	    $buildData->{'content'} .= "     type(".$directive."Wrapper), pointer :: wrapper\n";
	    $buildData->{'content'} .= "     call c_f_pointer(wrapperC,wrapper)\n";
	    $buildData->{'content'} .= "     ".$_->{'name'}."_C=wrapper\%wrappedObject\%".$_->{'name'}."(".$argumentList.")\n";
	    $buildData->{'content'} .= "     return\n";
	    $buildData->{'content'} .= "   end function ".$_->{'name'}."_C\n\n";
	}
    }
    # Insert post-contains implementation code.
    foreach my $class ( @{$buildData->{$directive}->{'classes'}} ) {
	my $unitDepth     = 0;
	my $containsFound = 0;
	open(my $classFile,$class->{'file'});
	until ( eof($classFile) ) {
	    &Fortran_Utils::Get_Fortran_Line($classFile,my $rawLine, my $processedLine, my $bufferedComments);
	    foreach my $unitType ( keys(%unitOpeners) ) {
		++$unitDepth
		    if ( $processedLine =~ m/$unitOpeners{$unitType}->{"regEx"}/i );
		--$unitDepth
		    if ( $processedLine =~ m/$unitClosers{$unitType}->{"regEx"}/i );
	    }
	    $buildData->{'content'} .= $rawLine
		if ( $containsFound == 1 );
	    $containsFound = 1
		if ( $processedLine =~ m/^\s*contains\s*$/i && $unitDepth == 0 );
	}
	close($classFile);
    }
    # Generate C-bindings here if required.
    if ( @methodsCBound ){
	# Iterate over methods and generate the necessary code.
	my $externCode;
	my $classCode;
	my $methodCode;
	foreach ( @methodsCBound ) {
	    my $type;
	    switch ( $_->{'type'} ) {
		case ( "double precision" ) {
		    $type = "double";
		}
		else {
		    die("Galacticus::Build::Functions::Functions_Generate_Output: type unsupported for C++-binding");
		}
	    }
	    my $separator     = "";
	    my $fullSeparator = "";
	    my $argumentList  = "";
	    my $variableList  = "";
	    my $fullList      = "";
	    if ( $_->{'pass'} eq "yes" ) {
		$argumentList .= $separator."void*";
		$variableList .= $separator."fortranSelf";
		$separator     = ",";
	    }
	    my @arguments;
	    if ( exists($_->{'argument'}) ) {
		if ( UNIVERSAL::isa($_->{'argument'},"ARRAY") ) {
		    push(@arguments,@{$_->{'argument'}});
		} else {
		    push(@arguments,  $_->{'argument'} );
		}
	    }
	    foreach my $argument ( @arguments ) {
		foreach my $intrinsic ( keys(%intrinsicDeclarations) ) {
		    my $declarator = $intrinsicDeclarations{$intrinsic};
		    if ( my @matches = $argument =~ m/$declarator->{'regEx'}/ ) {
			my $intrinsicName =                          $declarator->{'intrinsic' }    ;
			my $type          =                 $matches[$declarator->{'type'      }-1] ;
			my $attributeList =                 $matches[$declarator->{'attributes'}-1] ;
			$attributeList =~ s/^\s*,?\s*//;
			$attributeList =~ s/\s*$//;
			my @attributes = &Fortran_Utils::Extract_Variables($attributeList, keepQualifiers => 1, removeSpaces => 1);
			foreach my $attribute ( @attributes ) {
			    die("Galacticus::Build::Functions::Functions_Generate_Output:  attribute not supported for C++-binding")
				unless ( $attribute eq "intent(in)" );
			}
			my @variables     = split(/\s*,\s*/,$matches[$declarator->{'variables' }-1]);
			die("Galacticus::Build::Functions::Functions_Generate_Output: non-standard kinds are not supported for C++-binding")
			    if ( defined($type) );
			my $cType;
			switch ( $intrinsicName ) {
			    case ( "double precision" ) {
				$cType = "double";
			    }
			    else {
				die("Galacticus::Build::Functions::Functions_Generate_Output: type not supported for C++-binding");
			    }
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
	    $externCode .= " ".$type." ".$_->{'name'}."_C(".$argumentList.");\n";
	    my $classArgumentList = $argumentList;
	    $classArgumentList =~ s/^void\*,?//
		if ( $_->{'pass'} eq "yes" );
	    $classCode  .= " ".$type." ".$_->{'name'}."(".$classArgumentList.");\n";
	    # Build the method.
	    $methodCode .= $type." ".$directive."Class::".$_->{'name'}." (".$fullList.") {\n";
	    $methodCode .= " return ".$_->{'name'}."_C(".$variableList.");\n";
	    $methodCode .= "}\n\n";
	}
	my $cBindings;
	$cBindings  = "// Generated automatically by Galacticus::Build::Function\n";
	$cBindings .= "//  From: ".$buildData->{'fileName'}."\n";
	$cBindings .= "//  Time: ".$now."\n\n";
	# Generate external linkage for creator, destructor, and method functions.
	$cBindings .= "extern \"C\"\n";
	$cBindings .= "{\n";
	$cBindings .= " void* ".$directive."();\n";
	$cBindings .= " void ".$directive."Destructor(void*);\n";
	$cBindings .= $externCode;
	$cBindings .= "}\n\n";
	# Create a class for this object.
	$cBindings .= "class ".$directive."Class {\n";
	$cBindings .= "  void *fortranSelf;\n";
	$cBindings .= " public:\n";
	$cBindings .= " ".$directive."Class ();\n";
	$cBindings .= " ~".$directive."Class ();\n";
	$cBindings .= $classCode;
	$cBindings .= "};\n\n";	
	# Create a creator.
	$cBindings .= $directive."Class::".$directive."Class () {\n";
	$cBindings .= " fortranSelf=".$directive."();\n";
	$cBindings .= "};\n\n";
	# Create a destructor.
	$cBindings .= $directive."Class::~".$directive."Class () {\n";
	$cBindings .= " ".$directive."Destructor(fortranSelf);\n";
	$cBindings .= "};\n\n";
	# Create methods.
	$cBindings .= $methodCode;
	open(cHndl,">work/build/".$directive.".h");
	print cHndl $cBindings;
	close(cHndl);
    }
    # Generate documentation.
    my $documentation = "\\subsubsection{".$buildData->{'descriptiveName'}."}\\label{sec:methods".ucfirst($directive)."}\n\n";
    $documentation   .= "Additional implementations for ".lc($buildData->{'descriptiveName'})." are added using the {\\tt ".$directive."} class.\n";
    $documentation   .= "The implementation should be placed in a file containing the directive:\n";
    $documentation   .= "\\begin{verbatim}\n";
    $documentation   .= "!# <".$directive." name=\"".$directive."MyImplementation\">\n";
    $documentation   .= "!# <description>A short description of the implementation.</description>\n";
    $documentation   .= "!# </".$directive.">\n";
    $documentation   .= "\\end{verbatim}\n";
    $documentation   .= "where {\\tt MyImplementation} is an appropriate name for the implemention. This file should be treated as a regular Fortran module, but without the initial {\\tt module} and final {\\tt end module} lines. That is, it may contain {\\tt use} statements and variable declarations prior to the {\\tt contains} line, and should contain all functions required by the implementation after that line. Function names should begin with {\\tt ".$directive."MyImplementation}. The file \\emph{must} define a type that extends the {\\tt ".$directive."Class} class (or extends another type which is itself an extension of the {\\tt ".$directive."Class} class), containing any data needed by the implementation along with type-bound functions required by the implementation. The following type-bound functions are required (unless inherited from the parent type):\n";
    $documentation   .= "\\begin{description}\n";
    # Create functions.
    foreach my $method ( @methods ) {
	$documentation   .= "\\item[{\\tt ".$method->{'name'}."}] ".$method->{'description'}." Must have the following interface:\n";
	$documentation   .= "\\begin{verbatim}\n";
	# Insert arguments.
	my @arguments;
	if ( exists($method->{'argument'}) ) {
	    if ( UNIVERSAL::isa($method->{'argument'},"ARRAY") ) {
		push(@arguments,@{$method->{'argument'}});
	    } else {
		push(@arguments,  $method->{'argument'} );
	    }
	}
	unshift(@arguments,"class(".$directive."Class), intent(inout) :: self");
	my $argumentList = "";
	my $separator    = "";
	my @argumentDefinitions;
	foreach my $argument ( @arguments ) {
	    if ( $argument =~ $Fortran_Utils::variableDeclarationRegEx ) {
		my $intrinsic     = $1;
		my $type          = $2;
		my $attributeList = $3;
		my $variableList  = $4;
		my @variables  = &Fortran_Utils::Extract_Variables($variableList,keepQualifiers => 1);
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
	$documentation .= "   ".$type.$category." myImplementation".ucfirst($method->{'name'})."(self";
	$documentation .= ",".$argumentList
	    unless ( $argumentList eq "" );
	$documentation .= ")\n";
	$documentation .= &Fortran_Utils::Format_Variable_Defintions(\@argumentDefinitions);
	$documentation .= "   end ".$type.$category." myImplementation".ucfirst($method->{'name'})."\n";
	$documentation .= "\\end{verbatim}\n\n";
    }
    $documentation   .= "\\end{description}\n\n";


    $documentation   .= "Existing implementations are:\n";
    $documentation   .= "\\begin{description}\n";
    foreach my $class ( @{$buildData->{$directive}->{'classes'}} ) {
	$documentation   .= "\\item[{\\tt ".$class->{'name'}."}] ".$class->{'description'};
	$documentation   .= " \\iflabelexists{phys:".$directive.":".$class->{'name'}."}{See \\S\\ref{phys:".$directive.":".$class->{'name'}."}.}{}\n";
    }
    $documentation   .= "\\end{description}\n\n";
	
    system("mkdir -p doc/methods");
    open(my $docHndl,">doc/methods/".$directive.".tex");
    print $docHndl $documentation;
    close($docHndl);
}

sub Functions_Modules_Generate_Output {
    # Generate output for a "functionModules" directive.
    my $buildData = shift;

    # Assert that we have a file name and directive present.
    die("Galacticus::Build::Function::Functions_Parse_Directive: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );
    die("Galacticus::Build::Function::Functions_Parse_Directive: no directive present")
	unless ( exists($buildData->{'directive'}) );

    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;

    # Extract the directive.
    my $directive = $buildData->{'directive'};

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::Function\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "!  Time: ".$now."\n\n";

    # Insert module use statements from implementation code.
    my %modules;
    foreach my $class ( @{$buildData->{$directive}->{'functionModules'}} ) {
	my $unitDepth = 0;
	open(my $classFile,$class->{'file'});
	until ( eof($classFile) ) {
	    &Fortran_Utils::Get_Fortran_Line($classFile,my $rawLine, my $processedLine, my $bufferedComments);
	    if ( $processedLine =~ m/^\s*use\s+([a-zA-Z0-9_\s:\,]+)/ ){
		$modules{$1} = $1;
	    }
	    last
		if ( $processedLine =~ m/^\s*contains\s*$/i && $unitDepth == 0 );
	}
	close($classFile);
    }
    $buildData->{'content'} .= "use ".$_
	foreach ( sort(keys(%modules)) );
}

1;
