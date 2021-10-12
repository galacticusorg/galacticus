#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Build::SourceTree;
use Galacticus::Build::SourceTree::Parse::Declarations;
use List::ExtraUtils;
use XML::Simple;
use Data::Dumper;
use Sort::Topo;
use Text::Template 'fill_in_string';
use Scalar::Util qw(reftype);

# Build interfaces for libgalacticus.
# Andrew Benson (28-September-2021)

# Initialize a structure which will hold the generated code.
my $code;
$code  ->{'units'} = [];

# Initialize a structure which will hold the Python interfaces.
my $python;
$python->{'c_lib'} = [];

# Get an XML parser.
my $xml = new XML::Simple;

# Get directive locations.
my $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");

# Get state storable information .
my $stateStorables     = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml"    );

# Get the list of functionClasses that should be compiled into the library.
my $libraryFunctionClasses = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/source/libraryClasses.xml");

# Augment function class list with required information.
foreach my $functionClass ( &List::ExtraUtils::hashList($libraryFunctionClasses->{'classes'}, keyAs => "name") ) {
    $functionClass->{'module'} = $stateStorables->{'functionClasses'}->{$functionClass->{'name'}."Class"}->{'module'};
}

# Iterate over all files containing function class definitions.
foreach my $fileName ( @{$directiveLocations->{'functionClass'}->{'file'}} ) {
    my $tree   = &Galacticus::Build::SourceTree::ParseFile($fileName);
    my $node_  = $tree;
    my $depth  = 0;
    while ( $node_ ) {
	my $node = $node_;
	$node_ = &Galacticus::Build::SourceTree::Walk_Tree($node_,\$depth);
	# Skip anything that isn't a functionClass node.
	next
	    unless ( $node->{'type'} eq "functionClass" );
	next
	    unless ( grep {$node->{'directive'}->{'name'} eq $_} keys(%{$libraryFunctionClasses->{'classes'}}) );
	# Get the corresponding functionClass.
	my $functionClass = $libraryFunctionClasses->{'classes'}->{$node->{'directive'}->{'name'}};
	if ( exists($node->{'directive'}->{'method'}->{'name'}) && ! reftype($node->{'directive'}->{'method'}->{'name'}) ) {
	    # Only a single method is defined for this class. Turn it back into a hash reference.
	    $functionClass->{'methods'}->{$node->{'directive'}->{'method'}->{'name'}} = $node->{'directive'}->{'method'};	
	} else {
	    $functionClass->{'methods'}                                               = $node->{'directive'}->{'method'};	
	}
	# Find all implementations of this class.
	my $classID = 0;
	foreach my $fileNameImplementation ( &List::ExtraUtils::as_array($directiveLocations->{$functionClass->{'name'}}->{'file'}) ) {
	    my $treeImplementation   = &Galacticus::Build::SourceTree::ParseFile($fileNameImplementation);
	    my $nodeImplementation_  = $treeImplementation;
	    my $depthImplementation  = 0;
	    my $nameImplementation;
	    my $nameConstructor;
	    my @argumentsConstructor;
	    while ( $nodeImplementation_ ) {
		my $nodeImplementation = $nodeImplementation_;
		$nodeImplementation_ = &Galacticus::Build::SourceTree::Walk_Tree($nodeImplementation_,\$depthImplementation);
		if ( $nodeImplementation->{'type'} eq $functionClass->{'name'} ) {
		    # Implementation node found here.
		    $nameImplementation = $nodeImplementation->{'directive'}->{'name'};
		} elsif ( defined($nameImplementation) && $nodeImplementation->{'type'} eq "interface" && $nodeImplementation->{'name'} eq $nameImplementation ) {
		    # Constructor interface found here.
		    my $nodeInterface = $nodeImplementation->{'firstChild'};
		    while ( $nodeInterface ) {
			if ( $nodeInterface->{'type'} eq "moduleProcedure" ) {
			    my @constructorsInternal = grep {$_ =~ m/Internal$/} @{$nodeInterface->{'names'}};
			    $nameConstructor = $constructorsInternal[0]
				if ( scalar(@constructorsInternal) == 1 );
			}
			$nodeInterface = $nodeInterface->{'sibling'};
		    }
		} elsif ( defined($nameConstructor) && $nodeImplementation->{'type'} eq "function" && $nodeImplementation->{'name'} eq $nameConstructor ) {
		    # Internal constructor found. Extract declarations.
		    if ( $nodeImplementation->{'opener'} =~ m/^\s*function\s+$nameConstructor\s*\(([^\)]+)\)/ ) {
			@argumentsConstructor = map {{name => $_}} split(/\s*,\s*/,$1);
		    }		    
		    my $nodeConstructor = $nodeImplementation->{'firstChild'};
		    while ( $nodeConstructor ) {
			if ( $nodeConstructor->{'type'} eq "declaration" ) {
			    foreach my $declaration ( @{$nodeConstructor->{'declarations'}} ) {
				foreach my $variable ( @{$declaration->{'variables'}} ) {
				    foreach my $argument ( @argumentsConstructor ) {
					if ( lc($argument->{'name'}) eq $variable ) {
					    $argument->{$_} = $declaration->{$_}
					       foreach ( "intrinsic", "type", "attributes" );
					}
				    }
				}
			    }
			}
			$nodeConstructor = $nodeConstructor->{'sibling'};
		    }

		}
	    }
	    die("Unable to find implementation of '".$functionClass->{'name'}."' in '".$fileNameImplementation."'")
		unless ( defined($nameImplementation) );
	    # If no internal constructor was found, use the default constructor.
	    $nameConstructor      = $nameImplementation
		unless ( defined($nameConstructor   ) );
	    ++$classID;
	    my $implementation =
	    {
		name      => $nameImplementation    ,
		classID   => $classID               ,
		fileName  => $fileNameImplementation,
		arguments => @argumentsConstructor ? \@argumentsConstructor : []
	    };
	    push(
		@{$functionClass->{'implementations'}},
		$implementation
		);
	}
	# Add Python parent class.
	&interfacesPythonClasses(      $python,$functionClass                        );	
	# Add pointer get functions.
	&interfacesPointerGet   ($code        ,$functionClass                        );
	# Add constructors.
	&interfacesConstructors ($code,$python,$functionClass,$libraryFunctionClasses);
	# Add interfaces to all methods.
	&interfacesMethods      ($code,$python,$functionClass                        );
	# Add a destructor.
	&interfacesDestructor   ($code,$python,$functionClass                        );	
    }
}

# Append the initialization code. This is necessary to cause libgfortran to be initialized at run time.
my $libraryInitializer = fill_in_string(<<'CODE', PACKAGE => 'code');
program libGalacticusInit
end program libGalacticusInit
CODE
push(
    @{$code->{'units'}},
    $libraryInitializer
    );

# Serialize the code.
open(my $output,">",$ENV{'BUILDPATH'}."/libgalacticus.Inc");
print $output join("\n",@{$code->{'units'}})."\n";
close($output);

# Generate code to initialize the Python interface,
$python->{'units'}->{'init'}->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'ext');
from ctypes import *
# Load the shared library into ctypes.
libname = "./libgalacticus.so"
c_lib = CDLL(libname)
CODE
foreach my $clibFunction ( @{$python->{'c_lib'}} ) {
    $python->{'units'}->{'init'}->{'content'} .= "c_lib.".$clibFunction->{'name'}.".restype  = "  .            $clibFunction->{'restype' }  ."\n"
	if ( defined($clibFunction->{'restype' }) );
    $python->{'units'}->{'init'}->{'content'} .= "c_lib.".$clibFunction->{'name'}.".argtypes = [ ".join(", ",@{$clibFunction->{'argtypes'}})." ]\n"
	if ( defined($clibFunction->{'argtypes'}) );
}
$python->{'units'}->{'init'}->{'indent'} = 0;

# Serialize the Python interfaces.
my %dependencies;
foreach my $unitName ( keys(%{$python->{'units'}}) ) {
    next
	unless ( exists($python->{'units'}->{$unitName}->{'dependencies'}) );
    push(@{$dependencies{$unitName}},@{$python->{'units'}->{$unitName}->{'dependencies'}});
}
my @pythonNamesUnordered = sort(keys(%{$python->{'units'}}));
my @pythonNamesOrdered   = &Sort::Topo::sort(\@pythonNamesUnordered,\%dependencies);
my @pythonStackOrdered   = map {$python->{'units'}->{$_}} @pythonNamesOrdered;
open(my $pythonOutput,">python/galacticus.py");
while ( scalar(@pythonStackOrdered) > 0 ) {
    my $unit = pop(@pythonStackOrdered);
    push(@pythonStackOrdered,map {$_->{'indent'} = $unit->{'indent'}+1; $_} @{$unit->{'subUnits'}})
	if ( exists($unit->{'subUnits'}) );
    my $indent = "    " x $unit->{'indent'};
    open(my $code,"<",\$unit->{'content'});
    while ( my $codeLine = <$code> ) {
	print $pythonOutput $indent.$codeLine;
    }
    close($code);
    print $pythonOutput "\n";
}

exit;

sub interfacesPointerGet {
    # Build functions which return pointers to the specific type.
    my $code                   = shift();
    $ext::functionClass        = shift();
    @ext::functionClassSymbols = ( $ext::functionClass->{'name'}."Class" );
    push(@ext::functionClassSymbols,map {$_->{'name'}} @{$ext::functionClass->{'implementations'}});
    my $function = fill_in_string(<<'CODE', PACKAGE => 'ext');
function {$functionClass->{'name'}}GetPtr({$functionClass->{'name'}}_,classID)
  use, intrinsic :: ISO_C_Binding               , only : c_ptr                          , c_int, c_f_pointer
  use            :: Galacticus_Error            , only : Galacticus_Error_Report
  use            :: {$functionClass->{'module'}}, only : {join(", ",@functionClassSymbols)}
  implicit none
  class({$functionClass->{'name'}}Class), pointer :: {$functionClass->{'name'}}GetPtr
  type(c_ptr), intent(in   ) :: {$functionClass->{'name'}}_
  integer(c_int), intent(in   ) :: classID
{join("\n",map {"  type(".$_->{'name'}."), pointer :: ".$_->{'name'}."_"} @{$functionClass->{'implementations'}})}

  select case (classID)
CODE
    foreach $ext::implementation ( @{$ext::functionClass->{'implementations'}} ) {
	$function .= fill_in_string(<<'CODE', PACKAGE => 'ext');
  case ({$implementation->{'classID'}})
     call c_f_pointer({$functionClass->{'name'}}_,{$implementation->{'name'}}_)
     {$functionClass->{'name'}}GetPtr => {$implementation->{'name'}}_
CODE
    }
    $function .= fill_in_string(<<'CODE', PACKAGE => 'ext');
  case default
     {$functionClass->{'name'}}GetPtr => null()
     call Galacticus_Error_Report('unknown classID'//\{introspection:location\})
  end select
  return
end function {$functionClass->{'name'}}GetPtr
CODE
    push(
	@{$code->{'units'}},
	$function
	);
}

sub interfacesConstructors {
    # Build interfaces to constructors.
    my $code                   = shift();
    my $python                 = shift();
    $ext::functionClass        = shift();
    my $libraryFunctionClasses = shift();
    foreach $ext::implementation ( @{$ext::functionClass->{'implementations'}} ) {
	# Reset all generated code fragments.
	undef($ext::declarations              );
	undef($ext::constructorArguments      );
	undef(@ext::pythonConstructorArguments);
	undef($ext::dereference               );
	undef(@ext::functionArguments         );
	undef(@ext::pythonFunctionArguments   );
	undef(%ext::isoCBindingSymbols        );
	undef($ext::moduleUses                );
	my @clibArgTypes;
	# Initialize the hash of ISO_C_Binding symbols that we must import. "c_ptr" and "c_loc" are always needed.
	$ext::isoCBindingSymbols{'c_ptr'} = 1;
	$ext::isoCBindingSymbols{'c_loc'} = 1;
	# Initialize a hash of which "GetPtr" functions have already been declared.
	my %getPtrClasses;
	# Add the "self" argument to the Python constructor arguments.
	push(@ext::pythonConstructorArguments,"self");
	# Iterate over each argument, building the argument list, declarations, pointer dereferencings, and constructor arguments.
	foreach my $argument ( @{$ext::implementation->{'arguments'}} ) {
	    # Flag if this is the last argument.
	    my $isLast = $argument == ${$ext::implementation->{'arguments'}}[-1];
	    # Push onto the function argument list.
	    push(@ext::functionArguments      ,$argument->{'name'});
	    push(@ext::pythonFunctionArguments,$argument->{'name'});
	    # Build the C-interoperable declaration.
	    my $cType;
	    my @cAttributes;
	    # By default the argument name to pass to the constructor is just the argument name passed to this interface function.
	    $argument->{'passName'} = $argument->{'name'};
	    # Map types.
	    my $declarations = "";
	    if ( $argument->{'intrinsic'} eq "double precision" ) {
		$cType                               = "real(c_double)";
		$ext::isoCBindingSymbols{'c_double'} = 1;
		push(@clibArgTypes,"c_double");
	    } elsif ( $argument->{'intrinsic'} eq "integer" ) {
		unless ( defined($argument->{'type'}) ) {
		    $cType                            = "integer(c_int)";
		    $ext::isoCBindingSymbols{'c_int'} = 1;
		    push(@clibArgTypes,"c_int");
		}
	    } elsif ( $argument->{'intrinsic'} eq "logical" ) {
		$cType                               = "logical(c_bool)";
		$ext::isoCBindingSymbols{'c_bool'} = 1;
		push(@clibArgTypes,"c_bool");
		$ext::dereference          .= $argument->{'name'}."_=logical(".$argument->{'name'}.")\n";
		$argument->{'passName'}    .= "_";
		$declarations              .= "  logical :: ".$argument->{'name'}."_\n";
	    } elsif ( $argument->{'intrinsic'} eq "character" ) {
		$cType                               = "character(c_char)";
		$ext::isoCBindingSymbols{'c_char'} = 1;
		push(@clibArgTypes,"c_char_p");
		push(@cAttributes,"dimension(*)");
		$argument->{'passName'}     = "char(String_C_to_Fortran(".$argument->{'name'}."))";
		$ext::moduleUses->{'String_Handling'}->{'String_C_to_Fortran'} = 1;
		$ext::moduleUses->{'ISO_Varying_String'}->{'char'} = 1;
	    } elsif ( $argument->{'intrinsic'} eq "class" ) {
		$cType                            = "type(c_ptr)";
		$ext::isoCBindingSymbols{'c_ptr'} = 1;
		push(@clibArgTypes,"c_void_p");
		# Check for a functionClass argument.
		if ( my @matchedClass = grep {$argument->{'type'} eq $_."Class"} keys(%{$libraryFunctionClasses->{'classes'}}) ) {
		    # Add a classID argument for this functionClass argument, and send the object pointer and ID separately to the Python interface.
		    push(@ext::functionArguments,$argument->{'name'}."_ID");
		    push(@clibArgTypes,"c_int");
		    $ext::pythonFunctionArguments[-1] .= "._glcObj";
		    push(@ext::pythonFunctionArguments,$argument->{'name'}."._classID");
		    $declarations                     .= "  integer(c_int), value :: ".$argument->{'name'}."_ID\n";
		    $ext::isoCBindingSymbols{'c_int'}  = 1;
		    # Add code to deference this to a pointer.
		    (my $className = $argument->{'type'}) =~ s/Class$//;
		    $declarations              .= "  class(".$argument->{'type'}."), pointer :: ".$className."GetPtr\n"
			unless ( exists($getPtrClasses{$className}) );
		    $declarations              .= "  class(".$argument->{'type'}."), pointer :: ".$argument->{'name'}."_\n";
		    $getPtrClasses{$className}  = 1;
		    $ext::dereference          .= $argument->{'name'}."_ => ".$className."GetPtr(".$argument->{'name'}.",".$argument->{'name'}."_ID)\n";
		    $argument->{'passName'}    .= "_";
		    my $moduleName              = $libraryFunctionClasses->{'classes'}->{$matchedClass[0]}->{'module'};
		    $ext::moduleUses->{$moduleName}->{$argument->{'type'}} = 1;
		}
	    } else {
		die "unsupported type '".$argument->{'intrinsic'}."' in constructor arguemnt '".$argument->{'name'}."' for class '".$ext::implementation->{'name'}."'";
	    }
	    # Determine intent.
	    if ( ( $cType eq "type(c_ptr)" || grep {$_ =~ m/intent\s*\(\s*in\s*\)/} @{$argument->{'attributes'}} ) && ! grep {$_ =~ m/^dimension/} @cAttributes ) {
		push(@cAttributes,"value");
	    } elsif ( grep {$_ =~ m/^dimension/} @cAttributes ) {
		# Array object - no need to add "value" attribute.
	    } else {
		die("unsupported intent '".join(", ",@{$argument->{'attributes'}})."'");
	    }
	    # Add this argument to the Python constructor.
	    my $pythonConstructorArgument = $argument->{'name'};
	    $pythonConstructorArgument .= "=None"
		if ( grep {$_ eq "optional"} @{$argument->{'attributes'}} );
	    push(@ext::pythonConstructorArguments,$pythonConstructorArgument);
	    # Add a declaration for this argument, and add it into the call to the constructor function.
	    $ext::declarations         .= "  ".$cType.(@cAttributes ? ", " : "").join(", ",@cAttributes)." :: ".$argument->{'name'}."\n".$declarations;
	    $ext::constructorArguments .= "  &amp; ".$argument->{'name'}."=".$argument->{'passName'}.($isLast ? " &amp;" : ", &amp;\n");
	}
	# Build code for module uses.
	$ext::moduleUseCode = "";
	foreach my $moduleName ( sort(keys(%{$ext::moduleUses})) ) {
	    $ext::moduleUseCode .= "use :: ".$moduleName.", only : ".join(", ",sort(keys(%{$ext::moduleUses->{$moduleName}})))."\n";
	}
	# Finally, build the interface constructor function.
	my $constructor .= fill_in_string(<<'CODE', PACKAGE => 'ext');
function {$implementation->{'name'}}PY({join(",",@functionArguments)}) bind(c,name='{$implementation->{'name'}}PY')
  use, intrinsic :: ISO_C_Binding               , only : {join(", ",keys(%isoCBindingSymbols))}
  use            :: {$functionClass->{'module'}}, only : {$implementation->{'name'}}
{$moduleUseCode}
  implicit none
  type(c_ptr                      )          :: {$implementation->{'name'}}PY
  type({$implementation->{'name'}}), pointer :: self
{$declarations}

{$dereference}
  allocate(self)
  !![
  <referenceConstruct object="self">
   <constructor>
    {$implementation->{'name'}}( &amp;
{$constructorArguments}
     &amp;                     )
   </constructor>
  </referenceConstruct>
  !!]
  {$implementation->{'name'}}PY=c_loc(self)
  return
end function {$implementation->{'name'}}PY
CODE
	push(
	    @{$code->{'units'}},
	    $constructor
	);
	# Add library interface descriptors.
	push(
	    @{$python->{'c_lib'}},
	    {
		name     => $ext::implementation->{'name'}."PY",
		restype  => "c_void_p",
		argtypes => \@clibArgTypes
	    }
	    );
    # Add a constructor to the Python class.
    my $pythonConstructor = fill_in_string(<<'CODE', PACKAGE => 'ext');
# Constructor
def __init__({join(",",@pythonConstructorArguments)}):
    # Assign class ID so relevant pointers can be constructed on the Fortran side.
    self._classID = {$implementation->{'classID'}}
    self._glcObj = c_lib.{$implementation->{'name'}}PY({join(",",@pythonFunctionArguments)})
CODE
    push(
	@{$python->{'units'}->{$ext::implementation->{'name'}}->{'subUnits'}},
	{
	    content => $pythonConstructor
	}
	);
    }
}

sub interfacesDestructor {
    # Build interfaces to destructors.
    my $code            = shift();
    my $python          = shift();
    $ext::functionClass = shift();
    my $destructor = fill_in_string(<<'CODE', PACKAGE => 'ext');
subroutine {$functionClass->{'name'}}DestructorPY(self,classID) bind(c,name='{$functionClass->{'name'}}DestructorPY')
  use, intrinsic :: ISO_C_Binding               , only : c_ptr                          , c_int
  use            :: {$functionClass->{'module'}}, only : {$functionClass->{'name'}}Class
  implicit none
  type   (c_ptr                          ), value  , intent(in   ) :: self
  integer(c_int                          ), value  , intent(in   ) :: classID
  class  ({$functionClass->{'name'}}Class), pointer                :: self_  , {$functionClass->{'name'}}GetPtr

  self_ => {$functionClass->{'name'}}GetPtr(self,classID)
  !![
  <objectDestructor name="self_"/>
  !!]
  return
end subroutine {$functionClass->{'name'}}DestructorPY
CODE
    push(
	@{$code->{'units'}},
	$destructor
	);
    # Add c_lib interface.
    push(
	@{$python->{'c_lib'}},
	{
	    name     => $ext::functionClass->{'name'}."DestructorPY",
	    restype  => undef(),
	    argtypes => [ "c_void_p", "c_int" ]
	}
	);
    # Add a destructor to the Python class.
    my $pythonDestructor = fill_in_string(<<'CODE', PACKAGE => 'ext');
# Destructor
def __del__(self):
    c_lib.{$functionClass->{'name'}}DestructorPY(self._glcObj,self._classID)
CODE
    push(
	@{$python->{'units'}->{$ext::functionClass->{'name'}}->{'subUnits'}},
	{
	    content => $pythonDestructor
	}
	);
}

sub interfacesMethods {
    # Build interface functions to class methods.
    my $code            = shift();
    my $python          = shift();
    $ext::functionClass = shift();
    foreach $ext::method ( &List::ExtraUtils::hashList($ext::functionClass->{'methods'},keyAs => 'name') ) {
	# Reset all generated code fragments.
	undef($ext::procedure               );
	undef($ext::declarations            );
	undef(@ext::interfaceArguments      );
	undef(@ext::pythonInterfaceArguments);
	undef($ext::reassignments           );
	undef(@ext::methodNames             );
	undef(@ext::methodArguments         );
	undef(@ext::pythonMethodArguments   );
	undef(%ext::isoCBindingSymbols      );
	undef($ext::resultConversionOpen    );
	undef($ext::resultConversionClose   );
	undef(@ext::clibArgTypes            );
	my $clibResType;
	# Initialize the hash of ISO_C_Binding symbols that we must import. "c_ptr" and "c_int" are always needed.
	$ext::isoCBindingSymbols{'c_ptr'} = 1;
	$ext::isoCBindingSymbols{'c_int'} = 1;
	# Parse arguments.
	my @arguments;
	if ( exists($ext::method->{'argument'}) ) {
	    foreach my $argument ( &List::ExtraUtils::as_array($ext::method->{'argument'}) ) {
		push(@arguments,&Galacticus::Build::SourceTree::Parse::Declarations::parseDeclaration($argument));
	    }
	}
	# Add function declaration.
	$ext::procedure             = "function";
	$ext::resultConversionOpen  = "";
	$ext::resultConversionClose = "";
	my $functionCType;
	if ( $ext::method->{'type'} eq "double precision" ) {
	    $functionCType                       = "real(c_double)";
	    $clibResType                         = "c_double";
	    $ext::isoCBindingSymbols{'c_double'} = 1;
	} elsif ( $ext::method->{'type'} eq "logical" ) {
	    $functionCType                       = "logical(c_bool)";
	    $clibResType                         = "c_bool";
	    $ext::isoCBindingSymbols{'c_bool'}   = 1;
	    $ext::resultConversionOpen           = "logical(";
	    $ext::resultConversionClose          = ",kind=c_bool)";
	} elsif ( $ext::method->{'type'} eq "void" ) {
	    	$ext::procedure = "subroutine";
	} else {
	    die("unsupported type '".$ext::method->{'type'}."'");
	}
	$ext::declarations  .= $functionCType." :: ".$ext::functionClass->{'name'}.ucfirst($ext::method->{'name'})."PY\n"
	    unless ( $ext::method->{'type'} eq "void" );
	# Add self argument, declaration, and dereference.
	push(@ext::interfaceArguments      ,"self"         ,"selfClassID"  );
	push(@ext::pythonInterfaceArguments,"self._glcObj" ,"self._classID");
	push(@ext::pythonMethodArguments   ,"self"                         );
	$ext::declarations  .= "type   (c_ptr                                 ), value  , intent(in   ) :: self\n";
	$ext::declarations  .= "integer(c_int                                 ), value  , intent(in   ) :: selfClassID\n";
	push(@ext::clibArgTypes,"c_void_p","c_int");
	$ext::declarations  .= "class  (".$ext::functionClass->{'name'}."Class), pointer                :: self_, ".$ext::functionClass->{'name'}."GetPtr\n";
	$ext::reassignments .= "self_ => ".$ext::functionClass->{'name'}."GetPtr(self,selfClassID)\n";
	# Iterate over arguments, adding them to the lists, adding declarations, and any reassignments.
	$ext::countOptional = 0;
	foreach my $argument ( @arguments ) {
	    push(@ext::interfaceArguments      ,@{$argument->{'variableNames'}});
	    push(@ext::pythonInterfaceArguments,@{$argument->{'variableNames'}});
	    my $cType;
	    my @cAttributes;
	    # Determine pass by. Scalar arguments that are 'intent(in   )' and not optional can be passed by value, others must be passed by reference.
	    my $passBy =
		(
		 (  grep {$_ =~ m/intent\s*\(\s*in\s*\)/} @{$argument->{'attributes'}})
		 &&
		 (! grep {$_ =~ m/dimension/            } @{$argument->{'attributes'}})
		 &&
		 (! grep {$_ eq "optional"              } @{$argument->{'attributes'}})
		)
		?
		"value"
		:
		"ref";
	    # Determine if argument is optional.
	    my $isOptional = grep {$_ eq "optional"} @{$argument->{'attributes'}};
	    $ext::countOptional += scalar(@{$argument->{'variableNames'}})
		if ( $isOptional );
	    # Generate code for this argument.
	    if ( $argument->{'intrinsic'} eq "double precision" ) {
		$cType                                = "real(c_double)";
		$ext::isoCBindingSymbols{'c_double'}  = 1;
		push(@ext::methodArguments,@{$argument->{'variableNames'}});
		push(@ext::methodNames    ,@{$argument->{'variableNames'}});
		my $suffix = $isOptional ? "=None" : "";
		push(@ext::pythonMethodArguments,map {$_.$suffix} @{$argument->{'variableNames'}});
		push(@ext::clibArgTypes,("c_double") x scalar(@{$argument->{'variableNames'}}));
	    } elsif ( $argument->{'intrinsic'} eq "integer" ) {
		if ( defined($argument->{'type'}) ) {
		    die("unsupported integer type '".$argument->{'type'}."'");
		} else {
		    # This is a standard integer, corresponding to a c_int.
		    $cType                             = "integer(c_int)";
		    $ext::isoCBindingSymbols{'c_int'}  = 1;
		    push(@ext::methodArguments,@{$argument->{'variableNames'}});   
		    push(@ext::methodNames    ,@{$argument->{'variableNames'}});   
		    my $suffix = $isOptional ? "=None" : "";
		    push(@ext::pythonMethodArguments,map {$_.$suffix} @{$argument->{'variableNames'}});
		    push(@ext::clibArgTypes,("c_int") x scalar(@{$argument->{'variableNames'}}));
		}
	    } elsif ( $argument->{'intrinsic'} eq "logical" ) {
		$cType                                = "logical(c_bool)";
		$ext::isoCBindingSymbols{'c_bool'  }  = 1;
		$ext::declarations                   .= "logical :: ".join(", ",map {$_."_"} @{$argument->{'variableNames'}})."\n";
		push(@ext::methodArguments,map {$_."_"} @{$argument->{'variableNames'}});
		push(@ext::methodNames    ,             @{$argument->{'variableNames'}});
		my $suffix = $isOptional ? "=None" : "";
		push(@ext::pythonMethodArguments,map {$_.$suffix} @{$argument->{'variableNames'}});
		push(@ext::clibArgTypes,("c_bool") x scalar(@{$argument->{'variableNames'}}));
		foreach my $variableName ( @{$argument->{'variableNames'}} ) {
		    $ext::reassignments              .=  ($isOptional ? "if (present(".$variableName.")) " : "").$variableName."_=".$variableName."\n";
		}
	    } else {
		print "unsupported type '".$argument->{'intrinsic'}."'\n";
	    }
	    # Add attributes for arguments.
	    push(@cAttributes,"value")
		if ( $passBy eq "value" );
	    push(@cAttributes,"optional")
		if ( $isOptional );
	    $ext::declarations .= $cType.(@cAttributes ? ", " : "").join(", ",@cAttributes)." :: ".join(", ",@{$argument->{'variableNames'}})."\n";
	    
	}
	# Generate the function.
	my $function = fill_in_string(<<'CODE', PACKAGE => 'ext');
{$procedure} {$functionClass->{'name'}}{ucfirst($method->{'name'})}PY({join(",",@interfaceArguments)}) bind(c,name='{$functionClass->{'name'}}{ucfirst($method->{'name'})}PY')
  use, intrinsic :: ISO_C_Binding               , only : {join(", ",keys(%isoCBindingSymbols))}
  use            :: {$functionClass->{'module'}}, only : {$functionClass->{'name'}}Class
  implicit none
{$declarations}

{$reassignments}
CODE
	if ( $ext::countOptional > 0 ) {
	    my $formatBinary = "%.".$ext::countOptional."b";
	    for(my $i=0;$i<2**$ext::countOptional;++$i) {
		@ext::state = split(//,sprintf($formatBinary,$i));
		$ext::condition = ($i > 0 ? "else " : "")."if (".join(" .and. ",map {($ext::state[$_-1] == 0 ? ".not." : "")."present(".$ext::methodNames[-$ext::countOptional-1+$_].")"} 1..$ext::countOptional).") then";
		$function .= fill_in_string(<<'CODE', PACKAGE => 'ext');
{$condition}
    {$method->{'type'} eq "void" ? "call " : $functionClass->{'name'}.ucfirst($method->{'name'})."PY="}{$resultConversionOpen}self_%{$method->{'name'}}({join(",",map {$_ < scalar(@methodArguments)-$countOptional ? $methodArguments[$_] : ($state[$_-scalar(@methodArguments)+$countOptional] == 1 ? $methodNames[$_]."=".$methodArguments[$_] : ())} 0..scalar(@methodArguments)-1)}){$resultConversionClose}
CODE
	    }
	    $function .= "else\n";
        }
	$function .= fill_in_string(<<'CODE', PACKAGE => 'ext');
    {$method->{'type'} eq "void" ? "call " : $functionClass->{'name'}.ucfirst($method->{'name'})."PY="}{$resultConversionOpen}self_%{$method->{'name'}}({join(",",@methodArguments)}){$resultConversionClose}
CODE
	if ( $ext::countOptional > 0 ) {
	    $function .= "end if\n";
        }
	$function .= fill_in_string(<<'CODE', PACKAGE => 'ext');
  return
end {$procedure} {$functionClass->{'name'}}{ucfirst($method->{'name'})}PY
CODE
	push(
	    @{$code->{'units'}},
	    $function
	    );
	# Add library interface descriptors.
	my @clibArgTypesNonOptional = @ext::clibArgTypes[0..scalar(@ext::clibArgTypes)-1-$ext::countOptional];
	push(
	    @{$python->{'c_lib'}},
	    {
		name     => $ext::functionClass->{'name'}.ucfirst($ext::method->{'name'})."PY",
		restype  => $clibResType,
		argtypes => \@clibArgTypesNonOptional
	    }
	    );
	# Add a method to the Python class.
	my $pythonMethod = fill_in_string(<<'CODE', PACKAGE => 'ext');
def {$method->{'name'}}({join(",",@pythonMethodArguments)}):
CODE
	if ( $ext::countOptional == 0 ) {
	    $pythonMethod .= fill_in_string(<<'CODE', PACKAGE => 'ext');
    return c_lib.{$functionClass->{'name'}}{ucfirst($method->{'name'})}PY({join(",",@pythonInterfaceArguments)})
CODE
	} else {
	    my $formatBinary = "%.".$ext::countOptional."b";
	    for(my $i=0;$i<2**$ext::countOptional;++$i) {
		@ext::state = split(//,sprintf($formatBinary,$i));
		$ext::condition = "    ".($i > 0 ? "el" : "")."if ".join(" and ",map {($ext::state[$_-1] == 0 ? "not " : "")."".$ext::pythonInterfaceArguments[-$ext::countOptional-1+$_]} 1..$ext::countOptional).":";
		$pythonMethod .= fill_in_string(<<'CODE', PACKAGE => 'ext');
{$condition}
        return c_lib.{$functionClass->{'name'}}{ucfirst($method->{'name'})}PY({join(",",map {$_ < scalar(@pythonInterfaceArguments)-$countOptional ? $pythonInterfaceArguments[$_] : ($state[$_-scalar(@pythonInterfaceArguments)+$countOptional] == 1 ? "byref(".$clibArgTypes[$_]."(".$pythonInterfaceArguments[$_]."))" : "None")} 0..scalar(@pythonInterfaceArguments)-1)})
CODE
	    }
	}
	push(
	    @{$python->{'units'}->{$ext::functionClass->{'name'}}->{'subUnits'}},
	    {
		content => $pythonMethod
	    }
	    );	
    }
}

sub interfacesPythonClasses {
    # Build Python parent class.
    my $python          = shift();
    $ext::functionClass = shift();

    # Add the parent class.
    my $class = fill_in_string(<<'CODE', PACKAGE => 'ext');
class {$functionClass->{'name'}}:

    # Constructor
    def __init__(self):
        # Assign class ID to a negative number - indicating this is not a concrete class.
        self._classID = -1
CODE
    $python->{'units'}->{$ext::functionClass->{'name'}}->{'content'     } = $class;
    $python->{'units'}->{$ext::functionClass->{'name'}}->{'indent'      } = 0;
    $python->{'units'}->{$ext::functionClass->{'name'}}->{'dependencies'} = [ "init" ];

    # Add child classes.
    foreach $ext::implementation ( @{$ext::functionClass->{'implementations'}} ) {
	my $class = fill_in_string(<<'CODE', PACKAGE => 'ext');
class {$implementation->{'name'}}({$functionClass->{'name'}}):
CODE
	$python->{'units'}->{$ext::implementation->{'name'}}->{'content'     } = $class;
	$python->{'units'}->{$ext::implementation->{'name'}}->{'indent'      } = 0;
	$python->{'units'}->{$ext::implementation->{'name'}}->{'dependencies'} = [ $ext::functionClass->{'name'} ];	
    }
}
