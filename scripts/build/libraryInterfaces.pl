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
$code->{'main'} = [];

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
    my $tree       = &Galacticus::Build::SourceTree::ParseFile($fileName);
    my $node_      = $tree;
    my $depth      = 0;
    my @moduleUses;
    while ( $node_ ) {
	my $node = $node_;
	$node_ = &Galacticus::Build::SourceTree::Walk_Tree($node_,\$depth);
	# Extract any module uses.
	push(@moduleUses,$node->{'moduleUse'})
	    if ( $node->{'type'} eq "moduleUse" );
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
	$functionClass->{'moduleUses'} = \@moduleUses;
	# Find the parent-child class dependencies.
	my $extensions;
	my $moduleUsesImplementations;
	foreach my $fileNameImplementation ( &List::ExtraUtils::as_array($directiveLocations->{$functionClass->{'name'}}->{'file'}) ) {
	    my $treeImplementation     = &Galacticus::Build::SourceTree::ParseFile($fileNameImplementation);
	    my $depthImplementation    = 0;
	    my $nodeImplementation     = $treeImplementation;
	    my $nameImplementation;
	    my @moduleUsesImplementation;
	    while ( $nodeImplementation ) {
		if ( $nodeImplementation->{'type'} eq $functionClass->{'name'} ) {
		    # Find the name of the implementation.
		    $nameImplementation = $nodeImplementation->{'directive'}->{'name'};
		} elsif ( $nodeImplementation->{'type'} eq "type" && $nodeImplementation->{'name'} eq $nameImplementation ) {
		    # Find the extended class.
		    if ( $nodeImplementation->{'opener'} =~ m/,\s*extends\s*\(\s*([a-zA-Z0-9_]+)\s*\)/ ) {
			$extensions->{$nodeImplementation->{'name'}} = $1;
		    }
		} elsif ( $nodeImplementation->{'type'} eq "moduleUse" ) {
		    # Module use statements, extract for later reference.
		    push(@moduleUsesImplementation,$nodeImplementation->{'moduleUse'});
		}
		
		$nodeImplementation = &Galacticus::Build::SourceTree::Walk_Tree($nodeImplementation,\$depthImplementation);
	    }
	    @{$moduleUsesImplementations->{$nameImplementation}} = @moduleUsesImplementation
		unless (
		    exists($functionClass->{$nameImplementation}                      )
		    &&
		    exists($functionClass->{$nameImplementation}->{'exclude'}         )
		    &&
		           $functionClass->{$nameImplementation}->{'exclude'} eq "yes"
		);
	}	
	# Find all implementations of this class.
	my $classID = 0;
	foreach my $fileNameImplementation ( &List::ExtraUtils::as_array($directiveLocations->{$functionClass->{'name'}}->{'file'}) ) {
	    my $treeImplementation     = &Galacticus::Build::SourceTree::ParseFile($fileNameImplementation);
	    my $nodeImplementation_    = $treeImplementation;
	    my $depthImplementation    = 0;
	    my $abstractImplementation;
	    my $nameImplementation;
	    my $nameConstructor;
	    my @argumentsConstructor;
	    my @moduleUsesImplementation;
	    while ( $nodeImplementation_ ) {
		my $nodeImplementation = $nodeImplementation_;
		$nodeImplementation_ = &Galacticus::Build::SourceTree::Walk_Tree($nodeImplementation_,\$depthImplementation);
		if ( $nodeImplementation->{'type'} eq $functionClass->{'name'} ) {
		    # Implementation node found here.
		    $nameImplementation = $nodeImplementation->{'directive'}->{'name'};
		    $abstractImplementation = exists($nodeImplementation->{'directive'}->{'abstract'}) && $nodeImplementation->{'directive'}->{'abstract'} eq "yes";
		} elsif ( defined($nameImplementation) && $nodeImplementation->{'type'} eq "interface" && defined($nodeImplementation->{'name'}) && $nodeImplementation->{'name'} eq $nameImplementation ) {
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
		    if ( $nodeImplementation->{'opener'} =~ m/^\s*(recursive\s+)*function\s+$nameConstructor\s*\(([^\)]+)\)/ ) {
			@argumentsConstructor = map {{name => $_}} split(/\s*,\s*/,$2);
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
		} elsif ( $nodeImplementation->{'type'} eq "moduleUse" ) {
		    # Module use statements, extract for later reference.
		    push(@moduleUsesImplementation,$nodeImplementation->{'moduleUse'});
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
		name       => $nameImplementation       ,
		classID    => $classID                  ,
		fileName   => $fileNameImplementation   ,
		moduleUses => \@moduleUsesImplementation,
		arguments  => @argumentsConstructor ? \@argumentsConstructor : []
	    };
	    push(
		@{$functionClass->{'implementations'}},
		$implementation
		)
		unless (
		    $abstractImplementation
		    ||
		    (
		     exists($functionClass->{$nameImplementation}                      )
		     &&
		     exists($functionClass->{$nameImplementation}->{'exclude'}         )
		     &&
		            $functionClass->{$nameImplementation}->{'exclude'} eq "yes"
		    )
		);
	}
	# Add Python parent class.
	&interfacesPythonClasses(      $python,$functionClass                                                               );	
	# Add pointer get functions.
	&interfacesPointerGet   ($code        ,$functionClass                                                               );
	# Add constructors.
	&interfacesConstructors ($code,$python,$functionClass,$libraryFunctionClasses,$extensions,$moduleUsesImplementations);
	# Add interfaces to all methods.
	&interfacesMethods      ($code,$python,$functionClass                        ,$extensions,$moduleUsesImplementations);
	# Add a destructor.
	&interfacesDestructor   ($code,$python,$functionClass                                                               );
    }
}

# Append the initialization code. This consists of a function that will be called on module import to initialize various things in
# Galacticus, plus a Fortran "program" unit. This latter is necessary to cause libgfortran to be initialized at run time.
my $libraryInitializer = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine libGalacticusInitL() bind(c,name='libGalacticusInitL')
  use:: Events_Hooks, only : eventsHooksInitialize
  use :: IO_HDF5    , only : ioHDF5AccessInitialize

  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Initialize HDF5 library acess lock.
  call ioHDF5AccessInitialize()
end subroutine libGalacticusInitL

program libGalacticusInit
end program libGalacticusInit
CODE
push(
    @{$code->{'main'}},
    $libraryInitializer
    );

# Serialize the code.
system("mkdir -p ".$ENV{'BUILDPATH'}."/libgalacticus");
foreach my $className ( sort(keys(%{$code})) ) {
    my $outputFileName = $ENV{'BUILDPATH'}."/".($className eq "main" ? "libgalacticus.Inc" : "libgalacticus/".$className.".F90");
    open(my $output,">",$outputFileName);
    print $output join("\n",@{$code->{$className}})."\n";
    close($output);
}

# Generate code to initialize the Python interface,
$python->{'units'}->{'init'}->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'ext');
from ctypes import *
# Load the shared library into ctypes.
import os
# Load the shared library into ctypes.
cwd = os.getcwd()
libname = os.path.join(cwd,"galacticus/lib/libgalacticus.so")
c_lib = CDLL(libname)
c_lib.libGalacticusInitL()
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
open(my $pythonOutput,">","./galacticus.py");
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
close($pythonOutput);

exit;

sub interfacesPointerGet {
    # Build functions which return pointers to the specific type.
    my $code                   = shift();
    $ext::functionClass        = shift();
    @ext::functionClassSymbols = ( $ext::functionClass->{'name'}."Class" );
    push(@ext::functionClassSymbols,map {$_->{'name'}} @{$ext::functionClass->{'implementations'}});
    my $function = fill_in_string(<<'CODE', PACKAGE => 'ext');
function {$functionClass->{'name'}}GetPtr({$functionClass->{'name'}}_,classID)
  use, intrinsic :: ISO_C_Binding               , only : c_ptr                             , c_int, c_f_pointer
  use            :: Error                       , only : Error_Report
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
     call Error_Report('unknown classID'//\{introspection:location\})
  end select
  return
end function {$functionClass->{'name'}}GetPtr
CODE
    push(
	@{$code->{$ext::functionClass->{'name'}}},
	$function
	);
}

sub interfacesConstructors {
    # Build interfaces to constructors.
    my $code                      = shift();
    my $python                    = shift();
    $ext::functionClass           = shift();
    my $libraryFunctionClasses    = shift();
    my $extensions                = shift();
    my $moduleUsesImplementations = shift();
    foreach $ext::implementation ( @{$ext::functionClass->{'implementations'}} ) {
	# Extract the list of arguments and process to determine how interfaces should be built.
	my @argumentList = @{$ext::implementation->{'arguments'}};
	@argumentList    = &assignCTypes             (\@argumentList                                                                                );
	@argumentList    = &assignCAttributes        (\@argumentList                                                                                );
	@argumentList    = &buildPythonReassignments (\@argumentList                                                                                );
	@argumentList    = &buildFortranReassignments(\@argumentList,$ext::functionClass,$ext::implementation,$extensions,$moduleUsesImplementations);
	# Construct pre- and post-arguments content for the call from Fortran to Galacticus.
	my $preArguments  .= fill_in_string(<<'CODE', PACKAGE => 'ext');
  !![
  <referenceConstruct object="self">
   <constructor>
    {$implementation->{'name'}}( &amp;
CODE
	my $postArguments .= fill_in_string(<<'CODE', PACKAGE => 'ext');
     &amp;                     )
   </constructor>
  </referenceConstruct>
  !!]
CODE
	# Finally, build the interface constructor function.
	$ext::isoCBindingImport = &isoCBindingImport   (\@argumentList,'c_ptr'      ,'c_loc'               );
	@ext::functionArguments = &fortranArgList      (\@argumentList                                     );
	$ext::declarations      = &fortranDeclarations (\@argumentList                                     );
	$ext::reassignments     = &fortranReassignments(\@argumentList                                     );	
	$ext::moduleUseCode     = &fortranModuleUses   (\@argumentList                                     );
	$ext::callCode          = &fortranCallCode     (\@argumentList,$preArguments,$postArguments,"&amp;");
	my $constructor .= fill_in_string(<<'CODE', PACKAGE => 'ext');
function {$implementation->{'name'}}L({join(",",@functionArguments)}) bind(c,name='{$implementation->{'name'}}L')
  use :: {$functionClass->{'module'}}, only : {$implementation->{'name'}}
{$isoCBindingImport}
{$moduleUseCode}
  implicit none
  type(c_ptr                      )          :: {$implementation->{'name'}}L
  type({$implementation->{'name'}}), pointer :: self
{$declarations}

{$reassignments}
  allocate(self)
{$callCode}
  {$implementation->{'name'}}L=c_loc(self)
  return
end function {$implementation->{'name'}}L
CODE
	push(
	    @{$code->{$ext::functionClass->{'name'}}},
	    $constructor
	);
	# Add library interface descriptors.
	my @clibArgTypes = &ctypesArgTypes(\@argumentList);
	push(
	    @{$python->{'c_lib'}},
	    {
		name     => $ext::implementation->{'name'}."L",
		restype  => "c_void_p",
		argtypes => \@clibArgTypes
	    }
	    );
	# Add a constructor to the Python class.
	@ext::pythonConstructorArguments = &pythonArgList      (\@argumentList                                                           );
	$ext::pythonReassignments        = &pythonReassignments(\@argumentList                                                           );
	$ext::pythonCallCode             = &pythonCallCode     (\@argumentList,"self._glcObj = c_lib.".$ext::implementation->{'name'}."L");
	my $pythonConstructor = fill_in_string(<<'CODE', PACKAGE => 'ext');
# Constructor
def __init__({join(",",@pythonConstructorArguments)}):
    # Assign class ID so relevant pointers can be constructed on the Fortran side.
    self._classID = {$implementation->{'classID'}}
{$pythonReassignments}
{$pythonCallCode}
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
subroutine {$functionClass->{'name'}}DestructorL(self,classID) bind(c,name='{$functionClass->{'name'}}DestructorL')
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
end subroutine {$functionClass->{'name'}}DestructorL
CODE
    push(
	@{$code->{$ext::functionClass->{'name'}}},
	$destructor
	);
    # Add c_lib interface.
    push(
	@{$python->{'c_lib'}},
	{
	    name     => $ext::functionClass->{'name'}."DestructorL",
	    restype  => undef(),
	    argtypes => [ "c_void_p", "c_int" ]
	}
	);
    # Add a destructor to the Python class.
    my $pythonDestructor = fill_in_string(<<'CODE', PACKAGE => 'ext');
# Destructor
def __del__(self):
    c_lib.{$functionClass->{'name'}}DestructorL(self._glcObj,self._classID)
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
    my $extensions                = shift();
    my $moduleUsesImplementations = shift();
    foreach $ext::method ( &List::ExtraUtils::hashList($ext::functionClass->{'methods'},keyAs => 'name') ) {
	# Construct function declaration.
	my %isoCBindingSymbols;
	$ext::procedure             = "function";
	$ext::resultConversionOpen  = "";
	$ext::resultConversionClose = "";
	my $functionCType;
	my $clibResType;
	if ( $ext::method->{'type'} eq "double precision" ) {
	    $functionCType                       = "real(c_double)";
	    $clibResType                         = "c_double";
	    $isoCBindingSymbols{'c_double'} = 1;
	} elsif ( $ext::method->{'type'} eq "logical" ) {
	    $functionCType                       = "logical(c_bool)";
	    $clibResType                         = "c_bool";
	    $ext::resultConversionOpen           = "logical(";
	    $ext::resultConversionClose          = ",kind=c_bool)";
	    $isoCBindingSymbols{'c_bool'} = 1;
	} elsif ( $ext::method->{'type'} eq "void" ) {
	    	$ext::procedure = "subroutine";
	} else {
	    print "unsupported type '".$ext::method->{'type'}."'\n";
	    delete($ext::functionClass->{'methods'}->{$ext::method});
	    next;
	}
	$ext::functionDeclaration  = $ext::method->{'type'} eq "void" ? "" : $functionCType." :: ".$ext::functionClass->{'name'}.ucfirst($ext::method->{'name'})."L\n";
	# Create a list of arguments.
	my @argumentList;
	# Add the self argument.
	push(
	    @argumentList,
	    {
		intrinsic     => "class",
		type          => $ext::functionClass->{'name'}."Class",
		attributes    => [ "intent(inout)" ],
		name          => "self"
	    }
	    );
	# Add all other arguments.
	if ( exists($ext::method->{'argument'}) ) {
	    foreach my $argument ( &List::ExtraUtils::as_array($ext::method->{'argument'}) ) {
		my $declaration = &Galacticus::Build::SourceTree::Parse::Declarations::parseDeclaration($argument);
		foreach my $name ( @{$declaration->{'variableNames'}} ) {
		    my $declarationSingle =
		    {
			intrinsic     => $declaration->{'intrinsic' },
			type          => $declaration->{'type'      },
			attributes    => $declaration->{'attributes'},
			name          => $name
		    };
		    push(
			@argumentList     ,
			$declarationSingle
			);
		}
	    }
	}
	# Process argument list.
	@argumentList    = &assignCTypes             (\@argumentList                                                                   );
	@argumentList    = &assignCAttributes        (\@argumentList                                                                   );
	@argumentList    = &buildPythonReassignments (\@argumentList                                                                   );
	@argumentList    = &buildFortranReassignments(\@argumentList,$ext::functionClass,undef(),$extensions,$moduleUsesImplementations);
	# Generate the pre- and post-arguments content of the function call.
	my $preArguments  = ($ext::method->{'type'} eq "void" ? "call " : $ext::functionClass->{'name'}.ucfirst($ext::method->{'name'})."L=").$ext::resultConversionOpen."self_%".$ext::method->{'name'}."( &\n";
	my $postArguments = "& )".$ext::resultConversionClose."\n";
	# Generate the function.
	$ext::isoCBindingImport = &isoCBindingImport   (\@argumentList,keys(%isoCBindingSymbols)                   );
	@ext::functionArguments = &fortranArgList      (\@argumentList                                             );
	$ext::declarations      = &fortranDeclarations (\@argumentList                                             );
	$ext::reassignments     = &fortranReassignments(\@argumentList                                             );	
	$ext::moduleUseCode     = &fortranModuleUses   (\@argumentList                                             );
	$ext::callCode          = &fortranCallCode     (\@argumentList,$preArguments            ,$postArguments,"&");
	my $function = fill_in_string(<<'CODE', PACKAGE => 'ext');
{$procedure} {$functionClass->{'name'}}{ucfirst($method->{'name'})}L({join(",",@functionArguments)}) bind(c,name='{$functionClass->{'name'}}{ucfirst($method->{'name'})}L')
  use :: {$functionClass->{'module'}}, only : {$functionClass->{'name'}}Class
{$isoCBindingImport}
{$moduleUseCode}
  implicit none
{$functionDeclaration}
{$declarations}

{$reassignments}
{$callCode}
  return
end {$procedure} {$functionClass->{'name'}}{ucfirst($method->{'name'})}L
CODE
	push(
	    @{$code->{$ext::functionClass->{'name'}}},
	    $function
	    );
	# Add library interface descriptors.
	my @clibArgTypes = &ctypesArgTypes(\@argumentList);
	push(
	    @{$python->{'c_lib'}},
	    {
		name     => $ext::functionClass->{'name'}.ucfirst($ext::method->{'name'})."L",
		restype  => $clibResType,
		argtypes => \@clibArgTypes
	    }
	    );
	# Add a method to the Python class.
	@ext::pythonMethodArguments = &pythonArgList (\@argumentList                                                                                  );
	$ext::pythonCallCode        = &pythonCallCode(\@argumentList,"return c_lib.".$ext::functionClass->{'name'}.ucfirst($ext::method->{'name'})."L");
	my $pythonMethod = fill_in_string(<<'CODE', PACKAGE => 'ext');
def {$method->{'name'}}({join(",",@pythonMethodArguments)}):
{$pythonCallCode}
CODE
	push(
	    @{$python->{'units'}->{$ext::functionClass->{'name'}}->{'subUnits'}},
	    {
		content => $pythonMethod
	    }
	    );	
    }
}

sub interfacesPythonClasses {
    # Build Python parent class. This corresponds to a Galacticus "functionClass" class.
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
    
    # Add child classes. These correspond to each implementation of the Galacticus "functionClass".
    foreach $ext::implementation ( @{$ext::functionClass->{'implementations'}} ) {
	my $class = fill_in_string(<<'CODE', PACKAGE => 'ext');
class {$implementation->{'name'}}({$functionClass->{'name'}}):
CODE
	$python->{'units'}->{$ext::implementation->{'name'}}->{'content'     } = $class;
	$python->{'units'}->{$ext::implementation->{'name'}}->{'indent'      } = 0;
	$python->{'units'}->{$ext::implementation->{'name'}}->{'dependencies'} = [ $ext::functionClass->{'name'} ];	
    }
}

sub assignCTypes {
    # Assign appropriate C-types for each argument in the list.
    my @argumentList = @{shift()};
    my @argumentListNew;
    while ( @argumentList ) {
	my $argument = pop(@argumentList);
	# Detect if this argument is optional.
	$argument->{'isOptional'} = grep {$_ eq "optional"} @{$argument->{'attributes'}};
	# Mark the argument as present in all argument lists.
	$argument->{'fortran'   }->{'isPresent'} = 1;
	$argument->{'python'    }->{'isPresent'} = 1;
	$argument->{'galacticus'}->{'isPresent'} = 1;
	# Mark the argument as not a functionClass object by default.
	$argument->{'isFunctionClass'} = 0;
	# Determine the type of the argument.
	if ( $argument->{'intrinsic'} eq "double precision" ) {
	    # Double precision type is represented as a C double.
	    $argument->{'ctypes' }->{'type'} = "c_double";
	    $argument->{'fortran'}->{'type'} = "real(c_double)";
	} elsif ( $argument->{'intrinsic'} eq "integer"     ) {
	    # Integer type - check the kind.
	    if ( defined($argument->{'type'}) ) {
		die("Unknown integer kind '".$argument->{'type'}."' for argument '".$argument->{'name'}."'");
	    } else {
		# No kind was specified - use a C int.
		$argument->{'ctypes' }->{'type'} = "c_int";
		$argument->{'fortran'}->{'type'} = "integer(c_int)";
	    }
	} elsif ( $argument->{'intrinsic'} eq "logical" ) {
	    # Logical type is represented as a C bool.
	    $argument->{'ctypes' }->{'type'} = "c_bool";
	    $argument->{'fortran'}->{'type'} = "logical(c_bool)";
	} elsif ( $argument->{'intrinsic'} eq "character" ) {
	    # Character type is represented as a C char pointer.
	    $argument->{'ctypes' }->{'type'} = "c_char_p";
	    $argument->{'fortran'}->{'type'} = "character(c_char)";
	} elsif ( $argument->{'intrinsic'} eq "type" ) {
	    # Derived-type argument. Check the type.
	    if ( $argument->{'type'} eq "varying_string" ) {
		# Varying strings are represented as C char pointers.
		$argument->{'ctypes' }->{'type'} = "c_char_p";
		$argument->{'fortran'}->{'type'} = "character(c_char)";
	    } elsif ( $argument->{'type'} =~ m/^enumeration[a-z0-9_]+type/i ) {
		# Enumerations are represented as C ints.
		$argument->{'ctypes' }->{'type'} = "c_int";
		$argument->{'fortran'}->{'type'} = "integer(c_int)";
	    } else {
		# Other types are represented as opaque (void) pointers.
		$argument->{'ctypes' }->{'type'} = "c_void_p";
		$argument->{'fortran'}->{'type'} = "type(c_ptr)";
	    }
	} elsif ( $argument->{'intrinsic'} eq "class" ) {
	    # Derived-type class arguments are represented as opaque (void) pointers.
	    $argument->{'ctypes' }->{'type'} = "c_void_p";
	    $argument->{'fortran'}->{'type'} = "type(c_ptr)";
	    # Check for a functionClass argument.
	    if ( my @matchedClass = grep {$argument->{'type'} eq $_."Class"} keys(%{$libraryFunctionClasses->{'classes'}}) ) {
		# Mark as a functionClass argument.
		$argument->{'isFunctionClass'} = 1;
		# If this is the "self" argument, it should not be passed to Galacticus (as we use the method bound to "self" instead).
		$argument->{'galacticus'}->{'isPresent'} = 0
		    if ( $argument->{'name'} eq "self" );
		# For functionClass types we must add an extra argument which passes an integer specifying the concrete type of
		# the class. This argument will not be present in the Python function arguments.
		my $argumentNew =
		         {
			     name       => $argument->{'name'}."_ID",
			     intrinsic  => "integer"                ,
			     type       => undef()                  ,
			     attributes => [ "intent(in   )" ],
			     ctypes     =>
			                   {
					       type      => "c_int"
					   },
			     fortran    =>
			                   {
					       type      => "integer(c_int)",
					       isPresent => 1
					   },
			     python     =>
			                   {
					       isPresent => 0
					   },
			     galacticus =>
			                   {
					       isPresent => 0
					   }
			 };
		if ( $argument->{'isOptional'} ) {
		    push(@{$argumentNew->{'attributes'}},"optional");
		    $argumentNew->{'isOptional'}              = 1;
		    $argumentNew->{'python'    }->{'present'} = $argument->{'name'};
		}
		unshift(@argumentListNew,$argumentNew);
	    }
	} else {
	    die("Unsupported argument type '".$argument->{'type'}."'");
	}
	# Put this argument onto the new list.
	unshift(@argumentListNew,$argument);
    }
    return @argumentListNew;
}

sub assignCAttributes {
    # Assign appropriate attributes for each argument in the list.
    my @argumentList = @{shift()};
    my @argumentListNew;
    while ( @argumentList ) {
	# Get the next argument from the list.
	my $argument = pop(@argumentList);
	# Build attributes.
	## Add "optional" attribute for any optional argument.
 	push(@{$argument->{'fortran'}->{'attributes'}},"optional")
	    if ( $argument->{'isOptional'}            );
	## Transfer any dimension attributes.
	push(@{$argument->{'fortran'}->{'attributes'}},grep {$_ =~ m/^dimension/ || $_ =~ m/^allocatable/} @{$argument->{'attributes'}});
	## String argument passed as type "c_char_p" must be "dimension(*)" in Fortran.
 	push(@{$argument->{'fortran'}->{'attributes'}},"dimension(*)")
	    if ( $argument->{'ctypes'}->{'type'} eq "c_char_p" );
	## Determine whether to pass by value or by reference.
	$argument->{'passBy'    } =
	    (
	     # Can pass by value if the C-type is a pointer type, or if intent is "in".
	     (
	      $argument->{'ctypes'}->{'type'} =~ m/_p$/
	      ||
	      ( grep {$_ =~ m/intent\s*\(\s*in\s*\)/} @{$argument             ->{'attributes'}})
	     )
	     # Optional arguments can not be passed by value.
	     &&
	      !                                         $argument             ->{'isOptional'}
	     # Non-scalar arguments can not be passed by value.
	     &&
	     (! grep {$_ =~ m/^dimension/           } @{$argument->{'fortran'}->{'attributes'}})
	    )
	    ?
	    "value"
	    :
	    "reference";
	## Add "value" attribute for any argument passed by value.
 	push(@{$argument->{'fortran'}->{'attributes'}},"value")
	    if ( $argument->{'passBy'} eq "value" );
	# Record whether ctypes should pass a pointer. This is the case if we pass by reference, except for c_char_p types which
	# are passed as a pointer value from Python, but received as a pointer reference by Fortran.
	$argument->{'ctypes'}->{'pointer'} =
	    $argument->{'passBy'}           eq "reference"
	    &&
	    $argument->{'ctypes'}->{'type'} ne "c_char_p" ;
	# Put this argument onto the new list.
	unshift(@argumentListNew,$argument);	
    }
    return @argumentListNew;
}

sub buildPythonReassignments {
    # Geneate reassignments of Python arguments to allow passing between languages.
    my @argumentList = @{shift()};
    my @argumentListNew;
    while ( @argumentList ) {
	# Get the next argument from the list.
	my $argument = pop(@argumentList);
	# For optional functionClass arguments we must use temporary variables within the Python function.
	if ( $argument->{'isFunctionClass'} ) {
	    # Get the class ID argument.
	    my $argumentID = shift(@argumentListNew);
	    if ( $argument->{'isOptional'} ) {
		# Set the name to use when passing this argument, and for the following argument (which will be the class ID number).
		$argument  ->{'python'}->{'passAs'} = $argument->{'name'}."_glcObj" ;
		$argumentID->{'python'}->{'passAs'} = $argument->{'name'}."_classID";
		# Create the reassignment code. If the argument is present, simply extract the opaque pointer and class ID into the
		# variable that will be passed. Otherwise, set those to "None".
		$ext::argumentName = $argument->{'name'};
		$argument->{'python'}->{'reassignment'} = fill_in_string(<<'CODE', PACKAGE => 'ext');
    if {$argumentName}:
        {$argumentName}_glcObj ={$argumentName}._glcObj
        {$argumentName}_classID={$argumentName}._classID
    else:
        {$argumentName}_glcObj =None
        {$argumentName}_classID=None
CODE
	    } else {
		# Set the name to use when passing this argument, and for the following argument (which will be the class ID
		# number). Argument is non-optional, so we can access the relevant properties directly.
		$argument  ->{'python'}->{'passAs'} = $argument->{'name'}."._glcObj" ;
		$argumentID->{'python'}->{'passAs'} = $argument->{'name'}."._classID";
	    }
	    # Push the class ID argument back onto the list.
	    unshift(@argumentListNew,$argumentID);
	}
	# Put this argument onto the new list.
	unshift(@argumentListNew,$argument);
    }
    return @argumentListNew;
}

sub buildFortranReassignments {
    # Geneate reassignments of Fortran arguments to allow passing between languages.
    my @argumentList              = @{shift()};
    my $functionClass             =   shift() ;
    my $implementation            =   shift() ;
    my $extensions                =   shift() ;
    my $moduleUsesImplementations =   shift() ;
    my @argumentListNew;
    while ( @argumentList ) {
	# Get the next argument from the list.
	my $argument = pop(@argumentList);
	# Determine any reassignments.
	if ( $argument->{'intrinsic'} eq "logical" ) {
	    # Logical arguments are passed in as type "c_bool". They must be reassigned to a regular "logical".
	    $argument->{'fortran'}->{'reassignment'} = ($argument->{'isOptional'} ? "if (present(".$argument->{'name'}."))" : "").$argument->{'name'}."_=logical(".$argument->{'name'}.")\n";
	    $argument->{'fortran'}->{'declarations'} = "logical :: ".$argument->{'name'}."_\n";
	    $argument->{'fortran'}->{'passAs'      } = $argument->{'name'}."_";
	} elsif ( $argument->{'intrinsic'} eq "character" ) {
	    # Character arguments are passed as type c_char_p. We can convert it to varying_string via the String_C_to_Fortran() function.
	    $argument->{'fortran'}->{'modules'}->{'String_Handling'   }->{'String_C_to_Fortran'} = 1;
	    $argument->{'fortran'}->{'modules'}->{'ISO_Varying_String'}->{'char'               } = 1;
	    $argument->{'fortran'}->{'passAs' }                                                  = "char(String_C_to_Fortran(".$argument->{'name'}."))"
	} elsif ( $argument->{'intrinsic'} eq "type" ) {
	    if ( $argument->{'type'} eq "varying_string" ) {
		# varying_string arguments are passed as type c_char_p. We can convert it to varying_string via the String_C_to_Fortran() function.
		$argument->{'fortran'}->{'modules'}->{'String_Handling'}->{'String_C_to_Fortran'} = 1;
		$argument->{'fortran'}->{'passAs' }                                               = "String_C_to_Fortran(".$argument->{'name'}.")"
	    } elsif ( $argument->{'type'} =~ m/^enumeration[a-z0-9_]+type$/i ) {
		# Enumerations are passed as type c_int and must be set into the appropriate enumeration type.
		$argument  ->{'fortran'}->{'declarations'} = "type(".$argument->{'type'}.") :: ".$argument->{'name'}."_\n";
		$argument  ->{'fortran'}->{'passAs'      } = $argument->{'name'}."_";
		$argument  ->{'fortran'}->{'reassignment'} = ((grep {$_ eq "optional"} @{$argument->{'attributes'}}) ? "if (present(".$argument->{'name'}.")) " : "").$argument->{'name'}."_%ID=".$argument->{'name'}."\n";
		# Search for any module import of this enumeration type.
		my $importModule;
		if ( defined($implementation) ) {

		    my $className = $implementation->{'name'};
		    while ( defined($className) && ! defined($importModule) ) {
			foreach my $useBlock ( @{$moduleUsesImplementations->{$className}} ) {
			    foreach my $module ( keys(%{$useBlock}) ) {
				if ( grep {$_ eq $argument->{'type'}} keys(%{$useBlock->{$module}->{'only'}}) ) {
				    $importModule = $module;
				}
			    }
			}
			$className = exists($extensions->{$className}) ? $extensions->{$className} : undef();
		    }
		}
		unless ( defined($importModule) ) {
		    foreach my $useBlock ( @{$functionClass->{'moduleUses'}} ) {
			foreach my $module ( keys(%{$useBlock}) ) {
			    if ( grep {$_ eq $argument->{'type'}} keys(%{$useBlock->{$module}->{'only'}}) ) {
				$importModule = $module;
			    }
			}
		    }
		}
		unless ( defined($importModule) ) {
		    $importModule = $functionClass->{'module'};
		}
		$argument  ->{'fortran'}->{'modules'}->{$importModule}->{$argument->{'type'}}  = 1
		    if ( defined($importModule) );
	    } elsif ( $argument->{'type'} eq "treeNode" ) {
		# treeNode arguments are passed as type c_ptr and must be dereferenced.
		$argument  ->{'fortran'}->{'declarations'      }                                      = "type(".$argument->{'type'}."), pointer :: ".$argument->{'name'}."_\n";
		$argument  ->{'fortran'}->{'passAs'            }                                      = $argument->{'name'}."_";
		$argument  ->{'fortran'}->{'reassignment'      }                                      = "call c_f_pointer(".$argument->{'name'}.",".$argument->{'name'}."_)\n";
		$argument  ->{'fortran'}->{'reassignment'      }                                      = "if (present(".$argument->{'name'}.")) then\n ".$argument  ->{'fortran'}->{'reassignment'}."else\n ".$argument->{'name'}."_ => null()\nend if\n"
		    if ( $argument->{'isOptional'} );
		@{$argument->{'fortran'}->{'isoCBindingSymbols'}                                    } = ( "c_f_pointer" );
		$argument  ->{'fortran'}->{'modules'           }->{'Galacticus_Nodes'}->{'treeNode'}  = 1;
	    } else {
		# Derived types (other than "varying_string") are passed as opaque pointers, and must be dereferenced.
		die("no functionClass defined for this function")
		    unless ( defined($functionClass) );
		$argument  ->{'fortran'}->{'declarations'      }                                                       = "type(".$argument->{'type'}."), pointer :: ".$argument->{'name'}."_\n";
		$argument  ->{'fortran'}->{'passAs'            }                                                       = $argument->{'name'}."_";
		$argument  ->{'fortran'}->{'reassignment'      }                                                       = ($argument->{'isOptional'} ? "if (present(".$argument->{'name'}.")) " : "")."call c_f_pointer(".$argument->{'name'}.",".$argument->{'name'}."_)\n";
		@{$argument->{'fortran'}->{'isoCBindingSymbols'}                                                     } = ( "c_f_pointer" );
		$argument  ->{'fortran'}->{'modules'           }->{$functionClass->{'module'}}->{$argument->{'type'}}  = 1;
	    }
	} elsif ( $argument->{'isFunctionClass'} ) {
	    # functionClass arguments must be dereferenced via a specialized function.
	    (my $className = $argument->{'type'}) =~ s/Class$//;
	    my @matchedClass = grep {$className eq $_} keys(%{$libraryFunctionClasses->{'classes'}});
	    die("unable to locate matching functionClass")
		unless ( scalar(@matchedClass) == 1 );
	    my $moduleName = $libraryFunctionClasses->{'classes'}->{$matchedClass[0]}->{'module'};
	    $argument->{'fortran'}->{'declarations' }                                        = "class(".$argument->{'type'}."), pointer :: ".$argument->{'name'}."_\n";
	    $argument->{'fortran'}->{'reassignment' }                                        = ($argument->{'isOptional'} ? "if (present(".$argument->{'name'}.")) " : "").$argument->{'name'}."_ => ".$className."GetPtr(".$argument->{'name'}.",".$argument->{'name'}."_ID)\n";
	    $argument->{'fortran'}->{'passAs'       }                                        = $argument->{'name'}."_";
	    $argument->{'fortran'}->{'functionClass'}                                        = $className; 
	    $argument->{'fortran'}->{'modules'      }->{$moduleName}->{$argument->{'type'}}  = 1;
	} elsif ( $argument->{'intrinsic'} eq "class" && $argument->{'type'} eq "*" ) {
	    # Unlimited polymorphic argument. We must find the concrete type to use for this argument.
	    die("no functionClass defined for this function")
		unless ( defined($functionClass ) );
	    die("no implementation defined for this function")
		unless ( defined($implementation) );
	    die("no concrete type information provided for argument '".$argument->{'name'}."' of class '".$implementation->{'name'}."'")
	     	unless ( exists($libraryFunctionClasses->{'classes'}->{$functionClass->{'name'}}->{$implementation->{'name'}}->{'constructor'}->{'argument'}) );
	    my @concreteTypes = grep {$_->{'name'} eq $argument->{'name'}} &List::ExtraUtils::as_array($libraryFunctionClasses->{'classes'}->{$functionClass->{'name'}}->{$implementation->{'name'}}->{'constructor'}->{'argument'});
	    $argument  ->{'fortran'}->{'modules'           }->{$concreteTypes[0]->{'module'}}->{$concreteTypes[0]->{'type'}}  = 1;
	    $argument  ->{'fortran'}->{'declarations'      }                                                                  = "type(".$concreteTypes[0]->{'type'}."), pointer :: ".$argument->{'name'}."_\n";
	    $argument  ->{'fortran'}->{'reassignment'      }                                                                  = ($argument->{'isOptional'} ? "if (present(".$argument->{'name'}."))" : "")."call c_f_pointer(".$argument->{'name'}.",".$argument->{'name'}."_)\n";
	    $argument  ->{'fortran'}->{'passAs'            }                                                                  = $argument->{'name'}."_";
	    @{$argument->{'fortran'}->{'isoCBindingSymbols'}                                                                } = ( "c_f_pointer" );
	}
	# Put this argument onto the new list.
	unshift(@argumentListNew,$argument);
    }
    return @argumentListNew;
}

sub ctypesArgTypes {
    # Geneate an ".argtypes" definition for a function.
    my @argumentList = @{shift()};
    my @argTypes;
    foreach my $argument ( @argumentList ) {
	# Wrap the type inside a "POINTER" function if it should be passed by reference.
	my $ctype =
	               $argument->{'ctypes'}->{'pointer'}
	    ?
	    "POINTER(".$argument->{'ctypes'}->{'type'   }.")"
	    :
	               $argument->{'ctypes'}->{'type'   }    ;
	# Add the type.
	push(
	    @argTypes,
	    $ctype
	    );
    }
    return @argTypes;
}

sub pythonArgList {
    # Generate an argument list for the Python function.
    my @argumentList = @{shift()};
    # Initialize the argument list with "self".
    my @argList = (scalar(@argumentList) > 0 && $argumentList[0]->{'name'} eq "self") ? () : ( "self" );
    # Add all other arguments.
    my $firstOptionalFound = 0;
    foreach my $argument ( @argumentList ) {
	# Skip arguments not present in the Python function.
	next
	    unless ( $argument->{'python'}->{'isPresent'} );
	# Record when we see the first optional argument.
	$firstOptionalFound = 1
	    if ( $argument->{'isOptional'} );
	# Assign the argument name.
	my $pythonArg = $argument->{'name'};
	# For all arguments including and after the first optional argument we must provide a default (of "None").
	$pythonArg .= "=None"
	    if ( $firstOptionalFound );
	# Push this argument onto our list.
	push(
	    @argList,
	    $pythonArg
	    );
    }
    return @argList;
}

sub pythonReassignments {
    # Generate reassignments of arguments in the Python function.
    my @argumentList = @{shift()};
    # Build reassignments.
    return join("\n",map {exists($_->{'python'}->{'reassignment'}) ? $_->{'python'}->{'reassignment'} : ()} @argumentList );
}

sub pythonCallCode {
    # Generate an argument list for the Python call.
    my @argumentList = @{shift()};
    my $call         =   shift() ;
    # Initialize the code;
    my $code;
    # Determine if any optional arguments are present.
    my $hasOptionals = grep {$_->{'isOptional'}} @argumentList;
    # If there are no optional arguments simply pass all arguments by name directly.
    if ( ! $hasOptionals ) {
	# No optional arguments - pass all arguments directly by name. No need for any type conversion or explicit "byref()" as
	# the ctypes.arglist specifications will handle this automatically for us.
	$code =
	    "    ".
	    $call.
	    "(".
	    join
	    (
	     ",",
	     map
	     {
		 $_->{'fortran'}->{'isPresent'}                                                 # Skip arguments not present in the Fortran function.
		 ?
	         exists($_->{'python'}->{'passAs'}) ? $_->{'python'}->{'passAs'} : $_->{'name'} # Use the "passAs" name if defined, otherwise the regular argument name.
		 :
	         ()
	     }
	     @argumentList
	    ).
	    ")";
    } else {
	# Some optional arguments are present.
	my %presentNames;
	foreach my $argument ( @argumentList ) {
	    next
		unless ( $argument->{'isOptional'} );
	    # Assign default argument names to check for argument presence.
	    $argument->{'python'}->{'present'} = $argument->{'name'}
	        unless ( exists($argument->{'python'}->{'present'}) );
	    # Accumulate a dictionary of argument names for presence checks.
	    ++$presentNames{$argument->{'python'}->{'present'}};
	}
	# Determine names of optional arguments.
	my @optionalNames  = sort(keys(%presentNames));
	my $countOptionals = scalar(@optionalNames);
	# Iterate over all possible combinations of present parameters.
	my $formatBinary = "%.".$countOptionals."b";
	for(my $i=0;$i<2**$countOptionals;++$i) {
	    my @state = split(//,sprintf($formatBinary,$i));
	    my %present;
	    for(my $i=0;$i<scalar(@optionalNames);++$i) {
		$present{$optionalNames[$i]} = $state[$i];
	    }
	    # Build the condition for this state.
	    $code .= "    ".($i > 0 ? "el" : "")."if ".join(" and ",map {($present{$_} ? "" : "not ").$_} @optionalNames).":\n";
	    # Add all arguments.
	    my @argList;
	    my $firstOptionalFound = 0;
	    foreach my $argument ( @argumentList ) {
	    	# Skip arguments not present in the Fortran function.
	    	next
	    	    unless ( $argument->{'fortran'}->{'isPresent'} );
	    	# Record when we see the first optional argument.
	    	$firstOptionalFound = 1
	    	    if ( $argument->{'isOptional'} );
	    	# Assign the argument name. This is just the actual argument name, unless a specific "pass as" name has been specified for
	    	# the Python function.
	    	my $callArg = exists($argument->{'python'}->{'passAs'}) ? $argument->{'python'}->{'passAs'} : $argument->{'name'};
		# If we are at or past the first optional argument then we must perform explicit conversion of the argument as
		# ctypes does not know how to do it for us.
		if ( $firstOptionalFound ) {
		    # Add explicit conversion to the relevant ctype.
		    $callArg = $argument->{'ctypes' }->{'type'}."(".$callArg.")";
		    # Add explicit passing by reference is needed.
		    $callArg = "byref(".$callArg.")"
			if ( $argument->{'passBy'} eq "reference" );
		}
		# Check the state for this argument, and replace it with "None" if not present.
		$callArg = "None"
		    unless ( ! $argument->{'isOptional'} || $present{$argument->{'python'}->{'present'}} );
	    	# Push this argument onto our list.
	    	push(
	    	    @argList,
	    	    $callArg
	    	    );
	    }
	    $code .= "        ".$call."(".join(",",@argList).")\n";
	}
    }
    return $code;
}

sub fortranArgList {
    # Generate an argument list for the Fortran function.
    my @argumentList = @{shift()};
    # Create the argument list.
    my @argList = map {$_->{'fortran'}->{'isPresent'} ? $_->{'name'} : ()} @argumentList;
    return @argList;
}

sub fortranDeclarations {
    # Generate variable declarations list for the Fortran function.
    my @argumentList = @{shift()};
    # Initialize the declaration code.
    my $code = "";
    # Initialize list of functionClasses needed.
    my %functionClasses;
    # Iterate over all arguments.
    foreach my $argument ( @argumentList ) {
	# Add a declaration for this argument.
	my @attributes = @{$argument->{'fortran'}->{'attributes'}};
	$code .= "  ".$argument->{'fortran'}->{'type'}.(@attributes ? ", " : "").join(", ",@attributes)." :: ".$argument->{'name'}."\n";
	# Add any extra declarations needed for this argument.
	$code .= $argument->{'fortran'}->{'declarations'}
	    if ( exists($argument->{'fortran'}->{'declarations'}) );
	# Accumulate any functionClass needed.
	$functionClasses{$argument->{'fortran'}->{'functionClass'}} = 1
	    if ( exists($argument->{'fortran'}->{'functionClass'}) );
    }
    # Add interfaces for any functionClass pointer dereferencing functions.
    foreach my $functionClass ( sort(keys(%functionClasses)) ) {
	$code .="interface\n";
	$code .=" function ".$functionClass."GetPtr(ptr_,classID)\n";
	$code .="  import c_int, c_ptr, ".$functionClass."Class\n";
	$code .="  class(".$functionClass."Class), pointer :: ".$functionClass."GetPtr\n";
	$code .="  type   (c_ptr), intent(in   ) :: ptr_\n";
	$code .="  integer(c_int), intent(in   ) :: classID\n";
	$code .=" end function ".$functionClass."GetPtr\n";
	$code .="end interface\n";
    }
    # Return the declarations.
    return $code;
}

sub fortranReassignments {
    # Generate variable reassignments for the Fortran function.
    my @argumentList = @{shift()};
    # Generate the code.
    my $code = join("",map {exists($_->{'fortran'}->{'reassignment'}) ? $_->{'fortran'}->{'reassignment'} : ()} @argumentList);
    return $code;
}

sub fortranModuleUses {
    # Generate module use statements for the Fortran function.
    my @argumentList = @{shift()};

    # Initialize module data structure.
    my $modules;
    # Iterate over all arguments accumulating a dictionary of required symbols.
    foreach my $argument ( @argumentList ) {
	next
	    unless ( exists($argument->{'fortran'}->{'modules'}) );
	foreach my $moduleName ( keys(%{$argument->{'fortran'}->{'modules'}}) ) {
	    foreach my $symbolName ( keys(%{$argument->{'fortran'}->{'modules'}->{$moduleName}}) ) {
		++$modules->{$moduleName}->{$symbolName};
	    }
	}
    }
    # Generate the code.
    my $code = "";
    foreach my $moduleName ( keys(%{$modules}) ) {
	$code .= "use :: ".$moduleName.", only : ".join(", ",sort(keys(%{$modules->{$moduleName}})))."\n";
    }
    return $code;
}

sub fortranCallCode {
    # Generate an argument list for the Fortran call.
    my @argumentList  = @{shift()};
    my $preArguments  =   shift() ;
    my $postArguments =   shift() ;
    my $continuation  =   shift() ;
    # Initialize the code;
    my $code;
    # Determine if any optional arguments are present.
    my $hasOptionals = grep {$_->{'isOptional'}} @argumentList;
    # If there are no optional arguments simply pass all arguments by name directly.
    my $codeUnconditional =
	$preArguments.
	$continuation." ".
	join
	(
	 ",",
	 map
	 {
	     $_->{'galacticus'}->{'isPresent'}                                                # Skip arguments not present in the Galacticus function.
	     ?
	         exists($_->{'fortran'}->{'passAs'}) ? $_->{'fortran'}->{'passAs'} : $_->{'name'} # Use the "passAs" name if defined, otherwise the regular argument name.
	     :
	         ()
	 }
	 @argumentList
	).
	$continuation."\n".
	$postArguments;
    if ( ! $hasOptionals ) {
	# No optional arguments - pass all arguments directly.
	$code = $codeUnconditional;
    } else {
	# Some optional arguments are present.
	my %presentNames;
	foreach my $argument ( @argumentList ) {
	    next
		unless ( $argument->{'isOptional'} && $argument->{'galacticus'}->{'isPresent'} );
	    # Accumulate a dictionary of argument names for presence checks.
	    ++$presentNames{$argument->{'name'}};
	}
	# Determine names of optional arguments.
	my @optionalNames  = sort(keys(%presentNames));
	my $countOptionals = scalar(@optionalNames);
	# Iterate over all possible combinations of present parameters.
	my $formatBinary = "%.".$countOptionals."b";
	for(my $i=0;$i<2**$countOptionals;++$i) {
	    my @state = split(//,sprintf($formatBinary,$i));
	    my %present;
	    for(my $i=0;$i<scalar(@optionalNames);++$i) {
		$present{$optionalNames[$i]} = $state[$i];
	    }
	    # Build the condition for this state.
	    $code .= ($i > 0 ? "else " : "")."if (".join(" .and. ",map {($present{$_} ? "" : ".not.")."present(".$_.")"} @optionalNames).") then\n";
	    # Add all arguments.
	    my @argList;
	    foreach my $argument ( @argumentList ) {
	    	# Skip arguments not present in the Galacticus function, or which do not have their state set.
	    	next
	    	    unless ( $argument->{'galacticus'}->{'isPresent'} && ( ! $argument->{'isOptional'} || $present{$argument->{'name'}} ) );
	    	# Assign the argument name. This is just the actual argument name, unless a specific "pass as" name has been specified for
	    	# the Fortran function.
	    	my $callArg = $argument->{'name'}."=".(exists($argument->{'fortran'}->{'passAs'}) ? $argument->{'fortran'}->{'passAs'} : $argument->{'name'});
	    	# Push this argument onto our list.
	    	push(
	    	    @argList,
	    	    $callArg
	    	    );
	    }
	    $code .= $preArguments.$continuation." ".join(",",@argList)." ".$continuation."\n".$postArguments;
	}
	# Include the unconditional code in an "else" block. This will never be called but avoids compiler warnings about the function result possibly being unset.
	$code .= "else\n";
	$code .= $codeUnconditional;
	$code .= "end if\n";
    }
    return $code;
}

sub isoCBindingImport {
    # Generate code to import symbols from ISO_C_Binding.
    my @argumentList = @{shift()};
    my @extraSymbols = @_;
    # Initialize the list of required symbols with any passed directly to this function.
    my %symbols;
    $symbols{$_} = 1
	foreach ( @extraSymbols );
    # Add symbols from all arguments.
    foreach my $argument ( @argumentList ) {
	# Extract the ISO_C_Binding symbol from the Fortran type.
	if ( $argument->{'fortran'}->{'type'} =~ m/\(([a-z_]+)\)/ ) {
	    $symbols{$1} = 1;
	} else {
	    die("can not determine ISO_C_Binding symbol");
	}
	# Add any additional symbols needed for this argument.
	if ( exists($argument->{'fortran'}->{'isoCBindingSymbols'}) ) {
	    $symbols{$_} = 1
		foreach ( @{$argument->{'fortran'}->{'isoCBindingSymbols'}} );
	}
    }
    return "   use, intrinsic :: ISO_C_Binding, only : ".join(", ",sort(keys(%symbols)))."\n";
}
