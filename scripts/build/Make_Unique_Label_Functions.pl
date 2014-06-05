#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;
use File::Slurp;
use Digest::MD5 qw(md5_base64);
require Fortran::Utils;
require Galacticus::Build::Directives;
require Galacticus::Build::Make;
require List::ExtraUtils;
require File::Match;

# Construct a unique string that defines the operation of a given module. A function is constructed for every instance of a
# "uniqueLabel" directive embedded in the source code. This function is constructed as follows:
#
# 1) Find all files on which the file with the embedded "uniqueLabel" depends. (These are the "dependent files".)
#
# 2) Scan those files for instances of input parameters being read in.
#
# 3) Append the name of that parameter and its value (or default value) to the return value of the function on the following conditions:
#
#  a) For parameters not part of a Galacticus method, always append the parameter;
#
#  b) For parameters part of an old-style Galacticus method, append the parameter only if the specific implementation of the
#    method to which it belongs is selected at run time;
#
#  c) For parameters part of a new-style Galacticus method, append the parameter only if the specific implementation:
#
#     i) is directly used by any  dependent file;
#
#    ii) is used indirectly by virtue of being extended by some other implementation that is directly or indirectly used by any
#        dependent file;
#
#   iii) if any dependent file uses the default implementation, and the appropriate implementation to which the parameter belongs
#        is selected at run time.
#
# 4) Additionally, when appending parameters from a file, also append an MD5 digest of the file (with included files included
#    inline, comments and leading/trailing whitespace stripped) if the appropriate option is selected in the function
#    call. Include in this digest and XML or HDF5 file in the data/ directory that the file accesses. Also append an MD5 digest of
#    any other files explicitly requested in the "uniqueLabel" directive.
#
# Andrew Benson (10-April-2012)

# Get the source directory.
die('Usage: Make_Unique_Label_Functions.pl <galacticusDirectory>')
    unless ( scalar(@ARGV) == 1 );
(my $galacticusDirectory = $ARGV[0]) =~ s/\/$//;
my $sourceDirectory = $galacticusDirectory."/source";

# Parse the code directive locations file.
my $codeDirectiveLocationsXML = new XML::Simple;
my $codeDirectiveLocations    =$codeDirectiveLocationsXML->XMLin($galacticusDirectory."/work/build/Code_Directive_Locations.xml");

# Open the output file.
open(my $functionHandle,">".$galacticusDirectory."/work/build/utility.input_parameters.unique_labels.inc"             );
open(my $scopeHandle,">".$galacticusDirectory."/work/build/utility.input_parameters.unique_labels.visibilities.inc");

# Store for default parameter values.
my %defaultValues;

# Stores for the structure of new-style method implementations, and which files implement those methods.
my $methodStructure;
my $methodFiles;

# Extract C-include directories from the Makefile.
my $cFlags = $ENV{'GALACTICUS_CFLAGS'};
open(my $make,$galacticusPath."Makefile");
while ( my $line = <$make> ) {
    if ( $line =~ m/CFLAGS\s*\+??=\s*(.*)$/ ) {
	$cFlags .= " ".$1;
    }
}
close($make);
# Extract include directories from flags.
my @includeDirectories;
while ( $cFlags =~ m/\-I(\S+)/ ) {
    my $path = $1;
    $path =~ s/^\.\//$galacticusPath/;
    push(@includeDirectories,$path);
    $cFlags =~ s/\-I(\S+)//;
}

# Iterate over all files containing unique label directives.
foreach my $fileName ( @{$codeDirectiveLocations->{'uniqueLabel'}->{'file'}} ) {
    # Get the equivalent object file.
    (my $objectFile = $fileName) =~ s/\.F90$/.o/;
    # Get the name of the dependencies for this file.
    (my $depFile = $objectFile) =~ s/source\/(.*)\.o$/work\/build\/$1.d/;
    # Extract the function name and any parameters to ignore.
    my ($labelFunction, %ignoreParameters, @ignorePatterns, @hashFiles);
    foreach my $uniqueLabel ( &Directives::Extract_Directives($fileName,"uniqueLabel") ) {
	$labelFunction = $uniqueLabel->{'function'};
	map { $ignoreParameters{$_} = 1 } &ExtraUtils::as_array($uniqueLabel->{'ignore'});
	&ExtraUtils::smart_push(\@ignorePatterns,$uniqueLabel->{'ignoreRegex'});
	@hashFiles = &ExtraUtils::as_array($uniqueLabel->{'hashFile'});
    }
    # Begin creating the function for the definition.
    my $definitionCode =
<<CODE;
function $labelFunction(includeVersion,includeBuild,includeSourceDigest,asHash,parameters)
  implicit none
  type   (varying_string    )                          :: $labelFunction
  logical                    , intent(in   ), optional :: includeVersion,includeBuild,includeSourceDigest,asHash
  type   (varying_string    )                          :: parameterValue
  type   (inputParameterList), intent(  out), optional :: parameters

  $labelFunction=''
CODE
    # Scan dependencies for default parameter values, while also accumulating a full list of dependecy files.
    my @directives;
    my @dependencyFileList;
    my @dependencyFileStack = map {s/\.\/work\/build\/(.*)\.o$/$1/; $_;} split("\n",read_file($depFile));
    while ( scalar(@dependencyFileStack) > 0 ) {
	my $depName = pop(@dependencyFileStack);
	push(@dependencyFileList,$depName);
	$depName =~ s/\.\/work\/build\/(.*)\.o$/$1/;       
	# Scan the file for default parameter values.
	my $sourceFile = $sourceDirectory."/".$depName.".F90";
	unless ( $depName eq "utility.input_parameters" ) {
	    # Extract default values for all parameters defined in this file.
	    map {
		$_->{'submatches'}->[1] =~ s/^\s*//;
		$_->{'submatches'}->[1] =~ s/\s*$//;
		$_->{'submatches'}->[1] =~ s/^'//;
		$_->{'submatches'}->[1] =~ s/'$//;
		$_->{'submatches'}->[1] =~ s/'/''/;
		$defaultValues{$_->{'submatches'}->[0]} = $_->{'submatches'}->[1];
	    } &Fortran_Utils::Get_Matching_Lines($sourceFile,qr/Get_Input_Parameter\s*\(\s*'([^']*)'.*defaultValue\s*=\s*(.*)[,\)]/);
	    # Locate and process any new-style method definitions.
	    foreach my $include ( &Directives::Extract_Directives($sourceFile,"include",conditions => {type => "function"}) ) {
		# Add any files implementing this method to the list of dependencies.
		push(
		    @dependencyFileStack,
		    map
		    {$_ =~ s/^.*$sourceDirectory\/(.*)\.F90/$1/; $_;}
		    &ExtraUtils::as_array($codeDirectiveLocations->{$include->{'directive'}}->{'file'})
		    );
	        # Extract and store the default implementation for this method.
		$defaultValues{$include->{'directive'}."Method"} = $include->{'default'};
		# Add this method to the list of method directives to watch.
		push(@directives,$include->{'directive'});
 	    }
	}		
    }
    # Construct dependencies between new-style method implementations.
    foreach my $directive ( @directives ) {
	$methodStructure->{$directive}->{$directive."Class"}->{'extends'} = undef();
	foreach my $sourceFile ( &ExtraUtils::as_array($codeDirectiveLocations->{$directive}->{'file'}) ) {
	    foreach my $classDeclaration ( &Fortran_Utils::Get_Matching_Lines($sourceFile,$Fortran_Utils::classDeclarationRegEx) ) {
		my $type      = $classDeclaration->{'submatches'}->[3];
		my $extends   = $classDeclaration->{'submatches'}->[1];
		(my $leafName = $sourceFile) =~ s/.*source\///;
		$methodStructure->{$directive}->{$type}->{'extends'} = $extends;
		$methodFiles->{$leafName} = {implementation => $type, directive => $directive};
	    }
	}
    }
    # Build a list of which implementations of new-style methods are actually being used.
    foreach my $directive ( @directives ) {
	foreach my $dependencyFile ( @dependencyFileList ) {
	    my $sourceFile = $sourceDirectory."/".$dependencyFile.".F90";
	    # Detect the specific implementation of any new-style method used.
	    my @directiveNamesUsed;
	    # Construct the selection of all non-base-class implementations of this directive.
	    my $allowed = join("|",grep {$_ ne $directive."Class"} keys(%{$methodStructure->{$directive}}));
	    # Match instantiations of the form:
	    #  type(cosmologyParametersSimple) :: abcd
	    unless ( grep {$_ =~ m/$sourceFile/} &ExtraUtils::as_array($codeDirectiveLocations->{$directive}->{'file'}) ) {
		push(@directiveNamesUsed,$_->{'submatches'}->[0])
		    foreach ( &Fortran_Utils::Get_Matching_Lines($sourceFile,qr/^\s*type\s*\(\s*($allowed)\s*\).*::/i) );
	    }
	    # Match instantiations of the form:
	    #  abcd => cosmologyParametersSimple(1,2,3,4)
	    push(@directiveNamesUsed,$_->{'submatches'}->[0])
		foreach ( &Fortran_Utils::Get_Matching_Lines($sourceFile,qr/^\s*[a-zA-Z0-9_]+\s*=>\s*($allowed)\s*\(.*\)\s*$/i) );
	    # Match instantiations of the form:
	    #  abcd => cosmologyParameters('simple')
	    push(@directiveNamesUsed,$directive.ucfirst($_->{'submatches'}->[0]))
		foreach ( &Fortran_Utils::Get_Matching_Lines($sourceFile,qr/^\s*[a-zA-Z0-9_]+\s*=>\s*$directive\s*\(\s*[\'\"]\s*(.*?)\s*[\'\"]\s*\)\s*$/i) );
	    # Match instantiations of the form:
	    #  abcd => cosmologyParameters()
	    push(@directiveNamesUsed,$directive."Class")
		if ( &Fortran_Utils::Get_Matching_Lines($sourceFile,qr/^\s*[a-zA-Z0-9_]+\s*=>\s*$directive\s*\(\s*\)\s*$/i) );
	    # Record which implementations were used.
	    $methodStructure->{$directive}->{$_}->{'isUsed'} = 1
		foreach ( @directiveNamesUsed );
	}   
	# Where an implementation extends another implementation, ensure that the extended implementation is
	# used if the extension is used.
	foreach my $directive ( @directives ) {
	    foreach ( keys(%{$methodStructure->{$directive}}) ) {
		if ( exists($methodStructure->{$directive}->{$_}->{'isUsed'}) ) {
		    my $parentImplementation = $methodStructure->{$directive}->{$_};
		    while ( defined($parentImplementation->{'extends'}) ) {
			my $parentImplementationName      = $parentImplementation->{'extends'};
			$parentImplementation             = $methodStructure->{$directive}->{$parentImplementationName};
			$parentImplementation->{'isUsed'} = 1;
		    }
		}
	    }
	}
    }
    # Process the list of dependencies.
    foreach my $dependencyFile ( @dependencyFileList ) {
	# Get the name of the associated module.
	my $moduleName    = &Make::Module_Name($dependencyFile,default => "self");
	# Initialize the code that will be used to label this module.
	my $moduleCode    = "  ".$labelFunction."=".$labelFunction."//'::".$moduleName."'\n";
	# Scan this file for parameters.
	my $sourceFile = $sourceDirectory."/".$dependencyFile.".F90";
	if ( -e $sourceFile ) {
	    # Initialize.
	    my $hasParameters = 0;
	    my $methodParameter;
	    my $methodValue;
	    # Extract old-style method activation parameter names and values.
	    map 
	    {$methodParameter = $_->{'submatches'}->[0]."Method"; $methodValue = $_->{'submatches'}->[1];} 
	    &Fortran_Utils::Get_Matching_Lines($sourceFile,qr/^\s*if\s*\(\s*([a-zA-Z0-9_]+)Method\s*==\s*\'([a-zA-Z0-9_\-\+]+)\'\s*\)/);
	    # Extract new-style method activation parameter names and values.
	    foreach my $directive ( @directives ) {
		map 
		{$methodParameter = $directive."Method"; ($methodValue = $_->{'name'}) =~ s/^$directive//; $methodValue = lcfirst($methodValue);} 
		&Directives::Extract_Directives($sourceFile,$directive);
	    }
	    # Extract any input parameters from this file.
	    foreach my $directive ( &Directives::Extract_Directives($sourceFile,"inputParameter",comment => qr/^\s*(\!|\/\/)@/) ) {
		# Use parameter regEx as name if no name is defined.
		$directive->{'name'} = $directive->{'regEx'}
		    unless ( exists($directive->{'name'}) );
		# Ignore parameters that match an ignored name or regex.
		unless ( exists($ignoreParameters{$directive->{'name'}}) || map {$directive->{'name'} =~ m/$_/} @ignorePatterns ) {
		    $hasParameters = 1;
		    # Handle special cases of parameter names that depend on other directives.
		    if ( $directive->{'name'} =~ m/\(\#([a-zA-Z0-9]+)\->([a-z]+)\)/ ) {
			my $dependentDirectiveName = $1;
			my $dependentPropertyName  = $2;
			foreach my $dependentFileName ( @{$codeDirectiveLocations->{$dependentDirectiveName}->{'file'}} ) {
			    foreach my $dependentDirective ( &Directives::Extract_Directives($dependentFileName,$dependentDirectiveName) ) {
				my $dependentProperty = $dependentDirective->{$dependentPropertyName};
				(my $parameterName = $directive->{'name'}) =~ s/\(\#[a-zA-Z0-9]+\->[a-z]+\)/$dependentProperty/;
				$moduleCode .=
<<CODE;
  if (Input_Parameter_Is_Present('$parameterName')) then
   call Get_Input_Parameter_VarString('$parameterName',parameterValue,writeOutput=.false.)
   $labelFunction=$labelFunction//'#$parameterName\['//parameterValue//']'
   if (present(parameters)) then
     call parameters%add("$parameterName",char(parameterValue))
   end if
  end if
CODE
			    }
			}
		    } else {
			# Determine the default value for this parameter.
			my $defaultValue = exists($defaultValues{$directive->{'name'}}) ? $defaultValues{$directive->{'name'}} : "";
			$moduleCode .=
<<CODE;
  call Get_Input_Parameter_VarString('$directive->{'name'}',parameterValue,defaultValue='$defaultValue',writeOutput=.false.)
  $labelFunction=$labelFunction//'#$directive->{'name'}\['//parameterValue//']'
  if (present(parameters)) then
    call parameters%add("$directive->{'name'}",char(parameterValue))
  end if
CODE
		    }
		}
	    }
	    # Add a source MD5 digest for this file.
	    my $ctx = Digest::MD5->new();
	    $ctx->add(&Fortran_Utils::read_file($sourceFile,state => "raw", followIncludes => 1, includeLocations => [ "../source", "../work/build" ], stripRegEx => qr/^\s*![^\#\@].*$/, stripLeading => 1, stripTrailing => 1));
	    # Search for use on any files from the data directory by this source file.
	    &Hash_Data_Files(
		$ctx,
		map 
		{$_->{'submatches'}->[0]} 
		&Fortran_Utils::Get_Matching_Lines($sourceFile,qr/[\"\'](data\/[a-zA-Z0-9_\.\-\/]+\.(xml|hdf5))[\"\']/)
		);
	    my $sourceCodeDigest = $ctx->b64digest();
	    $moduleCode .= "  if (present(includeSourceDigest)) then\n";
	    $moduleCode .= "    if (includeSourceDigest) ".$labelFunction."=".$labelFunction."//'_source:".$dependencyFile."[".$sourceCodeDigest."]'\n";
	    $moduleCode .= "  end if\n";
	    # Assume that we will not test the value of the method parameter by default.
	    my $testMethodParameter = 0;
	    # Test if a method parameter is defined.
	    if ( defined($methodParameter) ) {
		# A method parameter is defined, determine if this file implements a new-style method.
		(my $leafName = $sourceFile) =~ s/.*source\///;
		if ( exists($methodFiles->{$leafName}) ) {	
		    # It does implement a new-style method, extract the method and implementation.
		    my $directive          = $methodFiles->{$leafName}->{'directive'     };
		    my $implementationName = $methodFiles->{$leafName}->{'implementation'};
		    # We must test the value of the method parameter iff the specific implementation is not used, but the default
		    # implementation is used.
		    $testMethodParameter = 1
			if (
			      exists($methodStructure->{$directive}->{$directive."Class" }->{'isUsed'})
			    &&
			    ! exists($methodStructure->{$directive}->{$implementationName}->{'isUsed'})
			);
		    # If neither the default or specific implementation is used, we do not want to include the parameters from
		    # this implementation.
		    $hasParameters = 0
			unless ( 
			    exists($methodStructure->{$directive}->{$directive."Class" }->{'isUsed'}) 
			    ||
			    exists($methodStructure->{$directive}->{$implementationName}->{'isUsed'})
			);
		} else {
		    # This file implements an old-style method - always test the value of the method parameter in this case.
		    $testMethodParameter = 1;
		}
	    }
	    if ( $hasParameters == 1 ) {
		if ( $testMethodParameter == 1 ) {
		    my $defaultValue  = exists($defaultValues{$methodParameter}) ? $defaultValues{$methodParameter} : "";
 		    $definitionCode  .=
<<CODE;
  call Get_Input_Parameter_VarString('$methodParameter',parameterValue,defaultValue='$defaultValue',writeOutput=.false.)
  if (present(parameters)) then
    call parameters%add("$methodParameter",char(parameterValue))
  end if
  if (parameterValue == '$methodValue') then
$moduleCode  end if
CODE
		} else {
		    $definitionCode .= $moduleCode;
		}
	    }
	} else {
	    # A matching Fortran source file was not found. In this case, look for a matching file of other type, and include a
	    # suitable digest in the function.
	    foreach my $suffix ( ".c", ".cpp" ) {
		my $sourceFile = $sourceDirectory."/".$dependencyFile.$suffix;
		if ( -e $sourceFile ) {
		    my $ctx = Digest::MD5->new();
		    my $sourceHandle;
		    if ( $suffix eq ".c" || $suffix eq ".cpp" ) {
			# For C and C++ files, run them through the preprocessor to have any include files included.
			my $includeOptions = join(" ",map {"-I".$_} @includeDirectories);
			open($sourceHandle,"cpp -I".$galacticusDirectory."/source/ -I".$galacticusDirectory."/work/build/ ".$includeOptions." ".$sourceFile."|");
		    } else {
			# For other files, simply read the file directly.
			open($sourceHandle,       $sourceFile    );
		    }
		    $ctx->addfile($sourceHandle);
		    close($sourceHandle);
		    # Search for use on any files from the data directory by this source file.
		    &Hash_Data_Files(
			$ctx,
			map 
			{$_->{'submatches'}->[0]} 
			&File_Match::Get_Matching_Lines($sourceFile,qr/\"(data\/[a-zA-Z0-9_\.\-\/]+\.(xml|hdf5))\"/)
			);
		    my $sourceCodeDigest = $ctx->b64digest();
		    $definitionCode .=
<<CODE;
  if (present(includeSourceDigest)) then
    if (includeSourceDigest) $labelFunction=$labelFunction//'_source:$dependencyFile\[$sourceCodeDigest\]'
  end if
CODE
		}
	    }
	}
    }  
    # Add a hash of any other files which we were explicitly instructed to include as dependencies.
    if ( scalar(@hashFiles) > 0 ) {
	my $ctx = Digest::MD5->new();
	&Hash_Data_Files($ctx,@hashFiles);
	my $sourceCodeDigest = $ctx->b64digest();
	$definitionCode .=
<<CODE;
  if (present(includeSourceDigest)) then
    if (includeSourceDigest) $labelFunction=$labelFunction//'_extraFiles\[$sourceCodeDigest\]'
  end if
CODE
    }
    # Finish the definition code.
    $definitionCode .=
<<CODE;
  if (present(includeVersion    )) then
    if (includeVersion     ) $labelFunction=$labelFunction//'_'//Galacticus_Version     ()
  end if
  if (present(includeBuild       )) then
    if (includeBuild       ) $labelFunction=$labelFunction//'_'//Galacticus_Build_String()
  end if
  if (present(asHash)) then
    if (asHash) $labelFunction=Hash_MD5($labelFunction)
  end if
  return
end function $labelFunction

CODE
    # Write the function definition to file.
    print $functionHandle $definitionCode;
    print $scopeHandle "public :: ".$labelFunction."\n";
}
# Close the output files.
close($functionHandle);
close($scopeHandle);
# Done.
exit;

sub Hash_Data_Files {
    # Run a supplied list of files through a supplied MD5 hash object.
    my $ctx   = shift;
    my @files = @_;
    foreach ( @files ) {
    	# Run each data file through the MD5 hash.
    	my $dataFileName = $galacticusDirectory."/".$_;
    	if ( -e $dataFileName ) {
    	    open(my $dataHandle,$dataFileName);
    	    $ctx->addfile($dataHandle);
    	    close($dataHandle);
    	}
    }
}
