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
use Scalar::Util 'reftype';
use Data::Dumper;
require Fortran::Utils;

# Construct a unique string that defines the operation of a given module.
# Andrew Benson (10-April-2012)

# Get the source directory.
die('Usage: Make_Unique_Label_Functions.pl <sourceDir>') unless ( scalar(@ARGV) == 1 );
my $sourceDirectory = $ARGV[0]."/source";

# Parse the code directive locations file.
my $codeDirectiveLocationsXML = new XML::Simple;
my $codeDirectiveLocations    =$codeDirectiveLocationsXML->XMLin($ARGV[0]."/work/build/Code_Directive_Locations.xml");

# Open the output file.
open(oHndl,">./work/build/utility.input_parameters.unique_labels.inc"             );
open(vHndl,">./work/build/utility.input_parameters.unique_labels.visibilities.inc");

# Store for default parameter values.
my %defaultValues;

# Set of directives that we're watching.
my @directives;

# Scan the source directory for source files.
opendir(dHndl,$sourceDirectory);
while ( my $fileName = readdir(dHndl) ) {
    # Initialize.
    my $processFile = 0;
    my $labelFunction;
    my %ignoreParameters;
    my @ignorePatterns;

    # Select Fortran source files.
    if ( $fileName =~ m/\.F90$/ ) {
	my $xmlBuffer;
	open(iHndl,$sourceDirectory."/".$fileName);
	while ( my $xmlLine = <iHndl> ) {
	    if ( $xmlLine =~ m/^\s*\!#(.*)/ ) {
		my $xmlLine = $1;
		$xmlBuffer = "" if ( $xmlLine =~ m/^\s*<(uniqueLabel)>\s*$/ );
		$xmlBuffer .= $xmlLine;
		if ( $xmlLine =~ m/^\s*<\/uniqueLabel>\s*$/ ) {
		    # Parse the XML.
		    my $xml = new XML::Simple;
		    my $uniqueLabel = $xml->XMLin($xmlBuffer);
		    $labelFunction = $uniqueLabel->{'function'};
		    if ( exists($uniqueLabel->{'ignore'}) && reftype($uniqueLabel->{'ignore'}) ) {
			my @ignores;
			if ( reftype($uniqueLabel->{'ignore'}) eq "ARRAY" ) {
			    @ignores = @{$uniqueLabel->{'ignore'}};
			} else {
			    push(@ignores,$uniqueLabel->{'ignore'});
			}
			foreach my $ignore ( @ignores ) {
			    $ignoreParameters{$ignore} = 1;
			}
		    }
		    if ( exists($uniqueLabel->{'ignoreRegex'}) ) {
			if ( UNIVERSAL::isa($uniqueLabel->{'ignoreRegex'},"ARRAY") ) {
			    @ignorePatterns = @{$uniqueLabel->{'ignoreRegex'}};
			} else {
			    push(@ignorePatterns,$uniqueLabel->{'ignoreRegex'});
			}
		    }
		    $processFile = 1;
		}
	    }
	}
	close(iHndl);
    }

    # Process the file if necessary.
    if ( $processFile == 1 ) {
	# Initialize the definition code.
	my $definitionCode;
	
	# Get the equivalent object file.
	my $objectFile = "./work/build/".$fileName;
	$objectFile =~ s/\.F90$/.o/;

	# Get the name of the module supplied by this file.
	my $moduleFile = $objectFile;
	$moduleFile =~ s/\.o$/.m/;
	my $depFile = $objectFile;
	$depFile =~ s/\.o$/.d/;
	open(iHndl,$moduleFile);
	unless ( eof(iHndl) ) {
	    my $selfName = <iHndl>;
	    chomp($selfName);
	    close(iHndl);
	    $selfName =~ s/\.\/work\/build\/(.*)\.mod$/$1/;
	    
	    # Begin creating the function for the definition.
	    $definitionCode .= "function ".$labelFunction."(includeVersion,asHash)\n";
	    $definitionCode .= "  implicit none\n";
	    $definitionCode .= "  type(varying_string)                       :: ".$labelFunction."\n";
	    $definitionCode .= "  logical             , intent(in), optional :: includeVersion,asHash\n";
	    $definitionCode .= "  type(varying_string)                       :: parameterValue\n";
	    $definitionCode .= "  ".$labelFunction."=''\n";
	    
	    # Build a stack of dependency files.
	    my @initialDependencyFiles;
	    open(iHndl,$depFile);
	    while ( my $depName = <iHndl> ) {
		chomp($depName);
		$depName =~ s/\.\/work\/build\/(.*)\.o$/$1/;
		push(@initialDependencyFiles,$depName);
	    }
	    # Scan dependencies for default parameter values.
	    my @depFiles = @initialDependencyFiles;
	    while ( scalar(@depFiles) > 0 ) {
		my $depName = pop(@depFiles);
		chomp($depName);
		$depName =~ s/\.\/work\/build\/(.*)\.o$/$1/;       
		# Scan the file for default parameter values.
		my $sourceFile = $sourceDirectory."/".$depName.".F90";
		unless ( $depName eq "utility.input_parameters" ) {
		    my $fileHandle;
		    my $methodXmlBuffer;
		    open($fileHandle,$sourceFile);
		    until ( eof($fileHandle) ) {
			# Grab the next Fortran line.
			my $rawLine;
			my $processedLine;
			my $bufferedComments;
			&Fortran_Utils::Get_Fortran_Line($fileHandle,$rawLine,$processedLine,$bufferedComments);
			if ( $processedLine =~ m/Get_Input_Parameter/i ) {
			    if ( $processedLine =~ m/defaultValue\s*=\s*(.*)[,\)]/i ) {
				my $defaultValue = $1;
				$defaultValue =~ s/^\s*'//;
				$defaultValue =~ s/'\s*$//;
				$defaultValue =~ s/'/''/;
				$defaultValue =~ s/^\s*//;
				$defaultValue =~ s/\s*$//;
				if ( $processedLine =~ m/Get_Input_Parameter\s*\(\s*'([^']*)'/i ) {
				    my $parameterName = $1;
				    $defaultValues{$parameterName} = $defaultValue;
				}
			    }
			}
			# Add files contains new-style method implementations to the stack.
			if ( $bufferedComments =~ m/^\#(.*)/ ) {
			    my $xmlLine = $1;
			    $methodXmlBuffer = "" if ( $xmlLine =~ m/^\s*<include[\s>]/ );
			    $methodXmlBuffer .= $xmlLine;
			    if ( $xmlLine =~ m/^\s*<\/include>\s*$/ ) {
				# Parse the XML.
				my $xml = new XML::Simple;
				my $method = $xml->XMLin($methodXmlBuffer);
				if ( $method->{'type'} eq "function" ) {
				    # Find the files which implement this function.
				    my @implementations;
				    if ( UNIVERSAL::isa($codeDirectiveLocations->{$method->{'directive'}}->{'file'},"ARRAY") ) {
					push(@implementations,@{$codeDirectiveLocations->{$method->{'directive'}}->{'file'}});
				    } else {
					push(@implementations,  $codeDirectiveLocations->{$method->{'directive'}}->{'file'} );
				    }
				    $_ =~ s/^.*$sourceDirectory\/(.*)\.F90/$1/
					foreach ( @implementations );
				    push(@depFiles,@implementations);
				}
			    }
			}
		    }
		    close($fileHandle);
		}		
	    }
	    close(iHndl);
	    # Process the list of dependencies.
	    @depFiles = @initialDependencyFiles;
	    while ( scalar(@depFiles) > 0 ) {
		# Pop a dependency file off the stack.
		my $depName = pop(@depFiles);	    
		# Get the name of the associated module.
		my $moduleFile = "./work/build/".$depName.".m";
		my $moduleName = $depName;
		open(jHndl,$moduleFile);
		unless ( eof(jHndl) ) {
		    $moduleName = <jHndl>;
		    chomp($moduleName);
		    $moduleName =~ s/\.\/work\/build\/(.*)\.mod$/$1/;
		}
		close(jHndl);
		my $moduleCode    = "  ".$labelFunction."=".$labelFunction."//'::".$moduleName."'\n";	
		my $hasParameters = 0;
		# Scan this file for parameters.
		my $methodParameter;
		my $methodValue;
		my $xmlBuffer;
		my $methodXmlBuffer;
		my $sourceFile = $sourceDirectory."/".$depName.".F90";
		if ( -e $sourceFile ) {
		    open(sFile,$sourceFile);
		    while ( my $line = <sFile> ) {
			# Find old-style method activations.
			if ( $line =~ m/^\s*if\s*\(\s*([a-zA-Z0-9_]+)Method\s*==\s*\'([a-zA-Z0-9_\-\+]+)\'\s*\)/ ) {
			    $methodParameter = $1."Method";
			    $methodValue     = $2;
			}
			# Find new-style method activations.
			foreach my $directive ( @directives ) {
			    if ( $line =~ m/\s*(\!|\/\/)\#\s*<$directive\s+name\s*=\s*\"$directive(.*)\"\s*\/>/ ) {
				$methodParameter = $directive."Method";
				$methodValue     = lcfirst($2);
			    }
			}
			# Add files contains new-style method implementations to the stack.
			if ( $line =~ m/^\s*(\!|\/\/)\#(.*)/ ) {
			    my $xmlLine = $2;
			    $methodXmlBuffer = "" if ( $xmlLine =~ m/^\s*<include[\s>]/ );
			    $methodXmlBuffer .= $xmlLine;
			    if ( $xmlLine =~ m/^\s*<\/include>\s*$/ ) {
				# Parse the XML.
				my $xml = new XML::Simple;
				my $method = $xml->XMLin($methodXmlBuffer);
				if ( $method->{'type'} eq "function" ) {
				    # Find the files which implement this function.
				    my @implementations;
				    if ( UNIVERSAL::isa($codeDirectiveLocations->{$method->{'directive'}}->{'file'},"ARRAY") ) {
					push(@implementations,@{$codeDirectiveLocations->{$method->{'directive'}}->{'file'}});
				    } else {
					push(@implementations,  $codeDirectiveLocations->{$method->{'directive'}}->{'file'} );
				    }
				    $_ =~ s/^.*$sourceDirectory\/(.*)\.F90/$1/
					foreach ( @implementations );
				    push(@depFiles,@implementations);
				    push(@directives,$method->{'directive'});
				    $defaultValues{$method->{'directive'}."Method"} = $method->{'default'};
				}
			    }
			}
			# Find XML blobs.
			if ( $line =~ m/^\s*(\!|\/\/)@(.*)/ ) {
			    my $xmlLine = $2;
			    $xmlBuffer = "" if ( $xmlLine =~ m/^\s*<(inputParameter)>\s*$/ );
			    $xmlBuffer .= $xmlLine;
			    if ( $xmlLine =~ m/^\s*<\/inputParameter>\s*$/ ) {
				# Parse the XML.
				my $xml = new XML::Simple;
				my $inputParameter = $xml->XMLin($xmlBuffer);
				my $ignore = 0;
				$ignore = 1
				    if ( exists($ignoreParameters{$inputParameter->{'name'}}) );
				foreach ( @ignorePatterns ) {
				    $ignore = 1
					if ( $inputParameter->{'name'} =~ m/$_/ );
				}
				unless ( $ignore == 1 ) {
				    $moduleCode .= "  call Get_Input_Parameter_VarString('".$inputParameter->{'name'}."',parameterValue";
				    my $defaultValue = "";
				    $defaultValue = $defaultValues{$inputParameter->{'name'}}
				    if ( exists($defaultValues{$inputParameter->{'name'}}) );
				    $moduleCode .= ",defaultValue='".$defaultValue."'";
				    $moduleCode .= ",writeOutput=.false.)\n";
				    $moduleCode .= "  ".$labelFunction."=".$labelFunction."//'#".$inputParameter->{'name'}."['//parameterValue//']'\n";
				    $hasParameters = 1;
				}
			    }
			}
		    }
		    close(sFile);
		}
		if ( $hasParameters == 1 ) {
		    if ( defined($methodParameter) ) {
			$definitionCode .= "  call Get_Input_Parameter_VarString('".$methodParameter."',parameterValue";
			my $defaultValue = "";
			$defaultValue = $defaultValues{$methodParameter}
			if ( exists($defaultValues{$methodParameter}) );
			$definitionCode .= ",defaultValue='".$defaultValue."'";
			$definitionCode .= ",writeOutput=.false.)\n";
		    }
		    $definitionCode .= "  if (parameterValue == '".$methodValue."') then\n"
			if ( defined($methodParameter) );
		    $definitionCode .= $moduleCode;
		    $definitionCode .= "  end if\n"
			if ( defined($methodParameter) );	
		}
	    }
	    close(iHndl);
	    
	    # Finish the definition code.
	    $definitionCode .= "  if (present(includeVersion)) then\n";
	    $definitionCode .= "    if (includeVersion) ".$labelFunction."=".$labelFunction."//'_'//Galacticus_Version()\n";
	    $definitionCode .= "  end if\n";
	    $definitionCode .= "  if (present(asHash)) then\n";
	    $definitionCode .= "    if (asHash) ".$labelFunction."=Hash_MD5(".$labelFunction.")\n";
	    $definitionCode .= "  end if\n";
	    $definitionCode .= "  return\n";
	    $definitionCode .= "end function ".$labelFunction."\n\n";
	    
	    # Write the function definition to file.
	    print oHndl $definitionCode;
	    print vHndl "public :: ".$labelFunction."\n";
	} else {
	    close(iHndl);
	}	
    }
}
closedir(dHndl);

# Close the output files.
close(oHndl);
close(vHndl);

exit;
