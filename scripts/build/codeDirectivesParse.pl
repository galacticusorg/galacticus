#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Fcntl qw(SEEK_SET);
use XML::Simple;
use Data::Dumper;
use File::Changes;
use Galacticus::Build::Directives;
use List::ExtraUtils;
use List::Uniq qw(uniq);
use Storable;

# Scans source code for "!![...!!]" directives and generates a Makefile.
# Andrew Benson (09-September-2016)

# Read command line arguments.
die "Usage: codeDirectivesParse.pl <installDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $installDirectoryName = $ARGV[0];
# Specify verbosity.
my $verbosity            = 0;
# Get an XML parser.
my $xml                  = new XML::Simple();
# Initialize data structure to hold per-file information.
my $directivesPerFile;
my $havePerFile = -e $ENV{'BUILDPATH'}."/codeDirectives.blob";
my $updateTime;
if ( $havePerFile ) {
    $directivesPerFile = retrieve($ENV{'BUILDPATH'}."/codeDirectives.blob");
    $updateTime        = -M       $ENV{'BUILDPATH'}."/codeDirectives.blob" ;
}
# Open the source directory and get a list of source files to process.
my $sourceDirectoryName = $installDirectoryName."/source";
opendir(my $sourceDirectory,$sourceDirectoryName) or die "Can't open the source directory: #!";
my @sourceFileNames = grep {$_ =~ m/\.(f|f90|c|cpp|h)$/i && $_ !~ m/^\.\#/} readdir($sourceDirectory);
closedir($sourceDirectory);

# Look for changes in source files. Force a rescan of all files if anything has changed.
my $forceRescan = 0;
my @fileIdentifiers;
foreach my $sourceFileName ( @sourceFileNames ) {
    (my $fileIdentifier = $sourceDirectoryName."/".$sourceFileName) =~ s/\//_/g;
    push(@fileIdentifiers,$fileIdentifier);
}
# Check for new files.
foreach my $fileIdentifier ( @fileIdentifiers ) {
    unless ( exists($directivesPerFile->{$fileIdentifier}) ) {
	$forceRescan = 1;
    }
}
# Check for removed files.
foreach my $fileIdentifier ( keys(%{$directivesPerFile}) ) {
    unless ( grep {$_ eq $fileIdentifier} @fileIdentifiers ) {
	$forceRescan = 1;
    }
}

# Iterate over source files.
foreach my $fileName ( @sourceFileNames ) {

    # Add the file to the stack of filenames to process.
    my @fileNames = ( $sourceDirectoryName."/".$fileName );    
    (my $fileIdentifier = $sourceDirectoryName."/".$fileName) =~ s/\//_/g;
    $fileIdentifier =~ s/^\._??//;

    # Check if file is updated. If it is not, skip processing it. If it is, remove previous record of uses and rescan.
    my $rescan = 1;
    if ( $havePerFile && exists($directivesPerFile->{$fileIdentifier}) ) {
	$rescan = 0
	    unless ( grep {-M $_ < $updateTime} &List::ExtraUtils::as_array($directivesPerFile->{$fileIdentifier}->{'files'}) );
    }
    if ( $rescan || $forceRescan ) {
	delete($directivesPerFile->{$fileIdentifier})
    	    if ( $havePerFile && exists($directivesPerFile->{$fileIdentifier}) );
	push(@{$directivesPerFile->{$fileIdentifier}->{'files'}},$sourceDirectoryName."/".$fileName);
	# Process files until none remain.
	while ( scalar(@fileNames) > 0 ) {
	    my $filePathName = pop(@fileNames);
	    # Push any included files onto the stack.
	    open(my $fileHandle,$filePathName) 
		or die "codeDirectivesParse.pl: can not open input file: #!";
	    my @includedFiles = 
		map
	    {
		$_ =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/
		    ?
		    do {(my $includeFileName = $sourceDirectoryName."/".$1) =~ s/\.inc$/.Inc/; -e $includeFileName ? $includeFileName : ()}
		:
		    ()
	    }
	    <$fileHandle>;
	    close($fileHandle);
	    push(@fileNames,@includedFiles);
	    push(@{$directivesPerFile->{$fileIdentifier}->{'files'}},@includedFiles);
	    # Get all directives in the file.
	    foreach my $directive ( &Galacticus::Build::Directives::Extract_Directives($filePathName,"*", setRootElementType => 1) ) {
		# Act on the directive. "Include" directives are handled separately from other directives.
		if ( $directive->{'rootElementType'} eq "include" ) {
		    # Store the source file name for this directive.
		    $directive->{'source'} = $filePathName;
		    # Determine the file name to be included, store it, then remove the include statement from the directive.
		    ($directive->{'fileName'} = $ENV{'BUILDPATH'}."/".$1) =~ s/\.inc$/\.Inc/
			if ( $directive->{'content'} =~ m/^\s*\#??include\s*["'<](.+)["'>]/im );
		    delete($directive->{'content'});
		    # Create an entry for this directive in the list of include directives.
		    my $xmlOutput                        = new XML::Simple( NoAttr => 1, RootName => $directive->{'rootElementType'} );
		    my $directiveName                    = (exists($directive->{'name'}) ? $directive->{'name'} : $directive->{'directive'}).".".$directive->{'type'};
		    $directivesPerFile->{$fileIdentifier}->{'includeDirectives'}->{$directiveName} = 
		    {
			source   => $filePathName,
			fileName => $directive->{'fileName'},
			xml      => $xmlOutput->XMLout($directive)
		    };
		} else {
		    # For non-include directives, simply record the file which originated the directive.
		    push(@{$directivesPerFile->{$fileIdentifier}->{'nonIncludeDirectives'}->{$directive->{'rootElementType'}}->{'files'}},$sourceDirectoryName."/".$fileName);
		    if ( $directive->{'rootElementType'} eq "functionClass" ) {
			# For functionClass directives, separately store the name of the file originating the directive, and add
			# implicit directives originating from the preprocessed source file.
			(my $fileName = $ENV{'BUILDPATH'}."/".$fileName) =~ s/\.F90$/.p.F90/;
			$directivesPerFile->{$fileIdentifier}->{'functionClasses'}->{$directive->{'name'}} = $fileName;
			&addImplicitDirectives($directive,$directivesPerFile->{$fileIdentifier}->{'nonIncludeDirectives'},$fileIdentifier,$fileName,$fileName);
		    }
		}
	    }
	}
    }
}
# Reduce over files.
my $includeDirectives;
my $nonIncludeDirectives;
my $functionClasses;
foreach my $fileIdentifier ( keys(%{$directivesPerFile}) ) {
    if ( exists($directivesPerFile->{$fileIdentifier}->{'includeDirectives'}) ) {
	$includeDirectives->{$_} = $directivesPerFile->{$fileIdentifier}->{'includeDirectives'}->{$_}
            foreach ( keys(%{$directivesPerFile->{$fileIdentifier}->{'includeDirectives'}}) );
    }
    if ( exists($directivesPerFile->{$fileIdentifier}->{'nonIncludeDirectives'}) ) {
	foreach my $directive ( keys(%{$directivesPerFile->{$fileIdentifier}->{'nonIncludeDirectives'}}) ) {
	    foreach my $content ( 'files', 'dependency' ) {
		push(@{$nonIncludeDirectives->{$directive}->{$content}},@{$directivesPerFile->{$fileIdentifier}->{'nonIncludeDirectives'}->{$directive}->{$content}})
		    if ( exists($directivesPerFile->{$fileIdentifier}->{'nonIncludeDirectives'}->{$directive}->{$content}) );
	    }
	}
    }
    if ( exists($directivesPerFile->{$fileIdentifier}->{'functionClasses'}) ) {
	$functionClasses->{$_} = $directivesPerFile->{$fileIdentifier}->{'functionClasses'}->{$_}
            foreach ( keys(%{$directivesPerFile->{$fileIdentifier}->{'functionClasses'}}) );
    }
}
# Uniqueify.
foreach my $directive ( keys(%{$nonIncludeDirectives}) ) {
    foreach my $content ( 'files', 'dependency' ) {
	@{$nonIncludeDirectives->{$directive}->{$content}} = uniq(sort(@{$nonIncludeDirectives->{$directive}->{$content}}))
	    if ( exists($nonIncludeDirectives->{$directive}->{$content}) );
    }
}
# Output a file listing which files contain which directives.
my $outputDirectives;
@{$outputDirectives->{$_}->{'file'}} = @{$nonIncludeDirectives->{$_}->{'files'}}
    foreach ( keys(%{$nonIncludeDirectives}) );
my $xmlOutput = new XML::Simple( NoAttr => 1, RootName => "directives" );
open(my $directiveLocationsFile,">".$ENV{'BUILDPATH'}."/directiveLocations.xml.tmp");
print $directiveLocationsFile $xmlOutput->XMLout($outputDirectives);
close($directiveLocationsFile);
&File::Changes::Update($ENV{'BUILDPATH'}."/directiveLocations.xml",$ENV{'BUILDPATH'}."/directiveLocations.xml.tmp");
# Output the Makefile and XML files containing the directives.
open(my $directivesMakefile,">".$ENV{'BUILDPATH'}."/Makefile_Directives");
foreach my $directive ( sort(keys(%{$includeDirectives})) ) {
    # Construct the unpreprocessed include file name.
    (my $fileName = $includeDirectives->{$directive}->{'fileName'}) =~ s/\.inc$/\.Inc/;
    # Build a list of extra dependencies.
    my @extraDependencies;
    # For "function" actions, add the file containing the directive as an additional dependency
    # as these files are simply copied into the include file as part of the include file construction.
    push(@extraDependencies,sort(keys(%{$nonIncludeDirectives->{$1}->{'files'}})))
	if ( $directive =~ m/^([a-zA-Z0-9_]+)\.function$/ );
    # Add on any other dependencies.
    push(@extraDependencies,@{$nonIncludeDirectives->{$1}->{'dependency'}})
	if ( $directive =~ m/^([a-zA-Z0-9_]+)\.(moduleUse|functionCall)$/ &&  exists($nonIncludeDirectives->{$1}->{'dependency'}) );
    # Output the Makefile rule. Add an explicit dependence on "hdf5FCInterop.dat" to ensure that HDF5 type interoperability is
    # determined before attempting to build the file. Similarly add explicit dependence on "openMPCriticalSections.xml" so that
    # OpenMP critical sections can be referenced by name.
    print $directivesMakefile $fileName.".up: ".$ENV{'BUILDPATH'}."/".$directive.".xml ".join(" ",@extraDependencies)." \$(BUILDPATH)/hdf5FCInterop.dat \$(BUILDPATH)/openMPCriticalSections.xml\n";
    print $directivesMakefile "\t./scripts/build/buildCode.pl ".$installDirectoryName." ".$ENV{'BUILDPATH'}."/".$directive.".xml\n";
    print $directivesMakefile $fileName.": ".$fileName.".up\n";
    print $directivesMakefile "\n";
    # Output the directive itself.
    open(my $directiveFile,">".$ENV{'BUILDPATH'}."/".$directive.".xml.tmp");
    print $directiveFile $includeDirectives->{$directive}->{'xml'};
    close($directiveFile);
    &File::Changes::Update($ENV{'BUILDPATH'}."/".$directive.".xml",$ENV{'BUILDPATH'}."/".$directive.".xml.tmp");
}
# Add additional dependencies for object files of source files that contain functionClass directives. These source files get other
# source files incorporated into them via the source tree preprocessor.
foreach my $directiveName ( sort(keys(%{$functionClasses})) ) {
    print $directivesMakefile $functionClasses->{$directiveName}.".up: ".join(" ",sort(@{$nonIncludeDirectives->{$directiveName}->{'files'}}))."\n\n";
}
# Include explicit dependencies for Makefile_Use_Dependencies to ensure that module dependencies get rebuilt
# after these directive include files are constructed.
print $directivesMakefile $ENV{'BUILDPATH'}."/Makefile_Use_Dependencies: ".join(" ",sort(map {(my $fileName = $includeDirectives->{$_}->{'fileName'}) =~ s/\.inc$/\.Inc/; $fileName} keys(%{$includeDirectives})))."\n\n";
# Include a rule for including Makefile_Component_Includes. This has to go here since Makefile_Component_Includes depends on
# objects.nodes.components.Inc for which Makefile_Directive contains the build rule.
print $directivesMakefile "-include ".$ENV{'BUILDPATH'}."/Makefile_Component_Includes\n";
print $directivesMakefile $ENV{'BUILDPATH'}."/Makefile_Component_Includes: ".$ENV{'BUILDPATH'}."/objects.nodes.components.Inc\n\n";
close($directivesMakefile);
# Output the per file directive data.
store($directivesPerFile,$ENV{'BUILDPATH'}."/codeDirectives.blob.tmp");
&File::Changes::Update($ENV{'BUILDPATH'}."/codeDirectives.blob",$ENV{'BUILDPATH'}."/codeDirectives.blob.tmp");
exit;

sub addImplicitDirectives {
    # Add implicit dependencies required by certain directives.
    my $directive          = shift();
    my $implicitDirectives = shift();
    my $fileIdentifier     = shift();
    my $fileName           = shift();
    my $filePathName       = shift();
    my %implicitDirectives =
	(
	 stateful         => 
	 {
	     always => $directive->{'rootElementType'} eq "functionClass"           ,
	     tasks  => [ "galacticusStateRetrieveTask", "galacticusStateStoreTask" ]
	 },
	 functionClassDestroy => 
	 {
	     always => $directive->{'rootElementType'} eq "functionClass" && ! exists($directive->{'functionClassDestroy'}),
	     tasks  => [ "functionClassDestroyTask"                                ]
	 }
	);
    foreach my $implicitDirective ( &List::ExtraUtils::hashList(\%implicitDirectives, keyAs => "directive") ) {
	if ( ( exists($directive->{$implicitDirective->{'directive'}}) && $directive->{$implicitDirective->{'directive'}} eq "yes" ) || $implicitDirective->{'always'} ) {
	    foreach my $implicitDependency ( @{$implicitDirective->{'tasks'}} ) {
		push(@{$directivesPerFile->{$fileIdentifier}->{'nonIncludeDirectives'}->{$implicitDependency}->{'files'     }},$filePathName);
		push(@{$directivesPerFile->{$fileIdentifier}->{'nonIncludeDirectives'}->{$implicitDependency}->{'dependency'}},$fileName    );
	    }			    
	}
    }
}
