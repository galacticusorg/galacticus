#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Fcntl qw(SEEK_SET);
use XML::Simple;
use Data::Dumper;
use File::Changes;
use Galacticus::Build::Directives;

# Scans source code for "!#" directives and generates a Makefile.
# Andrew Benson (09-September-2016)

# Read command line arguments.
die "Usage: codeDirectivesParse.pl <installDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $installDirectoryName = $ARGV[0];
# Specify verbosity.
my $verbosity            = 0;
# Get an XML parser.
my $xml                  = new XML::Simple();
# Initialize hashes.
my $includeDirectives;
my $nonIncludeDirectives;
my $functionClasses;
# Open the source directory and get a list of source files to process.
my $sourceDirectoryName = $installDirectoryName."/source";
opendir(my $sourceDirectory,$sourceDirectoryName) or die "Can't open the source directory: #!";
my @sourceFileNames = grep {$_ =~ m/\.(f|f90|c|cpp|h)$/i && $_ !~ m/^\.\#/} readdir($sourceDirectory);
closedir($sourceDirectory);
# Iterate over source files.
foreach my $fileName ( @sourceFileNames ) {
    # Add the file to the stack of filenames to process.
    my @fileNames = ( $sourceDirectoryName."/".$fileName );    
    # Process files until none remain.
    while ( scalar(@fileNames) > 0 ) {
	my $filePathName = pop(@fileNames);
	# Push any included files onto the stack.
	open(my $fileHandle,$filePathName) 
	    or die "codeDirectivesParse.pl: can not open input file: #!";
	push(
	    @fileNames,
	    map
	    {
		$_ =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/
	        ?
	        do {(my $includeFileName = $sourceDirectoryName."/".$1) =~ s/\.inc$/.Inc/; -e $includeFileName ? $includeFileName : ()}
		:
	        ()
	    }
	    <$fileHandle>
	    );
	close($fileHandle);
	# Get all directives in the file.	
	foreach my $directive ( &Galacticus::Build::Directives::Extract_Directives($filePathName,"*", comment => qr/^\s*(!|\/\/)\#\s+/, setRootElementType => 1) ) {
	    # Act on the directive. "Include" directives are handled separately from other directives.
	    if ( $directive->{'rootElementType'} eq "include" ) {
		# Store the source file name for this directive.
		$directive->{'source'} = $filePathName;
		# Determine the file name to be included, store it, then remove the include statement from the directive.
		($directive->{'fileName'} = $ENV{'BUILDPATH'}."/".$1) =~ s/\.inc$/\.Inc/
		    if ( $directive->{'content'} =~ m/^\s*\#??include\s*["'<](.+)["'>]/i );
		delete($directive->{'content'});
		# Create an entry for this directive in the list of include directives.
		my $xmlOutput                        = new XML::Simple( NoAttr => 1, RootName => $directive->{'rootElementType'} );
		my $directiveName                    = (exists($directive->{'name'}) ? $directive->{'name'} : $directive->{'directive'}).".".$directive->{'type'};
		$includeDirectives->{$directiveName} = 
		{
		    source   => $filePathName,
		    fileName => $directive->{'fileName'},
		    xml      => $xmlOutput->XMLout($directive)
		};
		# Add implicit directives for function directives.
		&addImplicitDirectives($directive,$nonIncludeDirectives,$fileName,$sourceDirectoryName."/".$fileName)
		    if ( $directive->{'type'} eq "function" );
	    } else {
		# For non-include directives, simply record the file which originated the directive.
		$nonIncludeDirectives->{$directive->{'rootElementType'}}->{'files'}->{$sourceDirectoryName."/".$fileName} = 1;
		if ( $directive->{'rootElementType'} eq "functionClass" ) {
		    # For functionClass directives, separately store the name of the file originating the directive, and add
		    # implicit directives originating from the preprocessed source file.
		    (my $fileName = $ENV{'BUILDPATH'}."/".$fileName) =~ s/\.F90$/.p.F90/;
		    $functionClasses->{$directive->{'name'}} = $fileName;
		    &addImplicitDirectives($directive,$nonIncludeDirectives,$fileName,$fileName);
		}
	    }
	}
    }
}
# Output a file listing which files contain which directives.
my $outputDirectives;
@{$outputDirectives->{$_}->{'file'}} = keys(%{$nonIncludeDirectives->{$_}->{'files'}})
    foreach ( keys(%{$nonIncludeDirectives}) );
my $xmlOutput = new XML::Simple( NoAttr => 1, RootName => "directives" );
open(my $directiveLocationsFile,">".$ENV{'BUILDPATH'}."/directiveLocations.xml");
print $directiveLocationsFile $xmlOutput->XMLout($outputDirectives);
close($directiveLocationsFile);
# Output the Makefile and XML files containing the directives.
open(my $directivesMakefile,">".$ENV{'BUILDPATH'}."/Makefile_Directives");
foreach my $directive ( keys(%{$includeDirectives}) ) {
    # Construct the unpreprocessed include file name.
    (my $fileName = $includeDirectives->{$directive}->{'fileName'}) =~ s/\.inc$/\.Inc/;
    # Build a list of extra dependencies.
    my @extraDependencies;
    # For "function" actions, add the file containing the directive as an additional dependency
    # as these files are simply copied into the include file as part of the include file construction.
    push(@extraDependencies,keys(%{$nonIncludeDirectives->{$1}->{'files'}}))
	if ( $directive =~ m/^([a-zA-Z0-9_]+)\.function$/ );
    # Add on any other dependencies.
    push(@extraDependencies,keys(%{$nonIncludeDirectives->{$1}->{'dependency'}}))
	if ( $directive =~ m/^([a-zA-Z0-9_]+)\.(moduleUse|functionCall)$/ &&  exists($nonIncludeDirectives->{$1}->{'dependency'}) );
    # Output the Makefile rule. Add an explicit dependence on "hdf5FCInterop.dat" to ensure that HDF% type interoperability is
    # determined before attempting to build the file.
    print $directivesMakefile $fileName.": ".$ENV{'BUILDPATH'}."/".$directive.".xml ".join(" ",@extraDependencies)." \$(BUILDPATH)/hdf5FCInterop.dat\n";
    print $directivesMakefile "\t./scripts/build/buildCode.pl ".$installDirectoryName." ".$ENV{'BUILDPATH'}."/".$directive.".xml\n";
    print $directivesMakefile "\n";
    # Output the directive itself.
    open(my $directiveFile,">".$ENV{'BUILDPATH'}."/".$directive.".xml.tmp");
    print $directiveFile $includeDirectives->{$directive}->{'xml'};
    close($directiveFile);
    &File::Changes::Update($installDirectoryName."/".$ENV{'BUILDPATH'}."/".$directive.".xml",$installDirectoryName."/".$ENV{'BUILDPATH'}."/".$directive.".xml.tmp");
}
# Add additional dependencies for object files of source files that contain functionClass directives. These source files get other
# source files incorporated into them via the source tree preprocessor.
foreach my $directiveName ( keys(%{$functionClasses}) ) {
    print $directivesMakefile $functionClasses->{$directiveName}.": ".join(" ",keys(%{$nonIncludeDirectives->{$directiveName}->{'files'}}))."\n\n";
    # Include explicit dependencies for Makefile_Use_Dependencies to ensure that module dependencies get rebuilt
    # after these directive include files are constructed.
    print $directivesMakefile $ENV{'BUILDPATH'}."/Makefile_Use_Dependencies: ".join(" ",map {(my $fileName = $includeDirectives->{$_}->{'fileName'}) =~ s/\.inc$/\.Inc/; $fileName} keys(%{$includeDirectives}))."\n\n";
}
# Include a rule for including Makefile_Component_Includes. This has to go here since Makefile_Component_Includes depends on
# objects.nodes.components.Inc for which Makefile_Directive contains the build rule.
print $directivesMakefile "-include ".$ENV{'BUILDPATH'}."/Makefile_Component_Includes\n";
print $directivesMakefile $ENV{'BUILDPATH'}."/Makefile_Component_Includes: ".$ENV{'BUILDPATH'}."/objects.nodes.components.Inc\n\n";
close($directivesMakefile);
exit;

sub addImplicitDirectives {
    # Add implicit dependencies required by certain directives.
    my $directive          = shift();
    my $implicitDirectives = shift();
    my $fileName           = shift();
    my $filePathName       = shift();
    my %implicitDirectives =
	(
	 stateful         => [ "galacticusStateRetrieveTask", "galacticusStateStoreTask", "galacticusStateSnapshotTask" ],
	 calculationReset => [ "calculationResetTask"                                                                   ]
	);
    foreach my $implicitDirective ( keys(%implicitDirectives) ) {
	if ( exists($directive->{$implicitDirective}) && $directive->{$implicitDirective} eq "yes" ) {
	    foreach my $implicitDependency ( @{$implicitDirectives{$implicitDirective}} ) {
		$nonIncludeDirectives->{$implicitDependency}->{'files'     }->{$filePathName} = 1;
		$nonIncludeDirectives->{$implicitDependency}->{'dependency'}->{$fileName    } = 1;
	    }			    
	}
    }
}
