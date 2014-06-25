#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use Fcntl qw(SEEK_SET);
use XML::Simple;
use Data::Dumper;
require System::Redirect;

# Scans source code for "!#" directives and generates a Makefile.

# Define the source directory.
die "Usage: Code_Directive_Parser.pl sourcedir"
    unless ( scalar(@ARGV) == 1 );
my $sourcedir = $ARGV[0];
my @sourcedirs = ( $sourcedir."/source" );

# Specify verbosity.
my $verbosity = 0;

# Create XML object.
my $xml = new XML::Simple;

# Initialize hashes.
my %includeDirectives;
my $otherDirectives;

# Open the source directory.
foreach my $srcdir ( @sourcedirs ) {
    opendir(my $indir,$srcdir) or die "Can't open the source directory: #!";
    while (my $fname = readdir $indir) {	
	if (
	    ( $fname =~ m/\.[fF](90)??t??$/ || $fname =~ m/\.c(pp)??$/ || $fname =~ m/\.h$/ )
	    && $fname !~ m/^\.\#/
	    ) {
	    my$fullname = $srcdir."/".$fname;

	    # Add the file to the list of filenames to process.
	    my @fileNames;
	    my @filePositions;
	    unshift(@fileNames,$fullname);
	    unshift(@filePositions,-1);
	    
	    # Process files until none remain.
	    while ( scalar(@fileNames) > 0 ) {
		
		# Open the file.
		open(my $fileHandle,$fileNames[0]) or die "Can't open input file: #!";
		seek($fileHandle,$filePositions[0],SEEK_SET) unless ( $filePositions[0] == -1 );

		while (my $line = <$fileHandle>) {
		    my $lineNumber = $.;

		    # Detect include files, and recurse into them.
		    if ( $line =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ ) {
			my $includeFile = $srcdir."/".$1;
			$includeFile =~ s/\.inc$/.Inc/;
			if ( -e $includeFile ) {
			    $filePositions[0] = tell($fileHandle);
			    unshift(@fileNames,$includeFile);
			    unshift(@filePositions,-1);
			    last;
			}
		    }

		    if (
			$line =~ m/^\s*!\#\s+(<\s*([a-zA-Z]+)+.*>)\s*$/
			|| $line =~ m/^\s*\/\/\#\s+(<\s*([a-zA-Z]+)+.*>)\s*$/
			) {
			my $xmlCode = $1."\n";
			my $xmlTag  = $2;
			# Read ahead until a matching close tag is found.
			unless ( $xmlCode =~  m/\/>/ ) {
			    my $nextLine = "";
			    until ( $nextLine =~ m/<\/$xmlTag>/ || eof($fileHandle) ) {
				$nextLine = <$fileHandle>;
				$nextLine =~ s/^\s*!\#\s+//;
				$nextLine =~ s/^\s*\/\/\#\s+//;
				$xmlCode .= $nextLine;
			    }
			}
			my $data = eval{$xml->XMLin($xmlCode)};
			die("Code_Directive_Parser.pl failed in ".$fullname." at line ".$lineNumber." with message:\n".$@) if ($@);
			if ( $verbosity == 1 ) {
			    print "$fname : $xmlCode\n";
			    print Dumper($data);
			}
			# Act on the directive.
			my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>$xmlTag);
			if ( $xmlTag eq "include" ) {
			    my $fileName;
			    if ( ${$data}{'content'} =~ m/^\s*\#??include\s*["'<](.+)["'>]/i ) {
				($fileName = "work/build/".$1) =~ s/\.inc$/\.Inc/;
				${$data}{'fileName'} = $fileName;
			    }
			    delete(${$data}{'content'});
			    my $directive = ${$data}{'directive'};
			    $directive = ${$data}{'name'}
			    if ( exists(${$data}{'name'}) );
			    $directive .= ".".${$data}{'type'};
			    ${$includeDirectives{$directive}}{'source'  } = $fileNames[0];
			    ${$includeDirectives{$directive}}{'fileName'} = $fileName;
			    ${$includeDirectives{$directive}}{'xml'     } = $xmlOutput->XMLout($data);
			} else {
			    $otherDirectives->{$xmlTag}->{$srcdir."/".$fname} = 1;
			}
		    }
		}
		if ( eof($fileHandle) ) {
		    shift(@fileNames);
		    shift(@filePositions);
		}
		close($fileHandle);
	    }
	}

    }
    closedir($indir);
}

# Output a file listing which files contain which directives.
my $outputDirectives;
foreach my $xmlTag ( keys(%{$otherDirectives}) ) {
    my @fileNames;
    foreach my $fileName ( keys(%{$otherDirectives->{$xmlTag}} ) ){
	push(@fileNames,$fileName);
    }
    @{$outputDirectives->{$xmlTag}->{'file'}} = @fileNames;
}
my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"directives");
open(directiveHndl,">./work/build/Code_Directive_Locations.xml");
print directiveHndl $xmlOutput->XMLout($outputDirectives);
close(directiveHndl);

# Output the Makefile
open(makefileHndl,">./work/build/Makefile_Directives");
foreach my $directive ( keys(%includeDirectives) ) {
    (my $fileName = ${$includeDirectives{$directive}}{'fileName'}) =~ s/\.inc$/\.Inc/;
    my $extraDependencies = "";
    # For "function" actions, add the file containing the directive as an additional dependency
    # as these files are simply copied into the include file as part of the include file construction.
    if ( $directive =~ m/^([a-zA-Z0-9_]+)\.function$/ ) {
	my $name = $1;
	my @fileNames = keys(%{$otherDirectives->{$name}});
	$extraDependencies .= " ".join(" ",@fileNames);
    }
    print makefileHndl $fileName.": ./work/build/".$directive.".xml".$extraDependencies."\n";
    print makefileHndl "\t./scripts/build/Build_Include_File.pl ".$sourcedir." ./work/build/".$directive.".xml\n";
    print makefileHndl "\n";
    open(xmlHndl,">./work/build/".$directive.".xml.tmp");
    print xmlHndl ${$includeDirectives{$directive}}{'xml'};
    close(xmlHndl);
    &SystemRedirect::tofile("diff -q  $sourcedir/work/build/".$directive.".xml.tmp $sourcedir/work/build/".$directive.".xml","/dev/null");
    if ( $? == 0 ) {
	system("rm -f $sourcedir/work/build/".$directive.".xml.tmp");
    } else {
	system("mv $sourcedir/work/build/".$directive.".xml.tmp $sourcedir/work/build/".$directive.".xml");
    }
}
# Include a rule for including Makefile_Component_Includes. This has to go here since Makefile_Component_Includes depends on
# objects.nodes.components.Inc for which Makefile_Directive contains the build rule.
print makefileHndl "-include ./work/build/Makefile_Component_Includes\n";
print makefileHndl "./work/build/Makefile_Component_Includes: ./work/build/objects.nodes.components.Inc\n\n";
close(makefileHndl);
exit;
