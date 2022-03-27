#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Data::Dumper;

# Build a Makefile with dependencies for libgalacticus.
# Andrew Benson (25-March-2022)

# Get an XML parser.
my $xml = new XML::Simple;

# Get the list of functionClasses that should be compiled into the library.
my $libraryFunctionClasses = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/source/libraryClasses.xml");

# Generate the Makefile rules.
my $rules;
$code::buildPath = $ENV{'BUILDPATH'}.($ENV{'BUILDPATH'} =~ m/\/$/ ? "" : "/");
foreach my $functionClass ( &List::ExtraUtils::hashList($libraryFunctionClasses->{'classes'}, keyAs => "name") ) {
    $code::name  = $functionClass->{'name'};
    $rules      .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$buildPath}libgalacticus/{$name}.F90:
	./scripts/build/libraryInterfaces.pl
{$buildPath}libgalacticus/{$name}.p.F90.up: {$buildPath}libgalacticus/{$name}.F90 {$buildPath}hdf5FCInterop.dat {$buildPath}openMPCriticalSections.xml
	./scripts/build/preprocess.pl {$buildPath}libgalacticus/{$name}.F90 {$buildPath}libgalacticus/{$name}.p.F90
{$buildPath}libgalacticus/{$name}.p.F90 : | {$buildPath}libgalacticus/{$name}.p.F90.up
	@true
{$buildPath}libgalacticus/{$name}.o : {$buildPath}libgalacticus/{$name}.p.F90 {$buildPath}libgalacticus/{$name}.d Makefile
	@mkdir -p {$buildPath}/moduleBuild
	$(FCCOMPILER) -c {$buildPath}libgalacticus/{$name}.p.F90 -o {$buildPath}libgalacticus/{$name}.o $(FCFLAGS) 2>&1 | ./scripts/build/postprocess.pl {$buildPath}libgalacticus/{$name}.p.F90

CODE
}

# Add rules for building module dependency files.
$code::dependencies = join(" ",map {$code::buildPath."libgalacticus/".$_.".d"} sort(keys(%{$libraryFunctionClasses->{'classes'}})));
$rules .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$buildPath}libgalacticus_classes.d: {$dependencies}
	@rm -f {$buildPath}libgalacticus_classes.d~
CODE
foreach my $functionClass ( &List::ExtraUtils::hashList($libraryFunctionClasses->{'classes'}, keyAs => "name") ) {
    $code::name  = $functionClass->{'name'};
    $rules      .= fill_in_string(<<'CODE', PACKAGE => 'code');
	@cat {$buildPath}libgalacticus/{$name}.d >> {$buildPath}libgalacticus_classes.d~
CODE
}
$rules .= fill_in_string(<<'CODE', PACKAGE => 'code');
	@sort -u {$buildPath}libgalacticus_classes.d~ -o {$buildPath}libgalacticus_classes.d~
	@if cmp -s {$buildPath}libgalacticus_classes.d {$buildPath}libgalacticus_classes.d~ ; then \
	 rm {$buildPath}libgalacticus_classes.d~ ; \
	else \
	 mv {$buildPath}libgalacticus_classes.d~ {$buildPath}libgalacticus_classes.d ; \
	fi
CODE

# Add object file dependencies.
$code::dependencies = join(" ",map {$code::buildPath."libgalacticus/".$_.".o"} sort(keys(%{$libraryFunctionClasses->{'classes'}})));
$rules .= fill_in_string(<<'CODE', PACKAGE => 'code');
libgalacticus.so: {$dependencies}
CODE

# Create the Makefile.
open(my $makefile,">",$code::buildPath."Makefile_Library_Dependencies");
print $makefile $rules;
close($makefile);

exit;
