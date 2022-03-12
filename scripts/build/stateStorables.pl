#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Galacticus::Build::Directives;
use Galacticus::Build::SourceTree;
use List::ExtraUtils;
use Data::Dumper;
use Storable;
use Fortran::Utils;
use File::Changes;

# Builds a list of classes which support store/restore of their state.
# Andrew Benson (10-April-2018)

# Read command line arguments.
die "Usage: stateStoreables.pl <installDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $installDirectoryName = $ARGV[0];
# Get an XML parser.
my $xml                  = new XML::Simple();
# Parse the directive locations file.
my $directiveLocations   = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
# Initialize data structure to hold per-file information.
my $storablesPerFile;
my $havePerFile = -e $ENV{'BUILDPATH'}."/stateStorables.blob";
my $updateTime;
if ( $havePerFile ) {
    $storablesPerFile = retrieve($ENV{'BUILDPATH'}."/stateStorables.blob");
    $updateTime       = -M       $ENV{'BUILDPATH'}."/stateStorables.blob" ;
}
# Find all files which contain functionClass objects - these all support state store/restore.
foreach my $functionClassFileName ( &List::ExtraUtils::as_array($directiveLocations->{'functionClass'}->{'file'}) ) {
    (my $fileIdentifier = $functionClassFileName) =~ s/\//_/g;
    $fileIdentifier =~ s/^\._??//;
    # Extract a functionClass directives from this file.
    unless ( $havePerFile && exists($storablesPerFile->{$fileIdentifier}) && -M $functionClassFileName > $updateTime  ) {
	delete($storablesPerFile->{$fileIdentifier});
	# Look for modules in this file.
	my @modules;
	open(my $code,$functionClassFileName);
	do {
	    my ($rawLine, $processedLine, $bufferedComments);
	    &Fortran::Utils::Get_Fortran_Line($code,$rawLine,$processedLine,$bufferedComments);
	    if ( my @matches = ( $processedLine =~ $Fortran::Utils::unitOpeners{'module'}->{'regEx'} ) ) {
		push(@modules,$matches[$Fortran::Utils::unitOpeners{'module'}->{'unitName'}]);
	    }
	} until ( eof($code) );
	close($code);
	foreach my $functionClass ( &Galacticus::Build::Directives::Extract_Directives($functionClassFileName,'functionClass') ) {
	    my $functionClassDescriptor;
	    $functionClassDescriptor->{'name'  } = $functionClass->{'name'}."Class";
	    $functionClassDescriptor->{'module'} = $modules[0]
		if ( scalar(@modules) == 1 );
	    push(@{$storablesPerFile->{$fileIdentifier}->{'functionClasses'        }},$functionClassDescriptor);
	    push(@{$storablesPerFile->{$fileIdentifier}->{'functionClassDirectives'}},$functionClass->{'name'});
	}
    }
    # Iterate over functionClasses in this file.
    foreach my $functionClassName ( @{$storablesPerFile->{$fileIdentifier}->{'functionClassDirectives'}} ) {
	# Find all files which contain instances of this functionClass.
	foreach my $instanceFileName ( &List::ExtraUtils::as_array($directiveLocations->{$functionClassName}->{'file'}) ) {
	    (my $instanceFileIdentifier = $instanceFileName) =~ s/\//_/g;
	    $instanceFileIdentifier =~ s/^\._??//;
	    next
		if ( $havePerFile && exists($storablesPerFile->{$instanceFileIdentifier}) && -M $instanceFileName > $updateTime  );
	    delete($storablesPerFile->{$instanceFileIdentifier});
	    # Extract the relevant directives from this file.
	    foreach my $instance ( &Galacticus::Build::Directives::Extract_Directives($instanceFileName,$functionClassName) ) {
		push(@{$storablesPerFile->{$instanceFileIdentifier}->{'functionClassInstances'}},$instance->{'name'});
	    }
	}
    }   
}
# Find all files which contain functionClassType objects - these all support state store/restore.
foreach my $functionClassFileName ( &List::ExtraUtils::as_array($directiveLocations->{'functionClassType'}->{'file'}) ) {
    (my $fileIdentifier = $functionClassFileName) =~ s/\//_/g;
    $fileIdentifier =~ s/^\._??//;
    next
	if ( $havePerFile && exists($storablesPerFile->{$fileIdentifier}) && -M $functionClassFileName > $updateTime  );
    delete($storablesPerFile->{$fileIdentifier});
    # Extract a functionClassType directives from this file.
    push(@{$storablesPerFile->{$fileIdentifier}->{'functionClassTypes'}},map {{name => $_->{'name'}, file => $functionClassFileName}} &Galacticus::Build::Directives::Extract_Directives($functionClassFileName,'functionClassType'));
}
# Find all files which contain stateStorable objects - these explicitly support state store/restore.
foreach my $stateStorableFileName ( &List::ExtraUtils::as_array($directiveLocations->{'stateStorable'}->{'file'}) ) {
    (my $fileIdentifier = $stateStorableFileName) =~ s/\//_/g;
    $fileIdentifier =~ s/^\._??//;
    next
	if ( $havePerFile && exists($storablesPerFile->{$fileIdentifier}) && -M $stateStorableFileName > $updateTime  );
    delete($storablesPerFile->{$fileIdentifier});
    # Parse the source of this file.
    my $tree = &Galacticus::Build::SourceTree::ParseFile($stateStorableFileName);
    # Walk the tree.
    my $node           = $tree;
    my $depth          = 0    ;
    my @directiveNodes        ;
    my %classes               ;
    while ( $node ) {
	# Capture stateStorable directives.
	push(@directiveNodes,$node)
	    if ( $node->{'type'} eq "stateStorable" );
	# Capture derived type definitions.
	if ( $node->{'type'} eq "type" ) {
	    # Parse class openers to find dependencies.
	    if ( $node->{'opener'} =~ m/^\s*type\s*(,\s*(abstract)\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*([a-z0-9_]+)\s*$/i ) {
		my $type     = $5;
		my $extends  = $3;
		my $abstract = defined($2);
		$classes{$type} =
		{
		     node     => $node   ,
		     extends  => $extends,
		     abstract => $abstract
		};
	    } else {
		die("stateStorables.pl: unable to parse type definition opener");
	    }
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
    # Process any stateStorable directives.
    foreach my $directiveNode ( @directiveNodes ) {
	my $directive = $directiveNode->{'directive'};
	my $moduleName = "";
	my $parentNode = $directiveNode;
	while ( $parentNode ) {
	    if ( $parentNode->{'type'} eq "module" && $parentNode->{'opener'} =~ m/^module\s+([a-zA-Z0-9_]+)/ ) {
		$moduleName = $1;
		last;
	    }
	    $parentNode = $parentNode->{'parent'};
	}
	# Scan all known classes, finding all which derive from the base class.
	foreach my $className ( sort(keys(%classes)) ) {
	    my $matches         = 0;
	    my $parentClassName = $className;
	    while ( defined($parentClassName) ) {
		if ( $parentClassName eq $directive->{'class'} ) {
		    $matches = 1;
		    last;
		}
		$parentClassName = $classes{$parentClassName}->{'extends'};
	    }
	    push(@{$storablesPerFile->{$fileIdentifier}->{'stateStorables'}},{type => $className, class => $directive->{'class'}, module => $moduleName})
		if ( $matches );
	}
    }
}
# Sort results.
my $stateStorables;
@{$stateStorables->{'functionClasses'       }} = sort {$a->{'name'} cmp $b->{'name'}} map {exists($storablesPerFile->{$_}->{'functionClasses'       }) ? @{$storablesPerFile->{$_}->{'functionClasses'       }} : ()} keys(%{$storablesPerFile});
@{$stateStorables->{'functionClassTypes'    }} = sort                                 map {exists($storablesPerFile->{$_}->{'functionClassTypes'    }) ? @{$storablesPerFile->{$_}->{'functionClassTypes'    }} : ()} keys(%{$storablesPerFile});
@{$stateStorables->{'functionClassInstances'}} = sort                                 map {exists($storablesPerFile->{$_}->{'functionClassInstances'}) ? @{$storablesPerFile->{$_}->{'functionClassInstances'}} : ()} keys(%{$storablesPerFile});
@{$stateStorables->{'stateStorables'        }} = sort {$a->{'type'} cmp $b->{'type'}} map {exists($storablesPerFile->{$_}->{'stateStorables'        }) ? @{$storablesPerFile->{$_}->{'stateStorables'        }} : ()} keys(%{$storablesPerFile});
# Output the results.
open(my $outputFile,">".$ENV{'BUILDPATH'}."/stateStorables.xml.tmp");
print $outputFile $xml->XMLout($stateStorables, RootName => "storables");
close($outputFile);
&File::Changes::Update($ENV{'BUILDPATH'}."/stateStorables.xml",$ENV{'BUILDPATH'}."/stateStorables.xml.tmp");
# Output the per file module use data.
store($storablesPerFile,$ENV{'BUILDPATH'}."/stateStorables.blob.tmp");
&File::Changes::Update($ENV{'BUILDPATH'}."/stateStorables.blob",$ENV{'BUILDPATH'}."/stateStorables.blob.tmp");

exit 0;
