#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use XML::Simple;
use System::Redirect;

# Locate all allocatable arrays in the code base and determine their type and dimensionality.

# Define the source directory
die "Usage: Find_Allocatable_Arrays.pl <sourcedir>"
    unless ( scalar(@ARGV) == 1 );
my $sourcedir = $ARGV[0];

# Build a list of source directories.
my @sourcedirs = ( $sourcedir."/source" );
my @bases      = ( "" );

# Open the source directorys and scan
my $ibase = -1;
my %type_dim;
foreach my $srcdir ( @sourcedirs ) {
    ++$ibase;
    my $base = $bases[$ibase];
    opendir(my $indir,$srcdir) or die "Can't open the source directory: #!";
    while (my $fname = readdir $indir) {
	
	if ( ( ( ( lc($fname) =~ m/\.f(90)??$/ ) && ! -e $srcdir."/".$fname."t" ) || ( lc($fname) =~ m/\.f90t$/ ) ) && lc($fname) !~ m/^\.\#/ ) {
	    my $fullname = "$srcdir/$fname";
	    my @scanfiles = ( $fullname );
	    while ( scalar(@scanfiles) > 0 ) {
		my $fullname = pop(@scanfiles);
		if ( ! -e $fullname ) {
		    my $leaf;
		    if ( $fullname =~ m/\/([\w\.]+?)$/ ) {
			$leaf = $1;
		    } else {
			$leaf=$fullname;
		    }
		    system("make $leaf");
		}
		open(my $infile,$fullname) or die "Can't open input file: $fullname";
		while (my $line = <$infile>) {
		    if ( $line =~ m/,\s*allocatable\s*[,:]/i && $line =~ m/\(\s*:/i && $line =~ m/^\s*([a-zA-Z0-9_\s]+)(\((len|kind)=[\sa-z0-9_]+\))?\s*,/i && $line !~ m/\(\s*kind\s*=\s*HID/ && $line !~ m/\(\s*kind\s*=\s*c_size_t/ && $line !~ m/\(\s*kind\s*=\s*c_char/ ) {
			while ( $line =~ m/&\s*$/ ) {
			    $line =~ s/&\s*$//;
			    my $tline = <$infile>;
			    $tline =~ s/^\s*&//;
			    $line .= $tline;
			}
			# Found a line declaring allocatable arrays.
			$line =~ m/^\s*([a-zA-Z0-9_\s]+)/i;
			my $vtype = $1;
			$vtype =~ s/\s*$//;
			my $kind;
			if ( $line =~ m/\(kind=([a-zA-Z0-9_]+)\)/i ) {
			    $kind = $1;
			} else {
			    $kind = "";
			}
			if ( $line =~ m/dimension\(([:,]+)\)/ ) {
			    my $dimstr = $1;
			    my $dimension = () = $dimstr =~ /:/g;
			    my $type = lc($vtype).":".$dimension.":".$kind;
			    ++$type_dim{$type};
			} else {
			    while ( $line =~ m/\(([:,]+)\)/ ) {
				my $dimstr = $1;
				my $dimension = () = $dimstr =~ /:/g;
				$line =~ s/\([:,]+\)//;
				my $type = lc($vtype).":".$dimension.":".$kind;
				++$type_dim{$type};
			    }
			}
		    }
		}
	    }
	}
    }
}

# Create auxilliary file containing the list of dimensions and lengths for character arrays.
my $iAllocatable = -1;
my %allocatables;
foreach my $type ( sort keys %type_dim ) {
    my @entries = split(/:/,$type);
    ++$iAllocatable;
    ${${$allocatables{'allocatable'}}[$iAllocatable]}{"type"}      =                        $entries[0]     ;
    ${${$allocatables{'allocatable'}}[$iAllocatable]}{"dimension"} =                        $entries[1]     ;
    ${${$allocatables{'allocatable'}}[$iAllocatable]}{"kind"}      = defined($entries[2]) ? $entries[2] : "";
}

# Create XML object.
my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"allocatables");
open(outHndl,">".$ENV{'BUILDPATH'}."/Allocatable_Arrays.xml.tmp");
print outHndl $xmlOutput->XMLout(\%allocatables);
close(outHndl);

# Replace original file only if it has changed.
if ( -e $sourcedir."/".$ENV{'BUILDPATH'}."/Allocatable_Arrays.xml" ) {
    &System::Redirect::tofile("diff -q  ".$sourcedir."/".$ENV{'BUILDPATH'}."/Allocatable_Arrays.xml.tmp ".$sourcedir."/".$ENV{'BUILDPATH'}."/Allocatable_Arrays.xml","/dev/null");
    if ( $? == 0 ) {
	system("rm -f ".$sourcedir."/".$ENV{'BUILDPATH'}."/Allocatable_Arrays.xml.tmp");
    } else {
	system("mv ".$sourcedir."/".$ENV{'BUILDPATH'}."/Allocatable_Arrays.xml.tmp ".$sourcedir."/".$ENV{'BUILDPATH'}."/Allocatable_Arrays.xml");
    }
} else {
    system("mv ".$sourcedir."/".$ENV{'BUILDPATH'}."/Allocatable_Arrays.xml.tmp ".$sourcedir."/".$ENV{'BUILDPATH'}."/Allocatable_Arrays.xml");
}

exit;
