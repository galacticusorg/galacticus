#!/usr/bin/env perl
use lib './perl';
use System::Redirect;
use XML::Simple;

# Locate all allocatable arrays in the code base and determine their type and dimensionality.

# Define the source directory
if ( $#ARGV != 0 ) {die "Usage: Find_Allocatable_Arrays.pl <sourcedir>"};
my $sourcedir = $ARGV[0];

# Build a list of source directories.
$sourcedirs[0] = $sourcedir."/source";
$bases[0] = "";
if ( -e $sourcedir."/Source_Codes" ) {
    opendir(sdir,"$sourcedir/Source_Codes");
    while ( $line = readdir sdir ) {
	$line =~ s/\n//;
	if ( -d $sourcedir."/Source_Codes/".$line ) {
	    unless ( $line =~ m/^\.+$/ ) {
		$sourcedirs[++$#sourcedirs] = $sourcedir."/Source_Codes/".$line;
		$bases[++$#bases] = "Source_Codes/".$line."/";
	    }
	}
    }
    closedir(sdir);
}

# Open the source directorys and scan
$ibase=-1;
foreach $srcdir ( @sourcedirs ) {
    ++$ibase;
    $base = $bases[$ibase];
    opendir(indir,$srcdir) or die "Can't open the source directory: #!";
    while (my $fname = readdir indir) {
	
	if ( ( ( ( lc($fname) =~ m/\.f(90)??$/ ) && ! -e $srcdir."/".$fname."t" ) || ( lc($fname) =~ m/\.f90t$/ ) ) && lc($fname) !~ m/^\.\#/ ) {
	    my $fullname = "$srcdir/$fname";
	    $scanfiles[++$#scanfiles] = $fullname;
	    while ( $#scanfiles >= 0 ) {
		$fullname = $scanfiles[$#scanfiles];
		--$#scanfiles;
		if ( ! -e $fullname ) {
		    if ( $fullname =~ m/\/([\w\.]+?)$/ ) {
			$leaf = $1;
		    } else {
			$leaf=$fullname;
		    }
		    system("$make $leaf");
		}
		open(infile,$fullname) or die "Can't open input file: $fullname";
		while (my $line = <infile>) {
		    if ( $line =~ m/,\s*allocatable\s*[,:]/i && $line =~ m/\(\s*:/i && $line =~ m/^\s*([a-zA-Z0-9_\s]+)(\((len|kind)=[\sa-z0-9_]+\))?\s*,/i && $line !~ m/\(\s*kind\s*=\s*HID/ ) {
			while ( $line =~ m/&\s*$/ ) {
			    $line =~ s/&\s*$//;
			    $tline = <infile>;
			    $tline =~ s/^\s*&//;
			    $line .= $tline;
			}
			# Found a line declaring allocatable arrays.
			$line =~ m/^\s*([a-zA-Z0-9_\s]+)/i;
			$vtype = $1;
			$vtype =~ s/\s*$//;
			if ( $line =~ m/\(kind=([a-zA-Z0-9_]+)\)/i ) {
			    $kind = $1;
			} else {
			    $kind = "";
			}
			if ( $line =~ m/dimension\(([:,]+)\)/ ) {
			    $dimstr = $1;
			    $dimension = () = $dimstr =~ /:/g;
			    $type = lc($vtype).":".$dimension.":".$kind;
			    ++$type_dim{$type};
			} else {
			    while ( $line =~ m/\(([:,]+)\)/ ) {
				$dimstr = $1;
				$dimension = () = $dimstr =~ /:/g;
				$line =~ s/\([:,]+\)//;
				$type = lc($vtype).":".$dimension.":".$kind;
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
$iAllocatable = -1;
for $type ( sort keys %type_dim ) {
    @entries = split(/:/,$type);
    ++$iAllocatable;
    ${${$allocatables{'allocatable'}}[$iAllocatable]}{"type"}      = $entries[0];
    ${${$allocatables{'allocatable'}}[$iAllocatable]}{"dimension"} = $entries[1];
    ${${$allocatables{'allocatable'}}[$iAllocatable]}{"kind"}      = $entries[2];
}

# Create XML object.
$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"allocatables");
open(outHndl,">./work/build/Allocatable_Arrays.xml.tmp");
print outHndl $xmlOutput->XMLout(\%allocatables);
close(outHndl);

# Replace original file only if it has changed.
if ( -e $sourcedir."/work/build/Allocatable_Arrays.xml" ) {
    &SystemRedirect::tofile("diff -q  $sourcedir/work/build/Allocatable_Arrays.xml.tmp $sourcedir/work/build/Allocatable_Arrays.xml","/dev/null");
    if ( $? == 0 ) {
	system("rm -f $sourcedir/work/build/Allocatable_Arrays.xml.tmp");
    } else {
	system("mv $sourcedir/work/build/Allocatable_Arrays.xml.tmp $sourcedir/work/build/Allocatable_Arrays.xml");
    }
} else {
    system("mv $sourcedir/work/build/Allocatable_Arrays.xml.tmp $sourcedir/work/build/Allocatable_Arrays.xml");
}

exit;
