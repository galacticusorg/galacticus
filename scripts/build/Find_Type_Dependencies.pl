#!/usr/bin/env perl
use lib './perl';
use System::Redirect;
use XML::Simple;

# Find which modules define all derived type variables.

# Define the source directory
if ( $#ARGV != 0 ) {die "Usage: Find_Type_Dependencies.pl <sourcedir>"};
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
		    if ( $line !~ m/^\s*module\s+procedure\s/i && $line =~ m/^\s*module\s+([a-zA-Z0-9_]+)/i ) {$moduleName = $1};
		    if ( $line =~ m/^\s*type\s+([a-zA-Z0-9_]+)\s*(!.*)??$/i ) {
			$type = lc($1);
			$dependencies{$type} = $moduleName;
		    }
		}
	    }
	}
    }
}

# Create XML object.
$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"dependencies");
open(outHndl,">./work/build/Type_Definitions.xml.tmp");
print outHndl $xmlOutput->XMLout(\%dependencies);
close(outHndl);

# Replace original file only if it has changed.
if ( -e $sourcedir."/work/build/Type_Definitions.xml" ) {
    &SystemRedirect::tofile("diff -q  $sourcedir/work/build/Type_Definitions.xml.tmp $sourcedir/work/build/Type_Definitions.xml","/dev/null");
    if ( $? == 0 ) {
	system("rm -f $sourcedir/work/build/Type_Definitions.xml.tmp");
    } else {
	system("mv $sourcedir/work/build/Type_Definitions.xml.tmp $sourcedir/work/build/Type_Definitions.xml");
    }
} else {
    system("mv $sourcedir/work/build/Type_Definitions.xml.tmp $sourcedir/work/build/Type_Definitions.xml");
}

exit;
