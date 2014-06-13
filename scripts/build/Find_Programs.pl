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

# Locate files which contain programs and append to a list of executables.

# Define the source directory
die "Usage: Find_Programs.pl sourcedir"
    unless ( scalar(@ARGV) == 1 );
my $sourcedir = $ARGV[0];

#
# Open an output file
open(my $outfile,">$sourcedir/work/build/Makefile_All_Execs");

#
# Specify work directory.
my $workDir = "/work/build/";

# Build a list of source directories.
my @sourcedirs = ( $sourcedir."/source" );
my @bases      = ( ""                   );

#
# Open the source directorys and scan
#
my @exes;
my $ibase = -1;
foreach my $srcdir ( @sourcedirs ) {
    ++$ibase;
    my $base = $bases[$ibase];
    opendir(my $indir,$srcdir) or die "Can't open the source directory: #!";
    while (my $fname = readdir $indir) {

	if ( $fname =~ m/\.[fF](90)?t?$/ && lc($fname) !~ m/^\.\#/ ) {	    
	    my $pname = $fname;
	    my $fullname = "$srcdir/$fname";
	    my $doesio = 0;
	    open(my $infile,$fullname) or die "Can't open input file: #!";
	    my @fileinfo = stat $infile;
	    
	    my $hasuses = 0;
	    my $oname = $fname;
	    $oname =~ s/\.[fF](90)?t?$/\.o/;
	    my @incfiles;
	    my @modfiles;
	    my $exclude = 0;
	    while (my $line = <$infile>) {
		$exclude = 1
		    if ( $line =~ m/^\s*!\/\s+exclude/ );
		if ( $line =~ m/^\s*program\s/i ) {
		    my $ename = $fname;
		    $ename =~ s/\.[fF](90)?t?$/\.exe/;
		    my $do_entry = 0;
		    if ( $ibase == 0 ) {
			push(@exes,$ename)
			    unless ( $exclude == 1 );
			if ( ! -e $sourcedir."/Extensions/Sources/".$fname ) {$do_entry = 1};
		    } else {
			if ( ! -e $sourcedir."/$fname" && $exclude == 0 ) {push(@exes,$ename)};
			$do_entry = 1;
		    }
		    if ( $do_entry == 1 ) {
			my $ofile = $fname;
			$ofile =~ s/\.[fF](90)?t?$/\.o/;
			my $dfile = $fname;
			$dfile =~ s/\.[fF](90)?t?$/\.d/;
			my $root = $fname;
			$root =~ s/\.[fF](90)?t?$//;
			my $eleaf = $root.".exe";
			print $outfile "$root.exe: .$workDir$base$ofile .$workDir$base$dfile \$(MAKE_DEPS)\n";
			print $outfile "\t\$(FCCOMPILER) `cat .$workDir$base$dfile` -o $root.exe \$(FCFLAGS) `scripts/build/Library_Dependencies.pl $root.exe \$(FCFLAGS)`\n";
			print $outfile "\t./scripts/build/Find_Executable_Size.pl $root.exe .$workDir$root.size\n";
			print $outfile "\t./scripts/build/Find_Parameter_Dependencies.pl `pwd` $root.exe\n\n";
		    }
                }
	    }
	    close($infile);
	}
    }
    closedir($indir);
}

print $outfile "\nall_exes = ".join(" ",@exes)."\n";
close($outfile);

exit;

