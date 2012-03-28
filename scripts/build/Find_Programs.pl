#!/usr/bin/env perl

# Locate files which contain programs and append to a list of executables.

# Define the source directory
if ( $#ARGV != 0 ) {die "Usage: Find_Programs.pl sourcedir"};
my $sourcedir = $ARGV[0];

#
# Open an output file
open(outfile,">$sourcedir/work/build/Makefile_All_Execs");

#
# Specify work directory.
$workDir = "/work/build/";

# Build a list of source directories.
$sourcedirs[0] = $sourcedir."/source";
$bases[0] = "";

#
# Open the source directorys and scan
#
$ibase=-1;
foreach $srcdir ( @sourcedirs ) {
    ++$ibase;
    $base = $bases[$ibase];
    opendir(indir,$srcdir) or die "Can't open the source directory: #!";
    while (my $fname = readdir indir) {

	if ( $fname =~ m/\.[fF](90)?t?$/ && lc($fname) !~ m/^\.\#/ ) {	    
	    my $pname = $fname;
	    my $fullname = "$srcdir/$fname";
	    my $doesio = 0;
	    open(infile,$fullname) or die "Can't open input file: #!";
	    @fileinfo=stat infile;
	    
	    my $hasuses = 0;
	    my $oname = $fname;
	    $oname =~ s/\.[fF](90)?t?$/\.o/;
	    @incfiles = ();
	    @modfiles = ();
	    while (my $line = <infile>) {
		if ( $line =~ m/^\s*program\s/i ) {
		    $ename = $fname;
		    $ename =~ s/\.[fF](90)?t?$/\.exe/;
		    $do_entry = 0;
		    if ( $ibase == 0 ) {
			@exes[++$#exes] = $ename;
			if ( ! -e $sourcedir."/Extensions/Sources/".$fname ) {$do_entry = 1};
		    } else {
			if ( ! -e $sourcedir."/$fname" ) {@exes[++$#exes] = $ename};
			$do_entry = 1;
		    }
		    if ( $do_entry == 1 ) {
			$ofile = $fname;
			$ofile =~ s/\.[fF](90)?t?$/\.o/;
			$dfile = $fname;
			$dfile =~ s/\.[fF](90)?t?$/\.d/;
			$root = $fname;
			$root =~ s/\.[fF](90)?t?$//;
			$eleaf = $root.".exe";
			print outfile "$root.exe: .$workDir$base$ofile .$workDir$base$dfile \$(MAKE_DEPS)\n";
			print outfile "\t\$(F03COMPILER) `cat .$workDir$base$dfile` -o $root.exe \$(F03FLAGS) \$(LIBS)\n";
			print outfile "\t./scripts/build/Find_Executable_Size.pl $root.exe .$workDir$root.size\n";
			print outfile "\t./scripts/build/Find_Parameter_Dependencies.pl $root.exe\n\n";
		    }
                }
	    }
	    close(infile);
	}
    }
    closedir(indir);
}

print outfile "\nall_exes = @exes\n";
close(outfile);

exit;

