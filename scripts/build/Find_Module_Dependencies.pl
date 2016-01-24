#!/usr/bin/env perl
use strict;
use warnings;

# Locate source files which create modules

# Define the source directory
die "Usage: Find_Module_Dependencies.pl sourcedir"
    unless ( scalar(@ARGV) == 1 );
my $sourcedir = $ARGV[0];

#
# Open an output file
open(my $outfile,">".$sourcedir."/".$ENV{'BUILDPATH'}."/Makefile_Module_Deps");

# Specify a work directory.
my $workDir = $ENV{'BUILDPATH'}."/";

#
# Build a list of source directories,
my @sourcedirs = ( $sourcedir."/source" );
my @bases      = ( "" );
if ( -e $sourcedirs[0] ) {
    opendir(my $sdir,$sourcedirs[0]);
    while ( my $line = readdir $sdir ) {
	$line =~ s/\n//;
	if ( -d $sourcedirs[0]."/".$line ) {
	    unless ( $line =~ m/^\.+$/ ) {
		push(@sourcedirs,$sourcedirs[0]."/".$line    );
		push(@bases     ,                   $line."/");
		system("mkdir -p ".$ENV{'BUILDPATH'}."/".$line);
	    }
	}
    }
    closedir($sdir);
}

#
# Open the source directory
#
my $ibase=-1;
foreach my $srcdir ( @sourcedirs ) {
    ++$ibase;
    my $base = $bases[$ibase];
    opendir(my $indir,$srcdir) or die "Can't open the source directory:";
    while (my $fname = readdir $indir) {
	if ( ( ( ( lc($fname) =~ m/\.f(90)??$/ ) && ! -e $srcdir."/".$fname."t" ) || ( lc($fname) =~ m/\.f90t$/ ) ) && lc($fname) !~ m/^\.\#/ ) {
	    my $pname = $fname;
	    my $fullname = "$srcdir/$fname";
	    my $doesio = 0;
	    my $oname = $fname;
	    $oname =~ s/\.f90$/\.o/;
	    $oname =~ s/\.F90$/\.o/;
	    $oname =~ s/\.f90t$/\.o/;
	    $oname =~ s/\.F90t$/\.o/;
	    $oname =~ s/\.f$/\.o/;
	    $oname =~ s/\.F$/\.o/;

	    my @generics;
	    my @allmods;
	    my @scanfiles = ( $fullname );
	    while ( scalar(@scanfiles) > 0 ) {
		$fullname = pop(@scanfiles);
		open(my $infile,$fullname) or die "Can't open input file: $fullname";
		
		while (my $line = <$infile>) {
		    # Locate lines which contain a Forpedo "#definetype"
		    if ( $line =~ m/^#definetype\s+([a-zA-Z0-9_]+)\s+([a-zA-Z0-9_]+)\s/i ) {
			my $label = lc($1);
			my $type = lc($2);
			$generics[$label] .= " ".$type;
		    }

                    # Locate any lines which use the "module" statement and extract the name of the file they use
		    $line =~ s/!.*$//;
		    $line =~ s/\s+$//;
		    if ( lc($line) =~ m/^\s*module\s/ && lc($line) !~ m/^\s*module\s+procedure/ ) {
			# If necessary, generate list of non-generic name replacements to swap out (for Forpedo).
			my @swaplist;
			if ( $line =~ m/<([a-zA-Z0-9_]+)>/i ) {
			    # Generic type module.
			    my $label = lc($1);
			    my $swap = $generics[$label];
			    $swap =~ s/^\s+//;
			    @swaplist = split(/\s+/, $swap);
			} else {
			    # Standard module.
			    @swaplist = ( "null" ); 
			}
			foreach my $swap ( @swaplist ) {
			    my $tline = lc($line);
			    $tline =~ s/<[a-zA-Z0-9_]+>/$swap/;
			    my $startpos = index($tline,"module")+7;
			    my $sublen = length($tline)-$startpos;
			    my $incfile = substr($tline,$startpos,$sublen).".mod";
			    $incfile =~ s/\r//g;
			    push(@allmods,$incfile);
			    print $outfile $workDir.lc($incfile),": ",$workDir,$base,$oname,"\n";
			    print $outfile "\t\@if [ ! -f ",$workDir.lc($incfile)," ]; then \\\n";
			    print $outfile "\t  rm $workDir$base$oname ; \\\n";
			    print $outfile "\t  \$(MAKE) $workDir$base$oname ; \\\n";
			    print $outfile "\tfi\n\n";
			    my $dname = $oname;
			    $dname =~ s/.o$/.d/;
			    print $outfile $workDir,lc($incfile),".d: ",$workDir,$base,$dname,"\n";
			    print $outfile "\t\@echo ",$workDir.$base,$oname," > ",$workDir,lc($incfile),".d\n";
			    print $outfile "\t\@cat ",$workDir.$base,$dname," >> ",$workDir,lc($incfile),".d\n";
			    
			    # Create rule for making a *.mod.gv file which is used in building GraphViz descriptions of source
			    # file dependencies.
			    my $gvname = $pname.".gv";
			    my $modname = $incfile;
			    $modname =~ s/.mod$//;
			    print $outfile $workDir,lc($incfile),".gv: ",$workDir,$base,$dname," ",$workDir,$base,$gvname."\n";
			    print $outfile "\t\@echo ",$base,$pname." > ",$workDir,lc($incfile).".gv\n";

			}
		    }
		    if ( $line =~ m/include\s+\'(\w+)\'/i ) {push(@scanfiles,$sourcedir."/".$1);};
		}
		close($infile);

		# Create a rule for the module list file if needed.
		if ( @allmods ) {
		    my $mname = $fname;
		    $mname =~ s/\.[fFt90]+$/\.m/;
		    print $outfile "$workDir$base$mname:\n";
		    my $first = 1;
		    foreach my $mod ( @allmods ) {
			if ( $first == 1 ) {
			    $first = 0;
			    print $outfile "\t\@echo $workDir$mod > $workDir$base$mname\n";
			} else {
			    print $outfile "\t\@echo $workDir$mod >> $workDir$base$mname\n";
			}
		    }
		    print $outfile "\n";
		}

	    }
	}
    }
    closedir($indir);
}
