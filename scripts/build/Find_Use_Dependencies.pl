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
unshift(@INC,$galacticusPath."perl"); 
use List::Uniq ':all';
use Data::Dumper;

# Locate source files which have dependencies on modules.

# Define the source directory
die "Usage: Find_Use_Dependencies.pl sourcedir"
    unless ( scalar(@ARGV) == 1 || scalar(@ARGV) == 2 );
my $sourcedir = $ARGV[0];
my $make = "make";
$make = $ARGV[1]
    if ( scalar(@ARGV) == 2 );

# Specify work directory.
my $workDir = "/work/build/";

# List of modules to ignore (as they're external to the source code).
my %ignoreList = (
    "omp_lib" => 1,
    "hdf5" => 1,
    "h5tb" => 1,
    "h5lt" => 1,
    "h5global" => 1,
    "h5fortran_types" => 1,
    "fox_common" => 1,
    "fox_dom" => 1,
    "fox_wxml" => 1,
    "fox_utils" => 1,
    "fgsl" => 1
    );

# Modules that require a library to be linked.
my %moduleLibararies = (
    fftw3      => "fftw3",
    fgsl       => "fgsl_gfortran",
    fox_common => "FoX_common",
    fox_dom    => "FoX_dom",
    fox_wxml   => "FoX_wxml",
    fox_utils  => "FoX_utils",
    hdf5       => "hdf5_fortran"
    );
my %includeLibararies = (
    crypt      => "crypt"
    );

# Open the compiler options file and find preprocessor flags.
my @preprocs;
foreach my $makefile ( "Makefile" ) {
    open(my $ophndl,$makefile);
    while ( my $line = <$ophndl> ) {
        if ( $line =~ m/^\s*FCFLAGS\s*\+??=/ ) {
            while ( $line =~ s/\s\-D([0-9A-Z]+)\s// ) {
                push(@preprocs,$1);
            }
        }
    }
    close($ophndl);
}
my $environmentOptions = $ENV{"GALACTICUS_FLAGS"};
if ( defined($environmentOptions) ) {
    while ( $environmentOptions =~ s/\s\-D([0-9A-Z]+)\s// ) {
	push(@preprocs,$1);
    }
    while ( $environmentOptions =~ s/\s\-D([0-9A-Z]+)$// ) {
	push(@preprocs,$1);
    }
}

# Open an output file
open(my $outfile,">$sourcedir/work/build/Makefile_Use_Deps");

#
# Build a list of source directories.
my @sourcedirs  = ( $sourcedir."/source" );
my @bases       = ( "" );
if ( -e $sourcedirs[0] ) {
    opendir(my $sdir,$sourcedirs[0]);
    while ( my $line = readdir $sdir ) {
	$line =~ s/\n//;
	if ( -d $sourcedirs[0]."/".$line ) {
	    unless ( $line =~ m/^\.+$/ ) {
		push(@sourcedirs,$sourcedirs[0]."/".$line    );
		push(@bases     ,                   $line."/");
	    }
	}
    }
    closedir($sdir);
}

#
# Open the source directorys and scan
#
my @scanfiles;
my $ibase = -1;
foreach my $srcdir ( @sourcedirs ) {
    ++$ibase;
    my $base = $bases[$ibase];
    opendir(my $indir,$srcdir) or die "Can't open the source directory: #!";
    while (my $fname = readdir $indir) {

	if (
	    (
	     ( ( lc($fname) =~ m/\.f(90)??$/ ) && ! -e $srcdir."/".$fname."t" )
	     || ( lc($fname) =~ m/\.f90t$/ )
	     || ( lc($fname) =~ m/\.c(pp)??$/ )
	     || ( lc($fname) =~ m/\.inc$/ )
	    )
	    && lc($fname) !~ m/^\.\#/ 
	    ) {	    
	    my $pname = $fname;
	    my $fullname = "$srcdir/$fname";
	    my $doesio = 0;
	    my $hasuses = 0;
	    my $oname = $fname;
	    $oname =~ s/\.f90$/\.o/;
	    $oname =~ s/\.F90$/\.o/;
	    $oname =~ s/\.f90t$/\.o/;
	    $oname =~ s/\.F90t$/\.o/;
	    $oname =~ s/\.f$/\.o/;
	    $oname =~ s/\.F$/\.o/;
	    $oname =~ s/\.c(pp)??$/\.o/;
	    my @incfiles;
	    my @modfiles;
	    my @extra_includes;
	    my %libraryDependencies;

	    push(@scanfiles,$fullname);
	    while ( scalar(@scanfiles) > 0 ) {
		$fullname = pop(@scanfiles);
		if ( ! -e $fullname ) {
		    my $leaf;
		    if ( $fullname =~ m/\/([\w\.]+?)$/ ) {
			$leaf = $1;
		    } else {
			$leaf=$fullname;
		    }
		    system("make $leaf");
		}
		my @preproc_stack_name;
		my @preproc_stack_state;
		$libraryDependencies{"stdc++"} = 1
		    if ( lc($fname) =~ m/\.cpp$/ );
		my $active = 1;
		open(my $infile,$fullname) or die "Can't open input file: $fullname";
		while (my $line = <$infile>) {
		    if ( $line =~ m/^\s*\!;\s*([a-zA-Z0-9_]+)\s*$/ ) {
			$libraryDependencies{$1} = 1;	
		    }
		    if ( $line =~ m/^\s*\#include\s+<([a-zA-Z0-9_]+)\.h>/ ) {
			my $includeFile = $1;
			$libraryDependencies{$includeLibararies{lc($includeFile)}} = 1
			    if ( exists($includeLibararies{lc($includeFile)}) );
		    }
		    if ( $line =~ m/^\#/ ) {
			if ( $line =~ m/^\#ifdef\s+([0-9A-Z_]+)\s*$/ ) {
			    my $this_preproc = $1;
			    push(@preproc_stack_name ,$this_preproc);
			    push(@preproc_stack_state,1            );
			}
			if ( $line =~ m/^\#ifndef\s+([0-9A-Z_]+)\s*$/ ) {
			    my $this_preproc = $1;
			    push(@preproc_stack_name ,$this_preproc);
			    push(@preproc_stack_state,0            );
			}
			if ( $line =~ m/^\#endif\s*$/ ) {
			    pop(@preproc_stack_name );
			    pop(@preproc_stack_state);
			}
			if ( $line =~ m/^\#else\s*$/ ) {
			    $preproc_stack_state[-1] = 1-$preproc_stack_state[-1];
			}
			$active = 1;
			for(my $i=0;$i<scalar(@preproc_stack_state);++$i) {
			    my $this_active;
			    if ( $preproc_stack_state[$i] == 1 ) {
				$this_active = 0;
				foreach my $preproc ( @preprocs ) {
				    if ( $preproc eq $preproc_stack_name[$i] ) {$this_active = 1};
				}
			    } else {
				$this_active = 1;
				foreach my $preproc ( @preprocs ) {
				    if ( $preproc eq $preproc_stack_name[$i] ) {$this_active = 0};
				}
			    }
			    if ( $this_active == 0 ) {$active=0};
			}
		    }
		    if ( $active == 1 ) {
# Locate any lines which use the "use" statement and extract the name of the file they use. (We ignore the "hdf5" module as it is external.)
			if ( $line =~ m/^\s*use\s+([a-zA-Z0-9_]+)/i ) {
			    my $incfile = $1;
			    $libraryDependencies{$moduleLibararies{lc($incfile)}} = 1
				if ( exists($moduleLibararies{lc($incfile)}) );
			    unless ( $incfile =~ m/^hdf5$/i || exists($ignoreList{lc($incfile)}) ) {
				$incfile .= ".mod ";
				$incfile =~ s/\r//g;
				@incfiles = ( @incfiles, $workDir.lc($incfile));
				$hasuses = 1;
			    }
			}
			if ( $line =~ m/^\s*\!:\s*(.*)$/ || $line =~ m/^\s*\/\/:\s*(.*)$/ ) {
			    my $includes = $1;
			    my @includes = split(/\s+/,$1);
			    push(@extra_includes,@includes);
			    $hasuses = 1;
			}
			if ( $line =~ m/^\s*module / ) {
			    my $startpos = index($line,"module")+7;
			    my $sublen = length($line)-$startpos-1;
			    my $incfile = substr($line,$startpos,$sublen).".mod ";
			    @modfiles = ( @modfiles, lc($incfile));
			}
			if ( $line =~ m/^\s*include\s+(\'|\")([\w\.\-]+)(\'|\")/i ) {
			    my $ifile = $2;
			    (my $Ifile = $ifile) =~ s/\.inc$/.Inc/;
			    if ( -e $sourcedir."/source/".$Ifile ) {
				push(@scanfiles,$sourcedir."/source/".$Ifile);
			    } elsif ( -e $sourcedir."/work/build/".$ifile ) {
				push(@scanfiles,$sourcedir."/work/build/".$ifile);
			    }
			}
		    }
		}
		close($infile);
	    }

	    # Output library file rule.
	    unless ( $oname =~ m/\.Inc$/ ) {
		my $lname = $oname;
		$lname =~ s/.o$/.fl/;
		print $outfile ".".$workDir.$base.$lname,":\n";
		if ( scalar(keys(%libraryDependencies)) > 0 ) {
		    my $direct = ">";
		    foreach my $library ( keys(%libraryDependencies) ) {
			print $outfile "\t\@echo ".$library." ".$direct." .$workDir$base$lname\n";
			$direct = ">>";
		    }
		} else {
		    print $outfile "\t\@touch .$workDir$base$lname\n";
		}
	    }
	    
	    # Process output for files which had use statements
	    if ( $hasuses == 1 ) {
		# Sort the list of used files
		my @sortedinc = sort @incfiles;
		# Remove any duplicate entries
		@sortedinc = uniq(@sortedinc);
		
		for (my $i = 0; $i < scalar(@sortedinc); $i += 1) {
		    foreach my $item (@modfiles) {
			if ( $sortedinc[$i] eq $item ) {
			    if ( $i < scalar(@sortedinc)-1 ) {
				for (my $j = $i; $j < scalar(@sortedinc)-1; $j += 1) {
				    $sortedinc[$j] = $sortedinc[$j+1];
				}
			    }
			    $i -= 1;
			    pop(@sortedinc);
			    last;
			}
		    }
		}

                # Output the dependencies
		if (scalar(@sortedinc) > 0 || scalar(@extra_includes) > 0) {
		    print $outfile ".".$workDir.$base.$oname,": ";
		    if ( scalar(@sortedinc) > 0 ) {print $outfile ".",join(".",@sortedinc)};
		    print $outfile " Makefile";
		    foreach my $extra_include ( @extra_includes ) {
			print $outfile " $extra_include";
		    }
		    print $outfile "\n";
		    for (my $i = 0; $i < scalar(@sortedinc); $i += 1) {
			$sortedinc[$i] =~ s/.mod\s$/.mod.d /;
		    }
		    my $dname = $oname;
		    $dname =~ s/.o$/.d/;
		    $dname =~ s/.Inc$/.d/;
		    print $outfile ".".$workDir.$base.$dname,": ";
		    if ( scalar(@sortedinc) > 0 ) {print $outfile ".",join(".",@sortedinc)};
		    foreach my $extra_include ( @extra_includes ) {
			(my $dFile = $extra_include) =~ s/\.o$/.d/;
			print $outfile " ".$dFile;
		    }
		    print $outfile "\n";
		    print $outfile "\t\@echo .$workDir$base$oname > .$workDir$base$dname\n";
		    foreach my $extra_include ( @extra_includes ) {
			(my $dFile = $extra_include) =~ s/\.o$/.d/;
			if ( $extra_include =~ m/\// ) {
			    print $outfile "\t\@cat $dFile >> .$workDir$base$dname\n";
			} else {
			    print $outfile "\t\@cat .$workDir$dFile >> .$workDir$base$dname\n";
			}
		    }
		    foreach my $item (@sortedinc) {
			$item =~ s/\s+$//;
			print $outfile "\t\@cat .$item >> .$workDir$base$dname\n";
		    }
		    print $outfile "\t\@sort -u .$workDir$base$dname -o .$workDir$base$dname\n\n";

		    # Create rules for making dependency trees with GraphViz.
		    for (my $i = 0; $i < scalar(@sortedinc); $i += 1) {
			$sortedinc[$i] =~ s/\.d$/.gv /;
		    }
		    my $gvname = $pname.".gv";
		    $dname = $oname;
		    $dname =~ s/.o$/.d/;
		    print $outfile ".".$workDir.$base.$gvname,": .".$workDir.$base.$dname." ";
		    if ( scalar(@sortedinc) > 0 ) {print $outfile ".",join(".",@sortedinc)};
		    print $outfile "\n";
		    print $outfile "\t\@echo \\\"$pname\\\" > .$workDir$base$gvname\n";
		    foreach my $extra_include ( @extra_includes ) {
			(my $dFile = $extra_include) =~ s/\.o$/.d/;
			if ( $extra_include =~ m/\// ) {
			    print $outfile "\t\@awk '{print \"\\\"".$pname."\\\" -> \\\"\"\$\$1\"\\\"\"}' $dFile >> .$workDir$base$gvname\n";
			} else {
			    print $outfile "\t\@awk '{print \"\\\"".$pname."\\\" -> \\\"\"\$\$1\"\\\"\"}' .$workDir$dFile >> .$workDir$base$gvname\n";
			}
		    }
		    foreach my $item (@sortedinc) {
			$item =~ s/\s+$//;
			print $outfile "\t\@awk '{print \"\\\"".$pname."\\\" -> \\\"\"\$\$1\"\\\"\"}' .$item >> .$workDir$base$gvname\n";
			print $outfile "\t\@cat `awk '{print \".".$workDir.$base."\"\$\$1\".gv\"}' .$item` >> .$workDir$base$gvname\n";
		    }
		    print $outfile "\t\@sort -u .$workDir$base$gvname -o .$workDir$base$gvname\n\n";
		}
	    }
	}
    }
    closedir($indir);
}
