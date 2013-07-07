#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 

# Locate source files which have dependencies on include files

# Define the source directory
if ( $#ARGV != 0 ) {die "Usage: Find_Include_Dependencies.pl sourcedir"};
my $sourcedir = $ARGV[0];

#
# Open an output file
open(my $outfile,">$sourcedir/work/build/Makefile_Include_Deps");
#
# Build a list of source directories.
my @sourcedirs = ( $sourcedir."/source" );
my @bases      = ( ""                   );
if ( -e $sourcedir."/Source_Codes" ) {
    opendir(my $sdir,"$sourcedir/Source_Codes");
    while ( my $line = readdir $sdir ) {
	$line =~ s/\n//;
	if ( -d $sourcedir."/Source_Codes/".$line ) {
	    unless ( $line =~ m/^\.+$/ ) {
		$sourcedirs[++$#sourcedirs] = $sourcedir."/Source_Codes/".$line;
		$bases[++$#bases] = "Source_Codes/".$line."/";
	    }
	}
    }
    closedir($sdir);
}

# Array of all auto-generated include files.
my @All_Auto_Includes;

#
# Open the source directory
#
my $ibase = -1;
my %nonAutoInclude;
foreach my $srcdir ( @sourcedirs ) {
    ++$ibase;
    my $base = $bases[$ibase];
    opendir(my $indir,$srcdir) or die "Can't open the source directory: #!";
    while (my $fname = readdir $indir) {
	
	if ( ( lc($fname) =~ m/\.f(90)??t??$/ || lc($fname) =~ m/\.inc$/ || lc($fname) =~ m/\.c(pp)??$/ || lc($fname) =~ m/\.h??$/ ) && lc($fname) !~ m/^\.\#/ ) {
	    my $pname = $fname;
	    my $fullname = "$srcdir/$fname";
	    my $doesio = 0;
	    open(my $infile,$fullname) or die "Can't open input file: #!";
	    my @fileinfo = stat $infile;
	    
	    my $hasincludes = 0;
	    my $oname = $fname;
	    $oname =~ s/\.F90$/\.o/;
	    $oname =~ s/\.f90t$/\.o/;
	    $oname =~ s/\.F90t$/\.o/;
	    $oname =~ s/\.f$/\.o/;
	    $oname =~ s/\.F$/\.o/;
	    $oname =~ s/\.c(pp)$/\.o/;
	    $oname =~ s/\.h$/\.o/;
	    my @incfiles;
	    while (my $line = <$infile>) {
		# Locate any lines which use the "include" statemalent and extract the name of the file they include
		if ( $line =~ m/^\s*#??include\s*['"](.+\.(inc|h)\d*)['"]/ ) {
		    my $incfile = $1;
		    $nonAutoInclude{$incfile} = 1 if ( $line =~ m/! NO_USES/ );
		    $incfile .= " ";
		    if ( $hasincludes == 0 ) {$hasincludes = 1};
		    @incfiles = ( @incfiles, $incfile);
		}
	    }
	    # Process output for files which had include statements
	    if ( $hasincludes == 1 ) {
		# Sort the list of included files
		my @sortedinc = sort @incfiles;
		# Remove any duplicate entries
		if ($#sortedinc > 0) {
		    for (my $i = 1; $i <= $#sortedinc; $i += 1) {
			if ($sortedinc[$i] eq $sortedinc[$i-1]) {
			    if ( $i < $#sortedinc ) {
				for (my $j = $i; $j < $#sortedinc; $j += 1) {
				    $sortedinc[$j] = $sortedinc[$j+1];
				}
			    }
			    $i -= 1;
			    $#sortedinc -= 1;
			}
		    }
		}
		# Output the dependencies
		print $outfile "./work/build/".$base,$oname,":";
		foreach my $inc ( @sortedinc ) {
		    $inc =~ s/\s*$//;
		    (my $Iinc = $inc) =~ s/\.inc$/\.Inc/;
		    my $ext_inc = $srcdir."/".$inc;
		    my $ext_Iinc = $srcdir."/".$Iinc;
		    if ( -e $ext_Iinc ) {
			if ( $ibase == 0 ) {
			    print $outfile " ./work/build/$inc";
			    if ( $inc =~ m/\.inc$/ && ! exists($nonAutoInclude{$inc}) ) {$All_Auto_Includes[++$#All_Auto_Includes]=$Iinc};
			} else {
			    print $outfile " ./work/build/$srcdir/$inc";
			    if ( $inc =~ m/\.inc$/ && ! exists($nonAutoInclude{$inc}) ) {$All_Auto_Includes[++$#All_Auto_Includes]=$srcdir."/".$Iinc};			}
		    } elsif ( -e $ext_inc ) {
			if ( $ibase == 0 ) {
			    print $outfile " ./work/build/$inc";
			    if ( $inc =~ m/\.inc$/ && ! exists($nonAutoInclude{$inc}) ) {$All_Auto_Includes[++$#All_Auto_Includes]=$inc};
			} else {
			    print $outfile " ./work/build/$srcdir/$inc";
			    if ( $inc =~ m/\.inc$/ && ! exists($nonAutoInclude{$inc}) ) {$All_Auto_Includes[++$#All_Auto_Includes]=$srcdir."/".$inc};
			}
		    } else {
			print $outfile " ./work/build/$inc";
			if ( $inc =~ m/\.inc$/ && ! exists($nonAutoInclude{$inc}) ) {$All_Auto_Includes[++$#All_Auto_Includes]=$inc};
		    }
		}
		print $outfile "\n\n";
	    }
	    close($infile);
	}
    }
    closedir($indir);
}
print $outfile "\n./work/build/Makefile_Use_Deps: ./work/build/".join(" ./work/build/",@All_Auto_Includes)."\n";
close($outfile);
