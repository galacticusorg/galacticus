#!/usr/bin/env perl
use strict;
use warnings;

# Make a copy of all arguments.
my @arguments = @ARGV;

# Find the executable.
my $executable;
for(my $i=0;$i<scalar(@arguments);++$i) {
    $executable = $arguments[$i+1]
	if ( $arguments[$i] eq "-o" );
}
die('unable to determine executable')
    unless ( defined($executable) );
print "Executable is '".$executable."'\n";

# Get the compile command.
my $compileCommand;
while ( scalar(@arguments) > 0 ) {
    my $arg = shift(@arguments);
    if ( $arg =~ m/^`/ ) {
	my $toExpand = $arg;
	while ( $toExpand !~ m/`$/ ) {
	    $arg = shift(@arguments);
	    $toExpand .= " ".$arg;
	}
	$toExpand =~ s/`//g;
	my $expanded = `$toExpand`;
	$expanded =~ s/\n/ /g;
	$compileCommand .= " ".$expanded;
    } else {
	$compileCommand .= " ".$arg;
    }
}

# Find dynamic libraries being linked.
my $isGCC      = 0;
my $isGFortran = 0;
my $isGPP      = 0;
my @mvLibs        ;
my $mvDirName     ;
open(my $otool,"otool -L ".$executable." |");
while ( my $line = <$otool> ) {
    my @columns = split(" ",$line);
    # Skip unless this is a dylib or an so. MacOS shared libraries should have a .dylib extension, but some libraries that we
    # install (e.q. qhull) install with a .so extension, so we need to handle that case also.
    next
	unless ( $columns[0] =~ m/\.dylib$/ || $columns[0] =~ m/\.so[0-9\.]+$/ );
    # Check for existance of the corresponding static library.
    my $dynamicName = $columns[0];
    (my $libraryName = $dynamicName) =~ s/^.*\/lib([a-zA-Z0-9_\-\+]+)\..*/$1/;
    my $libraryNameOriginal = $libraryName;
    my $staticName;
    if ( $columns[0] =~ m/\.dylib$/ ) {
	($staticName = $dynamicName) =~ s/(\.\d+)?\.dylib$/.a/;
    } elsif ( $columns[0] =~ m/\.so[0-9\.]+$/ ) {
	($staticName = $dynamicName) =~ s/\.so[0-9\.]+$/.a/;
    }
    $libraryName = "qhullstatic_r"
	if ( $libraryName eq "qhull_r" );
    print "Looking for static library for '".$libraryName."'\n";
    if      ( $libraryName =~ m/^gcc/      ) {
	$isGCC      = 1;
	print " -> Can use gcc compiler option\n";
    } elsif ( $libraryName =~ m/^gfortran/ ) {
	$isGFortran = 1;
	print " -> Can use gfortran compiler option\n";
    } elsif ( $libraryName =~ m/^stdc\+\+/ ) {
	$isGPP      = 1;
	print " -> Can use g++ compiler option\n";
    } elsif ( -e $staticName ) {
	print " -> Found static library at '".$staticName."'\n";
	$libraryName =~ s/\+/\\\+/g;
	if ( $compileCommand =~ m/\-l$libraryNameOriginal/ ) {
	    $compileCommand =~ s/\-l$libraryNameOriginal/$staticName/;
	} else {
	    $compileCommand .= " ".$staticName;
	    if ( $libraryName eq "quadmath" ) {
		($mvDirName = $dynamicName) =~ s/(.*)\/.*/$1/;
		opendir(my $libDir,$mvDirName);
		while ( my $fileName = readdir($libDir) ) {
		    if ( $fileName =~ m/^libquadmath.*\.dylib/ ) {
			push(@mvLibs,$fileName);
		    }
		}
		close($libDir);
		my $mvCommand = "sudo -- sh -c '".join("; ",map {"mv ".$mvDirName."/".$_." ".$mvDirName."/".$_."~"} @mvLibs)."'";
		print "Must move dylibs temporarily (requires sudo):\n";
		system($mvCommand);
	    }
	}
    } else {
	# Look in other possible locations.
	my $found     = 0;
	my @locations = ( "/usr/local/lib" );
	push(@locations,split(/:/,$ENV{'LD_LIBRARY_PATH'}))
	    if ( exists($ENV{'LD_LIBRARY_PATH'}) );
	foreach my $location ( @locations ) {
	    my $fileName = $location."/lib".$libraryName.".a";
	    if ( -e $fileName ) {
		print " -> Found static library at '".$fileName."'\n";
		if ( $compileCommand =~ m/\-l$libraryNameOriginal/ ) {
		    $compileCommand =~ s/\-l$libraryNameOriginal/$staticName/;
		} else {
		    $compileCommand .= " ".$fileName;
		}
		$found = 1;
		last;
	    }
	}
	unless ( $found ) {
	    (my $path = $staticName) =~ s/\/[^\/]+$//;
	    print " -> No static library found\n";
	}
    }
}
close($otool);
# Add static GCC flags.
$compileCommand .= " -static-libgfortran"
    if ( $isGFortran );
$compileCommand .= " -static-libgcc"
    if ( $isGCC      );
$compileCommand .= " -static-libstdc++"
    if ( $isGPP      );
# Compile the static binary.
print "Relinking with: ".$compileCommand."\n";
system($compileCommand);
# Restore dylibs.
if ( scalar(@mvLibs) > 0 ) {
    my $mvCommand = "sudo -- sh -c '".join("; ",map {"mv ".$mvDirName."/".$_."~ ".$mvDirName."/".$_} @mvLibs)."'";
    print "Must restore temporarily moved dylibs (requires sudo):\n";
    system($mvCommand);
}
exit;
