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
my $compileCommand = join(" ",@arguments);

# Find dynamic libraries being linked.
my $isGCC      = 0;
my $isGFortran = 0;
my $isGPP      = 0;
my @mvLibs        ;
my $mvDirName     ;
open(my $otool,"otool -L ".$executable." |");
while ( my $line = <$otool> ) {
    my @columns = split(" ",$line);
    # Skip unless this is a dylib.
    next
	unless ( $columns[0] =~ m/\.dylib$/ );
    # Check for existance of the corresponding static library.
    my $dynamicName = $columns[0];
    (my $libraryName = $dynamicName) =~ s/^.*\/lib([a-zA-Z0-9_\-\+]+)\..*/$1/;
    (my $staticName  = $dynamicName) =~ s/(\.\d+)?\.dylib$/.a/;
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
	if ( $compileCommand =~ m/\-l$libraryName/ ) {
	    $compileCommand =~ s/\-l$libraryName/$staticName/;
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
	(my $path = $staticName) =~ s/\/[^\/]+$//;
	print " -> No static library found, looked in '".$path."'\n";
	system("ls ".$path);
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
