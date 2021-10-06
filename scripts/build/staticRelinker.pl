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
    (my $staticName = $dynamicName) =~ s/(\.\d+)?\.dylib$/.a/;
    if ( $libraryName =~ m/^gcc*/ ) {
	$isGCC = 1;
    } elsif ( $libraryName =~ m/^gfortran/ ) {
	$isGFortran = 1;
    } elsif ( -e $staticName ) {
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
   }
}
close($otool);
# Add static GCC flags.
$compileCommand .= " -static-libgfortran"
    if ( $isGFortran );
$compileCommand .= " -static-libgcc"
    if ( $isGCC      );
# Compile the static binary.
system($compileCommand);
# Restore dylibs.
if ( scalar(@mvLibs) > 0 ) {
    my $mvCommand = "sudo -- sh -c '".join("; ",map {"mv ".$mvDirName."/".$_."~ ".$mvDirName."/".$_} @mvLibs)."'";
    print "Must restore temprarily moved dylibs (requires sudo):\n";
    system($mvCommand);
}
exit;
