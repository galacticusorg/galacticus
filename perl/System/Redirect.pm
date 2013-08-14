package SystemRedirect;
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

sub tofile {
    (my $command, my $file) = @_;
    
    # Save current standard output and standard errer    
    open OLDOUT, '>&', \*STDOUT or die "Can't duplicate STDOUT: $!";
    open OLDERR, '>&', \*STDERR or die "Can't duplicate STDERR: $!";
    
    # Open new standard output and error
    open(STDOUT, ">$file") || die "Can't redirect stdout";
    open(STDERR, ">&STDOUT") || die "Can't dup stdout";
    
    select(STDERR); $| = 1;     # make unbuffered
    select(STDOUT); $| = 1;     # make unbuffered
    
    # Run the system command.
    system($command);
    my $result = $?;
    
    # Close the output file
    close(STDOUT);
    close(STDERR);
    
    # Restore standard output and error
    open STDOUT, '>&', \*OLDOUT or die "Cannot duplicate OLDOUT: $!";
    open STDERR, '>&', \*OLDERR or die "Cannot duplicate OLDERR: $!";

    return $result;
}

1;
