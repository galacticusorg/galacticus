package SystemRedirect;
  my $status = 1;
  $status;

sub tofile {
    ($command, $file) = @_;

    # Save current standard output and standard errer    
    open(OLDOUT, ">&STDOUT");
    open(OLDERR, ">&STDERR");
    
    # Open new standard output and error
    open(STDOUT, ">$file") || die "Can't redirect stdout";
    open(STDERR, ">&STDOUT") || die "Can't dup stdout";
    
    select(STDERR); $| = 1;     # make unbuffered
    select(STDOUT); $| = 1;     # make unbuffered
    
    # Run the system command.
    system($command);
    $result = $?;
    
    # Close the output file
    close(STDOUT);
    close(STDERR);
    
    # Restore standard output and error
    open(STDOUT, ">&OLDOUT");
    open(STDERR, ">&OLDERR");

    return $result;
}
