#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
 $ENV{"GALACTICUS_ROOT_V092"} = getcwd()."/";
}
unshift(@INC,$galacticusPath."perl"); 
use threads;

# A simple wrapper script which launches multiple threads to process the input list of files through Grasil in parallel.
# Andrew Benson (19-July-2011)

# Specify default options.
my $cpuLimit = 3600;

# Get a list of the file root names to run through Grasil.
my @grasilFilesRoots = @ARGV;
while ( $grasilFilesRoots[0] =~ m/^\-\-/ ) {
    my $optionName  = shift(@grasilFilesRoots);
    my $optionValue = shift(@grasilFilesRoots);
    $cpuLimit = $optionValue
	if ( $optionName eq "--cpuLimit" );
}

# Ensure that we have the Grasil executable and data files.
my $grasilPath = $galacticusPath."/aux/Grasil/";
system("mkdir -p ".$grasilPath);
system("wget http://adlibitum.oat.ts.astro.it/silva/grasil/download/gslib.tar.gz -O ".$grasilPath."/gslib.tar.gz")
    unless ( -e $grasilPath."/gslib.tar.gz" );
die("Grasil_Wrapper.pl: unable to download gslib.tar.gz")
    unless ( $? == 0 );
unless ( -e $grasilPath."/TAU96.OF" ) {
    system("cd ".$grasilPath."; tar xvfz gslib.tar.gz");
    die("Grasil_Wrapper.pl: unable to extract gslib.tar.gz")
	unless ( $? == 0 );
    unlink($grasilPath."/grasil")
	if ( -e $grasilPath."/grasil" );
}
system("wget http://adlibitum.oat.ts.astro.it/silva/grasil/SSP_zip/grasil.tar.gz -O ".$grasilPath."/grasil.tar.gz")
    unless ( -e $grasilPath."/grasil.tar.gz" );
die("Grasil_Wrapper.pl: unable to download grasil.tar.gz")
    unless ( $? == 0 );
unless ( -e $grasilPath."/template_grasil.sf" ) {
    system("cd ".$grasilPath."; tar xvfz grasil.tar.gz");
    die("Grasil_Wrapper.pl: unable to extract gslib.tar.gz")
	unless ( $? == 0 );
    unlink($grasilPath."/grasil")
	if ( -e $grasilPath."/grasil" );
}
system("wget http://users.obs.carnegiescience.edu/abenson/galacticus/tools/grasil -O ".$grasilPath."/grasil; chmod u=wrx ".$grasilPath."/grasil")
    unless ( -e $grasilPath."/grasil" );
die("Grasil_Wrapper.pl: unable to download grasil")
    unless ( $? == 0 );

# Create an array of threads.
my @grasilThreads;

if ( scalar(@grasilFilesRoots) > 1 ) {

    # Launch a thread for each given name.
    foreach my $grasilFilesRoot ( @grasilFilesRoots ) {
	if ( $^V lt v5.14.0 ) { 
	    $grasilThreads[++$#grasilThreads] = threads->create( sub {system("cd `dirname ".$grasilFilesRoot."`; ulimit -t ".$cpuLimit."; ".$grasilPath."/grasil `basename ".$grasilFilesRoot."` &> `basename ".$grasilFilesRoot.".log`");});
	} else {
	    $grasilThreads[++$#grasilThreads] = threads->create({'context' => 'void'}, sub {system("cd `dirname ".$grasilFilesRoot."`; ulimit -t ".$cpuLimit."; ".$grasilPath."/grasil `basename ".$grasilFilesRoot."` &> `basename ".$grasilFilesRoot.".log`");});
	}
    }
    
    # Wait for threads to finish.
    foreach my $grasilThread ( @grasilThreads ) {
	$grasilThread->join();
    }

} else {

    # Run a single instance without launching a thread.
    system("cd `dirname ".$grasilFilesRoots[0]."`; ulimit -t ".$cpuLimit."; ".$grasilPath."/grasil `basename ".$grasilFilesRoots[0]."` &> `basename ".$grasilFilesRoots[0].".log`");

}

exit;
