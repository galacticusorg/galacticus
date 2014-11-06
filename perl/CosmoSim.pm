package CosmoSim;
use strict;
use warnings;
use Cwd;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = cwd()."/";
}
unshift(@INC,$galacticusPath."perl"); 
use WWW::Curl::Easy;
use WWW::Curl::Form;
use XML::Simple;

sub query {
    # Run a query on the CosmoSim server and download results.
    my $queryString = shift;
    my $fileName    = shift;
    # Parse the Galacticus config file if it is present.
    my $sqlUser;
    my $sqlPassword;
    if ( -e $galacticusPath."/galacticusConfig.xml" ) {
	my $xml    = new XML::Simple();
	my $config = $xml->XMLin($galacticusPath."/galacticusConfig.xml");
	if ( exists($config->{'cosmosimDB'}->{'host'}) ) {
	    my %hosts;
	    if ( exists($config->{'cosmosimDB'}->{'host'}->{'name'}) ) {
		$hosts{'default'} = $config->{'cosmosimDB'}->{'host'};
	    } else {
		%hosts = %{$config->{'cosmosimDB'}->{'host'}};
	    }
	    foreach ( keys(%hosts) ) {
		if ( $_ eq $ENV{'HOSTNAME'} || $_ eq "default" ) {
		    $sqlUser     = $hosts{$_}->{'user'    }
		    if ( exists($hosts{$_}->{'user'    }) );
		    $sqlPassword = $hosts{$_}->{'password'}
		    if ( exists($hosts{$_}->{'password'}) );
		    if ( exists($hosts{$_}->{'passwordFrom'}) ) {
			if ( $hosts{$_}->{'passwordFrom'} eq "input" ) {
			    $sqlPassword = <>;
			    chomp($sqlPassword);
			}
		    }
		}
	    }
	}
    }
    die("generateHalosNBody: CosmoSim database username and password must be defined")
	unless ( defined($sqlUser) && defined($sqlPassword) );
    # Get CUrl objects.
    my $xml      = new XML::Simple    ();
    my $curlPost = new WWW::Curl::Easy();
    my $curlGet  = new WWW::Curl::Easy();
    # Create the query job.
    my $createText;
    my $date       = DateTime->now();
    my $createForm = new WWW::Curl::Form();
    $createForm->formadd("query",$queryString);
    $createForm->formadd("queue","long"      );
    $createForm->formadd("table",$date       );
    $curlPost->setopt(CURLOPT_URL           ,"http://www.cosmosim.org/uws/query");
    $curlPost->setopt(CURLOPT_HTTPPOST      ,$createForm                        );
    $curlPost->setopt(CURLOPT_FOLLOWLOCATION,1                                  );
    $curlPost->setopt(CURLOPT_USERPWD       ,$sqlUser.":".$sqlPassword          );
    $curlPost->setopt(CURLOPT_WRITEDATA     ,\$createText                       );
    print "generateHalosNBody: creating CosmoSim UWS job\n";
    die("generateHalosNBody(): failed to create job")
	unless ( $curlPost->perform() == 0 );
    my $query = $xml->XMLin($createText);
    # Submit the job.
    my $submitText;
    my $submitForm = new WWW::Curl::Form();
    $submitForm->formadd("phase","run");
    $curlPost->setopt(CURLOPT_URL      ,"http://www.cosmosim.org/uws/query/".$query->{'uws:jobId'});
    $curlPost->setopt(CURLOPT_HTTPPOST ,$submitForm                                               );
    $curlPost->setopt(CURLOPT_WRITEDATA,\$submitText                                              );
    print "generateHalosNBody: submitting CosmoSim UWS job\n";
    die("generateHalosNBody(): failed to submit job")
	unless ( $curlPost->perform() == 0 );
    # Check status.
    my $statusText;
    do {
	undef($statusText);
	sleep(10);
	$curlGet->setopt(CURLOPT_URL      ,"http://www.cosmosim.org/uws/query/".$query->{'uws:jobId'}."/phase");
	$curlGet->setopt(CURLOPT_USERPWD  ,$sqlUser.":".$sqlPassword                                          );
	$curlGet->setopt(CURLOPT_WRITEDATA,\$statusText                                                       );
	die("generateHalosNBody(): failed to check status")
	    unless ( $curlGet->perform() == 0 );
	print "generateHalosNBody: CosmoSim UWS job status is '".$statusText."'\n";
	die
	    if ( $statusText eq "ERROR" || $statusText eq "ABORTED" );
    }
    until ( $statusText eq "COMPLETED" );
    # Get results information.
    my $resultsText;
    $curlGet->setopt(CURLOPT_URL      ,"http://www.cosmosim.org/uws/query/".$query->{'uws:jobId'}."/results");
    $curlGet->setopt(CURLOPT_USERPWD  ,$sqlUser.":".$sqlPassword                                            );
    $curlGet->setopt(CURLOPT_WRITEDATA,\$resultsText                                                        );
    die("generateHalosNBody(): failed to get results information")
	unless ( $curlGet->perform() == 0 );
    my $results = $xml->XMLin($resultsText);
    # Download results.
    open(my $resultsFile,">".$fileName);
    $curlGet->setopt(CURLOPT_URL      , $results->{'uws:result'}->{'xlink:href'});
    $curlGet->setopt(CURLOPT_USERPWD  ,$sqlUser.":".$sqlPassword          );
    $curlGet->setopt(CURLOPT_WRITEDATA,$resultsFile                             );
    print "generateHalosNBody: downloading CosmoSim data\n";
    die("generateHalosNBody(): failed to get results information")
	unless ( $curlGet->perform() == 0 );
    close($resultsFile);
    # Delete the job.
    my $deleteText;
    $curlGet->setopt(CURLOPT_URL          , "http://www.cosmosim.org/uws/query/".$query->{'uws:jobId'});
    $curlGet->setopt(CURLOPT_USERPWD      ,$sqlUser.":".$sqlPassword                                  );
    $curlGet->setopt(CURLOPT_CUSTOMREQUEST, "DELETE"                                                  ); 
    $curlGet->setopt(CURLOPT_WRITEDATA    ,\$deleteText                                               );
    print "generateHalosNBody: deleting CosmoSim UWS job\n";
    die("generateHalosNBody(): failed to get results information")
	unless ( $curlGet->perform() == 0 );
}

1;
