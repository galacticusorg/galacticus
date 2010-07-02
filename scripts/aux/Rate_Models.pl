#!/usr/bin/env perl
use Image::Magick;
use XML::Simple;
use Data::Dumper;
use Term::ReadKey;
use File::Find;
use File::Copy;

# Script to display Galacticus model plots and request a user rating for them. Rating is stored to file.
# Andrew Benson (14-June-2010)

# Set default list of model directories to scan.
$parameters{'modelDirectories'} = "models";

# Set default name of plot to rate.
$parameters{'plotToRate'}       = "Star_Formation_History.pdf";

# Name of the corresponding fit.
$parameters{'fitName'}          = "Volume averaged star formation rate history";

# Name of the rater.
$parameters{'userName'}         = getlogin();

# Parse command line parameters.
for($iParameter=0;$iParameter<=$#ARGV;++$iParameter) {
    $parameter = $ARGV[$iParameter];
    if ( $parameter =~ m/\"/ ) {
	do {
	    ++$iParameter;
	    $parameter .= $ARGV[$iParameter];
	} until ( $parameter =~ m/\"$/ );
    }
    if ( $parameter =~ m/(\S+)=(.*)/ ) {
	$parameterName  = $1;
	$parameterValue = $2;
	$parameters{$parameterName} = $parameterValue;
    }
}

# Specify directory containing models.
@modelDirectories = split(/,/,$parameters{'modelDirectories'});

# Name of plot files to rate.
$plotToRate = $parameters{'plotToRate'};

# Name of fitting entry to rate.
$fitName = $parameters{'fitName'};

# Name of the rating user.
$userName = $parameters{'userName'};

# Scan for model directories.
find(\&processFile, @modelDirectories);

exit;

sub processFile {

    $fileName = $_;
    chomp($fileName);

    # Check if the file matches the name of the image file we're searching for.
    if ( $fileName eq $plotToRate.".bz2" ) {
	system("bunzip2 ".$plotToRate.".bz2");
	$fileName = $plotToRate;
	$reCompressFile = 1;
    } else {
	$reCompressFile = 0;
    }
    if ( $fileName eq $plotToRate ) {

        # Display the image.
	$modelName = $File::Find::dir;
	print "Displaying for model: ".$modelName."\n";
	$image = Image::Magick->new;
	open(IMAGE,$fileName);
	$image->Read(file=>\*IMAGE);
	close(IMAGE);
	$image->Display();
	if ( $reCompressFile == 1 ) {system("bzip ".$fileName)};

        # Get user rating.
	ReadMode 4; # Turn off controls keys
	$rating = -1;
	while ( $rating == -1 ) {
	    while (not defined ($key = ReadKey(-1))) {
		# No key yet
	    }
	    if ( $key =~ m/^[12345]$/ ) {
		$rating = $key;
	    } else {
		print "Please enter rating of 1-5.\n";
	    }
	}
	ReadMode 0;
	print " -> Rating: ".$rating."\n";

	# See if we can find a fit results file.
	if ( -e "galacticusFits.xml.bz2" ) {
	    system("bunzip2 galacticusFits.xml.bz2");
	    $reCompressFile = 1;
	} else {
	    $reCompressFile = 0;
	}
	if ( -e "galacticusFits.xml" ) {
	    $xml = new XML::Simple;
	    $data = $xml->XMLin("galacticusFits.xml");
	    $netRating      = 0;
	    $netRatingCount = 0;
	    foreach $fit ( keys(%{$data->{'galacticusFit'}}) ) {
		if ( $fit =~ m/$fitName/ ) {
		    $now = `date`;
		    chomp($now);
		    $data->{'galacticusFit'}->{$fit}->{'rating'}->{'value'} = $rating;
		    $data->{'galacticusFit'}->{$fit}->{'rating'}->{'user'}  = $userName;
		    $data->{'galacticusFit'}->{$fit}->{'rating'}->{'time'}  = $now;
		}
		if ( exists($data->{'galacticusFit'}->{$fit}->{'rating'}->{'value'}) && $fit ne "net" ) {
		    ++$netRatingCount;
		    $netRating += $data->{'galacticusFit'}->{$fit}->{'rating'}->{'value'};
		}
	    }
	    # Store the net rating for this model.
	    foreach $fit ( keys(%{$data->{'galacticusFit'}}) ) {
		if ( $fit eq "net" ) {
		    $data->{'galacticusFit'}->{$fit}->{'rating'}->{'value'} = $netRating/$netRatingCount;
		    $data->{'galacticusFit'}->{$fit}->{'rating'}->{'user'}  = $userName;
		    $data->{'galacticusFit'}->{$fit}->{'rating'}->{'time'}  = $now;
		}
	    }
	    move("galacticusFits.xml","galacticusFits.xml~");
	    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFits");
	    open(outHndl,">galacticusFits.xml");
	    print outHndl $xmlOutput->XMLout($data);
	    close(outHndl);
	    if ( $reCompressFile == 1 ) {system("bzip2 galacticusFits.xml")};
	}
    }
}
