#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use Clone qw(clone);
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Switch;
require Galacticus::Constraints::Parameters;
require Galacticus::HDF5;

# Run calculations to determine the model discrepancy arising from the use of Monte Carlo merger trees.
# Andrew Benson (16-November-2012)

# Get arguments.
die("Usage: monteCarloTrees.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
     make  => "yes"           ,
     trees => "millennium:MPA"
    );
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Parse the constraint config file.
my $xml    = new XML::Simple;
my $config = $xml->XMLin($configFile, KeyAttr => 0);

# Validate the config file.
die("monteCarloTrees.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'workDirectory' }) );
die("monteCarloTrees.pl: compilation must be specified in config file"   ) unless ( exists($config->{'compilation'   }) );
die("monteCarloTrees.pl: baseParameters must be specified in config file") unless ( exists($config->{'baseParameters'}) );

# Determine the scratch and work directories.
my $workDirectory    = $config->{'workDirectory'};
my $scratchDirectory = $config->{'workDirectory'};
$scratchDirectory    = $config->{'scratchDirectory'} if ( exists($config->{'scratchDirectory'}) );

# Create the work and scratch directories.
system("mkdir -p ".$config->{'workDirectory'});

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("monteCarloTrees.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Determine the number of subvolumes to obtain.
my $subVolumeCount;
$subVolumeCount    = $arguments{'subVolumeCount'}
    if ( exists($arguments{'subVolumeCount'}) );
die('monteCarloTrees.pl: subVolumeCount must be a power of 2')
    unless ( ! defined($subVolumeCount) || ($subVolumeCount & ($subVolumeCount-1)) == 0 );

# Get all required merger trees from the Millennium Simulation database.
my $millenniumDbCommand;
my $subvolumesMaximum;
switch ( $arguments{'trees'} ) {
    case ( 'millennium:MPA' ) {
        $subvolumesMaximum = 512;
	$millenniumDbCommand = $galacticusPath."/scripts/aux/Millennium_Trees_Grab.pl --table MPAHaloTrees..MHalo --traceParticles no";
	# Get Millennium database username and password if available.
	$millenniumDbCommand .= " --user "    .$arguments{'user'    } 
	if ( exists($arguments{'user'    }) );
	$millenniumDbCommand .= " --password ".$arguments{'password'} 
	if ( exists($arguments{'password'}) );
    }
    case ( 'MillGas' ) {
	$subvolumesMaximum = 0;
    }
}
if ( $subvolumesMaximum == 0 ) {
    die("monteCarloTrees.pl: no subvolumes available for these trees")
	if ( defined($subVolumeCount) );
} else {
    die('monteCarloTrees.pl: subVolumeCount must be between 4 and '.$subvolumesMaximum.' inclusive')
	if ($subVolumeCount < 4 || $subVolumeCount > $subvolumesMaximum);
}

# Determine the path to which merger trees should be stored.
my $treeDirectory     = $workDirectory;
if ( -e $galacticusPath."/galacticusConfig.xml" ) {
    my $xml    = new XML::Simple;
    my $config = $xml->XMLin($galacticusPath."/galacticusConfig.xml");
    switch ( $arguments{'trees'} ) {
	case ( 'millennium:MPA' ) {
	    if ( exists($config->{'millenniumDB'}->{'host'}) ) {
		foreach ( keys(%{$config->{'millenniumDB'}->{'host'}}) ) {
		    $treeDirectory = $config->{'millenniumDB'}->{'host'}->{$_}->{'treePath'}
		    if ( ( $ENV{'HOSTNAME'} =~ m/$_/ || $_ eq "default" ) && exists($config->{'millenniumDB'}->{'host'}->{$_}->{'treePath'}) );
		}
	    }
	    $treeDirectory .= "/MPAHalo"
		if ( -e $treeDirectory."/MPAHalo" );
	}
	case ( 'MillGas' ) {
	    if ( exists($config->{'millGas'}->{'host'}) ) {
		foreach ( keys(%{$config->{'millGas'}->{'host'}}) ) {
		    $treeDirectory = $config->{'millGas'}->{'host'}->{$_}->{'treePath'}
		    if ( ( $ENV{'HOSTNAME'} =~ m/$_/ || $_ eq "default" ) && exists($config->{'millGas'}->{'host'}->{$_}->{'treePath'}) );
		}
	    }
	    $treeDirectory .= "/treedir_976"
		if ( -e $treeDirectory."/treedir_976" );	   
	}
    }
}

# Retrieve merger trees from the Millennium database.
if ( defined($millenniumDbCommand) ) {
    system("cd ".$galacticusPath."; make Millennium_Merger_Tree_File_Maker.exe");
    my $failuresEncountered     = 0;
    my $subVolumeRetrievedCount = -1;
    while ( $subVolumeRetrievedCount < $subVolumeCount-1 ) {
	++$subVolumeRetrievedCount;
	# Check if we already have this file.
	unless ( -e $treeDirectory."/subvolume".$subVolumeRetrievedCount.".hdf5" ) {
	    # Retrieve the raw data.
	    my $thisMillenniumDbCommand = $millenniumDbCommand." --select \"root.fileNr = ".$subVolumeRetrievedCount." and root.snapNum = 63\" --output ".$treeDirectory."/subvolume".$subVolumeRetrievedCount.".csv";
	    system($thisMillenniumDbCommand);
	    unless ( $? == 0 ) {
		print "monteCarloTrees.pl: failed to retrieve subvolume ".$subVolumeRetrievedCount." from Millennium database\n";
		++$failuresEncountered;
	    }
	    # Convert to Galacticus format.
	    system($galacticusPath."/Millennium_Merger_Tree_File_Maker.exe ".$treeDirectory."/subvolume".$subVolumeRetrievedCount.".csv none ".$treeDirectory."/subvolume".$subVolumeRetrievedCount.".hdf5 galacticus 1");
	    unless ( $? == 0 ) {
		print "monteCarloTrees.pl: failed to convert subvolume ".$subVolumeRetrievedCount." to Galacticus format\n";
		++$failuresEncountered;
	    } else {
		unlink($treeDirectory."/subvolume".$subVolumeRetrievedCount.".csv");
	    }
	}
    }
    die("monteCarloTrees.pl: failures encountered when downloading from Millennium database")
	unless ( $failuresEncountered == 0 );
}

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'compilation'},$config->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Ensure we have redshift zero included in the outputs.
my @outputRedshifts = split(/\s+/,$parameters->{'parameter'}->{'outputRedshifts'}->{'value'});
my $gotRedshiftZero = 0;
foreach ( @outputRedshifts ) {
$gotRedshiftZero = 1
    if ( $_ == 0.0 );
}
push(@outputRedshifts,0.0)
    unless ( $gotRedshiftZero == 1 );
$parameters->{'parameter'}->{'outputRedshifts'}->{'value'} = join(" ",@outputRedshifts);

# Switch off thread locking.
$parameters->{'parameter'}->{'treeEvolveThreadLock'}->{'value'} = "false";

# Initialize a stack for PBS models.
my @pbsStack;

# Define a path name for the models.
(my $modelsDirectory = $workDirectory."/modelDiscrepancyNew/monteCarloTrees.noJumpNoPromoteNoMonotonic.".$arguments{'trees'}."/") =~ s/:/_/g;

# Specify models to run.
my @models;
my $randomSeed = 826;
my $haveSubvolumes = 1;
$haveSubvolumes = 0
    if ( $subvolumesMaximum == 0 );
$subVolumeCount = 1
    unless ( defined($subVolumeCount) );
for(my $iSubvolume=0;$iSubvolume<$subVolumeCount;++$iSubvolume) {
    ++$randomSeed;
    my $subVolumeSuffix = "";
    $subVolumeSuffix = $iSubvolume
	if ( $haveSubvolumes == 1 );
    my @monteCarloParameters =
	(
	 # Switch to using Monte-Carlo trees.
	 {
	     name  => "mergerTreeConstructMethod",
	     value => "build"
	 },
	 {
	     name  => "mergerTreeBuildTreesHaloMassDistribution",
	     value => "read"
	 },
	 {
	     name  => "mergerTreeBuildTreeMassesFile",
	     value => $modelsDirectory."nBody".$subVolumeSuffix."/treeMasses.hdf5"
	 },
	 # Match the mass resolution and timestepping of the N-body merger trees.
	 {
	     name  => "mergerTreesBuildMassResolutionMethod",
	     value => "fixed"
	 },
	 {
	     name  => "mergerTreeRegridTimes",
	     value => "true"
	 },
	 {
	     name  => "haloMassFunctionSimpleSystematicAlpha",
	     value => 0.0
	 },
	 {
	     name  => "haloMassFunctionSimpleSystematicBeta",
	     value => 0.0
	 },
	 {
	     name  => "powerSpectrumIndex",
	     value => 1.0
	 },
	 {
	     name  => "powerSpectrumReferenceWavenumber",
	     value => 1.0
	 },
	 {
	     name  => "powerSpectrumRunning",
	     value => 0.0
	 },
	 {
	     name  => "transferFunctionMethod",
	     value => "file"
	 },
	 # Set mass function analysis to use Poisson model appropriate for N-body mass function sampling.
	 {
	     name  => "analysisMassFunctionCovarianceModel",
	     value => "Poisson"
	 },
	 # Set an initial random number seed.
	 {
	     name  => "randomSeed",
	     value => $randomSeed
	 }
	);
    # Adjust the number of trees to run if specified.
    push(
	@monteCarloParameters,
	{
	    name  => "mergerTreeBuildTreesPerDecade",
	    value => $arguments{'treesPerDecade'}
	}
	)
	if ( exists($arguments{'treesPerDecade'}) );
    # Add tree specific settings.
    switch ( $arguments{'trees'} ) {
	case ( 'millennium:MPA' ) {
	    push(
		@monteCarloParameters,
		{
		    name  => "mergerTreeBuildMassResolutionFixed",
		    value => 2.3578789e10
		},
		{
		    name  => "mergerTreeRegridCount",
		    value => 60
		},
		{
		    name  => "mergerTreeRegridSpacing",
		    value => "millennium"
		},
		# Set cosmological parameters to match the Millennium Simulation.
		{
		    name  => "H_0",
		    value => 73.0
		},
		{
		    name  => "Omega_Matter",
		    value => 0.25
		},
		{
		    name  => "Omega_DE",
		    value => 0.75
		},
		{
		    name  => "Omega_b",
		    value => 0.0455
		},
		{
		    name  => "sigma_8",
		    value => 0.9
		},
		{
		    name  => "transferFunctionFile",
		    value => "data/largeScaleStructure/transferFunctionMillenniumSimulation.xml"
		}
		);
	}
	case ( 'MillGas' ) {
	    push(
		@monteCarloParameters,
		{
		    name  => "mergerTreeBuildMassResolutionFixed",
		    value => 2.6602117120e10
		},
		{
		    name  => "mergerTreeRegridCount",
		    value => 873
		},
		{
		    name  => "mergerTreeRegridSpacing",
		    value => "read"
		},
		{
		    name  => "mergerTreeRegridRedshifts",
		    value => "0.000000000000000 0.00122199999168515 0.00244900002144277 0.00367700005881488 0.00491000013425946 0.00614499999210238 0.00738300010561943 0.00862399954348803 0.00986799970269203 0.0111149996519089 0.0123650003224611 0.0136190000921488 0.0148750003427267 0.01613499969244 0.0173979997634888 0.0173979997634888 0.0186640005558729 0.0199330002069473 0.0212479997426271 0.0225680004805326 0.0238899998366833 0.0252170003950596 0.026545999571681 0.027878999710083 0.0292160008102655 0.0305560007691383 0.0318990014493465 0.0332470014691353 0.034596998244524 0.0359509997069836 0.0373089984059334 0.0386699996888638 0.0400349982082844 0.0414029993116856 0.0428170002996922 0.044234998524189 0.0456550009548664 0.0470810011029243 0.0485100001096725 0.0499439984560013 0.0513809993863106 0.0528220012784004 0.0542659983038902 0.0557149983942509 0.0571690015494823 0.0586259998381138 0.0600869990885258 0.0615510009229183 0.0630199983716011 0.0644930005073547 0.0660089999437332 0.0675309970974922 0.0690559968352318 0.0705860033631325 0.0721189975738525 0.0736579969525337 0.0752009972929955 0.0767470002174377 0.0782999992370605 0.0798550024628639 0.0814170017838478 0.0829809978604317 0.0845519974827766 0.0861250013113022 0.0877050012350082 0.0892879962921143 0.0909119993448257 0.0925429984927177 0.0941770002245903 0.0958169996738434 0.097461998462677 0.0991109982132912 0.100764997303486 0.102425001561642 0.104088999330997 0.105759002268314 0.107432998716831 0.109113000333309 0.110798001289368 0.112488001585007 0.114183999598026 0.115883000195026 0.117623001337051 0.119368001818657 0.121119000017643 0.12287399917841 0.124635003507137 0.126403003931046 0.128175005316734 0.129953995347023 0.131738007068634 0.133526995778084 0.135322004556656 0.137123003602028 0.1389289945364 0.140741005539894 0.142559006810188 0.144382998347282 0.146243005990982 0.148111000657082 0.149983003735542 0.151862993836403 0.153747007250786 0.155638992786407 0.157536000013351 0.159438997507095 0.1613499969244 0.163265004754066 0.165188997983932 0.16711699962616 0.169053003191948 0.170993998646736 0.172942996025085 0.174897998571396 0.176886007189751 0.178882002830505 0.180884003639221 0.182894006371498 0.184909999370575 0.186931997537613 0.188962996006012 0.190999001264572 0.193042993545532 0.195094004273415 0.197152003645897 0.19921800494194 0.201288998126984 0.203367993235588 0.205455005168915 0.207549005746841 0.209674000740051 0.211805000901222 0.213945999741554 0.216094002127647 0.21824799478054 0.220411002635956 0.222581997513771 0.224758997559547 0.226945996284485 0.229139998555183 0.231342002749443 0.233550995588303 0.235769003629684 0.237994998693466 0.240226998925209 0.242468997836113 0.2447379976511 0.247016996145248 0.249301999807358 0.251596003770828 0.253899991512299 0.256209999322891 0.258531004190445 0.26085901260376 0.263195991516113 0.26554200053215 0.267895996570587 0.270258992910385 0.272632002830505 0.275012999773026 0.277404010295868 0.279801994562149 0.282225996255875 0.284657001495361 0.287099003791809 0.289548993110657 0.292010009288788 0.294477999210358 0.29695799946785 0.299445986747742 0.301943987607956 0.304452985525131 0.306968986988068 0.309496998786926 0.312032997608185 0.314579993486404 0.317135989665985 0.319703012704849 0.322290003299713 0.324885994195938 0.327493995428085 0.330110996961594 0.332738012075424 0.335375994443893 0.338023990392685 0.340683996677399 0.343353003263474 0.346033006906509 0.348722994327545 0.351424008607864 0.354135990142822 0.356860995292664 0.359595000743866 0.362340003252029 0.365101993083954 0.367873013019562 0.370656996965408 0.373450994491577 0.376257985830307 0.37907400727272 0.381904989480972 0.384743988513947 0.387598007917404 0.3904629945755 0.393337994813919 0.39622700214386 0.399125993251801 0.402038991451263 0.404962003231049 0.407898992300034 0.410847008228302 0.41380500793457 0.416777998209 0.419761002063751 0.422758013010025 0.425767004489899 0.428790003061295 0.431825995445251 0.434873014688492 0.437934011220932 0.441006988286972 0.444094985723495 0.447194993495941 0.450309991836548 0.453435987234116 0.456577003002167 0.459722012281418 0.46288201212883 0.466053992509842 0.469242006540298 0.472442001104355 0.475657999515533 0.478884994983673 0.482129007577896 0.485385000705719 0.488656014204025 0.491943001747131 0.49524199962616 0.498558014631271 0.501887023448944 0.505232989788055 0.508590996265411 0.511949002742767 0.515323996543884 0.518710970878601 0.522113978862762 0.525532007217407 0.528964996337891 0.532416999340057 0.535880982875824 0.539361000061035 0.542859971523285 0.54637199640274 0.549899995326996 0.553444027900696 0.557003974914551 0.560582995414734 0.564176976680756 0.567762017250061 0.571362972259521 0.574984014034271 0.578619003295898 0.582270979881287 0.585940003395081 0.589625000953674 0.59332799911499 0.597051024436951 0.60078901052475 0.60454398393631 0.608317017555237 0.612107992172241 0.61591899394989 0.619745016098022 0.623589992523193 0.627418994903564 0.63126802444458 0.635133028030396 0.639015972614288 0.642920970916748 0.646840989589691 0.650780975818634 0.654742002487183 0.658719003200531 0.662715017795563 0.666730999946594 0.670768022537231 0.674822986125946 0.678897023200989 0.682995021343231 0.687108993530273 0.691200017929077 0.695312023162842 0.699440002441406 0.703592002391815 0.707764029502869 0.711956024169922 0.716170012950897 0.720400989055634 0.724654972553253 0.72893100976944 0.733228027820587 0.737546980381012 0.741886973381042 0.746245980262756 0.750630021095276 0.755035996437073 0.759407997131348 0.763801991939545 0.768215000629425 0.772653996944427 0.777113974094391 0.781597971916199 0.786103010177612 0.790632009506226 0.795180976390839 0.799755990505219 0.804354012012482 0.808975994586945 0.813620984554291 0.818288028240204 0.822980999946594 0.82769900560379 0.832373976707458 0.837069988250732 0.841794013977051 0.846538007259369 0.851310014724731 0.856104016304016 0.8609259724617 0.865769028663635 0.870640993118286 0.875539004802704 0.880459010601044 0.885408997535706 0.89038097858429 0.895381987094879 0.900407016277313 0.905462026596069 0.910461008548737 0.915485024452209 0.920535981655121 0.92561399936676 0.930718004703522 0.935850024223328 0.941008985042572 0.946195006370544 0.951409995555878 0.956651985645294 0.961923003196716 0.967221975326538 0.972549974918365 0.977907001972198 0.983292996883392 0.988708019256592 0.994058012962341 0.999432027339935 1.00484001636505 1.01027297973633 1.01573896408081 1.02123498916626 1.02675700187683 1.03231394290924 1.03790104389191 1.04351496696472 1.04916405677795 1.05483996868134 1.06055104732513 1.06629502773285 1.07206594944 1.07787501811981 1.08359801769257 1.08935296535492 1.09514498710632 1.10096395015717 1.1068160533905 1.11269998550415 1.11862194538116 1.12457203865051 1.13055598735809 1.13657903671265 1.14263105392456 1.148717045784 1.15483796596527 1.16099894046783 1.16718995571136 1.17341697216034 1.17954695224762 1.18571603298187 1.1919150352478 1.19815003871918 1.20441997051239 1.21072602272034 1.21707403659821 1.22345304489136 1.22986805438995 1.23632597923279 1.24281597137451 1.2493439912796 1.25591099262238 1.26251995563507 1.26916396617889 1.27584600448608 1.2824170589447 1.2890260219574 1.2956680059433 1.30235397815704 1.30908000469208 1.31584405899048 1.32264304161072 1.32948696613312 1.33637201786041 1.34329199790955 1.35025894641876 1.35726702213287 1.36431801319122 1.37140500545502 1.3785400390625 1.38571798801422 1.39276194572449 1.39984798431396 1.40697598457336 1.41414105892181 1.42135405540466 1.42861104011536 1.43591105937958 1.44325494766235 1.45064401626587 1.45807695388794 1.46555602550507 1.47307395935059 1.48064494132996 1.48826205730438 1.495924949646 1.5036369562149 1.5111939907074 1.51879703998566 1.52644598484039 1.53413498401642 1.54187703132629 1.54966795444489 1.55750596523285 1.5653920173645 1.57332694530487 1.58131098747253 1.58934497833252 1.59742295742035 1.60555803775787 1.61374402046204 1.62198102474213 1.63027095794678 1.63838303089142 1.64653897285461 1.6547520160675 1.66301703453064 1.67133295536041 1.67969405651093 1.68811404705048 1.69658803939819 1.70511496067047 1.71369695663452 1.72232604026794 1.73101699352264 1.73976397514343 1.74856698513031 1.75741899013519 1.76633596420288 1.77504897117615 1.78381705284119 1.79263305664062 1.80151295661926 1.81044900417328 1.81944298744202 1.82849395275116 1.83759605884552 1.84676396846771 1.85599195957184 1.8652800321579 1.87462902069092 1.88403904438019 1.89350199699402 1.90303599834442 1.91263294219971 1.92199397087097 1.9314169883728 1.94089996814728 1.950443983078 1.96006000041962 1.9697300195694 1.97946298122406 1.98925995826721 1.99912095069885 2.0090479850769 2.01904106140137 2.02910995483398 2.0392370223999 2.04943108558655 2.05969500541687 2.07002711296082 2.08009791374207 2.09023499488831 2.10043907165527 2.11070990562439 2.12104988098145 2.13145899772644 2.1419370174408 2.1524760723114 2.16309595108032 2.17378711700439 2.18455100059509 2.19538807868958 2.2063000202179 2.21728610992432 2.22834706306458 2.23948502540588 2.25031995773315 2.26123809814453 2.27221894264221 2.28327393531799 2.29440402984619 2.30562090873718 2.31690406799316 2.3282630443573 2.33970093727112 2.3512179851532 2.3628249168396 2.37450194358826 2.38626098632812 2.39810109138489 2.41003608703613 2.42204403877258 2.43371200561523 2.44547200202942 2.45730090141296 2.46921110153198 2.48121500015259 2.49329090118408 2.505450963974 2.51769590377808 2.53003907203674 2.54245591163635 2.55496191978455 2.56756806373596 2.58025097846985 2.59302496910095 2.60590410232544 2.61886191368103 2.63143801689148 2.64410209655762 2.65685510635376 2.66969799995422 2.68264389038086 2.69566893577576 2.70878601074219 2.72199702262878 2.73530197143555 2.74870204925537 2.76219892501831 2.77579307556152 2.78950095176697 2.80329394340515 2.81718707084656 2.83118295669556 2.84474897384644 2.85841202735901 2.87217211723328 2.88603091239929 2.89998888969421 2.91404795646667 2.92820811271667 2.94247102737427 2.95683908462524 2.97131109237671 2.98588991165161 3.00057601928711 3.01537108421326 3.03027510643005 3.04529094696045 3.0604190826416 3.0750629901886 3.08981204032898 3.10466909408569 3.11963391304016 3.13470888137817 3.14989399909973 3.1651918888092 3.18060207366943 3.19610905647278 3.21175003051758 3.22750806808472 3.24338388442993 3.2593789100647 3.27549600601196 3.29173493385315 3.30809807777405 3.32391309738159 3.3398449420929 3.35589408874512 3.37206292152405 3.38835191726685 3.40476202964783 3.42129707336426 3.43797492980957 3.45476007461548 3.47167205810547 3.48871302604675 3.50588488578796 3.52318811416626 3.54062509536743 3.55819702148438 3.57590508460999 3.59299206733704 3.61022806167603 3.62757301330566 3.64507007598877 3.66267895698547 3.68044304847717 3.69832110404968 3.7163360118866 3.73451209068298 3.75280690193176 3.77126598358154 3.78984594345093 3.80859398841858 3.82746601104736 3.846510887146 3.86568307876587 3.88417196273804 3.90280199050903 3.9215989112854 3.94051599502563 3.95957899093628 3.97879004478455 3.99815106391907 4.01766204833984 4.03735208511353 4.05717086791992 4.07714700698853 4.09728193283081 4.11757612228394 4.13806009292603 4.15868091583252 4.17946910858154 4.19950485229492 4.21966981887817 4.24001884460449 4.26049995422363 4.28114223480225 4.30197381973267 4.32297086715698 4.34410715103149 4.36541080474854 4.38691520690918 4.40859222412109 4.43041515350342 4.45244407653809 4.47462177276611 4.49698209762573 4.51955604553223 4.54124116897583 4.56312799453735 4.58515691757202 4.60736083984375 4.62977504730225 4.65233612060547 4.67511081695557 4.69803810119629 4.72115087509155 4.74448490142822 4.76797723770142 4.79169607162476 4.81557607650757 4.83965492248535 4.86396789550781 4.88844919204712 4.91198205947876 4.93570423126221 4.95958185195923 4.98368883132935 5.00799083709717 5.0324912071228 5.0571551322937 5.08205890655518 5.10716915130615 5.13244915008545 5.15797710418701 5.18371820449829 5.20967578887939 5.23585176467896 5.26221084594727 5.28883409500122 5.31436681747437 5.34010887145996 5.36606121063232 5.39222717285156 5.41860914230347 5.44520902633667 5.4720311164856 5.49907684326172 5.5263500213623 5.55385303497314 5.58158922195435 5.60956001281738 5.63777017593384 5.6662220954895 5.69491910934448 5.72386407852173 5.75160217285156 5.77956914901733 5.80776882171631 5.83620500564575 5.86483192443848 5.89374685287476 5.92290687561035 5.9523138999939 5.98197317123413 6.01188516616821 6.04205513000488 6.07243585586548 6.10312986373901 6.13409185409546 6.16532611846924 6.19683313369751 6.22700023651123 6.25736808776855 6.28804492950439 6.31898307800293 6.35013103485107 6.38159894943237 6.41333818435669 6.44529581069946 6.4775857925415 6.5101580619812 6.54301404953003 6.57610177993774 6.60953903198242 6.6432728767395 6.67724895477295 6.71158695220947 6.74437379837036 6.77744102478027 6.81085205078125 6.84449100494385 6.87841987609863 6.91264390945435 6.94722986221313 6.98205614089966 7.01718902587891 7.05263185501099 7.08845520019531 7.12453317642212 7.16093397140503 7.19772911071777 7.23478984832764 7.27218818664551 7.3079252243042 7.34390211105347 7.38026285171509 7.41694211959839 7.45387077331543 7.49119901657104 7.52885723114014 7.5668511390686 7.60511112213135 7.64378881454468 7.68281698226929 7.72212219238281 7.76186084747314 7.80196523666382 7.84235811233521 7.88320398330688 7.922119140625 7.96137619018555 8.00106239318848 8.04101943969727 8.08133220672607 8.12200736999512 8.16304683685303 8.20445919036865 8.2463321685791 8.28850078582764 8.33105659484863 8.37400436401367 8.41734886169434 8.46118640899658 8.50534152984619 8.54991340637207 8.59241771697998 8.63521099090576 8.67848110198975 8.72204685211182 8.7661018371582 8.81046199798584 8.85532379150391 8.90059757232666 8.94619083404541 8.99230575561523 9.03874969482422 9.08572864532471 9.13304710388184 9.18091487884521 9.22913265228271 9.27791500091553 9.32428646087646 9.3710765838623 9.41829490661621 9.46594429016113 9.51391983032227 9.56245136260986 9.61143112182617 9.66086673736572 9.71076679229736 9.76113510131836 9.81198024749756 9.86330699920654 9.91500473022461 9.96731662750244 10.020133972168 10.0734605789185 10.1240892410278 10.1753072738647 10.2268733978271 10.2789163589478 10.3315725326538 10.384593963623 10.4381141662598 10.4922714233398 10.5468111038208 10.6018695831299 10.6574573516846 10.7137174606323 10.7703828811646 10.8276014328003 10.88551902771 10.9438638687134 10.9993276596069 11.0551643371582 11.1116695404053 11.1685590744019 11.2261352539062 11.284107208252 11.3427839279175 11.402024269104 11.4616804122925 11.5220699310303 11.5828895568848 11.6444625854492 11.7064800262451 11.7692718505859 11.8325233459473 11.8965702056885 11.9572277069092 12.0184602737427 12.0802736282349 12.1426773071289 12.2056789398193 12.269287109375 12.3335113525391 12.3983602523804 12.4638433456421 12.5299692153931 12.5967483520508 12.6641893386841 12.7323026657104 12.8010988235474 12.8705883026123 12.9407796859741 13.0071716308594 13.0743970870972 13.1420707702637 13.2103986740112 13.2795944213867 13.3492612838745 13.4196109771729"
		},
		 # Force output of trees to allow us to extract final time halo masses.
		 {
		     name  => "mergerTreeOutput",
		     value => "true"
		 },
		# Set cosmological parameters to match the Millennium Simulation.
		{
		    name  => "H_0",
		    value => 70.4
		},
		{
		    name  => "Omega_Matter",
		    value => 0.272
		},
		{
		    name  => "Omega_DE",
		    value => 0.728
		},
		{
		    name  => "Omega_b",
		    value => 0.0455
		},
		{
		    name  => "sigma_8",
		    value => 0.810
		},
		{
		    name  => "transferFunctionFile",
		    value => "data/largeScaleStructure/transferFunctionMillGas.xml"
		},

	 {
	     name => "mergerTreeComputeConditionalMassFunction",
	     value => "true"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionParentMassCount",
	     value => "10"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionParentMassMinimum",
	     value => "1.0e10"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionParentMassMaximum",
	     value => "1.0e15"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionMassRatioCount",
	     value => "30"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionMassRatioMinimum",
	     value => "1.0e-5"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionMassRatioMaximum",
	     value => "1.0e+1"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionParentRedshifts",
	     value => "0.0 0.0 0.0 0.0 0.0"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionProgenitorRedshifts",
	     value => "0.019933 0.508591 0.988708 2.070027 3.865683"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth",
	     value => "4"
	 },


		);
	}
    }
    # Push the model to the stack.
    my %monteCarloModel =
	(
	 label      => "monteCarlo".$subVolumeSuffix,
	 dependsOn  => $modelsDirectory."nBody".$subVolumeSuffix."/treeMasses.hdf5",
	 parameters => \@monteCarloParameters
	);
    push(@models,\%monteCarloModel);
}
for(my $iSubvolume=0;$iSubvolume<$subVolumeCount;++$iSubvolume) {
    my $subVolumeSuffix = "";
    $subVolumeSuffix = $iSubvolume
	if ( $haveSubvolumes == 1 );
    my @nBodyParameters =
	(
		 # Switch to using N-body trees.
		 {
		  name  => "mergerTreeConstructMethod",
		  value => "read"
		 },
		 # Force output of trees to allow us to extract final time halo masses.
		 {
		     name  => "mergerTreeOutput",
		     value => "true"
		 },
		 # Prune branches of the tree to keep only those with at least 20 particles.
		 {
		     name  => "mergerTreePruneBranches",
		     value => "true"
		 },
	         # Set mass function analysis to use Poisson model appropriate for N-body mass function sampling.
	         {
	             name  => "analysisMassFunctionCovarianceModel",
	             value => "Poisson"
	         },
		 # Set merger tree reading options.
		 {
		     name  => "mergerTreeReadTreeIndexToRootNodeIndex",
		     value => "true"
		 },
		 {
		     name  => "allTreesExistAtFinalTime",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadOutputTimeSnapTolerance",
		     value => 1.0e-3
		 },
		 {
		     name  => "mergerTreeReadPresetPositions",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetOrbits",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetMergerTimes",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetMergerNodes",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetSubhaloMasses",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetSpins",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetScaleRadii",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadAllowBranchJumps",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadAllowSubhaloPromotions",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeEnforceMonotonicGrowth",
		     value => "false"
		 },
		 # {
		 #     name  => "accretionHalosSimpleAccreteNewGrowthOnly",
		 #     value => "true"
		 # },
		 # {
		 #     name  => "treeNodeMethodBasic",
		 #     value => "standardTracking"
		 # },
		 {
		     name  => "powerSpectrumIndex",
		     value => 1.0
		 },
		 {
		     name  => "powerSpectrumReferenceWavenumber",
		     value => 1.0
		 },
		 {
		     name  => "powerSpectrumRunning",
		     value => 0.0
		 },
		 {
		     name  => "transferFunctionMethod",
		     value => "file"
		 },
		 # Set an initial random number seed.
		 {
		     name  => "randomSeed",
		     value => $randomSeed
		 },

	 {
	     name => "mergerTreeComputeConditionalMassFunction",
	     value => "true"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionParentMassCount",
	     value => "10"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionParentMassMinimum",
	     value => "1.0e10"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionParentMassMaximum",
	     value => "1.0e15"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionMassRatioCount",
	     value => "30"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionMassRatioMinimum",
	     value => "1.0e-5"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionMassRatioMaximum",
	     value => "1.0e+1"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionParentRedshifts",
	     value => "0.0 0.0 0.0 0.0 0.0"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionProgenitorRedshifts",
	     value => "0.019933 0.508591 0.988708 2.070027 3.865683"
	 },
	 {
	     name => "mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth",
	     value => "4"
	 },

	);
    # Add tree specific settings.
    switch ( $arguments{'trees'} ) {
	case ( 'millennium:MPA' ) {
	    push(
		@nBodyParameters,
		 {
		     name  => "mergerTreeReadFileName",
		     value => $treeDirectory."/subvolume".$iSubvolume.".hdf5"
		 },
		 # Set cosmological parameters.
		 {
		     name  => "H_0",
		     value => 73.0
		 },
		 {
		     name  => "Omega_Matter",
		     value => 0.25
		 },
		 {
		     name  => "Omega_DE",
		     value => 0.75
		 },
		 {
		     name  => "Omega_b",
		     value => 0.0455
		 },
		 {
		     name  => "sigma_8",
		     value => 0.9
		 },
		 {
		     name  => "transferFunctionFile",
		     value => "data/largeScaleStructure/transferFunctionMillenniumSimulation.xml"
		 },
		 {
		     name  => "mergerTreePruningMassThreshold",
		     value => 2.3578789e10
		 }
		);
	}
	case ( 'MillGas' ) {
	    push(
		@nBodyParameters,
		 {
		     name  => "mergerTreeReadFileName",
		     value => $treeDirectory."/tree_976.glc.hdf5"
		 },
		 # Set cosmological parameters.
		 {
		     name  => "H_0",
		     value => 70.4
		 },
		 {
		     name  => "Omega_Matter",
		     value => 0.272
		 },
		 {
		     name  => "Omega_DE",
		     value => 0.728
		 },
		 {
		     name  => "Omega_b",
		     value => 0.0455
		 },
		 {
		     name  => "sigma_8",
		     value => 0.810
		 },
		 {
		     name  => "transferFunctionFile",
		     value => "data/largeScaleStructure/transferFunctionMillGas.xml"
		 },
		 {
		     name  => "mergerTreePruningMassThreshold",
		     value => 2.6602117120e10
		 }
		);
	}
    }
    push(
	@models,
	{
	    label      => "nBody".$subVolumeSuffix,
	    parameters => \@nBodyParameters
	}
	);
}

# Iterate over models.
foreach my $model ( @models ) {
    # Specify the output name.
    (my $safeLabel      = $model->{'label'}) =~ s/:/_/g;
    my $modelDirectory = $modelsDirectory.$safeLabel;
    system("mkdir -p ".$modelDirectory);
    my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
    push(
	@{$model->{'parameters'}},
	{
	    name  => "galacticusOutputFileName",
	    value => $galacticusFileName
	}
	);
    my $newParameters = clone($parameters);
    $newParameters->{'parameter'}->{$_->{'name'}}->{'value'} = $_->{'value'}
    foreach ( @{$model->{'parameters'}} );
    # Run the model.
    unless ( -e $galacticusFileName ) {
	# Generate the parameter file.
	my $parameterFileName = $modelDirectory."/parameters.xml";
	&Parameters::Output($newParameters,$parameterFileName);
	# Create a batch script for PBS.
	my $batchScriptFileName = $modelDirectory."/launch.pbs";
	open(oHndl,">".$batchScriptFileName);
	print oHndl "#!/bin/bash\n";
	print oHndl "#PBS -N monteCarloTrees".$model->{'label'}."\n";
	print oHndl "#PBS -l walltime=3:00:00\n";
	print oHndl "#PBS -l mem=4gb\n";
	print oHndl "#PBS -l nodes=1:ppn=12\n";
	print oHndl "#PBS -j oe\n";
	print oHndl "#PBS -o ".$modelDirectory."/launch.log\n";
	print oHndl "#PBS -V\n";
	print oHndl "cd \$PBS_O_WORKDIR\n";
	print oHndl "export LD_LIBRARY_PATH=/home/abenson/Galacticus/Tools/lib:/home/abenson/Galacticus/Tools/lib64:\$LD_LIBRARY_PATH\n";
	print oHndl "export PATH=/home/abenson/Galacticus/Tools/bin:\$PATH\n";
	print oHndl "export GFORTRAN_ERROR_DUMPCORE=YES\n";
	print oHndl "ulimit -t unlimited\n";
	print oHndl "ulimit -c unlimited\n";
	print oHndl "export OMP_NUM_THREADS=12\n";
	print oHndl "mpirun --bynode -np 1 Galacticus.exe ".$parameterFileName."\n";
	print oHndl $galacticusPath."/constraints/modelDiscrepancy/monteCarloTreesExtractHaloMasses.pl ".$modelDirectory."\n"
	    if ( $model->{'label'} =~ m/^nBody/ );
	foreach my $constraint ( @constraints ) {
	    # Parse the definition file.
	    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
	    # Insert code to run the analysis code.
	    my $analysisCode = $constraintDefinition->{'analysis'};
	    print oHndl $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".xml\n";
	}
	close(oHndl);
	# Queue the calculation.
	my %pbsData = 
	    (
	     script    => $batchScriptFileName,
	     model     => $modelDirectory,
	     label     => $model->{'label'}
	    );
	$pbsData{'dependsOn'} = $model->{'dependsOn'}
	   if ( exists( $model->{'dependsOn'}) );
	push(@pbsStack,\%pbsData);   
    }
}
# Send jobs to PBS.
&PBS_Submit(@pbsStack)
    if ( scalar(@pbsStack) > 0 );

# Iterate over constraints.
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Iterate over Monte Carlo models.
    my $monteCarloMass;
    my $monteCarloMassFunction;
    my $monteCarloCovariance;
    for(my $iSubvolume=0;$iSubvolume<$subVolumeCount;++$iSubvolume) {
	# Construct the subvolume suffix.
	my $subVolumeSuffix = "";
	$subVolumeSuffix = $iSubvolume
	    if ( $haveSubvolumes == 1 );
	# Locate the Monte Carlo model results.
	my $monteCarloResultFileName   = $modelsDirectory."monteCarlo".$subVolumeSuffix."/".$constraintDefinition->{'label'}.".xml";
	# Read the results.
	my $monteCarloResult           = $xml->XMLin($monteCarloResultFileName);
	# Extract the results.
	my $thisMonteCarloMass         = pdl @{$monteCarloResult->{'x'         }};
	my $thisMonteCarloMassFunction = pdl @{$monteCarloResult->{'y'         }};
	my $thisMonteCarloCovariance   = pdl @{$monteCarloResult->{'covariance'}};
	my $ySize                      = nelem($thisMonteCarloMassFunction);
	$thisMonteCarloCovariance      = reshape($thisMonteCarloCovariance,$ySize,$ySize);
	unless ( defined($monteCarloMass) ) {
	    $monteCarloMass         = $thisMonteCarloMass;
	    $monteCarloMassFunction = pdl zeroes($ySize       );
	    $monteCarloCovariance   = pdl zeroes($ySize,$ySize);
	}
	$monteCarloMassFunction += $thisMonteCarloMassFunction;
	$monteCarloCovariance   += $thisMonteCarloCovariance  ;
    }
    if ( $haveSubvolumes == 1 ) {
	$monteCarloMassFunction *=  $subvolumesMaximum/$subVolumeCount    ;
	$monteCarloCovariance   *= ($subvolumesMaximum/$subVolumeCount)**2;
    }
    # Iterate over N-body models.
    my $nBodyMass;
    my @nBodyMassFunctions;
    my @nBodyMassFunctionCovariances;
    for(my $iSubvolume=0;$iSubvolume<$subVolumeCount;++$iSubvolume) {
	# Construct the subvolume suffix.
	my $subVolumeSuffix = "";
	$subVolumeSuffix = $iSubvolume
	    if ( $haveSubvolumes == 1 );
	# Locate the N-body model results.
	my $nbodyResultFileName = $modelsDirectory."nBody".$subVolumeSuffix."/".$constraintDefinition->{'label'}.".xml";
	# Read the results.
	my $nbodyResult         = $xml->XMLin($nbodyResultFileName);
	# Extract the results.
	$nBodyMass                        = pdl @{$nbodyResult->{'x'         }};
	my $thisMassFunction              = pdl @{$nbodyResult->{'y'         }};
	my $thisCovariance                = pdl @{$nbodyResult->{'covariance'}};
	my $ySize                         = nelem($thisMassFunction);
	$thisCovariance                   = reshape($thisCovariance,$ySize,$ySize);
	push(@nBodyMassFunctions          ,$thisMassFunction);
	push(@nBodyMassFunctionCovariances,$thisCovariance);
    }
    # Construct samples
    my @covariances;
    my $combineCount = 1;
    my $nBodyMassFunction;
    if ( $haveSubvolumes == 1 ) {
	while ( $combineCount <= $subVolumeCount ) {
	    # Combine mass functions.
	    my @combinedMassFunctions;
	    my $subVolumeBegin = 0;
	    while ( $subVolumeBegin+$combineCount-1 < $subVolumeCount ) {
		my $massFunction = pdl zeroes(nelem($nBodyMass));
		for(my $subVolume=$subVolumeBegin;$subVolume<$subVolumeBegin+$combineCount;++$subVolume) {
		    $massFunction += $nBodyMassFunctions[$subVolume];
		}
		$massFunction      *= $subvolumesMaximum/$combineCount;
		push(@combinedMassFunctions,$massFunction);
		$nBodyMassFunction  = $massFunction
		    if ( $combineCount == $subVolumeCount );
		$subVolumeBegin    += $combineCount
	    }
	    # Measure the covariance.
	    if ( $combineCount < $subVolumeCount ) {
		my $covariance = pdl zeroes(nelem($nBodyMass),nelem($nBodyMass));
		for(my $i=0;$i<nelem($nBodyMass);++$i) {
		    for(my $j=0;$j<nelem($nBodyMass);++$j) {
			my $iMean = pdl 0.0;
			my $jMean = pdl 0.0;
			foreach my $combinedMassFunction ( @combinedMassFunctions ) {
			    $iMean                   += $combinedMassFunction->($i)                            ;
			    $jMean                   +=                             $combinedMassFunction->($j);
			    $covariance->(($i),($j)) += $combinedMassFunction->($i)*$combinedMassFunction->($j);
			}
			$covariance->(($i),($j)) .= ($covariance->(($i),($j))-$iMean*$jMean/scalar(@combinedMassFunctions))/(scalar(@combinedMassFunctions)-1);
		    }
		}
		push(@covariances,$covariance);
	    }
	    # Double the number of subvolumes to combined.
	    $combineCount *= 2;
	}
    } else {
	$nBodyMassFunction = $nBodyMassFunctions[0];	
    }
    # Find the multiplicative between these two models.
    (my $nonZero, my $zero)                      = which_both($monteCarloMassFunction > 0.0);
    my $modelDiscrepancyMultiplicative           = $nBodyMassFunction->copy();
    $modelDiscrepancyMultiplicative->($nonZero) /= $monteCarloMassFunction->($nonZero);
    $modelDiscrepancyMultiplicative->($zero   ) .= 1.0;
    # Find the covariance remaining in the complete set of subvolumes.
    my $modelDiscrepancyCovariance;
    if ( $haveSubvolumes == 1 ) {
	my $extrapolateFrom                       = scalar(@covariances)-3;
	$extrapolateFrom                          = 0
	    if ( $extrapolateFrom < 0 );
	$modelDiscrepancyCovariance               = $covariances[$extrapolateFrom]/(2**(scalar(@covariances)-$extrapolateFrom));
    } else {
	$modelDiscrepancyCovariance               = $nBodyMassFunctionCovariances[0];
    }
    # Compute the covariance.
    my $modelDiscrepancyCovarianceMultiplicative = 
	 $monteCarloCovariance      *outer($nBodyMassFunction/$monteCarloMassFunction**2,$nBodyMassFunction/$monteCarloMassFunction**2)
	+$modelDiscrepancyCovariance*outer(               1.0/$monteCarloMassFunction   ,               1.0/$monteCarloMassFunction   );
    # Output the model discrepancy to file.
    my $outputFile = new PDL::IO::HDF5(">".$modelsDirectory."discrepancy".ucfirst($constraintDefinition->{'label'}).".hdf5");
    $outputFile->dataset('multiplicative'          )->set($modelDiscrepancyMultiplicative          );
    $outputFile->dataset('multiplicativeCovariance')->set($modelDiscrepancyCovarianceMultiplicative);
    $outputFile->attrSet(
     	description => "Model discrepancy for ".$constraintDefinition->{'name'}." due to use of Monte Carlo merger trees."
     	);
}

exit;

sub PBS_Submit {
    # Submit jobs to PBS and wait for them to finish.
    my @pbsStack = @_;
    my %pbsJobs;
    # Determine maximum number allowed in queue at once.
    my $jobMaximum = 10;
    $jobMaximum = $arguments{'pbsJobMaximum'}
    if ( exists($arguments{'pbsJobMaximum'}) );
    # Submit jobs and wait.
    print "Waiting for PBS jobs to finish...\n";
    while ( scalar(keys %pbsJobs) > 0 || scalar(@pbsStack) > 0 ) {
	# Find all PBS jobs that are running.
	my %runningPBSJobs;
	undef(%runningPBSJobs);
	open(pHndl,"qstat -f|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^Job\sId:\s+(\S+)/ ) {$runningPBSJobs{$1} = 1};
	}
	close(pHndl);
	foreach my $jobID ( keys(%pbsJobs) ) {
	    unless ( exists($runningPBSJobs{$jobID}) ) {
		print "PBS job ".$jobID." has finished.\n";
		# Remove the job ID from the list of active PBS jobs.
		delete($pbsJobs{$jobID});
	    }
	}
	# If fewer than maximum number of jobs are in the queue, pop one off the stack.
	while ( scalar(@pbsStack) > 0 && scalar(keys %pbsJobs) < $jobMaximum ) {
	    my $pbsData     = pop(@pbsStack);
	    my $batchScript = $pbsData->{'script'};
	    if ( exists($pbsData->{'dependsOn'}) && ! -e $pbsData->{'dependsOn'} ) {
		unshift(@pbsStack,$pbsData);
		sleep 5;
		last;
	    } else {
		# Submit the PBS job.
		open(pHndl,"qsub ".$batchScript."|");
		my $jobID = "";
		while ( my $line = <pHndl> ) {
		    if ( $line =~ m/^(\d+\S+)/ ) {$jobID = $1};
		}
		close(pHndl);	    
		# Add the job number to the active job hash.
		unless ( $jobID eq "" ) {
		    $pbsJobs{$jobID} = $pbsData;
		}
	    }
	}
	sleep 5;
    }
}
