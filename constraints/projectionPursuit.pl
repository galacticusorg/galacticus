#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::MatrixOps;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
use Data::Dumper;
use Galacticus::Options;
use Galacticus::Launch::PBS;
use Galacticus::Constraints::Parameters;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Compute projections of the posterior which give the strongest constraints.
# Andrew Benson (27-May-2016)

# Get arguments.
die("Usage: projectionPursuit.pl <configFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Process options.
my %options = 
    (
     make                => "yes",
     sampleFrom          =>    -1,
     sampleCount         =>    -1,
     chain               => "all",
     outliers            =>    "",
     contributionMinimum =>  0.05,
     eigenVectorsRetain  =>     5,
     tableFile           => "projection.tex",
     modelsDirectory     => "projectionModels"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::Parse_Config($configFile);

# Validate the config file.
die("projectionPursuit.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'likelihood'}->{'workDirectory' }) );

# Find maximum likelihood model parameters.
print "Finding maximum likelihood model...\n";
my $parametersMaximumLikelihood = &Galacticus::Constraints::Parameters::Maximum_Likelihood_Vector($config,\%options);

# Generate a sample of model parameters.
print "Building matrix of samples from model posterior...\n";
my $parametersMatrix            = &Galacticus::Constraints::Parameters::Sample_Matrix            ($config,\%options);
print "...found ".$parametersMatrix->dim(1)." samples\n";

# For each parameter, apply mappings and normalize sampled states by the variance in the prior.
print "Mapping and normalizing...\n";
my @parameters;
my $j = -1;
my $parameterVariances = pdl [];
for(my $i=0;$i<scalar(@{$config->{'parameters'}->{'parameter'}});++$i) {
    # Look for active parameters.
    if ( exists($config->{'parameters'}->{'parameter'}->[$i]->{'prior'}) ) {
	# Get the parameter.
	++$j;
	my $parameter = $config->{'parameters'}->{'parameter'}->[$i];
	# Apply any mapping.
	if ( exists($parameter->{'mapping'}) ) {
	    if ( $parameter->{'mapping'}->{'type'} eq "linear" ) {
		# Nothing to do.
	    } elsif ( $parameter->{'mapping'}->{'type'} eq "logarithmic" ) {
		# Apply logarithmic mapping.
		$parametersMatrix->(($i),:) .= log($parametersMatrix->(($i),:));
		$parametersMaximumLikelihood->(($i)) .= log($parametersMaximumLikelihood->(($i)));
		foreach my $limit ( "minimum", "maximum" ) {
		    $parameter->{'prior'}->{'distribution'}->{$limit} = log($parameter->{'prior'}->{'distribution'}->{$limit})
			if ( exists($parameter->{'prior'}->{'distribution'}->{$limit}) );
		}		
	    } else {
		die("Unknown mapping '".$parameter->{'mapping'}->{'type'}."'");
	    }
	}
	# Compute the variance of the prior and normalize by the root-variance.
	my $variance;
	if ( $parameter->{'prior'}->{'distribution'}->{'type'} eq "uniform" ) {
	    $variance = 
		+(
		  +$parameter->{'prior'}->{'distribution'}->{'maximum'}
		  -$parameter->{'prior'}->{'distribution'}->{'minimum'}
		 )**2
		/12.0;
	} elsif ( $parameter->{'prior'}->{'distribution'}->{'type'} eq "normal" ) {
	    my $limitLower = exists($parameter->{'prior'}->{'distribution'}->{'minimum'}) ? ($parameter->{'prior'}->{'distribution'}->{'minimum'}-$parameter->{'prior'}->{'distribution'}->{'mean'})/sqrt($parameter->{'prior'}->{'distribution'}->{'variance'}) : -1.0e6;
	    my $limitUpper = exists($parameter->{'prior'}->{'distribution'}->{'maximum'}) ?( $parameter->{'prior'}->{'distribution'}->{'maximum'}-$parameter->{'prior'}->{'distribution'}->{'mean'})/sqrt($parameter->{'prior'}->{'distribution'}->{'variance'}) : +1.0e6;
	    $variance = $parameter->{'prior'}->{'distribution'}->{'variance'}*(($limitLower*exp(-0.5*$limitLower**2)-$limitUpper*exp(-0.5*$limitUpper**2))/sqrt(2.0*PI)-0.5*(erf($limitLower/sqrt(2.0))-erf($limitUpper/sqrt(2.0))));
	} else {
	    die("Unknown prior '".$parameter->{'prior'}->{'distribution'}->{'type'}."'");
	}
	$parametersMatrix->(($j),:) /= sqrt($variance);
	$parameterVariances = $parameterVariances->append($variance);
	# Subtract the mean.
	my $mean = $parametersMatrix->(($j),:)->avg();
	$parametersMatrix->(($j),:) -= $mean;
	# Push the parameter onto the active list.
	push(@parameters,$parameter);
    }
}

# Compute eigenvectors of the sample matrix.
print "Finding eigenvectors...\n";
my $covarianceMatrix = transpose($parametersMatrix) x $parametersMatrix;
$covarianceMatrix /= $parametersMatrix->dim(1)-1;
(my $eigenVectors, my $eigenValues) = eigens_sym($covarianceMatrix);

# Generate index into sorted eigenvalues.
my $eigenValuesRank = $eigenValues->qsorti();

# Iterate over eigenvalues from smallest to largest. Generate a LaTeX table showing the dominant contributions to the lowest
# variance eigenvectors.
my $componentsMaximum = 0;
my @tableBody;
for(my $i=0;$i<nelem($eigenValues) && $i<$options{'eigenVectorsRetain'};++$i) {
    # Begin LaTeX table entry for this eigenvector.
    my $tableLine = "\$P_".$i."\$";
    # Generate an index into ranked contributions to eigenvector magnitude.
    my $j                = $eigenValuesRank->(($i)  );
    my $vector           = $eigenVectors   ->(($j),:);
    my $vectorSquared    = $vector**2;
    my $vectorMagnitude  = $vectorSquared->sum();
    $vectorSquared      /= $vectorMagnitude;
    my $parameterRank    = $vectorSquared->qsorti();   
    my $componentsCount  = 0; 
    for(my $k=nelem($parameterRank)-1;$k>=0;--$k) {
	my $l           = $parameterRank->(($k));
	my $label       = exists($parameters[$l]->{'label'}) ? "\$".$parameters[$l]->{'label'}."\$" : "{\tt ".$parameters[$l]->{'name'}."}";
	my $coefficient = sprintf("%4.2f",$vector->(($l)));
	$coefficient    =~ s/\./&.&/;
	$coefficient    = "+".$coefficient
	    if ( $vector->(($l)) > 0.0 );
	if ( $vectorSquared->(($l)) > $options{'contributionMinimum'} ) {
	    $tableLine .= " & ".$coefficient." & ".$label;
	    ++$componentsCount;
	}
    }
    # Finsih LaTeX table entry for this eigenvector.
    push(
	@tableBody,
	{
	    componentsCount => $componentsCount,
	    content         => $tableLine
	}
	);
    $componentsMaximum = $componentsCount
	if ( $componentsCount > $componentsMaximum );
}
# Construct the table.
open(my $tableFile,">".$options{'tableFile'});
print $tableFile "\\begin{tabular}{l".("r\@{}c\@{}l\@{ }l" x $componentsMaximum)."} \\\\\n";
print $tableFile "\\hline\n";
print $tableFile "{\\bf Vector} & \\multicolumn{".(4*$componentsMaximum)."}{c}{\\bf Components} \\\\\n";
print $tableFile "\\hline\n";
foreach ( @tableBody ) {
    print $tableFile $_->{'content'};
    if ( $_->{'componentsCount'} < $componentsMaximum ) {
	my $padding = " &&&" x ($componentsMaximum-$_->{'componentsCount'});
	print $tableFile $padding;
    }
    print $tableFile " \\\\\n";
}
print $tableFile "\\hline\n";
print $tableFile "\\end{tabular}\n";
close($tableFile);

# Iterate over eigenvalues from smallest to largest. Also generate parameter sets for models perturb the maximum likelihood model
# along these eigenvalues.
my $xml = new XML::Simple;
my @pbsJobs;
for(my $i=0;$i<nelem($eigenValues) && $i<$options{'eigenVectorsRetain'};++$i) {
    my $j = $eigenValuesRank->(($i));
    # Iterate over number of "sigma" perturbation to the model.
    for(my $n=-3;$n<=+3;++$n) {
	# Construct the perturbed parameters.
	my $perturbedModel = 
	    +$parametersMaximumLikelihood
	    +$n
	    *sqrt($eigenValues ->(($j)  ))
	    *     $eigenVectors->(($j),:) 
	    *sqrt($parameterVariances    );
	# Unmap the parameters.
	for(my $k=0;$k<scalar(@parameters);++$k) {
	    if ( exists($parameters[$k]->{'mapping'}) ) {
		if ( $parameters[$k]->{'mapping'}->{'type'} eq "linear" ) {
		    # Nothing to do.
		} elsif ( $parameters[$k]->{'mapping'}->{'type'} eq "logarithmic" ) {
		    # Remove logarithmic mapping.
		    $perturbedModel->(($k)) .= exp($perturbedModel->(($k)));
		} else {
		    die("Unknown mapping '".$parameters[$k]->{'mapping'}->{'type'}."'");
		}
	    }
	}
	# Generate model parameters, and a job to run the model.
	my $modelDirectory = $config->{'likelihood'}->{'workDirectory'}."/".$options{'modelsDirectory'}."/vector".$i."_perturb".$n."/";
	unless ( -e $modelDirectory."galacticus.hdf5" ) {
	    system("mkdir -p ".$modelDirectory);
	    my $perturbedParameters = &Galacticus::Constraints::Parameters::Convert_Parameter_Vector_To_Galacticus($config,$perturbedModel);
	    &Galacticus::Constraints::Parameters::Apply_Command_Line_Parameters($perturbedParameters,\%options);
	    $perturbedParameters->{'galacticusOutputFileName'}->{'value'} = $modelDirectory."galacticus.hdf5";
	    open(my $parameterFile,">".$modelDirectory."parameters.xml");
	    print $parameterFile $xml->XMLout($perturbedParameters, rootName => "parameters");
	    close($parameterFile);
	    my $command = "mpirun --bynode -np 1 Galacticus.exe ".$modelDirectory."parameters.xml\n";
	    my %job =
		(
		 launchFile => $modelDirectory."/launch.pbs",
		 label      => "projection_".$i."_".$n,
		 logFile    => $modelDirectory."/launch.log",
		 command    => $command
		);
	    foreach ( 'ppn', 'walltime', 'memory' ) {
		$job{$_} = $options{$_}
		if ( exists($options{$_}) );
	    }
	    # Queue the calculation.
	    push(@pbsJobs,\%job);
	}
     }
}
# Send jobs to PBS.
&Galacticus::Launch::PBS::SubmitJobs(\%options,@pbsJobs);

exit;
