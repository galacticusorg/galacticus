# Provide tools for validation purposes.

package Galacticus::Validation;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use JSON::PP;
use Git;

sub extract {
    my $fileName = shift();
    my $name     = shift();
    my $suffix   = shift();
    my @likelihoods;
    my @results;
    my $model    = new PDL::IO::HDF5($fileName);
    my $analyses = $model->group('analyses');
    foreach my $analysisName ( $analyses->groups() ) {
	my $analysisGroup = $analyses->group($analysisName);
	my @attributeNames = $analysisGroup->attrs();
	my $attributes =
	{
	    name        => $analysisName,
	    xAxisLabel  => "x",
	    yAxisLabel  => "y",
	    xAxisIsLog  =>  0 ,
	    yAxisIsLog  =>  1 ,
	    targetLabel => "data"
	};
	($attributes->{$_}) = $analysisGroup->attrGet($_)
	    foreach ( @attributeNames );
	(my $logLikelihood) = $analysisGroup->attrGet('logLikelihood');
	print $analysisName."\t".$logLikelihood."\n";
	push(
	    @likelihoods,
	    {
		name  => $name." - Likelihood - ".$analysisName,
		unit  => "-logℒ"                                                 ,
		value => abs($logLikelihood->sclr())
	    }
	    );
	# Skip cases for which we have no "type" specified.
	unless ( exists($attributes->{'type'}) ) {
	    print "Warning: analysis '".$analysisName."' has no 'type' attribute, so it can not be processed.\n";
	} elsif ( $attributes->{'type'} eq "function1D" ) {
	    # Simple 1D function - will be shown as an x-y scatter plot.
	    # Validate attributes.
	    my @datasetNames =
		(
		 {name => 'xDataset'         , required => 1},
		 {name => 'yDataset'         , required => 1},
		 {name => 'yDatasetTarget'   , required => 0},
		 {name => 'yCovariance'      , required => 0},
		 {name => 'yCovarianceTarget', required => 0},
		 {name => 'yErrorLower'      , required => 0},
		 {name => 'yErrorUpper'      , required => 0},
		 {name => 'yErrorLowerTarget', required => 0},
		 {name => 'yErrorUpperTarget', required => 0},
		);
	    foreach my $attribute ( @datasetNames ) {
		die("Error: attribute '".$attribute->{'name'}."' is missing from analysis '".$analysisName."' but is required.")
		    unless ( ! $attribute->{'required'} || grep {$_ eq $attribute->{'name'}} keys(%{$attributes}) );
	    }
	    # Read the datasets.
	    my $data;
	    foreach my $dataset ( @datasetNames ) {
		if ( grep {$_ eq $dataset->{'name'}} keys(%{$attributes}) ) {
		    (my $analysisDatasetName) = $analysisGroup->attrGet($dataset->{'name'});
		    unless ( grep {$_ eq $analysisDatasetName} $analysisGroup->datasets() ) {
			print "Analysis: ".$analysisName."\n";
			print "Generic name: ".$dataset->{'name'}."\n";
			print "Actual name: ".$analysisDatasetName."\n";
			print "Available datasets:\n";
			print join("\n",map {"\t".$_} $analysisGroup->datasets())."\n";
			die("failed to find dataset");
		    }
		    $data->{$dataset->{'name'}} = $analysisGroup->dataset($analysisDatasetName)->get();
		    unless ( defined($data->{$dataset->{'name'}}) ) {
			print "Analysis: ".$analysisName."\n";
			print "Generic name: ".$dataset->{'name'}."\n";
			print "Actual name: ".$analysisDatasetName."\n";
			print "Available datasets:\n";
			print join("\n",map {"\t".$_} $analysisGroup->datasets())."\n";
			system("h5dump -A -g analyses/".$analysisName." ".$fileName);
			die("failed to read dataset");
		    }
		}
	    }
	    # Extract errors.
	    if ( exists($data->{'yCovariance'      }) ) {
		$data->{'yError'      } = $data->{'yCovariance'      }->diagonal(0,1)->sqrt();
		delete($data->{'yCovariance'      });
	    }
	    if ( exists($data->{'yCovarianceTarget'}) ) {
		$data->{'yErrorTarget'} = $data->{'yCovarianceTarget'}->diagonal(0,1)->sqrt();
		delete($data->{'yCovarianceTarget'});
	    }
	    # Store results.
	    my $result;
	    foreach my $attributeName ( keys(%{$attributes}) ) {
		if ( $attributeName =~ m/^[xy]AxisIsLog$/ || $attributeName eq "logLikelihood" ) {
		    $result->{'attributes'}->{$attributeName} = $attributes->{$attributeName}->sclr();
		} else {
		    # LaTeX conversions.
		    $attributes->{$attributeName} =~ s/\$//g;
		    $attributes->{$attributeName} =~ s/\\mathrm\{([^\}]+)\}/$1/g;
		    $attributes->{$attributeName} =~ s/\\hbox\{([^\}]+)\}/$1/g;
		    $attributes->{$attributeName} =~ s/\\odot/☉/g;
		    $attributes->{$attributeName} =~ s/\\langle/⟨/g;
		    $attributes->{$attributeName} =~ s/\\rangle/⟩/g;
		    $result->{'attributes'}->{$attributeName} = $attributes->{$attributeName};
		}
	    }
	    foreach my $dataName ( keys(%{$data}) ) {
		@{$result->{'data'}->{$dataName}} = $data->{$dataName}->list();
	    }
	    push(
		@results,
		$result
		);
	}
    }
    {
	# Write benchmark results.
	my $json = JSON::PP->new()->pretty()->encode(\@likelihoods);
	open(my $reportFile,">","outputs/validate_".$suffix.".json");
	print $reportFile $json;
	close($reportFile);
    }
    {
	# Interface with git.
	my $repo         = Git->repository(Directory => $ENV{'GALACTICUS_EXEC_PATH'});
	my $lastRevision = $repo->command_oneline( [ 'rev-list', '--all' ], STDERR => 0 );
	(my $authorName  = $repo->command_oneline( [ 'show', '-s', '--format="%an"', $lastRevision ], STDERR => 0 )) =~ s/"//g;
	(my $authorEmail = $repo->command_oneline( [ 'show', '-s', '--format="%ae"', $lastRevision ], STDERR => 0 )) =~ s/"//g;
	(my $authorDate  = $repo->command_oneline( [ 'show', '-s', '--format="%aD"', $lastRevision ], STDERR => 0 )) =~ s/"//g;
	(my $message     = $repo->command_oneline( [ 'show', '-s', '--format="%s"' , $lastRevision ], STDERR => 0 )) =~ s/"//g;

	# Write results.
	my $output;
	$output =
	{
	    repoUrl => "https://github.com/galacticusorg/galacticus",
	    commit =>
	    {
		author =>
		{
		    name => $authorName,
		    email => $authorEmail
		},
		id        => $lastRevision,
		message   => $message,
		timestamp => $authorDate,
		url       => "https://github.com/galacticusorg/galacticus/commit/".$lastRevision
	    },
	    results => \@results
	};
	my $json = JSON::PP->new()->pretty()->encode($output);
	open(my $reportFile,">","outputs/results_".$suffix.".json");
	print $reportFile "window.ANALYSES_DATA = ";
	print $reportFile $json;
	close($reportFile);
    }
}

1;
