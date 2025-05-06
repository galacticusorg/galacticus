# Contains a Perl module which implements various useful functions for working with simulations used as constraints.

package Galacticus::Constraints::Simulations;
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use PDL;
use Storable qw(dclone);
use Exporter 'import';
use Data::Dumper;
our @EXPORT_OK = qw(iterate selectSimulations matchSelection);

sub iterate {
    # Construct a list that can be iterated over to visit each simulation of interest.
    # The list can be truncated after any stage (simulation, group, etc.).
    # Various useful properties (e.g. Hubble parameter, particle mass, are automatically added to the data structure as needed).
    my $simulations    =   shift() ;
    my %options        = %{shift()};
    my (%optionsExtra) = @_
	if ( $#_ >= 1 );
    # Validate the 'stopAfter' option.
    $optionsExtra{'stopAfter'} = "redshift"
	unless ( exists($optionsExtra{'stopAfter'}) );
    my @allowedStopAfter = ( 'suite', 'group', 'resolution', 'simulation', 'realization', 'redshift' );
    die('unrecognized "stopAfter" option - allowed choices are: '.join(", ",map {'"'.$_.'"'} @allowedStopAfter))
	unless ( grep {$_ eq $optionsExtra{'stopAfter'}} @allowedStopAfter );
    # Build the list of selections if we do not yet have it.
    @{$simulations->{'selections'}} = &selectSimulations(\%options)
	unless ( exists($simulations->{'selections'}) );
    # Construct the list.
    my @simulationList;
    # Iterate over suites.
    foreach my $suite ( &List::ExtraUtils::hashList($simulations->{'suite'}, keyAs => "name") ) {
	next
	    unless ( &matchSelection($simulations->{'selections'},$suite->{'name'}) );
	# Find the cosmological parameters for this suite.
	unless ( exists($suite->{'cosmology'}) ) {
	    my $xml = new XML::Simple();
	    my $cosmologyParameters = $xml->XMLin($options{'pipelinePath'}."cosmology_".$suite->{'name'}.".xml");
	    $suite->{'cosmology'}->{$_} = $cosmologyParameters->{'cosmologyParameters'}->{$_}->{'value'}
	        foreach ( 'HubbleConstant', 'OmegaMatter', 'OmegaDarkEnergy', 'OmegaBaryon' );
	    $suite->{'cosmology'}->{'hubbleConstant'} = $suite->{'cosmology'}->{'HubbleConstant'}/100.0;
	}
	# Push to the list and move on if we are to stop after the "suite" stage.  
	if ( $optionsExtra{'stopAfter'} eq "suite" ) {
	    push(@simulationList,{suite => $suite});
	    next;
	}
 	# Iterate over groups in the suite.
	foreach my $group ( &List::ExtraUtils::hashList($suite->{'group'}, keyAs => "name") ) {
	    next
		unless ( &matchSelection($simulations->{'selections'},$suite->{'name'},$group->{'name'}) );
	    # Find metadata.
	    unless ( exists($group->{'metaData'}) ) {
		my $xml = new XML::Simple();
		my $simulationParameters = $xml->XMLin($options{'pipelinePath'}."simulation_".$suite->{'name'}."_".$group->{'name'}.".xml");
		$group->{'metaData'}->{'reference'} = $simulationParameters->{'simulation'}->{'simulationReference'}->{'value'};
		$group->{'metaData'}->{'url'      } = $simulationParameters->{'simulation'}->{'simulationURL'      }->{'value'};
	    }
	    # Find best available resolution.
	    my @resolutions = sort map {(my $resolution = $_->{'name'}) =~ s/^resolutionX//; $resolution*1.0} &List::ExtraUtils::hashList($group->{'resolution'}, keyAs => "name");
	    $group->{'resolutionBest'} = "resolutionX".$resolutions[-1];
	    # Push to the list and move on if we are to stop after the "group" stage.  	    
	    if ( $optionsExtra{'stopAfter'} eq "group" ) {
		push(@simulationList,{suite => $suite, group => $group});
		next;
	    }
	    # Iterate over resolutions in the group.
	    foreach my $resolution ( &List::ExtraUtils::hashList($group->{'resolution'}, keyAs => "name") ) {
		next
		    unless ( &matchSelection($simulations->{'selections'},$suite->{'name'},$group->{'name'},$resolution->{'name'}) );
		# Find the number of subvolumes for this resolution.
		unless ( exists($resolution->{'subvolumes'}) ) {
		    my $xml = new XML::Simple();
		    my $simulationParameters = $xml->XMLin($options{'pipelinePath'}."simulation_".$suite->{'name'}."_".$group->{'name'}.".xml");	
		    $resolution->{'subvolumes'}   = $simulationParameters->{'simulation'}->{'subvolumes'}->{$resolution->{'name'}}->{'value'};
		}
		# Find the particle mass for this resolution.
		unless ( exists($resolution->{'massParticle'}) ) {
		    my $xml = new XML::Simple();
		    my $simulationParameters = $xml->XMLin($options{'pipelinePath'}."simulation_".$suite->{'name'}."_".$group->{'name'}.".xml");
		    my $massParticle         = $simulationParameters->{'simulation'}->{'massParticle'}->{$resolution->{'name'}}->{'value'};
		    # Handle any Hubble parameter scaling.
		    if ( $massParticle =~ m/^=/ ) {
			$massParticle =~ s/^=//;
			$massParticle =~ s/\[cosmologyParameters\/HubbleConstant\]/$suite->{'cosmology'}->{'HubbleConstant'}/;
			$massParticle = eval($massParticle);
		    }
		    $resolution->{'massParticle'} = $massParticle;
		}
		# Extract lists of redshifts.
		my $expansionFactors = pdl split(" ",$resolution->{'haloMassFunction'}->{'expansionFactors'}->{'value'});
		@{$resolution->{'expansionFactors'}} = list(    $expansionFactors    );
		@{$resolution->{'redshifts'       }} = list(1.0/$expansionFactors-1.0);
		# Extract lists of snapshots.
		@{$resolution->{'snapshots'       }} = split(" ",$resolution->{'haloMassFunction'}->{'snapshots'}->{'value'})
		if ( exists($resolution->{'haloMassFunction'}->{'snapshots'}) );
		# Push to the list and move on if we are to stop after the "resolution" stage.  	    
		if ( $optionsExtra{'stopAfter'} eq "resolution" ) {
		    push(@simulationList,{suite => $suite, group => $group, resolution => $resolution});
		    next;
		}
		# Iterate over simulations in the group.
		foreach my $simulation ( &List::ExtraUtils::hashList($resolution->{'simulation'}, keyAs => "name") ) {
		    next
			unless ( &matchSelection($simulations->{'selections'},$suite->{'name'},$group->{'name'},$resolution->{'name'},$simulation->{'name'}) );
		    # Extract lists of realizations.
		    @{$simulation->{'realizationsList'}} = exists($simulation->{'realizations'}) ? split(" ",$simulation->{'realizations'}->{'value'}) : ( "only" );
		    # Push to the list and move on if we are to stop after the "simulation" stage.  	    
		    if ( $optionsExtra{'stopAfter'} eq "simulation" ) {
			push(@simulationList,{suite => $suite, group => $group, resolution => $resolution, simulation => $simulation});
			next;
		    }		
		    # Iterate over realizations in the simulation.
		    foreach my $realization ( @{$simulation->{'realizationsList'}} ) {
			# Find best available resolution for this realization.
			my $isBestResolution = 1;
			(my $resolutionNumerical = $resolution->{'name'}) =~ s/^resolutionX//;
			$resolutionNumerical *= 1.0;
			foreach my $resolutionOther ( &List::ExtraUtils::hashList($group->{'resolution'}, keyAs => "name") ) {
			    (my $resolutionOtherNumerical = $resolutionOther->{'name'}) =~ s/^resolutionX//;
			    $resolutionOtherNumerical *= 1.0;
			    foreach my $simulationOther ( &List::ExtraUtils::hashList($resolutionOther->{'simulation'}, keyAs => "name") ) {
				next
				    unless ( $simulationOther->{'name'} eq $simulation->{'name'} );
				# Extract lists of realizations.
				my @realizationsListOther = exists($simulationOther->{'realizations'}) ? split(" ",$simulationOther->{'realizations'}->{'value'}) : ( "only" );
				my $realizationIsPresent = grep {$_ eq $realization} @realizationsListOther;
				next
				    unless ( $realizationIsPresent );
				$isBestResolution = 0
				    if ( $resolutionOtherNumerical > $resolutionNumerical );
			    }
			}
			next
			    unless ( &matchSelection($simulations->{'selections'},$suite->{'name'},$group->{'name'},$resolution->{'name'},$simulation->{'name'},$realization,$isBestResolution) );
			# Push to the list and move on if we are to stop after the "realization" stage.  	    
			if ( $optionsExtra{'stopAfter'} eq "realization" ) {
			    push(@simulationList,{suite => $suite, group => $group, resolution => $resolution, simulation => $simulation, realization => $realization});
			    next;
			}
			# Iterate over redshifts in the simulation.
			foreach my $redshift ( @{$resolution->{'redshifts'}} ) {
			    my $redshiftShort = sprintf("%5.3f",$redshift);
			    next
				unless ( &matchSelection($simulations->{'selections'},$suite->{'name'},$group->{'name'},$resolution->{'name'},$simulation->{'name'},$realization,$isBestResolution,$redshiftShort) );
			    # Begin constructing the complete entry that will be pushed to the list.
			    my $entry = {suite => $suite, group => $group, resolution => $resolution, simulation => $simulation, redshift => $redshiftShort, realization => $realization};
			    # Set the target data file.
			    $entry->{'fileTargetData'} = 
				$ENV{'GALACTICUS_DATA_PATH'}             .
				"/static/darkMatter/haloMassFunction_"   .
			            $entry->{'suite'      }->{'name'}."_".
			            $entry->{'group'      }->{'name'}."_".
			            $entry->{'resolution' }->{'name'}."_".
			            $entry->{'simulation' }->{'name'}."_".
			            $entry->{'realization'}          ."_".
				"z".$entry->{'redshift'   }              .
				".hdf5";
			    # Extract data from the target data file.
			    die("target data file '".$entry->{'fileTargetData'}."' does not exist")
				unless ( -e $entry->{'fileTargetData'} );
			    my $dataTarget       = new PDL::IO::HDF5($entry->{'fileTargetData'});
			    my $simulationTarget = $dataTarget->group('simulation0001');
			    # Extract properties of the primary halo from the target data file.
			    ($entry->{'massPrimary'}) = $simulationTarget->attrGet('massPrimary')
				if ( $entry->{'suite'}->{'limitMassMaximum'}->{'value'} eq "primaryFraction" );
			    # Extract properties of the environment, if needed.
			    if ( $entry->{'suite'}->{'includeEnvironment'}->{'value'} eq "true" ) {
				foreach my $attributeName ( 'massEnvironment', 'overdensityEnvironment' ) {
				    die("Attribute '".$attributeName."' is missing in file ".$entry->{'fileTargetData'})
					unless ( grep {$_ eq $attributeName} $simulationTarget->attrs() );
				    ($entry->{'environment'}->{$attributeName}) = $simulationTarget->attrGet($attributeName);
				}
			    }
			    # Push to the list and move on if we are to stop after the "redshift" stage.  	    
			    if ( $optionsExtra{'stopAfter'} eq "redshift" ) {
				push(@simulationList,$entry);
				next;
			    }
			}
		    }
		}
	    }
	}
    }
    # Re-order to spread models with the same power spectrum class. This optimizes load balancing when recalculation of σ(M) is
    # required, allowing different MPI processes to begin with models that have different power spectra and so will store these to
    # different files.
    if ( exists( $options{'reOrder'}) && $options{'reOrder'} eq "yes" && $optionsExtra{'stopAfter'} eq "redshift" ) {
	my %powerSpectrumClasses;
	foreach my $entry ( @simulationList ) {
	    push(@{$powerSpectrumClasses{$entry->{'simulation'}->{'powerSpectrumClass'}->{'value'}}},$entry);
	}
	my @simulationListReordered = ();
	while ( scalar(keys(%powerSpectrumClasses)) > 0 ) {
	    foreach my $powerSpectrumClass ( sort(keys(%powerSpectrumClasses)) ) {
		push(@simulationListReordered,pop(@{$powerSpectrumClasses{$powerSpectrumClass}}));
		delete($powerSpectrumClasses{$powerSpectrumClass})
		    if ( scalar(@{$powerSpectrumClasses{$powerSpectrumClass}}) == 0 );
	    }
	}
	@simulationList = @simulationListReordered;
	foreach my $entry ( @simulationList ) {
	    print $entry->{'simulation'}->{'powerSpectrumClass'}->{'value'}."\n";
	}
    }
    # Return the iterable list that we have constructed.
    return @simulationList;
}

sub selectSimulations {
    # Construct a list of selected simulations.
    my %options = %{shift()};
    my @selections;
    foreach my $selection ( &List::ExtraUtils::as_array($options{'select'}) ) {
	# Split the selection into sub-selections.
	my @subselections = split(/::/,$selection);
	foreach my $subselection ( @subselections ) {
	    my @splitSelection = split(/,/,$subselection);
	    $subselection = \@splitSelection;
	}
	push(
	    @selections,
	    {
		suite       => $subselections[0],
		group       => $subselections[1],
		resolution  => $subselections[2],
		simulation  => $subselections[3],
		realization => $subselections[4],
		redshift    => $subselections[5]
	    }
	    );
    }
    die("No entries matched the provided selections")
	if ( exists($options{'select'}) && scalar(@selections) == 0 );
    return @selections;
}

sub matchSelection {
    # Test if the current simulation matches a selection.
    my @selections       = @{shift()};
    my $suite            =   shift() ;
    my $group            =   shift() ;
    my $resolution       =   shift() ;
    my $simulation       =   shift() ;
    my $realization      =   shift() ;
    my $isBestResolution =   shift() ;
    my $redshift         =   shift() ;
    # If there are no selections, everything is a match.
    return 1
	if ( scalar(@selections) == 0 );
    foreach my $selection ( @selections ) {
	next unless ( ! defined($selection->{'suite'      }) || ! defined($suite      ) || grep {$_ eq $suite       || $_ eq "*"} @{$selection->{'suite'      }} );
	next unless ( ! defined($selection->{'group'      }) || ! defined($group      ) || grep {$_ eq $group       || $_ eq "*"} @{$selection->{'group'      }} );
	my $resolutionMatches;
	if ( defined($selection->{'resolution'}) && defined($resolution) && defined($isBestResolution) ) {
	    foreach my $selectResolution ( @{$selection->{'resolution'}} ) {
		if ( $selectResolution eq "best" ) {
		    $resolutionMatches = $isBestResolution;
		} else {
		    $resolutionMatches = $selectResolution eq $resolution || $selectResolution eq "*";
		}
	    }
	} else {
	    $resolutionMatches = 1;
	}
	next unless ( $resolutionMatches );
	next unless ( ! defined($selection->{'simulation' }) || ! defined($simulation ) || grep {$_ eq $simulation  || $_ eq "*"} @{$selection->{'simulation' }} );
	next unless ( ! defined($selection->{'redshift'   }) || ! defined($redshift   ) || grep {$_ eq $redshift    || $_ eq "*"} @{$selection->{'redshift'   }} );
	next unless ( ! defined($selection->{'realization'}) || ! defined($realization) || grep {$_ eq $realization || $_ eq "*"} @{$selection->{'realization'}} );
	return 1;
    }
    return 0
}

1;
