# Contains a Perl module which implements dependency sorting of tasks.

package Galacticus::Build::Dependencies;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Sort::Topo;
use Scalar::Util 'reftype';
use Data::Dumper;

sub Dependency_Sort {
    # Generate an array of sorted names based on the provided hash containing dependency information.
    my $tasks     = shift;
    my $buildData = shift;

    # Initialize a hash of dependencies and sorting names.
    my %dependencies;
    my %sortNames;
    my %sortNamesReverse;

    # Generate empty lists of dependencies.
    foreach my $key ( keys(%{$tasks}) ) {
	# Get the name for sorting.
	my $sortName = $key;
	$sortName = $tasks->{$key}->{'sortName'}
	    if ( exists($tasks->{$key}->{'sortName'}) );
	$tasks->{$key}->{'sortName'} = $sortName;
	# Create an empty list of dependencies for this task.
	@{$dependencies{$sortName}} = ();
    }

    # Iterate over all keys in the provided hash.
    foreach my $key ( keys(%{$tasks}) ) {
	# Get the name for sorting.
	my $sortName = $key;
	$sortName = $tasks->{$key}->{'sortName'}
	    if ( exists($tasks->{$key}->{'sortName'}) );
     	# Check for "after" dependencies.
     	if ( exists($tasks->{$key}->{'after'}) ) {
     	    if ( reftype($tasks->{$key}->{'after'}) && reftype($tasks->{$key}->{'after'}) eq "ARRAY" ) {
     		foreach my $after ( @{$tasks->{$key}->{'after'}} ) {
		    my $afterSortName = $after;
		    $afterSortName = $tasks->{$after}->{'sortName'}
		        if ( exists($tasks->{$after}->{'sortName'}) );
     		    push(
     			@{$dependencies{$afterSortName}},
     			$sortName
     			);
     		}
     	    } else {
		my $after         = $tasks->{$key}->{'after'};
		my $afterSortName = $after;
		$afterSortName = $tasks->{$afterSortName}->{'sortName'}
		    if ( exists($tasks->{$afterSortName}->{'sortName'}) );
     		push(
     		    @{$dependencies{$afterSortName}},
     		    $sortName
     		    );
     	    }
     	}
     	# Check for "before" dependencies.
     	if ( exists($tasks->{$key}->{'before'}) ) {
     	    my @dependencyList;
     	    if ( reftype($tasks->{$key}->{'before'}) && reftype($tasks->{$key}->{'before'}) eq "ARRAY" ) {
     		@dependencyList = @{$tasks->{$key}->{'before'}};
     	    } else {
     		@dependencyList = ( $tasks->{$key}->{'before'} );
     	    }
	    foreach my $before ( @dependencyList ) {
		my $beforeSortName = $before;
		$beforeSortName = $tasks->{$before}->{'sortName'}
		        if ( exists($tasks->{$before}->{'sortName'}) );
     	    push(@{$dependencies{$sortName}},$beforeSortName);
	    }
     	}
	# Store the sort names and a reverse lookup hash.
	$sortNames{$key} = $tasks->{$key}->{'sortName'};
	push(@{$sortNamesReverse{$tasks->{$key}->{'sortName'}}},$key);
    }

    # Scan through dependency keys, searching for regular expressions.
    foreach my $key ( keys(%dependencies) ) {
	if ( $key =~ m/^re:/ ) {
	    (my $regEx = $key) =~ s/^re://;
	    # Iterate over dependency keys again, searching for matches.
	    foreach my $matchKey ( keys(%dependencies) ) {
		unless ( $key eq $matchKey ) {
		    if ( $matchKey =~ m/$regEx/ ) {
			push(@{$dependencies{$matchKey}},@{$dependencies{$key}});
		    }
		}
	    }
	    delete($dependencies{$key});
	}
    }
    # Iterate over dependency keys. Look for dependencies in terms of regular expressions.
    foreach my $key ( keys(%dependencies) ) {
	foreach my $dependency ( @{$dependencies{$key}} ) {
	    if ( $dependency =~ m/^re:/ ) {
		(my $regEx = $dependency) =~ s/^re://;
		# Iterate through dependency keys and find matches.
		foreach my $matchKey ( keys(%dependencies) ) {
		    push(@{$dependencies{$key}},$matchKey)
			if ( $matchKey =~ m/$regEx/ );
		}
	    }
	}
	@{$dependencies{$key}} = grep(!/^re:/,@{$dependencies{$key}});
    }

    # Create a list of names for sorting.
    my @unsortedUnits = keys(%dependencies);
   
    # Perform initial alphanumerical sort.
    my @presortedUnits = sort(@unsortedUnits);

    # Perform dependency sort.
    my @sortedSortNames = &Sort::Topo::sort(\@presortedUnits,\%dependencies);
    my @sortedUnits;
    foreach my $sortName ( @sortedSortNames ) {
	push(@sortedUnits,@{$sortNamesReverse{$sortName}});
    }

    # Store the sorted lists.
    @{$buildData->{'sortNames'}} = @sortedSortNames;
    @{$buildData->{'unitNames'}} = @sortedUnits;

}

1;
