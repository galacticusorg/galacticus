# Contains a Perl module which implements various useful functionality for handling covariance matrices in
# Galacticus when fitting to constraints.

package Covariances;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::MatrixOps;
use PDL::Matrix;
use PDL::LinearAlgebra;
use Data::Dumper;

sub LogNormalCovariance {
    # Given a vector, y, and its covariance matrix, C, compute the covariance in log(y). C_ln = ln(1+C/y*y) where "*" indicates
    # the outer product.
    my $y                     = shift;
    my $C                     = shift;
    # Compute the outer product. Where the covariance is zero, set the outer product to unity to avoid division by zero errors.
    my $outerY                = outer($y,$y);
    $outerY    ->($C==0.0;?) .= 1.0;
    # Compute the ratio of covariance to outer product.
    my $r                     = $C/$outerY;
    # Comptue the covariance in log(y). Where the covariance is small, use a Taylor expansion to maintain accuracy.
    my $eps                   = pdl 1.0e-9;
    my $CLog                  = log(1.0+$r);
    $CLog->(abs($r)<$eps;?)  .= $r->(abs($r)<$eps;?)*(1.0+$r->(abs($r)<$eps;?)*(-1.0/2.0+$r->(abs($r)<$eps;?)*1.0/3.0));
    return $CLog;
}

sub SVDInvert {
    # Invert a covariance matrix using Singular Value Decomposition.
    my $C = shift;
    # Get any options.
    my %options;
    (%options) = @_
	if ( scalar(@_) > 0 );
    # Set default options.
    $options{'errorTolerant'} = 0
	unless ( exists($options{'errorTolerant'}) );
    # Do the Singular Value Decomposition.
    (my $r1, my $s, my $r2)                             = svd($C);
    # Invert the matrix.
    my $nonZeroSingularValues                           = which($s > 0.0);
    my $sInverse                                        = zeroes($r1);
    $sInverse->diagonal(0,1)->($nonZeroSingularValues) .= 1.0/$s->($nonZeroSingularValues);
    my $CInverse                                        = $r1 x $sInverse x transpose($r2);
    # Force the inverse covariance matrix to be semi-positive definite.
    my $CInverseSPD                                     = &MakeSemiPositiveDefinite($CInverse);
    # Compute the log of the determinant of the covariance matrix.
    my $logDeterminant                                  = sum(log($s->($nonZeroSingularValues)));
    # Perform sanity checks on the inverse covariance matrix and its determinant.
    (my $eigenVectors, my $eigenValues)                 = eigens_sym($CInverseSPD);
    print "SVDInvert: inverse covariance matrix is not semi-positive definite\n"
	unless ( all(         $eigenValues >= 0.0)  ); 
    unless (     isfinite($logDeterminant)  ) {
	print "SVDInvert: covariance matrix determinant failed\n";
	die
	    unless ( $options{'errorTolerant'} == 1 );
    }
    unless ( all(isfinite($CInverseSPD   )) ) {
	print "SVDInvert: covariance matrix inversion failed\n";
	die
	    unless ( $options{'errorTolerant'} == 1 );
    }
    return $CInverseSPD, $logDeterminant;
}

sub MakeSemiPositiveDefinite {
    # Force a matrix to be semi-positive definite by decomposing it into its eigenvectors, setting any negative eigenvalues to
    # zero, and reconstructing the original matric from the eigenvectors and (modified) eigenvalues.
    my $C                               = shift;
    # Decompose into eigenvectors and eigenvalues.
    (my $eigenVectors, my $eigenValues) = eigens_sym($C);
    # Test for non-semi-positive definiteness.
    if ( any($eigenValues < 0.0) ) {
	# Force any negative eigenvalues to zero.
	$eigenValues->($eigenValues<0.0;?) .= 0.0;
	# Reconstruct the matrix using these modified eigenvalues.
	my $eigenValuesMatrix               = zeroes($C);
	$eigenValuesMatrix->diagonal(0,1)  .= $eigenValues;
	my $CSPD                            = $eigenVectors x $eigenValuesMatrix x transpose($eigenVectors);
	return $CSPD;
    } else {
	# Matrix is already semi-positive definite, so return it unchanged.
	return $C;
    }
}

sub ComputeLikelihood {
    # Compute the likelihood given two y-vectors and their covariance matrix.
    my $y1                            = shift;
    my $y2                            = shift;
    my $C                             = shift;
    # Get any options.
    my %options;
    (%options) = @_
	if ( scalar(@_) > 0 );
    # Find the difference between the vectors.
    my $d                             = $y1-$y2;
    # Where upper limits are present, truncate the difference vector.
    if ( exists($options{'upperLimits'}) ) {
	my $limitTruncate = which($d->($options{'upperLimits'}) > 0.0);
	$d->($options{'upperLimits'})->($limitTruncate) .= 0.0;
    }
    # Invert the covariance matrix.
    (my $CInverse,my $logDeterminant) = &SVDInvert($C);
    # Construct the likelihood.
    my $vCv                           = $d x $CInverse x transpose($d);
    die("ComputeLikelihood: inverse covariance matrix is not semi-positive definite")
	unless ( $vCv->((0),(0)) >= 0.0 );
    my $logLikelihoodLog              = -0.5*$vCv->((0),(0))-0.5*nelem($y1)*log(2.0*3.1415927)-0.5*$logDeterminant;
    # Optionally return the logarithm of the determinant.
    ${$options{'determinant'}} = $logDeterminant
	if ( exists($options{'determinant'}) );
    # Return the likelihood.
    return $logLikelihoodLog->sclr();
}

1;
