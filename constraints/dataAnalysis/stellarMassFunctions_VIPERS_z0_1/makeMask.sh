#!/bin/sh

# Combine mangle polygon files for surveyed areas and masked regions supplied by I. Davidzon into
# a single mangle file. Follows the procedure described at https://jila.colorado.edu/~ajsh/mangle/index.html#complement

# Andrew Benson (13-June-2014)

# Define work directory.
work=constraints/dataAnalysis/stellarMassFunctions_VIPERS_z0_1
# Define mangle directory.
mangle=aux/mangle/bin
# Create an all sky file.
if [ ! -f $work/allsky.pol ]; then
 cat > $work/allsky.pol <<EOF
1 polygons
snapped
polygon 0 ( 0 caps, 1 weight, 0 pixel, 12.566370614359172954 str):
EOF
fi
# Pixelize the all sky file.
if [ ! -f $work/allsky.ply ]; then
    $mangle/pixelize -Ps0,11 $work/allsky.pol $work/allsky.ply
fi
# Balkanize the survey regions.
if [ ! -f $work/maskCombinedB.ply ]; then
    $mangle/balkanize                   \
	$work/manglemasks/mask_*.mangle \
	$work/maskCombinedB.ply
fi
# Unify the surveyed regions.
if [ ! -f $work/maskCombinedBU.ply ]; then
    $mangle/unify                      \
	$work/maskCombinedB.ply        \
	$work/maskCombinedBU.ply
fi
# Holeize the survey regions.
if [ ! -f $work/maskCombinedBUH.ply ]; then
    sed -r s/" 1 weight"/" 0 weight"/ \
	$work/maskCombinedBU.ply >    \
	$work/maskCombinedBUH.ply
fi
# Balkanize all sky with holeized survey regions.
if [ ! -f $work/allskyB.ply ]; then
    $mangle/balkanize             \
	$work/allsky.ply          \
	$work/maskCombinedBUH.ply \
	$work/allskyB.ply
fi
# Unify the balkanized all sky with holeized survery regions.
if [ ! -f $work/allskyBU.ply ]; then
    $mangle/unify          \
	$work/allskyB.ply  \
	$work/allskyBU.ply
fi
# Remove zero weight and zero area poygons from the balkanized all sky with holeized survery regions.
if [ ! -f $work/allskyBU2.ply ]; then
    $mangle/poly2poly -k1.0e-16 -j1.0e-16  \
	$work/allskyBU.ply                 \
	$work/allskyBU2.ply
fi
# Create hole all sky with survey regions removed.
if [ ! -f $work/allskyBU2H.ply ]; then
    sed -r s/" 1 weight"/" 0 weight"/ \
	$work/allskyBU2.ply >          \
	$work/allskyBU2H.ply
fi
# Balkanize photometric holes with hole all sky with survey regions removed.
for field in W1 W4; do
    if [ ! -f $work/photo$field.B.ply ]; then
	$mangle/balkanize                                \
	    $work/samhain_mangle/photo_$field.reg.mangle \
	    $work/allskyBU2H.ply                         \
	    $work/photo$field.B.ply
    fi
done
# Unify photometric holes in survey regions.
for field in W1 W4; do
    if [ ! -f $work/photo$field.BU.ply ]; then
	$mangle/unify                \
	    $work/photo$field.B.ply  \
	    $work/photo$field.BU.ply
    fi
done
# Remove zero weight and zero area polygons frmo the photometric holes in survey regions.
for field in W1 W4; do
    if [ ! -f $work/photo$field.BU2.ply ]; then
	$mangle/poly2poly -k1.0e-16 -j1.0e-16  \
	    $work/photo$field.BU.ply           \
	    $work/photo$field.BU2.ply
    fi
done

exit
