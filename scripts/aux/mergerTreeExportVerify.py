#!/usr/bin/env python3
import numpy as np
import h5py
import sys
import argparse

# Test script which verifies that merger trees exported from Galacticus are consistent in structure (an key attributes of nodes)
# with the original imported trees.
# Andrew Benson (16-December-2025)

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='mergerTreeExportVerify.py',description='Verify that an exported merger tree matches the original.')
parser.add_argument('originalTreeFileName',help='The file name of the original merger tree file.')
parser.add_argument('exportedTreeFileName',help='The file name of the exported merger tree file.')
args = parser.parse_args()

# Set status.
success = True

# Define SI units needed for conversions.
massSolar  = 1.98910e+30
megaParsec = 3.08568e+22
kilo       = 1.00000e+03

# Define the merger tree files to be compared.
trees = {}
trees['original'] = {"file": h5py.File(args.originalTreeFileName,'r')}
trees['new'     ] = {"file": h5py.File(args.exportedTreeFileName,'r')}

# Read required propteries from the trees.
for tree in trees:
    for property in ( 'nodeIndex', 'hostIndex', 'descendantIndex', 'nodeMass', 'scaleRadius', 'position', 'velocity' ):
        trees[tree][property] = trees[tree]['file']['forestHalos'][property][:]
    # Derive redshift/expansion factor as needed.
    if tree == 'original':
        trees[tree]['expansionFactor'] = trees[tree]['file']['forestHalos']['expansionFactor'][:]
        trees[tree]['redshift'       ] = 1.0/trees[tree]['expansionFactor']-1.0
    else:
        trees[tree]['redshift'       ] = trees[tree]['file']['forestHalos']['redshift'][:]
        trees[tree]['expansionFactor'] = 1.0/(1.0+trees[tree]['redshift'])
    # Get unit (and cosmology) conversion factors.
    for unit in ( 'mass', 'length', 'velocity' ):
        for attribute in ( 'ScaleFactorExponent', 'HubbleExponent', 'UnitsInSI' ):
            trees[tree][unit+attribute] = trees[tree]['file']['units'].attrs[unit+attribute]
    trees[tree]['HubbleParam'] = trees[tree]['file']['cosmology'].attrs['HubbleParam']
    # Perform unit conversion to Galacticus standards.
    trees[tree]['nodeMass'   ] *= trees[tree]['HubbleParam']**trees[tree][  'massHubbleExponent']*trees[tree]['expansionFactor']**trees[tree][  'massScaleFactorExponent']*(trees[tree][  'massUnitsInSI']/massSolar )
    trees[tree]['scaleRadius'] *= trees[tree]['HubbleParam']**trees[tree]['lengthHubbleExponent']*trees[tree]['expansionFactor']**trees[tree]['lengthScaleFactorExponent']*(trees[tree]['lengthUnitsInSI']/megaParsec)
    for i in range(3):
        trees[tree]['position'][:,i] *= trees[tree]['HubbleParam']**trees[tree][  'lengthHubbleExponent']*trees[tree]['expansionFactor']**trees[tree][  'lengthScaleFactorExponent']*(trees[tree][  'lengthUnitsInSI']/megaParsec)
        trees[tree]['velocity'][:,i] *= trees[tree]['HubbleParam']**trees[tree]['velocityHubbleExponent']*trees[tree]['expansionFactor']**trees[tree]['velocityScaleFactorExponent']*(trees[tree]['velocityUnitsInSI']/kilo      )

# Nodes which were cloned in Galacticus are output with an offset in their node index (to avoid conflicting with the original
# node). Determine the clone offset.
nodeIndexMaximum = np.max(trees['original']['nodeIndex'])
nodeIndexOffset  = int(10**np.ceil(np.log10(nodeIndexMaximum)))

# Iterate through nodes in the new tree.
for i in range(len(trees['new']['nodeIndex'])):
    if i % 1000 == 0:
        print(f"Processing node {i} of {len(trees['new']['nodeIndex'])}")
    # Find the corresponding node in the original tree.
    matchedNode = trees['original']['nodeIndex'] == trees['new']['nodeIndex'][i]
    if np.count_nonzero(matchedNode) == 0:
        # No match - check if this is a clone.
        if trees['new']['nodeIndex'][i] > nodeIndexOffset:
            # Clone.
            continue
        else:
            print("Unmatched index")
            sys.exit(1)
    elif np.count_nonzero(matchedNode) > 1:
        # Multiple matches.
        print("Multiple matches found")
        sys.exit(1)
    else:
        # A single match was found - proceed to verify it is correct.
        #
        # Only single level hierarchies are supported by the mergerTreeConstructorRead class.  Therefore, first move to the
        # earliest progenitor, then follow this forward, moving to the ultimate host at each step, to trace the actual host that
        # will have been assigned by Galacticus.
        #
        # Array index of the original matched node.
        indexOriginal = np.where(matchedNode)[0][0]
        # While this is a subhalo, trace.
        while trees['original']['hostIndex'][indexOriginal] != trees['original']['nodeIndex'][indexOriginal]:
            # Find the ultimate host of the current traced node.
            host = np.where(trees['original']['nodeIndex'] == trees['original']['hostIndex'][indexOriginal])[0]
            while trees['original']['hostIndex'][host][0] != trees['original']['nodeIndex'][host][0]:
                host = np.where(trees['original']['nodeIndex'] == trees['original']['hostIndex'][host])[0]
            # Find the descendant of the traced node.
            descendant = np.where(trees['original']['nodeIndex'] == trees['original']['descendantIndex'][indexOriginal])[0]
            # If we have both a descendant and a host, and the redshift is prior to the original node (i.e. we have already
            # stepped back to a progenitor), then check for branch jumps.
            if trees['original']['descendantIndex'][host] >= 0 and trees['original']['hostIndex'][descendant] >= 0 and trees['original']['redshift'][indexOriginal] > trees['original']['redshift'][matchedNode][0]:
                # Find the ultimate host's descendant node.
                hostDescendant = np.where(trees['original']['nodeIndex'] == trees['original']['descendantIndex'][host      ])[0][0]
                # Find the descendant node's ultimate host.
                descendantHost = np.where(trees['original']['nodeIndex'] == trees['original'][      'hostIndex'][descendant])[0][0]
                while trees['original']['hostIndex'][descendantHost] != trees['original']['nodeIndex'][descendantHost]:
                    descendantHost = np.where(trees['original']['nodeIndex'] == trees['original']['hostIndex'][descendantHost])[0]
                # If these two do not match it indicates a branch jump - which Galacticus will respect. Therefore, the host node
                # is just the ultimate host of the descendant. Set that an exit tracing.
                if hostDescendant != descendantHost:
                    indexOriginal = descendantHost
                    break
            # Move the the most massive progenitor node.
            progenitors = trees['original']['descendantIndex'] == trees['original']['nodeIndex'][indexOriginal]
            if np.count_nonzero(progenitors) == 0:
                # No progenitors found - exit tracing.
                break
            mostMassiveProgenitor = np.argmax(trees['original']['nodeMass'][progenitors])
            indexOriginal         = np.where(progenitors)[mostMassiveProgenitor][0]
        # Having traced back through progenitors (or until the most recent branch jump), now trace forward, moving to the ultimate
        # host at each step.
        while trees['original']['redshift'][indexOriginal] >= trees['original']['redshift'][matchedNode][0]:
            # Move to the ultimate host.
            while trees['original']['hostIndex'][indexOriginal] != trees['original']['nodeIndex'][indexOriginal]:
                hostNext      = trees['original']['nodeIndex'] == trees['original']['hostIndex'][indexOriginal]
                indexOriginal = np.where(hostNext)[0]
            if trees['original']['redshift'][indexOriginal] == trees['original']['redshift'][matchedNode][0]:
                # Exit if the redshift of the original node is reached.
                break
            else:
                # Move to the descendant node.
                descendant    = trees['original']['nodeIndex'] == trees['original']['descendantIndex'][indexOriginal]
                indexOriginal = np.where(descendant)[0]
        # Look for mismatches between the properties of the exported node and the original. Here we append reasons for failure to
        # a list.
        failures = []
        if trees['new']['descendantIndex'][i] != trees['original']['descendantIndex'][matchedNode  ][0]:
            failures.append("descendant indices do not match")
        if trees['new'][      'hostIndex'][i] != trees['original'][      'hostIndex'][indexOriginal]:
            failures.append("host indices do not match")
        if not np.isclose (trees['new']['redshift'][i  ],trees['original']['redshift'][matchedNode  ][0],atol=1.0e-3):
            failures.append("redshifts do not match")
        if not np.isclose (trees['new']['nodeMass'][i  ],trees['original']['nodeMass'][matchedNode  ][0],rtol=1.0e-3):
            failures.append("node masses do not match")
        if not np.allclose(trees['new']['position'][i,:],trees['original']['position'][matchedNode,:][0],rtol=1.0e-3):
            failures.append("positions do not match")
        if not np.allclose(trees['new']['velocity'][i,:],trees['original']['velocity'][matchedNode,:][0],rtol=1.0e-3):
            failures.append("velocities do not match")
        # Tests for isolated halos only - subhalos do not track the scale radius, this is computed internally for them.
        if trees['new']['nodeIndex'][i] == trees['new']['hostIndex'][i]:
            if not np.isclose(trees['new']['scaleRadius'][i],trees['original']['scaleRadius'][matchedNode][0],rtol=1.0e-3):
                failures.append("scale radii do not match")
        # Report on any failures.
        if len(failures) > 0:
            success = False
            print(f"Node {trees['new']['nodeIndex'][i]} failed because:")
            for failure in failures:
                print(f"   --> {failure}")
            print(f"   Descendant index: {trees['new']['descendantIndex'][i  ]} {trees['original']['descendantIndex'][matchedNode  ][0]}")
            print(f"   Host index:       {trees['new']['hostIndex'      ][i  ]} {trees['original']['hostIndex'      ][matchedNode  ][0]}")
            print(f"   Host index [top]: {trees['new']['hostIndex'      ][i  ]} {trees['original']['hostIndex'      ][indexOriginal]   }")
            print(f"   Redshift:         {trees['new']['redshift'       ][i  ]} {trees['original']['redshift'       ][matchedNode  ][0]}")
            print(f"   Node mass:        {trees['new']['nodeMass'       ][i  ]} {trees['original']['nodeMass'       ][matchedNode  ][0]}")
            print(f"   Scale radius:     {trees['new']['scaleRadius'    ][i  ]} {trees['original']['scaleRadius'    ][matchedNode  ][0]}")
            print(f"   Position [x]:     {trees['new']['position'       ][i,0]} {trees['original']['position'       ][matchedNode,0][0]}")
            print(f"   Position [y]:     {trees['new']['position'       ][i,1]} {trees['original']['position'       ][matchedNode,1][0]}")
            print(f"   Position [z]:     {trees['new']['position'       ][i,2]} {trees['original']['position'       ][matchedNode,2][0]}")
            print(f"   Velocity [x]:     {trees['new']['velocity'       ][i,0]} {trees['original']['velocity'       ][matchedNode,0][0]}")
            print(f"   Velocity [y]:     {trees['new']['velocity'       ][i,1]} {trees['original']['velocity'       ][matchedNode,1][0]}")
            print(f"   Velocity [z]:     {trees['new']['velocity'       ][i,2]} {trees['original']['velocity'       ][matchedNode,2][0]}")

# Report status.
if not success:
    print("Failures found")
    sys.exit(1)
