#!/usr/bin/env python3
# Analyze a variety of cosmological N-body simulations to extract statistics of interest.
# Python port of constraints/pipelines/darkMatter/simulationsAnalyze.pl
# Andrew Benson (ported to Python 2026)

import argparse
import math
import os
import re
import shutil
import sys
from datetime import datetime

import h5py
import lxml.etree as ET
import numpy as np

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))
import queueManager
from Galacticus.Constraints.Simulations import iterate


# ---------------------------------------------------------------------------
# XML parsing helpers
# ---------------------------------------------------------------------------

def _element_to_dict_keyed(el, keyed_tags):
    """Recursively convert an lxml element to a nested dict.

    Elements whose tag is in keyed_tags are converted to dicts keyed by their
    'name' attribute, matching XML::Simple's KeyAttr behaviour. All other
    children with a unique tag become plain dict values; repeating tags become
    lists. An element's own XML attributes are included in the result dict.
    """
    result = dict(el.attrib)
    children_by_tag = {}
    for child in el:
        children_by_tag.setdefault(child.tag, []).append(child)
    for tag, children in children_by_tag.items():
        if tag in keyed_tags:
            result[tag] = {
                child.get('name'): _element_to_dict_keyed(child, keyed_tags)
                for child in children
            }
        elif len(children) == 1:
            result[tag] = _element_to_dict_keyed(children[0], keyed_tags)
        else:
            result[tag] = [_element_to_dict_keyed(c, keyed_tags) for c in children]
    return result


def parse_simulations_xml(path):
    """Parse simulations.xml into a nested dict matching XML::Simple KeyAttr output.

    suite, group, resolution, and simulation elements are keyed by their
    'name' attribute in the result dict.
    """
    keyed_tags = {'suite', 'group', 'resolution', 'simulation'}
    tree = ET.parse(path)
    root = tree.getroot()
    result = {}
    suites = root.findall('suite')
    if suites:
        result['suite'] = {
            el.get('name'): _element_to_dict_keyed(el, keyed_tags)
            for el in suites
        }
    return result


def _parse_param_xml(path):
    """Return the root element of a Galacticus parameter XML file."""
    return ET.parse(path).getroot()


def _write_param_xml(root, path):
    """Serialise a Galacticus parameter XML element tree to file."""
    ET.indent(root)
    tree = ET.ElementTree(root)
    tree.write(path, xml_declaration=True, encoding='utf-8')


def _find_or_create(parent, tag):
    """Return the child element with the given tag, creating it if absent."""
    el = parent.find(tag)
    if el is None:
        el = ET.SubElement(parent, tag)
    return el


# ---------------------------------------------------------------------------
# Job submission helpers
# ---------------------------------------------------------------------------

def _translate_job(job):
    """Translate Perl-style job dict keys to Python queueManager keys."""
    j = dict(job)
    # OMP thread count
    if 'ompThreads' in j:
        j['countOpenMPThreads'] = j.pop('ompThreads')
    # Memory: Perl uses 'mem', Python manager uses 'memory'
    if 'mem' in j:
        j['memory'] = j.pop('mem')
    # Log file: map to both output and error
    if 'logFile' in j:
        log = j.pop('logFile')
        j.setdefault('logOutput', log)
        j.setdefault('logError',  log)
    # ppn maps to tasksPerNode
    if 'ppn' in j:
        j['tasksPerNode'] = j.pop('ppn')
    return j


def _submit_jobs(manager, jobs):
    """Submit a list of Perl-style job dicts via the Python queue manager."""
    if jobs:
        manager.submitJobs([_translate_job(j) for j in jobs])


# ---------------------------------------------------------------------------
# Active-step resolution
# ---------------------------------------------------------------------------

def _build_active_steps(active_analyses):
    """Return a dict mapping step IDs to the analyses that require them."""
    active_steps = {}
    if 'haloMassFunction' in active_analyses:
        for step in ('identifyAlwaysIsolated', 'extractHalos', 'massFunctions'):
            active_steps.setdefault(step, []).append('haloMassFunction')
    if 'subhaloStatistics' in active_analyses:
        for step in ('identifyAlwaysIsolated', 'extractSubhalos', 'subhaloFunctions'):
            active_steps.setdefault(step, []).append('subhaloStatistics')
    return active_steps


# ---------------------------------------------------------------------------
# Pre/post-process dispatch
# ---------------------------------------------------------------------------

def _run_hooks(hook_key, step_id, entries, suites_cfg, manager, options):
    """Run all preprocess or postprocess hooks for a step in iteration order.

    Mirrors the Perl while(workDone) loop: iterates through entries calling
    the i-th hook function until no entry has an i-th hook, then submits jobs.
    """
    iteration = 0
    while True:
        jobs = []
        work_done = False
        for entry in entries:
            suite_name = entry['suite']['name']
            hooks = (suites_cfg
                     .get(suite_name, {})
                     .get('steps', {})
                     .get(step_id, {})
                     .get(hook_key, []))
            if iteration < len(hooks):
                work_done = True
                hooks[iteration](entry, jobs, options)
        _submit_jobs(manager, jobs)
        if not work_done:
            break
        iteration += 1


# ---------------------------------------------------------------------------
# Symphony / COZMIC processor stubs  (implemented in later steps)
# ---------------------------------------------------------------------------

def symphony_process_identify_always_isolated(entry, expansion_factor, parameters, options):
    """Set cosmology, extra import columns, and insert physicalToComoving operator."""
    realization = entry['realization']
    host_halo_ids = entry['simulation'].get('hostHaloIDs', {})
    if realization not in host_halo_ids:
        raise RuntimeError(
            f"cannot find host halo ID for {entry['suite']['name']}; "
            f"{entry['group']['name']}; {entry['resolution']['name']}; "
            f"{entry['simulation']['name']}; {realization}"
        )
    entry['hostHaloID'] = host_halo_ids[realization]

    # Append extra columns to the importer if not already present.
    importer     = parameters.find('nbodyImporter')
    current_cols = importer.find('readColumns').get('value', '').split()
    for prop in ('X', 'Y', 'Z', 'Rvir', 'rs', 'Vmax'):
        if prop not in current_cols:
            current_cols.append(prop)
    importer.find('readColumns').set('value', ' '.join(current_cols))

    # Remove position/radius/velocity properties from the delete-properties operator
    # (they must be kept so downstream operators can use them).
    nb_ops    = parameters.find('nbodyOperator').findall('nbodyOperator')
    keep_out  = {'position', 'radiusScale', 'radiusVirial', 'velocityMaximum'}
    to_delete = [p for p in nb_ops[3].find('propertyNames').get('value', '').split()
                 if p not in keep_out]
    nb_ops[3].find('propertyNames').set('value', ' '.join(to_delete))

    # Splice a physicalToComoving operator at index 2 (before the former ops[2]).
    new_op = ET.Element('nbodyOperator', value='physicalToComoving')
    parameters.find('nbodyOperator').insert(2, new_op)


def symphony_preprocess_extract_halos_locate(entry, jobs, options):
    """Identify the primary progenitor halo via zoomInExtract.py.

    Builds a single job whose command accumulates one zoomInExtract.py call per
    epoch whose primaryHalo_*.xml file does not yet exist.
    """
    realization   = entry['realization']
    host_halo_ids = entry['simulation'].get('hostHaloIDs', {})
    if realization not in host_halo_ids:
        raise RuntimeError(
            f"cannot find host halo ID for {entry['suite']['name']}; "
            f"{entry['group']['name']}; {entry['resolution']['name']}; "
            f"{entry['simulation']['name']}; {realization}"
        )
    host_halo_id = host_halo_ids[realization]

    job = None
    for epoch in entry['resolution']['epochs']:
        rl           = epoch['redshiftLabel']
        primary_file = entry['path'] + f'primaryHalo_{rl}.xml'
        if os.path.exists(primary_file):
            continue
        if job is None:
            mem = '16G'
            if entry['group'     ]['name'] == 'Group':
                mem = '32G'
            if entry['resolution']['name'] == 'resolutionX64':
                mem = '64G'
            s, g, sim, r = (entry['suite']['name'], entry['group']['name'],
                            entry['simulation']['name'], realization)
            job = {
                'command':    '',
                'launchFile': entry['path'] + f'zoomInExtract_{s}_{g}_{sim}_{r}.sh',
                'logFile':    entry['path'] + f'zoomInExtract_{s}_{g}_{sim}_{r}.log',
                'label':                      f'zoomInExtract_{s}_{g}_{sim}_{r}',
                'ppn':        1,
                'ompThreads': 1,
                'nodes':      1,
                'mem':        mem,
                'walltime':   '8:00:00',
                'mpi':        'no',
            }
        hmf = entry['resolution'].get('haloMassFunction', {})
        job['command'] += (
            options['pipelinePath'] + 'zoomInExtract.py'
            + ' ' + entry['path']
            + ' ' + primary_file
            + ' ' + str(epoch['expansionFactor'])
            + ' ' + str(entry['suite']['cosmology']['hubbleConstant'])
            + ' ' + str(entry['resolution']['massParticle'])
            + ' ' + hmf.get('massHostLogMin', {}).get('value', '')
            + ' ' + hmf.get('massHostLogMax', {}).get('value', '')
            + ' ' + str(host_halo_id)
            + '\n'
        )

    if job is not None:
        jobs.append(job)


def symphony_preprocess_extract_halos_uncontaminated(entry, jobs, options):
    """Find the uncontaminated high-resolution region around the primary halo.

    Reads primaryHalo_*.xml to cache the halo centre on entry, then submits a
    Galacticus.exe job to measure the uncontaminated radius if not done yet.
    """
    omp_threads = options.get('ompThreadsResolved', 1)
    for epoch in entry['resolution']['epochs']:
        rl   = epoch['redshiftLabel']
        a    = epoch['expansionFactor']
        af_label     = f'sphericalOrigin:a{a:.3f}'
        primary_file = entry['path'] + f'primaryHalo_{rl}.xml'

        # Cache halo centre coordinates and central mass from the primaryHalo XML.
        primary_root = ET.parse(primary_file).getroot()
        entry[af_label]      = (primary_root.get('x') + ' '
                                + primary_root.get('y') + ' '
                                + primary_root.get('z'))
        entry['massCentral'] = float(primary_root.get('mc'))

        uncontaminated_file = entry['path'] + f'uncontaminated_{rl}.hdf5'
        if os.path.exists(uncontaminated_file):
            continue

        # Find the snapshot index matching this expansion factor.
        expansion_factors = entry['resolution'].get('expansionFactors', [])
        snapshots         = entry['resolution'].get('snapshots', [])
        snapshot = next(
            (snapshots[i] for i, af in enumerate(expansion_factors) if af == a),
            None,
        )

        root = _parse_param_xml(options['pipelinePath'] + 'zoomInSelectUncontaminated.xml')
        root.find('nbodyImporter').find('fileName').set(
            'value', entry['path'] + f'snapshots/snapshot_{snapshot}')
        nb_ops = root.find('nbodyOperator').findall('nbodyOperator')
        nb_ops[0].find('point'   ).set('value', entry[af_label])
        nb_ops[1].find('fileName').set('value', uncontaminated_file)
        _find_or_create(root, 'outputFileName').set(
            'value', entry['path'] + f'uncontaminatedExtract_{rl}.hdf5')

        param_file = entry['path'] + f'uncontaminated_{rl}.xml'
        _write_param_xml(root, param_file)

        s, g, sim, r = (entry['suite']['name'], entry['group']['name'],
                        entry['simulation']['name'], entry['realization'])
        mem = '196G' if entry['resolution']['name'] == 'resolutionX64' else '32G'
        jobs.append({
            'command':    os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe ' + param_file,
            'launchFile': entry['path'] + f'uncontaminatedExtract_{s}_{g}_{sim}_{r}_{rl}.sh',
            'logFile':    entry['path'] + f'uncontaminatedExtract_{s}_{g}_{sim}_{r}_{rl}.log',
            'label':                      f'uncontaminatedExtract_{s}_{g}_{sim}_{r}_{rl}',
            'ppn':        omp_threads,
            'ompThreads': omp_threads,
            'nodes':      1,
            'mem':        mem,
            'walltime':   '8:00:00',
            'mpi':        'no',
        })


def symphony_process_extract_halos(entry, expansion_factor, parameters, options):
    """Set the spherical extraction region and splice distance/filter operators."""
    redshift      = 1.0 / expansion_factor - 1.0
    rl            = f'z{redshift:.3f}'
    af_label      = f'sphericalOrigin:a{expansion_factor:.3f}'
    primary_file  = entry['path'] + f'primaryHalo_{rl}.xml'
    unconts_file  = entry['path'] + f'uncontaminated_{rl}.hdf5'

    print(f"Processing halo extraction for {entry['suite']['name']} : "
          f"{entry['group']['name']} : {entry['simulation']['name']} : "
          f"{entry['realization']}")

    # Read the uncontaminated radius from the HDF5 output.
    with h5py.File(unconts_file, 'r') as hdf:
        radius_uncontaminated = float(
            hdf['Snapshot00001/HaloCatalog'].attrs['radiusUncontaminated'])

    # Record the uncontaminated radius in the primaryHalo XML so later steps can use it.
    primary_root = ET.parse(primary_file).getroot()
    box_size     = primary_root.get('l')
    primary_root.set('ur', str(radius_uncontaminated))
    with open(primary_file, 'w') as fh:
        fh.write(ET.tostring(primary_root, encoding='unicode'))

    # Append virial/scale radii and position to the importer property list.
    importer = parameters.find('nbodyImporter')
    props    = importer.find('properties').get('value', '')
    importer.find('properties').set('value', props + ' radiusScale radiusVirial position')

    # Extend the delete-properties operator to also remove distanceFromPoint.
    nb_ops = parameters.find('nbodyOperator').findall('nbodyOperator')
    del_props = nb_ops[1].find('propertyNames').get('value', '')
    nb_ops[1].find('propertyNames').set('value', del_props + ' distanceFromPoint')

    # Store the spherical radii on entry for use by the postprocess hooks.
    entry.setdefault(rl, {})['sphericalRadiusMinimum'] = 0.0
    entry[rl]['sphericalRadiusMaximum'] = radius_uncontaminated

    # Splice three operators at index 1: distanceFromPoint, filterProperties, setBoxSize.
    parent  = parameters.find('nbodyOperator')
    dist_op = ET.Element('nbodyOperator', value='distanceFromPoint')
    ET.SubElement(dist_op, 'point', value=entry[af_label])

    filt_op = ET.Element('nbodyOperator', value='filterProperties')
    ET.SubElement(filt_op, 'propertyNames', value='distanceFromPoint')
    ET.SubElement(filt_op, 'rangeLow',  value=str(entry[rl]['sphericalRadiusMinimum']))
    ET.SubElement(filt_op, 'rangeHigh', value=str(entry[rl]['sphericalRadiusMaximum']))

    box_op = ET.Element('nbodyOperator', value='setBoxSize')
    ET.SubElement(box_op, 'boxSize', value=str(box_size))

    parent.insert(1, dist_op)
    parent.insert(2, filt_op)
    parent.insert(3, box_op)


def symphony_process_extract_subhalos(entry, expansion_factor, parameters, options):
    """Filter subhalos to those within the host virial radius."""
    redshift     = 1.0 / expansion_factor - 1.0
    rl           = f'z{redshift:.3f}'
    af_label     = f'sphericalOrigin:a{expansion_factor:.3f}'
    primary_file = entry['path'] + f'primaryHalo_{rl}.xml'

    print(f"Processing subhalo extraction for {entry['suite']['name']} : "
          f"{entry['group']['name']} : {entry['simulation']['name']} : "
          f"{entry['realization']}")

    # Read virial radius of host from primaryHalo XML.
    primary_root  = ET.parse(primary_file).getroot()
    virial_radius = primary_root.get('r')

    # Append position to the importer property list.
    importer = parameters.find('nbodyImporter')
    props    = importer.find('properties').get('value', '')
    importer.find('properties').set('value', props + ' position')

    # Extend the delete-properties operator to also remove distanceFromPoint.
    nb_ops    = parameters.find('nbodyOperator').findall('nbodyOperator')
    del_props = nb_ops[1].find('propertyNames').get('value', '')
    nb_ops[1].find('propertyNames').set('value', del_props + ' distanceFromPoint')

    # Splice two operators at index 1: distanceFromPoint and filterProperties.
    parent  = parameters.find('nbodyOperator')
    dist_op = ET.Element('nbodyOperator', value='distanceFromPoint')
    ET.SubElement(dist_op, 'point', value=entry[af_label])

    filt_op = ET.Element('nbodyOperator', value='filterProperties')
    ET.SubElement(filt_op, 'propertyNames', value='distanceFromPoint')
    ET.SubElement(filt_op, 'rangeLow',  value='0.0')
    ET.SubElement(filt_op, 'rangeHigh', value=str(virial_radius))

    parent.insert(1, dist_op)
    parent.insert(2, filt_op)


def symphony_process_subhalo_functions(entry, expansion_factor, parameters, options):
    """Set host position, mass, and mass/radius/Vmax bin ranges for subhalo functions."""
    redshift    = 1.0 / expansion_factor - 1.0
    rl          = f'z{redshift:.3f}'
    af_label    = f'sphericalOrigin:a{expansion_factor:.3f}'

    # Scan all realizations to find the maximum primary halo mass (for range limits).
    mass_primary_maximum = 0.0
    all_realizations     = entry['simulation'].get('realizations', {}).get('value', 'only').split()
    for realization in all_realizations:
        path         = entry['path'].replace(f'/{entry["realization"]}/', f'/{realization}/')
        primary_file = path + f'primaryHalo_{rl}.xml'
        if not os.path.exists(primary_file):
            continue
        primary_root = ET.parse(primary_file).getroot()
        m            = float(primary_root.get('m'))
        if m > mass_primary_maximum:
            mass_primary_maximum = m

    # Read properties of the primary halo for this realization.
    primary_root  = ET.parse(entry['path'] + f'primaryHalo_{rl}.xml').getroot()
    mass_primary  = float(primary_root.get('m'))
    virial_radius = float(primary_root.get('r'))

    # Particle count thresholds.
    mass_count_particles_minimum      = 300
    structure_count_particles_minimum = 2000
    mass_particle                     = entry['resolution']['massParticle']

    nb_ops = parameters.find('nbodyOperator').findall('nbodyOperator')

    # Set host position (operator [1]).
    nb_ops[1].find('point').set('value', entry[af_label])

    # Subhalo mass function (operator [2]).
    mf_count_per_decade = 3
    mass_ratio_minimum  = 10.0 ** (
        (math.floor(
            math.log10(mass_count_particles_minimum * mass_particle / mass_primary_maximum)
            * mf_count_per_decade
        ) + 1.5) / mf_count_per_decade
    )
    nb_ops[2].find('massHost'        ).set('value', str(mass_primary))
    nb_ops[2].find('massRatioMinimum').set('value', str(mass_ratio_minimum))
    nb_ops[2].find('massRatioMaximum').set('value', '1.0')
    nb_ops[2].find('massCountPerDecade').set('value', str(mf_count_per_decade))

    # Subhalo radial function (operator [3]).
    nb_ops[3].find('radiusVirialHost').set('value', str(virial_radius))
    nb_ops[3].find('massMinimum'     ).set('value', str(mass_count_particles_minimum * mass_particle))

    # Subhalo Vmax function (operator [4]).
    vf_count_per_decade = 5
    mass_minimum_vmax   = 10.0 ** (
        (math.floor(
            math.log10(structure_count_particles_minimum * mass_particle)
            * vf_count_per_decade
        ) + 1.5) / vf_count_per_decade
    )
    mass_maximum_vmax = 10.0 ** (math.floor(math.log10(mass_primary_maximum)) + 1.0)
    nb_ops[4].find('massMinimum'      ).set('value', str(mass_minimum_vmax))
    nb_ops[4].find('massMaximum'      ).set('value', str(mass_maximum_vmax))
    nb_ops[4].find('massCountPerDecade').set('value', str(vf_count_per_decade))


def symphony_postprocess_select_in_sphere(entry, jobs, options):
    """Select all N-body particles within the uncontaminated sphere of interest."""
    omp_threads       = options.get('ompThreadsResolved', 1)
    expansion_factors = entry['resolution'].get('expansionFactors', [])
    snapshots         = entry['resolution'].get('snapshots', [])

    # Parse the template XML once and mutate it per epoch.
    root = _parse_param_xml(options['pipelinePath'] + 'zoomInSelectInSphere.xml')
    nb_ops = root.find('nbodyOperator').findall('nbodyOperator')

    for epoch in entry['resolution']['epochs']:
        rl       = epoch['redshiftLabel']
        a        = epoch['expansionFactor']
        af_label = f'sphericalOrigin:a{a:.3f}'

        snapshot = next(
            (snapshots[i] for i, af in enumerate(expansion_factors) if af == a),
            None,
        )
        output_file = entry['path'] + f'selectedParticles_{rl}.hdf5'

        _find_or_create(root, 'outputFileName').set(
            'value', entry['path'] + f'selectInSphereGLC_{rl}.hdf5')
        root.find('nbodyImporter/fileName').set(
            'value', entry['path'] + f'snapshots/snapshot_{snapshot}')
        nb_ops[0].find('point'   ).set('value', entry[af_label])
        nb_ops[1].find('rangeLow' ).set('value', str(entry[rl]['sphericalRadiusMinimum']))
        nb_ops[1].find('rangeHigh').set('value', str(entry[rl]['sphericalRadiusMaximum']))
        nb_ops[4].find('fileName').set('value', output_file)
        nb_ops[4].find('redshift').set('value', str(epoch['redshift']))

        param_file = entry['path'] + f'selectInSphere_{rl}.xml'
        _write_param_xml(root, param_file)

        if os.path.exists(output_file):
            continue

        mem = '16G'
        if entry['group'     ]['name'] == 'Group':
            mem = '32G'
        if entry['resolution']['name'] == 'resolutionX64':
            mem = '196G'
        jobs.append({
            'command':    os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe ' + param_file,
            'launchFile': entry['path'] + f'selectInSphere_{rl}.sh',
            'logFile':    entry['path'] + f'selectInSphere_{rl}.log',
            'label':                      f'selectInSphere_{rl}',
            'ppn':        omp_threads,
            'ompThreads': omp_threads,
            'nodes':      1,
            'mem':        mem,
            'walltime':   '8:00:00',
            'mpi':        'no',
        })


def symphony_postprocess_select_in_ics(entry, jobs, options):
    """Identify the selected particles in the initial conditions snapshot."""
    omp_threads = options.get('ompThreadsResolved', 1)

    # Read the IC redshift from the MUSIC config file.
    music_file = entry['path'] + 'music.conf'
    if not os.path.exists(music_file):
        raise FileNotFoundError(f'missing file: {music_file}')
    with open(music_file) as fh:
        for line in fh:
            m = re.match(r'^zstart\s*=\s*([\d.]+)', line)
            if m:
                entry['redshiftICs'] = float(m.group(1))

    # Parse the template XML once and mutate it per epoch.
    root   = _parse_param_xml(options['pipelinePath'] + 'zoomInSelectInICs.xml')
    nb_ops = root.find('nbodyOperator').findall('nbodyOperator')

    for epoch in entry['resolution']['epochs']:
        rl          = epoch['redshiftLabel']
        output_file = entry['path'] + f'selectedParticles_{rl}_ICs.hdf5'

        root.find('nbodyImporter/fileName').set(
            'value', entry['path'] + 'ic/ic_gadget_dist')
        nb_ops[0].find('idSelectionFileName').set(
            'value', entry['path'] + f'selectedParticles_{rl}.hdf5')
        nb_ops[1].find('fileName').set('value', output_file)
        nb_ops[1].find('redshift').set('value', f'{entry["redshiftICs"]:.3f}')
        _find_or_create(root, 'outputFileName').set(
            'value', entry['path'] + f'selectInICs_{rl}.hdf5')

        param_file = entry['path'] + f'selectInICs_{rl}.xml'
        _write_param_xml(root, param_file)

        if os.path.exists(output_file):
            continue

        mem = '16G'
        if entry['group']['name'] == 'Group':
            mem = '32G'
        if entry['group']['name'] == 'MilkyWay' and entry['resolution']['name'] == 'resolutionX8':
            mem = '32G'
        if entry['group']['name'] == 'MilkyWay' and entry['resolution']['name'] == 'resolutionX64':
            mem = '196G'
        jobs.append({
            'command':    os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe ' + param_file,
            'launchFile': entry['path'] + f'selectInICs_{rl}.sh',
            'logFile':    entry['path'] + f'selectInICs_{rl}.log',
            'label':                      f'selectInICs_{rl}',
            'ppn':        omp_threads,
            'ompThreads': omp_threads,
            'nodes':      1,
            'mem':        mem,
            'walltime':   '8:00:00',
            'mpi':        'no',
        })


def symphony_postprocess_analyze(entry, jobs, options):
    """Measure mass and overdensity of the selected Lagrangian region."""
    omp_threads = options.get('ompThreadsResolved', 1)

    # Parse the template XML once and mutate it per epoch.
    root      = _parse_param_xml(options['pipelinePath'] + 'zoomInAnalyze.xml')
    nb_ops    = root.find('nbodyOperator').findall('nbodyOperator')
    importers = root.find('nbodyImporter').findall('nbodyImporter')

    for epoch in entry['resolution']['epochs']:
        rl = epoch['redshiftLabel']

        _find_or_create(root, 'outputFileName').set(
            'value', entry['path'] + f'environment_{rl}.hdf5')
        importers[0].find('fileName').set(
            'value', entry['path'] + f'selectedParticles_{rl}_ICs.hdf5')
        importers[1].find('fileName').set(
            'value', entry['path'] + 'ic/ic_gadget_dist')
        nb_ops[0].find('values').set('value', f'{entry["redshiftICs"]:.3f}')
        # Operator [5] contains a nested nbodyOperator element.
        nb_ops[5].find('nbodyOperator/fileName').set(
            'value', entry['path'] + f'allParticles_{rl}_ICs.hdf5')
        nb_ops[5].find('nbodyOperator/redshift').set(
            'value', f'{entry["redshiftICs"]:.3f}')

        param_file = entry['path'] + f'analyze_{rl}.xml'
        _write_param_xml(root, param_file)

        if os.path.exists(entry['path'] + f'environment_{rl}:MPI0000.hdf5'):
            continue

        mem = '48G'
        if entry['resolution']['name'] == 'resolutionX8':
            mem = '64G'
        if entry['resolution']['name'] == 'resolutionX64':
            mem = '256G'
        jobs.append({
            'command':    os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe ' + param_file,
            'launchFile': entry['path'] + f'analyze_{rl}.sh',
            'logFile':    entry['path'] + f'analyze_{rl}.log',
            'label':                      f'analyze_{rl}',
            'ppn':        omp_threads,
            'ompThreads': omp_threads,
            'nodes':      1,
            'mem':        mem,
            'walltime':   '8:00:00',
            'mpi':        'no',
        })


def symphony_postprocess_set_volume(entry, jobs, options):
    """Set the effective Lagrangian box size for each epoch."""
    G = 4.3011827419096073e-9
    omega_matter   = entry['suite']['cosmology']['OmegaMatter']
    hubble_const   = entry['suite']['cosmology']['HubbleConstant']
    density_mean   = 3.0 * omega_matter * hubble_const**2 / (8.0 * math.pi * G)

    for epoch in entry['resolution']['epochs']:
        rl = epoch['redshiftLabel']

        selected_file = entry['path'] + f'selectedParticles_{rl}.hdf5'
        with h5py.File(selected_file, 'r') as hdf:
            halo_cat = hdf['Snapshot00001/HaloCatalog']
            if 'massTotal' not in halo_cat.attrs:
                raise KeyError(f"'massTotal' attribute missing in '{selected_file}'")
            mass = float(halo_cat.attrs['massTotal'])

        primary_file = entry['path'] + f'primaryHalo_{rl}.xml'
        primary_root = ET.parse(primary_file).getroot()
        mass_primary = float(primary_root.get('m'))

        if mass_primary >= mass:
            raise ValueError(f"Primary halo mass ({mass_primary}) exceeds region mass ({mass})")

        box_size     = ((mass - mass_primary) / density_mean) ** (1.0 / 3.0)
        radius_region = (3.0 * mass / density_mean / (4.0 * math.pi)) ** (1.0 / 3.0)

        entry.setdefault(rl, {})['radiusRegion'] = radius_region
        entry[rl]['massRegion']  = mass
        entry[rl]['massPrimary'] = mass_primary

        halos_file = entry['path'] + f'nonFlyby_{rl}_subVolume0_0_0.hdf5'
        with h5py.File(halos_file, 'r+') as hdf:
            hdf['Snapshot00001/HaloCatalog'].attrs['boxSize']   = box_size
            hdf['SimulationProperties'].attrs['boxSize']         = box_size


def symphony_postprocess_mass_function(entry, jobs, options):
    """Store region mass, radius, and overdensity into the halo mass function file."""
    for epoch in entry['resolution']['epochs']:
        rl = epoch['redshiftLabel']

        env_file = entry['path'] + f'environment_{rl}:MPI0000.hdf5'
        with h5py.File(env_file, 'r') as hdf:
            if 'simulation0002' not in hdf:
                raise KeyError(f"'simulation0002' group missing in '{env_file}'")
            sim_grp = hdf['simulation0002']
            for attr in ('massTotal', 'convexHullOverdensity'):
                if attr not in sim_grp.attrs:
                    raise KeyError(f"'{attr}' attribute missing in '{env_file}'")
            mass       = float(sim_grp.attrs['massTotal'])
            overdensity = float(sim_grp.attrs['convexHullOverdensity'])

        hmf_file = entry['path'] + f'haloMassFunction_{rl}:MPI0000.hdf5'
        with h5py.File(hmf_file, 'r+') as hdf:
            sim_grp = hdf.require_group('simulation0001')
            sim_grp.attrs['massPrimary']            = entry[rl]['massPrimary']
            sim_grp.attrs['massRegion']             = entry[rl]['massRegion']
            sim_grp.attrs['radiusRegion']           = entry[rl]['radiusRegion']
            sim_grp.attrs['overdensityEnvironment'] = overdensity
            sim_grp.attrs['massEnvironment']        = mass


# ---------------------------------------------------------------------------
# Suite configuration
# ---------------------------------------------------------------------------

SUITES = {
    'MDPL': {
        'analyses': ['haloMassFunction'],
        'steps':    {},
    },
    'Symphony': {
        'analyses': ['haloMassFunction', 'subhaloStatistics'],
        'steps': {
            'identifyAlwaysIsolated': {
                'processParameters': symphony_process_identify_always_isolated,
            },
            'extractHalos': {
                'preprocess': [
                    symphony_preprocess_extract_halos_locate,
                    symphony_preprocess_extract_halos_uncontaminated,
                ],
                'processParameters': symphony_process_extract_halos,
                'postprocess': [
                    symphony_postprocess_select_in_sphere,
                    symphony_postprocess_select_in_ics,
                    symphony_postprocess_analyze,
                    symphony_postprocess_set_volume,
                ],
            },
            'extractSubhalos': {
                'preprocess': [
                    symphony_preprocess_extract_halos_locate,
                    symphony_preprocess_extract_halos_uncontaminated,
                ],
                'processParameters': symphony_process_extract_subhalos,
            },
            'massFunctions': {
                'postprocess': [symphony_postprocess_mass_function],
            },
            'subhaloFunctions': {
                'processParameters': symphony_process_subhalo_functions,
            },
        },
    },
}
# COZMIC uses exactly the same processing pipeline as Symphony.
SUITES['COZMIC'] = SUITES['Symphony']


# ---------------------------------------------------------------------------
# Step functions
# ---------------------------------------------------------------------------

def step_identify_always_isolated(entries, suites_cfg, active_steps, manager, options, omp_threads):
    """Identify always-isolated halos for each simulation subvolume."""
    active_analyses = set(options['analyses'].split(','))
    jobs = []
    for entry in entries:
        suite_name   = entry['suite']['name']
        suite_cfg    = suites_cfg.get(suite_name, {})
        suite_analyses = suite_cfg.get('analyses', [])
        if not any(a in active_analyses for a in suite_analyses):
            continue

        # Build epoch list (expansionFactor, redshift, redshiftLabel) for this resolution.
        entry['resolution']['epochs'] = [
            {
                'expansionFactor': a,
                'redshift':        z,
                'redshiftLabel':   f'z{z:.3f}',
            }
            for a, z in zip(
                entry['resolution']['expansionFactors'],
                entry['resolution']['redshifts'],
            )
        ]

        # Halo mass limits: minimum is two decades above particle mass.
        entry['massMinimum'] = 10.0 ** (int(math.log10(entry['resolution']['massParticle'])) + 2)
        entry['massMaximum'] = 1.0e16

        # Canonical path for this simulation / realization.
        entry['path'] = (
            options['simulationDataPath']
            + suite_name                      + '/'
            + entry['group'     ]['name']     + '/'
            + entry['resolution']['name']     + '/'
            + entry['simulation']['name']     + '/'
            + entry['realization']            + '/'
        )

        n_sub = int(entry['resolution']['subvolumes'])
        proc  = suite_cfg.get('steps', {}).get('identifyAlwaysIsolated', {}).get('processParameters')

        for i in range(n_sub):
            for j in range(n_sub):
                for k in range(n_sub):
                    if not os.path.exists(entry['path'] + f'tree_{i}_{j}_{k}.dat'):
                        continue

                    root = _parse_param_xml(options['pipelinePath'] + 'identifyAlwaysIsolated.xml')
                    _find_or_create(root, 'outputFileName').set(
                        'value', entry['path'] + f'identifyAlwaysIsolatedGLC_{i}_{j}_{k}.hdf5')
                    root.find('nbodyImporter/fileName').set(
                        'value', entry['path'] + f'tree_{i}_{j}_{k}.dat')
                    nb_ops = root.find('nbodyOperator').findall('nbodyOperator')
                    nb_ops[4].find('fileName').set(
                        'value', entry['path'] + f'alwaysIsolated_subVolume{i}_{j}_{k}.hdf5')
                    cosmo = root.find('cosmologyParameters')
                    for key in ('HubbleConstant', 'OmegaMatter', 'OmegaDarkEnergy', 'OmegaBaryon'):
                        cosmo.find(key).set('value', str(entry['suite']['cosmology'][key]))
                    if proc:
                        proc(entry, None, root, options)

                    param_file = entry['path'] + f'identifyAlwaysIsolated_{i}_{j}_{k}.xml'
                    _write_param_xml(root, param_file)

                    if os.path.exists(entry['path'] + f'alwaysIsolated_subVolume{i}_{j}_{k}.hdf5'):
                        continue

                    mem = '64G' if entry['resolution']['name'] == 'resolutionX64' else '8G'
                    jobs.append({
                        'command':    os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe ' + param_file,
                        'launchFile': entry['path'] + f'identifyAlwaysIsolated_{i}_{j}_{k}.sh',
                        'logFile':    entry['path'] + f'identifyAlwaysIsolated_{i}_{j}_{k}.log',
                        'label':                      f'identifyAlwaysIsolated_{i}_{j}_{k}',
                        'ppn':        omp_threads,
                        'ompThreads': omp_threads,
                        'nodes':      1,
                        'mem':        mem,
                        'walltime':   '8:00:00',
                        'mpi':        'no',
                    })

    _submit_jobs(manager, jobs)


def step_extract_halos(entries, suites_cfg, active_steps, manager, options, omp_threads):
    """Extract non-flyby halo snapshots at each expansion factor."""
    jobs = []
    for entry in entries:
        suite_name = entry['suite']['name']
        n_sub      = int(entry['resolution']['subvolumes'])
        proc       = suites_cfg.get(suite_name, {}).get('steps', {}).get('extractHalos', {}).get('processParameters')

        for i in range(n_sub):
            for j in range(n_sub):
                for k in range(n_sub):
                    if not os.path.exists(entry['path'] + f'tree_{i}_{j}_{k}.dat'):
                        continue

                    for epoch in entry['resolution']['epochs']:
                        a       = epoch['expansionFactor']
                        af_low  = (1.0 - 5.0e-4) * a
                        af_high = (1.0 + 5.0e-4) * a
                        rl      = epoch['redshiftLabel']

                        root = _parse_param_xml(options['pipelinePath'] + 'extractHalosSnapshot.xml')
                        _find_or_create(root, 'outputFileName').set(
                            'value', entry['path'] + f'alwaysIsolated_subVolumeGLC{i}_{j}_{k}.hdf5')
                        root.find('nbodyImporter/fileName').set(
                            'value', entry['path'] + f'alwaysIsolated_subVolume{i}_{j}_{k}.hdf5')
                        root.find('nbodyImporter/properties').set(
                            'value', 'particleID isFlyby expansionFactor massVirial hostedRootID')
                        nb_ops = root.find('nbodyOperator').findall('nbodyOperator')
                        nb_ops[0].find('propertyNames').set('value', 'isFlyby expansionFactor')
                        nb_ops[0].find('rangeLow' ).set('value', f'0 {af_low}')
                        nb_ops[0].find('rangeHigh').set('value', f'0 {af_high}')
                        nb_ops[1].find('propertyNames').set('value', 'isFlyby expansionFactor')
                        nb_ops[2].find('fileName').set(
                            'value', entry['path'] + f'nonFlyby_{rl}_subVolume{i}_{j}_{k}.hdf5')
                        nb_ops[2].find('redshift').set('value', str(epoch['redshift']))
                        if proc:
                            proc(entry, a, root, options)

                        param_file = entry['path'] + f'identifyNonFlyby_{rl}_{i}_{j}_{k}.xml'
                        _write_param_xml(root, param_file)

                        if os.path.exists(entry['path'] + f'nonFlyby_{rl}_subVolume{i}_{j}_{k}.hdf5'):
                            continue

                        mem = '64G' if entry['resolution']['name'] == 'resolutionX64' else '8G'
                        jobs.append({
                            'command':    os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe ' + param_file,
                            'launchFile': entry['path'] + f'identifyNonFlyby_{rl}_{i}_{j}_{k}.sh',
                            'logFile':    entry['path'] + f'identifyNonFlyby_{rl}_{i}_{j}_{k}.log',
                            'label':                      f'identifyNonFlyby_{rl}_{i}_{j}_{k}',
                            'ppn':        1,
                            'ompThreads': 1,
                            'nodes':      1,
                            'mem':        mem,
                            'walltime':   '8:00:00',
                            'mpi':        'no',
                        })

    _submit_jobs(manager, jobs)


def step_extract_subhalos(entries, suites_cfg, active_steps, manager, options, omp_threads):
    """Extract subhalo snapshots (halos with a parent) at each expansion factor."""
    jobs = []
    for entry in entries:
        suite_name = entry['suite']['name']
        n_sub      = int(entry['resolution']['subvolumes'])
        proc       = suites_cfg.get(suite_name, {}).get('steps', {}).get('extractSubhalos', {}).get('processParameters')

        for i in range(n_sub):
            for j in range(n_sub):
                for k in range(n_sub):
                    if not os.path.exists(entry['path'] + f'tree_{i}_{j}_{k}.dat'):
                        continue

                    for epoch in entry['resolution']['epochs']:
                        a       = epoch['expansionFactor']
                        af_low  = (1.0 - 5.0e-4) * a
                        af_high = (1.0 + 5.0e-4) * a
                        rl      = epoch['redshiftLabel']

                        root = _parse_param_xml(options['pipelinePath'] + 'extractSubhalosSnapshot.xml')
                        _find_or_create(root, 'outputFileName').set(
                            'value', entry['path'] + f'subhalos_subVolumeGLC{i}_{j}_{k}.hdf5')
                        root.find('nbodyImporter/fileName').set(
                            'value', entry['path'] + f'alwaysIsolated_subVolume{i}_{j}_{k}.hdf5')
                        root.find('nbodyImporter/properties').set(
                            'value', 'particleID isFlyby expansionFactor massVirial velocityMaximum')
                        nb_ops = root.find('nbodyOperator').findall('nbodyOperator')
                        nb_ops[0].find('propertyNames').set('value', 'isFlyby expansionFactor')
                        nb_ops[0].find('rangeLow' ).set('value', f'1 {af_low}')
                        nb_ops[0].find('rangeHigh').set('value', f'1 {af_high}')
                        nb_ops[1].find('propertyNames').set('value', 'isFlyby expansionFactor')
                        nb_ops[2].find('fileName').set(
                            'value', entry['path'] + f'subhalos_{rl}_subVolume{i}_{j}_{k}.hdf5')
                        nb_ops[2].find('redshift').set('value', str(epoch['redshift']))
                        if proc:
                            proc(entry, a, root, options)

                        param_file = entry['path'] + f'identifySubhalos_{rl}_{i}_{j}_{k}.xml'
                        _write_param_xml(root, param_file)

                        if os.path.exists(entry['path'] + f'subhalos_{rl}_subVolume{i}_{j}_{k}.hdf5'):
                            continue

                        mem = '64G' if entry['resolution']['name'] == 'resolutionX64' else '8G'
                        jobs.append({
                            'command':    os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe ' + param_file,
                            'launchFile': entry['path'] + f'identifySubhalos_{rl}_{i}_{j}_{k}.sh',
                            'logFile':    entry['path'] + f'identifySubhalos_{rl}_{i}_{j}_{k}.log',
                            'label':                      f'identifySubhalos_{rl}_{i}_{j}_{k}',
                            'ppn':        1,
                            'ompThreads': 1,
                            'nodes':      1,
                            'mem':        mem,
                            'walltime':   '8:00:00',
                            'mpi':        'no',
                        })

    _submit_jobs(manager, jobs)


def step_mass_functions(entries, suites_cfg, active_steps, manager, options, omp_threads):
    """Compute halo mass functions by combining all subvolumes for each epoch."""
    jobs = []
    for entry in entries:
        suite_name = entry['suite']['name']
        n_sub      = int(entry['resolution']['subvolumes'])
        proc       = suites_cfg.get(suite_name, {}).get('steps', {}).get('massFunctions', {}).get('processParameters')

        for epoch in entry['resolution']['epochs']:
            rl = epoch['redshiftLabel']

            # Collect importers for each subvolume that has a tree file.
            importers = [
                entry['path'] + f'nonFlyby_{rl}_subVolume{i}_{j}_{k}.hdf5'
                for i in range(n_sub)
                for j in range(n_sub)
                for k in range(n_sub)
                if os.path.exists(entry['path'] + f'tree_{i}_{j}_{k}.dat')
            ]

            if os.path.exists(entry['path'] + f'haloMassFunction_{rl}:MPI0000.hdf5'):
                continue

            root = _parse_param_xml(options['pipelinePath'] + 'haloMassFunctionCompute.xml')

            # Replace the merge importer's children with one child per subvolume.
            merge_el = root.find('nbodyImporter')
            for child in merge_el.findall('nbodyImporter'):
                merge_el.remove(child)
            for fname in importers:
                imp_el = ET.SubElement(merge_el, 'nbodyImporter', value='IRATE')
                ET.SubElement(imp_el, 'fileName',   value=fname)
                ET.SubElement(imp_el, 'snapshot',   value='1')
                ET.SubElement(imp_el, 'properties', value='massVirial')

            _find_or_create(root, 'outputFileName').set(
                'value', entry['path'] + f'haloMassFunction_{rl}.hdf5')
            nb_ops = root.find('nbodyOperator').findall('nbodyOperator')
            nb_ops[0].find('values').set('value', str(entry['resolution']['massParticle']))
            nb_ops[1].find('simulationReference').set('value', entry['group']['metaData']['reference'])
            nb_ops[1].find('simulationURL'      ).set('value', entry['group']['metaData']['url'      ])
            nb_ops[1].find('massMinimum'        ).set('value', str(entry['massMinimum']))
            nb_ops[1].find('massMaximum'        ).set('value', str(entry['massMaximum']))
            nb_ops[1].find('description'        ).set('value',
                f'Halo mass function of non-flyby halos for the '
                f'{entry["suite"]["name"]} {entry["group"]["name"]} '
                f'{entry["simulation"]["name"]} {rl} simulation')
            if proc:
                proc(entry, epoch['expansionFactor'], root, options)

            param_file = entry['path'] + f'haloMassFunction_{rl}.xml'
            _write_param_xml(root, param_file)
            jobs.append({
                'command':    os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe ' + param_file,
                'launchFile': entry['path'] + f'haloMassFunction_{rl}.sh',
                'logFile':    entry['path'] + f'haloMassFunction_{rl}.log',
                'label':                      f'haloMassFunction_{rl}',
                'ppn':        omp_threads,
                'ompThreads': omp_threads,
                'nodes':      1,
                'mem':        '8G',
                'walltime':   '8:00:00',
                'mpi':        'no',
            })

    _submit_jobs(manager, jobs)


def step_subhalo_functions(entries, suites_cfg, active_steps, manager, options, omp_threads):
    """Compute subhalo mass, radial, and velocity-maximum functions for each epoch."""
    jobs = []
    for entry in entries:
        suite_name = entry['suite']['name']
        n_sub      = int(entry['resolution']['subvolumes'])
        proc       = suites_cfg.get(suite_name, {}).get('steps', {}).get('subhaloFunctions', {}).get('processParameters')

        for epoch in entry['resolution']['epochs']:
            rl = epoch['redshiftLabel']

            # Collect importers for each subvolume that has a tree file.
            importers = [
                entry['path'] + f'subhalos_{rl}_subVolume{i}_{j}_{k}.hdf5'
                for i in range(n_sub)
                for j in range(n_sub)
                for k in range(n_sub)
                if os.path.exists(entry['path'] + f'tree_{i}_{j}_{k}.dat')
            ]

            if os.path.exists(entry['path'] + f'subhaloFunctions_{rl}:MPI0000.hdf5'):
                continue

            root = _parse_param_xml(options['pipelinePath'] + 'subhaloFunctionsCompute.xml')

            # Replace the merge importer's children with one child per subvolume.
            merge_el = root.find('nbodyImporter')
            for child in merge_el.findall('nbodyImporter'):
                merge_el.remove(child)
            for fname in importers:
                imp_el = ET.SubElement(merge_el, 'nbodyImporter', value='IRATE')
                ET.SubElement(imp_el, 'fileName',   value=fname)
                ET.SubElement(imp_el, 'snapshot',   value='1')
                ET.SubElement(imp_el, 'properties', value='position massVirial velocityMaximum')

            _find_or_create(root, 'outputFileName').set(
                'value', entry['path'] + f'subhaloFunctions_{rl}.hdf5')
            nb_ops = root.find('nbodyOperator').findall('nbodyOperator')
            nb_ops[0].find('values').set('value', str(entry['resolution']['massParticle']))
            ref = entry['group']['metaData']['reference']
            url = entry['group']['metaData']['url']
            suite_grp_sim = (f'{entry["suite"]["name"]} {entry["group"]["name"]} '
                             f'{entry["simulation"]["name"]} {rl} simulation')
            nb_ops[2].find('simulationReference').set('value', ref)
            nb_ops[2].find('simulationURL'      ).set('value', url)
            nb_ops[2].find('description'        ).set('value',
                f'Subhalo mass function of non-flyby halos for the {suite_grp_sim}')
            nb_ops[3].find('simulationReference').set('value', ref)
            nb_ops[3].find('simulationURL'      ).set('value', url)
            nb_ops[3].find('description'        ).set('value',
                f'Subhalo radial distribution function of non-flyby halos for the {suite_grp_sim}')
            nb_ops[4].find('simulationReference').set('value', ref)
            nb_ops[4].find('simulationURL'      ).set('value', url)
            nb_ops[4].find('description'        ).set('value',
                r'Subhalo $V_\mathrm{max}$ function of non-flyby halos for the ' + suite_grp_sim)
            if proc:
                proc(entry, epoch['expansionFactor'], root, options)

            param_file = entry['path'] + f'subhaloFunctions_{rl}.xml'
            _write_param_xml(root, param_file)
            jobs.append({
                'command':    os.environ['GALACTICUS_EXEC_PATH'] + '/Galacticus.exe ' + param_file,
                'launchFile': entry['path'] + f'subhaloFunctions_{rl}.sh',
                'logFile':    entry['path'] + f'subhaloFunctions_{rl}.log',
                'label':                      f'subhaloFunctions_{rl}',
                'ppn':        omp_threads,
                'ompThreads': omp_threads,
                'nodes':      1,
                'mem':        '8G',
                'walltime':   '8:00:00',
                'mpi':        'no',
            })

    _submit_jobs(manager, jobs)


STEP_FUNCTIONS = {
    'identifyAlwaysIsolated': step_identify_always_isolated,
    'extractHalos':           step_extract_halos,
    'extractSubhalos':        step_extract_subhalos,
    'massFunctions':          step_mass_functions,
    'subhaloFunctions':       step_subhalo_functions,
}

STEP_IDS = [
    'identifyAlwaysIsolated',
    'extractHalos',
    'extractSubhalos',
    'massFunctions',
    'subhaloFunctions',
]


# ---------------------------------------------------------------------------
# Results storage stub  (implemented in later steps)
# ---------------------------------------------------------------------------

def store_results(active_analyses, entries, options):
    """Copy halo mass functions and aggregate subhalo statistics to the data repo."""
    data_path = os.environ.get('GALACTICUS_DATA_PATH', '')

    for step_id in active_analyses:

        if step_id == 'haloMassFunction':
            for entry in entries:
                suite_name = entry['suite']['name']
                group_name = entry['group']['name']
                res_name   = entry['resolution']['name']
                sim_name   = entry['simulation']['name']
                realization = entry['realization']
                for epoch in entry['resolution']['epochs']:
                    rl = epoch['redshiftLabel']
                    src = entry['path'] + f'haloMassFunction_{rl}:MPI0000.hdf5'
                    dst = (data_path + '/static/darkMatter/haloMassFunction_'
                           + suite_name + '_' + group_name + '_' + res_name + '_'
                           + sim_name + '_' + realization + '_' + rl + '.hdf5')
                    shutil.copy(src, dst)

        elif step_id == 'subhaloStatistics':
            subhalo_stats = {}

            for entry in entries:
                suite_name  = entry['suite']['name']
                group_name  = entry['group']['name']
                res_name    = entry['resolution']['name']
                sim_name    = entry['simulation']['name']
                reference   = entry['group'].get('metaData', {}).get('reference', '')
                url         = entry['group'].get('metaData', {}).get('url', '')
                # Shorten reference: keep only "(Author YYYY)" portion.
                ref_short = re.sub(r'^([^\)]+)\s+\((\d+).*\)', r'(\1 \2)', reference)

                for epoch in entry['resolution']['epochs']:
                    rl    = epoch['redshiftLabel']
                    label = f'{suite_name}_{group_name}_{res_name}_{sim_name}_{rl}'

                    stat = subhalo_stats.setdefault(label, {
                        'redshift':           epoch['redshift'],
                        'reference':          reference,
                        'referenceURL':       url,
                        'label':              suite_name + ' ' + ref_short,
                        'countRealizations':  0,
                    })
                    stat['countRealizations'] += 1
                    n = stat['countRealizations']

                    func_file = entry['path'] + f'subhaloFunctions_{rl}:MPI0000.hdf5'
                    with h5py.File(func_file, 'r') as hdf:
                        sim_grp = hdf['simulation0001']

                        # --- Subhalo mass function ---
                        smf_grp   = sim_grp['subhaloMassFunction']
                        count     = smf_grp['count'][:]
                        mf        = smf_grp['massFunction'][:]
                        mass_ratio = smf_grp['massRatio'][:]
                        mass_host  = float(smf_grp.attrs['massHost'])
                        if n == 1:
                            stat['smf'] = {
                                'count':      np.zeros_like(mf),
                                'massFunction': np.zeros_like(mf),
                                'massRatio':  mass_ratio.copy(),
                                'massHost':   np.array([]),
                            }
                        else:
                            if not np.array_equal(mass_ratio, stat['smf']['massRatio']):
                                raise ValueError(f'mass ratios do not align for {label}')
                        stat['smf']['count']        += count.astype(float)
                        stat['smf']['massFunction'] += mf
                        stat['smf']['massHost']      = np.append(stat['smf']['massHost'], mass_host)

                        # --- Subhalo radius function ---
                        srf_grp  = sim_grp['subhaloRadiusFunction']
                        count_r  = srf_grp['count'][:]
                        radial   = srf_grp['radialDistribution'][:]
                        rad_ratio = srf_grp['radiusRatio'][:]
                        mass_min  = float(srf_grp.attrs['massMinimum'])
                        if n == 1:
                            stat['srf'] = {
                                'count':             np.zeros_like(radial),
                                'radialDistribution': np.zeros_like(radial),
                                'radiusRatio':       rad_ratio.copy(),
                                'massMinimum':       mass_min,
                            }
                        else:
                            if not np.array_equal(rad_ratio, stat['srf']['radiusRatio']):
                                raise ValueError(f'radius ratios do not align for {label}')
                        stat['srf']['count']             += count_r.astype(float)
                        stat['srf']['radialDistribution'] += radial

                        # --- Subhalo Vmax function ---
                        svf_grp   = sim_grp['subhaloVelocityMaximumMeanFunction']
                        count_v   = svf_grp['count'][:]
                        vmax_mean = svf_grp['velocityMaximumMean'][:]
                        vmax_err  = svf_grp['velocityMaximumMeanError'][:]
                        mass_v    = svf_grp['mass'][:]
                        if n == 1:
                            stat['svf'] = {
                                'count':                    np.zeros_like(vmax_mean),
                                'countRealizations':        np.zeros_like(vmax_mean),
                                'velocityMaximumMean':      np.zeros_like(vmax_mean),
                                'velocityMaximumMeanSquared': np.zeros_like(vmax_mean),
                                'velocityMaximumMeanVariance': np.zeros_like(vmax_mean),
                                'mass':                     mass_v.copy(),
                            }
                        else:
                            if not np.array_equal(mass_v, stat['svf']['mass']):
                                raise ValueError(f'mass arrays do not align for {label}')
                        nonzero = np.where(vmax_mean > 0.0)[0]
                        stat['svf']['count']                     += count_v.astype(float)
                        if nonzero.size > 0:
                            stat['svf']['countRealizations'][nonzero] += 1
                        stat['svf']['velocityMaximumMean']       += vmax_mean
                        stat['svf']['velocityMaximumMeanSquared'] += vmax_mean ** 2
                        stat['svf']['velocityMaximumMeanVariance'] += vmax_err ** 2

            # Generate output for each label.
            for label, stat in subhalo_stats.items():
                n = stat['countRealizations']
                stat['smf']['count']        /= n
                stat['smf']['massFunction'] /= n
                stat['srf']['count']             /= n
                stat['srf']['radialDistribution'] /= n

                svf = stat['svf']
                nonzero_cr = np.where(svf['countRealizations'] > 0)[0]
                if nonzero_cr.size > 0:
                    svf['velocityMaximumMean'][nonzero_cr]       /= svf['countRealizations'][nonzero_cr]
                    svf['velocityMaximumMeanSquared'][nonzero_cr] /= svf['countRealizations'][nonzero_cr]
                    svf['velocityMaximumMeanVariance'][nonzero_cr] /= svf['countRealizations'][nonzero_cr]
                svf['velocityMaximumMeanError'] = np.sqrt(
                    svf['velocityMaximumMeanVariance']
                    + svf['velocityMaximumMeanSquared']
                    - svf['velocityMaximumMean'] ** 2
                )
                nonzero_count = np.where(svf['count'] > 0)[0]
                if nonzero_count.size > 0:
                    svf['velocityMaximumMeanError'][nonzero_count] /= np.sqrt(svf['count'][nonzero_count])

                # Workaround: replicate length-1 massHost arrays.
                mass_host = stat['smf']['massHost']
                if mass_host.size == 1:
                    mass_host = np.append(mass_host, mass_host)
                tree_weights = np.ones_like(mass_host)

                host_file = data_path + f'/static/darkMatter/hostHaloMasses_{label}.hdf5'
                with h5py.File(host_file, 'w') as hdf:
                    hdf.create_dataset('treeRootMass', data=mass_host,    maxshape=(None,))
                    hdf.create_dataset('treeWeight',   data=tree_weights, maxshape=(None,))

                z_str = f'{stat["redshift"]:.3f}'
                store_file = data_path + f'/static/darkMatter/subhaloDistributions_{label}.hdf5'
                with h5py.File(store_file, 'w') as hdf:
                    hdf.attrs['label']        = stat['label']
                    hdf.attrs['redshift']     = stat['redshift']
                    hdf.attrs['reference']    = stat['reference']
                    hdf.attrs['referenceURL'] = stat['referenceURL']
                    hdf.attrs['provenance']   = (
                        'Computed from Symphony Rockstar `tree_?_?_?.dat` files at '
                        + datetime.now().isoformat() + '.'
                    )

                    smf_grp = hdf.create_group('massFunction')
                    smf_grp.attrs['selection'] = f'Subhalos within the host virial radius at z={z_str}.'
                    smf_grp.create_dataset('massRatio',         data=stat['smf']['massRatio'],    maxshape=(None,))
                    smf_grp.create_dataset('massFunction',      data=stat['smf']['count'],        maxshape=(None,))
                    smf_grp.create_dataset('massFunctionError', data=np.sqrt(stat['smf']['count'] / n), maxshape=(None,))

                    srf_grp = hdf.create_group('radialDistribution')
                    srf_grp.attrs['selection']  = f'Subhalos within the host virial radius at z={z_str}.'
                    srf_grp.attrs['massMinimum'] = stat['srf']['massMinimum']
                    srf_grp.create_dataset('radiusFractional',        data=stat['srf']['radiusRatio'],        maxshape=(None,))
                    srf_grp.create_dataset('radialDistribution',      data=stat['srf']['count'],              maxshape=(None,))
                    srf_grp.create_dataset('radialDistributionError', data=np.sqrt(stat['srf']['count'] / n), maxshape=(None,))

                    svf_grp = hdf.create_group('velocityMaximum')
                    svf_grp.attrs['selection'] = f'Subhalos within the host virial radius at z={z_str}.'
                    svf_grp.create_dataset('mass',                     data=svf['mass'],                     maxshape=(None,))
                    svf_grp.create_dataset('velocityMaximumMean',      data=svf['velocityMaximumMean'],      maxshape=(None,))
                    svf_grp.create_dataset('velocityMaximumMeanError', data=svf['velocityMaximumMeanError'], maxshape=(None,))


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Analyze cosmological N-body simulations to extract statistics of interest.'
    )
    parser.add_argument('--simulationDataPath', required=True,
                        help='Path to the simulation data directory')
    parser.add_argument('--pipelinePath', required=True,
                        help='Path to pipeline configuration files')
    parser.add_argument('--analyses', default='haloMassFunction,subhaloStatistics',
                        help='Comma-separated list of analyses to perform')
    parser.add_argument('--ompThreads', default='max',
                        help='Number of OpenMP threads, or "max" to use the queue config value')
    parser.add_argument('--jobMaximum', type=int, default=64,
                        help='Maximum number of concurrent jobs')
    parser.add_argument('--waitOnSubmit', type=int, default=5,
                        help='Seconds to wait between job submissions')
    parser.add_argument('--waitOnActive', type=int, default=30,
                        help='Seconds to wait when polling active jobs')
    parser.add_argument('--partition', default=None,
                        help='SLURM partition override')
    parser.add_argument('--select', default=None,
                        help='Simulation selection filter (suite::group::resolution::...)')
    parser.add_argument('--reOrder', default='no',
                        help='Re-order simulations by power spectrum class ("yes"/"no")')
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = _parse_args()

    # Normalise paths to end with '/'.
    for attr in ('simulationDataPath', 'pipelinePath'):
        val = getattr(args, attr)
        if not val.endswith('/'):
            setattr(args, attr, val + '/')

    # Build the options dict that Simulations.iterate() expects.
    options = vars(args)

    # Determine which analyses and steps are active.
    active_analyses = {a: True for a in args.analyses.split(',')}
    active_steps    = _build_active_steps(active_analyses)

    # Initialise the queue manager.
    manager = queueManager.factory(args)

    # Resolve the number of OMP threads and store in options so symphony hooks can read it.
    if args.ompThreads == 'max':
        omp_threads = manager.options.get('tasksPerNode', 1)
    else:
        omp_threads = int(args.ompThreads)
    options['ompThreadsResolved'] = omp_threads

    # Parse the simulation definition file.
    simulations = parse_simulations_xml(args.pipelinePath + 'simulations.xml')

    # Build the list of simulation entries to process (up to realization level).
    entries = iterate(simulations, options, stop_after='realization')

    # Execute each active step.
    for step_id in STEP_IDS:
        if step_id not in active_steps:
            continue

        # Preprocess hooks.
        _run_hooks('preprocess', step_id, entries, SUITES, manager, options)

        # Main step function.
        STEP_FUNCTIONS[step_id](entries, SUITES, active_steps, manager, options, omp_threads)

        # Postprocess hooks.
        _run_hooks('postprocess', step_id, entries, SUITES, manager, options)

    # Store results to the datasets repository.
    store_results(active_analyses, entries, options)


if __name__ == '__main__':
    main()
