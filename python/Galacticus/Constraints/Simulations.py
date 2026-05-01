# Provides utilities for iterating over simulation data structures used as constraints.
# Python port of perl/Galacticus/Constraints/Simulations.pm
# Andrew Benson (ported to Python 2026)

import os
import re
import math

import numpy as np
import lxml.etree as ET
import h5py

from List.ExtraUtils import hash_list, as_array
from XML.Utils import xml_to_dict

__all__ = ['select_simulations', 'match_selection', 'iterate', 'parse_simulations_xml']


def select_simulations(options):
    """Parse options['select'] into a list of selection filter dicts.

    Input format: "suite::group::resolution::simulation::realization::redshift"
    where multiple values at each level are comma-separated.
    Returns a list of dicts, each with keys 'suite', 'group', 'resolution',
    'simulation', 'realization', 'redshift' — values are lists of strings or None.
    """
    levels = ('suite', 'group', 'resolution', 'simulation', 'realization', 'redshift')
    selections = []
    for selection_str in as_array(options.get('select')):
        parts = selection_str.split('::')
        sel = {}
        for i, level in enumerate(levels):
            if i < len(parts):
                sel[level] = parts[i].split(',') if parts[i] else None
            else:
                sel[level] = None
        selections.append(sel)
    if 'select' in options and options['select'] is not None and len(selections) == 0:
        raise ValueError("No entries matched the provided selections")
    return selections


def match_selection(selections, suite=None, group=None, resolution=None,
                    simulation=None, realization=None,
                    is_best_resolution=False, redshift=None):
    """Return True if the given hierarchy matches any selection filter.

    An empty selections list matches everything. '*' is a wildcard. For the
    resolution level, the special value 'best' matches only when
    is_best_resolution is True.
    """
    if not selections:
        return True
    for sel in selections:
        # Suite level
        if sel.get('suite') is not None and suite is not None:
            if not any(v == suite or v == '*' for v in sel['suite']):
                continue
        # Group level
        if sel.get('group') is not None and group is not None:
            if not any(v == group or v == '*' for v in sel['group']):
                continue
        # Resolution level — supports 'best' keyword
        if sel.get('resolution') is not None and resolution is not None:
            resolution_matches = False
            for v in sel['resolution']:
                if v == 'best':
                    if is_best_resolution:
                        resolution_matches = True
                elif v == resolution or v == '*':
                    resolution_matches = True
            if not resolution_matches:
                continue
        # Simulation level
        if sel.get('simulation') is not None and simulation is not None:
            if not any(v == simulation or v == '*' for v in sel['simulation']):
                continue
        # Realization level
        if sel.get('realization') is not None and realization is not None:
            if not any(v == realization or v == '*' for v in sel['realization']):
                continue
        # Redshift level
        if sel.get('redshift') is not None and redshift is not None:
            if not any(v == redshift or v == '*' for v in sel['redshift']):
                continue
        return True
    return False


def iterate(simulations, options, stop_after='redshift'):
    """Build and return a list of simulation entry dicts.

    Traverses suite → group → resolution → simulation → realization → redshift,
    applying match_selection() at each level. Lazily loads cosmology and
    simulation metadata from XML sidecars in options['pipelinePath'].

    stop_after controls how deep to iterate:
        'suite', 'group', 'resolution', 'simulation', 'realization', 'redshift'

    Returns a list of dicts. Keys present depend on stop_after level; at the
    deepest level ('redshift') each dict also contains fileTargetData, and
    conditionally massPrimary and environment.
    """
    allowed_stop_after = ('suite', 'group', 'resolution', 'simulation', 'realization', 'redshift')
    if stop_after not in allowed_stop_after:
        raise ValueError(
            f'Unrecognized "stop_after" option "{stop_after}" — '
            f'allowed: {", ".join(allowed_stop_after)}'
        )

    # Build selection list once and cache on the simulations dict.
    if 'selections' not in simulations:
        simulations['selections'] = select_simulations(options)
    selections = simulations['selections']

    pipeline_path = options.get('pipelinePath', '')

    simulation_list = []

    for suite in hash_list(simulations.get('suite', {}), key_as='name'):
        suite_name = suite['name']
        if not match_selection(selections, suite=suite_name):
            continue

        # Lazily load cosmology for this suite.
        if 'cosmology' not in suite:
            _load_cosmology(suite, pipeline_path)

        if stop_after == 'suite':
            simulation_list.append({'suite': suite})
            continue

        for group in hash_list(suite.get('group', {}), key_as='name'):
            group_name = group['name']
            if not match_selection(selections, suite=suite_name, group=group_name):
                continue

            # Lazily load simulation metadata (reference, URL) for this group.
            if 'metaData' not in group:
                _load_group_metadata(suite_name, group_name, group, pipeline_path)

            # Determine the best (highest-numbered) resolution in this group.
            res_numbers = []
            for res in hash_list(group.get('resolution', {}), key_as='name'):
                num = _resolution_number(res['name'])
                if num is not None:
                    res_numbers.append(num)
            best_res_number = max(res_numbers) if res_numbers else None
            group['resolutionBest'] = f'resolutionX{int(best_res_number)}' if best_res_number is not None else None

            if stop_after == 'group':
                simulation_list.append({'suite': suite, 'group': group})
                continue

            for resolution in hash_list(group.get('resolution', {}), key_as='name'):
                res_name = resolution['name']
                if not match_selection(selections, suite=suite_name, group=group_name,
                                       resolution=res_name):
                    continue

                # Lazily load subvolumes, particle mass, and epoch lists.
                if 'subvolumes' not in resolution:
                    _load_resolution_data(suite_name, group_name, res_name, resolution, suite, pipeline_path)

                if stop_after == 'resolution':
                    simulation_list.append({'suite': suite, 'group': group, 'resolution': resolution})
                    continue

                for simulation in hash_list(resolution.get('simulation', {}), key_as='name'):
                    sim_name = simulation['name']
                    if not match_selection(selections, suite=suite_name, group=group_name,
                                           resolution=res_name, simulation=sim_name):
                        continue

                    # Build realization list.
                    if 'realizationsList' not in simulation:
                        if 'realizations' in simulation:
                            simulation['realizationsList'] = simulation['realizations']['value'].split()
                        else:
                            simulation['realizationsList'] = ['only']

                    if stop_after == 'simulation':
                        simulation_list.append({
                            'suite': suite, 'group': group,
                            'resolution': resolution, 'simulation': simulation,
                        })
                        continue

                    res_number = _resolution_number(res_name)

                    for realization in simulation['realizationsList']:
                        # Determine whether this is the best resolution for this realization.
                        is_best = True
                        if res_number is not None:
                            for other_res in hash_list(group.get('resolution', {}), key_as='name'):
                                other_num = _resolution_number(other_res['name'])
                                if other_num is None or other_num <= res_number:
                                    continue
                                for other_sim in hash_list(other_res.get('simulation', {}), key_as='name'):
                                    if other_sim['name'] != sim_name:
                                        continue
                                    if 'realizations' in other_sim:
                                        other_realizations = other_sim['realizations']['value'].split()
                                    else:
                                        other_realizations = ['only']
                                    if realization in other_realizations:
                                        is_best = False

                        if not match_selection(selections, suite=suite_name, group=group_name,
                                               resolution=res_name, simulation=sim_name,
                                               realization=realization,
                                               is_best_resolution=is_best):
                            continue

                        if stop_after == 'realization':
                            simulation_list.append({
                                'suite': suite, 'group': group,
                                'resolution': resolution, 'simulation': simulation,
                                'realization': realization,
                            })
                            continue

                        # Iterate over redshifts.
                        for redshift in resolution.get('redshifts', []):
                            redshift_label = f'{redshift:.3f}'
                            if not match_selection(selections, suite=suite_name, group=group_name,
                                                   resolution=res_name, simulation=sim_name,
                                                   realization=realization,
                                                   is_best_resolution=is_best,
                                                   redshift=redshift_label):
                                continue

                            entry = {
                                'suite':       suite,
                                'group':       group,
                                'resolution':  resolution,
                                'simulation':  simulation,
                                'realization': realization,
                                'redshift':    redshift_label,
                            }

                            # Set target data file path.
                            entry['fileTargetData'] = (
                                os.environ.get('GALACTICUS_DATA_PATH', '') +
                                '/static/darkMatter/haloMassFunction_' +
                                suite_name + '_' +
                                group_name + '_' +
                                res_name + '_' +
                                sim_name + '_' +
                                realization + '_' +
                                'z' + redshift_label +
                                '.hdf5'
                            )

                            if not os.path.exists(entry['fileTargetData']):
                                raise FileNotFoundError(
                                    f"Target data file '{entry['fileTargetData']}' does not exist"
                                )

                            # Extract optional attributes from target data file.
                            with h5py.File(entry['fileTargetData'], 'r') as hdf:
                                sim_grp = hdf['simulation0001']
                                limit_str = suite.get('limitMassMaximum', {}).get('value', '')
                                if 'primaryFraction' in limit_str.split(':'):
                                    entry['massPrimary'] = float(sim_grp.attrs['massPrimary'])
                                if suite.get('includeEnvironment', {}).get('value') == 'true':
                                    entry['environment'] = {}
                                    for attr in ('massEnvironment', 'overdensityEnvironment'):
                                        if attr not in sim_grp.attrs:
                                            raise KeyError(
                                                f"Attribute '{attr}' missing in {entry['fileTargetData']}"
                                            )
                                        entry['environment'][attr] = float(sim_grp.attrs[attr])

                            if stop_after == 'redshift':
                                simulation_list.append(entry)

    # Re-order to spread models with the same power spectrum class for MPI load balancing.
    if options.get('reOrder') == 'yes' and stop_after == 'redshift':
        ps_classes = {}
        for entry in simulation_list:
            ps_class = entry['simulation'].get('powerSpectrumClass', {}).get('value', '')
            ps_classes.setdefault(ps_class, []).append(entry)
        reordered = []
        while ps_classes:
            for ps_class in sorted(list(ps_classes.keys())):
                reordered.append(ps_classes[ps_class].pop())
                if not ps_classes[ps_class]:
                    del ps_classes[ps_class]
        simulation_list = reordered
        for entry in simulation_list:
            print(entry['simulation'].get('powerSpectrumClass', {}).get('value', ''))

    return simulation_list


def parse_simulations_xml(path):
    """Parse simulations.xml into a nested dict matching XML::Simple KeyAttr output.

    suite, group, resolution, and simulation elements are keyed by their
    'name' attribute in the result dict.
    """
    keyed_tags = {'suite', 'group', 'resolution', 'simulation'}
    root = ET.parse(path).getroot()
    result = {}
    suites = root.findall('suite')
    if suites:
        result['suite'] = {
            el.get('name'): xml_to_dict(el, keyed_tags=keyed_tags)
            for el in suites
        }
    return result


# Perl-compatible camelCase aliases for public API.
selectSimulations    = select_simulations
matchSelection       = match_selection
parseSimulationsXml  = parse_simulations_xml


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _resolution_number(name):
    """Extract the numeric suffix from a resolution name like 'resolutionX8'."""
    m = re.match(r'^resolutionX(.+)$', name)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return None
    return None


def _load_cosmology(suite, pipeline_path):
    """Parse cosmology_<suite>.xml and store values in suite['cosmology']."""
    path = os.path.join(pipeline_path, f"cosmology_{suite['name']}.xml")
    tree = ET.parse(path)
    root = tree.getroot()
    cosmo_el = root if root.tag == 'cosmologyParameters' else root.find('.//cosmologyParameters')
    suite['cosmology'] = {}
    for key in ('HubbleConstant', 'OmegaMatter', 'OmegaDarkEnergy', 'OmegaBaryon'):
        el = cosmo_el.find(key) if cosmo_el is not None else root.find(f'.//{key}')
        if el is not None:
            value = el.get('value')
            if value is None:
                raise ValueError(f"cosmology_{suite['name']}.xml: <{key}> element has no 'value' attribute")
            suite['cosmology'][key] = float(value)
    if 'HubbleConstant' in suite['cosmology']:
        suite['cosmology']['hubbleConstant'] = suite['cosmology']['HubbleConstant'] / 100.0


def _load_group_metadata(suite_name, group_name, group, pipeline_path):
    """Parse simulation_<suite>_<group>.xml and store metadata in group['metaData']."""
    path = os.path.join(pipeline_path, f'simulation_{suite_name}_{group_name}.xml')
    tree = ET.parse(path)
    root = tree.getroot()
    sim_el = root if root.tag == 'simulation' else root.find('.//simulation')
    group['metaData'] = {}
    for field, xml_name in (('reference', 'simulationReference'), ('url', 'simulationURL')):
        el = sim_el.find(xml_name) if sim_el is not None else root.find(f'.//{xml_name}')
        if el is not None:
            group['metaData'][field] = el.get('value', '')


def _load_resolution_data(suite_name, group_name, res_name, resolution, suite, pipeline_path):
    """Parse simulation_<suite>_<group>.xml and populate resolution fields."""
    path = os.path.join(pipeline_path, f'simulation_{suite_name}_{group_name}.xml')
    tree = ET.parse(path)
    root = tree.getroot()

    # Subvolumes count.
    sv_el = root.find(f'.//subvolumes/{res_name}')
    if sv_el is not None:
        resolution['subvolumes'] = int(sv_el.get('value'))

    # Particle mass (may contain a Hubble-scaled expression).
    mp_el = root.find(f'.//massParticle/{res_name}')
    if mp_el is not None:
        mass_particle_str = mp_el.get('value', '')
        if mass_particle_str.startswith('='):
            expr = mass_particle_str[1:]
            if 'HubbleConstant' not in suite['cosmology']:
                raise ValueError(
                    f"massParticle/{res_name} expression references "
                    f"[cosmologyParameters/HubbleConstant] but no HubbleConstant "
                    f"is defined in suite '{suite['name']}'"
                )
            hubble = suite['cosmology']['HubbleConstant']
            expr = re.sub(r'\[cosmologyParameters/HubbleConstant\]', str(hubble), expr)
            # eval() is intentional: expression comes from a trusted pipeline XML file.
            resolution['massParticle'] = float(eval(expr))  # noqa: S307
        else:
            resolution['massParticle'] = float(mass_particle_str)

    # Expansion factors and derived redshifts.
    hmf = resolution.get('haloMassFunction', {})
    ef_str = hmf.get('expansionFactors', {}).get('value', '')
    if ef_str:
        a_vals = np.array([float(x) for x in ef_str.split()])
        resolution['expansionFactors'] = a_vals.tolist()
        resolution['redshifts'] = (1.0 / a_vals - 1.0).tolist()

    # Optional snapshot identifiers.
    snap_str = hmf.get('snapshots', {}).get('value', '')
    if snap_str:
        resolution['snapshots'] = snap_str.split()
