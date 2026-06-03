#!/usr/bin/env python3
# Parse an MCMC config file and generate a LaTeX table describing the priors applied to each parameter.
# Andrew Benson (15-April-2024)

import argparse
import math
import re
import sys
import xml.etree.ElementTree as ET


def _fmt(value, fmt):
    try:
        return fmt % float(value)
    except (ValueError, TypeError):
        return str(value)


def main():
    parser = argparse.ArgumentParser(
        description='Generate a LaTeX table of MCMC parameter priors from a config file.'
    )
    parser.add_argument('configFileName')
    parser.add_argument('--outputFile', default='priors.tex')
    parser.add_argument('--format',     default='%+.2f')
    parser.add_argument('--groupBefore', action='append', default=[], metavar='PARAM LABEL')
    args = parser.parse_args()

    # Parse parameter group labels.
    group_labels = {}
    for opt in args.groupBefore:
        m = re.match(r'^([a-zA-Z0-9:]+)\s+(.*)', opt)
        if not m:
            sys.exit('invalid parameter name')
        group_labels[m.group(1)] = m.group(2)

    # Read the config file.
    tree = ET.parse(args.configFileName)
    root = tree.getroot()
    sim  = root.find('.//posteriorSampleSimulation')

    with open(args.outputFile, 'w') as out:
        out.write('\\begin{tabular}{ll}\n')
        out.write('\\hline\n')
        out.write('\\textbf{Parameter} & \\textbf{Prior} \\\\\n')
        out.write('\\hline\n')
        for i, param_elem in enumerate(sim.findall('modelParameter')):
            name_val = param_elem.find('name').get('value')
            label_el = param_elem.find('label')
            name     = ('$' + label_el.get('value') + '$') if label_el is not None else name_val

            prior_el   = param_elem.find('distributionFunction1DPrior')
            prior_type = prior_el.get('value')

            if prior_type == 'uniform':
                prior_name       = 'Uniform'
                prior_parameters = [
                    float(prior_el.find('limitLower').get('value')),
                    float(prior_el.find('limitUpper').get('value')),
                ]
            elif prior_type == 'logUniform':
                m = re.match(r'^\$(.+)\$$', name)
                if m:
                    name = r'$\log_{10}(' + m.group(1) + r')$'
                else:
                    name = r'$\log_{10}(\hbox{' + name + r'})$'
                prior_name       = 'Uniform'
                prior_parameters = [
                    math.log10(float(prior_el.find('limitLower').get('value'))),
                    math.log10(float(prior_el.find('limitUpper').get('value'))),
                ]
            elif prior_type == 'normal':
                prior_name       = 'Normal'
                prior_parameters = [
                    float(prior_el.find('mean').get('value')),
                    float(prior_el.find('variance').get('value')),
                ]
                lower_el = prior_el.find('limitLower')
                upper_el = prior_el.find('limitUpper')
                if lower_el is not None or upper_el is not None:
                    prior_parameters.append(
                        float(lower_el.get('value')) if lower_el is not None else r'$-\infty$'
                    )
                    prior_parameters.append(
                        float(upper_el.get('value')) if upper_el is not None else r'$+\infty$'
                    )
            elif prior_type == 'logNormal':
                prior_name       = 'Lognormal'
                prior_parameters = [
                    float(prior_el.find('mean').get('value')),
                    float(prior_el.find('variance').get('value')),
                ]
            else:
                sys.exit(f'unknown prior: {prior_type}')

            # Write any parameter group label.
            if name_val in group_labels:
                if i > 0:
                    out.write('\\hline\n')
                out.write('\\multicolumn{2}{c}{\\emph{' + group_labels[name_val] + '}} \\\\\n')

            # Write the parameter prior.
            formatted = ','.join(_fmt(p, args.format) for p in prior_parameters)
            out.write(name + ' & \\hbox{' + prior_name + '}(' + formatted + ') \\\\\n')

        out.write('\\hline\n')
        out.write('\\end{tabular}\n')


if __name__ == '__main__':
    main()
