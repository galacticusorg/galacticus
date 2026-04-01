#!/usr/bin/env python3
import argparse
import colorsys
import re
import sys
from datetime import datetime

# Parse profiling information from a build log file.
# Andrew Benson (21-June-2021)

# Generates an HTML representation of the tasks being performed during each second of compilation.

# A color bar is shown next to each task indicating the period over which it was running.

# The color of the bar changes from red to green to indicate the number of simultaneous tasks
# running, such that red indicates regions where build parallelism is minimal.

parser = argparse.ArgumentParser(description="Parse profiling information from a build log file.")
parser.add_argument('buildLogFile',  help="Path to the build log file")
parser.add_argument('profileFile',   help="Path to the output HTML profile file")
parser.add_argument('--durationMinimum', type=float, default=1,
                    help="The shortest duration task (in seconds) to include in the profile")
args = parser.parse_args()


def gradient_color(fraction, gradient):
    """Interpolate a color along a gradient and return a CSS hex string."""
    fracs = gradient['fraction']
    # Linear interpolation between the two control points.
    t = (fraction - fracs[0]) / (fracs[1] - fracs[0]) if fracs[1] != fracs[0] else 0.0
    t = max(0.0, min(1.0, t))
    hue = gradient['hue'][0]        + t * (gradient['hue'][1]        - gradient['hue'][0])
    sat = gradient['saturation'][0] + t * (gradient['saturation'][1] - gradient['saturation'][0])
    val = gradient['value'][0]      + t * (gradient['value'][1]      - gradient['value'][0])
    # colorsys uses hue in [0, 1]; convert from degrees.
    r, g, b = colorsys.hsv_to_rgb(hue / 360.0, sat, val)
    r = int(round(r * 255))
    g = int(round(g * 255))
    b = int(round(b * 255))
    a = 255
    return f"#{r:02x}{g:02x}{b:02x}{a:02x}"


def create_cell(cell, out):
    """Emit JavaScript for a run of accumulated identical cells."""
    if cell['count'] == 0:
        return
    if cell['type'] == 'empty':
        out.write(f"    chartContent += '<td/>'.repeat({cell['count']})\n")
    else:
        out.write(f"    chartContent += '<td style=\"background: {cell['type']}\"/>'.repeat({cell['count']})\n")
    cell['type']  = None
    cell['count'] = 0


# Parse profiling information from the file.
tasks        = []
time_earliest = None
time_latest   = None

with open(args.buildLogFile) as f:
    for line in f:
        m = re.match(r'^\+\+Task: \{([\s\d\+\-:]+)\|([\s\d\+\-:]+)\} \'(.*)', line)
        if not m:
            continue
        start_str = m.group(1).strip()
        stop_str  = m.group(2).strip()
        command   = m.group(3)
        start_time = datetime.fromisoformat(start_str)
        stop_time  = datetime.fromisoformat(stop_str)
        duration   = (stop_time - start_time).total_seconds()
        # Process only tasks above the minimum duration.
        if duration < args.durationMinimum:
            continue
        # Find earliest and latest times in the build.
        if time_earliest is None or start_time < time_earliest:
            time_earliest = start_time
        if time_latest is None or stop_time > time_latest:
            time_latest = stop_time
        # Simplify the command to a more human-readable form where possible.
        elements = command.split()
        if len(elements) > 1 and elements[0] == "./scripts/build/preprocess.pl":
            elements[1] = elements[1].replace("source/", "")
            command = elements[1] + " (preprocess)"
        elif len(elements) > 4 and elements[0] == "gfortran":
            elements[4] = elements[4].replace("./work/build/", "")
            command = elements[4] + " (compile)"
        elif len(elements) > 2 and elements[0] == "./scripts/build/sourceDigests.pl":
            elements[2] = elements[2].replace("./work/build/", "").replace("'", "")
            command = elements[2] + " (source digests)"
        elif len(elements) > 2 and elements[0] == "./scripts/build/parameterDependencies.pl":
            elements[2] = elements[2].replace("./work/build/", "").replace("'", "")
            command = elements[2] + " (parameter dependencies)"
        elif len(elements) > 2 and elements[0] == "./scripts/build/buildCode.pl":
            elements[2] = elements[2].replace("./work/build/", "").replace("'", "")
            command = elements[2] + " (build)"
        elif len(elements) > 8 and elements[1] == "-MRegexp::Common":
            elements[8] = elements[8].replace("work/build/", "")
            command = elements[8] + " (cpp)"
        else:
            if elements:
                mb = re.match(r'\./scripts/build/(.*)\.pl', elements[0])
                if mb:
                    command = mb.group(1)
        tasks.append({
            'description': command,
            'startTime':   start_time,
            'endTime':     stop_time,
        })

if not tasks:
    print("No tasks found in build log.", file=sys.stderr)
    raise SystemExit(1)

# Find the maximum time (in seconds) for the build.
time_maximum = 0
for task in tasks:
    task['start'] = int((task['startTime'] - time_earliest).total_seconds())
    task['end']   = int((task['endTime']   - time_earliest).total_seconds())
    if task['end'] > time_maximum:
        time_maximum = task['end']

# Find the number of threads running during each second of the build.
thread_count_maximum = 0
thread_count = [0] * (time_maximum + 1)
# Use a sweep-line / difference-array approach to compute thread_count in O(tasks + seconds).
# diff[i] represents the change in thread count at time i.
diff = [0] * (time_maximum + 2)
for task in tasks:
    start = task['start']
    end = task['end']
    if 0 <= start <= time_maximum:
        diff[start] += 1
    if end + 1 <= time_maximum:
        diff[end + 1] -= 1
current = 0
for i in range(time_maximum + 1):
    current += diff[i]
    thread_count[i] = current
    if thread_count[i] > thread_count_maximum:
        thread_count_maximum = thread_count[i]

# Compute a cost for each task. This is the sum of the inverse of the number of threads running
# each second.
cost_maximum = 0.0
cost_total   = 0.0
task_maximum = None
for task in tasks:
    task['cost']          = 0.0
    task['isMaximumCost'] = False
    for i in range(task['start'], task['end'] + 1):
        task['cost'] += 1.0 / thread_count[i]
    cost_total += task['cost']
    if task['cost'] > cost_maximum:
        cost_maximum = task['cost']
        task_maximum = task
if task_maximum is not None:
    task_maximum['isMaximumCost'] = True

# Create output.
gradient = {
    'fraction':   [0.0,   1.0],
    'hue':        [0.0, 120.0],
    'saturation': [1.0,   1.0],
    'value':      [1.0,   1.0],
}

with open(args.profileFile, 'w') as out:
    out.write("<html>\n")
    out.write(" <head>\n")
    out.write(" <title>Galacticus Build Profile</title>\n")
    out.write("  <style>\n")
    out.write("   body { font:16px Calibri;}\n")
    out.write("   td, th {\n")
    out.write("           margin: 0;\n")
    out.write("           border: 1px;\n")
    out.write("           height: 1em;\n")
    out.write("           white-space: nowrap;\n")
    out.write("           border-top-width: 0px;\n")
    out.write("           width: 1px;\n")
    out.write("          }\n")
    out.write("   div {\n")
    out.write("        width: 500px;\n")
    out.write("        overflow-x: scroll;\n")
    out.write("        margin-left: 22em;\n")
    out.write("        overflow-y: visible;\n")
    out.write("        padding-bottom: 1px;\n")
    out.write("       }\n")
    out.write("   .headcol {\n")
    out.write("             position: absolute;\n")
    out.write("             width: 0em;\n")
    out.write("             left: 0;\n")
    out.write("             top: auto;\n")
    out.write("             border-top-width: 1px;\n")
    out.write("             margin-top: -1px;\n")
    out.write("            }\n")
    out.write("   .headcolbold {\n")
    out.write("             position: absolute;\n")
    out.write("             width: 0em;\n")
    out.write("             left: 0;\n")
    out.write("             top: auto;\n")
    out.write("             border-top-width: 1px;\n")
    out.write("             margin-top: -1px;\n")
    out.write("             color: red;\n")
    out.write("            }\n")
    out.write("  </style>\n")
    out.write(" </head>\n")

    # Table of build profile.
    out.write("<body>\n")
    out.write(f"Total compile time = {time_maximum} seconds<p>\n")
    out.write(
        f"Build profile (all tasks which ran for {args.durationMinimum} seconds or longer: "
        f"{100.0 * cost_total / time_maximum:5.2f}% of total time)<br>\n"
    )
    out.write("<div>\n")
    out.write("<table id=\"chart\" style=\"border-spacing: 0; border-collapse: separate; border-top: 1px\">\n")
    out.write("</table>\n")
    out.write("</div><p>\n")

    # Table of task costs.
    tasks_ordered = sorted(tasks, key=lambda t: t['cost'], reverse=True)
    out.write("Task costs (ordered)<br>\n")
    out.write("<table>\n")
    for task in tasks_ordered:
        out.write(
            f"<tr><td>{task['description']}</td>"
            f"<td>{task['cost']:7.2f}</td>"
            f"<td>({100.0 * task['cost'] / time_maximum:5.2f}%)</td></tr>\n"
        )
    out.write("</table>\n")
    out.write("  <script type=\"text/javascript\">\n")
    out.write("   function makeChart() {\n")
    out.write("    var chartContent = ''\n")
    for task in tasks:
        cls = "bold" if task['isMaximumCost'] else ""
        description = task['description'].replace("'", "\\'")
        out.write(f"    chartContent += '<tr><th class=\"headcol{cls}\">{description}</th>'\n")
        cell = {'count': 0, 'type': None}
        for i in range(time_maximum + 1):
            if i < task['start']:
                if cell['type'] != 'empty':
                    create_cell(cell, out)
                cell['type'] = 'empty'
                cell['count'] += 1
            elif i > task['end']:
                # Missing columns at the end of the table row are ignored.
                create_cell(cell, out)
            else:
                fraction = thread_count[i] / thread_count_maximum
                color    = gradient_color(fraction, gradient)
                if cell['type'] != color:
                    create_cell(cell, out)
                cell['type'] = color
                cell['count'] += 1
        create_cell(cell, out)
        out.write("     chartContent += '</tr>'\n")
    out.write("    document.getElementById('chart').innerHTML = chartContent\n")
    out.write("   }\n")
    out.write("   window.onLoad = makeChart()\n")
    out.write("  </script>\n")
    out.write("</body>\n")
    out.write("</html>\n")
