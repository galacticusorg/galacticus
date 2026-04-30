# Cross-module registry of build-handler hooks.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Hooks.pm: a module-level dict (`module_hooks`,
# `our %moduleHooks` in Perl) that handler modules populate at import time.
# Each entry maps a build `type` (e.g. `"component"`) to a dict of
# `validate` / `parse` / `generate` callables that buildCode.py invokes
# during a build.

module_hooks = {}
