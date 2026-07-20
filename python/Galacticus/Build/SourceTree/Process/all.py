"""Import every source-tree Process submodule, registering all process hooks.

Importing a Process submodule registers its hook with the SourceTree
pipeline. Every driver that runs ``process_tree`` must import the
FULL hook set -- a missing import in one driver but not another would
silently generate different code from the same source. Import this module
instead of listing the submodules individually:

    import Galacticus.Build.SourceTree.Process.all  # noqa: F401

When adding a new Process submodule, add it here (and nowhere else).
"""

import Galacticus.Build.SourceTree.Process.AddMetaProperty          # noqa: F401
import Galacticus.Build.SourceTree.Process.Allocate                 # noqa: F401
import Galacticus.Build.SourceTree.Process.ComponentBuilder         # noqa: F401
import Galacticus.Build.SourceTree.Process.ConditionalCall          # noqa: F401
import Galacticus.Build.SourceTree.Process.Constants                # noqa: F401
import Galacticus.Build.SourceTree.Process.Constructors             # noqa: F401
import Galacticus.Build.SourceTree.Process.DebugHDF5                # noqa: F401
import Galacticus.Build.SourceTree.Process.DebugMPI                 # noqa: F401
import Galacticus.Build.SourceTree.Process.DeepCopyActions          # noqa: F401
import Galacticus.Build.SourceTree.Process.DeepCopyFinalize         # noqa: F401
import Galacticus.Build.SourceTree.Process.DeepCopyReset            # noqa: F401
import Galacticus.Build.SourceTree.Process.Dependencies             # noqa: F401
import Galacticus.Build.SourceTree.Process.Enumeration              # noqa: F401
import Galacticus.Build.SourceTree.Process.EventHooks               # noqa: F401
import Galacticus.Build.SourceTree.Process.EventHooksStatic         # noqa: F401
import Galacticus.Build.SourceTree.Process.ForEach                  # noqa: F401
import Galacticus.Build.SourceTree.Process.FunctionClass            # noqa: F401
import Galacticus.Build.SourceTree.Process.FunctionsGlobal          # noqa: F401
import Galacticus.Build.SourceTree.Process.Generics                 # noqa: F401
import Galacticus.Build.SourceTree.Process.HDF5FCInterop            # noqa: F401
import Galacticus.Build.SourceTree.Process.InputParameter           # noqa: F401
import Galacticus.Build.SourceTree.Process.InputParametersValidate  # noqa: F401
import Galacticus.Build.SourceTree.Process.MetaPropertyDatabase     # noqa: F401
import Galacticus.Build.SourceTree.Process.NonProcessed             # noqa: F401
import Galacticus.Build.SourceTree.Process.ObjectBuilder            # noqa: F401
import Galacticus.Build.SourceTree.Process.OptionalArgument         # noqa: F401
import Galacticus.Build.SourceTree.Process.ParameterMigration       # noqa: F401
import Galacticus.Build.SourceTree.Process.ProfileOpenMP            # noqa: F401
import Galacticus.Build.SourceTree.Process.SourceDigest             # noqa: F401
import Galacticus.Build.SourceTree.Process.SourceIntrospection      # noqa: F401
import Galacticus.Build.SourceTree.Process.StateStorable            # noqa: F401
import Galacticus.Build.SourceTree.Process.StateStore               # noqa: F401
import Galacticus.Build.SourceTree.Process.ThreadSafeIO             # noqa: F401
