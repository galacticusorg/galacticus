.. _manual-sec-editorSetup:

Editor Setup
============

This page describes how to set up an editor for working on Galacticus. The
recommended editor is `Visual Studio Code <https://code.visualstudio.com/>`__,
for which the repository ships a ready-made workspace configuration and a
companion extension. Other editors can be used, but the guidance below on the
embedded languages and on Fortran auto-formatting applies to any editor.

.. _manual-sec-editorVSCode:

Visual Studio Code
------------------

Opening the Galacticus checkout as a workspace folder picks up the configuration
committed under ``.vscode/`` and the language-server configuration in the
repository-root ``.fortls`` file. No manual setup of these is required.

Recommended extensions
~~~~~~~~~~~~~~~~~~~~~~~~

The workspace recommends the extensions below (VSCode offers to install them
when the folder is first opened, via ``.vscode/extensions.json``):

* `Modern Fortran <https://marketplace.visualstudio.com/items?itemName=fortran-lang.linter-gfortran>`__
  (``fortran-lang.linter-gfortran``) — Fortran syntax highlighting, linting, and
  integration with the ``fortls`` language server.
* `Python <https://marketplace.visualstudio.com/items?itemName=ms-python.python>`__
  and `Pylance <https://marketplace.visualstudio.com/items?itemName=ms-python.vscode-pylance>`__
  — for the Python infrastructure under ``python/`` and ``scripts/``.
* `XML <https://marketplace.visualstudio.com/items?itemName=redhat.vscode-xml>`__
  (``redhat.vscode-xml``) — XML editing, plus real-time validation and
  autocompletion of parameter files against the generated
  ``schema/parameters.xsd``. The shipped ``.vscode/settings.json`` wires this up;
  see :ref:`editor-parameter-validation` in the user guide for what the schema
  checks (and what it deliberately does not).
* `YAML <https://marketplace.visualstudio.com/items?itemName=redhat.vscode-yaml>`__
  (``redhat.vscode-yaml``) — for the workflow and configuration files.
* `reStructuredText Syntax highlighting <https://marketplace.visualstudio.com/items?itemName=trond-snekvik.simple-rst>`__
  (``trond-snekvik.simple-rst``) — provides the RST grammar used to highlight the
  embedded docstrings (see below).
* `C/C++ <https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools>`__
  (``ms-vscode.cpptools``) — for the handful of C sources.
* `Galacticus <https://marketplace.visualstudio.com/items?itemName=GalacticusOrg.galacticus-code>`__
  (``GalacticusOrg.galacticus-code``) — Galacticus-specific support, described in
  :ref:`manual-sec-editorExtension`.

.. _manual-sec-editorLanguageServer:

Fortran language server
-----------------------

The Modern Fortran extension uses the `fortls <https://fortls.fortran-lang.org/>`__
language server to provide go-to-definition, hover, signature help, and
completion. Install it (together with the ``findent`` indenter) via the
``editor`` optional-dependencies group of the Python package:

.. code-block:: bash

   pip install -e .[editor]

The repository-root ``.fortls`` file points the server at the ``source/``
directory and excludes the build directories. Note that Galacticus generates a
large amount of Fortran at build time (the ``!![ ... !!]`` directives are
expanded by the build scripts, and additional include files are produced under
``work/build``); symbols that only exist in generated code will therefore not be
resolved by the language server, which indexes the un-preprocessed source only.

.. _manual-sec-editorIndentation:

Indentation and formatting
--------------------------

Indentation while typing is handled by the editor's built-in auto-indent. The
Galacticus source is indented with entity bodies (``module`` / ``program`` /
``subroutine`` / ``function``) at two spaces and block constructs (``do`` /
``if`` / ``select`` and similar) at three.

.. warning::

   **Do not run a whole-file Fortran formatter on Galacticus source.** Formatters
   such as ``findent`` and ``fprettify`` do not understand the embedded, non-Fortran
   directive blocks (:ref:`manual-sec-editorEmbedded`): they re-indent the raw
   ``<name>`` / ``<method>`` lines inside ``!![ ... !!]`` blocks and relocate the
   ``!$`` OpenMP sentinels, corrupting the file. For this reason whole-file
   formatting is disabled in the committed ``.vscode/settings.json``.

If you want to re-indent a hand-selected region of *pure* Fortran (containing no
embedded directive or docstring blocks), you can run ``findent`` on that
selection only. The flags that reproduce the Galacticus style are:

.. code-block:: bash

   findent -i3 -r2 -m2 --openmp=0

where ``-i3`` sets the indent for block constructs, ``-r2`` the indent for
procedure bodies, ``-m2`` the indent for module bodies, and ``--openmp=0``
leaves the ``!$`` sentinels untouched.

.. _manual-sec-editorEmbedded:

Embedded languages in Fortran source
------------------------------------

Galacticus ``.F90`` files embed two non-Fortran languages inside ``!!`` comment
markers, both expanded by the build/documentation tooling:

* ``!![ ... !!]`` — **XML** directive blocks (for example ``<functionClass>``,
  class implementation registrations, ``<inputParameter>``, and many others).
* ``!!{RST ... !!}`` — **reStructuredText** docstrings, including inline
  ``:math:`` (LaTeX), ``:term:``, and ``:ref:`` roles. These docstrings are
  extracted at documentation-build time by ``scripts/doc/extractDocsRST.py`` to
  produce the physics pages of this manual.

The :ref:`Galacticus extension <manual-sec-editorExtension>` highlights both of
these blocks with the appropriate language grammar.

.. _manual-sec-editorExtension:

The Galacticus extension
------------------------

The ``GalacticusOrg.galacticus-code`` extension
(`source <https://github.com/galacticusorg/galacticus-code>`__) adds
Galacticus-specific support:

* **Embedded syntax highlighting** of the XML and RST blocks described in
  :ref:`manual-sec-editorEmbedded`.
* **Open Documentation for functionClass** — a command (also on the ``.F90``
  editor right-click menu, and in the Command Palette as
  ``Galacticus: Open Documentation for functionClass``, or the keybinding
  ``Ctrl+K Ctrl+G`` / ``Cmd+K Cmd+G`` while a ``.F90`` file is focused) that opens
  the online documentation for the class defined in the current file. It
  recognises both a ``functionClass`` base class and a concrete implementation,
  and deep-links to the family page at the corresponding ``physics-<name>``
  anchor. The base URL is configurable via the ``galacticus.docsBaseUrl`` setting.
