# Contributing to Galacticus

Thank you for considering contributing to Galacticus! We welcome contributions from developers, scientists, and researchers who are passionate about galaxy formation physics and computational methods. Whether you're fixing bugs, adding features, improving documentation, or using AI tools to accelerate your development, your work helps advance scientific discovery.

## Types of Contributions

We welcome many types of contributions:

- **Bug reports** - Help us identify and fix issues
- **Features & enhancements** - New functionality and improvements to existing code
- **Documentation** - Updates to comments, docstrings, wiki pages, and guides
- **Performance improvements** - Optimizations and efficiency gains
- **AI-assisted contributions** - Code developed with tools like GitHub Copilot, Claude, ChatGPT, etc.

*All contributions, regardless of source, must be thoroughly verified by a human developer before submission.*

## Getting Started

### Prerequisites

Before developing, ensure you have:

- A modern Fortran compiler (`gfortran` ≥ 16; earlier versions will not compile Galacticus)
- `make`
- HDF5, FFTW3, and GSL libraries
- Python 3 (≥ 3.9)

**Tip:** Use [GitHub Codespaces](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=204087682) for a pre-configured environment with all dependencies installed.

For detailed setup instructions, see the [README](README.md).

### Building & Testing

1. Clone the repository:
   ```bash
   git clone https://github.com/galacticusorg/galacticus.git
   cd galacticus
   ```

2. Install Python dependencies (puts the modules under `python/` on the import path and pulls in numpy / scipy / h5py / lxml / matplotlib / …):
   ```bash
   pip install -e '.[test]'
   ```
   The `[test]` extra also installs `pytest`. Drop it if you don't need to run the unit tests.

3. Build Galacticus:
   ```bash
   make -jN Galacticus.exe
   ```
   Replace `N` with the number of CPU cores to use.

4. Run the quick test:
   ```bash
   ./Galacticus.exe parameters/quickTest.xml
   ```

5. (Optional) Run the Python unit tests:
   ```bash
   export GALACTICUS_EXEC_PATH=`pwd`
   python -m pytest -v
   ```
   These cover the build-system modules under `python/` and the documentation tooling under `scripts/doc/`. They do **not** run the model regression tests.

6. (Optional) Run the model regression suite — the tests CI relies on to check the compiled model:
   ```bash
   python3 testSuite/test-all.py
   ```
   This runs every `testSuite/test-*.py` script (logs are written to `testSuite/outputs/`); you can also run one directly, e.g. `python3 testSuite/test-Python-interface.py`. Most of these require a built `Galacticus.exe` and the run-time datasets (`GALACTICUS_DATA_PATH`).

See the [Testing and Continuous Integration](https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/continuous-integration.html) guide for full details on the test suites and the CI pipeline.

### Installing the git hooks

Galacticus uses git hooks (kept in the separate [`gitHooks`](https://github.com/galacticusorg/gitHooks) repository) to run pre-commit checks and to enforce the commit-message format. They are recommended for all contributors and are **required** if you commit directly to the repository. The simplest way to install them is to clone the hooks repository and point your Galacticus clone at it:

```bash
git clone https://github.com/galacticusorg/gitHooks.git
git config core.hooksPath /path/to/gitHooks
```

The hooks provide:

- **`commit-msg`** — enforces [Conventional Commits](https://www.conventionalcommits.org) format (see [Submitting a Pull Request](#submitting-a-pull-request) below for the allowed types).
- **`pre-commit`** — runs static checks on staged files (Fortran static analysis; `.bib`, XML, YAML, and Python validation; docstring and spell checks; leftover-debug detection) and, if [`claude`](https://docs.claude.com/en/docs/claude-code/overview) is installed, an automated review.
- **`prepare-commit-msg`** — pre-fills a suggested Conventional-Commits message (when `claude` is installed).
- **`pre-push`** — asks you to confirm before pushing directly to `master`.

See the [`gitHooks` README](https://github.com/galacticusorg/gitHooks) for more detail.

For more help, see the [README troubleshooting section](README.md#troubleshooting).

## Development Workflow

### Creating a Branch

For new features or bug fixes, create a branch:

```bash
git checkout -b feature/your-feature-name
```

For simple changes (e.g. fixing a typo) you can commit directly to `master`; for anything larger, use a branch and open a pull request.

### Making Changes

When making changes:

1. **Follow code conventions** - See the [Coding](https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/coding.html) documentation for detailed style guidelines, naming conventions, and component patterns
2. **Test locally** - Build and test your changes before submitting
3. **Add yourself as a contributor** (see [Contributor Attribution](#contributor-attribution) below)
4. **Document your changes** - Update comments and documentation as needed

### Testing Before Submission

Before submitting a pull request, ensure:

- [ ] The code builds successfully: `make -jN Galacticus.exe`
- [ ] The quick test passes: `./Galacticus.exe parameters/quickTest.xml`
- [ ] Any new code you added is tested (or explain why not in the PR)
- [ ] Your changes don't break existing tests or functionality
- [ ] Documentation and comments are clear and updated

## Submitting a Pull Request

When your changes are ready:

1. **Push your branch** to your fork:
   ```bash
   git push origin your-branch-name
   ```

2. **Create a pull request** - Use the [PR template](.github/pull_request_template.md), which asks for:
   - A summary of your changes
   - Type of change (bug fix, feature, documentation, etc.)
   - How to test your changes
   - Confirmation that you've tested locally

3. **Follow Conventional Commits** - The `commit-msg` hook enforces the [Conventional Commits](https://www.conventionalcommits.org) format `type(scope): summary`. Use clear, descriptive messages:
   ```
   fix: resolve memory leak in component initialization
   feat: add new parameter for galactic winds model
   docs: update building instructions for gfortran 16
   ```
   The allowed types are: `fix`, `feat`, `build`, `chore`, `ci`, `docs`, `style`, `test`, `refactor`, `perf`, `revert`, and `clean`. For a breaking change, append `!` after the type/scope (e.g. `feat!: …`) and include a `BREAKING CHANGE:` footer.

4. **Address CI/CD checks** - Automated checks run on your PR (see [Testing and Continuous Integration](https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/continuous-integration.html) for what runs and when). If any checks fail, review the error messages and update your code accordingly.

5. **Respond to code review** - A maintainer will review your PR and may request changes. This is normal and helps maintain code quality!

## Contributor Attribution

We use inline markers to track who contributed to each part of the code. This is automatically extracted to generate contributor lists.

### Fortran Contributions

Add a comment at the top of the file or function:

```fortran
module myModule
  !+ This module was created by Jane Developer
  !+ Contributions to this file made by: Jane Developer, Bob Smith
  
  implicit none
  
contains
  
  subroutine mySubroutine()
    !+ Implementation by Jane Developer with assistance from Claude
    ! ... rest of code
  end subroutine mySubroutine
  
end module myModule
```

The [extractContributors.py](scripts/doc/extractContributors.py) script automatically generates contributor lists from these markers.

## AI-Assisted Contributions

We welcome contributions developed with AI tools like **GitHub Copilot, Claude, ChatGPT**, and similar services. The key requirement is **human verification**.

### What This Means

AI is an excellent development tool that can help you:

- Generate boilerplate code faster
- Suggest implementations for routine tasks
- Refactor existing code
- Write documentation and comments

However, **you are responsible** for:

1. **Reviewing the code** - Ensure it's correct, efficient, and follows project conventions
2. **Testing thoroughly** - AI-generated code needs the same level of testing as human-written code
3. **Understanding the logic** - You should be able to explain every part of your code
4. **Catching issues** - Look for edge cases, performance problems, and style inconsistencies that AI may have missed

### Typical Workflow

1. **Use the AI tool** - Prompt Claude, Copilot, etc. to help with your implementation
2. **Review what it generated** - Read through all the code carefully
3. **Test it locally** - Build and run your changes, including edge cases
4. **Fix any issues** - Update the code based on your testing and review
5. **Add attribution comments** - Document that AI tools assisted you (optional but appreciated)
6. **Submit your PR** - You're now the verified human author of this contribution

### Attribution

If you used AI tools, add a comment noting the assistance:

```fortran
subroutine doSomething()
  !+ This subroutine was generated with Claude and thoroughly tested by Jane Developer
  ! ... implementation
end subroutine doSomething
```

This helps future maintainers understand the code's origin and appreciate your use of modern development tools.

### Why This Matters

Human verification ensures:

- **Correctness** - The code actually does what it's supposed to do
- **Quality** - It meets project standards and best practices
- **Maintainability** - Future developers can understand and modify it
- **Responsibility** - You stand behind the code you submit

This is **not a distrust of AI**—it's professional responsibility. Just as you'd review code from any source before merging it, AI-generated code needs human eyes. In fact, using AI tools responsibly can help you write better code faster.

## Getting Help

- **Questions?** Ask in the [discussion forum](https://github.com/galacticusorg/galacticus/discussions)
- **Found a bug?** Open an [issue](https://github.com/galacticusorg/galacticus/issues) with details about your system and the error
- **Need more details?** See the [developer guide](https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/index.html) on ReadTheDocs
- **Development docs?** Check the [Development](https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/development.html) documentation for build system details and the [Coding](https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/coding.html) documentation for code conventions

## License

By contributing to Galacticus, you agree that your contributions will be licensed under the same [GPL-3.0 license](COPYING) as the project. This ensures that the code remains open and available to the scientific community.

---

**Thank you for contributing to Galacticus!** We appreciate your time and effort in helping advance this project.
