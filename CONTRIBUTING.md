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

- A modern Fortran compiler (e.g., `gfortran` ≥ 11)
- `make`
- HDF5, FFTW3, and GSL libraries
- Perl and Python 3

**Tip:** Use [GitHub Codespaces](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=204087682) for a pre-configured environment with all dependencies installed.

For detailed setup instructions, see the [README](README.md).

### Building & Testing

1. Clone the repository:
   ```bash
   git clone https://github.com/galacticusorg/galacticus.git
   cd galacticus
   ```

2. Build Galacticus:
   ```bash
   make -jN Galacticus.exe
   ```
   Replace `N` with the number of CPU cores to use.

3. Run the quick test:
   ```bash
   ./Galacticus.exe parameters/quickTest.xml
   ```

For more help, see the [README troubleshooting section](README.md#troubleshooting).

## Development Workflow

### Creating a Branch

For new features or bug fixes, create a branch:

```bash
git checkout -b feature/your-feature-name
```

For simple changes, you can work on `master` and follow the [simple changes workflow](https://github.com/galacticusorg/galacticus/wiki/Contributing#making-simple-changes).

### Making Changes

When making changes:

1. **Follow code conventions** - See the [Coding](https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Development.pdf#coding) documentation for detailed style guidelines, naming conventions, and component patterns
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

3. **Follow Conventional Commits** - Use clear, descriptive commit messages:
   ```
   fix: resolve memory leak in component initialization
   feat: add new parameter for galactic winds model
   docs: update building instructions for Fortran 11
   ```
   See the [wiki](https://github.com/galacticusorg/galacticus/wiki/Contributing#commit-messages) for details.

4. **Address CI/CD checks** - Our automated tests will run on your PR. If any checks fail, review the error messages and update your code accordingly.

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

### Perl Contributions

Add comments in `.pm` files:

```perl
# Module created by Jane Developer
# Contributions: Jane Developer, Bob Smith
# This subroutine was generated with Claude and verified by Jane Developer

sub my_subroutine {
  # ... code
}
```

The [Extract_Contributors.py](scripts/doc/Extract_Contributors.py) script automatically generates contributor lists from these markers.

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
- **Need more details?** See the comprehensive [wiki](https://github.com/galacticusorg/galacticus/wiki/Contributing)
- **Development docs?** Check the [Development](https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Development.pdf#development) documentation for build system details and the [Coding](https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Development.pdf#coding) documentation for code conventions

## License

By contributing to Galacticus, you agree that your contributions will be licensed under the same [GPL-3.0 license](COPYING) as the project. This ensures that the code remains open and available to the scientific community.

---

**Thank you for contributing to Galacticus!** We appreciate your time and effort in helping advance this project.
