Compile-Time Error Messages
===========================

If you are developing Galacticus (adding new features, or modifying existing ones) it's very likely that you'll eventually run into an error while trying to compile Galacticus. This page lists common types of such compile-time errors, along with guidance on what they mean and how to fix them.

For *run-time* errors (i.e. things that go wrong when you are actually running the compiled executable, ``Galacticus.exe``) you should refer to :doc:`run-time-errors`. For other problems (Galacticus is not performing as you expect, or doing what you want) see the general :doc:`index`.

To use the information on this page, search for an error message that most closely resembles that which you receive when you try to compile your code. Note that, throughout this page, we make use of the commonly used `placeholder names <https://en.wikipedia.org/wiki/Foobar>`_ ``foo``, ``bar``, and ``baz`` to represent the names of variables and functions. Of course, your error message will contain the names of the *actual* variables and functions in Galacticus.

Understanding the general content of error messages
---------------------------------------------------

Errors from the ``gfortran`` compiler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An error message originating from the ``gfortran`` compiler will typically have a form like this:

.. code-block:: text

   source/foo.F90; line 76 [preprocessed line 86]; code 34

      86 |   call bar(baz)
         |                                  1
   Error: Symbol 'baz' at (1) has no IMPLICIT type
   make: *** [work/build/foo.o] Error 1

The first line:

.. code-block:: text

   source/foo.F90; line 76 [preprocessed line 86]; code 34

tells you where the error occurred. In this case, the error occurs in file ``source/foo.F90``, and at (or close to) line ``76`` of that file. Also reported is the error code (``34`` here), which is typically not too useful. You'll also see that this line reports ``preprocessed line 86``. Galacticus first pre-processes any Fortran file from the ``source/`` directory (manipulating the code in a variety of ways), and emits it to a similarly named file in ``work/build/`` with an extension ``.p.F90`` (so, ``source/foo.F90`` is pre-processed to become ``work/build/foo.p.F90``). This part of the error message tells you the corresponding line causing the error in the pre-processed file - which is sometimes useful in tracking down errors that arise as a result of that pre-processing.

The middle part of the error message:

.. code-block:: text

      86 |   call bar(baz)
         |                                  1
   Error: Symbol 'baz' at (1) has no IMPLICIT type

shows you the line of code that caused the error, along with the error message itself. Here the compiler is telling us that we're trying to use a symbol (a "symbol" refers to either a variable or a function name) ``baz`` of which it is unaware.

The last line:

.. code-block:: text

   make: *** [work/build/foo.o] Error 1

is simply the ``make`` command telling you that, because an error occurred (in trying to make ``work/build/foo.o`` in this case) it is giving up.

Errors from the pre-processor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may also see errors that originate from Galacticus' pre-processor itself.  These appear as Python tracebacks ending with a raised exception; an illustrative example is:

.. code-block:: text

   Traceback (most recent call last):
     ...
     File "./python/Galacticus/Build/Baz.py", line 354, in foo
       raise ValueError("no matching unit closer was found for opener 'function bar(a,b,c) result(self)'")
   ValueError: no matching unit closer was found for opener 'function bar(a,b,c) result(self)'

These can be identified by the ``python/Galacticus/Build`` path in the traceback, which corresponds to the pre-processor.  In the above, ``foo`` indicates the name of the preprocessor function which identified the error, and ``./python/Galacticus/Build/Baz.py`` indicates the file in which that function can be found.  The error message itself is, in this case, ``no matching unit closer was found for opener 'function bar(a,b,c) result(self)'`` (which means that a ``function`` statement in a source file had no matching ``end function`` line).

Common error messages
----------------------

Missing (or misnamed) function closer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   ValueError: Build_Children: no matching unit closer was found for opener 'function foo(a,b,c) result(self)'

(raised from ``python/Galacticus/Build/SourceTree.py``)

This error means that the preprocessor could not find a matching closer line in a source file. An *opener* is a line which begins a function (or subroutine, or module, or other code unit), such as ``function foo`` in the above. Every such opener must have a matching *closer*, in this case ``end function foo``.

If this error occurs then either:

#. The closer line is entirely missing;
#. The name of the function is incorrect, e.g. ``end function foob`` would not match ``function foo``; or
#. The name of of the function is missing from the closer, e.g. ``end function`` does not match ``function foo``.

To fix this error, check that you have a matching ``end`` line in the source file, and that the name of the function/subroutine/module in that ``end`` line matches that in the opener.

Derived-type must be abstract
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   Error: Derived-type ‘foo’ declared at (1) must be ABSTRACT because ‘bar’ is DEFERRED and not overridden

Generally means that the named method, ``bar``, in the derived-type (i.e. ``class``) named ``foo`` has not been given a concrete implementation in the class - so the class can not be a concrete class. The parent class has marked this method as ``deferred`` which means that it provides no implementation of the method, and so can not be a concrete class, and must instead be ``abstract`` (a class which is just a template for other classes, and can not actually be used itself). Any child class of this parent must therefore either provide an implementation of this method, or must itself be marked a being ``abstract``.

Commonly occurs if you forget to implement a function for the method, implemented the function (``baz``) but forget to add a suitable

.. code-block:: fortran

   procedure :: bar => baz

line in the ``class`` definition, or misspelled ``bar`` in this line.

To fix this problem check that you *have* a function in your class implementing this method, that you have the appropriate ``procedure`` line in the ``class`` definition, and that you have spelled the names of the method and function correctly.

Result mismatches
~~~~~~~~~~~~~~~~~~

.. code-block:: text

   procedure :: bar => fooBar; Error: Result mismatch for the overriding procedure ‘bar’ at (1): Rank mismatch in function result (2/1)

This error is telling you that the result returned by the function ``fooBar`` has a different rank than the function associated with method ``bar`` in the parent class. In this specific error ``Rank mismatch in function result (2/1)`` means that our function ``fooBar`` is returning a 2D array, while the parent classes function returned a 1D array. Ranks of function results of method must always match between parent and child classes.

To fix this, change the rank of the function result in the function ``fooBar`` to match that in the parent. For example, in the above, the function result may have been defined using something like:

.. code-block:: fortran

   double precision, allocatable, dimension(:,:) :: fooBar

but should be:

.. code-block:: fortran

   double precision, allocatable, dimension(:) :: fooBar

to match the rank from the parent class.

Private components
~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   Error: Component ‘foo’ at (1) is a PRIVATE component of ‘bar’

This is often a confusing error message, which can be misleading. Usually, the cause is that you are calling the constructor for some class, but the arguments you are passing to it do not match those of any of the available constructors for that class. (The compiler then assumes that you're trying to default constructor, which simply sets all member variables of the object directly - but you are not! Since the content of the class is marked ``private`` these member variables can not be set, resulting in an error.)

Usually the cause of this error is a mismatch in the arguments you are providing to the constructor call, and those of the actual constructor that you intended to call. The mismatch could be in the number of arguments (i.e. missing one or more arguments, or having an extra argument), an argument with the wrong type, wrong dimension, of other attribute mismatch.

To fix this error, find the constructor function that you intended to call and check carefully that the arguments you are providing in the call to the constructor match those that the constructor function is expecting.

``Parse_Directives`` failed to parse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   Traceback (most recent call last):
     ...
     File "python/Galacticus/Build/Directives.py", line 134, in _parse_xml_block
       raise RuntimeError(
   RuntimeError: extract_directives: failed parsing XML while extracting directive 'functionClass' from 'source/star_formation/rates/disks/_class.F90': mismatched tag: line 12, column 2
   XML content was:
   <functionClass>
   <name>starFormationRateDisks</name>
   <descriptiveName>Rates for star formation in disks.</descriptiveName>
   <description>Class providing models of rates of star formation in disks.</description>
   <default>intgrtdSurfaceDensity</default>
   <method name="rate" >
   <description>Returns the rate (in units of $\mathrm{M}_\odot$ Gyr$^{-1}$) for star formation in the disk component of {\normalfont \ttfamily node}.</description>
   <type>double precision</type>
   <pass>yes</pass>
   <selfTarget>yes</selfTarget>
   <argument>type(treeNode), intent(inout), target :: node</argument>
   </functionClass>

This is an error message from the pre-processor. The pre-processor was attempting to parse a block of XML embedded in a Galacticus source file, but parsing that XML failed. The error message will include the XML that was being parsed (in the above, the section starting ``<functionClass>``), and, below this, the error message reported by the XML parser. In this particular case, the XML parser reports ``mismatched tag at line 12, column 2``. Looking at this line in the block of XML above we see that we have a closing ``</functionClass>`` tag, but this does not match the currently open tag, ``<method>`` - that means we are missing a closing ``</method>`` tag here.

To fix this type of problem, use the error message from the XML parser to help you understand why the XML is syntactically incorrect, and then correct it.

Missing keyword name
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   foo=bar(compute=.false.,fast=.false.,baz)
      |                                 1
   Error: Missing keyword name in actual argument list at (1)

Here we make a function call and make use of named arguments (e.g. ``compute=.false.`` here - we are naming the argument ``compute`` of the function and assigning it a value, ``=.false.``, which allows us to put the arguments in any order, not just the order that they appear in the function itself). In this case, after the first named argument, *all* arguments must be named. So, the above should be changed to:

.. code-block:: fortran

   foo=bar(compute=.false.,fast=.false.,darkMatterHaloScale_=baz)

for example (assuming that the argument you want to pass ``baz`` to was named ``darkMatterHaloScale_`` in the function itself).

Not a member of the ``functionClass`` structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   Galacticus::Build::SourceTree::Process::Constructor::Process_Constructors():meta:26
    463 |  self%foo=foo
      |             1
   Error: ‘foo’ at (1) is not a member of the ‘functionclass’ structure

This error occurs because the object ``self`` is a member of a class that does not contain a variable called ``foo``. This type of error often occurs in constructors, particularly if you include a directive:

.. code-block:: text

   !![
   <constructorAssign variables="foo"/>
   !!]

which says to take the argument ``foo`` to this constructor and assign it to a variable named ``foo`` in the object (``self``) being constructed. We typically do this as the value of ``foo`` is something that we later want to use in one of the class' methods.

The solution to this is to add a variable named ``foo`` to the type definition for this class (with the same type and rank as the argument ``foo`` in the constructor).
