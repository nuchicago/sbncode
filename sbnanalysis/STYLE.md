Style Guide
===========
A uniform code style within the project helps improve readability, making
it easier to debug, review, and contribute to the code. The style guidelines
that follow are suggestions which will lead to uniformity across the
sbnanalysis project. In some cases, these are purely aesthetic, with no
functional impact beyond improved readability through consistency.
In other cases (with respect to globals and `using`, for example), they are
intended to make code more flexible and portable. Developers are asked to
follow these guidelines to the extent possible. Requirements that developers
*must* follow in order to contribute code are indicated explicitly.

General Guidelines
------------------
* All projects, plus all classes, functions, and variables, *must* be
  documented, with user documentation and Doxygen-style comments
* Global variables *must* be avoided
* The standard ROOT random generator (`gRandom`) *must* be used for
  random numbers
* Most code *should* fit within an 80 character wide terminal window
* Most functions *should* fit on a page
* Brand-new language features (e.g. lambda functions) *should* be avoided
  unless there is a significant benefit (these often make code harder to
  follow).
* Code blocks *should* be indented with two spaces, except the main body
  of code within the `namespace` block.

Documentation
-------------
Every analysis *must* include at minimum a detailed `README.md` file
which explains how files are organized, what the analysis is, and how to
run it. It should also identify the authors and contributors.

### Code ###
The entire `sbnanalysis` project *must* have 100% Doxygen coverage, with
no warnings reported when running `doxygen doc/Doxyfile` in the
top level directory. This means that every class, function, and member
variable requires documentation (including units where relevant!).
Exhaustive documentation is available at the Doxygen homepage, but this
section provides a few essentials.

In general, all documentation which will be built for browsing goes in the
header file, and the only documentation in the source file is in-code comments
with details one would need to understand or modify the code, but are not
necessary in order to *use* it.

### Files ###
Header files *must* include a block describing the contents (without
specifically documenting classes within), and provide a list of authors or
maintainers.

```c++
/**
 * \file ExampleFile.h
 *
 * An example file.
 *
 * This is an longer description of this file that is so long that it even
 * span multiple lines.
 *
 * Author: A. Coder <emailaddress@domain>, 2018/01/31
 */
```

### Classes ###
Classes *must* be documented with a description and documentation for all
members. Note that there is a specific comment format for Doxygen.

```c++
/**
 * \class ExampleClass
 * \brief An example class
 *
 * This class is just an example, but this would be a detailed description
 * of what it does. You can put code examples or images here, see the
 * Doxygen docs!
 */
class ExampleClass {
public:
  /** Constructor. */
  ExampleClass();

  /** Destructor. */
  virtual ~ExampleClass();

  /**
   * Divide the arguments.
   *
   * This function performs division on the arguments, using the formula:
   *
   *     OUTPUT = a / b
   *
   * \param a The numerator
   * \param b The demoninator
   * \returns The quotient
   */
  double Divide(double a, double b);

protected:
  int fNumApples;  //!< The number of apples
};
```

Headers and Forward Declarations
--------------------------------
Wherever possible, use forward declarations in headers rather than including
headers. For example, when declaring a class that references a pointer to
a class type, you can use:

```c++
class SomeType;

class MyClass {
public:
  SomeType* var;
};
```

That is, you don't need to `#include <SomeType.h>`. Forward declarations can
help clarify code and speed up compilation time. You will still need to include
the header where `SomeType` is defined in the source file where `SomeType` is
actually used.

Indeed, every file *must* include exactly what it needs, and not rely on
includes in included files.

Namespaces
----------
All code for a given analysis *must* live in a namespace matching the name
of the `ana` subdirectory. For example, `ExampleAnalysis` code is defined in
`ana::ExampleAnalysis`.

```c++
namespace ana {
  namespace AnalysisName {

// ...

  }  // namespace AnalysisName
}  // namespace ana
```

Note the convention of not indenting the main body of code, to improve
readability.

Feel free to define any namespaces you like within the analysis namespace.

It is recommended to *not* use the `using namespace` directive in source files,
as it can make code harder to understand. It should never be used in headers.

Header Guards
-------------
Header guards *must* be used in all header files, and use the following format:

```c++
#ifndef __sbnanalysis_ana_AnalysisName_ClassName__
#define __sbnanalysis_ana_AnalysisName_ClassName__
// ...
#endif  // __sbnanalysis_ana_AnalysisName_ClassName__
```

Note that these header guards should be the very first and last thing in the
header file. Source (cxx) files do not need guards.

Naming
------
* Class member variables should be prefixed with `f` and use CamelCase,
  like `fSomeVariable`. *Exception:* Cases where users will actually
  interact with the variable name, for example classes written to a ROOT file.
  Here, `words_with_underscores` is suggested.
* Function names should use CamelCase, following the convention of ROOT.
* Identifiers in `enum`s should start with `k`.
* `lowerCamelCase` is suggested for local variables

C++ Style
---------
The following style is recommended for C++ (see also code in the `core`
directory):

* Classes: See examples above.
* Functions:
```c++
void FunctionName(int argument) {
  // ...
}
```
* Long lines:
```c++
auto const& mcshower_list = \
  *(ev.getValidHandle<std::vector<sim::MCShower> >(mcshower_tag));
```
* Loops:
```c++
for (size_t i=0; i<n; i++) {
  // ...
}
```
* Conditionals:
```c++
if (2 + 2 == 5 || true) {
  // ...
}
```

Python Style
------------
Python code should follow PEP-8 standards. Python 3 is strongly recommended.

JSON
----
JSON configuration files for a processor should store all settings in an object
keyed by the processor name. The top level is reserved for framework-level
parameters, e.g. the output filename.

Example:

```
{
    "OutputFilename": "something.root",
    "MyProcessor": {
        "mysetting": 42
    }
}
```

All key names should be valid C++ identifiers. (JSON allows arbitrary strings,
but using valid identifiers is recommended here).
