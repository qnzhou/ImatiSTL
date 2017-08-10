* Command line usage

- "make ImatiSTL" compiles the library only.

- "make test" compiles a test binary.

- "make" compiles both.

- "make clean" to clean objects and targets.

- "make debug" for debugging symbols.

- "make release" for optimization.



* Makefile mantainance and usage

- For the library, add new source directories and new filenames of source code modules to Makefile.ImatiSTL

- You can use Makefile.test as a reference for compiling projects that include ImatiSTL on UNIX.

- See https://github.com/ggrocca/eggmake for further information on how to use eggmake.



* Compile options

By default only static libraries are compiled. To compile dynamic
libraries, uncomment dynamic library support in Makefile.ImatiSTL.

As it is, the Makefile provided produces a library compiled for 64bit architectures
in "Fast" mode (see Readme). To change these default settings just edit the Makefile.ImatiSTL.


* Code practices for UNIX support

- Linux is case-sensitive, always include files using the correct case.

- Remember to use proper const modifiers for reference paramaters in functions.
