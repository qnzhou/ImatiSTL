----------------------------
ImatiSTL - Version 4.2                                                          
----------------------------

by Marco Attene
                                                                    
Consiglio Nazionale delle Ricerche                                        
Istituto di Matematica Applicata e Tecnologie Informatiche                
Sezione di Genova                                                         
IMATI-GE / CNR                                                            
                                                                         
ImatiSTL is a C++ library for applications that need to finely post-process raw polygon meshes represented by STL files.
It is particularly useful in 3D printing scenarios: the repairing functions provided by ImatiSTL bridge the gap between printable meshes and real-world meshes with diverse potential defects and flaws (e.g. degenerate triangles, self-intersections, surface holes, ...).
ImatiSTL provides both functions for local mesh repairing and methods for global remeshing. Local repairing is useful when the available data must remain exactly as it is in the input STL file, whereas global remeshing may fix even highly corrupted meshes at the cost of a small distortion in the resulting fixed mesh.

ImatiSTL implements an edge-based data structure with all its fundamental functionalities (i.e., file I/O, mesh construction/destruction, traversal).
It includes support for reading and writing the following file formats:
STL (http://www.sdsc.edu/tmf/Stl-specs/stl.html)
OFF (http://shape.cs.princeton.edu/benchmark/documentation/off_format.html)
PLY (http://www.cs.unc.edu/~geom/Powerplant/Ply.doc)
and partially:
IV 2.1, VRML 1.0, VRML 2.0, OBJ.

ImatiSTL is a sort of evolution of an older library called JMeshLib, but it is not backward-compatible.
In ImatiSTL V4.0, most of the old JMeshLib functionalities are part of a module called TMesh.

This package provides pre-compiled static libraries for Windows only (in lib/).
For Linux you must compile the source tree yourself.

See the comments within the source files for details or use doxygen to
produce documentation in a more readable format.

-------------------
Citation policy
--------------------
You are free to use ImatiSTL for research purposes according to the licensing terms specified at the end of this document.
If your research produces publications, please cite the following paper that describes the underlying arithmetics in ImatiSTL:

> Marco Attene. ImatiSTL - Fast and Reliable Mesh Processing with a Hybrid Kernel.
  LNCS Transactions on Computational Science, Springer (to appear).


-------------------
System Requirements
--------------------

ImatiSTL 4.2 has been tested on 32 and 64 bit PCs running:
 - Microsoft Windows OS with MSVC 12.0 (Visual C++ 2013)

Older versions were successfully compiled on Ubuntu Linux with standard development tools (gcc/g++).

ImatiSTL can exploit both MPIR and CGAL.
MPIR is used to deal with multiprecision and rational coordinates.
Precompiled .lib and .dll files are provided for Windows systems
in 'mpir32' and 'mpir64'. If you are using another system you need
to download MPIR at http://www.mpir.org/ and compile the library yourself.

Use of MPIR is required from version 2.0 to 3.1-1.
From version 3.2 use of MPIR is optional.

CGAL is used to add lazy evaluation to MPIR functionalities.
This possibility is enabled from version 4.0 on.
Lazy evaluation makes the code faster but more memory demanding.

-------------------
Building the tree
-------------------

On WINDOWS: double-click on vc12\imatistl.vcxproj and press F7.
Select the proper combination of Configuration/Platform based on your needs.
Configurations include:
'Fast' = top speed / lowest memory footprint / non-robust (self-contained / no dependancies)
'Hybrid' = slowest / intermediate memory footprint / robust (uses MPIR)
'Lazy' = intermediate speed / high memory footprint / robust (uses CGAL)

Platforms include:
Win32
x64

On Linux/Unix:
The "make" folder contains basic Makefiles to build the library in Fast mode and 64bits.
Otherwise you must create your own Makefile or edit the examples provided,
possibly inspired by the corresponding .vcxproj files.

Compilation options through preprocessor definitions:
1) IS64BITPLATFORM - set to compile a 64bit-compliant library
2) EXTENSIBLE_TMESH - set to enable polymorphism
3) USE_HYBRID_KERNEL - set to enable the hybrid kernel. If set, you need to have Mpir installed.
4) USE_CGAL_LAZYNT - set to enable CGAL's lazy evaluation. If set, you need to have CGAL installed.
5) _HAVE_PPL - set if your development environment supports Microsoft's Parallel Patterns Library

If you set USE_CGAL_LAZYNT make sure that USE_HYBRID_KERNEL is set as well.

-------------------
Using the library
-------------------

--- include path must include:
$(IMATISTL_HOME)/include/Kernel
$(IMATISTL_HOME)/include/TMesh
$(IMATISTL_HOME)/include/ImatiSTL
$(IMATISTL_HOME)/mpirXX (only if USE_HYBRID_KERNEL, XX is either 32 or 64)
CGAL\auxiliary\gmp\include (only if USE_CGAL_LAZYNT)
CGAL\build\include (only if USE_CGAL_LAZYNT)
boost_root (only if USE_CGAL_LAZYNT)
CGAL\include (only if USE_CGAL_LAZYNT)

--- preprocessor definitions
must be the same used to compile the library (see above)

--- library path must include:
$(IMATISTL_HOME)/lib (if compiled for 32bit)
$(IMATISTL_HOME)/lib64 (if compiled for 64bit)
$(IMATISTL_HOME)/mpirXX (only if USE_HYBRID_KERNEL, XX is either 32 or 64)
CGAL\build\lib\Release (only if USE_CGAL_LAZYNT)

--- static libraries to be linked
imatistl_XXX.lib (where, XXX can be either Fast, Hybrid or Lazy, possibly with 64 extension)
mpirXX.lib (only if USE_HYBRID_KERNEL, XX is either 32 or 64)
CGAL_Core-vc120-mt-YYY.lib (only if USE_CGAL_LAZYNT, where YYY is the version of CGAL you are using)
mpir.lib (only if USE_HYBRID_KERNEL alone. Do not include this when using USE_CGAL_LAZYNT)

--- other options to be set
Support for OpenMP must be active

---------
Copyright
---------

ImatiSTL is

Copyright(C) 2013: IMATI-GE / CNR                                         
                                                                          
All rights reserved.                                                      
                                                                          
This program is dual-licensed as follows:

(1) You may use ImatiSTL as free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published 
by the Free Software Foundation; either version 3 of the License, or     
(at your option) any later version.                                      
In this case the program is distributed in the hope that it will be      
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            
GNU General Public License (http://www.gnu.org/licenses/gpl.txt)         
for more details.                                                        
                                                                         
(2) You may use ImatiSTL as part of a commercial software. In this case a
proper agreement must be reached with the Authors and with IMATI-GE/CNR  
based on a proper licensing contract.