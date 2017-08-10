/****************************************************************************
* ImatiSTL                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2013: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is dual-licensed as follows:                                 *
*                                                                           *
* (1) You may use ImatiSTL as free software; you can redistribute it and/or *
* modify it under the terms of the GNU General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or      *
* (at your option) any later version.                                       *
* In this case the program is distributed in the hope that it will be       *
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
* (2) You may use ImatiSTL as part of a commercial software. In this case a *
* proper agreement must be reached with the Authors and with IMATI-GE/CNR   *
* based on a proper licensing contract.                                     *
*                                                                           *
****************************************************************************/

#ifndef _TRIMESH_H
#define _TRIMESH_H

#include "tmesh.h"

namespace IMATI_STL
{

//! TriMesh

//! Extends Basic_TMesh with additional methods.
//! This class does not add fields, thus it is safe to type-cast a Basic_TMesh to a TriMesh.

class TriMesh : public Basic_TMesh
{
	public:

		/////////////////////////////////////////////////////////////////////////////
		//
		// Constructors (inherited from Basic_TMesh)
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Empty triangulation. Should be used only prior to a call to load().
		TriMesh() : Basic_TMesh() {}

		//! Pre-defined triangulation. Currently, only "triangle" and "tetrahedron"
		//! are recognized.
		TriMesh(const char *c) : Basic_TMesh(c) {}

		//! Clones an existing Trianglation.
		TriMesh(const Basic_TMesh *b, const bool clone_info = false) : Basic_TMesh(b, clone_info) {}

		//! Clones an existing connected component.
		TriMesh(const Triangle *t, const bool keep_ref = false) : Basic_TMesh(t, keep_ref) {}

		//! Returns true only if object is a basic TriMesh. All the reimplementations must return false.
		TMESH_VIRTUAL bool isBaseType() const { return true; }



		/////////////////////////////////////////////////////////////////////////////
		//
		// Global Operations
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Cuts longest edges until all of them are shorter than 'epsilon'. O(NlogN).
		//! If maxnumv != 0, the refinement stops as soon as the mesh has max_numv vertices.
		int epsilonSample(coord squared_epsilon, bool only_on_selection =false, int max_numv =0);

		//! Uniformly re-samples the surface using 'numver' vertices and 'ns' iterations. 
		int uniformRemesh(int ns, int numver);

		/////////////////////////////////////////////////////////////////////////////
		//
		// Advanced geometrical fixes (Implemented in "MESH_STRUCTURE/
		// cutIntersections.cpp" and "MESH_STRUCTURE/computeOuterHull.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! This function resolves all the intersections using a mixed arithmetic
		//! kernel. Resulting nonmanifold edges are stored as multiple copies of boundary edges.
		//! Returns true on success - theoretically it should never fail, but in a very few
		//! cases it does (probably due to some unconsidered configs that I cannot figure out !).
		//! If 'timeout_secs' is nonzero, the function stops with failure after the specified time.
		bool cutIntersections(bool print_debug_info =false, int timeout_secs =0);

		//! This function resolves all the intersections and transforms the mesh
		//! into its outer hull.
		//! The outer hull is the set of closed surfaces that bound solid regions.
		//! Possible remaining open surfaces are returned as a new mesh.
		//! If 'timeout_secs' is nonzero, the function stops with NULL after the specified time.
		TriMesh *computeOuterHull(bool print_debug_info = false, int timeout_secs = 0);

		//! Displace the surface outward by an offset distance.
		//! Result may have degenerate and self-intersecting triangles if the input has either
		//! flat or concave dihedral angles.
		void toRawOffset(const coord& offset_distance);

		//! Turn the surface into a thin solid by creating an inverted copy.
		//! Possible initial boundaries are turned into triangles strips.
		TriMesh *toThinShell(coord& thickness);


		//! Remeshes the entire mesh based on a uniform sampling.
		//! New vertices are at the intersections of the original mesh with
		//! a voxel grid whose size (i.e. num. voxels per side) is specified
		//! by 'gridsize'. The mesh topology is completely rebuilt.
		//! The resulting mesh may have nearly-degenerate elements and other
		//! sorts of defects that one may need to remove using 'meshclean()'.
		//! If 'simplify_result' is true, coplanar vertices are merged.
		TMESH_VIRTUAL void latticeRemesh(UINT16 gridsize = 128, bool simplify_result = false);

	};

} //namespace T_MESH

#endif //_TRIMESH_H

