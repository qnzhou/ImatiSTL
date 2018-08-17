/****************************************************************************
* TMesh                                                                  *
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
* (1) You may use TMesh as free software; you can redistribute it and/or *
* modify it under the terms of the GNU General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or      *
* (at your option) any later version.                                       *
* In this case the program is distributed in the hope that it will be       *
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
* (2) You may use TMesh as part of a commercial software. In this case a *
* proper agreement must be reached with the Authors and with IMATI-GE/CNR   *
* based on a proper licensing contract.                                     *
*                                                                           *
****************************************************************************/

#ifndef SIMPLIFICATION_H
#define SIMPLIFICATION_H

#include "tin.h"
#include "edgeHeap.h"

using namespace T_MESH;

///////////////////////////////////////////////////
//
// Quadric Error matrix
//
// This class stores the elements of the QEM
// associated to a vertex. The matrix is symmetric,
// so only 10 elements need to be stored explicitly.
//
// The field 'compute_optimal_point' is used to
// switch to the optimal point computation. If not
// set, the best among the mid and the two
// end-points is used as representative for the
// collapse.
//
// The field 'use_check_collapse' is used to run
// the inversion-test prior to a contraction. If
// the test fails, the contraction is assigned an
// infinite cost.
// If 'fast_check_collapse' is FALSE, then the
// algorithm checks for creation of degenerate or
// overlapping triangles too.
//
////////////////////////////////////////////////////

class Q_matrix : public SymMatrix4x4
{
 public:

 static char compute_optimal_point;	 // Read above.
 static char use_check_collapse;	 // Read above.
 static char fast_check_collapse;	 // Read above.

 Q_matrix() : SymMatrix4x4() {};	 // Empty constructor
 Q_matrix(const SymMatrix4x4& q) : SymMatrix4x4(q) {};	// Copy
 Q_matrix(Vertex *);			 // Init with incident planes

 void addPlane(Triangle *);		 // Add one plane

 Point getOptimalPoint(Edge *);		 // Return the new position
 Point bestAmongMidAndEndpoints(Edge *); // Return best of v1, v2 or mid
 int checkCollapse(Edge *, Point *);	 // Check for mesh inversion
 coord getError(Point *, Edge * =NULL); // Compute the error at a point
};

#endif // SIMPLIFICATION_H
