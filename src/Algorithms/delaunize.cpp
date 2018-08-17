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

#include "tin.h"
#include "edgeHeap.h"

using namespace T_MESH;

//coord squaredHeightCostFunction(Edge *e)
//{
//	if (e->isOnBoundary()) return DBL_MAX; // Must be not on boundary
//	Vertex *ov1 = e->t1->oppositeVertex(e);
//	Vertex *ov2 = e->t2->oppositeVertex(e);
//	if (e->v1->exactOrientation(e->v2, ov1, ov2) != 0) return DBL_MAX; // Must be flat
//
//	if (!(!e->v1->exactSameSideOnPlane(e->v2, ov1, ov2) && ov1->exactMisalignment(e->v1, ov2) && ov1->exactMisalignment(e->v2, ov2)))
//		return DBL_MAX; // Must be properly convex
//
//	coord e2 = e->squaredLength();
//	coord h1 = Point::squaredTriangleArea3D(e->v1, e->v2, ov1) / e2;
//	coord h2 = Point::squaredTriangleArea3D(e->v1, e->v2, ov2) / e2;
//
//	return MIN(h1, h2);
//}

coord delaunayMinAngleSqSin(Edge *e) { return e->delaunayMinAngleSquaredSin(); }

bool Basic_TMesh::delaunizeFlatAreas()
{
	if (E.numels()<1) return false;		// Check for bad input

	Edge *e;
//	edgeHeap *eh = new edgeHeap(E, &squaredHeightCostFunction);		// Allocate the priority queue
	edgeHeap *eh = new edgeHeap(E, &delaunayMinAngleSqSin);		// Allocate the priority queue

	TMesh::begin_progress();

	while (!eh->isEmpty())	// Enter only if there are edges to be processed
	{
		e = eh->popHead();		// Pick the first one (which is the one with minimum cost)
		coord prev_cost = edgeHeap::getEdgeCost(e);
		if (prev_cost == DBL_MAX) break;

		if (e->swap())
		{
//			if (squaredHeightCostFunction(e) <= prev_cost) e->swap(1);
			if (delaunayMinAngleSqSin(e) <= prev_cost) e->swap(1);
			else
			{
				eh->update(e);
				eh->update(e->t1->nextEdge(e));
				eh->update(e->t1->prevEdge(e));
				eh->update(e->t2->nextEdge(e));
				eh->update(e->t2->prevEdge(e));
			}
		}
	}

	TMesh::end_progress();
	delete eh;

	return true;
}


