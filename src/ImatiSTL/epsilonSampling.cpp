/****************************************************************************
* ImatiSTLExt                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
*                                                                           *
* Copyright(C) 2006: IMATI-GE / CNR                                         *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include "imatistl.h"
#include "edgeHeap.h"

namespace IMATI_STL
{

coord squaredEdgeLength(Edge *e) { return e->squaredLength(); }

int TriMesh::epsilonSample(coord epsilon, bool only_on_selection, int max_numv)
{
 Node *n;
 Edge *e;
 Point p;
 int bdr, missv;
 coord msl=0.0;

 if (max_numv == 0) max_numv = INT_MAX;

 ImatiSTL::begin_progress();

 dynamicEdgeHeap *eh = new dynamicEdgeHeap(E.numels(), &squaredEdgeLength);

 if (only_on_selection)
 {
	 Triangle *t;
	 FOREACHEDGE(e, n) UNMARK_BIT(e, 6);
	 FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) { MARK_BIT(t->e1, 6); MARK_BIT(t->e2, 6); MARK_BIT(t->e3, 6); }
	 FOREACHEDGE(e, n) if (IS_BIT(e, 6))
	 {
		 UNMARK_BIT(e, 6);
		 if (e->squaredLength() > epsilon) eh->push(e);
	 }
	 missv = eh->getnum();
	 printf("subdividing %d selected edges...\n", missv);
	 while (!eh->isEmpty() && V.numels()<max_numv)
	 {
		 e = eh->popHead();
		 ImatiSTL::report_progress("%f to %f             ", sqrt(TMESH_TO_DOUBLE(epsilon)), e->length());
		 bdr = e->isOnBoundary();
		 p = e->getMidPoint();
		 splitEdge(e, &p, true);
		 if (e->squaredLength() > epsilon) eh->push(e);
		 e = ((Edge *)E.head()->data);
		 if ((((e->t1 != NULL && IS_VISITED(e->t1)) || (e->t2 != NULL && IS_VISITED(e->t2)))) && e->squaredLength() > epsilon) eh->push(e);
		 e = ((Edge *)E.head()->next()->data);
		 if ((((e->t1 != NULL && IS_VISITED(e->t1)) || (e->t2 != NULL && IS_VISITED(e->t2)))) && e->squaredLength() > epsilon) eh->push(e);
		 if (!bdr)
		 {
			 e = ((Edge *)E.head()->next()->next()->data);
			 if ((((e->t1 != NULL && IS_VISITED(e->t1)) || (e->t2 != NULL && IS_VISITED(e->t2)))) && e->squaredLength() > epsilon) eh->push(e);
		 }
	 }
 }
 else
 {
	 FOREACHEDGE(e, n) if (e->squaredLength() > epsilon) eh->push(e);

	 missv = eh->getnum();
	 while (!eh->isEmpty() && V.numels()<max_numv)
	 {
		 e = eh->popHead();
		 ImatiSTL::report_progress("%f to %f             ", sqrt(TMESH_TO_DOUBLE(epsilon)), e->length());
		 bdr = e->isOnBoundary();
		 p = e->getMidPoint();
		 splitEdge(e, &p);
		 if (e->squaredLength() > epsilon) eh->push(e);
		 e = ((Edge *)E.head()->data);
		 if (e->squaredLength() > epsilon) eh->push(e);
		 e = ((Edge *)E.head()->next()->data);
		 if (e->squaredLength() > epsilon) eh->push(e);
		 if (!bdr)
		 {
			 e = ((Edge *)E.head()->next()->next()->data);
			 if (e->squaredLength() > epsilon) eh->push(e);
		 }
	 }
 }
 
 delete(eh);
 ImatiSTL::end_progress();

 return 1;
}

} //namespace IMATI_STL
