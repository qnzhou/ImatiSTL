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

#include "imatistl.h"

namespace IMATI_STL
{

	void TriMesh::toRawOffset(const coord& offset)
	{
		Node *n;
		Vertex *v;
		FOREACHVERTEX(v, n) v->info = new Point(v->getNormal());
		FOREACHVERTEX(v, n) {
			v->setValue((*v) + ((*((Point *)v->info))*offset));
			delete ((Point *)v->info);
			v->info = NULL;
		}
	}

	TriMesh *TriMesh::toThinShell(coord& thickness)
	{
		TriMesh *reverse = (TriMesh *)newObject(this);
		Node *n = E.head(), *m = reverse->E.head();
		Edge *ep, *er;
		while (m != NULL)
		{
			ep = (Edge *)n->data;
			er = (Edge *)m->data;
			ep->info = er; er->info = ep;
			n = n->next();
			m = m->next();
		}

		reverse->flipNormals();
		moveMeshElements(reverse);
		deselectTriangles();
		toRawOffset(thickness);
		Triangle *t;
		FOREACHEDGE(ep, n) if (ep->isOnBoundary())
		{
			bridgeBoundaries(ep, (Edge *)ep->info);
			t = (Triangle *)(T.head()->data); MARK_VISIT(t);
			t = (Triangle *)(T.head()->next()->data); MARK_VISIT(t);
		}

		Vertex *v;
		Point *p, nor;
		FOREACHVERTEX(v, n) v->info = new Point();
		FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
		{
			nor = t->getNormal();
			p = (Point *)t->v1()->info; (*p) += nor;
			p = (Point *)t->v2()->info; (*p) += nor;
			p = (Point *)t->v3()->info; (*p) += nor;
			UNMARK_VISIT(t);
		}

		FOREACHVERTEX(v, n)
		{
			p = (Point *)v->info;
			if (!p->isNull())
			{
				p->normalize();
				v->setValue((*v) + (*p)*thickness);
			}
			delete (p);
			v->info = NULL;
		}

		return this;
	}


} //namespace IMATI_STL
