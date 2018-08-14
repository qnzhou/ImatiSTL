/****************************************************************************
* TMesh                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2012: IMATI-GE / CNR                                         *
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

#include "tmesh.h"

using namespace T_MESH;

//////////////////// Loop's Subdivision Scheme //////////////////

class loopSplit
{
 public:
 Edge *e;
 Point p;

 loopSplit(Edge *_e, int md)
 {
  e = _e;
  if (e->t1 == NULL || e->t2 == NULL || md) p = ((*e->v1)+(*e->v2))/2.0;
  else
  {
   Vertex *ov1 = e->t1->oppositeVertex(e);
   Vertex *ov2 = e->t2->oppositeVertex(e);
   p = ((((*e->v1)+(*e->v2))*3.0)+((*ov1)+(*ov2)))/8.0;
  }
 }
};

double subsurfbeta_loop(int k)
{
 double beta = (cos((2.0*M_PI)/((double)k))/4.0)+(3.0/8.0);
 return ((5.0/8.0)-(beta*beta))/((double)k);
}


void loopRelaxOriginal(Vertex *v)
{
 Node *n;
 Edge *e;
 List *ve;
 Point np;
 int k;
 double beta;

 if (v->isOnBoundary())
 {
  ve = v->VE();
  np = (*(((Edge *)ve->head()->data)->oppositeVertex(v)));
  np += (*(((Edge *)ve->tail()->data)->oppositeVertex(v)));
  np = (((*v)*6.0)+np)/8.0;
  delete(ve);
 }
 else
 {
  ve = v->VE();
  k = ve->numels();
  beta = subsurfbeta_loop(k);
  FOREACHVEEDGE(ve, e, n) np += (*(e->oppositeVertex(v)));
  np = ((*v)*(1.0-k*beta))+(np*beta);
  delete(ve);
 }

 v->info = new Point(&np);
}


//// Performs one subdivision step using Loop's scheme     ////
//// If 'md' is set the method does not alter the geometry ////

void Basic_TMesh::loopSubdivision(bool midpoint)
{
 List nvs;
 Edge *e, *e1, *e2;
 Node *n;
 Vertex *v;
 Triangle *t;
 loopSplit *ls;
 int k, detected_sharp=0, is_selection=0;
 if (!midpoint) FOREACHVERTEX(v, n) loopRelaxOriginal(v);

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) {is_selection=1; break;}
 if (is_selection) {FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) {MARK_BIT(t->e1, 3); MARK_BIT(t->e2, 3); MARK_BIT(t->e3, 3);}}
 else FOREACHEDGE(e, n) MARK_BIT(e, 3);

 FOREACHEDGE(e, n) if (IS_BIT(e, 3))
 {
  MARK_BIT(e->v1, 3); MARK_BIT(e->v2, 3);
  if (!midpoint && IS_SHARPEDGE(e)) detected_sharp = 1;
  nvs.appendHead(new loopSplit(e, midpoint));
 }

 if (detected_sharp)
  TMesh::warning("loopSubdivision: Crease-preservation is not supported.\n");

 FOREACHNODE(nvs, n)
 {
  ls = ((loopSplit *)n->data);
  k = ls->e->isOnBoundary();
  v = splitEdge(ls->e, &(ls->p));
  e1 = (Edge *)E.head()->data;
  e2 = (Edge *)E.head()->next()->data;
  MARK_VISIT2(v); MARK_VISIT2(e1); if (!k) MARK_VISIT2(e2);

  if (ls->e->t2)
  {
   if (IS_VISITED(ls->e->t2)) {t=(Triangle *)T.head()->data; MARK_VISIT(t);}
   if (ls->e->t1 && IS_VISITED(ls->e->t1)) {t=(Triangle *)T.head()->next()->data; MARK_VISIT(t);}
  }
  else if (IS_VISITED(ls->e->t1)) {t=(Triangle *)T.head()->data; MARK_VISIT(t);}

  if (IS_SHARPEDGE(ls->e))
  {
   if (k) TAG_SHARPEDGE(e2);
   else TAG_SHARPEDGE((Edge *)E.head()->next()->next()->data);
  }
 }

 nvs.freeNodes();

 FOREACHEDGE(e, n)
  if (IS_VISITED2(e))
  {
   UNMARK_VISIT2(e);
   if ((IS_VISITED2(e->v1) && !IS_VISITED2(e->v2)) ||
       (!IS_VISITED2(e->v1) && IS_VISITED2(e->v2)))
    if (e->swap(1) && (!IS_VISITED2(e->v1) || !IS_VISITED2(e->v2))) e->swap(1);
  }

 FOREACHVERTEX(v, n) if (!IS_VISITED2(v) && !midpoint && IS_BIT(v, 3))
 {
  v->setValue((Point *)v->info);
  delete((Point *)v->info);
  v->info = NULL;
 }
 FOREACHVERTEX(v, n) {UNMARK_VISIT2(v); UNMARK_BIT(v, 3);}

 if (detected_sharp)
  TMesh::warning("loopSubdivision: Tagged sharp edges have been smoothed.\n");

 FOREACHEDGE(e, n) UNMARK_BIT(e, 3);
}


////////////////// SQRT(3) Subdivision Scheme //////////////////////

void sqrt3RelaxOriginal(Vertex *v)
{
	Node *n;
	Vertex *w;
	Point sn;
	List *vv = v->VV();
	double val = vv->numels();

	FOREACHVVVERTEX(vv, w, n) sn += (*w);
	delete(vv);

	double an = (4.0 - 2.0*cos((2 * M_PI) / val)) / 9.0;
	sn = ((*v)*(1.0 - an)) + (sn*(an / val));
	v->info = new Point(&sn);
}

void Basic_TMesh::sqrt3Subdivision()
{
	Node *n;
	Triangle *t, *y, *b;
	Edge *e;
	Vertex *v;
	Point c;
	int detected_sharp = 0, is_selection = 0;

	FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) { is_selection = 1; break; }
	if (is_selection) { FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) { MARK_VISIT2(t->e1); MARK_VISIT2(t->e2); MARK_VISIT2(t->e3); } } else FOREACHEDGE(e, n) MARK_VISIT2(e);

	FOREACHEDGE(e, n) if (IS_VISITED2(e))
	{
		MARK_BIT(e->v1, 3); MARK_BIT(e->v2, 3);
		if (e->isOnBoundary())
		{
			TMesh::warning("sqrt3Subdivision: Meshes with boundary are not supported.\n");
			FOREACHEDGE(e, n) UNMARK_VISIT2(e);
			FOREACHVERTEX(v, n) UNMARK_BIT(v, 3);
			return;
		} else if (IS_SHARPEDGE(e)) detected_sharp = 1;
	}

	if (detected_sharp)
		TMesh::warning("sqrt3Subdivision: Feature-preservation is not supported.\n");


	FOREACHVERTEX(v, n) if (IS_BIT(v, 3)) sqrt3RelaxOriginal(v);

	if (is_selection)
	{
		FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
		{
			v = splitTriangle(t, &(c = t->getCenter()));
			b = t->oppositeEdge(v)->oppositeTriangle(t); if (!IS_VISITED(b)) UNMARK_VISIT(t);
			y = ((Triangle *)T.head()->data); b = y->oppositeEdge(v)->oppositeTriangle(y); if (IS_VISITED(b)) MARK_VISIT(y);
			y = ((Triangle *)T.head()->next()->data); b = y->oppositeEdge(v)->oppositeTriangle(y); if (IS_VISITED(b)) MARK_VISIT(y);
		}
	} else
	{
		FOREACHTRIANGLE(t, n) splitTriangle(t, &(c = t->getCenter()), 1);
	}

	FOREACHEDGE(e, n) if (IS_VISITED2(e)) { e->swap(1); UNMARK_VISIT2(e); }

	for (n = V.tail(); n != NULL; n = n->prev())
	{
		v = (Vertex *)n->data;
		if (IS_BIT(v, 3))
		{
			UNMARK_BIT(v, 3);
			v->setValue((Point *)v->info);
			delete((Point *)v->info);
			v->info = NULL;
		}
	}

	if (detected_sharp)
		TMesh::warning("sqrt3Subdivision: Tagged sharp edges have been smoothed.\n");
}









void mbs_ComputeClosedButterflyCoefs(double cf[], double *cf_center, int sz)
{
	int i;
	if (sz == 6)
	{
		cf[0] = 1.0 / 2.0;
		cf[1] = 1.0 / 8.0;
		cf[2] = -1.0 / 8.0;
		cf[3] = 0.0;
		cf[4] = -1.0 / 8.0;
		cf[5] = 1.0 / 8.0;
		*cf_center = 1.0 / 2.0;
	} else if (sz == 3)
	{
		cf[0] = 5.0 / 12.0;
		cf[1] = -1.0 / 12.0;
		cf[2] = -1.0 / 12.0;
		*cf_center = 3.0 / 4.0;
	} else if (sz == 4)
	{
		cf[0] = 3.0 / 8.0;
		cf[1] = 0.0;
		cf[2] = -1.0 / 8.0;
		cf[3] = 0.0;
		*cf_center = 3.0 / 4.0;
	} else
	{
		for (i = 0; i < sz; i++)
			cf[i] = (0.25 + cos(2 * M_PI*i / ((double)sz)) +
			0.5*cos(4 * M_PI*i / ((double)sz))) / ((double)sz);
		*cf_center = 3.0 / 4.0;
	}
}

void mbs_ComputeOpenButterflyCoefs(double cf[], double *cf_center, int i, int sz)
{
	int j;

	double thetak = M_PI / ((double)sz);
	double k = double(sz);
	double lambda1 = 0.5;
	double lambda2 = 0.25;

	cf[0] = 0.25*cos(i*thetak) -
		(1 / k)*lambda2*sin(2 * thetak)*sin(2 * thetak*i) / (cos(thetak) - cos(2 * thetak));
	cf[sz] = -cf[0];
	*cf_center =
		1.0 - (2 / k)*lambda1*(sin(thetak)*sin(i*thetak)) / (1 - cos(thetak));
	for (j = 1; j <= sz - 1; j++)
		cf[j] = (2 / k)*(lambda1*sin(i*thetak)*sin(j*thetak) +
		lambda2*sin(2 * i*thetak)*sin(2 * j*thetak));
	double sum = 0.0;
	for (j = 0; j <= sz; j++)
		sum += cf[j];
	sum += *cf_center;
}

Point mbs_applyCoefs(List *ve, double cf[], double cfc, int order)
{
	int i = 0;
	Node *n;
	Edge *e;
	Vertex *v0 = ((Edge *)ve->head()->data)->commonVertex(((Edge *)ve->tail()->data));
	Vertex *v;
	Point disp, p = (*v0)*cfc;

	FOREACHVEEDGE(ve, e, n)
	{
		v = e->oppositeVertex(v0);
		disp = ((*v)*(cf[i]));
		p = p + disp;
		i++;
	}

	return p;
}


Point mbs_sharpedge_rule(Edge *e)
{
	List *ve1 = e->v1->VE(), *ve2 = e->v2->VE();
	Edge *f, *oe1 = NULL, *oe2 = NULL;
	Vertex *ov1, *ov2;
	Point op;
	int ne1 = 0, ne2 = 0;
	Node *n;
	FOREACHVEEDGE(ve1, f, n) if (IS_SHARPEDGE(f) && f != e)
	{
		oe1 = f; ne1++;
	}
	FOREACHVEEDGE(ve2, f, n) if (IS_SHARPEDGE(f) && f != e)
	{
		oe2 = f; ne2++;
	}
	delete(ve1); delete(ve2);

	if (ne1 == 1 && ne2 == 1)
	{
		ov1 = oe1->oppositeVertex(e->v1);
		ov2 = oe2->oppositeVertex(e->v2);
		return ((((*(e->v1)) + (*(e->v2))) * 9) - ((*ov1) + (*ov2))) / 16;
	} else if (ne1 == 1)
	{
		ov1 = oe1->oppositeVertex(e->v1);
		op = (*(e->v2)) - (*(e->v1));
		op = (*ov1) + (op*((((op*((*(e->v1)) - (*ov1))) * 2) / (op*op)) + 1));
		return ((((*(e->v1)) + (*(e->v2))) * 9) - ((*ov1) + (op))) / 16;
	} else if (ne2 == 1)
	{
		ov2 = oe2->oppositeVertex(e->v2);
		op = (*(e->v1)) - (*(e->v2));
		op = (*ov2) + (op*((((op*((*(e->v2)) - (*ov2))) * 2) / (op*op)) + 1));
		return ((((*(e->v1)) + (*(e->v2))) * 9) - ((*ov2) + (op))) / 16;
	}

	return e->getMidPoint();
}


class modbutSplit
{
public:
	Edge *e;
	Point p;

	modbutSplit(Edge *_e)
	{
		e = _e;
		if (IS_SHARPEDGE(e))
		{
			p = mbs_sharpedge_rule(e);
			return;
		}
		Edge *f;
		e->v1->e0 = e->v2->e0 = e;
		List *vrTail = e->v2->VT();
		List *vrHead = e->v1->VT();
		List *veTail = e->v2->VE();
		List *veHead = e->v1->VE();
		Node *n;

		int tailSize = vrTail->numels();
		int headSize = vrHead->numels();
		bool tailClosed = (vrTail->numels() == veTail->numels());
		bool headClosed = (vrHead->numels() == veHead->numels());

		double *cfTail = new double[tailSize + 1];
		double *cfHead = new double[headSize + 1];
		double cfTail_center;
		double cfHead_center;

		int tailEdge = 0, headEdge = 0;
		FOREACHVEEDGE(veTail, f, n) if (f == e) break; else tailEdge++;
		FOREACHVEEDGE(veHead, f, n) if (f == e) break; else headEdge++;

		Vertex *ov1, *ov2, *ov3;

		if (tailClosed && headClosed && ((tailSize == 6 && headSize == 6) || (tailSize != 6 && headSize != 6)))
		{
			mbs_ComputeClosedButterflyCoefs(cfTail, &cfTail_center, tailSize);
			mbs_ComputeClosedButterflyCoefs(cfHead, &cfHead_center, headSize);
			p = (mbs_applyCoefs(veTail, cfTail, cfTail_center, 0) + mbs_applyCoefs(veHead, cfHead, cfHead_center, 1)) / 2.0;
		} else if (tailClosed && headClosed && (tailSize != 6))
		{
			mbs_ComputeClosedButterflyCoefs(cfTail, &cfTail_center, tailSize);
			p = mbs_applyCoefs(veTail, cfTail, cfTail_center, 0);
		} else if (tailClosed && headClosed && (headSize != 6))
		{
			mbs_ComputeClosedButterflyCoefs(cfHead, &cfHead_center, headSize);
			p = mbs_applyCoefs(veHead, cfHead, cfHead_center, 1);
		} else if (e->isOnBoundary())
		{
			ov1 = e->v1->prevOnBoundary();
			ov2 = e->v2->nextOnBoundary();
			if (ov1 == e->v2) ov1 = e->v1->nextOnBoundary();
			if (ov2 == e->v1) ov2 = e->v2->prevOnBoundary();
			p = ((((*(e->v1)) + (*(e->v2))) * 9) - ((*ov1) + (*ov2))) / 16;
		} else if (!tailClosed && (tailSize != 3) && headClosed && (headSize != 6))
		{
			mbs_ComputeOpenButterflyCoefs(cfTail, &cfTail_center, tailEdge, tailSize);
			mbs_ComputeClosedButterflyCoefs(cfHead, &cfHead_center, headSize);
			p = (mbs_applyCoefs(veTail, cfTail, cfTail_center, 0) + mbs_applyCoefs(veHead, cfHead, cfHead_center, 1)) / 2.0;
		} else if (!headClosed && (headSize != 3) && tailClosed && (tailSize != 6))
		{
			mbs_ComputeOpenButterflyCoefs(cfHead, &cfHead_center, headEdge, headSize);
			mbs_ComputeClosedButterflyCoefs(cfTail, &cfTail_center, tailSize);
			p = (mbs_applyCoefs(veTail, cfTail, cfTail_center, 0) + mbs_applyCoefs(veHead, cfHead, cfHead_center, 0)) / 2.0;
		} else if (!headClosed && (headSize != 3) && !tailClosed && (tailSize != 3))
		{
			mbs_ComputeOpenButterflyCoefs(cfHead, &cfHead_center, headEdge, headSize);
			mbs_ComputeOpenButterflyCoefs(cfTail, &cfTail_center, tailEdge, tailSize);
			p = (mbs_applyCoefs(veTail, cfTail, cfTail_center, 0) + mbs_applyCoefs(veHead, cfHead, cfHead_center, 0)) / 2.0;
		} else if (!tailClosed && (tailSize != 3) && ((headClosed && (headSize == 6)) || (!headClosed && (headSize == 3))))
		{
			mbs_ComputeOpenButterflyCoefs(cfTail, &cfTail_center, tailEdge, tailSize);
			p = mbs_applyCoefs(veTail, cfTail, cfTail_center, 0);
		} else if (!headClosed && (headSize != 3) && ((tailClosed && (tailSize == 6)) || (!tailClosed && (tailSize == 3))))
		{
			mbs_ComputeOpenButterflyCoefs(cfHead, &cfHead_center, headEdge, headSize);
			p = mbs_applyCoefs(veHead, cfHead, cfHead_center, 0);
		} else if (tailClosed && (tailSize != 6) && (/*(headClosed && (headSize == 6)) || */(!headClosed && (headSize == 3))))
		{
			mbs_ComputeClosedButterflyCoefs(cfTail, &cfTail_center, tailSize);
			p = mbs_applyCoefs(veTail, cfTail, cfTail_center, 0);
		} else if (headClosed && (headSize != 6) && (/*(tailClosed && (tailSize == 6)) || */(!tailClosed && (tailSize == 3))))
		{
			mbs_ComputeClosedButterflyCoefs(cfHead, &cfHead_center, headSize);
			p = mbs_applyCoefs(veHead, cfHead, cfHead_center, 1);

		} else if (!tailClosed && !headClosed)
		{
			if (tailEdge == 1 && headEdge == 2)
			{
				ov1 = ((Edge *)veHead->head()->data)->oppositeVertex(e->v1);
				ov2 = ((Edge *)veHead->head()->next()->data)->oppositeVertex(e->v1);
				ov3 = ((Edge *)veTail->tail()->data)->oppositeVertex(e->v2);
				p = (((*(e->v1)) + (*(e->v2)))*0.5) + ((*ov2)*0.25) - ((*ov1)*0.125) - ((*ov2)*0.125);
			} else if (tailEdge == 2 && headEdge == 1)
			{
				ov1 = ((Edge *)veTail->head()->data)->oppositeVertex(e->v2);
				ov2 = ((Edge *)veTail->head()->next()->data)->oppositeVertex(e->v2);
				ov3 = ((Edge *)veHead->tail()->data)->oppositeVertex(e->v1);
				p = (((*(e->v1)) + (*(e->v2)))*0.5) + ((*ov2)*0.25) - ((*ov1)*0.125) - ((*ov2)*0.125);
			} else if ((tailEdge == 1 && headEdge == 1) || (tailEdge == 2 && headEdge == 2))
			{
				p = (((*(e->v1)) + (*(e->v2)))*0.5);
			}
		} else if (!tailClosed && headClosed)
		{
			if (tailEdge == 1)
			{
				ov1 = ((Edge *)veTail->head()->data)->oppositeVertex(e->v2);
				ov2 = ((Edge *)veTail->tail()->data)->oppositeVertex(e->v2);
				f = ((Edge *)veHead->head()->next()->data);
				ov3 = f->leftTriangle(e->v1)->oppositeVertex(f);
				p = (*ov1) - ((*ov2) + (*ov3));
				f = ((Edge *)veHead->tail()->data);
				ov1 = f->rightTriangle(e->v1)->oppositeVertex(f);
				ov2 = f->oppositeVertex(e->v1);
				p = (p + ((*ov2) * 3)) - ((*ov1) * 2);
				p = p + ((*(e->v2)) * 6) + ((*(e->v1)) * 10);
				p = p / 16;
			} else
			{
				ov1 = ((Edge *)veTail->tail()->data)->oppositeVertex(e->v2);
				ov2 = ((Edge *)veTail->head()->data)->oppositeVertex(e->v2);
				f = ((Edge *)veHead->tail()->data);
				ov3 = f->rightTriangle(e->v1)->oppositeVertex(f);
				p = (*ov1) - ((*ov2) + (*ov3));
				f = ((Edge *)veHead->head()->next()->data);
				ov1 = f->leftTriangle(e->v1)->oppositeVertex(f);
				ov2 = f->oppositeVertex(e->v1);
				p = (p + ((*ov2) * 3)) - ((*ov1) * 2);
				p = p + ((*(e->v2)) * 6) + ((*(e->v1)) * 10);
				p = p / 16;
			}
		} else if (tailClosed && !headClosed)
		{
			if (headEdge == 1)
			{
				ov1 = ((Edge *)veHead->head()->data)->oppositeVertex(e->v1);
				ov2 = ((Edge *)veHead->tail()->data)->oppositeVertex(e->v1);
				f = ((Edge *)veTail->head()->next()->data);
				ov3 = f->leftTriangle(e->v2)->oppositeVertex(f);
				p = (*ov1) - ((*ov2) + (*ov3));
				f = ((Edge *)veTail->tail()->data);
				ov1 = f->rightTriangle(e->v2)->oppositeVertex(f);
				ov2 = f->oppositeVertex(e->v2);
				p = (p + ((*ov2) * 3)) - ((*ov1) * 2);
				p = p + ((*(e->v1)) * 6) + ((*(e->v2)) * 10);
				p = p / 16;
			} else
			{
				ov1 = ((Edge *)veHead->tail()->data)->oppositeVertex(e->v1);
				ov2 = ((Edge *)veHead->head()->data)->oppositeVertex(e->v1);
				f = ((Edge *)veTail->tail()->data);
				ov3 = f->rightTriangle(e->v2)->oppositeVertex(f);
				p = (*ov1) - ((*ov2) + (*ov3));
				f = ((Edge *)veTail->head()->next()->data);
				ov1 = f->leftTriangle(e->v2)->oppositeVertex(f);
				ov2 = f->oppositeVertex(e->v2);
				p = (p + ((*ov2) * 3)) - ((*ov1) * 2);
				p = p + ((*(e->v1)) * 6) + ((*(e->v2)) * 10);
				p = p / 16;
			}
		} else
		{
			TMesh::error("ModButSubdivide: Unexpected configuration !\n");
		}

		delete(vrTail);
		delete(vrHead);
		delete(veTail);
		delete(veHead);
		delete(cfTail);
		delete(cfHead);
	}

};

Basic_TMesh *mbs_cutoffSubtin(Vertex *v, Edge *e, Edge **out)
{
	Edge *f, *e1, *e2, *e3;
	Vertex *w, *v1, *v2;
	Triangle *t, *nt;
	List nts, nes, *vt;
	Basic_TMesh *tin = new Basic_TMesh;
	Node *n;

	f = e;
	do
	{
		w = f->oppositeVertex(v);
		t = f->leftTriangle(v);
		if (t == NULL) break;
		if (!IS_VISITED2(t)) { nts.appendHead(t); MARK_VISIT2(t); }
		f = t->oppositeEdge(w);
		if (f == e) break;
	} while (!IS_SHARPEDGE(f));

	f = e;
	do
	{
		w = f->oppositeVertex(v);
		t = f->rightTriangle(v);
		if (t == NULL) break;
		if (!IS_VISITED2(t)) { nts.appendHead(t); MARK_VISIT2(t); }
		f = t->oppositeEdge(w);
		if (f == e) break;
	} while (!IS_SHARPEDGE(f));

	w = e->oppositeVertex(v);
	vt = w->VT();
	FOREACHVTTRIANGLE(vt, t, n)
	if (!IS_VISITED2(t)) { nts.appendHead(t); MARK_VISIT2(t); }
	delete(vt);

	while (nts.numels())
	{
		t = (Triangle *)nts.popHead();
		UNMARK_VISIT2(t);
		if (t->e1->info == NULL) { t->e1->info = e1 = new Edge(t->e1->v1, t->e1->v2); nes.appendHead(t->e1); } else e1 = (Edge *)t->e1->info;
		if (t->e2->info == NULL) { t->e2->info = e2 = new Edge(t->e2->v1, t->e2->v2); nes.appendHead(t->e2); } else e2 = (Edge *)t->e2->info;
		if (t->e3->info == NULL) { t->e3->info = e3 = new Edge(t->e3->v1, t->e3->v2); nes.appendHead(t->e3); } else e3 = (Edge *)t->e3->info;
#ifdef USE_PER_TRIANGLE_COLORS
		nt = tin->newTriangle(e1, e2, e3, t->getColor());
#else
		nt = tin->newTriangle(e1, e2, e3);
#endif
		tin->T.appendHead(nt);
		if (t->e1->t1 == t) e1->t1 = nt; else e1->t2 = nt;
		if (t->e2->t1 == t) e2->t1 = nt; else e2->t2 = nt;
		if (t->e3->t1 == t) e3->t1 = nt; else e3->t2 = nt;
		if (IS_SHARPEDGE(t->e1)) TAG_SHARPEDGE(e1);
		if (IS_SHARPEDGE(t->e2)) TAG_SHARPEDGE(e2);
		if (IS_SHARPEDGE(t->e3)) TAG_SHARPEDGE(e3);
	}

	FOREACHVEEDGE((&(nes)), f, n)
	{
		if (f == e) *out = ((Edge *)f->info);
		if (f->v1->info == NULL)
		{
			f->v1->info = v1 = new Vertex(f->v1); tin->V.appendHead(v1); if (IS_SHARPEDGE(f->v1)) TAG_SHARPEDGE(v1);
		} else v1 = (Vertex *)f->v1->info;
		if (f->v2->info == NULL)
		{
			f->v2->info = v2 = new Vertex(f->v2); tin->V.appendHead(v2); if (IS_SHARPEDGE(f->v2)) TAG_SHARPEDGE(v2);
		} else v2 = (Vertex *)f->v2->info;
		f = (Edge *)f->info;
		tin->E.appendHead(f);
		f->v1 = v1; f->v2 = v2; v1->e0 = v2->e0 = f;
	}

	while (nes.numels()) { f = (Edge *)nes.popHead(); f->v1->info = f->v2->info = f->info = NULL; }

	return tin;
}


//// Performs one subdivision step using the modified butterfly scheme ////

void Basic_TMesh::modbutSubdivision()
{
	List nvs;
	Edge *e, *e1, *e2;
	Node *n;
	Vertex *v;
	Triangle *t;
	modbutSplit *ls;
	int k, bde, is_selection = 0;

	FOREACHEDGE(e, n) if (IS_SHARPEDGE(e))
	{
		TAG_SHARPEDGE(e->v1); TAG_SHARPEDGE(e->v2);
	}

	FOREACHVERTEX(v, n) v->info = 0;

	FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) { is_selection = 1; break; }
	if (is_selection) { FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) { MARK_BIT(t->e1, 3); MARK_BIT(t->e2, 3); MARK_BIT(t->e3, 3); } } else FOREACHEDGE(e, n) MARK_BIT(e, 3);

	FOREACHEDGE(e, n) if (IS_BIT(e, 3))
	{
		if (!IS_SHARPEDGE(e))
		{
			Basic_TMesh *st1 = NULL, *st2 = NULL;
			e2 = e;
			if (IS_SHARPEDGE(e->v1)) { st1 = mbs_cutoffSubtin(e->v1, e, &e1); e2 = e1; }
			if (IS_SHARPEDGE(e->v2))
			{
				if (st1 != NULL) st2 = mbs_cutoffSubtin(e1->v2, e1, &e2);
				else st2 = mbs_cutoffSubtin(e->v2, e, &e2);
			}
			nvs.appendHead(new modbutSplit(e2));
			((modbutSplit *)nvs.head()->data)->e = e;
			if (st1 != NULL) delete(st1);
			if (st2 != NULL) delete(st2);
		} else nvs.appendHead(new modbutSplit(e));
	}

	FOREACHNODE(nvs, n)
	{
		ls = ((modbutSplit *)n->data);
		k = (ls->e->isOnBoundary());
		bde = IS_SHARPEDGE(ls->e);
		v = splitEdge(ls->e, &(ls->p));
		e1 = (Edge *)E.head()->data;
		e2 = (Edge *)E.head()->next()->data;
		MARK_VISIT2(v); MARK_VISIT2(e1); if (!k) MARK_VISIT2(e2);

		if (ls->e->t2)
		{
			if (IS_VISITED(ls->e->t2)) { t = (Triangle *)T.head()->data; MARK_VISIT(t); }
			if (ls->e->t1 && IS_VISITED(ls->e->t1)) { t = (Triangle *)T.head()->next()->data; MARK_VISIT(t); }
		} else if (IS_VISITED(ls->e->t1)) { t = (Triangle *)T.head()->data; MARK_VISIT(t); }

		if (bde) { e = (Edge *)E.head()->next()->next()->data; TAG_SHARPEDGE(e); }
	}

	nvs.freeNodes();

	FOREACHEDGE(e, n)
	if (IS_VISITED2(e))
	{
		UNMARK_VISIT2(e);
		if ((IS_VISITED2(e->v1) && !IS_VISITED2(e->v2)) ||
			(!IS_VISITED2(e->v1) && IS_VISITED2(e->v2)))
		if (e->swap(1) && (!IS_VISITED2(e->v1) || !IS_VISITED2(e->v2))) e->swap(1);
	}

	FOREACHVERTEX(v, n) { UNMARK_VISIT2(v); v->info = NULL; }
	FOREACHEDGE(e, n) UNMARK_BIT(e, 3);
}

