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

#include "tin.h"
#include <stdlib.h>

using namespace T_MESH;

#define ES_GREY			0
#define ES_RED			8
#define ES_GREEN		16
#define ES_BLUE			32
#define ES_BROWN		64
#define ES_PAINT(e, c)	((e)->mask = (((unsigned char)(c)) | ((e)->mask & 7)))
#define ES_COLOR(e)	((e)->mask & ((unsigned char)248))

//////////////////////////////////////////////////////////////////////////
//		FILTERS							//
//////////////////////////////////////////////////////////////////////////

void filter0(List *el, double eps)
{
 Node *n;
 Edge *e;
 FOREACHVEEDGE(el, e, n) ES_PAINT(e, ((e->curvature() < eps)?(ES_BROWN):(ES_GREY)));
}

int filter1_single(Vertex *v)
{
 List *ve = v->VE();
 Node *n;
 Edge *e;
 FOREACHVEEDGE(ve, e, n) if (ES_COLOR(e) != ES_BROWN) {delete(ve); return 0;}
 delete(ve);
 return 1;
}

void filter1(List *vl)
{
 Node *n;
 Vertex *v;
 FOREACHVVVERTEX(vl, v, n) ES_PAINT(v, ((filter1_single(v))?(ES_RED):(ES_GREY)));
}

int filter2_single(Triangle *t)
{
 Vertex *v1 = t->v1(), *v2 = t->v2(), *v3 = t->v3();
 return (ES_COLOR(v1) == ES_RED || ES_COLOR(v2) == ES_RED || ES_COLOR(v3) == ES_RED);
}

void filter2(List *tl)
{
 Node *n;
 Triangle *t;
 FOREACHVTTRIANGLE(tl, t, n) ES_PAINT(t, ((filter2_single(t))?(ES_RED):(ES_GREY)));
}

void filter3_single(Triangle *t)
{
 if (ES_COLOR(t) != ES_RED) return;
 List triList;
 Triangle *s;

 triList.appendHead(t);
 while(triList.numels() > 0)
 {
  t = (Triangle *)triList.popHead();
  if ((s = t->t1()) != NULL && ES_COLOR(s) != ES_RED && ES_COLOR(t->e1) == ES_BROWN)
   {triList.appendHead(s); ES_PAINT(s, ES_RED);}
  if ((s = t->t2()) != NULL && ES_COLOR(s) != ES_RED && ES_COLOR(t->e2) == ES_BROWN)
   {triList.appendHead(s); ES_PAINT(s, ES_RED);}
  if ((s = t->t3()) != NULL && ES_COLOR(s) != ES_RED && ES_COLOR(t->e3) == ES_BROWN)
   {triList.appendHead(s); ES_PAINT(s, ES_RED);}
 }
}

void filter3(List *tl)
{
 Node *n;
 Triangle *t;
 FOREACHVTTRIANGLE(tl, t, n) filter3_single(t);
}

void filter4(List *tl)
{
 Node *n;
 Triangle *t;

 FOREACHVTTRIANGLE(tl, t, n) if (ES_COLOR(t) == ES_RED)
 {
  ES_PAINT(t->v1(), ES_RED); ES_PAINT(t->v2(), ES_RED); ES_PAINT(t->v3(), ES_RED);
  ES_PAINT(t->e1, ES_RED); ES_PAINT(t->e2, ES_RED); ES_PAINT(t->e3, ES_RED);
 } 
}

void filter5(List *el)
{
 Node *n;
 Edge *e;
 FOREACHVEEDGE(el, e, n)
  if (ES_COLOR(e) != ES_RED && !e->isOnBoundary() &&
      ES_COLOR(e->v1) == ES_RED && ES_COLOR(e->v2) == ES_RED)
       ES_PAINT(e, ES_BLUE);
}

void filter6(List *tl)
{
 Node *n;
 Triangle *t;
 FOREACHVTTRIANGLE(tl, t, n)
  if (ES_COLOR(t->e1) == ES_BLUE && ES_COLOR(t->e2) == ES_BLUE && ES_COLOR(t->e3) == ES_BLUE)
   ES_PAINT(t, ES_GREEN);
}


//////////////////////////////////////////////////////////////////////////
//		SUBDIVIDERS						//
//////////////////////////////////////////////////////////////////////////

// The pair (chamfer edge-split point)
class splitPoint
{
 public:
 Point p;	// Splitting point
 Edge *e;	// Edge to be split

 splitPoint(Point _p, Edge *_e) {p=_p; e=_e;}
};

// The pair (corner triangle-split point)
class splitCorner
{
 public:
 Point p;	// Splitting point
 Triangle *t;	// Triangle to be split

 splitCorner(Point _p, Triangle *_t) {p=_p; t=_t;}
};

Point weightedNormal(Vertex *v)
{
 Point nor;
 Node *n;
 Triangle *t;
 List *vt = v->VT();
 FOREACHVTTRIANGLE(vt, t, n) if (ES_COLOR(t) == ES_RED) nor = nor+(t->getNormal()*t->getAngle(v));
 delete(vt);
 if (nor.isNull()) return v->getNormal();
 nor.normalize();

 return nor;
}

double maxDistanceFromEdgeVertices(Point *p, Edge *e)
{
 double d1 = p->distance(e->v1);
 double d2 = p->distance(e->v1);
 return MAX(d1,d2);
}

// Computes wether the edge has to be split and, if so, where it needs to be split.
Point getPointToSplit(Edge *e, double len)
{
 if (ES_COLOR(e) != ES_BLUE) return INFINITE_POINT; // Point(COORD_MAX, COORD_MAX, COORD_MAX);

 Point pr1, pr2, cp1, cp2, sp;

 Point nor1 = weightedNormal(e->v1);
 Point nor2 = weightedNormal(e->v2);

 pr2 = nor1&nor2;
 pr1 = pr2&e->toVector();
 if ((nor1&pr1)*(pr1&nor2) <= 0) return e->getMidPoint();
 pr1 = Point(((nor1*(*e->v1))), ((nor2*(*e->v2))), 0.0).linearSystem(nor1, nor2, pr2);
 pr2 = pr1 + (pr2*len);
 if (!e->v1->closestPoints(e->v2, &pr1, &pr2, &cp1, &cp2)) return e->getMidPoint();
 sp = cp2;

 if (len != 0 && maxDistanceFromEdgeVertices(&sp, e) > len) return e->getMidPoint();

 return sp;
}


double maxDistanceFromTriangleVertices(Point *p, Triangle *t)
{
 double d1 = p->distance(t->v1());
 double d2 = p->distance(t->v2());
 double d3 = p->distance(t->v3());
 return MAX(d1,MAX(d2,d3));
}

// Computes wether the triangle has to be split and, if so, where it needs to be split.
Point getTriangleToSplit(Triangle *s, double len)
{
 if (ES_COLOR(s) != ES_GREEN) return INFINITE_POINT; // Point(COORD_MAX, COORD_MAX, COORD_MAX);

 Vertex *v1 = s->v1(), *v2 = s->v2(), *v3 = s->v3();
 Point nor1 = weightedNormal(v1);
 Point nor2 = weightedNormal(v2);
 Point nor3 = weightedNormal(v3);
 Point sp = Point((nor1*(*v1)), (nor2*(*v2)), (nor3*(*v3))).linearSystem(nor1, nor2, nor3);

 // sp not too far from the triangle
 if (len != 0 && maxDistanceFromTriangleVertices(&sp, s) > len) return s->getCenter();

 return sp;
}


//////////////////////////////////////////////////////////////////////////
//		EDGE_SHARPENER						//
//////////////////////////////////////////////////////////////////////////

//#define USE_MARK_FOR_SUBDIVISION

int Basic_TMesh::featureRecover(double ang, double edge_mult)
{
 Triangle *t;
 Edge *e;
 Vertex *v;
 Node *n;
 Point p;
 List splitPoints;
#ifdef USE_MARK_FOR_SUBDIVISION
 List fvs;
#endif

 TMesh::begin_progress();
 if (ang == 0)
 {
  FOREACHEDGE(e, n) ang += e->curvature();
  ang *= 1.5;
  ang /= E.numels();
  TMesh::info("threshold angle: %f (automatic setting)\n",ang);
 }
 else TMesh::info("threshold angle: %f (user defined)\n",ang);

 filter0(&E, ang); TMesh::report_progress("5 %% done   ");
 filter1(&V); TMesh::report_progress("5 %% done   ");
 filter2(&T); TMesh::report_progress("10 %% done   ");
 filter3(&T); TMesh::report_progress("15 %% done   ");
 filter4(&T); TMesh::report_progress("20 %% done   ");
 filter5(&E); TMesh::report_progress("25 %% done   ");
 filter6(&T); TMesh::report_progress("30 %% done   ");

 if (edge_mult < 0) return 0;

#ifdef USE_MARK_FOR_SUBDIVISION
 FOREACHEDGE(e, n)
  if (!e->isOnBoundary() && ES_COLOR(e->t1) == ES_RED && 
	ES_COLOR(e->t2) == ES_RED && (e->curvature() >= ang))
   {fvs.appendHead(e->v1); fvs.appendHead(e->v2);}
 FOREACHTRIANGLE(t, n)
 {
  if (ES_COLOR(t->e1) == ES_BLUE && ES_COLOR(t->e2) == ES_RED && ES_COLOR(t->e3) == ES_RED)
   fvs.appendHead(t->oppositeVertex(t->e1));
  if (ES_COLOR(t->e2) == ES_BLUE && ES_COLOR(t->e3) == ES_RED && ES_COLOR(t->e1) == ES_RED)
   fvs.appendHead(t->oppositeVertex(t->e2));
  if (ES_COLOR(t->e3) == ES_BLUE && ES_COLOR(t->e1) == ES_RED && ES_COLOR(t->e2) == ES_RED)
   fvs.appendHead(t->oppositeVertex(t->e3));
 }
#endif

 double ael, mel = 0.0;
 FOREACHEDGE(e, n) if ((ael = e->length()) > mel) mel=ael;
 TMesh::report_progress("35 %% done   ");

 FOREACHTRIANGLE(t, n)
 {
  p = getTriangleToSplit(t, edge_mult*mel);
  if (p.x != DBL_MAX) splitPoints.appendHead(new splitCorner(p, t));
 }
 TMesh::report_progress("45 %% done   ");

 // Perform the actual triangle split
 FOREACHNODE(splitPoints, n)
#ifdef USE_MARK_FOR_SUBDIVISION
  fvs.appendHead(splitTriangle((((splitCorner *)n->data)->t), &(((splitCorner *)n->data)->p)));
#else
  splitTriangle((((splitCorner *)n->data)->t), &(((splitCorner *)n->data)->p));
#endif
 // Free the triangle list
 splitPoints.freeNodes();
 TMesh::report_progress("50 %% done   ");

 FOREACHEDGE(e, n)
 {
  p = getPointToSplit(e, edge_mult*mel);
  if (p.x != DBL_MAX) splitPoints.appendHead(new splitPoint(p, e));
 }
 TMesh::report_progress("60 %% done   ");

 // Perform the actual split
 FOREACHNODE(splitPoints, n)
#ifdef USE_MARK_FOR_SUBDIVISION
  fvs.appendHead(splitEdge((((splitPoint *)n->data)->e), &(((splitPoint *)n->data)->p)));
#else
  splitEdge((((splitPoint *)n->data)->e), &(((splitPoint *)n->data)->p));
#endif

 // Free the edge list
 splitPoints.freeNodes();
 TMesh::report_progress("70 %% done   ");

 FOREACHTRIANGLE(t, n) {ES_PAINT(t, ES_GREY); UNMARK_VISIT(t);}
 TMesh::report_progress("80 %% done   ");
 FOREACHEDGE(e, n) ES_PAINT(e, ES_GREY);
 TMesh::report_progress("90 %% done   ");
 FOREACHVERTEX(v, n) ES_PAINT(v, ES_GREY);
 TMesh::report_progress("95 %% done   ");

#ifdef USE_MARK_FOR_SUBDIVISION
 for (n=fvs.head(); n != NULL; n=n->next)
 {
  v = (Vertex *)n->data;
  MARK_BIT(v,1);
 }
 // Actually, the following should mark only feature edges ...
 FOREACHEDGE(e, n) if (IS_BIT(e->v1, 1) && IS_BIT(e->v2, 1)) TAG_SHARPEDGE(e);
 for (n=fvs.head(); n != NULL; n=n->next)
 {
  v = (Vertex *)n->data;
  UNMARK_BIT(v,1);
 }
#endif
 TMesh::end_progress();

 return 1;
}
