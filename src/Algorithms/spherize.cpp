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
#include "simplification.h"


class VC_Cost
{
 public:
 int index;
 double cost;
 double (*costFunction)(Vertex *);

 VC_Cost(Vertex *, int, double (*)(Vertex *));
 double updateCost(Vertex *);
};

VC_Cost::VC_Cost(Vertex *v, int ind, double (*cf)(Vertex *))
{
 costFunction = cf;
 index = ind;
 updateCost(v);
}

double VC_Cost::updateCost(Vertex *v)
{
 cost = costFunction(v); return cost;
}


class vertexHeap : abstractHeap
{
 public:
 Vertex **vertices;

 vertexHeap(int, Vertex **);
 ~vertexHeap();

 inline void push(Vertex *v) {insert((void *)((VC_Cost *)v->info)->index);}
 inline Vertex *popHead() {return vertices[(int)removeHead()];}
 inline int isEmpty() {return (numels==0);}

 void remove(Vertex *);
 void update(Vertex *);
 int compare(const void *, const void *);
};

vertexHeap::vertexHeap(int n, Vertex **vtx) : abstractHeap(n)
{
 positions = new int[n+1];
 vertices = vtx;
}

vertexHeap::~vertexHeap() {delete(positions);}

void vertexHeap::remove(Vertex *v)
{
 VC_Cost *vcc = ((VC_Cost *)v->info);
 vcc->cost = -1;
 int i=positions[vcc->index];
 if (i) {upheap(i); removeHead();}
}

void vertexHeap::update(Vertex *v)
{
 VC_Cost *vcc = ((VC_Cost *)v->info);
 vcc->updateCost(v);
 if (!positions[vcc->index]) push(v);
 else downheap(upheap(positions[vcc->index]));
}

int vertexHeap::compare(const void *e1, const void *e2)
{
 Vertex *a = vertices[(int)e1];
 Vertex *b = vertices[(int)e2];
 double l1 = ((VC_Cost *)a->info)->cost;
 double l2 = ((VC_Cost *)b->info)->cost;
 if (l1 < l2) return -1;
 if (l2 < l1) return 1;

 return 0;
}

Point sprz_getCenterOfMass(Vertex *v2)
{
 Point np;
 Vertex *v;
 Node *n;
 List *vv = v2->VV();
 FOREACHVVVERTEX(vv, v, n) np = np+(*v);
 np = np/vv->numels();
 delete(vv);

 return np;
}

void sprz_setCoordinates(Vertex *v2)
{
 Point np = sprz_getCenterOfMass(v2);
 v2->setValue(&np);
}

double sprz_costFunction(Vertex *v)
{
 double t = v->totalAngle();
 return (t>=0)?(t):(DBL_MAX);
}

int Basic_TMesh::spherize(int numits)
{
 if (numits <= 0) return 0;
 Node *n;
 Vertex *v, *ov;
 List *vv;
 int i;

 Vertex **vertices = new Vertex *[V.numels()];
 vertexHeap vh(V.numels(), vertices);

 i=0; FOREACHVERTEX(v, n)
 {
  vertices[i] = v; v->info = new VC_Cost(v, i++, sprz_costFunction); vh.push(v);
 }

 i=0;
 TMesh::begin_progress();
 while (i++ < numits)
 {
  v = vh.popHead();
  sprz_setCoordinates(v);
  vh.update(v);
  vv = v->VV();
  FOREACHVVVERTEX(vv, ov, n) vh.update(ov);
  delete(vv);
  TMesh::report_progress("%d%% done ", (100 * i) / numits);
 }
 TMesh::end_progress();

 delete(vertices);
 FOREACHVERTEX(v, n) {delete((VC_Cost *)v->info); v->info = NULL;}

 return 1;
}
