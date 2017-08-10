/****************************************************************************
* JMeshExt                                                                  *
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

#include "tin.h"
#include <stdio.h>
#include <stdlib.h>

using namespace T_MESH;


Point sharpLaplacianDisplacement(Vertex *v)
{
 if (IS_SHARPEDGE(v)) return *v;

 List *ve = v->VE();
 Vertex *w;
 Node *m;
 Edge *e;
 Point np;
 int nse=0;

 FOREACHVEEDGE(ve, e, m)
  if (IS_SHARPEDGE(e) || e->isOnBoundary())
  {
   if (nse==0) np.setValue(e->oppositeVertex(v));
   else if (nse==1) np = np + (*(e->oppositeVertex(v)));
   else {delete(ve); return (*v);}
   nse++;
  }
  else if (!nse) {w = e->oppositeVertex(v); np = np+(*w);}

 if (!nse) np = (np)/(ve->numels());
 else if (nse == 1) np = (*v);
 else np = np/2;

 delete(ve);
 return np;
}

void Basic_TMesh::laplacianSmooth(int ns, coord l)
{
 Triangle *t;
 Edge *e;
 Vertex *v;
 Node *n;
 int i = 0, is_selection = 0, ins = ns;
 coord ln = coord(1.0)-l;
 Point np;

 FOREACHEDGE(e, n) UNMARK_VISIT(e);
 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
  {MARK_VISIT(t->e1); MARK_VISIT(t->e2); MARK_VISIT(t->e3);}
 FOREACHEDGE(e, n) if (IS_VISITED(e))
  {MARK_VISIT(e->v1); MARK_VISIT(e->v2); is_selection = 1;}

 List vts;
 FOREACHVERTEX(v, n) if (!is_selection || IS_VISITED(v)) vts.appendHead(v);

 double *xyz = (double *)malloc(sizeof(double)*vts.numels() * 3);
 if (xyz == NULL) {TMesh::warning("Not enough memory for vertex coordinates.\n"); return;}

 TMesh::begin_progress();
 for (; ns>0; ns--)
 {
  i=0;
  FOREACHVVVERTEX((&vts), v, n)
  {
   np = sharpLaplacianDisplacement(v);
   if (!(i%3000)) TMesh::report_progress("%d %% done - %d steps left",((i*33)/(vts.numels()) + (100*(ins-ns)))/ins, ns);
   xyz[i++] = TMESH_TO_DOUBLE(np.x*l + v->x*ln); xyz[i++] = TMESH_TO_DOUBLE(np.y*l + v->y*ln); xyz[i++] = TMESH_TO_DOUBLE(np.z*l + v->z*ln);
  }

  i=0;
  FOREACHVVVERTEX((&vts), v, n)
  {
   v->x = xyz[i++]; v->y = xyz[i++]; v->z = xyz[i++];
  }
 }
 TMesh::end_progress();
 free(xyz);

 FOREACHEDGE(e, n) UNMARK_VISIT(e);
 FOREACHVERTEX(v, n) UNMARK_VISIT(v);
}


/////////////////////////////////////////////////////////////////////////
//
// BILATERAL MESH DENOSING (only for whole mesh)
//
/////////////////////////////////////////////////////////////////////////

Point bilateralDenoisedPoint(Vertex *v, double sigma_c, double sigma_s)
{
 List *ve = v->VE();
 Vertex *q;
 Node *m;
 Edge *e;
 Point np;
 Point nor = v->getNormal();
 coord crdh;
 double w_c, w_s, t, h, sum = 0.0, normalizer = 0.0;

 FOREACHVEEDGE(ve, e, m)
 {
  q = e->oppositeVertex(v);
  t = ((*v)-(*q)).length();
  crdh = nor*((*v) - (*q)); h = TMESH_TO_DOUBLE(crdh);
  w_c = exp(-(t*t)/(2.0*sigma_c*sigma_c));
  w_s = exp(-(h*h)/(2.0*sigma_s*sigma_s));
  sum += w_c*w_s*h;
  normalizer += w_c*w_s;
 }

 np = (*v) - nor*(sum/normalizer);

 delete(ve);
 return np;
}

void Basic_TMesh::bilateralSmooth(int ns)
{
 Vertex *v;
 Node *n;
 int i;
 Point np, nor;

 Triangle *t;
 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
 {
  TMesh::warning("bilateralSmooth: Selections have no special treatment. Smoothing whole mesh.");
  break;
 }

 Edge *e;
 FOREACHEDGE(e, n) if (IS_SHARPEDGE(e))
 {
  TMesh::warning("bilateralSmooth: Sharp edges have no special treatment. Smoothing whole mesh.");
  break;
 }

 // Compute sigma_c as the average length of mesh edges
 double sigma_c = 0.0;
 FOREACHEDGE(e, n) sigma_c += e->length();
 sigma_c /= E.numels();

 // Compute all the vertex normals
 Point *nors = new Point[V.numels()];
 i=0; FOREACHVERTEX(v, n) {nors[i]=v->getNormal(); v->info = nors+i; i++;}

 // Compute average offset
 coord av_offset = 0.0;
 FOREACHEDGE(e, n)
 {
  np = (*(e->v2)) - (*(e->v1));
  nor = *((Point *)e->v2->info);
  av_offset += np*nor;
  np = (*(e->v1)) - (*(e->v2));
  nor = *((Point *)e->v1->info);
  av_offset += np*nor;
 }
 av_offset /= (E.numels()*2);

 // Compute sigma_s as standard deviation of offset
 coord ps, sigma_s = 0.0;
 FOREACHEDGE(e, n)
 {
  np = (*(e->v2)) - (*(e->v1));
  nor = *((Point *)e->v2->info);
  ps = av_offset - np*nor; sigma_s += ps*ps;
  np = (*(e->v1)) - (*(e->v2));
  nor = *((Point *)e->v1->info);
  ps = av_offset - np*nor; sigma_s += ps*ps;
 }
 sigma_s = sqrt(TMESH_TO_DOUBLE(sigma_s/(E.numels()*2)));

 delete [] nors;

 Point *dps = new Point[V.numels()];

 TMesh::begin_progress();
 for (; ns>0; ns--)
 {
  i=0; FOREACHVERTEX(v, n)
  {
   dps[i] = bilateralDenoisedPoint(v, sigma_c, TMESH_TO_DOUBLE(sigma_s));
   if (!(i%3000)) TMesh::report_progress("%d %% done - %d steps left",(i*100)/V.numels(), ns);
   i++;
  }
  i=0; FOREACHVERTEX(v, n) v->setValue(dps[i++]);
 }
 TMesh::end_progress();
 delete [] dps;
}
