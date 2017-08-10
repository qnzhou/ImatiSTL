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
#include <faceClustering.h>
#include <simplification.h>


////////////////////////////////////////////////////////////////////
//
// HIERARCHICAL FACE CLUSTERING
//
////////////////////////////////////////////////////////////////////

double facecNode::max_L_inf = DBL_MAX;

int clusterHeap::compare(const void *e1, const void *e2)
{
	clusterEdge *a = edges[(j_voidint)e1];
	clusterEdge *b = edges[(j_voidint)e2];
	double l1 = a->cost;
	double l2 = b->cost;

	if (l1 < l2) return -1;
	if (l1 > l2) return 1;
	return 0;
}


clusterEdge *clusterGraph::createEdge(graphNode *n1, graphNode *n2)
{
	Node *n;
	FOREACHNODE(n1->edges, n)
	if (((clusterEdge *)n->data)->hasNode(n2))
		return (clusterEdge *)n->data;

	if (curEdgeIndex >= maxNumEdges) return NULL;

	clusterEdge *ne = new clusterEdge(n1, n2, curEdgeIndex);
	edges.appendHead(ne);
	ces[curEdgeIndex] = ne;
	ne->cost = costFunction(ne->n1, ne->n2);
	ch->update(ne);
	curEdgeIndex++;

	return ne;
}


clusterEdge *clusterGraph::getFirstEdge()
{
	while (!ch->isEmpty() && ch->getFirst()->isUnlinked()) ch->popHead();
	return ch->getFirst();
}

int clusterGraph::collapseFirstEdge(void(*mergenodes)(const void *, const void *))
{
	clusterEdge *e;

	while ((e = ch->popHead()) != NULL) if (!e->isUnlinked()) break;
	if (e == NULL) return 0;

	if (mergenodes != NULL) mergenodes(e->n1, e->n2);

	graphNode *gn = e->n1;
	e->collapse();
	Node *n;
	clusterEdge *ne;
	FOREACHNODE(gn->edges, n)
	{
		ne = (clusterEdge *)n->data;
		ne->cost = costFunction(ne->n1, ne->n2);
		ch->update(ne);
	}
	return 1;
}

facecNode::facecNode(Triangle *t)
{
 triangles.appendHead(t);
 double *a = new double(t->area());
 areas.appendHead(a);
 tot_area = *a;
 sum_ctr = t->getCenter();
 Vertex *v1=t->v1(), *v2=t->v2(), *v3=t->v3();
 Cov_v = (SymMatrix3x3(v1->x, v1->y, v1->z)+SymMatrix3x3(v2->x, v2->y, v2->z)+SymMatrix3x3(v3->x, v3->y, v3->z))*(tot_area/3.0);
 sum_ctr *= tot_area;
}

facecNode::facecNode()
{
 tot_area = 0;
}


void facecNode::merge(const void *n1, const void *n2)
{
 facecNode *c1 = (facecNode *)n1;
 facecNode *c2 = (facecNode *)n2;

 c1->triangles.joinTailList(&(c2->triangles));
 c1->areas.joinTailList(&(c2->areas));
 c1->sum_ctr += c2->sum_ctr;
 c1->tot_area += c2->tot_area;
 c1->Cov_v += c2->Cov_v;
}


double facecNode::edgeCostFunction(const void *cn1, const void *cn2)
{
// The ACT_AREA_BIAS is a small coefficient that prevents the
// creation of too unbalanced trees on flat areas
 const double ACT_AREA_BIAS = 1.0e-12;
 facecNode *n1 = (facecNode *)cn1;
 facecNode *n2 = (facecNode *)cn2;
 double tot_area = n1->tot_area + n2->tot_area;

 double cost = fittingPlaneCost(cn1, cn2);
 if (cost < DBL_MAX) cost += tot_area*ACT_AREA_BIAS;
 return cost;
}


double facecNode::fittingPlaneCost(const void *cn1, const void *cn2)
{
 facecNode *n1 = (facecNode *)cn1;
 facecNode *n2 = (facecNode *)cn2;

 double lcost = 0;
 Node *n, *m;
 Triangle *t;
 facecNode *gn;
 double area, tot_area = (n1->tot_area+n2->tot_area);
 Point app = (n1->sum_ctr+n2->sum_ctr)/tot_area;

 SymMatrix3x3 TM = n1->Cov_v+n2->Cov_v;
 TM -= (SymMatrix3x3(app.x, app.y, app.z)*tot_area);

 double norx, nory, norz;
 TM.getMinEigenvector(&norx, &nory, &norz);
 Point nor(norx, nory, norz);
 if (nor.isNull()) return DBL_MAX;

 coord d, a = -(nor*app);

 for (gn=n1; gn!=NULL; gn=(gn==n1)?(n2):(NULL))
 {
  m=gn->areas.head(); FOREACHVTTRIANGLE((&(gn->triangles)), t, n)
  {
   area = (*((double *)m->data));
   d = ((t->getCenter())*(nor))+a; d *= d; lcost += (TMESH_TO_DOUBLE(d)*area);
   // Probably the following could be slightly optimized
   d = ((*(t->v1()))*(nor))+a; if (d*d > max_L_inf) return DBL_MAX;
   d = ((*(t->v2()))*(nor))+a; if (d*d > max_L_inf) return DBL_MAX;
   d = ((*(t->v3()))*(nor))+a; if (d*d > max_L_inf) return DBL_MAX;
   m=m->next();
  }
 }

 return lcost;
}



////////////////////////////////////////////////////////////////////
//
// LLOYD RELAXATION
//
////////////////////////////////////////////////////////////////////

lloydCluster::lloydCluster(Triangle *t, int index)
{
 id = index;
 triangles.appendHead(t);
}

lloydCluster::lloydCluster(lloydCluster *a, lloydCluster *b, int index)
{
 id = index;
 triangles.appendList(&(a->triangles));
 triangles.appendList(&(b->triangles));
 update(false);
}

lloydCluster::lloydCluster(facecNode *fc, int index)
{
 id = index;
 triangles.appendList(&(fc->triangles));
}



void lloydCluster::update(bool removetriangles)
{
 Triangle *t;
 Vertex *v1, *v2, *v3;
 Node *n;
 double area, totarea = 0;
 SymMatrix3x3 Cov_v;
 Point sum_ctr;

 FOREACHNODE(triangles, n)
 {
  t = (Triangle *)n->data;
  t->info = NULL;
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  area = t->area();
  totarea += area;
  Cov_v += ((SymMatrix3x3(v1->x, v1->y, v1->z) + SymMatrix3x3(v2->x, v2->y, v2->z) + SymMatrix3x3(v3->x, v3->y, v3->z))*(area/3.0));
  sum_ctr += t->getCenter()*area;
 }

 sum_ctr /= totarea;

 Cov_v -= (SymMatrix3x3(sum_ctr.x, sum_ctr.y, sum_ctr.z)*totarea);
 double norx, nory, norz;
 Cov_v.getMinEigenvector(&norx, &nory, &norz);
 proxy_normal.setValue(norx, nory, norz);
 if (proxy_normal.isNull()) TMesh::error("Unexpected zero-length proxy normal!\n");

 coord pdcrd = -(proxy_normal*sum_ctr);
 proxy_d = TMESH_TO_DOUBLE(pdcrd);

 double acost, mincost=DBL_MAX;
 seed=NULL;
 FOREACHNODE(triangles, n)
 {
  t = (Triangle *)n->data;
  if ((acost=getCost(t))<mincost) {seed=t; mincost=acost;}
 }
 if (seed==NULL) seed = (Triangle *)triangles.head()->data;
 seed->info = this;

 if (removetriangles)
 {
  triangles.removeNodes();
  triangles.appendHead(seed);
 }
}


double lloydCluster::getCost(Triangle *t)
{
 Point proxy_center = proxy_normal*(-proxy_d);
 coord d1 = ((*(t->v1())) - proxy_center)*proxy_normal; d1 = FABS(d1);
 coord d2 = ((*(t->v2())) - proxy_center)*proxy_normal; d2 = FABS(d2);
 coord d3 = ((*(t->v3())) - proxy_center)*proxy_normal; d3 = FABS(d3);
 coord d = d1*d1 + d2*d2 + d3*d3 + d1*d2 + d1*d3 + d2*d3;
 coord ret = d*t->area();
 return TMESH_TO_DOUBLE(ret);
}

double lloydCluster::getTotalCost()
{
 Node *n;
 double cost=0.0;
 FOREACHNODE(triangles, n) cost += getCost((Triangle *)n->data);
 return cost;
}


int triangleHeap::compare(const void *a, const void *b)
{
 double i1 = ((clusteredTriangle *)a)->cost;
 double i2 = ((clusteredTriangle *)b)->cost;

 if (i1 < i2) return -1;
 if (i2 < i1) return 1;

 return 0;
}



//////////////////////////////////////////////////////////////////////////
//
// Main function to perform Lloyd relaxation
//
//////////////////////////////////////////////////////////////////////////

void etm_lloyd_optimization(Basic_TMesh *tin, List *fcclusters)
{
 Node *n;
 unsigned int i;
 Triangle *t, *y;
 lloydCluster *lc;
 clusteredTriangle *ct;
 List clusters;

 // Create initial clusters
 i=0;
 for (n=fcclusters->head(); n!=NULL; n=n->next())
	 if (((facecNode *)n->data)->triangles.numels()) clusters.appendTail(new lloydCluster((facecNode *)n->data, i++));
 TMesh::info("Relaxing %d clusters\n",clusters.numels());

 // Cluster have been created. Now start Lloyd relaxation
 triangleHeap th(tin->T.numels()*3);
 int iterations = 0;

 TMesh::begin_progress();
 while(clusters.numels() > 1 && iterations < 10)
 {
  TMesh::report_progress("Iteration %d of %d",iterations++);
  // Update
  FOREACHNODE(clusters, n) ((lloydCluster *)n->data)->update();

  FOREACHNODE(clusters, n)
  {
   lc = ((lloydCluster *)n->data);
   t = lc->seed;
   y=t->t1(); if (y!=NULL && y->info == NULL) th.push(new clusteredTriangle(y, lc));
   y=t->t2(); if (y!=NULL && y->info == NULL) th.push(new clusteredTriangle(y, lc));
   y=t->t3(); if (y!=NULL && y->info == NULL) th.push(new clusteredTriangle(y, lc));
  }

  // Flood
  while (!th.isEmpty())
  {
   ct = (clusteredTriangle *)th.popHead();
   t = ct->t;
   if (t->info == NULL)
   {
    t->info = ct->l;
    ct->l->triangles.appendHead(ct->t);
    y=t->t1(); if (y!=NULL && y->info == NULL) th.push(new clusteredTriangle(y, ct->l));
    y=t->t2(); if (y!=NULL && y->info == NULL) th.push(new clusteredTriangle(y, ct->l));
    y=t->t3(); if (y!=NULL && y->info == NULL) th.push(new clusteredTriangle(y, ct->l));
   }
   delete(ct);
  }
 }
 TMesh::end_progress();

 FOREACHNODE(clusters, n) delete((lloydCluster *)n->data);
 clusters.removeNodes();
}

//////////////////////////////////////////////////////////////////////////
//
// MAIN FUNCTION TO TAG REGION BOUNDARIES
//
//////////////////////////////////////////////////////////////////////////

void Basic_TMesh::tagPlanarRegionsBoundaries(double max_distance)
{
 facecNode::max_L_inf = max_distance*max_distance;

 // Init face graph
 clusterGraph cg(E.numels(), &facecNode::edgeCostFunction);

 Triangle *t;
 Edge *e;
 Node *n;

 FOREACHTRIANGLE(t, n) t->info = cg.addNode(new facecNode(t));

 FOREACHEDGE(e, n) if (!e->isOnBoundary())
  cg.createEdge(((facecNode *)e->t1->info),((facecNode *)e->t2->info));

 // Perform the clustering
 int i=0;
 TMesh::begin_progress();
 while (cg.getLowestCost()!=DBL_MAX && cg.collapseFirstEdge(&facecNode::merge))
 {
	 if (!(i%1000)) TMesh::report_progress("%d %% done", (i*100)/T.numels());
	 i++;
 }
 TMesh::end_progress();

 // Tag triangles with cluster pointers
 facecNode *fn;
 Node *m;
 FOREACHNODE(cg.nodes, n)
 {
  fn = (facecNode *)n->data;
  FOREACHVTTRIANGLE((&(fn->triangles)), t, m) t->info = fn;
 }

// etm_lloyd_optimization(this, &cg.nodes);

 // Mark edge faces
 FOREACHEDGE(e, n)
	 if (!e->isOnBoundary() && e->t1->info != e->t2->info) TAG_SHARPEDGE(e);

 FOREACHTRIANGLE(t, n) t->info = NULL;
}
