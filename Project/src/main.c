/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

int main(void) {
  femGeo *theGeometry = geoGetGeometry();
  geoMeshRead("../../ProjectPreProcessor/data/mesh.txt");
  femProblem *theProblem = femElasticityRead(theGeometry, "../../ProjectPreProcessor/data/problem.txt");
  femElasticityPrint(theProblem);
  double *theSoluce = femElasticitySolve(theProblem);
  int nNodes = theGeometry->theNodes->nNodes;
  femSolutionWrite(nNodes, 2, theSoluce, "../../Project/data/UV.txt");
  femElasticityFree(theProblem);
  geoFree();
  return 0;
}
