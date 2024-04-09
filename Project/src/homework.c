
 #include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

void femElasticityAssembleElements(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
  int nLocal = theMesh->nLocalNode;
  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double gx = theProblem->gx;
  double gy = theProblem->gy;
  double **A = theSystem->A;
  double *B = theSystem->B;
  

  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }

    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      double xsi = theRule->xsi[iInteg];
      double eta = theRule->eta[iInteg];
      double weight = theRule->weight[iInteg];
      femDiscretePhi2(theSpace, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      double r = 0.0;
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];   // dxdxsi,dxdeta,dydxsi,dydeta: dérivées spatiales de x et y
        dxdeta += x[i] * dphideta[i];   // dphidxsi, dphideta, dphidxsi, dphideta : dérivées des fonctions de forme
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0)
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
      jac = fabs(jac);

      for (i = 0; i < theSpace->n; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      for (i = 0; i < theSpace->n; i++) {
        for (j = 0; j < theSpace->n; j++) {
          A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
          A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
          A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
          A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
        }
      }
      for (i = 0; i < theSpace->n; i++) {
        B[mapX[i]] += phi[i] * gx * rho * jac * weight;
        B[mapY[i]] += phi[i] * gy * rho * jac * weight;
      }
    }
  }
}






void femElasticityAssembleNeumann(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->ruleEdge;
  femDiscrete *theSpace = theProblem->spaceEdge;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;
  double x[2], y[2], phi[2];
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2];
  int nLocal = 2;
  double *B = theSystem->B;

  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;

    if(type != NEUMANN_X && type != NEUMANN_Y && type != NEUMANN_N && type != NEUMANN_T && type != NEUMANN_HYDROSTAT){
      continue;
    }
//////// Le but ici est de connaitre la hauteur totale de la stucture. 
//////// On va parcourir les positions 'y' des noeuds pour avoir la positions min et max.

    int iEl;//utiliser d'autres variables pour ne pas ecraser les variables de boucle
    int hauteur;// calcule de la hauteur du domaine en question
    int y_max = -1;//
    int y_min = -1;//
    int y_moy;//distance au centre d'une arete depuis l'origine
    int POS; // position correcte avec 0 au dessus du barrage, donc négative. Correspond à y dans rho*g*y
    int number_elem = theCondition->domain->nElem;
    for (int iEd = 0; iEd < theCondition->domain->nElem; iEd++) {
      iEl = theCondition->domain->elem[iEd];
      for (int k = 0; k < nLocal; k++) {
        map[k] = theEdges->elem[iEl * nLocal + k];
        y[k] = theNodes->Y[map[k]];
        if(y_max == -1 && y_min == -1){
          y_max = y[k];
          y_min = y[k];
        }
        if(y[k] > y_max){
          y_max = y[k];
        }
        if(y[k] < y_min){
          y_min = y[k];
        }
      }
    }

    hauteur = y_max - y_min;//
////////////////////////////////////////////
////// Maintenant on applique la condition de Neumann sur chaque morceaux de chaque domaine 
    
    for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
      iElem = theCondition->domain->elem[iEdge];
      for (j = 0; j < nLocal; j++) {
        map[j] = theEdges->elem[iElem * nLocal + j];
        x[j] = theNodes->X[map[j]];
        y[j] = theNodes->Y[map[j]];
        
      }
       
      int y_moy = (y[0]+y[1])/2; //
      POS = y_moy - hauteur;//

      double tx = x[1] - x[0];
      double ty = y[1] - y[0];
      double length = hypot(tx, ty);
      double jac = length / 2.0;
      
      double f_x = 0.0;
      double f_y = 0.0;
      if (type == NEUMANN_X) {
        f_x = value;
      }
      if (type == NEUMANN_Y) {
        f_y = value;
      }
      if (type == NEUMANN_N) {
        double nx =  ty / length;
        double ny = -tx / length;
        f_x = value * nx;
        f_y = value * ny;
      }
      if (type == NEUMANN_T) {
        f_x = value * tx/length;
        f_y = value * ty/length;
      }
      ///////
      if (type == NEUMANN_HYDROSTAT){
        f_x = value * POS;//
      }
      ///////


      //
      // A completer :-)
      // Attention, pour le normal tangent on calcule la normale (sortante) au SEGMENT, surtout PAS celle de constrainedNodes
      // Une petite aide pour le calcul de la normale :-)
      // double nx =  ty / length;
      // double ny = -tx / length;

      for (iInteg = 0; iInteg < theRule->n; iInteg++) {
        double xsi = theRule->xsi[iInteg];
        double weight = theRule->weight[iInteg];
        femDiscretePhi(theSpace, xsi, phi);
        for (i = 0; i < theSpace->n; i++) {
          B[2*map[i] + 0] += jac * weight * phi[i] * f_x;
          B[2*map[i] + 1] += jac * weight * phi[i] * f_y;
        }
      }
    }
  }
}

void femElasticityApplyDirichlet(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;

  for (int node = 0; node < theNodes->nNodes; node++) {
    femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
    if (theConstrainedNode->type == UNDEFINED)
      continue;
    femBoundaryType type = theConstrainedNode->type;

    if (type == DIRICHLET_X) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 0, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femFullSystemConstrain(theSystem, 2 * node + 0, value_x);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_y);
    }

    if (type == DIRICHLET_N) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer :-)
      femFullSystemConstrain(theSystem, 2 * node + 0, value * nx);
      femFullSystemConstrain(theSystem, 2 * node + 1, value * ny);
    }
    if (type == DIRICHLET_T) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;


      // A completer :-)
      double tx = -ny*value;
      double ty = nx*value;
      femFullSystemConstrain(theSystem, 2 * node + 0, tx);
      femFullSystemConstrain(theSystem, 2 * node + 1, ty);
    }
    if (type == DIRICHLET_NT) {
      double value_n = theConstrainedNode->value1;
      double value_t = theConstrainedNode->value2;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer :-)
      double tx = -ny;
      double ty = nx;
      femFullSystemConstrain(theSystem, 2 * node + 0, value_n * nx + value_t * tx); 
      femFullSystemConstrain(theSystem, 2 * node + 1, value_n * ny + value_t * ty);


    }
  }
}

double *femElasticitySolve(femProblem *theProblem) {
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  double *soluce = femFullSystemEliminate(theProblem->system);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}



void RenumberCuthill(femGeo *theGeometry) {
  int i;
  int *tab = malloc(sizeof(int) * theMesh->nodes->nNodes);

//avoir une liste avec les noeuds et leurs degres
femMesh *mesh = theGeometry->theElements;
int *elem_list = mesh->elem;

int list_size = mesh->nElem*mesh->nLocalNode*2
int *list = malloc(sizeof(int) * list_size);


  //étape 1 : créer tous les éléments i-j et j-i
for (int i = 0; i < (mesh->nElem); i++) {
  
    list[12*i]=malloc(2*sizeof(int));
    
    //on suppose nLocalNode = 3 sinon trop compliqué
    list[2*i*mesh->nLocalNode][0] = elem_list[i*mesh->nLocalNode];
    list[2*i*mesh->nLocalNode][1] = elem_list[i*mesh->nLocalNode+1];
    list[2*i*mesh->nLocalNode+1][0] = elem_list[i*mesh->nLocalNode];
    list[2*i*mesh->nLocalNode+1][1] = elem_list[i*mesh->nLocalNode+2];
    list[2*i*mesh->nLocalNode+2][0] = elem_list[i*mesh->nLocalNode+1];
    list[2*i*mesh->nLocalNode+2][1] = elem_list[i*mesh->nLocalNode];
    list[2*i*mesh->nLocalNode+3][0] = elem_list[i*mesh->nLocalNode+1];
    list[2*i*mesh->nLocalNode+3][1] = elem_list[i*mesh->nLocalNode+2];
    list[2*i*mesh->nLocalNode+4][0] = elem_list[i*mesh->nLocalNode+2];
    list[2*i*mesh->nLocalNode+4][1] = elem_list[i*mesh->nLocalNode];
    list[2*i*mesh->nLocalNode+5][0] = elem_list[i*mesh->nLocalNode+2];
    list[2*i*mesh->nLocalNode+5][1] = elem_list[i*mesh->nLocalNode+1];
  }
  //etape 2 : trier la liste 

  int compare(int *a, int *b) {

    if(a[0] == b[0])
      return (a[1] - b[1]);

    return (a[0] - b[0]);
  }

  qsort(list, list_size, sizeof(int*), compare)

  //etape 3 : supprimer les doublons   => bien faux mais hassoul

    int* new_list = malloc(n * sizeof(int*));
    int size = 0;
    new_list[j++] = list[0];
    for (int i = 1; i < list_size; i++) {
        if (list[i] != list[i - 1]) {
            new_list[size++] = list[i];
        }
    }
    

  //etape 4 : créer les vecteurs col et rptr 

  int *col = malloc(sizeof(int) * new_size); //new_size à instancier
  int *rptr = malloc(sizeof(int) * (theGeometry->theNodes->nNodes+1));
  rptr[0]=0;


  int num_col = 0;
  int num = 0;
  int noeud_actuel = list[0][0];

  for (int i = 0; i < theGeometry->theNodes->nNodes; i++) {
    
    while(list[num][0] == noeud_actuel){
      
      col[num_col++] = list[num++][1];}
    
    rptr[i] = num_col;
    noeud_actuel = list[num][0];
    


  }


  






//Instantiate an empty queue Q and empty array for permutation order of the objects R.
  int q = malloc(sizeof(int) todo);
  int r = malloc(sizeof(int) todo);
  int next_empty_q = 0;
  int next_empty_r = 0;
  int not_visited = malloc(sizeof(int) todo);

//S1: We first find the object with minimum degree whose index has not yet been added to R. 
//    Say, object corresponding to pth row has been identified as the object with a minimum degree. Add p to R.

  int min_index = todo;
  r[next_empty_r] = min_index;
  next_empty_r++;


//S2: As an index is added to R, and add all neighbors of the corresponding object at the index, 
//    in increasing order of degree, to Q. The neighbors are nodes with non-zero value amongst the 
//    non-diagonal elements in the pth row.

//S3: Extract the first node in Q, say C. If C has not been inserted in R, add it to R, add to Q the neighbors 
//    of C in increasing order of degree.
//S4: If Q is not empty, repeat S3.
//S5: If Q is empty, but there are objects in the matrix which have not been included in R, start from S1, once again. (This could happen if there are disjoint graphs)
//S6: Terminate this algorithm once all objects are included in R.
//S7: Finally, reverse the indices in R, i.e. (swap(R[i], R[P-i+1])).
  
  
  
 }