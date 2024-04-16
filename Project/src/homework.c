#include "fem.h"
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

int *DEGREE;
// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

void femElasticityAssembleElements(femProblem *theProblem) {
  femBandSystem *theSystem = theProblem->system;///ajout
  //femFullSystem *theSystem = theProblem->system;
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
        r += x[i]*phi[i];
      }
     
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0)
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
      jac = fabs(jac);

      for (i = 0; i < theSpace->n; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      if (theProblem->planarStrainStress == AXISYM){
        for (i = 0; i < theSpace->n; i++) {
          for (j = 0; j < theSpace->n; j++) {
            // on assemble la matrice de raideur qui sera symétrique donc que la partie triangulaire supérieure 
            A[mapX[i]][mapX[j]] += (mapX[j] >= mapX[i]) ? (dphidx[i] * a * r * dphidx[j] + dphidy[i] * c * r * dphidy[j] + phi[i] * ((b * dphidx[j]) + (a * phi[j] / r)) + dphidx[i] * b * phi[j]) * jac * weight : 0.0;
            A[mapX[i]][mapY[j]] += (mapY[j] >= mapX[i]) ? (dphidx[i] * b * r * dphidy[j] + dphidy[i] * c * r * dphidx[j] + phi[i] * b * dphidy[j]) * jac * weight : 0.0;
            A[mapY[i]][mapX[j]] += (mapX[j] >= mapY[i]) ? (dphidy[i] * b * r * dphidx[j] + dphidx[i] * c * r * dphidy[j] + dphidy[i] * b * phi[j]) * jac * weight : 0.0;
            A[mapY[i]][mapY[j]] += (mapY[j] >= mapY[i]) ? (dphidy[i] * a * r * dphidy[j] + dphidx[i] * c * r * dphidx[j]) * jac * weight : 0.0;
          }
        }
        for (i = 0; i < theSpace->n; i++) {
          B[mapX[i]] += phi[i] * gx * rho * jac * weight*r;
          B[mapY[i]] += phi[i] * gy * rho * jac * weight*r;
        }
      }else{
      for (i = 0; i < theSpace->n; i++) {
          for (j = 0; j < theSpace->n; j++) {
            // on assemble la matrice de raideur qui sera symétrique donc que la partie triangulaire supérieure 
            A[mapX[i]][mapX[j]] += (mapX[j] >= mapX[i]) ? (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight : 0.0;
            A[mapX[i]][mapY[j]] += (mapY[j] >= mapX[i]) ? (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight : 0.0;
            A[mapY[i]][mapX[j]] += (mapX[j] >= mapY[i]) ? (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight : 0.0;
            A[mapY[i]][mapY[j]] += (mapY[j] >= mapY[i]) ? (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight : 0.0;
          }
        }
      }
      
      for (i = 0; i < theSpace->n; i++) {
        B[mapX[i]] += phi[i] * gx * rho * jac * weight;
        B[mapY[i]] += phi[i] * gy * rho * jac * weight;
      }
    }
  }
}

//////// Le but ici est de connaitre la hauteur totale de la stucture. 
//////// On va parcourir les positions 'y' des noeuds pour avoir la positions min et max.
double calcul_hauteur(femProblem *theProblem, int iBnd){
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2];
  double x[2], y[2], phi[2];
  int nLocal = 2;
  int y_max = -1;
  int y_min = -1;
  int y_moy;
  int POS;
  femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
  int iEl;
  int hauteur;
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

    hauteur = y_max - y_min;
  
  return hauteur;
}



void femElasticityAssembleNeumann(femProblem *theProblem) {
  femBandSystem *theSystem = theProblem->system;///ajout
  //femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->ruleEdge;
  femDiscrete *theSpace = theProblem->spaceEdge;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;
  double x[2], y[2], phi[2];
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2];
  int nLocal = 2;
  double *B = theSystem->B;
  double pos;
  double y_moy;
  double hauteur;

  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;

    if(type != NEUMANN_X && type != NEUMANN_Y && type != NEUMANN_N && type != NEUMANN_T && type != NEUMANN_HYDROSTAT){
      continue;
    }
    /*
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
*/
////////////////////////////////////////////
////// Maintenant on applique la condition de Neumann sur chaque morceaux de chaque domaine 
    hauteur = calcul_hauteur(theProblem, iBnd);

    for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
      iElem = theCondition->domain->elem[iEdge];
      for (j = 0; j < nLocal; j++) {
        map[j] = theEdges->elem[iElem * nLocal + j];
        x[j] = theNodes->X[map[j]];
        y[j] = theNodes->Y[map[j]];
        
      }
       
      y_moy = (y[0]+y[1])/2; //
      pos = y_moy - hauteur;//

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
        f_x = value * pos;//
      }
      ///////


      //
      // A completer :-)
      // Attention, pour le normal tangent on calcule la normale (sortante) au SEGMENT, surtout PAS celle de constrainedNodes
      // Une petite aide pour le calcul de la normale :-)
      // double nx =  ty / length;
      // double ny = -tx / length;
      double r = 0.0;
      for (iInteg = 0; iInteg < theRule->n; iInteg++) {
        double xsi = theRule->xsi[iInteg];
        double weight = theRule->weight[iInteg];
        femDiscretePhi(theSpace, xsi, phi);
        for(i = 0; i < theSpace->n; i++){r = x[i]*phi[i];}
        if (theProblem->planarStrainStress == AXISYM){
          for (i = 0; i < theSpace->n; i++) {
            B[2*map[i] + 0] += jac * weight * phi[i] * f_x*r;
            B[2*map[i] + 1] += jac * weight * phi[i] * f_y*r;
          }
        }else{
          for (i = 0; i < theSpace->n; i++) {
            B[2*map[i] + 0] += jac * weight * phi[i] * f_x;
            B[2*map[i] + 1] += jac * weight * phi[i] * f_y;
        }

        }
        }
      }
    }
  }


void femElasticityApplyDirichlet(femProblem *theProblem) {
  femBandSystem *theSystem = theProblem->system;///ajout
  //femFullSystem *theSystem = theProblem->system;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;

  for (int node = 0; node < theNodes->nNodes; node++) {
    femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
    if (theConstrainedNode->type == UNDEFINED)
      continue;
    femBoundaryType type = theConstrainedNode->type;

    if (type == DIRICHLET_X) {
      double value = theConstrainedNode->value1;
      femBandSystemConstrain(theSystem, 2 * node + 0, value);///ajout
      //femFullSystemConstrain(theSystem, 2 * node + 0, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femBandSystemConstrain(theSystem, 2 * node + 1, value);///ajout
      //femFullSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femBandSystemConstrain(theSystem, 2 * node + 0, value_x);///ajout
      femBandSystemConstrain(theSystem, 2 * node + 1, value_y);///ajout
      //femFullSystemConstrain(theSystem, 2 * node + 0, value_x);
      //femFullSystemConstrain(theSystem, 2 * node + 1, value_y);
    }

    if (type == DIRICHLET_N) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer :-)
      femBandSystemConstrain(theSystem, 2 * node + 0, value * nx);///ajout
      femBandSystemConstrain(theSystem, 2 * node + 1, value * ny);///ajout
      //femFullSystemConstrain(theSystem, 2 * node + 0, value * nx);
      //femFullSystemConstrain(theSystem, 2 * node + 1, value * ny);
    }
    if (type == DIRICHLET_T) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;


      // A completer :-)
      double tx = -ny*value;
      double ty = nx*value;
      femBandSystemConstrain(theSystem, 2 * node + 0, tx);///ajout
      femBandSystemConstrain(theSystem, 2 * node + 1, ty);///ajout
      //femFullSystemConstrain(theSystem, 2 * node + 0, tx);
      //femFullSystemConstrain(theSystem, 2 * node + 1, ty);
    }
    if (type == DIRICHLET_NT) {
      double value_n = theConstrainedNode->value1;
      double value_t = theConstrainedNode->value2;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      // A completer :-)
      double tx = -ny;
      double ty = nx;
      femBandSystemConstrain(theSystem, 2 * node + 0, value_n * nx + value_t * tx);///ajout
      femBandSystemConstrain(theSystem, 2 * node + 1, value_n * ny + value_t * ty);///ajout
      //femFullSystemConstrain(theSystem, 2 * node + 0, value_n * nx + value_t * tx); 
      //femFullSystemConstrain(theSystem, 2 * node + 1, value_n * ny + value_t * ty);


    }
  }
}

# define MAX(a, b) ((a > b) ? (a) : (b)) // Pour comprendre les macros
# define MIN(a, b) ((a < b) ? (a) : (b)) // cf. notre précédent solutionnaire

int femMeshComputeBand(femMesh *theMesh) {

    int myMax, myMin, myBand, map[4];
    int nLocal = theMesh->nLocalNode;
    myBand = 0;

    for (int iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (int j = 0; j < nLocal; ++j)
            map[j] = theMesh->nodes->number[theMesh->elem[iElem * nLocal + j]];

        // On trouve le noeud maximum et minimum
        myMin = map[0];
        myMax = map[0];
        for (int j = 1; j < nLocal; j++) {
            myMax = MAX(map[j], myMax);
            myMin = MIN(map[j], myMin);
        }

        if (myBand < (myMax - myMin))
            myBand = myMax - myMin;
    }

    return (++myBand); // On incrémente de 1 myBand avant de le renvoyer (formule)
}

// assemblage de la  bande 
void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    for (int i = 0; i < nLoc; ++i){
        int row = map[i]; // on saute partout dans la matrice en prenant toutes les combinaisons des différentes indices
        // possibles des noeuds d'un element (pex pour 0,3,2 on a 9 possibilités).
        for (int j = 0; j < nLoc; ++j){
            int column = map[j];
            //on veut assembler la matrice mais que la partie supérieure donc on met le if ici pour s'en assurer
            if (column >= row){
                myBandSystem->A[map[i]][map[j]] += Aloc[i*nLoc + j];
            }
        }
        myBandSystem->B[map[i]] += Bloc[i];

    }
}

double *femBandSystemEliminate(femBandSystem *myBand) {
  double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    // on adapte la fonction de l'élimination de gauss et back substitution de fem.c pour qu'elle fonctionne avec une bande

    for (int k = 0; k < size; ++k){
        if (fabs(A[k][k]) <= 1e-8){ Error("Cannot eliminate with such a pivot.");}
        //on prend ici en compte la taille de bande 
        jend = (k+band < size) ? k+band : size; // on fait un min car quand on arrive à la fin de la matrice, la diagonale
        // se rapprochera de plus en plus de la dernière colonne et diminue la taille de la bande, il y a des risques de
        // débordement de mémoire.
        for (int i = k+1; i < jend; ++i) {
            factor = A[k][i] / A[k][k];
            for (int j = i; j < jend; ++j) {
                A[i][j] -= A[k][j] * factor;
            }
            B[i] -= B[k] * factor;
        }
    }
    for (int i = size-1; i >= 0; --i) {
        factor = 0;
        jend = (i+band < size) ? i+band : size;
        for (int j = i+1; j < jend; ++j) {
            //la matrice est symmétrique donc on échange les indices pour toujours travailler dans la partie supérieure
            factor += A[i][j] * B[j];
        }
        B[i] = (B[i] - factor)/A[i][i];
    }
    return(myBand->B);
}






 int compare(const void *a, const void *b){

    int x = *(int*)a;
    int y = *(int*)b;

    if(x[0] == y[0])
      return (x[1] - y[1]);
    return (x[0] - y[0]);
  }

  int compare2(const void *a, const void *b){
      int x = (int)a;
      int y = (int)b;
  return (DEGREE[x] - DEGREE[y]);
  }

int* RenumberCuthill(femGeo *theGeometry) {
  printf("début Renumber Cuthill\n");
  //avoir une liste avec les noeuds et leurs degres
  femMesh *mesh = theGeometry->theElements;
  int *elem_list = mesh->elem;

  int list_size = mesh->nElem*mesh->nLocalNode*2;
  int *list = malloc(sizeof(int) * list_size);

  for (int i = 0; i < list_size; i++) {
      list[i] = malloc(2*sizeof(int));
    }


    //étape 1 : créer tous les éléments i-j et j-i
  for (int i = 0; i < (mesh->nElem); i++) {
    
      
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



  qsort(list, list_size, sizeof(int*), compare);

  //etxpe 3 : supprimer les doublons   => bien faux mais hassoul

    int *list_whithout_duplicates = malloc(list_size * sizeof(int));
    //for (int i = 0; i < list_size; i++) {
    //    list_whithout_duplicates[i] = malloc(2*sizeof(int));
    //  }
    printf("avant list_whithout_duplicate\n");
        

    int size = 1;
    list_whithout_duplicates[0] = list[0];
    for (int i = 1; i < list_size; i++) {
        if (list[i][0] != list[i - 1][0] || list[i][1] != list[i - 1][1]) {
            list_whithout_duplicates[size++] = list[i];
        }
    }
    printf("après list_whithout_duplicate\n");
    

  //etape 4 : créer les vecteurs col et rptr 

  int *col = malloc(sizeof(int) * size); 
  int *rptr = malloc(sizeof(int) * (theGeometry->theNodes->nNodes+1));
  int *degree = malloc(sizeof(int) * (theGeometry->theNodes->nNodes));

  rptr[0]=0;


  int num_col = 0;
  int num = 0;
  int noeud_actuel = list_whithout_duplicates[0][0];

  printf("intermédiaire\n");
  
  for (int i = 0; i < theGeometry->theNodes->nNodes; i++) {
    degree[i]=0;
    while(num < size-1 && list_whithout_duplicates[num][0] == noeud_actuel ){


      col[num_col++] = list_whithout_duplicates[num++][1]; //size-1 me parait bizarre 
      degree[i]++;}
    
    rptr[i] = num_col;
    noeud_actuel = list_whithout_duplicates[num][0];
    


  }
  printf("millieu Renumber Cuthill\n");

//Instantiate an empty queue Q and empty array for permutation order of the objects R.
  int *q = malloc(sizeof(int)*theGeometry->theNodes->nNodes*theGeometry->theNodes->nNodes);
  int *r = malloc(sizeof(int)*theGeometry->theNodes->nNodes);
  int next_empty_q = 0;
  int first_filled_q = 0;
  int next_empty_r = 0;
  bool *not_visited = malloc(sizeof(int)*theGeometry->theNodes->nNodes);

//S1: We first find the object with minimum degree whose index has not yet been added to R. 
//    Say, object corresponding to pth row has been identified as the object with a minimum degree. Add p to R.

  for (int i = 0; i < theGeometry->theNodes->nNodes; i++) {
    not_visited[i] = true;
  }

  int actual_degree = degree[0];
  int minimum_node = 0;
  bool change = true;

    for (int i = 0; i< theGeometry->theNodes->nNodes; i++) {
      if(degree[i] < actual_degree && not_visited[i]){
        actual_degree = degree[i];
        minimum_node = i;
      }
    }

  while(change){
    printf("------------début du while------------\n");

   //next_empty_r<theGeometry->theNodes->nNodes &&

      r[next_empty_r++] = minimum_node;
      not_visited[minimum_node] = false;


    //S2: As an index is added to R, and add all neighbors of the corresponding object at the index, 
    //    in increasing order of degree, to Q. The neighbors are nodes with non-zero value amongst the 
    //    non-diagonal elements in the pth row.
      
        for (int i = rptr[minimum_node]; i < rptr[minimum_node+1]; i++) {
          q[next_empty_q++] = col[i];
        }
        printf("next_empty_q : %d, size : %d \n", next_empty_q, theGeometry->theNodes->nNodes * theGeometry->theNodes->nNodes);
        printf("avant compare2 Renumber Cuthill\n");

      DEGREE = degree;
      qsort(q, next_empty_q, sizeof(int), compare2);
      printf("après compare2 Renumber Cuthill\n");




    //S3: Extract the first node in Q, say C. If C has not been inserted in R, add it to R, add to Q the neighbors 
    //    of C in increasing order of degree.
    //S4: If Q is not empty, repeat S3.
    int C;
        while(first_filled_q < next_empty_q){
          C = q[first_filled_q++];
          if(not_visited[C]){
            r[next_empty_r++] = C;
            not_visited[C] = false;
            int scission = next_empty_q;
            for (int j = rptr[C]; j < rptr[C+1]; j++) {
              q[next_empty_q++] = col[j];
            }
            qsort(&q[scission], next_empty_q-scission, sizeof(int), compare2);//il faut sort à partir de scission à mon avis ce que j'ai mis est faux 
        }
      
    }
    //S5: If Q is empty, but there are objects in the matrix which have not been included in R, start from S1, 
    //    once again. (This could happen if there are disjoint graphs)
    //S6: Terminate this algorithm once all objects are included in R.
    change = false;
    actual_degree = 10000;
    for (int i = 0; i< theGeometry->theNodes->nNodes; i++) {
    if(degree[i] < actual_degree && not_visited[i]){
      actual_degree = degree[i];
      minimum_node = i;
      change = true;
      }
    }
  }


    //S7: Finally, reverse the indices in R, i.e. (swap(R[i], R[P-i+1])).
  
    int temp;
    for (int i = 0; i < theGeometry->theNodes->nNodes / 2; i++) {
        temp = r[i];
        r[i] = r[theGeometry->theNodes->nNodes - i - 1];
        r[theGeometry->theNodes->nNodes - i - 1] = temp;
    }

    printf("avant de free Renumber Cuthill\n");
    // section free (j'ai été un peu vite la dessus)   
    for (int i = 0; i < list_size; i++) {
        free(list[i]);
    }

    free(list);

    printf("millieu free Renumber Cuthill\n");
    //for (int i = 0; i < list_size; i++) {
    //    free(list_whithout_duplicates[i]);
    //}
    printf("après for free Cuthill\n");
    //free(list_whithout_duplicates);

    free(col);
    free(rptr);
    printf("trois quart free Renumber Cuthill\n");
    free(degree);
    free(q);
    //free(r);
    free(not_visited);



    for (int i = 0; i < theGeometry->theNodes->nNodes; i++) {
        theGeometry->theNodes->number[i] = r[i];
    }
    
    printf("fin Renumber Cuthill\n");
    return r; 

  
 }

  double *femElasticitySolve(femProblem *theProblem) {
  int* array = RenumberCuthill(theProblem->geometry);
    int size = theProblem->geometry->theNodes->nNodes;
    for (int i = 0; i < size; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");

  //femElasticityAssembleElements(theProblem);

  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  femBandSystemAssemble(theProblem->system, theProblem->system->A, theProblem->system->B, array, size);
  //double *soluce = femFullSystemEliminate(theProblem->system);
  double *soluce = femBandSystemEliminate(theProblem->system);///ajout
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}