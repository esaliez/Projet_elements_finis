#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

double geoSize(double x, double y) {

  femGeo *theGeometry = geoGetGeometry();
  return theGeometry->h * (1.0 - 0.5 * x);
}

void geoMeshGenerate(void) {
    femGeo* theGeometry = geoGetGeometry();



    int ierr;

    int* draw(double * t, int n, double X(double), double Y(double)) {
    int ierr;
    int idn, ido;
    int * idl = (int *) malloc((n +1) * sizeof(int));
    idl[0] = gmshModelOccAddPoint(X(t[0]), Y(t[0]), 0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    ido = idl[0];
    for (int i = 1; i < n; i++) {
        idn = gmshModelOccAddPoint(X(t[i]), Y(t[i]), 0, 0, -1, &ierr);
        ErrorGmsh(ierr);
        idl[i + 1] = gmshModelOccAddLine(ido, idn, -1, &ierr);
        ErrorGmsh(ierr);
        ido = idn;
    }
    idl[1] = idn;
    ErrorGmsh(ierr);
    return idl;
}

double X(double t) {
    int pas = 20;
    double * x = (double *) malloc((pas) * sizeof(double));
    x[0] = 0;
    x[1] = 2;
    x[2] = 0.5;
    x[3] = 0.5;
    x[4] = 0;
    int next=5;
    for (int i = next; i < pas; i++) {
        x[i] = 0;
    }
    return x[(int)t];
}
double Y(double t) {
    int pas = 20;
    double *y = (double *)malloc(pas * sizeof(double));
    y[0] = 0;
    y[1] = 0;
    y[2] = 1.5;
    y[3] = 2;
    y[4] = 2;
    int next = 5;
    for (int i = next; i < pas; i++) {
        y[i] = 2 - 2 * (i - next + 1) / (double)(pas - next + 1);
        printf("y[%d] = %f\n", i, y[i]);
    }

    return y[(int)t];
    
}

    int pas = 20;

    double * t = (double *) malloc((pas) * sizeof(double));
    for (int i = 0; i < pas; i++) {
        t[i] = i ;
    }
    
    int * idd = draw(t, pas, X, Y);
    int * ids = malloc(pas * sizeof(int));
    
    for (int i = 0; i < pas-1; i++) {
        ids[i + 1] = idd[i + 2];
    }
    printf("pass %d -- %d\n", idd[0], idd[1]);
    ids[0] = gmshModelOccAddLine(idd[0], idd[1], -1, &ierr);
    
    int curve = gmshModelOccAddCurveLoop(ids, pas, -1, &ierr);
    ErrorGmsh(ierr);
    int surface = gmshModelOccAddPlaneSurface(&curve, 1, -1, &ierr);
    ErrorGmsh(ierr);
    








 
//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    //geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
  if (theGeometry->elementType == FEM_QUAD) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
    gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  if (theGeometry->elementType == FEM_TRIANGLE) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

}

void geoMeshGenerateMshFile(const char *filename) {
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);
  return;
}

void geoMeshGenerateGeoFile(const char *filename) {
  femGeo *theGeometry = geoGetGeometry();
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);
  if (theGeometry->elementType == FEM_QUAD) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
    gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  if (theGeometry->elementType == FEM_TRIANGLE) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }
  return;}
