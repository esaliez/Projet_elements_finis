#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

double geoSize(double x, double y) {
  
  femGeo *theGeometry = geoGetGeometry();
  double h = theGeometry->h;
  int ref;
  ref = 5;
  double h1;
  if(x <= 15) {
    h = ref;
  }
  else if(x > 15 && x <= 50){
    h = ref + 2*ref*(x-15)/50;
    
  }
  else if(x > 50 && x <= 195){
    h = 12;
  }
  //return theGeometry->h * (1.0 - 0.5 * x);
  return h;
}

void geoMeshGenerate(void) {
  femGeo *theGeometry = geoGetGeometry();
  
  //theGeometry->h = Lx * 0.05;
  theGeometry->elementType = FEM_TRIANGLE;
  
  geoSetSizeCallback(geoSize);
  

  int ierr;
  //TRIANGLE
   int point1 = gmshModelOccAddPoint(0, 0, 0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int point2 = gmshModelOccAddPoint(195, 0, 0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int point3 = gmshModelOccAddPoint(40,220, 0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int point4 = gmshModelOccAddPoint(30, 245, 0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int point5 = gmshModelOccAddPoint(15, 285, 0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int point6 = gmshModelOccAddPoint(0, 285, 0, 0, -1, &ierr);
    ErrorGmsh(ierr);

    // Créer les lignes à partir des points
    int idLine1 = gmshModelOccAddLine(point1, point2, -1, &ierr);
    ErrorGmsh(ierr);
    int idLine2 = gmshModelOccAddLine(point2, point3, -1, &ierr);
    ErrorGmsh(ierr);
    int idLine3 = gmshModelOccAddLine(point3, point4, -1, &ierr);
    ErrorGmsh(ierr);
    int idLine4 = gmshModelOccAddLine(point4, point5, -1, &ierr);
    ErrorGmsh(ierr);
    int idLine5 = gmshModelOccAddLine(point5, point6, -1, &ierr);
    ErrorGmsh(ierr);
    int idLine6 = gmshModelOccAddLine(point6, point1, -1, &ierr);
    ErrorGmsh(ierr);

    // Créer une boucle de courbes
    int curveTags[] = {idLine1, idLine2, idLine3, idLine4, idLine5, idLine6};
    int idCurveLoop = gmshModelOccAddCurveLoop(curveTags, 6, -1, &ierr);
    ErrorGmsh(ierr);

    // Ajouter la boucle de courbes au modèle
    int surfaceTags[] = {idCurveLoop};
    int idSurface = gmshModelOccAddPlaneSurface(surfaceTags, 1, -1, &ierr);
    ErrorGmsh(ierr);
 
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

  return;
}



