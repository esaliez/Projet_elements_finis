#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH


//définir la géométrie
double geoSize(double x, double y) {
  
  femGeo *theGeometry = geoGetGeometry();
  double h = theGeometry->h;


  int h_min = 1; //taille minimale des triangles 3
  int h_max = 40; //taille maximale des triangles 6

  int b_X_1 = 2; //taille de la  zone grande précision coté gauche 20
  int b_X_2 = 5; //taille de la zone de transistion entre grande et petite précision 5
  int b_X_3 = 195; //taille maximale du barrage en x 
  
  int coeff_y = 1; //coefficient de mutiplication de la taille minimal des triangles en fonction de y
  int b_Y_1 = 285; //début de la zone de grande précision en y
  int b_Y_2 = 220; //zone intermédiaire en y (pas utilisée pour l'instant)
  int b_Y_3 = 0; //fin de zone en y (pas utilisée pour l'instant)

  //zone précise à gauche
  if(x <= b_X_1) {
    h = h_min;
  }
  
  //transisition entre la précision à gauche et la zone qui évolue selon y
  else if(x > b_X_1 && x <= b_X_2+b_X_1){
    h = h_min + (x - b_X_1)/(b_X_1);
  }
  //zone qui évolue selon y, grande précision en haut et faible en bas
  else if(x > b_X_1+b_X_2){
     h = h_min + coeff_y*(b_Y_1-y)/(b_Y_1/h_max);
  }
  
  h = 2;
  
  return h;
}


void geoMeshGenerate(void) {
  femGeo *theGeometry = geoGetGeometry();
  
 
  //theGeometry->elementType = FEM_QUAD;
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
    int point4 = gmshModelOccAddPoint(15, 285, 0, 0, -1, &ierr);
    ErrorGmsh(ierr);
    int point5 = gmshModelOccAddPoint(0, 285, 0, 0, -1, &ierr);
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
    int idLine5 = gmshModelOccAddLine(point5, point1, -1, &ierr);
    ErrorGmsh(ierr);

    // Créer une boucle de courbes
    int curveTags[] = {idLine1, idLine2, idLine3, idLine4, idLine5};
    int idCurveLoop = gmshModelOccAddCurveLoop(curveTags, 5, -1, &ierr);
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



