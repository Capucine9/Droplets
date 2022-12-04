class ImpliciteParticles {
  ArrayList<PVector> points = new ArrayList<PVector>();
  //ArrayList<PVector> dirs = new ArrayList<PVector>();
  ArrayList<PVector> velocity = new ArrayList<PVector>();
  //radius
  ArrayList<Float> radius = new ArrayList<Float>();
  //float [] radius = {100, 100}; //TODO : r[0]>=r[1]


  float tailleZ;
  float posX1;
  float posX2;
  float posY = height/2;
  float diff_hauteur = 10;


  /**
   * Creation des spheres
   **/
  ImpliciteParticles() {
    radius.add(6.0);
    radius.add(6.0);
    // Initialisation de variables
    tailleZ = width-radius.get(0);
    posX1 = radius.get(0);
    posX2 = width-radius.get(1);
    
    // creation de la gouttelle de gauche
    points.add(new PVector(posX1, posY, tailleZ/2));
    //dirs.add(new PVector(1, 0, 0));
    velocity.add( new PVector(4, 0, 0));
    
    // creation de la gouttelle de droite
    points.add(new PVector(posX2, posY - diff_hauteur, tailleZ/2));
    //dirs.add(new PVector(-1, 0, 0));
    velocity.add( new PVector(-4, 0, 0));
  }
  
  
  /**
   * Creer une nouvelle particule dans la liste "point" pour afficher une boule supplémentaire
   * dans le rendu. Les listes velocity et radius doivent contenir les nouvelles valeur pour la
   * nouvelle boule qui sera ajouté.
   **/
  void printMissingParticles () {
    // Si velocity et radius ne possède pas les données pour la nouvelle boule à ajouter
    if ( velocity.size() <= points.size() && radius.size() <= points.size() ) {
      System.err.println("printMissingParticles() : Les listes Velocity ou Radius ne possèdent pas les informations pour créer une nouvelles boule");
    }else{
      points.add(new PVector(points.get(0).x, points.get(0).y-10, points.get(0).z)); 
      N ++;
    }
  }
  

  /**
   * Creation du mouvement
   **/
  void nextStep() {
    for (int i = 0; i < N; i++) {
      PVector p = points.get(i);
      p.add(velocity.get(i));
    }
  }

  
  /**
   *
   **/
  float eval(float i, float j, float k) {
    PVector n = new PVector(i, j, k);
    float v = 0;
    for (int r=0; r<N; r++){
      v += f(points.get(r).dist(n)/(radius.get(r)*10));
    }
    return v;
  }


  /**
   * Evaluation es sommets (0 ou 1)
   **/
  int evalInt(float i, float j, float k) {
    return (eval(i, j, k) >= 0.00001) ? 1 : 0;
  }
  

  PVector interp(float x1, float y1, float z1, float x2, float y2, float z2) {
    return interpRec(x1, y1, z1, x2, y2, z2, REC_INTERP);
  }


  
  /**
   * Algo dichotomie
   **/
  PVector interpRec(float x1, float y1, float z1, float x2, float y2, float z2, int level) {
    PVector p1 = new PVector(x1, y1, z1);
    PVector p2 = new PVector(x2, y2, z2);
    PVector m = PVector.lerp(p1, p2, 0.5); //lerp permet de retrouver le milieu de p1p2 (milieu d'un bord d'un pixel) 
    float v1 = evalInt(x1, y1, z1);
    float v2 = evalInt(x2, y2, z2);
    if ((level == 0) || (v1 == v2)) 
      return m;
    if (evalInt(m.x, m.y, m.z) == v1) return interpRec(m.x, m.y, m.z, x2, y2, z2, level-1);
    return interpRec(x1, y1, z1, m.x, m.y, m.z, level-1);
    // on arrete la recursion si les points a droite et a gauche (coins d'un pixel) ont la meme valeur (meme couleur) et on recupere donc le milieu du cote du pixel
    //avec les 5 appels recursifs on arrive assez bien a ce rapprocher du point ou la valeur change
  } 
}
