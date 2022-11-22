class ImpliciteParticles {
  ArrayList<PVector> points = new ArrayList<PVector>();
  ArrayList<PVector> dirs = new ArrayList<PVector>();
  //velocity
  ArrayList<PVector> velocity = new ArrayList<PVector>();
  //float [] vel = {2, 8}; 

  float tailleZ = width-radius[0];
  float posX1 = radius[0];
  float posX2 = width-radius[1];
  float posY = height/2;
  float diff_hauteur = 75;


  /**
   * Creation des spheres
   **/
  ImpliciteParticles() {
    // creation de la gouttelle de gauche
    points.add(new PVector(posX1, posY, tailleZ/2));
    dirs.add(new PVector(1, 0, 0));
    velocity.add( new PVector(2, 0, 0));
    
    // creation de la gouttelle de droite
    points.add(new PVector(posX2, posY + diff_hauteur, tailleZ/2));
    dirs.add(new PVector(-1, 0, 0));
    velocity.add( new PVector(8, 0, 0));
  }
  

  /**
   * Creation du mouvement
   **/
  void nextStep() {
    for (int i = 0; i < N; i++) {
      PVector p = points.get(i);
      PVector d = dirs.get(i);
      PVector dd = d.copy();
      if (i==0) dd.mult(velocity.get(i).x);
      else dd.mult(velocity.get(i).x);
      p.add(dd);
    }
  }

  
  /**
   *
   **/
  float eval(float i, float j, float k) {
    PVector n = new PVector(i, j, k);
    float v = 0;
    for (int r=0; r<N; r++){
      v += f(points.get(r).dist(n)/radius[r]);
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
