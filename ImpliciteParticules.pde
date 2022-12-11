class ImpliciteParticles {
  ArrayList<PVector> points = new ArrayList<PVector>();
  //ArrayList<PVector> dirs = new ArrayList<PVector>();
  ArrayList<PVector> velocity = new ArrayList<PVector>();
  //radius
  ArrayList<Float> radius = new ArrayList<Float>();
  //float [] radius = {100, 100}; //TODO : r[0]>=r[1]

  int minimumY;
  int maximumY;
  int minimumZ;
  int maximumZ;
  float tailleZ;
  float posX1;
  float posX2;
  float posY = 0;
  float diff_hauteur = 10;
  float zoom = 20;
  int espace_begin = 300;
  

  /**
   * Creation des spheres
   **/
  ImpliciteParticles() {
    radius.add(6.0);
    radius.add(6.0);
    // Initialisation de variables
    tailleZ = width-radius.get(0);
    posX1 = 0- espace_begin;
    posX2 = 0+ espace_begin;
    
    // creation de la gouttelle de gauche
    points.add(new PVector(posX1, posY+diff_hauteur*0.5, tailleZ/2));
    //dirs.add(new PVector(1, 0, 0));
    velocity.add( new PVector(4, 0, 0));
   
    // creation de la gouttelle de droite
    points.add(new PVector(posX2, posY - diff_hauteur*0.5, tailleZ/2));
    //dirs.add(new PVector(-1, 0, 0));
    velocity.add( new PVector(-4, 0, 0));
    
    /** Stretch
      Radius = 6
      Diff_h = 10
      Vitesse = 4 / -4
     */

    // minimumY = (int)min(points.get(1).y-(radius.get(1)+1)*zoom, points.get(0).y-(radius.get(0)+1)*zoom);
    // maximumY = (int)max(points.get(1).y+(radius.get(1)+1)*zoom, points.get(0).y+(radius.get(0)+1)*zoom);  
  }
  
  
  /**
   * Creer une nouvelle particule dans la liste "point" pour afficher une boule supplémentaire
   * dans le rendu. Les listes velocity et radius doivent contenir les nouvelles valeur pour la
   * nouvelle boule qui sera ajouté.
   **/
  void printMissingParticles () {
    System.out.println("N = "+N);
    System.out.println("Nb points = "+points.size());
    System.out.println("Nb velocity = "+velocity.size());
    System.out.println("Nb radius = "+radius.size());

    float decalage_y = particles.diff_hauteur/(particles.velocity.size()-2);
    float y = decalage_y;

    for ( int i = 2; i < velocity.size(); i++ ) {
      // Si velocity et radius ne possède pas les données pour la nouvelle boule à ajouter
      if ( velocity.size() <= points.size() && radius.size() <= points.size() ) {
        System.err.println("printMissingParticles() : Les listes Velocity ou Radius ne possèdent pas les informations pour créer une nouvelles goutte");
      }else{
        points.add(new PVector(points.get(0).x, particles.points.get(1).y + y, points.get(0).z)); 
        y += decalage_y;
        N ++;
      }
    }
    
    
    System.out.println("N = "+N);
    System.out.println("Nb points = "+points.size());
    System.out.println("Nb velocity = "+velocity.size());
    System.out.println("Nb radius = "+radius.size());
  }
  

  /**
   * Creation du mouvement
   **/
  void nextStep() {
    for (int i = 0; i < N; i++) {
      if ( i < N ) {
        PVector p = points.get(i);
        p.add(velocity.get(i));
      }
    }

    // to have the origin at the center of the screen
    translate(width/2, height/2, 0);

    for ( int i = 0; i < N; i++ ) {
      translate(particles.points.get(i).x, particles.points.get(i).y*zoom, 0);
      sphere(particles.radius.get(i)*zoom);
      translate(-particles.points.get(i).x, -particles.points.get(i).y*zoom, 0);
    }
    
    translate(-width/2, -height/2, 0);

  }

  
  /**
   *
   **/
  float eval(float i, float j, float k) {
    PVector n = new PVector(i, j, k);
    float v = 0;
    for (int r=0; r<N; r++){
      v += f(points.get(r).dist(n)/(radius.get(r)*zoom));
    }
    return v;
  }


  /**
   * Evaluation es sommets (0 ou 1)
   **/
  int evalInt(float i, float j, float k) {
    return (eval(i, j, k) >= 0.1) ? 1 : 0;
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
