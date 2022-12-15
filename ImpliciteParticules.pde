class ImpliciteParticles {
  ArrayList<PVector> points = new ArrayList<PVector>();
  //ArrayList<PVector> dirs = new ArrayList<PVector>();
  ArrayList<PVector> velocity = new ArrayList<PVector>();
  //radius
  ArrayList<Float> radius = new ArrayList<Float>();
  //float [] radius = {100, 100}; //TODO : r[0]>=r[1]

  boolean isRuptured = false;
  boolean isCollision = false;

  float minimumY = height+1;
  float maximumY = -1;
  float tailleZ;
  float posX1;
  float posX2;
  float posY = 0;
  float zoom = 500;
  float zoom_tetra = 130;
  float diff_hauteur = 0;
  int espace_begin = 500;
  

  /**
   * Creation des spheres
   **/
  ImpliciteParticles(float dist) {
    //if (!Sphere)
    //  dist*=zoom_tetra;
    this.diff_hauteur = dist;
    System.out.println("diff_hauteur "+diff_hauteur);
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
     //diff_hauteur/= zoom_tetra;
     //points.get(0).y = posY+diff_hauteur*0.5;

    for (int i = 0; i < N-1; i++){
      if (minimumY > points.get(i).y-(radius.get(i)+1)*29)
          minimumY = points.get(i).y-(radius.get(i)+1)*29;
          
      if (maximumY < points.get(i).y+(radius.get(i)+1)*29)
          maximumY = points.get(i).y+(radius.get(i)+1)*29;
    }
  }
  
  
  /**
   * Creer une nouvelle particule dans la liste "point" pour afficher une boule supplémentaire
   * dans le rendu. Les listes velocity et radius doivent contenir les nouvelles valeur pour la
   * nouvelle boule qui sera ajouté.
   **/
  void printMissingParticles () {
    float decalage_y = ((particles.points.get(1).y+particles.radius.get(1)) - (particles.points.get(0).y-particles.radius.get(0)))/(particles.velocity.size()-2);
    float y = 0/2;
    System.out.println("decalage_y = "+decalage_y);
    System.out.println("diff = "+(points.get(0).y+radius.get(0)-points.get(1).y-radius.get(1)));
    System.out.println("P0 = "+points.get(0).y);
    System.out.println("P1 = "+points.get(1).y);
    System.out.println("R0 = "+radius.get(0));
    System.out.println("R1 = "+radius.get(1));
    for ( int i = 2; i < velocity.size(); i++ ) {
      // Si velocity et radius ne possède pas les données pour la nouvelle boule à ajouter
      if ( velocity.size() <= points.size() && radius.size() <= points.size() ) {
        System.err.println("printMissingParticles() : Les listes Velocity ou Radius ne possèdent pas les informations pour créer une nouvelles goutte");
      }else{
        if (particles.points.get(0).y == particles.points.get(1).y)
          points.add(new PVector(points.get(0).x, particles.points.get(0).y, points.get(0).z)); 
        else
          points.add(new PVector(points.get(0).x, particles.points.get(0).y-radius.get(0) + y, points.get(0).z)); 
        y += decalage_y;
        N ++;
      }
    }
    
    




  }
  

  /**
   * Creation du mouvement
   **/
  void nextStep() {
    for (int i = 0; i < N; i++) {
      if ( i < N ) {
        if ( (i < 2 || isRuptured ) && run ) {
          PVector p = points.get(i);
          p.add(velocity.get(i));
        }
      }
    }

    // to have the origin at the center of the screen
    //translate(width/2, height/2, 0);

    // for ( int i = 0; i < N; i++ ) {
    //   if ( isCollision && !isRuptured ) {
    //     float distance = resolv.norm(new PVector(points.get(0).x-points.get(1).x,points.get(0).y-points.get(1).y,0)) - radius.get(0) - radius.get(1);
    //     float decalage = (distance / (N-2);

    //     // cas impair, satelite du milieu
    //     if ( (N-2)%2 == 1 && i == round((N-2)/2) )  {

    //     }else if ( (N-2)%2 == 0 && i == round((N-2)/2) )  {



    //   }else{
    //     translate(particles.points.get(i).x, particles.points.get(i).y*zoom, 0);
    //     sphere(particles.radius.get(i)*zoom);
    //     translate(-particles.points.get(i).x, -particles.points.get(i).y*zoom, 0);
    //   }
    // }


    if ( isCollision && !isRuptured ) {
      for ( int i = 0; i < 2; i++ ) {
        if (print_sphere){
          translate(particles.points.get(i).x, particles.points.get(i).y*zoom, 0);
          sphere(particles.radius.get(i)*zoom);
          translate(-particles.points.get(i).x, -particles.points.get(i).y*zoom, 0);
        }
      }
      float demi_distance = (points.get(0).x-points.get(1).x-(radius.get(0)+radius.get(1))*0.5)*0.5;
      



      
      // cas impair
      if ( (N-2)%2 == 1 ) {
        float decalage = (demi_distance*2) / ((N-2)+1);
        int indice_milieu = (int) ((N+1)*0.5);
        
        // centre
        particles.points.get(indice_milieu).x = (particles.points.get(0).x - demi_distance);
        if (print_sphere) {
          translate(particles.points.get(indice_milieu).x, particles.points.get(indice_milieu).y*zoom, 0);
          sphere(particles.radius.get(indice_milieu)*zoom);
          translate(-particles.points.get(indice_milieu).x, -particles.points.get(indice_milieu).y*zoom, 0);
        }
        

        for ( int i = 1; i <= (N-3)*0.5; i++ ) {

          int indice_gauche = indice_milieu - i;
          int indice_droite = indice_milieu + i;
          int x_gauche = 0;
          int x_droite = 0;
          x_gauche = (int) (particles.points.get(1).x + demi_distance - i*decalage);
          x_droite = (int) (particles.points.get(0).x - demi_distance + i*decalage );
           
          
          particles.points.get(indice_gauche).x = x_gauche;
          particles.points.get(indice_droite).x = x_droite;
          if (print_sphere) {
            translate(particles.points.get(indice_gauche).x, particles.points.get(indice_gauche).y*zoom, 0);
            sphere(particles.radius.get(indice_gauche)*zoom);
            translate(-particles.points.get(indice_gauche).x, -particles.points.get(indice_gauche).y*zoom, 0);
            translate(particles.points.get(indice_droite).x, particles.points.get(indice_droite).y*zoom, 0);
            sphere(particles.radius.get(indice_droite)*zoom);
            translate(-particles.points.get(indice_droite).x, -particles.points.get(indice_droite).y*zoom, 0);
          }
        }

      // cas pair
      }else{
        float decalage = (demi_distance*2) / ((N-2));
        int indice_milieu = (int) ((N-2)*0.5) + 1;
        

        for ( int i = 0; i < ((N-2)*0.5); i++ ) {

          int indice_gauche = indice_milieu - i;
          int indice_droite = indice_milieu + (i+1);
          int x_gauche = 0;
          int x_droite = 0;
          if ( i == 0 ) {
            x_gauche = (int) (particles.points.get(1).x + demi_distance - 0.5*decalage);
            x_droite = (int) (particles.points.get(0).x - demi_distance );
          }else{
            x_gauche = (int) (particles.points.get(1).x + demi_distance - i*decalage - 0.5*decalage);
            x_droite = (int) (particles.points.get(0).x - demi_distance + i*decalage );
          }
          

          particles.points.get(indice_gauche).x = x_gauche;
          particles.points.get(indice_droite).x = x_droite;
          if (print_sphere) {
            translate(particles.points.get(indice_gauche).x, particles.points.get(indice_gauche).y*zoom, 0);
            sphere(particles.radius.get(indice_gauche)*zoom);
            translate(-particles.points.get(indice_gauche).x, -particles.points.get(indice_gauche).y*zoom, 0);
            translate(particles.points.get(indice_droite).x, particles.points.get(indice_droite).y*zoom, 0);
            sphere(particles.radius.get(indice_droite)*zoom);
            translate(-particles.points.get(indice_droite).x, -particles.points.get(indice_droite).y*zoom, 0);
          }

        }
      }
      float decalage = demi_distance / (N-3);


    }else{
      for ( int i = 0; i < N; i++ ) {
        if ( i < 2 ) {
          if (print_sphere) {  
            translate(particles.points.get(i).x, particles.points.get(i).y*zoom, 0);
            sphere(particles.radius.get(i)*zoom);
            translate(-particles.points.get(i).x, -particles.points.get(i).y*zoom, 0);
          }
        }else{
          if (print_sphere) {
            translate(particles.points.get(i).x, particles.points.get(i).y*zoom, 0);
            sphere(particles.radius.get(i)*zoom);
            translate(-particles.points.get(i).x, -particles.points.get(i).y*zoom, 0);
          }
        }
      }
    }
    
    //translate(-width/2, -height/2, 0);
    
  }

  
  /**
   *
   **/
  float eval(float i, float j, float k) {
    PVector n = new PVector(i, j, k);
    float v = 0;
    for (int r=0; r<2; r++){
      PVector tmp = new PVector(points.get(r).x, points.get(r).y*zoom_tetra, points.get(r).z);
      v += f(tmp.dist(n)/(radius.get(r)*zoom_tetra));
    }
    for (int r=2; r<N; r++){
      float zoom_sat = 1;
      if ( radius.get(r) < 0.3 ) zoom_sat = 3;
      PVector tmp = new PVector(points.get(r).x, points.get(r).y*zoom_tetra*zoom_sat, points.get(r).z);
      v += f(tmp.dist(n)/(radius.get(r)*zoom_tetra*zoom_sat));
    }
    return v;
  }


  /**
   * Evaluation es sommets (0 ou 1)
   **/
  int evalInt(float i, float j, float k) {
    return (eval(i, j, k) >= 0.0001) ? 1 : 0;
  }
  

  PVector interp(float x1, float y1, float z1, float x2, float y2, float z2) {
    return interpRec(x1, y1, z1, x2, y2, z2, REC_INTERP);
  }


  boolean isRupture () {
    boolean ret = false;

    if ( N > 2 && isCollision ) {
      float distance = 0;
      if ( print_sphere ) {
        distance = resolv.norm(new PVector(points.get(0).x-points.get(1).x,0,0)) - radius.get(0) - radius.get(1);
        if ( distance > radius.get(2)*2*(N-2)*zoom ) {
          isRuptured = true;
          System.err.println("Rupture");
        }
    
      }else{
        distance = resolv.norm(new PVector(points.get(0).x-points.get(1).x,0,0)) - radius.get(0) - radius.get(1);
        if ( distance > radius.get(2)*2*(N-2)*zoom_tetra ) {
          isRuptured = true;
          System.err.println("Rupture");
        }
      
      }

      
    }
    return ret;
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
