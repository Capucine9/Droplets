class Resolving_interaction{
  
  ImpliciteParticles particles;
  
<<<<<<< HEAD
  // calcul de la norme de la velocité relative
  float norm_vel_relative = norm(vel_relative(particles.velocity));
  
  // Weber Number (rapport entre les forces d'inertie et les forces de tension superficielle)
  float We = (2 * rho * pow(norm_vel_relative,2)) / sigma;
  
  //// Reflexive separation
  //EQ2
  float We_reflex = (3*(7*pow((1+pow(delta,3)),2/3)-4*(1+pow(delta,2)))*delta*pow((1+pow(delta,3)),2)) / (pow(delta,6)*eta_1+eta_2);
    
  
  
  
  
  
  void test(){
=======
>>>>>>> 02e51abccbc1f8644a103d6320445dc295c66a33
    // fraction de volume qui se chevauche
    float [] phi = {0.1, 0.2};
    
    float rho = 1; // densite (kg/L)
    float sigma = 73; // tension superficielle (mN/m)
    
    
    // Size ratio
    float delta = radius[1] / radius[0];
    
    
    float f_delta =           -1.23456;
    float norm_vel_relative = -1.23456;
    float We =                -1.23456;
    PVector u_dist =          new PVector(0,0,0);
    PVector u_dist2 =         new PVector(0,0,0);
    PVector x =               new PVector(0,0,0);
    float norm_x =            -1.23456;
    float X =                 -1.23456;
    float ksi =               -1.23456;
    float eta_1 =             -1.23456;
    float eta_2 =             -1.23456;
    
  void init() {
    
  
    f_delta = pow(delta,-3)-2.4*pow(delta,-2)+2.7*pow(delta,-1);
    
    
<<<<<<< HEAD
    
    
    
    
    // Impact (plus courte distance entre les deux centres des gouttelettes orthogonaux à la vitesse relative)
    PVector u_dist = vel_relative(particles.velocity).mult(1/norm_vel_relative);
    PVector u_dist2 = new PVector(u_dist.x * u_dist.x, u_dist.y * u_dist.y, u_dist.z * u_dist.z);
    PVector x = u_dist2.mult((particles.points.get(1).x - particles.points.get(0).x) - (particles.points.get(1).y - particles.points.get(0).y));
    float norm_x = norm(x) ;
    float X = norm_x / (radius[0] + radius[1]); 
=======
    // Weber Number (rapport entre les forces d'inertie et les forces de tension superficielle)
    // calcul de la norme de la velocité relative
    norm_vel_relative = norm(vel_relative(particles.velocity));
    
    We = (2 * rho * pow(norm_vel_relative,2)) / sigma;
    
    
    u_dist = vel_relative(particles.velocity).mult(1/norm_vel_relative);
     u_dist2 = new PVector(u_dist.x * u_dist.x, u_dist.y * u_dist.y, u_dist.z * u_dist.z);
    x = u_dist2.mult((particles.points.get(1).x - particles.points.get(0).x) - (particles.points.get(1).y - particles.points.get(0).y));
    norm_x = norm(x) ;
    X = norm_x / (radius[0] + radius[1]); 
>>>>>>> 02e51abccbc1f8644a103d6320445dc295c66a33
    
    ksi = 0.5*X*(1+delta);
    
<<<<<<< HEAD
    float eta_1 = 2*pow(1-ksi,2)*sqrt((1-pow(ksi,2)))-1;
    float eta_2 = 2*pow(delta-ksi,2)*sqrt(pow(delta,2)-pow(ksi,2))-pow(delta,3);
    
    
    
    
    
  
  
  
    // pour calcul de la nvelle vitesse de la gouttelette
    float  z; 
    // si stretching separation : 
    z = (X-sqrt((2.4*f_delta)/We)) / (1-sqrt((2.4*f_delta)/We));
    // si reflexive separation : 
    z = sqrt(1-We_reflex/We);  



    
    
    
    //EQ4 ??????????
    // volume du ligament liquide resultant
    // Reflexive separation
    // V_lig = volume(0) + volume(1);
    
    
    // Rayon de rupture du ligament
    float radius_bu = 1;
    // Rayon des satellites
    float radius_sat = 1.89* radius_bu;

  }
  
  
  int calcul_We(){
    // Reflexive separation
    if (We > We_reflex){return 3;}
    
    //// Stretching separation
    //EQ1
    float We_stretch = (4*pow(1+pow(delta,3),2) * sqrt(3*(1+delta)*(1-X)*(pow(delta,3)*phi[1]+phi[0]))) / (pow(delta,2)*((1+pow(delta,3))-(1-pow(X,2))*(phi[1]+pow(delta,3)*phi[0])));
    else if (We > We_stretch){return 2;} //TODO : PROBLEME LIGNE D'AVANT
    
    // Coalescence
    else {return 1;}
  }
  
  
  //TODO : EXECUTER METHODE EN FONCTION DU INT RESORTI PAR CALCUL_WE()
  
  
  // Reflexive separation
  void reflexive_separation(){
    // mise à jour de la vitesse des gouttelettes
    for (int i=0; i<1; i++){ //TODO : verifier si i<1 ou<2 (2 gouttelettes)
      velocity.get(i).x = calcul_velocity(i, velocity.get(i).x, velocity.get(1-i).x);
      velocity.get(i).y = calcul_velocity(i, velocity.get(i).y, velocity.get(1-i).y);
      velocity.get(i).z = calcul_velocity(i, velocity.get(i).z, velocity.get(1-i).z);
    }  
    // TODO : 
    // EQ8 r_sat  
    // n_sat
    
    // presence de satellites
    if (n_sat > 2){
      //vitesse des satellites
      for (int i=1; i<n_sat; i++){
        velocity_sat(i);
      }
      //rayons des gouttelettes = rayon du satellite (3 gouttelettes de meme rayon)
      for (int i=0; i<1; i++){
        radius[i]=r_sat;
      }
    }

  }
  
  // Stretching separation
  void stretching_separation(){
    // mise à jour de la vitesse des gouttelettes
    for (int i=0; i<1; i++){ //TODO : verifier si i<1 ou<2 (2 gouttelettes)
      velocity.get(i).x = calcul_velocity(i, velocity.get(i).x, velocity.get(1-i).x);
      velocity.get(i).y = calcul_velocity(i, velocity.get(i).y, velocity.get(1-i).y);
      velocity.get(i).z = calcul_velocity(i, velocity.get(i).z, velocity.get(1-i).z);
    }
    V_lig = V_ligament();
    // TODO : 
    // EQ8 r_sat 
    // n_sat
    
    // presence de satellites
    if ( n_sat > 0){
      //vitesse des satellites
      for (int i=1; i<n_sat; i++){
        velocity_sat(i);
      }
      
      //calcul des rayons des gouttelettes par conservation de masse
      for (int i=0; i<1; i++){
        float new_volume = volume(i) - v_interact(i);
        radius[i] = new_radius(new_volume);
      }
    } 
  }
  
  // Coalescence
  void coalescence(){
    //calcul du nouveau volume
    float new_volume = volume(0) + volume(1);
    //mise a jour de la vitesse (conservation de la quantité de mouvement [energie cinetique])
    velocity.get(0).x = (volume(0)*velocity.get(0).x + volume(1)*velocity.get(1).x) / new_volume;
    velocity.get(0).y = (volume(0)*velocity.get(0).y + volume(1)*velocity.get(1).y) / new_volume;
    velocity.get(0).z = (volume(0)*velocity.get(0).z + volume(1)*velocity.get(1).z) / new_volume;
    //mise a jour de la position (centre de masse)
    points.get(0).x = (volume(0)*points.get(0).x + volume(1)*points.get(1).x) / new_volume;
    points.get(0).y = (volume(0)*points.get(0).y + volume(1)*points.get(1).y) / new_volume;
    points.get(0).z = (volume(0)*points.get(0).z + volume(1)*points.get(1).z) / new_volume;
    // mise a jour du rayon (conservation de masse)
    radius[0] = new_radius(new_volume);
    radius[1] = 0;
  }
  
  
  
  
  
  
  
  
  // velocite relative
  PVector vel_relative(ArrayList<PVector> vel){
=======
    eta_1 = 2*pow(1-ksi,2)*pow((1-pow(ksi,2)),1/2)-1;
    eta_2 = 2*pow(delta-ksi,2)*pow((pow(delta,2)-pow(ksi,2)),(1/2))-pow(delta,3);
  }
  
  
  /**
   * Calcule la vitesse relative entre les deux sphere
   **/
  PVector vel_relative( ArrayList<PVector> vel){
>>>>>>> 02e51abccbc1f8644a103d6320445dc295c66a33
    PVector u_j = vel.get(0).copy();
    u_j.mult(-1);
    
    PVector u_ij = vel.get(1).copy();
    u_ij.add(u_j); 
    return u_ij;
  }
  
<<<<<<< HEAD
  // calcul de norme d'un vecteur
=======
  
  /**
   * Calcule la norme d'un vecteur
   **/
>>>>>>> 02e51abccbc1f8644a103d6320445dc295c66a33
  float norm( PVector vector){
    return pow(pow(vector.x,2)+pow(vector.y,2)+pow(vector.z,2),1/2);
  }
    
  // Actualisation de la velocite des gouttelettes
  //EQ3
  // k la gouttelette dont on calcul la nvelle vitesse
  void calcul_velocity(int k, float vel_k, float vel_l){
    float u = (pow(radius[k],3)*vel_k + pow(radius[1-k],3)*vel_l - pow(radius[1-k],3)*(vel_k-vel_l)*z) / (pow(radius[k],3)+pow(radius[1-k],3));
  }
  
  // volume dune sphere
  float volume(int k){
    return 4/3*PI*pow(radius[k],3); 
  }
  
  // rayon a partir dun volume
  float new_radius(float V){
    return pow(((V/PI)*3/4),0.33333333); //TODO : RACINE CUBIQUE..... AVEC TEST (pow(27,0.33333333)) CELA DONNE 3.0
  }
  
  // volume en interaction au niveau de la collision
  float v_interact(int k){
    return volume(k)*phi[k];
  }
  
  // calcul de la vitesse des satellites avec n le numero du sat dont on calcul la vitesse et nsat le nbr totale de sat
  //EQ9
  void velocity_sat(int n){
    //float u = V_lig_i/V_lig*u~_i+V_lig_j/V_lig*(u~j+(1-2n*(n_sat+1))*u~_ij)
  }
  
  // volume du ligament liquide resultant
  float V_ligament(){
    float tau = (1-X)*(1+delta);
    float E_dissip = 0.30; // pour l'energie cinetique initiale des gouttelettes en collision
    //EQ7
    float E_surface = 2*sigma*sqrt(PI*volume(0)*radius[0]*tau*(phi[0]+ pow(delta,3)* phi[1]));
    //// Volume des ligaments liquide
    //EQ6
    float E_stretch = 1/2* rho* pow(norm_vel_relative,2)* volume(0)* (pow(delta,3)/(pow((1+pow(delta,3)),2)))* ((1+pow(delta,3))- (1-pow(X,2))* ((phi[1]+ pow(delta,3)* phi[0])));
    //EQ5
    // coeff de separation de volume
    float C_sep = (E_stretch - E_surface - E_dissip) / (E_stretch + E_surface + E_dissip);
    
    // Stretching separation
    //EQ4
    // volume du ligament liquide resultant
    return C_sep* (v_interact(0)+ v_interact(1));
  }
    
}
