class Resolving_interaction{
    
  //ImpliciteParticles particles;
 
  float rho = 1; // densite (kg/L)
  float sigma = 73; // tension superficielle (mN/m)
  float k_1 = 11.5; //donnees
  float k_2 = 0.45;
  
  float [] phi = {-1.23456, -1.23456}; // fraction de volume qui se chevauche
  
  float delta =             -1.23456;; // Size ratio
  float f_delta =           -1.23456;
  float norm_vel_relative = -1.23456;
  float We =                -1.23456;
  float We_reflex =         -1.23456;
  float X =                 -1.23456;
  float ksi =               -1.23456;
  float eta_1 =             -1.23456;
  float eta_2 =             -1.23456;
  float V_lig =             -1.23456; // volume du ligament
  float Z =                 -1.23456; // pour calcul de la nvelle vitesse de la gouttelette
  float beta =              -1.23456; // pour calcul rayon rupture du ligament
  float r_0 =               -1.23456;
  float We_0 =              -1.23456;
  //float radius_bu =         -1.23456;
  float radius_bu =         2;
  float radius_sat =        -1.23456;
  float n_sat =             -1.23456;
   
    
  void init() {
    //particles = new ImpliciteParticles();
    
    delta = particles.radius.get(1) / particles.radius.get(0);
    f_delta = pow(delta,-3)-2.4*pow(delta,-2)+2.7*pow(delta,-1);
    
    // Impact (plus courte distance entre les deux centres des gouttelettes orthogonaux à la vitesse relative)
    PVector u_dist = vel_relative(particles.velocity).mult(1/norm_vel_relative);
    PVector u_dist2 = new PVector(u_dist.x * u_dist.x, u_dist.y * u_dist.y, u_dist.z * u_dist.z);
    PVector x = u_dist2.mult((particles.points.get(1).x - particles.points.get(0).x) - (particles.points.get(1).y - particles.points.get(0).y));
    float norm_x = norm(x) ;
    float X = norm_x / (particles.radius.get(0) + particles.radius.get(1)); 
    
    // fraction de volume qui se chevauche
    float hauteur_collision = (particles.points.get(1).y+particles.radius.get(1)) - (particles.points.get(0).y-particles.radius.get(0));
    for (int i=0; i<2; i++){
      phi[i] = volume_calotte(i, hauteur_collision)/volume(i);
    }
    
    ksi = 0.5*X*(1+delta);
    
    eta_1 = 2*pow(1-ksi,2)*sqrt((1-pow(ksi,2)))-1;
    eta_2 = 2*pow(delta-ksi,2)*sqrt(pow(delta,2)-pow(ksi,2))-pow(delta,3);
    

    // calcul de la norme de la velocité relative
    norm_vel_relative = norm(vel_relative(particles.velocity));
    
    // Weber Number (rapport entre les forces d'inertie et les forces de tension superficielle)
    We = (2 * rho * particles.radius.get(1) * pow(norm_vel_relative,2)) / sigma;
    
    //// Reflexive separation
    //EQ2
    We_reflex = (3*(7*pow((1+pow(delta,3)),2/3)-4*(1+pow(delta,2)))*delta*pow((1+pow(delta,3)),2)) / (pow(delta,6)*eta_1+eta_2);

    //// Instabilite des ligaments
    beta = 3/(4*sqrt(2))*(k_1*k_2);
    //TODO : déterminer r_0
    We_0 = 2* r_0 * rho/sigma * pow(norm_vel_relative,2); //avec r_0 rayon initial egal a longueur initiale du cylindre 

    // determine le type de collision
    //calcul_We();
    
  }
  
  // determine le type de collision et lance la fonction de mise a jour des parametres associé
  void calcul_We(){
    //// Stretching separation
    //EQ1
    float We_stretch = (4*pow(1+pow(delta,3),2) * sqrt(3*(1+delta)*(1-X)*(pow(delta,3)*phi[1]+phi[0]))) / (pow(delta,2)*((1+pow(delta,3))-(1-pow(X,2))*(phi[1]+pow(delta,3)*phi[0])));
    
    System.out.println("We "+ We);
    System.out.println("We_stretch "+ We_stretch);
                  
    // Reflexive separation
    if (We > We_reflex){
      System.out.println("Reflexive separation");
      reflexive_separation();
    }
    
    // Stretching separation
    else if (We > We_stretch){
      System.out.println("Stretching separation");
      stretching_separation();} //TODO : PROBLEME LIGNE D'AVANT
    
    // Coalescence
    else {
      System.out.println("Coalescence");
      coalescence();}
  }
  
  
  //TODO : EXECUTER METHODE EN FONCTION DE LA SORTIE DE CALCUL_WE()
  
  
  // Reflexive separation
  void reflexive_separation(){
    // mise à jour de la vitesse des gouttelettes
    Z = sqrt(1-We_reflex/We); 
    for (int i=0; i<2; i++){ 
      particles.velocity.get(i).x = calcul_velocity(i, particles.velocity.get(i).x, particles.velocity.get(1-i).x);
      particles.velocity.get(i).y = calcul_velocity(i, particles.velocity.get(i).y, particles.velocity.get(1-i).y);
      particles.velocity.get(i).z = calcul_velocity(i, particles.velocity.get(i).z, particles.velocity.get(1-i).z);
    }  
    V_lig = V_ligament();
    
    if (V_lig > 0){
      // TODO :
      // EQ8 r_sat 
      //radius_bu = equation_degre3(beta*sqrt(We_0), 1, 0, -1, r_0); choisir le bon parmis les 3
      // rayon des satellites
      radius_sat = 1.89 * radius_bu;
      particles.radius.add(radius_sat);
      // nombre de satellites par conservation de masse
      float volume_sat = volume(2);
      n_sat = V_lig/volume_sat;
      n_sat = n_sat - 2;
    }
    // presence de satellites
    if (n_sat > 0){
      // rayons des satellites
      for (int i=1; i<n_sat-1; i++){
        particles.radius.add(radius_sat);
      }
      // vitesse des satellites
      for (int i=1; i<n_sat; i++){
        particles.velocity.add(velocity_sat(i));
      }
      // rayons des gouttelettes = rayon du satellite (3 gouttelettes de meme rayon)
      for (int i=0; i<2; i++){
        particles.radius.set(i,radius_sat);
      }
    }
    else {
      particles.radius.remove(2);
    }
  }
  
  // Stretching separation
  void stretching_separation(){
    // mise à jour de la vitesse des gouttelettes
    Z = (X-sqrt((2.4*f_delta)/We)) / (1-sqrt((2.4*f_delta)/We));
    for (int i=0; i<2; i++){ 
      particles.velocity.get(i).x = calcul_velocity(i, particles.velocity.get(i).x, particles.velocity.get(1-i).x);
      particles.velocity.get(i).y = calcul_velocity(i, particles.velocity.get(i).y, particles.velocity.get(1-i).y);
      particles.velocity.get(i).z = calcul_velocity(i, particles.velocity.get(i).z, particles.velocity.get(1-i).z);
    }
    V_lig = V_ligament(); 
    
    if (V_lig > 0){
      // TODO :
      //r_0
      // EQ8 r_sat 
      // radius_bu = equation_degre3(beta*sqrt(We_0), 1, 0, -1, r_0)
      // rayon des satellites
      radius_sat = 1.89 * radius_bu;
      particles.radius.add(radius_sat);
      // nombre de satellites par conservation de masse
      float volume_sat = volume(2);
      n_sat = V_lig/volume_sat;
    }
    // presence de satellites
    if ( n_sat > 0){
      // rayon des satellites
      for (int i=1; i<n_sat-1; i++){
        particles.radius.add(radius_sat);
      }
      //vitesse des satellites
      for (int i=1; i<n_sat; i++){
        particles.velocity.add(velocity_sat(i));
      }
      
      //calcul des rayons des gouttelettes par conservation de masse
      for (int i=0; i<2; i++){
        float new_volume = volume(i) - v_interact(i);
        particles.radius.set(i,new_radius(new_volume));
      }
    } 
    else {
      particles.radius.remove(2);
    }
  }
  
  // Coalescence
  void coalescence(){
    //calcul du nouveau volume
    float new_volume = volume(0) + volume(1);
    //mise a jour de la vitesse (conservation de la quantité de mouvement [energie cinetique])
    particles.velocity.get(0).x = (volume(0)*particles.velocity.get(0).x + volume(1)*particles.velocity.get(1).x) / new_volume;
    particles.velocity.get(0).y = (volume(0)*particles.velocity.get(0).y + volume(1)*particles.velocity.get(1).y) / new_volume;
    particles.velocity.get(0).z = (volume(0)*particles.velocity.get(0).z + volume(1)*particles.velocity.get(1).z) / new_volume;
    //mise a jour de la position (centre de masse)
    particles.points.get(0).x = (volume(0)*particles.points.get(0).x + volume(1)*particles.points.get(1).x) / new_volume;
    particles.points.get(0).y = (volume(0)*particles.points.get(0).y + volume(1)*particles.points.get(1).y) / new_volume;
    particles.points.get(0).z = (volume(0)*particles.points.get(0).z + volume(1)*particles.points.get(1).z) / new_volume;
    // mise a jour du rayon (conservation de masse)
    particles.radius.set(0, new_radius(new_volume));
    particles.radius.set(1, 0.0);
  }
  
  
  /**
   * Calcule la vitesse relative entre les deux spheres
   **/
  PVector vel_relative( ArrayList<PVector> vel){
    PVector u_j = vel.get(0).copy();
    u_j.mult(-1);
    
    PVector u_ij = vel.get(1).copy();
    u_ij.add(u_j); 
    return u_ij;
  }
  
  /**
   * Calcule la norme d'un vecteur
   **/
  float norm( PVector vector){
    return sqrt(pow(vector.x,2)+pow(vector.y,2)+pow(vector.z,2)); //TODO : Verif formule
  }
  
  // Resolution d'une equation du troisieme degre
  ArrayList<Float> equation_degre3(float a, float b, float c, float d, float r_0){
    float q = (2*pow(b,3)-9*a*b*c+27*pow(a,2)*d)/(27*pow(a,3));
    float p = (3*a*c-pow(b,2))/(3*pow(a,2));
    float delta_1 = pow(q,2)+(4*pow(p,3))/27;
    float X1 = pow((-q-sqrt(delta_1))/2, 0.33333333) + pow((-q+sqrt(delta_1))/2, 0.33333333) - b/(3*a);
    float delta_2 = pow(b+a*X1,2) - 4*a*(c+(b+a*X1)*X1);
    float X2 = (-b-a*X1-sqrt(delta_2))/(2*a);
    float X3 = (-b-a*X1+sqrt(delta_2))/(2*a);
    
    float x1 = X1 * r_0; // TODO : X1 *=r_0 ?
    float x2 = X2 * r_0;
    float x3 = X3 * r_0;
    ArrayList<Float> r_bu = new ArrayList<Float>();
    r_bu.add(x1);
    r_bu.add(x2);
    r_bu.add(x3);
    return r_bu;
  }
  
  // Actualisation de la velocite des gouttelettes
  //EQ3
  float calcul_velocity(int k, float vel_k, float vel_l){
    return (pow(particles.radius.get(k),3)*vel_k + pow(particles.radius.get(1-k),3)*vel_l - pow(particles.radius.get(1-k),3)*(vel_k-vel_l)*Z) / (pow(particles.radius.get(k),3)+pow(particles.radius.get(1-k),3));
  }
  
  // volume dune sphere
  float volume(int k){
    return 4/3*PI*pow(particles.radius.get(k),3); 
  }
  
  // rayon a partir dun volume
  float new_radius(float V){
    return pow(((V/PI)*3/4),0.33333333); //TODO : RACINE CUBIQUE..... AVEC TEST (pow(27,0.33333333)) CELA DONNE 3.0
  }
  
  // calcul de phi 
  float volume_calotte(int k, float h){
    return (PI*pow(h,2)*(3*particles.radius.get(k)-h))/3;
  }
  
  // volume en interaction au niveau de la collision
  float v_interact(int k){
    return volume(k)*phi[k];
  }
  
  // calcul de la vitesse des satellites avec n le numero du sat dont on calcul la vitesse et nsat le nbr totale de sat
  //EQ9
  PVector velocity_sat(int n){
    float velx = v_interact(0)/V_lig * particles.velocity.get(0).x + v_interact(1)/V_lig * (particles.velocity.get(1).x+(1-2*n/(n_sat+1))*vel_relative(particles.velocity).x);
    float vely = v_interact(0)/V_lig * particles.velocity.get(0).y + v_interact(1)/V_lig * (particles.velocity.get(1).y+(1-2*n/(n_sat+1))*vel_relative(particles.velocity).y);
    float velz = v_interact(0)/V_lig * particles.velocity.get(0).z + v_interact(1)/V_lig * (particles.velocity.get(1).z+(1-2*n/(n_sat+1))*vel_relative(particles.velocity).z);
    return new PVector(velx, vely, velz);
  }
  
  // volume du ligament liquide resultant
  float V_ligament(){
    float tau = (1-X)*(1+delta);
    float E_dissip = 0.30; // pour l'energie cinetique initiale des gouttelettes en collision
    //EQ7
    float E_surface = 2*sigma*sqrt(PI*volume(0)*particles.radius.get(0)*tau*(phi[0]+ pow(delta,3)* phi[1]));
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
