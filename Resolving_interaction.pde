class Resolving_interaction{
  ImpliciteParticles particles;
  
  void test(){
    // fraction de volume qui se chevauche
    float [] phi = {0.1, 0.2};
    
    float rho = 1; // densite (kg/L)
    float sigma = 73; // tension superficielle (mN/m)
    
    // Size ratio
    float delta = radius[1] / radius[0];
  
    float f_delta = pow(delta,-3)-2.4*pow(delta,-2)+2.7*pow(delta,-1);
    
    
    // Weber Number (rapport entre les forces d'inertie et les forces de tension superficielle)
    // calcul de la norme de la velocité relative
    float norm_vel_relative = norm(vel_relative(particles.velocity));
    
    float We = (2 * rho * pow(norm_vel_relative,2)) / sigma;
    
    
    PVector u_dist = vel_relative(particles.velocity).mult(1/norm_vel_relative);
    PVector u_dist2 = new PVector(u_dist.x * u_dist.x, u_dist.y * u_dist.y, u_dist.z * u_dist.z);
    PVector x = u_dist2.mult((particles.points.get(1).x - particles.points.get(0).x) - (particles.points.get(1).y - particles.points.get(0).y));
    float norm_x = norm(x) ;
    float X = norm_x / (radius[0] + radius[1]); 
    
    float ksi = 0.5*X*(1+delta);
    
    float eta_1 = 2*pow(1-ksi,2)*pow((1-pow(ksi,2)),1/2)-1;
    float eta_2 = 2*pow(delta-ksi,2)*pow((pow(delta,2)-pow(ksi,2)),(1/2))-pow(delta,3);
  }
  
  
  // velocite relative
  PVector vel_relative( ArrayList<PVector> vel){
    PVector u_j = vel.get(0).copy();
    u_j.mult(-1);
    
    PVector u_ij = vel.get(1).copy();
    u_ij.add(u_j); 
    return u_ij;
  }
  
  float norm( PVector vector){
    return pow(pow(vector.x,2)+pow(vector.y,2)+pow(vector.z,2),1/2);
  }
    
}
