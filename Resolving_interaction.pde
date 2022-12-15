/**
 * Class in charge of resolving interaction between two droplets.
 **/
class Resolving_interaction{
 
	float rho = 1; 		// density (kg/L)
	float sigma = 73; 	// surface tension = 73 (mN/m)
	float k_1 = 11.5; 	
	float k_2 = 0.45;
  
  	float [] phi = {-1.23456, -1.23456}; // overlapping volume fraction
  
	float delta =             -1.23456; 	// Size ratio
	float f_delta =           -1.23456;
	float norm_vel_relative = -1.23456;
	float We =                -1.23456;
	float We_reflex =         -1.23456;
	double We_stretch =       -1.23456;
	float X =                 -1.23456;
	float ksi =               -1.23456;
	float eta_1 =             -1.23456;
	float eta_2 =             -1.23456;
	float V_lig =             -1.23456; 	// ligament volume
	float Z =                 -1.23456;
	float beta =              -1.23456;
	float r_0 =               -1.23456;
	float We_0 =              -1.23456;

	float radius_bu =         2;
	float radius_sat =        -1.23456;
	float n_sat =             -1.23456;
  
	float hauteur_collision;
   
    
	/**
	 * Initialize the value used to make the prediction
	 **/
  	void init() {
		delta = particles.radius.get(1) / particles.radius.get(0);
		f_delta = pow(delta,-3)-2.4*pow(delta,-2)+2.7*pow(delta,-1);
    
		// Impact (shortest distance between the center of droplets orthogonal to relative velocity)
		PVector u_dist = vel_relative(particles.velocity).mult(1/norm_vel_relative);
		PVector u_dist2 = new PVector(u_dist.x * u_dist.x, u_dist.y * u_dist.y, u_dist.z * u_dist.z);
		PVector x = u_dist2.mult((particles.points.get(1).x - particles.points.get(0).x) - (particles.points.get(1).y - particles.points.get(0).y));
		float norm_x = norm(x) ;
		
		// overlapping volume fraction
		hauteur_collision = (particles.points.get(1).y+particles.radius.get(1)) - (particles.points.get(0).y-particles.radius.get(0));
		X = particles.diff_hauteur / (particles.radius.get(0) + particles.radius.get(1)); 

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
		float terme = (7*pow((1+pow(delta,3)),0.66666666))-4*(1+pow(delta,2));
		We_reflex = (3*terme*delta*pow((1+pow(delta,3)),2)) / (pow(delta,6)*eta_1+eta_2);


		// Stretching separation
		//EQ1
		We_stretch = (4*pow(1+pow(delta,3),2) * sqrt(3*(1+delta)*(1-X)*(pow(delta,3)*phi[1]+phi[0]))) / (pow(delta,2)*((1+pow(delta,3))-(1-pow(X,2))*(phi[1]+pow(delta,3)*phi[0])));
		double f = 4*pow(1+pow(delta,3),2);
		double f1 = sqrt(3*(1+delta));
		double f2 = 1-X;
		double f3 = sqrt(pow(delta,3)*phi[1]+phi[0]);
		double s = (pow(delta,2)*((1+pow(delta,3))-(1-pow(X,2))*(phi[1]+pow(delta,3)*phi[0])));

		// Ligament instability
		beta = 3/(4*sqrt(2))*(k_1*k_2);
		r_0 = hauteur_collision * 0.5;
		We_0 = 2* r_0 * (rho/sigma) * pow(norm_vel_relative,2); //avec r_0 rayon initial egal a longueur initiale du cylindre 
	}
  


	/**
	 * Determine the type of collision et launch the associated function to update the environement
	 **/
	void calcul_We(){
		// REFLEXIVE SEPARATION
		if (We > We_reflex && We_reflex > 0 ){
			//System.err.println("Reflexive separation");
			reflexive_separation();
		}
		
		// STRETCHING SEPARATION
		else if (We > We_stretch){
			//System.err.println("Stretching separation");
			stretching_separation();
		}
		
		// COALESCENCE
		else {
			//System.err.println("Coalescence");
			coalescence();
		}
	}
  


	/**
	 * Return the name of the type of collision that will occur
   	 **/
  	String predict () {
		if (We > We_reflex && We_reflex > 0 )
			return "Reflexive separation";
		else if (We > We_stretch)
			return "Stretching separation";
		else 
			return "Coalescence";
		}
    
  

	/**
	 * Function to handle the reflexive separation case
	 **/
	void reflexive_separation(){
		// update droplets velocity
		Z = sqrt(1-(We_reflex/We));
		for (int i=0; i<2; i++){ 
			particles.velocity.get(i).x = calcul_velocity(i, particles.velocity.get(i).x, particles.velocity.get(1-i).x);
			particles.velocity.get(i).y = calcul_velocity(i, particles.velocity.get(i).y, particles.velocity.get(1-i).y);
			particles.velocity.get(i).z = calcul_velocity(i, particles.velocity.get(i).z, particles.velocity.get(1-i).z);
		}  
		particles.velocity.get(0).x = particles.velocity.get(0).x * -1;

		
		V_lig = volume(0)+volume(1);
    
		if (V_lig > 0){
			ArrayList<Float> resEq = equation_degre3(beta*sqrt(We_0), 1, 0, -1, r_0);

			float diff_0 = 10000;
			int i_tmp = 0;

			do {
				if ( r_0 > resEq.get(i_tmp) && r_0 - resEq.get(i_tmp) < diff_0 ) {
				diff_0 = r_0 - resEq.get(i_tmp);
				radius_bu = resEq.get(i_tmp);
				}
				i_tmp++;
			}
			while ( i_tmp < resEq.size());
      
			// satellite's radius
			radius_sat = 1.89 * radius_bu;
			particles.radius.add(radius_sat);
			// number of satellites with mass concervation
			float volume_sat = volume(2);
			n_sat = V_lig/volume_sat;
			n_sat = n_sat - 2;
      
			// if there is satellites
			if (n_sat >= 1) {
				for (int i=1; i<n_sat-1; i++){
					particles.radius.add(radius_sat);
				}
				for (int i=1; i<n_sat; i++){
					particles.velocity.add(velocity_sat(i));
				}
				// all droplets has the same radius
				for (int i=0; i<2; i++){
					particles.radius.set(i,radius_sat);
				}
			}else{
				particles.radius.remove(2);
			}
			particles.printMissingParticles();
		}
	}
  
  
  
  	/**
	 * Function to handle the stretching separation case
	 **/
	void stretching_separation(){
		// update droplets velocity
		Z = (X-sqrt((2.4*f_delta)/We)) / (1-sqrt((2.4*f_delta)/We));
		for (int i=0; i<2; i++){ 
			particles.velocity.get(i).x = calcul_velocity(i, particles.velocity.get(i).x, particles.velocity.get(1-i).x);
			particles.velocity.get(i).y = calcul_velocity(i, particles.velocity.get(i).y, particles.velocity.get(1-i).y);
			particles.velocity.get(i).z = calcul_velocity(i, particles.velocity.get(i).z, particles.velocity.get(1-i).z);
		}
		particles.velocity.get(0).x = particles.velocity.get(0).x * -1;


		V_lig = V_ligament();
    
		if (V_lig > 0){
			// EQ8 r_sat 
			ArrayList<Float> resEq = equation_degre3(beta*sqrt(We_0), 1, 0, -1, r_0);

			float diff_0 = 10000;
			int i_tmp = 0;
			do {
				if ( r_0 > resEq.get(i_tmp) && r_0 - resEq.get(i_tmp) < diff_0 ) {
					diff_0 = r_0 - resEq.get(i_tmp);
					radius_bu = resEq.get(i_tmp);
				}
				i_tmp++;
			}
			while ( i_tmp < resEq.size());

			// satellite radius
			radius_sat = 1.89 * radius_bu;
			particles.radius.add(radius_sat);
			// number of satellites with mass concervation
			float volume_sat = volume(2);

			n_sat = V_lig/volume_sat;
      
      		// if there is satellite
			if ( n_sat >= 1 ) {
				for (int i=1; i<=n_sat-1; i++){
					particles.radius.add(radius_sat);
				}
				for (int i=1; i<=n_sat; i++){
					particles.velocity.add(velocity_sat(i));
				}
				for (int i=0; i<2; i++){
					float new_volume = volume(i) - v_interact(i);
					particles.radius.set(i,new_radius(new_volume));
				}
			} else {
				particles.radius.remove(2);
			}
			particles.printMissingParticles();
		}
    }
  


	/**
	 * Function to handle the coalescence case
	 **/
  	void coalescence(){
		// compute the new volume
		float new_volume = volume(0) + volume(1);
		// velocity update (movement quantity conservation [kinetic energy])
		particles.velocity.get(0).x = (volume(0)*particles.velocity.get(0).x + volume(1)*particles.velocity.get(1).x) / new_volume;
		particles.velocity.get(0).y = (volume(0)*particles.velocity.get(0).y + volume(1)*particles.velocity.get(1).y) / new_volume;
		particles.velocity.get(0).z = (volume(0)*particles.velocity.get(0).z + volume(1)*particles.velocity.get(1).z) / new_volume;
		// position update (center of mass)
		particles.points.get(0).x = (volume(0)*particles.points.get(0).x + volume(1)*particles.points.get(1).x) / new_volume;
		particles.points.get(0).y = (volume(0)*particles.points.get(0).y + volume(1)*particles.points.get(1).y) / new_volume;
		particles.points.get(0).z = (volume(0)*particles.points.get(0).z + volume(1)*particles.points.get(1).z) / new_volume;
		// radius update (conservation of mass)
		particles.radius.set(0, new_radius(new_volume));
		particles.radius.set(1, 0.0);
    }
  

	
	/**
	 * Compute the relative velocity between the two main droplets
	 **/
	PVector vel_relative( ArrayList<PVector> vel){
		PVector u_j = vel.get(0).copy();
		u_j.mult(-1);
		
		PVector u_ij = vel.get(1).copy();
		u_ij.add(u_j); 
		return u_ij;
	}


  
	/**
	* Compute the norm of the vector passed by parameter
	**/
	float norm( PVector vector){
		return sqrt(pow(vector.x,2)+pow(vector.y,2)+pow(vector.z,2));
	}
  


  	/**
   	 * Resolution of the third degree equation
   	 **/
	ArrayList<Float> equation_degre3(float a, float b, float c, float d, float r_0){
		float q = (2*pow(b,3)-9*a*b*c+27*pow(a,2)*d)/(27*pow(a,3));
		float p = (3*a*c-pow(b,2))/(3*pow(a,2));
		float delta_1 = pow(q,2)+(4*pow(p,3))/27;
		float X1 = pow((-q-sqrt(delta_1))/2, 0.33333333) + pow((-q+sqrt(delta_1))/2, 0.33333333) - b/(3*a);
		float delta_2 = pow(b+a*X1,2) - 4*a*(c+(b+a*X1)*X1);
		float X2 = (-b-a*X1-sqrt(delta_2))/(2*a);
		float X3 = (-b-a*X1+sqrt(delta_2))/(2*a);
		
		float x1 = X1 * r_0;
		float x2 = X2 * r_0;
		float x3 = X3 * r_0;
		ArrayList<Float> r_bu = new ArrayList<Float>();
		if ( !Float.isNaN(x1) ) r_bu.add(x1);
		if ( !Float.isNaN(x2) ) r_bu.add(x2);
		if ( !Float.isNaN(x3) ) r_bu.add(x3);
    	return r_bu;
  	}
  


	/**
	 * Update the velocity of the two droplets in collision
	 * Equation n°3 of the paper
	 **/
	float calcul_velocity(int k, float vel_k, float vel_l){
			return (pow(particles.radius.get(k),3)*vel_k + pow(particles.radius.get(1-k),3)*vel_l - pow(particles.radius.get(1-k),3)*(vel_k-vel_l)*Z) / (pow(particles.radius.get(k),3)+pow(particles.radius.get(1-k),3));
	}
  
  

	/**
	 * Compute the volume of a droplet
	 **/
	float volume(int k){ 
		return ((4*PI)/3)*pow(particles.radius.get(k),3); 
	}
  


	/**
	 * Compute the radius of a droplet from his volume
	 **/
	float new_radius(float V){
		return pow(((V/PI)*0.75),0.33333333); //TODO : RACINE CUBIQUE..... AVEC TEST (pow(27,0.33333333)) CELA DONNE 3.0
	}
  
  

	/**
	 * Compute the volume of a part of the droplet
	 **/
	float volume_calotte(int k, float h){
		return (PI*pow(h,2)*(3*particles.radius.get(k)-h))/3;
	}
  

	/**
	 * Compute the volume involved in the collision of the two main droplets
	 **/
	float v_interact(int k){
		return volume(k)*phi[k];
	}



  	/**
	 * Compute the velocity of the satellites, with n his number
	 * Equation n°9 of the paper
	 **/
	PVector velocity_sat(int n){
		float velx = (v_interact(0)/V_lig) * particles.velocity.get(0).x + (v_interact(1)/V_lig) * (particles.velocity.get(1).x+(1-(2*n/(n_sat+1)))*vel_relative(particles.velocity).x);
		float vely = v_interact(0)/V_lig * particles.velocity.get(0).y + v_interact(1)/V_lig * (particles.velocity.get(1).y+(1-2*n/(n_sat+1))*vel_relative(particles.velocity).y);
		float velz = v_interact(0)/V_lig * particles.velocity.get(0).z + v_interact(1)/V_lig * (particles.velocity.get(1).z+(1-2*n/(n_sat+1))*vel_relative(particles.velocity).z);
		return new PVector(velx, vely, velz);
	}
  
  /**
   * Compute the volume of the ligament
   **/
	float V_ligament(){
		float tau = (1-X)*(1+delta);
		float E_dissip = 0.30; // pour l'energie cinetique initiale des gouttelettes en collision
		E_dissip = (0.5*volume(0)*pow(particles.velocity.get(0).x,2)) + (0.5*volume(1)*pow(particles.velocity.get(1).x,2));
		E_dissip += (0.5*volume(0)*pow(particles.velocity.get(0).y,2)) + (0.5*volume(1)*pow(particles.velocity.get(1).y,2));
		E_dissip += (0.5*volume(0)*pow(particles.velocity.get(0).z,2)) + (0.5*volume(1)*pow(particles.velocity.get(1).z,2));
		E_dissip *= 0.3;
		//Equation n°7 of the paper
		float E_surface = 2*sigma*sqrt(PI*volume(0)*particles.radius.get(0)*tau*(phi[0]+ pow(delta,3)* phi[1]));
		//Equation n°6 of the paper
		float E_stretch = 0.5 * rho* pow(norm_vel_relative,2)* volume(0)* (pow(delta,3)/pow(1+pow(delta,3),2))* ((1+pow(delta,3))- (1-pow(X,2))* (phi[1]+ pow(delta,3)* phi[0]) );
		float terme_1 = 0.5* rho* pow(norm_vel_relative,2)* volume(0);
		float terme_2 = (pow(delta,3)/pow(1+pow(delta,3),2));
		float terme_3 = ((1+pow(delta,3))- (1-pow(X,2))* (phi[1]+ pow(delta,3)* phi[0]));

		// Equation n°5 of the paper
		float C_sep = (E_stretch - E_surface - E_dissip) / (E_stretch + E_surface + E_dissip);
    
		// Equation n°4 of the paper
		return C_sep* (v_interact(0)+ v_interact(1));
	}
    
}
