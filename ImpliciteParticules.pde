/**
 * Class that manage the droplets.
 **/
class ImpliciteParticles {

	// the list of points in the simulation
  	ArrayList<PVector> points = new ArrayList<PVector>();
	
	// the list of velocities of the droplets
  	ArrayList<PVector> velocity = new ArrayList<PVector>();
  	
	// the list of radius of the droplets
	ArrayList<Float> radius = new ArrayList<Float>();


	// boolean to inform the state of the simulation
	boolean isRuptured = false;			// when the ligament has been broken
	boolean isCollision = false;		// when the droplets has been collided


	// variables linked to droplets' dimensions
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
	 * Iniatilize two droplets on the right and left side of the widonws, and update variables.
	 **/
	ImpliciteParticles(float dist) {
		
		// Initialisation de variables
		this.diff_hauteur = dist;
		posX1 = 0 - espace_begin;
		posX2 = 0 + espace_begin;
		tailleZ = width-6;

		// default init :
		//	- left droplet if smaller than the right one (always)
		//  - values by default are the smalest one in the IHM
		// left droplet
		radius.add(0.3);
    	points.add(new PVector(posX1, posY+diff_hauteur*0.5, tailleZ/2));
		velocity.add( new PVector(0.25, 0, 0));
		// right droplet
		radius.add(0.3);
		points.add(new PVector(posX2, posY - diff_hauteur*0.5, tailleZ/2));	
		velocity.add( new PVector(-0.25, 0, 0));


		// init the minimum and maximun Y to reduce calculation in the marching
		for (int i = 0; i < N-1; i++){
			if (minimumY > points.get(i).y-(6+1)*29)
				minimumY = points.get(i).y-(6+1)*29;
				
			if (maximumY < points.get(i).y+(6+1)*29)
				maximumY = points.get(i).y+(6+1)*29;
		}
	}
  

  
	/**
	 * Create a new droplet in the liste "point" to display an extra droplet in the simulation.
	 * The lists velocity and radius must already have the value for the new droplets to display.
	**/
	void printMissingParticles () {
		// compute the shift in Y to apply for each droplet to place them like a cascade
		float decalage_y = ((particles.points.get(1).y+particles.radius.get(1)) - (particles.points.get(0).y-particles.radius.get(0)))/(particles.velocity.size()-2);
		float y = 0;
	
		// for each new droplet (satelites)
		for ( int i = 2; i < velocity.size(); i++ ) {
			// prevent the lack of informations for the new droplets
			if ( velocity.size() <= points.size() && radius.size() <= points.size() ) {
				System.err.println("printMissingParticles() : Les listes Velocity ou Radius ne possèdent pas les informations pour créer une nouvelles goutte");
			}else{
				// case where there is no y difference for the two main droplet
				if (particles.points.get(0).y == particles.points.get(1).y)
					// then no shift to apply
					points.add(new PVector(points.get(0).x, particles.points.get(0).y, points.get(0).z)); 
				else
					points.add(new PVector(points.get(0).x, particles.points.get(0).y-radius.get(0) + y, points.get(0).z)); 
				y += decalage_y;
				// increase the number of droplet in the simulation
				N ++;
			}
		}
  	}
  


	/**
	 * Apply the mouvement for each droplet. If droplets has been collided but the ligament did not break, then the ligament is
	 * simulated by placing satelites in the space with a regular shift. If the marching display mode is disabled, then the
	 * droplets are printed with simple sphere
	 **/
	void nextStep() {
		// apply velocity on the droplets
		for (int i = 0; i < N; i++) {
			if ( i < N ) {
				if ( (i < 2 || isRuptured ) && run ) {
				PVector p = points.get(i);
				p.add(velocity.get(i));
				}
			}
		}

		// case where there is a ligament to display (collision done but no ligament rupture)
		if ( isCollision && !isRuptured ) {
			for ( int i = 0; i < 2; i++ ) {
				// if marching display is disabled
				if (print_sphere){
					translate(particles.points.get(i).x, particles.points.get(i).y*zoom, 0);
					sphere(particles.radius.get(i)*zoom);
					translate(-particles.points.get(i).x, -particles.points.get(i).y*zoom, 0);
				}
			}

			// distance between the two droplets divide by 2
			float demi_distance = (points.get(0).x-points.get(1).x-(radius.get(0)+radius.get(1))*0.5)*0.5;
		

			////////////////////////////////////////////////////////
			// CASE : number of satelite odd
			////////////////////////////////////////////////////////
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
			
				// other satelites
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


			////////////////////////////////////////////////////////
			// CASE : number of satelite even
			////////////////////////////////////////////////////////
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
			//float decalage = demi_distance / (N-3);


		// in case of no collision or broken ligament
		}else{
			// then display normaly the droplets
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
	}


  
	/**
	 * Method to evaluate at a cube of the marching tetrahedron division, the presence of a droplet to print. A multiplication
	 * is applied to the droplets coordinates to make them bigger is the simulation without change the environement (distance, size...)
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
			if ( radius.get(r) < 0.3 ) zoom_sat = 3;		// apply a zoom for too small satelites to see them
			PVector tmp = new PVector(points.get(r).x, points.get(r).y*zoom_tetra*zoom_sat, points.get(r).z);
			v += f(tmp.dist(n)/(radius.get(r)*zoom_tetra*zoom_sat));
		}
		return v;
	}



  	/**
	 * Evaluate the vertices (0 or 1)
   	 **/
	int evalInt(float i, float j, float k) {
		return (eval(i, j, k) >= 0.0001) ? 1 : 0;
	}
  


	/**
	 * Launch the recursivity
	 **/
  	PVector interp(float x1, float y1, float z1, float x2, float y2, float z2) {
    	return interpRec(x1, y1, z1, x2, y2, z2, REC_INTERP);
  	}



	/**
	 * Check of the ligament is broken or not. It is if the sum of diameter or each satelite is smaller
	 * than the space between the two main droplets to fill.
	 **/
  	boolean isRupture () {
    	boolean ret = false;

		if ( N > 2 && isCollision ) {
			float distance = 0;
			// is marching printing is disabled
			if ( print_sphere ) {
				distance = resolv.norm(new PVector(points.get(0).x-points.get(1).x,0,0)) - radius.get(0) - radius.get(1);
				if ( distance > radius.get(2)*2*(N-2)*zoom ) {
				isRuptured = true;
				System.err.println("Rupture");
				}
			
			// else marching printing is enabled
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



	PVector interpRec(float x1, float y1, float z1, float x2, float y2, float z2, int level) {
		PVector p1 = new PVector(x1, y1, z1);
		PVector p2 = new PVector(x2, y2, z2);	
		PVector m = PVector.lerp(p1, p2, 0.5);	// lerp to find the middle between p1-p2 (middle of a pixel)
		float v1 = evalInt(x1, y1, z1);
		float v2 = evalInt(x2, y2, z2);
		if ((level == 0) || (v1 == v2)) 
			return m;
		if (evalInt(m.x, m.y, m.z) == v1) return interpRec(m.x, m.y, m.z, x2, y2, z2, level-1);
			return interpRec(x1, y1, z1, m.x, m.y, m.z, level-1);
		// we stop the recusion if the points at the left and right (corner of a pixel) has the same value
	} 
}
