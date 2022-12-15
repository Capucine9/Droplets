/**
 * Main file of the project
 */


// the class which manage the interface of the simulation
IHM ihm;

// class to solve equations to predict droplet behaviour
Resolving_interaction resolv;

// boolean to let run the simaltion, or stop it
boolean run = false;
// boolean to inform the way of droplet printing
boolean print_sphere = false;



/**
 * Called after create the windows
 **/
void settings() {
  	size((int)(1280*1),(int)(720*1), P3D);
}



/**
 * Called at the creation od the windows. Inizialize variables, commands and the environment
 * to set up the simulation.
 **/
void setup() {
	// create and initialize the marching tetrahedron engine (for displaying)
  	initMarching();
	// initializing the solver
  	resolv = new Resolving_interaction();
  	resolv.init();
	// initialize the ihm
  	ihm = new IHM(this,resolv);
  
	// prevent no collision if the droplets are to espaced (in y)
  	if ( particles.radius.get(0)+particles.radius.get(1) <= particles.diff_hauteur)
    	System.err.println("No collision");
}



/**
 * Called each frame of the windows. Change droplets position according to their velocity,
 * check collision to ask the solver to predict the result, and display ihm and droplets.
 **/
void draw() {
  	background(75);
  
	// put the origin at the center of the screen (a little below)
	translate(width/2, height/1.7, -width/2);
	lights();
	noFill();
	fill(3, 161, 252);
	
	particles.nextStep();		//	move droplets (and launch basic display if marching is disabled)
	if ( !print_sphere )		// if marching printing is on
		frameMarching();		// display droplets with the marching algorithm
	stroke(30,180,255);
	
	fill(255);
	// cancel shift 
	translate(-width/2, -height/1.7, width/2);

  
	// Get value selected by the user in the IHM
	int distance_selected = Integer.parseInt(""+ihm.distance.getItem((int)ihm.distance.getValue()).get("text"));
	float size_0_selected = Float.parseFloat(""+ihm.L_diameter.getItem((int)ihm.L_diameter.getValue()).get("text"));
	float speed_0_selected = Float.parseFloat(""+ihm.L_velocity.getItem((int)ihm.L_velocity.getValue()).get("text"));
	float size_1_selected = Float.parseFloat(""+ihm.R_diameter.getItem((int)ihm.R_diameter.getValue()).get("text"));
	float speed_1_selected = Float.parseFloat(""+ihm.L_velocity.getItem((int)ihm.R_velocity.getValue()).get("text"));
	
	// if the value has change the value of a filed
	if ( distance_selected != ihm.diff_hauteur_selected ||
		size_0_selected != ihm.droplet_radius_selected[0] ||
		speed_0_selected != ihm.droplet_speed_selected[0] || 
		size_1_selected != ihm.droplet_radius_selected[1] ||
		speed_1_selected != ihm.droplet_speed_selected[1] ) {

		// then re-init the simulation
		ihm.diff_hauteur_selected = distance_selected;
		distance = (ihm.diff_hauteur_selected*0.01)*(size_0_selected+size_1_selected);
		initMarching();
    
		// save selected value
		if ( !particles.isCollision ) particles.radius.set(0, size_0_selected);
		if ( !particles.isCollision ) particles.velocity.set(0, new PVector(speed_0_selected,particles.velocity.get(0).y,0));
		ihm.droplet_radius_selected[0] = size_0_selected;
		ihm.droplet_speed_selected[0] = speed_0_selected;
		if ( !particles.isCollision ) particles.radius.set(1, size_1_selected);
		if ( !particles.isCollision ) particles.velocity.set(1, new PVector(-speed_1_selected,particles.velocity.get(1).y,0));
		ihm.droplet_radius_selected[1] = size_1_selected;
		ihm.droplet_speed_selected[1] = speed_1_selected;
	
		resolv = new Resolving_interaction();
		resolv.init(); 
		ihm.resolv = resolv;
	}
  	print_sphere = !ihm.checkbox.getState(0);
  

	// if the lecture mode is enable
	if ( run ) {    
		// if there is collision
		if ( particles.radius.get(0)+particles.radius.get(1) > particles.diff_hauteur) {
		if (particles.points.get(1).x - particles.points.get(0).x < 0 && particles.isCollision == false){
			particles.isCollision = true;
			resolv.calcul_We(); 		// prediction
		}
		// detect the moment where the ligament breaks
		if ( !particles.isRuptured && particles.isRupture() )
			particles.isRuptured = true;
		}
	}  
  
  // print the interface
  ihm.printInterface();
}



/**
 * Method executed when the user click on a mouse button
 **/
void mousePressed() {
	// at the right mouse button click, change the lecture mode
	if ( mouseButton == RIGHT ) run = !run;
}
