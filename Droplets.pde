//// Physically-Based Droplet Interaction


// Variables

// Parametres de la gouttelette i
float u_i; // velocite
float r_i; // rayon 
float phi_i; // fraction de volume qui se chevauche

// Parametres de la gouttelette j
float u_j; // velocite 
float r_j; // rayon
float phi_j; // fraction de volume qui se chevauche

// Parametres commums aux deux gouttelettes
float rho; // densite
float sigma; // tension superficielle
float ksi; // 0.5*X*(1+delta)
float eta_1; // 2*(1-ksi)^2*(1-ksi^2)^(1/2)-1
float eta_2; // 2*(delta-ksi)^2*(delta^2-ksi^2)^(1/2)-delta^3



float u_ij; // velocite relative (u_j - u_i)

// Weber Number (rapport entre les forces d'inertie et les forces de tension superficielle)
float We; // (2 * rho * ||u_ij||²) / sigma

// Impact (plus courte distance entre les deux centres des gouttelettes orthogonaux à la vitesse relative)
float X; // x / (r_i + r_j) = (||xij - xij^T * û * û||) / (r_i + r_j) avec xij = xj - xi et û = u_ij / ||u_ij||

// Size ratio
float delta; // r_j / r_i

float f_delta; // delta^(-3)-2.4*delta^(-2)+2.7*delta^(-1)


//// Stretching separation
//EQ1
float We_stretch; // (4*(1+delta^3)^2 *(3*(1+delta)*(1-X)*(delta^3*phi_j+phi_i))^(1/2)) / (delta^2*((1+delta^3)-(1-X^2)*(phi_j+delta^3*phi_i)))


//// Reflexive separation
//EQ2
float We_reflex; // (3*(7*(1+delta^3)^(2/3)-4*(1+delta^2))*delta*(1+delta^3)^2) / (delta^6*eta_1+eta_2)



// Actualisation de la velocite des gouttelettes
//valable pour k,l = i,j
//EQ3
float u; //^~k; // (r_k^3*u_k+r_l^3*u_l-r_l^3*(u_k-u_l)*z) / (r_k^3+r_l^3)
float  z; // si stretching separation : (X-((2.4*f_delta)/We)^(1/2)) / (1-((2.4*f_delta)/We)^(1/2))
// si reflexive separation : (1-We_reflex/We)^(1/2)


//// Volume des ligaments liquide
//EQ6
float E_stretch; // 1/2*rho*||u_ij||^2*V_i*(delta^3/((1+delta^3)^2))*((1+delta^3)-(1-X^2)*(phi_j+delta^3*phi_i))
float tau; // (1-X)*(1+delta)
//EQ7
float E_surface; // 2*sigma*(pi*V_i*r_i*tau*(phi_i+delta^3*phi_j))^(1*2)
float E_dissip; // 0.30 pour l'energie cinetique initiale des gouttelettes en collision
//EQ5
float C_sep; // coeff de separation de volume (E_stretch-E_surface-E_dissip) / (E_stretch+E_surface+E_dissip)

// Stretching separation
//EQ4
float V_lig; // volume du ligament liquide resultant C_sep*(Vinteract_i+Vinteract_j) avec Vinteract_k = V_k*phi_k
// V_k=4/3*pi*r_k^3 

// Reflexive separation
// V_lig = V_i+V_j;


//if V_lig >0 : ligaments et satellites


//// Instabilite des ligaments
float k_1; // 11.5
float k_2; // 0.45
float beta; // 3/(4*2^(1/2))*(k_1*k_2)
float We_0; //2*r_0*rho/sigma*||u_ij||^2 avec r_0 rayon initial egal a longueur initiale du cylindre 
float r_bu; // rayon de rupture du cylindre calculable avec :
//EQ8
// beta*We_0^(1/2)*(r_bu/r_0)^(7/2)+(r_bu/r_0)^2-1=0

float r_sat; // rayon des satellites 1.89*r_bu

//// Nombre de satellites
// Stretching separation

// Reflexive separation



//// Vitesse des satellites
// for n=1,....n_sat
//EQ9
// u_n=V_lig_i/V_lig*u~_i+V_lig_j/V_lig*(u~j+(1-2n*(n_sat+1))*u~_ij)


//5.2.3
//5.3.1






IHM ihm;


void settings() {
  size((int)(1280*1.3),(int)(720*1.3));
}



void setup() {
  
   ihm = new IHM(this);
  //p2 = controlP5.addDropdownList("myList-p2",220,100,100,120);
  //customize(p2); 
}


void draw() {
  background(75);
  
  
  ihm.printInterface();
  
  
  // get distance
  int variable = Integer.parseInt(""+ihm.distance.getItem((int)ihm.distance.getValue()).get("text"));
  System.out.println(variable);
  // get speed
  float variable2 = Float.parseFloat(""+ihm.speed_ddl.getItem((int)ihm.speed_ddl.getValue()).get("text"));
  System.out.println(variable2);
  // get size Ball 1
  variable2 = Float.parseFloat(""+ihm.L_diameter.getItem((int)ihm.L_diameter.getValue()).get("text"));
  System.out.println(variable2);
  // get speed Ball 1
  variable2 = Float.parseFloat(""+ihm.L_velocity.getItem((int)ihm.L_velocity.getValue()).get("text"));
  System.out.println(variable2);
  // get size Ball 2
  variable2 = Float.parseFloat(""+ihm.R_diameter.getItem((int)ihm.R_diameter.getValue()).get("text"));
  System.out.println(variable2);
  // get speed Ball 2
  variable2 = Float.parseFloat(""+ihm.R_velocity.getItem((int)ihm.R_velocity.getValue()).get("text"));
  System.out.println(variable2);
  
  // Structure (keys) of DropdownList.getItem()
  // view
  // color
  // name
  // text
  // state
  // value
}


void controlEvent(ControlEvent theEvent) {
  // DropdownList is of type ControlGroup.
  // A controlEvent will be triggered from inside the ControlGroup class.
  // therefore you need to check the originator of the Event with
  // if (theEvent.isGroup())
  // to avoid an error message thrown by controlP5.

  if (theEvent.isGroup()) {
    // check if the Event was triggered from a ControlGroup
    println("event from group : "+theEvent.getGroup().getValue()+" from "+theEvent.getGroup());
  } 
  else if (theEvent.isController()) {
    println("event from controller : "+theEvent.getController().getValue()+" from "+theEvent.getController());
  }
}
