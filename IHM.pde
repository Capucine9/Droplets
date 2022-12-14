import controlP5.*;

/**
 * Class that manage the interface of the simulation
 **/
class IHM {
  
  ControlP5 controlP5;
  Resolving_interaction resolv;
  
  DropdownList distance; 
  CheckBox checkbox;
  
  DropdownList L_diameter; 
  DropdownList R_diameter;
  DropdownList L_velocity; 
  DropdownList R_velocity; 
  
  int margeLeft = 10;
  
  int[] diffHeight = {0, 25, 50, 75, 80, 90};
  int diffHeight_selected = 0;
  
  float[] speed = {0.25, 0.5, 1, 2};
  int speed_selected = 0;
  
  float[] diameter = {0.3, 0.4, 0.5, 0.8};
  int L_diameter_selected = 0;
  int R_diameter_selected = 0;
  
  float[] velocity = {0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4, 6, 10, 15, 20};
  int L_velocity_selected = 0;
  int R_velocity_selected = 0;
  
  
  int level_height = 30;
  int level_side = 80;
  int level_diameter_velocity = 120;
  
  
  
  float diff_hauteur_selected = 0.0;
  float[] droplet_speed_selected = {0.0,0.0};
  float[] droplet_radius_selected = {0.0,0.0};
    
  
  IHM(PApplet app, Resolving_interaction r) {
    this.controlP5 = new ControlP5(app);
    this.resolv = r;
    
    // initialize the distance dropdownlist
    distance = controlP5.addDropdownList("...",550,level_height-15,100,300);
    distance.setFont(createFont("arial",15));
    customizeDistance();
    
    
    
    // customize diameter and velocity dropdownlist
    L_diameter = controlP5.addDropdownList("...  ",margeLeft*4,150,100,300);
    R_diameter = controlP5.addDropdownList("...   ", width-360+margeLeft*4,150,100,300);
    L_velocity = controlP5.addDropdownList("...    ",margeLeft*4+200,150,100,300);
    R_velocity = controlP5.addDropdownList("...     ", width-160+margeLeft*4,150,100,300);
    
    L_diameter.setFont(createFont("arial",15));
    R_diameter.setFont(createFont("arial",15));
    L_velocity.setFont(createFont("arial",15));
    R_velocity.setFont(createFont("arial",15));
    
    customizeDiameter(L_diameter);
    customizeDiameter(R_diameter);
    
    customizeVelocity(L_velocity);
    customizeVelocity(R_velocity);
    
    controlP5.setFont(createFont("arial",15));
    checkbox = controlP5.addCheckBox("checkBox")
                .setPosition(width-250, 10)
                .setSize(40, 40)
                .setItemsPerRow(3)
                .setSpacingColumn(30)
                .setSpacingRow(20)
                .addItem("Marching Tetrahedron", 0);
    
  }
  
  
  
  /**
   * Customize diameter dropdownlist
   **/
  void customizeDiameter(DropdownList ddl) {
    ddl.setBackgroundColor(color(190));
    ddl.setItemHeight(30);
    ddl.setBarHeight(40);
    for ( int i = 0; i < diameter.length; i++ ) {
       ddl.addItem(diameter[i]+"", i+1);
     }
    ddl.setColorBackground(color(60));
    ddl.setColorActive(color(255,128));
    ddl.setValue(0);
    ddl.setOpen(false);
    ddl.setLabel(diameter[0]+"");
  } 
  

  
  /**
   * Customize velocity dropdownlist
   **/
  void customizeVelocity(DropdownList ddl) {
    ddl.setBackgroundColor(color(190));
    ddl.setItemHeight(30);
    ddl.setBarHeight(40);
    for ( int i = 0; i < velocity.length; i++ ) {
       ddl.addItem(velocity[i]+"", i+1);
     }
    ddl.setColorBackground(color(60));
    ddl.setColorActive(color(255,128));
    ddl.setValue(0);
    ddl.setOpen(false);
    ddl.setLabel(velocity[0]+"");
  }
  
  

  /**
   * Customize distance dropdownlist
   **/
  void customizeDistance() {
    this.distance.setBackgroundColor(color(190));
    this.distance.setItemHeight(30);
    this.distance.setBarHeight(40);
    for ( int i = 0; i < diffHeight.length; i++ ) {
      this.distance.addItem(diffHeight[i]+"", i+1);
    }
    this.distance.setColorBackground(color(60));
    this.distance.setColorActive(color(255,128));
    this.distance.setValue(0);
    this.distance.setOpen(false);
    this.distance.setLabel(diffHeight[0]+"");
  }
  
  
  
  void printInterface() {  
    
    // ======================================================================================
    // Height difference between the two droplets
    // ======================================================================================
    textSize(25);
    text("Distance de hauteur des centres des gouttes (%) :", margeLeft, level_height);  // * 2.5 to avoid text in the windows bar
   
    
    // ======================================================================================
    // Left dorplet
    // ======================================================================================
    text("Goutte Gauche", margeLeft, level_side);
    text("Diametre (cm) :", margeLeft, level_diameter_velocity);
    text("Vitesse (m/s) :", margeLeft+200, level_diameter_velocity);
    
        
    // ======================================================================================
    // Right dorplet
    // ======================================================================================
    text("Goutte Droite", width-150, level_side);
    text("Diametre (cm) :", width-360, level_diameter_velocity);
    text("Vitesse (m/s) :", width-160, level_diameter_velocity);
    

    // ======================================================================================
    // Collision type
    // ======================================================================================
    try{
      text("Type : "+resolv.predict(),margeLeft, height-30);
    }catch(Exception e){}
  }
}
