/**
 * File that manage the marching tetrahedron 
 **/

// the resolution of the grid
final int rez = 15;
// number of droplets
int N = 2;
// number of recursive calls
final int REC_INTERP = 4;

// number of columns of the grid
int cols;
ImpliciteParticles particles;
Resolving_interaction interactions;

// X distance between the two main droplets, change by the IHM and use to save is value in ImpliciteParticles.pde
float distance = 0;


float f(float d) {
  return (d < 1) ? 8.0*pow(1 - d*d,2)/9 : 0;
}



/**
 * Init the algorithm
 **/
void initMarching() {
	N = 2;
	particles = new ImpliciteParticles(distance);
	cols = 1 + width / (2*rez);
}



/**
 * Triangle to print
 **/
void triangl(PVector v1, PVector v2, PVector v3) {
	beginShape();
	vertex(v1.x, v1.y, v1.z);
	vertex(v2.x, v2.y, v2.z);
	vertex(v3.x, v3.y, v3.z);
	endShape();
}



/**
 * Polygon to print
 **/
void poly(PVector v1, PVector v2, PVector v3, PVector v4) {
	beginShape();
	vertex(v1.x, v1.y, v1.z);
	vertex(v2.x, v2.y, v2.z);
	vertex(v3.x, v3.y, v3.z);
	endShape();
	beginShape();
	vertex(v1.x, v1.y, v1.z);
	vertex(v4.x, v4.y, v4.z);
	vertex(v3.x, v3.y, v3.z);
	endShape();
}



/**
 * Method that launch the course of the grid to print the droplets
 **/
void frameMarching() {
	for (int i = -cols+1; i < cols-1; i++) {
		for (int j = (int)particles.minimumY/rez; j < (int)particles.maximumY/rez; j++) {
			for (int k = (int)(particles.points.get(0).z-(particles.radius.get(0)+1)*particles.zoom_tetra)/rez; k <= (int)(particles.points.get(0).z+(particles.radius.get(0)+1)*particles.zoom_tetra)/rez; k++){
				
				float x = i * rez;
				float y = j * rez;
				float z = k * rez;
				
				PVector de = particles.interp(x, y, z+rez, x, y+rez, z); //de
				PVector dh = particles.interp(x, y, z+rez, x, y+rez, z+rez); //dh
				PVector dc = particles.interp(x, y, z+rez, x+rez, y, z+rez); //dc
				PVector ec = particles.interp(x, y+rez, z, x+rez, y, z+rez); //ec
				PVector he = particles.interp(x, y+rez, z+rez, x, y+rez, z); //he
				PVector hc = particles.interp(x, y+rez, z+rez, x+rez, y, z+rez); //hc
				
				PVector ab = particles.interp(x, y, z, x+rez, y, z); //ab
				PVector bc = particles.interp(x+rez, y, z, x+rez, y, z+rez); //bc
				PVector ac = particles.interp(x, y, z, x+rez, y, z+rez); //ac
				PVector ae = particles.interp(x, y, z, x, y+rez, z); //ae
				PVector eb = particles.interp(x, y+rez, z, x+rez, y, z); //eb
				
				PVector da = particles.interp(x, y, z+rez, x, y, z); //da
				
				PVector hg = particles.interp(x, y+rez, z+rez, x+rez, y+rez, z+rez); //hg
				PVector cg = particles.interp(x+rez, y, z+rez, x+rez, y+rez, z+rez); //cg
				PVector ge = particles.interp(x+rez, y+rez, z+rez, x, y+rez, z); //ge
				
				PVector bf = particles.interp(x+rez, y, z, x+rez, y+rez, z); //bf
				PVector ef = particles.interp(x, y+rez, z, x+rez, y+rez, z); //ef
				PVector fc = particles.interp(x+rez, y+rez, z, x+rez, y, z+rez); //fc
				
				PVector fg = particles.interp(x+rez, y+rez, z, x+rez, y+rez, z+rez); //fg
				
				// tetraedron evaluation
				int state1 = getState(  particles.evalInt(x, y+rez, z+rez), //cdeh
									particles.evalInt(x, y+rez, z),
									particles.evalInt(x, y, z+rez),
									particles.evalInt(x+rez, y, z+rez));
									
				int state2 = getState(  particles.evalInt(x, y+rez, z), //abce
									particles.evalInt(x+rez, y, z+rez),
									particles.evalInt(x+rez, y, z),
									particles.evalInt(x, y, z));
									
				int state3 = getState(  particles.evalInt(x, y+rez, z), //acde
									particles.evalInt(x, y, z+rez),
									particles.evalInt(x+rez, y, z+rez),
									particles.evalInt(x, y, z));
				
				int state4 = getState(  particles.evalInt(x, y+rez, z+rez), //cegh
									particles.evalInt(x+rez, y+rez, z+rez),
									particles.evalInt(x, y+rez, z),
									particles.evalInt(x+rez, y, z+rez));
									
				int state5 = getState(  particles.evalInt(x+rez, y+rez, z), //bcef
									particles.evalInt(x, y+rez, z),
									particles.evalInt(x+rez, y, z+rez),
									particles.evalInt(x+rez, y, z));
									
				int state6 = getState(  particles.evalInt(x+rez, y+rez, z+rez), //cefg
									particles.evalInt(x+rez, y+rez, z),
									particles.evalInt(x, y+rez, z),
									particles.evalInt(x+rez, y, z+rez));
								
				noStroke();
				
				// print triangle and polygon according to the state of the tetrahedron
				switch (state1) { //cdeh
					case 1: //c different
					case 14:
						triangl(dc, ec, hc);
						break;
					case 2: //d different
					case 13:
						triangl(de, dh, dc);
						break;
					case 3: //cd different
					case 12:
						poly(dh, de, ec, hc);
						break;
					case 4: //e different
					case 11:
						triangl(de, ec, he);
						break;
					case 5: //ce different
					case 10:
						poly(dc, hc, he, de);
						break;
					case 6: //de different
					case 9:
						poly(he, dh, dc, ec);
						break;
					case 7:  //h different
					case 8:
						triangl(dh, he, hc);
						break;
				}
				
				switch (state2) { //abce
					case 1: //a different
					case 14:
						triangl(ab, ac, ae);
						break;
					case 2: //b different
					case 13:
						triangl(ab, bc, eb);
						break;
					case 3: //ab different
					case 12:
						poly(ae, ac, bc, eb);
						break;
					case 4: //c different
					case 11:
						triangl(bc, ac, ec);
						break;
					case 5: //ac different
					case 10:
						poly(ae, ab, bc, ec);
						break;
					case 6: //bc different
					case 9:
						poly(eb, ab, ac, ec);
						break;
					case 7:  //e different
					case 8:
						triangl(ae, eb, ec);
						break;
				}
				
				switch (state3) { //acde
					case 1: //a different
					case 14:
						triangl(da, ac, ae);
						break;
					case 2: //c different
					case 13:
						triangl(ac, ec, dc);
						break;
					case 3: //ac different
					case 12:
						poly(da, ae, ec, dc);
						break;
					case 4: //d different
					case 11:
						triangl(da, dc, de);
						break;
					case 5: //ad different
					case 10:
						poly(de, dc, ac, ae);
						break;
					case 6: //cd different
					case 9:
						poly(da, de, ec, ac);
						break;
					case 7:  //e different
					case 8:
						triangl(ec, ae, de);
						break;
				}
				
				switch (state4) { //cegh
					case 1: //c different
					case 14:
						triangl(cg, ec, hc);
							break;
					case 2: //e different
					case 13:
						triangl(he, ge, ec);
						break;
					case 3: //ce different
					case 12:
						poly(he, ge, cg, hc);
						break;
					case 4: //g different
					case 11:
						triangl(hg, cg, ge);
						break;
					case 5: //cg different
					case 10:
						poly(hg, ge, ec, hc);
						break;
					case 6: //eg different
					case 9:
						poly(hg, cg, ec, he);
						break;
					case 7:  //h different
					case 8:
						triangl(hg, he, hc);
						break;
				}
				
				switch (state5) { //bcef
					case 1: //b different
					case 14:
						triangl(bc, eb, bf);
						break;
					case 2: //c different
					case 13:
						triangl(bc, fc, ec);
						break;
					case 3: //bc different
					case 12:
						poly(bf, eb, ec, fc);
						break;
					case 4: //e different
					case 11:
						triangl(ef, eb, ec);
						break;
					case 5: //be different
					case 10:
						poly(ec, ef, bf, bc);
						break;
					case 6: //ce different
					case 9:
						poly(ef, eb, bc, fc);
						break;
					case 7:  //f different
					case 8:
						triangl(fc, ef, bf);
						break;
				}
				
				switch (state6) { //cefg
					case 1: //c different
					case 14:
						triangl(ec, fc, cg);
						break;
					case 2: //e different
					case 13:
						triangl(ge, ef, ec);
						break;
					case 3: //ce different
					case 12:
						poly(ge, ef, fc, cg);
						break;
					case 4: //f different
					case 11:
						triangl(fg, ef, fc);
						break;
					case 5: //cf different
					case 10:
						poly(fg, ef, ec, cg);
						break;
					case 6: //ef different
					case 9:
						poly(ge, ec, fc, fg);
						break;
					case 7:  //g different
					case 8:
						triangl(fg, ge, cg);
						break;
				}
			}
		}
	}
}

/**
 * a, b, c and d egals to 0 or 1 according to the color
 **/
int getState(int a, int b, int c, int d) {
  	return a * 8 + b * 4  + c * 2 + d * 1;
}
