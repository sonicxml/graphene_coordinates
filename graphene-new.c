/* 
 * Tool to take an input of height and width and
 * generate a coordinate list of atoms in a sheet 
 * of graphene with aforementioned dimensions.
 *
 * Trevin Gandhi
 */

#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <stdbool.h>
#include "graphene.h"

// Difference Numbers: the hexagon won't always go right up to the 
// nanometers specified, so there are gaps created when translating. 
// Solving for fxdiff and fydiff eliminate those.
static unsigned int xdiff, ydiff;
static float fxdiff, fydiff;
// Hexagon Cut Horizontal
static unsigned int hexagonx;
// Hexagon Cut Vertical
static unsigned int hexagony;
// Height tracker
static double hunits = 0;
// Width tracker
static double wunits;
// Limits for the loops
static double ylimit, xlimit;
// Atom counters
static unsigned int acx = 0, acy = 0;
static unsigned int atoms = 0;
static unsigned int translations = 0;
static float distance;

int HexagonLoop () {
	// Add in height values
	hexagony = (y/hexy);
 	if ((hexagony * hexy) <= (y - 8.608292511)) {
		hexagony++;
	}
	// Add in width values
	hexagonx = (x/hexx);
	if ((hexagonx * hexx) <= (x - 5.629024292)) {
		hexagonx++;
	}
	// Open File
	FILE *file;
	file = fopen("graphenecoordinates.txt", "w");

	// Close the file
	fclose(file);

	return 0;
}

int RectLoop () {
	float transx, transy, i, iy = 1, ix = 1;
	unsigned int oppx, oppy;	
	bool cut;

	// Open File
	FILE *file;
	file = fopen("graphenecoordinates.txt", "w");
	
	// Convert to angstroms
	rectx *= 10;
	recty *= 10;
	recth *= 10;
	rectw *= 10;

	// Get Upper Left Y value of rectangle and Bottom Right X value of rectangle
	oppx = rectx + rectw;
	oppy = recty + recth;
	
	while (hunits <= ylimit) {
		for (wunits = (pointy>0?wleg:0); wunits <= xlimit; wunits += dwleg) {
			for (i = 0; i < 2; i += 2) {
				if ((hunits >= ((i * recth) + recty) && hunits <= ((i * recth) + oppy)) 
					&& (wunits >= (rectx) && wunits <= (oppx))) {
					cut = true;
					break;						
				} else {
					cut = false;
				}
			}
			if (!cut) {
				fprintf(file, "%f	%f	0\n", wunits, hunits);
				atoms++;
				if (ytrans == 0 && xtrans != 0) {
					for (ix = 1; ix <= xtrans; ix++) {
						transx = wunits + (ix * (xlimit - fxdiff));
						fprintf(file, "%f	%f	0\n", transx, hunits);
					}
				} else {
					for (iy = 1; iy <= ytrans; iy++) {
						transy = hunits + (iy * (ylimit - fydiff));
						for (ix = 1; ix <= xtrans; ix++) {
							transx = wunits + (ix * (xlimit - fxdiff));
							fprintf(file, "%f	%f	0\n", wunits, transy);
							fprintf(file, "%f	%f	0\n", transx, hunits);
							fprintf(file, "%f	%f	0\n", transx, transy);
						}
					}				
				}
			}
			
			if (hunits == 0)
				acx++;
		}
		
		acy++;
		if (hunits == 0) {
			pointy = -1;
			hunits += hleg;
		} else {
			pointy += sign(pointy);
			if ((abs(pointy)) > 2) 
				pointy = -sign(pointy);
			hunits += hleg * (abs(pointy) > 1 ? 2 : 1);
		}

	}

	
	// Close the file
	fclose(file);

	return 0;
}

int StandardLoop () {

	float transx, transy, xone, yone, xtwo, ytwo, distx, disty, iy = 1, ix = 1;

	// Open File
	FILE *file;
	file = fopen("graphenecoordinates.txt", "w");

	while (hunits <= ylimit) {
		for (wunits = (pointy>0?wleg:0); wunits <= xlimit; wunits += dwleg) {
			fprintf(file, "%f	%f	0\n", wunits, hunits);
			atoms++;
			if (ytrans == 0 && xtrans != 0) {
				for (ix = 1; ix <= xtrans; ix++) {
					if (wunits != 0) {
						// transx = (wunits) + (ix * x) - (ix * fxdiff);
						transx = wunits + (ix * (xlimit - fxdiff));			
						fprintf(file, "%f	%f	0\n", transx, hunits);
						translations++;
						if (wunits == wleg && hunits == 0) {
							xone = transx;
							yone = hunits;
						}
						if (wunits >= (xlimit - dwleg) && hunits == hleg && wunits <= (xlimit - dwleg + x) && ix == 1) {
							xtwo = wunits;
							ytwo = hunits;
						}
					}
				}			
			} else {
				for (iy = 1; iy <= ytrans; iy++) {
					transy = hunits + (iy * (ylimit - fydiff));
					for (ix = 1; ix <= xtrans; ix++) {
						transx = wunits + (ix * (xlimit - fxdiff));	
						fprintf(file, "%f	%f	0\n", wunits, transy);
						fprintf(file, "%f	%f	0\n", transx, hunits);
						fprintf(file, "%f	%f	0\n", transx, transy);
					}
				}				
			}
			
			if (hunits == 0)
				acx++;				
		}
	
		acy++;
		if (hunits == 0) {
			pointy = -1;
			hunits += hleg;
		} else {
			pointy += sign(pointy);
			if ((abs(pointy)) > 2) 
				pointy = -sign(pointy);
			hunits += hleg * (abs(pointy) > 1 ? 2 : 1);
		}

	}

	distx = xtwo - xone;
	disty = ytwo - yone;
	distx *= distx;
	disty *= disty;
	distance = distx + disty;
	distance = sqrtf(distance);
	
	// Close the file
	fclose(file);

	return 0;
}	
	
int main () {
	
	if (units) {
		// Since x and y are in nm, convert to angstroms
		x *= 10;
		y *= 10;
	}
	
	// Limits for the loops
	// Limits if x and y are hexagons
	// ylimit = (y * ((2 * hleg) + 1.42)) + 1;
	// xlimit = (x * dwleg) + 1;
	// Limits if x and y are nm
	// ylimit = y + (6 * hleg);
	// xlimit = x - dwleg;
	ylimit = y + 1; // - (5 * hleg);
	xlimit = x + 1;
	
	// xdiff = (xlimit / dwleg);
	// fxdiff = (xdiff + 0.5);
	// fxdiff = (xlimit - (fxdiff * dwleg));
	
	fxdiff = fmod(xlimit, dwleg);
	
	// ydiff = (ylimit / (5 * hleg));
	// fydiff = (ydiff + 0.25);
	// fydiff = (ylimit - (fydiff * (5 * hleg)));
	fydiff = 5 * hleg;
	fydiff = fmod(ylimit, fydiff);
	
	printf("Unit Cell X: %f\n", (xlimit - fxdiff));
	printf("Unit Cell Y: %f\n", (ylimit - fydiff));

	if (cuttype == 0) {
		StandardLoop();
	} else if (cuttype == 1) {
		RectLoop();
	} else {
		HexagonLoop();
	}

	if (ytrans != 0)
		acy /= (ytrans + 1);
		
	printf("Atoms along X Unit Cell: %d\n", acx);
	printf("Atoms along Y Unit Cell: %d\n", acy);
	printf("Atoms in Unit Cell: %d\n", atoms);
	printf("Translations: %d\n", translations);
	printf("Fxdiff: %f\n", fxdiff);
	printf("Fydiff: %f\n", fydiff);
	printf("Distance: %f\n", distance);
	
	return 0;
}
