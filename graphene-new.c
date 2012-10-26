/*
 * graphene-new.c
 * 
 * Tool to take an input of height and width and
 * generate a coordinate list of atoms in a sheet 
 * of graphene with aforementioned dimensions.
 * 
 * Copyright 2012 Trevin Gandhi <sonicxml@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
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
			
			// This only runs once, finds where antidots should be. 
			// TODO: Test to see if removing for loop doesn't have any
			// undesired effects, which it shouldn't have
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
			
			// Start Translations, if applicable
			if (ytrans == 0 && xtrans != 0) {
				for (ix = 1; ix <= xtrans; ix++) {
					if (wunits != 0) {
						transx = wunits + (ix * (xlimit - fxdiff));			
						fprintf(file, "%f	%f	0\n", transx, hunits);
						translations++;
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
	// Add 1 to x and y to allow for hexagons to complete
	// TODO: might want to make this just less than the distance across a hexagon
	ylimit = y + 1;
	xlimit = x + 1;
	
	fxdiff = fmod(xlimit, dwleg);
	
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
	
	return 0;
}
