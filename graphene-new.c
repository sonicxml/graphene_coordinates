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
// Tells what level of the hexagon it's creating
static signed int pointy = 1, holey = 1;
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

// Still a WIP
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

// Creates rectangular (including square) antidots
int RectLoop () {
	float transx, transy, iy = 1, ix = 1;
	unsigned int oppx, oppy, checker = 0;
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
	
	for (hunits = 0; hunits <= ylimit; hunits += (armchair) ? wleg : (hunits == 0) ? hleg : hleg * (abs(pointy) > 1 ? 2 : 1)) {
		for (wunits = ((pointy > 0) ? (armchair) ? hleg : wleg : 0); wunits <= xlimit; wunits += (armchair) ? (holey > 0) ? hleg : (dhleg - hleg) : dwleg) {
			
			// Check to see if there should be an antidot at that coordinate
			if ((hunits >= recty && hunits <= oppy) 
				&& (wunits >= rectx && wunits <= oppx)) {
				cut = true;
			} else {
				cut = false;
			}
			
			if (!cut) {
				fprintf(file, "%f	%f	0\n", wunits, hunits);
				atoms++;
				for (ix = 1; ix <= xtrans; ix++) {
						transx = wunits + (ix * (xlimit - fxdiff));
						fprintf(file, "%f	%f	0\n", transx, hunits);
						for (iy = 1; iy <= ytrans; iy++) {
							transy = hunits + (iy * (ylimit - fydiff));
							fprintf(file, "%f	%f	0\n", wunits, transy);
							fprintf(file, "%f	%f	0\n", transx, transy);
						}
				}
			}
			
			if (hunits == 0)
				acx++;
				
			if (armchair) {
				if (checker == 0 && pointy == 1) {
					holey = 1;
					checker = 1;
				} else if (checker == 0 && pointy == -1) {
					holey = -1;
					checker = 1;
				} else {
					holey += sign(holey);
					if ((abs(holey)) > 1) 
						holey = -sign(holey);
				}				
			}
		}
		
		acy++;
		if (hunits == 0) {
			pointy = -1;
		} else {
			pointy += sign(pointy);
			if ((abs(pointy)) > ((armchair)?1:2)) 
				pointy = -sign(pointy);
		}
	
	checker = 0;

	}

	
	// Close the file
	fclose(file);

	return 0;
}

// Creates a standard graphene sheet without antidots
int StandardLoop () {

	float transx, transy, iy = 1, ix = 1;
	int checker = 0;
	// Open File
	FILE *file;
	file = fopen("graphenecoordinates.txt", "w");

	for (hunits = 0; hunits <= ylimit; hunits += (armchair) ? wleg : (hunits == 0) ? hleg : hleg * (abs(pointy) > 1 ? 2 : 1)) {
		for (wunits = ((pointy > 0) ? (armchair) ? hleg : wleg : 0); wunits <= xlimit; wunits += (armchair) ? (holey > 0) ? hleg : (dhleg - hleg) : dwleg) {
			fprintf(file, "%f	%f	0\n", wunits, hunits);
			atoms++;
			
			// Start Translations, if applicable
			for (ix = 1; ix <= xtrans; ix++) {
				transx = wunits + (ix * (xlimit - fxdiff));
				fprintf(file, "%f	%f	0\n", transx, hunits);
				for (iy = 1; iy <= ytrans; iy++) {
					transy = hunits + (iy * (ylimit - fydiff));
					fprintf(file, "%f	%f	0\n", wunits, transy);
					fprintf(file, "%f	%f	0\n", transx, transy);
				}
			}
			
			if (hunits == 0)
				acx++;
				
			if (armchair) {
				if (checker == 0 && pointy == 1) {
					holey = 1;
					checker = 1;
				} else if (checker == 0 && pointy == -1) {
					holey = -1;
					checker = 1;
				} else {
					holey += sign(holey);
					if ((abs(holey)) > 1) 
						holey = -sign(holey);
				}				
			}
		}
		
		acy++;
		if (hunits == 0) {
			pointy = -1;
		} else {
			pointy += sign(pointy);
			if ((abs(pointy)) > ((armchair)?1:2)) 
				pointy = -sign(pointy);
		}
		
		checker = 0;

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
	ylimit = y;
	xlimit = x; 

	fxdiff = fmod(xlimit, ((armchair)?dhleg:dwleg));
	
	fydiff = ((armchair)?dwleg:(dhleg));
	fydiff = fmod(ylimit, fydiff);
	
	printf("Unit Cell X: %f\n", (xlimit - fxdiff));
	printf("Unit Cell Y: %f\n", (ylimit - fydiff));

	if (cuttype == 2) {
		HexagonLoop();
	} else if (cuttype == 1) {
		RectLoop();
	} else {
		StandardLoop();
	}

	if (ytrans != 0)
		acy /= (ytrans + 1);
		
	printf("Atoms along X Unit Cell: %d\n", acx);
	printf("Atoms along Y Unit Cell: %d\n", acy);
	printf("Atoms in Unit Cell: %d\n", atoms);
	printf("Fxdiff: %f\n", fxdiff);
	printf("Fydiff: %f\n", fydiff);
	
	return 0;
}
