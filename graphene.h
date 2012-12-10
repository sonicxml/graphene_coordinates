/*
 * 
 * Header file for the generation of graphene sheets
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

#ifndef GRAPHENE_H
#define GRAPHENE_H

//
// Graphene Measurements (in Zigzag orientation)
//

// Vertical leg: Hypotenuse (1.42 Angstroms) / 2 because it's a 30-60-90 triangle
#define hleg		0.71

// hleg * 4 (height of hexagon)
#define dhleg		2.84

// Horizontal leg: hleg * sqrt(3) (for the other leg of the triangle)
#define wleg 		1.229756073

// wleg * 2  (distance across hexagon)
#define dwleg		2.459512146


//
// Functions
//

// sign(x) function - returns -1 if x < 0, 0 if x == 0, and 1 if x > 0
#define sign(x) (( x > 0 ) - ( x < 0 ))


//
// General Configuration
//

// x and y are the dimensions of the unit cell
static unsigned int x = 4;
static unsigned int y = 4;
// 1 if above units are nm, 0 if angstroms
#define units		0
// Number of times to translate y-values
static int ytrans = 0;
// Number of times to translate x-values
#define xtrans		0
// Graphene orientation
// 0 = zigzag, 1 = armchair
static unsigned int armchair = 1;
// Which kind of cut to make
// 0 = none, 1 = rectangle, 2 = hexagon
#define cuttype		0


//
// Rectangle Cuts
//

// X and Y coordinates (nanometers) of bottom left corner of rectangle
static unsigned int rectx = 1;
static unsigned int recty = 1;
// Height of rectangle
static unsigned int recth = 2;
// Width of rectangle
static unsigned int rectw = 2;			


//
// Hexagon cuts
// N.B. - STILL A WIP
//

// TODO: Automatically calculate these values based on given {length, radius}
// {7,3} Hexagon Height (angstroms)
#define hexy		295.1414575
// {7,3} Hexagon Width (angstroms)
#define hexx		284
#endif
