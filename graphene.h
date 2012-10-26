/*
 * Trevin Gandhi - Header file for the generation of graphene sheets
 */


#ifndef GRAPHENE_H
#define GRAPHENE_H

//
// Graphene Measurements
//

// Vertical leg: Hypotenuse (1.42 Angstroms) / 2 because it's a 30-60-90 triangle
#define hleg		0.71

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
// General Config
//

// x and y are the dimensions of the unit cell
static unsigned int x = 40;
static unsigned int y = 80;
// 1 if above units are nm, 0 if angstroms
#define units		0
// Tells what level of the hexagon it's creating
static signed int pointy = 1;
// Number of times to translate y-values
static int ytrans = 0;
// Number of times to translate x-values
#define xtrans		3
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
//

// TODO: Automatically calculate these values based on given {length, radius}
// {7,3} Hexagon Height (angstroms)
#define hexy		295.1414575
// {7,3} Hexagon Width (angstroms)
#define hexx		284
#endif
