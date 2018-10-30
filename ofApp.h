#pragma once

#include "ofMain.h"
#include "ofNode.h"
#include "bacteria.hpp"
#include <vector>
#include <random>
#include <iostream>

// function pointer to an evaluation function
//typedef double (ofApp:: eval_func_t)(double *, unsigned int);

typedef struct {
    std::function<double(double *, unsigned int)> fitnessFunc; // the actual evaluation function
    bool isMin;                                                // whether or not we're trying to minimize or maximize
} probFunc;

class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);

	private:
		ofEasyCam cam;
		ofVboMesh mesh;
		ofSpherePrimitive sphere;
		bacteria visual;
		cell_t best;
		std::default_random_engine gen;
		std::uniform_real_distribution<double> domain;

	private:
		static double function(double * coords, unsigned int dim);
		static double function2(double * coords, unsigned int dim);
		
		const static unsigned int numFitnessFuncs = 2;
		probFunc fitnessFuncs[numFitnessFuncs];
		unsigned int fitnessFuncIndex = 0;
};
