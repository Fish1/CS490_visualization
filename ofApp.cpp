#include "ofApp.h"

#define MINIMUM -8
#define MAXIMUM 8
#define DIMENSION 2

// Fitness function
double ofApp::function(double * coords, unsigned int dim)
{
	double sum1 = 0;

	double sum2 = 0;

	for(unsigned int index = 1; index <= dim; ++index)
	{
		sum1 += pow((coords[index - 1] + pow(-1.0, index) * (index % 4)), 2);
		
		sum2 += pow(coords[index - 1], index);	
	}

	double result = (-1.0) * sqrt(sum1) + sin(sum2);

	return result;
}

void ofApp::setup()
{
	ofEnableDepthTest();

	// size is from -8 to 8
	const int size = 16;
	// how many vertices per 1 unit
	const int perUnit = 5;
	// square root of the number of vertices
	const int checks = perUnit * size;
	//size of spheres **bacteria
	const float width = 0.5f;
	//random number generator
	domain = std::uniform_real_distribution<double>(MINIMUM,MAXIMUM);

	// Create Verticies
	for(int z = 0; z < checks; ++z)
	{
		// the z position of the current vertex
		double currentZ = ((double)z / (double)perUnit) - ((double)size / 2.0);

		for(int x = 0; x < checks; ++x)
		{
			// the x position of the current vertex
			double currentX = ((double)x / (double)perUnit) - ((double)size / 2.0);

			// pass in these coordinates to the fitness function to get the y position
			double coord [] = {currentX, currentZ};

			// the y position of the current vertex
			double currentY = function(coord, 2);
			
			ofVec3f point(currentX, currentY, currentZ);
			mesh.addVertex(point);
		}
	}

	// Create indices

	for(unsigned int y = 0; y < checks - 1; ++y)
	{
		for(unsigned int x = 0; x < checks; ++x)
		{
			unsigned int current = x + checks * y;
			unsigned int below = x + checks * (y + 1);
			unsigned int left = (x - 1) + checks * y;
			unsigned int belowRight = (x + 1) + checks * (y + 1);

			if(x == 0)
			{
				mesh.addIndex(current);
				mesh.addIndex(below);
				mesh.addIndex(belowRight);	
			}
			else if(x == checks - 1)
			{
				mesh.addIndex(current);
				mesh.addIndex(left);
				mesh.addIndex(below);
			}
			else
			{
				mesh.addIndex(current);
				mesh.addIndex(below);
				mesh.addIndex(belowRight);
				
				mesh.addIndex(current);
				mesh.addIndex(left);
				mesh.addIndex(below);
			}
		}
	}

	//create spheres
	sphere.setRadius(width);


	best.fitness = -9999;

    /* Generate the initial population */
    for (int i = 0; i < visual.POP_SIZE; i++)
    {
        visual.population.at(i).pos = visual.genRandSol(DIMENSION);
        visual.population.at(i).fitness = 0.0;
        visual.population.at(i).health = 0.0;
    }


	// Initialize the camera closer to our graph
	cam.setTarget(glm::vec3(0.0f,-5.0f,0.0f));
	cam.setDistance(20.0f);
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){

	ofBackgroundGradient(ofColor(65,62,50), ofColor(25,22,10));	

	cam.begin();

	mesh.enableColors();
	ofSetColor(255,255,255);
	mesh.drawWireframe();
	mesh.disableColors();

	ofSetColor(255, 255, 0);

	/* Elimination/Dispersal Events */
    /* Swim about */
	visual.chemotaxisAndSwim(	DIMENSION, 
								visual.STEP_SIZE, 
								visual.ELDISP_STEPS, 
								visual.REPRO_STEPS, 
								visual.CHEMO_STEPS, 
								visual.SWIM_LEN, 
								visual.ELIM_PROB, 
								visual.ATTRACT_D, 
								visual.ATTRACT_W, 
								visual.REPEL_H, 
								visual.REPEL_W);

    /* Check for a new best */
    for (cell_t cell : visual.population)
		if (cell.fitness > best.fitness)
            best = cell;

	/* Randomly replace a cell at a new location */
	const double MAXPROB = 1.0;
	for (int cellNum = 0; cellNum < visual.population.size(); cellNum++)
	{
		double num = (double)rand() / ((double)RAND_MAX / (MAXPROB));
		if (num < visual.ELIM_PROB)
		{
			visual.population.at(cellNum).pos = visual.genRandSol(DIMENSION);
			visual.population.at(cellNum).health = 0.0;
			visual.population.at(cellNum).fitness = visual.evalFitness(visual.population.at(cellNum).pos);
		}
	}

	printf("Best: "); 
	visual.printVector(best.pos); printf("\n");
	printf("Fitness: %f\n", visual.evalFitness(best.pos));

	/**TEST CODE FOR DRAWING SPHERES**/
/*	double randX = domain(gen);
	double randZ = domain(gen);
	double random[] = {randX, randZ};
	double randY = function(random, 2);

	sphere.setPosition(randX, randY, randZ);

    sphere.draw();
*/
	cam.end();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
