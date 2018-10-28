#include "ofApp.h"

#define DIMENSION 2

unsigned frame_cnt = 0;

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
	domain = std::uniform_real_distribution<double>(visual.MIN_X, visual.MAX_X);

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



	// Initialize the camera closer to our graph
	cam.setTarget(glm::vec3(0.0f,-5.0f,0.0f));
	cam.setDistance(20.0f);
	//ofSetColor(255,255,0);
}

//--------------------------------------------------------------
void ofApp::update()
{
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
	{
		if (cell.fitness > best.fitness)
		{
            best = cell;

			printf("Best: "); 
			visual.printVector(best.pos); printf("\n");
			printf("Fitness: %f\n", visual.evalFitness(best.pos));
		}
	}

  	/* Randomly replace a cell at a new location */
	if(frame_cnt%(visual.REPRO_STEPS*visual.CHEMO_STEPS)==0)
	{
		const float MAXPROB = 1.0;
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


	}

	/*Kills bacteria with low health*/
	if(frame_cnt%visual.CHEMO_STEPS==0)
		visual.eliminatePop();
}

//--------------------------------------------------------------
void ofApp::draw(){

	static int red = 255;
	static int blue = 255;
	static int green = 255;

	frame_cnt++;
	ofBackgroundGradient(ofColor(65,62,50), ofColor(25,22,10));	

	cam.begin();

	mesh.enableColors();
	ofSetColor(255,255,255);
	mesh.drawWireframe();
	mesh.disableColors();

	//ofSetColor(255, 255, 0);
	if(frame_cnt%(visual.REPRO_STEPS*visual.CHEMO_STEPS)==0||frame_cnt==1)
	{
		re_roll:
		red = rand()%255;
		blue = rand()%255;
		green = rand()%255;
		if(red<10 && blue<10 && green<10)
			goto re_roll;
	}

	ofSetColor(red, green, blue);


	for(int i=0;i<visual.population.size();i++)
    {
    	ofDrawSphere(	glm::vec3(visual.population.at(i).pos[0],
    					ofApp::function(&visual.population.at(i).pos[0], DIMENSION),
    					visual.population.at(i).pos[1]), 0.2
    				);
    }

	cam.end();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

	if(key == 'z')
	{
		visual.POP_SIZE -= 1;
		if(visual.POP_SIZE == 0)
			visual.POP_SIZE = 1;

		std::cout << "POP_SIZE = " << visual.POP_SIZE << std::endl;

		visual.initializePopulation();
	}
		

	if(key == 'x')
	{
		visual.POP_SIZE += 1;

		std::cout << "POP_SIZE = " << visual.POP_SIZE << std::endl;
		
		visual.initializePopulation();
	}

	if(key == 'q')
	{
		visual.ATTRACT_D -= 0.01;

		if(visual.ATTRACT_D < 0.01)
			visual.ATTRACT_D = 0.01;

		std::cout << "ATTRACT_D = " << visual.ATTRACT_D << std::endl;
	}
	
	if(key == 'w')
	{
		visual.ATTRACT_D += 0.01;

		std::cout << "ATTRACT_D = " << visual.ATTRACT_D << std::endl;
	}

	if(key == 'a')
	{
		visual.ATTRACT_W -= 0.01;

		if(visual.ATTRACT_W < 0.01)
			visual.ATTRACT_W = 0.01;
		
		std::cout << "ATTRACT_W = " << visual.ATTRACT_W << std::endl;
	}
	
	if(key == 's')
	{
		visual.ATTRACT_W += 0.01;
		
		std::cout << "ATTRACT_W = " << visual.ATTRACT_W << std::endl;
	}

	if(key == 'e')
	{
		visual.REPEL_H -= 0.01;
		
		if(visual.REPEL_H < 0.01)
			visual.REPEL_H = 0.01;

		std::cout << "REPEL_H = " << visual.REPEL_H << std::endl;
	}

	if(key == 'r')
	{
		visual.REPEL_H += 0.01;
	
		std::cout << "REPEL_H = " << visual.REPEL_H << std::endl;
	}

	if(key == 'd')
	{
		visual.REPEL_W -= 0.01;

		if(visual.REPEL_W < 0.01)
		{
			visual.REPEL_W = 0.01;
		}

		std::cout << "REPEL_W = " << visual.REPEL_W << std::endl;
	}

	if(key == 'f')
	{
		visual.REPEL_W += 0.01;

		std::cout << "REPEL_W = " << visual.REPEL_W << std::endl;
	}
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
