#pragma once

#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>

const int MIN_X = -8.0;
const int MAX_X = 8.0;

class cell_t
{
public:
	std::vector<double> pos;
	double	          fitness;
	double	          health;
};



/* The cell interaction equation `g()` as described by the
 * clever algorithms textbook */


class BacteriaOptimization
{
private:
	const int	POP_SIZE = 20;       // population size - book uses 50
	const double STEP_SIZE = 0.1;	 // Same as book
	const int	ELDISP_STEPS = 5;    // elimination/dispersal events
	const int	REPRO_STEPS = 4;     // reproduction steps
	const int	CHEMO_STEPS = 350;    // chemotaxis/swim events - set as 7*ELDISP_STEPS
	const int	SWIM_LEN = 4;        // how long to swim?
	const double ELIM_PROB = 0.25;	// Probability of elimination
	const double ATTRACT_D = 0.1;	   // attraction coefficient
	const double ATTRACT_W = 0.2;	 // attraction weight?
	const double REPEL_H = 0.1;	   // repel coefficient
	const double REPEL_W = 10.0;	  // repel weight 
   
	const int n = 2;

	unsigned int updateCount = 0;

	unsigned int fitnessCalls;
 
	std::vector<cell_t> population;

	cell_t best; // best cell;

public:

	BacteriaOptimization()
	{
		best.fitness = -9999;

		for (int i = 0; i < POP_SIZE; i++)
		{
			population.push_back(cell_t());
			population.at(i).pos = genRandSol(n);
			population.at(i).fitness = 0.0;
			population.at(i).health = 0.0;
		}
	}


	std::vector<cell_t> & getPopulation()
	{
		return population;
	}

	cell_t & getBest()
	{
		return best;
	}

	unsigned int getFitnessCalls()
	{
		return fitnessCalls;	
	}

	void eliminatePopulation()
	{
		/* Sort by health 
		 * cells now sorted greatest health -> least health */
		std::sort(population.begin(), population.end(),
				[](cell_t a, cell_t b) {
					return a.health > b.health;
				}
		);

		/* Replace the bottom half (worse) of the cells
		 * with their top half (better) counterparts 
		 * 
		 * TODO - this could be played with/adjusted
		 * to calculate based on the total number of cells
		 * we want to replace !!! */
		for (int i = 0; i < population.size() / 2; i++)
		{
			population.at(i + (population.size()/2)) = population.at(i);
		}
	}

	void chemotaxisAndSwim(int n)
	{
		for( int cellNum = 0; cellNum < population.size(); cellNum++ )
		{
		// ------------------- TUMBLING ------------------------------
		// -----------------------------------------------------------

		//printf("cell num %d\n", cellNum);

		// calculate the current cell's fitness
		population.at(cellNum).fitness = evalFitness(population.at(cellNum).pos) + cellInteraction(population.at(cellNum));

		for (int stepNum = 0; stepNum < CHEMO_STEPS; stepNum++)
		{
			/* Create a temp cell and have it take a random step */
			cell_t tempCell;
			tempCell.health = population.at(cellNum).health;
			tempCell.pos.reserve(n);
			std::vector<double> dir = getRandDir(n);

			cell_t curCell = population.at(cellNum);
			for (int i = 0; i < n; i++)
			{
				//printf("Before tumble\n");
				//printf("Dir: "); printVector(dir); printf("\n");
				//printf("curcell pos: "); printVector(curCell.pos); printf("\n");
				//printf("tempcell pos size: %d\n", tempCell.pos.size());
				tempCell.pos.push_back( curCell.pos.at(i) + STEP_SIZE * dir.at(i) );
				//printf("After tumble\n");

				if (tempCell.pos.at(i) > MAX_X) tempCell.pos.at(i) = MAX_X;
				if (tempCell.pos.at(i) < MIN_X) tempCell.pos.at(i) = MIN_X;

			}

			tempCell.fitness = evalFitness(tempCell.pos) + cellInteraction(tempCell);
			/* Exit if we didn't find a better solution? 
			 * because we're MAXIMIZING a problem less is worse*/
			if (tempCell.fitness < population.at(cellNum).fitness) {
				stepNum = CHEMO_STEPS;
			}
			else {
				/* Otherwise the cell = the new cell, and add to the
				* overall health of the cell */
				population.at(cellNum) = tempCell;
				population.at(cellNum).health += tempCell.fitness;
			}
		}


		// -------------------- END TUMBLING ------------------------
		// ---------------------------------------------------------- 
		}
	}

	double cellInteraction(cell_t cell)
	{
		double sumAttract = 0.0;
		double sumRepel = 0.0;
		
		for (cell_t neighborCell : population)
		{
			double diff = 0.0;

			for (int i = 0; i < neighborCell.pos.size(); i++)
			{
				diff += pow((cell.pos.at(i) - neighborCell.pos.at(i)), 2.0);
			}

			sumAttract += -1.0 * ATTRACT_D * exp(-1.0 * ATTRACT_W * diff);
			sumRepel += REPEL_H * exp(-1.0 * REPEL_W * diff);
		}

		return sumAttract + sumRepel;
	}

	std::vector<double> getRandDir(int n)
	{
		std::vector<double> d(n);
		
		/* THESE CAN BE PLAYED WITH!!! */
		double MIN = -1.0;
		double MAX = 1.0;

		/* Since it's just a direction we want between 0 and 1,
		 * and the values will be scaled accordingly 
		 * 
		 * Right now step size is only 0.1 though, so might need to make
		 * higher... */
		for (int i = 0; i < n; i++)
		{
			double num = MIN + (double)rand() / ((double)RAND_MAX / (MAX-MIN));
			d.at(i) = num;
		}

		//printf("Random direction:\n");
		//printVector(d);
		//printf("\n");

		return d;
	}

	std::vector<double> genRandSol(int n)
	{
		std::vector<double> v(n);

		for (int i = 0; i < n; i++)
		{
			// start at -7 to 7 instead of -8 to 8
			v.at(i) = (MIN_X+1) + rand() % ((MAX_X-1) - (MIN_X+1) + 1);
		}

		return v;
	}

	void update()
	{
		//if(updateCount % CHEMO_STEPS)
		if(true)
		{
			/* Swim about */
			chemotaxisAndSwim(n);

			/* Check for a new best */
			for (cell_t cell : population)
			{
				/* -9999 for the initial cell */
				//if (best.fitness == -9999 || cell.fitness >= best.fitness)
				if (cell.fitness > best.fitness)
				{
					best = cell;
					printf("New Best: "); printVector(best.pos); printf("\n");
					printf("Fitness: %f\n", evalFitness(best.pos));
					std::cout << fitnessCalls << std::endl;
				}
			}
		}

		//if(updateCount % CHEMO_STEPS REPO_STEPS == 0)
		if(updateCount % CHEMO_STEPS == 0)
		{
		    eliminatePopulation();
		}

		//if(updateCount % ELDISP_STEPS == 0)
		if(updateCount % CHEMO_STEPS * REPRO_STEPS == 0)
		{
			/* Randomly replace a cell at a new location */
			const double MAXPROB = 1.0;
			for (int cellNum = 0; cellNum < population.size(); cellNum++)
			{
				double num = (double)rand() / ((double)RAND_MAX / (MAXPROB));
				if (num < ELIM_PROB) {
					population.at(cellNum).pos = genRandSol(n);
					population.at(cellNum).health = 0.0;
					population.at(cellNum).fitness = evalFitness(population.at(cellNum).pos);
				}
			}
		}

		++updateCount;
	}

	double evalFitness(const std::vector<double> & v)
	{
		++fitnessCalls;

		double sr = 0.0; // square root
		double sn = 0.0; // sin
		for (int i = 0; i < v.size(); i++)
		{
			/*****************************************/
			/**		     Square Root             **/
			/*****************************************/
			// x sub i
			double a = v.at(i);

			// (-1)^i (i % 4)
			// i+1 because formula states from i=1..N,
			// here i starts at 0
			a += pow(-1, i+1) * ((i+1) % 4);

			// sum to our current total, squared
			sr += pow(a, 2);

			/*****************************************/
			/**		         Sin                 **/
			/*****************************************/
			// x sub i
			a = v.at(i);

			// add to running total (x sub i) ^ i
			sn += pow(a, i+1);
		}

		// square root
		sr = sqrt(sr);
		// the square root is negated
		sr *= -1.0;

		// final value is sqrt(sum1) + sin(sum2)
		return sr + sin(sn);
	}

	void printVector(std::vector<double> v)
	{
		for (double x : v)
		{
			printf("%f, ", x);
		}
	}
};
