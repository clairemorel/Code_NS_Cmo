#ifndef PREDICTION_HPP
#define PREDICTION_HPP

// Petsc Libraries
#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

#include "bitpit.hpp"
#include <vector>
#include <cmath>

using namespace bitpit;

class Prediction
{
public:
	/*Constructor*/
	Prediction(PabloUniform& grid) : m_Velocity(2, vector<double>(grid.getGlobalNumOctants(),0.0)), m_u(grid.getGlobalNumOctants(),1.0), m_v(grid.getGlobalNumOctants(),0.0) {};

	vector<vector<double>> ComputeVirtualVelocity(PabloUniform& grid, double& dt, double& temps, int& Re, KSP& kspU, KSP& kspV, Mat& MatrixU, Mat& MatrixV, Vec& RHSU, Vec& RHSV, Vec& SolutionU, Vec& SolutionV);
	void ComputeRHS(PabloUniform& grid, double& dt, int& Re, Vec& RHS, vector<double> vect, int parameter);
	void computeMatrixDiamond(PabloUniform& grid, double& dt, int& Re, Mat& Mat_poisson, int parameter);
	vector<double> interp4p(PabloUniform& grid, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2, Octant*& Oct3);
	vector<double> interp3p(PabloUniform& grid, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2);

protected:
	vector<vector<double>> m_Velocity;
	vector<double> m_v;
	vector<double> m_u;

	double m_radius=1.0;
	int m_pen;
	double m_lambda=1.e8;
};


#endif