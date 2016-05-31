#ifndef HEATEQ_HPP
#define HEATEQ_HPP

// Petsc Libraries
#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

#include "bitpit.hpp"
#include <vector>
#include <cmath>

using namespace bitpit;


class HeatEq
{
public:
	/* Constructor */
	HeatEq(PabloUniform& grid) : m_T(grid.getGlobalNumOctants(),0.0), m_TBarre(50.0), m_TNext(grid.getGlobalNumOctants(),0.0) {};

	void computeMatrix(PabloUniform& grid, double& dt, Mat& Mat_poisson);
	void computeMatrixDiamond(PabloUniform& grid, double& dt, Mat& Mat_poisson);
	void computeVectorT(PabloUniform& grid, double& dt, Vec& RHS);
	vector<double> computeTNext(PabloUniform& grid, double& dt, double& temps, KSP& ksp, Mat& Mat_poisson, Vec& RHS, Vec& Solution);
	vector<double> interp4p(PabloUniform& grid, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2, Octant*& Oct3);
	vector<double> interp3p(PabloUniform& grid, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2);
	vector<int> getNeighbours(PabloUniform& grid, int i);

protected:
	vector<double> m_T, m_RHS;
	double m_TBarre;
	vector<double> m_TNext;

	double m_radius=1.0;
	int m_pen;
	double m_lambda=1.e8;

};

#endif