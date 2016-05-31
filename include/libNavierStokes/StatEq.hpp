#ifndef STATEQ_HPP
#define STATEQ_HPP

// Petsc Libraries
#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

#include "bitpit.hpp"
#include <vector>
#include <cmath>

using namespace bitpit;

class StatEq
{
public:
	/* Constructor */
	StatEq(PabloUniform& grid) : m_T(grid.getGlobalNumOctants(),0.0) {};

	void computeMatrix(PabloUniform& grid, Mat& Mat_poisson);
	void computeMatrixDiamond(PabloUniform& grid, Mat& Mat_poisson);
	void computeVectorT(PabloUniform& grid, Vec& RHS);
	vector<double> computeTNext(PabloUniform& grid, KSP& ksp, Mat& Mat_poisson, Vec& RHS, Vec& Solution);
	vector<double> interp4p(PabloUniform& grid, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2, Octant*& Oct3);
	vector<double> interp3p(PabloUniform& grid, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2);
	vector<double> barycenter2p(PabloUniform& grid, darray3 nodes, int owner, Octant*& Oct1);

protected:
	vector<double> m_T;
	int m_pen;
	double m_lambda=1.e16;
};

#endif