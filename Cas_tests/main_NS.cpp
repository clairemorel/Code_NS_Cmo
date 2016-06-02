static char help[] = "";

// Petsc Libraries
#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "bitpit.hpp"

// My Libraries
#include <libNavierStokes/Prediction.hpp>
#include <libNavierStokes/Export.hpp>

using namespace std;
using namespace bitpit;

int main(int argc, char** argv)
{
	PetscInitialize(&argc,&argv,(char *)0,help);
	PetscMPIInt size,rank;
	KSP kspU, kspV;
	Mat MatrixU, MatrixV;
	Vec U, UVirtual, V, VVirtual;

	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	/*------ Building of the grid ------*/
	double xBegin=-10.0, yBegin=-10.0, lengthDomain=20.0;
	int nbRefine=4;
	PabloUniform grid(xBegin,yBegin,0.0,lengthDomain);
	for (int i = 0; i < nbRefine; ++i){
		grid.adaptGlobalRefine();
	}

	//grid.setMarker(1,1);

	for (int i = 0; i < grid.getNumOctants(); ++i) {
		if (sqrt(pow(grid.getCenter(i)[0],2)+pow(grid.getCenter(i)[1],2))<5 && sqrt(pow(grid.getCenter(i)[0],2)+pow(grid.getCenter(i)[1],2))>2) {
			grid.setMarker(i,2);
		}
		if (sqrt(pow(grid.getCenter(i)[0],2)+pow(grid.getCenter(i)[1],2))<8 && sqrt(pow(grid.getCenter(i)[0],2)+pow(grid.getCenter(i)[1],2))>=5) grid.setMarker(i,1);
	}
	grid.adapt();

	grid.computeConnectivity();

	grid.loadBalance();

	grid.updateConnectivity();
	grid.computeIntersections();

	/*------ Create Data ------*/
	VecCreate(PETSC_COMM_WORLD, &U);
	VecSetType(U, VECSTANDARD);

	MatCreate(PETSC_COMM_WORLD, &MatrixU);
	MatSetType(MatrixU, MATMPIAIJ);
	MatCreate(PETSC_COMM_WORLD, &MatrixV);
	MatSetType(MatrixV, MATMPIAIJ);

	VecSetSizes(U,grid.getNumOctants(),PETSC_DECIDE);
	VecDuplicate(U, &UVirtual);
	VecDuplicate(U, &V);
	VecDuplicate(U, &VVirtual);
	MatSetSizes(MatrixU,grid.getNumOctants(),grid.getNumOctants(),PETSC_DETERMINE,PETSC_DETERMINE);
	MatMPIAIJSetPreallocation(MatrixU, 12, PETSC_NULL, 12, PETSC_NULL);
	MatSetSizes(MatrixV,grid.getNumOctants(),grid.getNumOctants(),PETSC_DETERMINE,PETSC_DETERMINE);
	MatMPIAIJSetPreallocation(MatrixV, 12, PETSC_NULL, 12, PETSC_NULL);

	KSPCreate(PETSC_COMM_WORLD, &kspU);
	KSPSetType(kspU,KSPBCGS);
	KSPSetTolerances(kspU,1.e-9,PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPCreate(PETSC_COMM_WORLD, &kspV);
	KSPSetType(kspV,KSPBCGS);
	KSPSetTolerances(kspV,1.e-9,PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

	/*------ Computing Prediction ------*/
	Prediction prediction(grid);
	vector<vector<double>> Velocity(2, vector<double>(grid.getGlobalNumOctants(), 0.));
	int Re=200;
	double dt=1.0, temps=dt;

	Velocity=prediction.ComputeVirtualVelocity(grid, dt, temps, Re, kspU, kspV, MatrixU, MatrixV, U, V, UVirtual, VVirtual);
	grid.writeTest("U", Velocity[0]);
	grid.writeTest("V", Velocity[1]);
	
	cout<<" Test 2 "<<'\n';

	VecDestroy(&U);
	VecDestroy(&UVirtual);
	VecDestroy(&V);
	VecDestroy(&VVirtual);
	MatDestroy(&MatrixU);
	MatDestroy(&MatrixV);
	KSPDestroy(&kspU);
	KSPDestroy(&kspV);
	PetscFinalize();
	return 0;
}