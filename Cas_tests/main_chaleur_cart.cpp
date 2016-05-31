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
#include <libNavierStokes/HeatEq.hpp>
#include <libNavierStokes/Export.hpp>

using namespace std;
using namespace bitpit;

int main(int argc, char** argv)
{
	/*------ Initialization of the Petsc communicator ------*/
 	PetscInitialize(&argc,&argv,(char *)0,help);
 	// Variables
	PetscMPIInt size,rank;
  KSP ksp;
  Mat Mat_Poisson;
  Vec Temp, TempNext;

	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	/*------ Building of the grid ------*/
  double xBegin=-10.0, yBegin=-10.0, lengthDomain=20.0;
  int nbRefine=6;
	PabloUniform gridCart(xBegin,yBegin,0.0,lengthDomain);
  for (int i = 0; i < nbRefine; ++i){
    gridCart.adaptGlobalRefine();
  }

  //gridCart.setMarker(1,1);

  for (int i = 0; i < gridCart.getNumOctants(); ++i) {
    if (sqrt(pow(gridCart.getCenter(i)[0],2)+pow(gridCart.getCenter(i)[1],2))<5 && sqrt(pow(gridCart.getCenter(i)[0],2)+pow(gridCart.getCenter(i)[1],2))>2) {
      gridCart.setMarker(i,2);
    }
    if (sqrt(pow(gridCart.getCenter(i)[0],2)+pow(gridCart.getCenter(i)[1],2))<10 && sqrt(pow(gridCart.getCenter(i)[0],2)+pow(gridCart.getCenter(i)[1],2))>=5)
      gridCart.setMarker(i,1);
  }
  gridCart.adapt();

  gridCart.computeConnectivity();

  gridCart.loadBalance();

  gridCart.updateConnectivity();
  gridCart.computeIntersections();


  /*------ Create Data ------*/
  VecCreate(PETSC_COMM_WORLD, &Temp);
  VecSetType(Temp, VECSTANDARD);

  MatCreate(PETSC_COMM_WORLD, &Mat_Poisson);
  MatSetType(Mat_Poisson , MATMPIAIJ);

  VecSetSizes(Temp,gridCart.getNumOctants(),PETSC_DECIDE);
  VecDuplicate(Temp, &TempNext);
  MatSetSizes(Mat_Poisson,gridCart.getNumOctants(),gridCart.getNumOctants(),PETSC_DETERMINE,PETSC_DETERMINE);
  MatMPIAIJSetPreallocation(Mat_Poisson, 12, PETSC_NULL, 12, PETSC_NULL);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp,KSPBCGS);
  KSPSetTolerances(ksp,1.e-9,PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

  /*------ Computing heat equation ------*/
  HeatEq heatEq(gridCart);
  double dt = 1, tmax=150, temps=dt;
  int iter = 0;
  vector<double> T;

  do {
    T=heatEq.computeTNext(gridCart,dt,temps, ksp, Mat_Poisson, Temp, TempNext);
    temps+=dt;
    iter+=1;
    gridCart.writeTest("grid_iter"+to_string(static_cast<unsigned long long>(iter)),T);

  } while (temps<=tmax);


  VecDestroy(&Temp);
  VecDestroy(&TempNext);
  MatDestroy(&Mat_Poisson);
  KSPDestroy(&ksp);
	PetscFinalize();

  return 0;
}