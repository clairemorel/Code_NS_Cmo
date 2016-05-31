static char help[] = "";

// Petsc Libraries
#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <libNavierStokes/StatEq.hpp>

#include "bitpit.hpp"

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
	//int nbRefine=6;
	//PabloUniform gridCart(xBegin,yBegin,0.0,lengthDomain);

	/*for (int i = 0; i < nbRefine; ++i){
		gridCart.adaptGlobalRefine();
	}

	gridCart.computeConnectivity();

	for (int i = 0; i < gridCart.getNumOctants(); ++i) {

    	if (sqrt(pow(gridCart.getCenter(i)[0],2)+pow(gridCart.getCenter(i)[1],2))<5 && sqrt(pow(gridCart.getCenter(i)[0],2)+pow(gridCart.getCenter(i)[1],2))>2) {
      		gridCart.setMarker(i,2);
    	}
    	if (sqrt(pow(gridCart.getCenter(i)[0],2)+pow(gridCart.getCenter(i)[1],2))<8 && sqrt(pow(gridCart.getCenter(i)[0],2)+pow(gridCart.getCenter(i)[1],2))>=5)
      	gridCart.setMarker(i,1);
	}*/

    PabloUniform gridCart(-0.5,-0.5,0,1);

	/*gridCart.loadBalance();
	gridCart.computeConnectivity();*/

    //PabloUniform gridCart(xBegin,yBegin,0.0,lengthDomain);

    int maxliv=6, nocts;
    for (int i = 0; i < maxliv; ++i) {
    	if (i<3) {
    		gridCart.adaptGlobalRefine();
    	}
    	else if (i==3) {
    		nocts=gridCart.getNumOctants();
    		for (int j = 0; j < nocts; ++j) {
    			if (j%4==2) {
    				gridCart.setMarker(j,1);
    			}
    		}
    		gridCart.adapt();
    		gridCart.updateConnectivity();
    	}
    	else {
    		gridCart.adaptGlobalRefine();
    	}
    	gridCart.adapt();
    	gridCart.updateConnectivity();
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

	/*------ Computing stationnary Heat Equation ------*/
	StatEq statEq(gridCart);
	vector<double> T, Texact(gridCart.getGlobalNumOctants(),0.0), error(gridCart.getGlobalNumOctants(),0.0);
	double normL2=0, normL1=0, normLinf=0, normL2inter2=0, normL1inter=0, normLinfinter=0, normL2inter=0, volumeInter=0, volumeTot=0, normL1inter2=0, normLinfinter2=0;

    T=statEq.computeTNext(gridCart, ksp, Mat_Poisson, Temp, TempNext);
    gridCart.writeTest("heat",T);

	for (int i = 0; i < gridCart.getNumOctants(); ++i){
		//Texact[i] = sin(gridCart.getCenter(i)[0]) + cos(gridCart.getCenter(i)[1]);
		Texact[i] = sin( pow(gridCart.getCenter(i)[0],2) + pow(gridCart.getCenter(i)[1],2) );
		error[i]=abs(T[i]-Texact[i]);
		normL1inter+=error[i]*gridCart.getVolume(i);
		normL2inter+=pow(error[i],2)*gridCart.getVolume(i);
		normLinfinter=max(normLinfinter,error[i]);//*gridCart.getVolume(i));
		volumeInter+=gridCart.getVolume(i);
	}

	gridCart.writeTest("heat_exact",Texact);
	gridCart.writeTest("error",error);

	MPI_Barrier(PETSC_COMM_WORLD);
	MPI_Reduce(&normL1inter, &normL1inter2, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	MPI_Reduce(&normL2inter, &normL2inter2, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
	MPI_Reduce(&normLinfinter, &normLinfinter2, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
	MPI_Reduce(&volumeInter, &volumeTot, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);

	normL2=sqrt(normL2inter2)/volumeTot;
	normL1=normL1inter2/volumeTot;
	normLinf=normLinfinter2;///volumeTot;

	PetscPrintf(PETSC_COMM_WORLD," Norm 1 of error : %g\n",normL1);
	PetscPrintf(PETSC_COMM_WORLD," Norm 2 of error : %g\n",normL2);
	PetscPrintf(PETSC_COMM_WORLD," Norm inf of error : %g\n",normLinf);


	VecDestroy(&Temp);
	VecDestroy(&TempNext);
	MatDestroy(&Mat_Poisson);
	KSPDestroy(&ksp);
	PetscFinalize();

	return 0;
}