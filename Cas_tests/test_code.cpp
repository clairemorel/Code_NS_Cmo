static char help[] = "";

// Petsc Libraries
#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "bitpit.hpp"

using namespace std;
using namespace bitpit;

int main(int argc, char** argv)
{
	PetscInitialize(&argc,&argv,(char *)0,help);
 	// Variables
	PetscMPIInt size,rank;

	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	/*------ Building of the grid ------*/
  	double xBegin=-10.0, yBegin=-10.0, lengthDomain=20.0;
  	int nbRefine=1;
	PabloUniform gridCart(xBegin,yBegin,0.0,lengthDomain);
  	for (int i = 0; i < nbRefine; ++i){
    	gridCart.adaptGlobalRefine();
  	}

  	gridCart.setMarker(1,1);

  	gridCart.adapt();

  	gridCart.computeConnectivity();

  	gridCart.loadBalance();

  	gridCart.updateConnectivity();
  	gridCart.computeIntersections();

  	/*------ Computing ------*/
  	double dt = 1, temps=dt;
  	// Variables
  	Octant *OctOwner, *localOwner, *Oct_n, *Oct_f2, *Oct_f1;
  	vector<uint32_t> owners(2,0);
  	Intersection *inter;
  	vector<uint> neigh_n,neigh_f1,neigh_f2;
	vector<bool> isghost_n,isghost_f1,isghost_f2;
	uint iface;
	PetscInt g_n2,g_n3,g_nn;

  	for (int i = 0; i < gridCart.getNumIntersections(); ++i){
		inter=gridCart.getIntersection(i);
		owners=gridCart.getOwners(inter);
		bool boundary=gridCart.getBound(inter);
		bool bound=0;
		bool face2=0;

		if (!boundary) {
			if (!gridCart.getIsGhost(inter)) OctOwner=gridCart.getOctant(owners[1]);
			else OctOwner=gridCart.getGhostOctant(owners[1]);

			if (gridCart.getLevel(OctOwner)!=gridCart.getLevel(owners[0])) {
				if (gridCart.getLevel(owners[0])>gridCart.getLevel(OctOwner)) {
					iface=gridCart.getFace(inter);
					localOwner=gridCart.getOctant(owners[0]);
				}
				else {
					iface=gridCart.getFace(inter);
					localOwner=OctOwner;
				}

				if (iface==0){
					gridCart.findNeighbours(localOwner,2,1,neigh_f1,isghost_f1);
					gridCart.findNeighbours(localOwner,3,1,neigh_f2,isghost_f2);
					gridCart.findNeighbours(localOwner,0,2,neigh_n,isghost_n);
					if (neigh_n.size()==0) {
						gridCart.findNeighbours(localOwner,2,2,neigh_n,isghost_n);
						face2=1;
						if (neigh_n.size()==0) bound=1;
					}
				}
				else if (iface==1){
					gridCart.findNeighbours(localOwner,2,1,neigh_f1,isghost_f1);
					gridCart.findNeighbours(localOwner,3,1,neigh_f2,isghost_f2);
					gridCart.findNeighbours(localOwner,1,2,neigh_n,isghost_n);
					if (neigh_n.size()==0) {
						gridCart.findNeighbours(localOwner,3,2,neigh_n,isghost_n);
						face2=1;
						if (neigh_n.size()==0) bound=1;
					}
				}
				else if (iface==2){
					gridCart.findNeighbours(localOwner,0,1,neigh_f1,isghost_f1);
					gridCart.findNeighbours(localOwner,1,1,neigh_f2,isghost_f2);
					gridCart.findNeighbours(localOwner,0,2,neigh_n,isghost_n);
					if (neigh_n.size()==0) {
						gridCart.findNeighbours(localOwner,1,2,neigh_n,isghost_n);
						face2=1;
						if (neigh_n.size()==0) bound=1;
					}
				}
				else if (iface==3){
					gridCart.findNeighbours(localOwner,0,1,neigh_f1,isghost_f1);
					gridCart.findNeighbours(localOwner,1,1,neigh_f2,isghost_f2);
					gridCart.findNeighbours(localOwner,2,2,neigh_n,isghost_n);
					if (neigh_n.size()==0) {
						gridCart.findNeighbours(localOwner,3,2,neigh_n,isghost_n);
						face2=1;
						if (neigh_n.size()==0) bound=1;
					}
				}

				// Là je me mets dans le cas où les intersections ne sont pas collées à une intersection avec le bord
				if (!bound) {
					cout<<"iface "<<iface<<'\n';
					cout<<"neigh_f2 "<<neigh_f2<<'\n';
					cout <<"owners0 "<<gridCart.getGlobalIdx(owners[0])<<'\n';
					cout <<"owners1 "<<gridCart.getGhostGlobalIdx(owners[1])<<'\n';
					cout <<"localOwner "<<gridCart.getGlobalIdx(localOwner)<<'\n';

					// je teste si mes voisins par les noeuds sont ghosts, je voudrais faire pareil avec les voisins par les faces mais c'est alors que je me rends compte que parfois j'ai un vecteur vide qd je calcule en //
					if (isghost_n[0]==1) {
						g_nn=gridCart.getGhostGlobalIdx(neigh_n[0]);
						Oct_n=gridCart.getGhostOctant(neigh_n[0]);
					}	
					else {
						g_nn=gridCart.getGlobalIdx(neigh_n[0]);
						Oct_n=gridCart.getOctant(neigh_n[0]);
					}
				}
			}
		}
	}
	return 0;
}