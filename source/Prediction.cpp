// My Libraries
#include <libNavierStokes/Prediction.hpp>

using namespace std;
using namespace bitpit;

vector<vector<double>> Prediction::ComputeVirtualVelocity(PabloUniform& grid, double& dt, double& temps, int& Re, KSP& kspU, KSP& kspV, Mat& MatrixU, Mat& MatrixV, Vec& RHSU, Vec& RHSV, Vec& SolutionU, Vec& SolutionV){

	if (temps==dt) {
		m_Velocity[0]=m_u;
		m_Velocity[1]=m_v;
	}

	ComputeRHS(grid, dt, Re, RHSU, m_u, 1);
	//ComputeRHS(grid, dt, Re, RHSV, m_v, 2);

	if (temps == dt) {
		computeMatrixDiamond(grid, dt, Re, MatrixU, 1);
		computeMatrixDiamond(grid, dt, Re, MatrixV, 2);

		KSPSetOperators(kspU, MatrixU, MatrixU);
  		KSPSetFromOptions(kspU);
  		KSPSetUp(kspU);
  		KSPSetOperators(kspV, MatrixV, MatrixV);
  		KSPSetFromOptions(kspV);
  		KSPSetUp(kspV);
	}

	/*KSPSolve(kspU, RHSU, SolutionU);
	KSPSolve(kspV, RHSV, SolutionV);

	PetscScalar *vecArray1, *vecArray2;
	VecGetArray(SolutionU, &vecArray1);
	VecGetArray(SolutionV, &vecArray2);

	for (int i = 0; i < grid.getNumOctants(); ++i) {
		m_u[i] = (double)PetscRealPart(vecArray1[i]);
		m_v[i] = (double)PetscRealPart(vecArray2[i]);
	}

	VecRestoreArray(SolutionU, &vecArray1);
	VecRestoreArray(SolutionV, &vecArray2);
*/
	m_Velocity[0]=m_u;
	m_Velocity[1]=m_v;

	/*cout<<"U "<<m_Velocity[0]<<'\n';
	cout<<"V "<<m_Velocity[1]<<'\n';*/

	return m_Velocity;
}

void Prediction::ComputeRHS(PabloUniform& grid, double& dt, int& Re, Vec& RHS, vector<double> vect, int parameter) {
	darray3 Position, node0, node1, node2, node3;
	int rankPrevious;
	Octant *OctPrevious;
	vector<uint32_t> neigh_n, neigh_f1, neigh_f2;
	vector<bool> isghost_n, isghost_f1, isghost_f2;
	bool NeighInterp;
	double centerX, centerY;
	vector<double> positionX(grid.getNumOctants(),0.), positionY(grid.getNumOctants(),0.);

	for (int i = 0; i < grid.getNumOctants(); ++i) vect[i] = 0;

	for (int i = 0; i < grid.getNumOctants(); ++i){
		bool boundary=grid.getBound(i, 0);
		NeighInterp = 0;

		if (boundary && parameter==1) vect[i] -= 2*grid.getVolume(i);


		// look for position
		positionX[i] = grid.getCenter(i)[0] - m_Velocity[0][i]*dt;
		positionY[i] = grid.getCenter(i)[1] - m_Velocity[1][i]*dt;
		Position[0] = positionX[i];
		Position[1] = positionY[i];

		OctPrevious = grid.getPointOwner(Position);

		if (OctPrevious == NULL){
			rankPrevious = grid.getPointOwnerRank(Position);
			if (rankPrevious==-1){
				/* code */
			}
			else {
				/* code*/
			}
		}
		else {
			centerX = grid.getCenter(OctPrevious)[0];	centerY = grid.getCenter(OctPrevious)[1];
			node0 = grid.getNode(OctPrevious, 0);	double deltaNode0 = sqrt(pow(centerX-node0[0],2) + pow(centerY-node0[1],2));
			node1 = grid.getNode(OctPrevious, 1);	double deltaNode1 = sqrt(pow(centerX-node1[0],2) + pow(centerY-node1[1],2));
			node2 = grid.getNode(OctPrevious, 2);	double deltaNode2 = sqrt(pow(centerX-node2[0],2) + pow(centerY-node2[1],2));
			node3 = grid.getNode(OctPrevious, 3);	double deltaNode3 = sqrt(pow(centerX-node3[0],2) + pow(centerY-node3[1],2));
			double mini1=min(deltaNode0, deltaNode1);	double mini2=min(mini1, deltaNode2);	double mini3=min(mini2, deltaNode3);
			if (mini3 == deltaNode0){
				grid.findNeighbours(OctPrevious, 0, 2, neigh_n, isghost_n);
				grid.findNeighbours(OctPrevious, 2, 1, neigh_f1, isghost_f1);
				grid.findNeighbours(OctPrevious, 0, 1, neigh_f2, isghost_f2);
				if (neigh_n.size()==0) NeighInterp=1;
			}
			else if (mini3 == deltaNode1) {
				grid.findNeighbours(OctPrevious, 1, 2, neigh_n, isghost_n);
				grid.findNeighbours(OctPrevious, 1, 1, neigh_f1, isghost_f1);
				grid.findNeighbours(OctPrevious, 2, 1, neigh_f2, isghost_f2);
				if (neigh_f2.size() == 2){
					neigh_f2[0] = neigh_f2[1];
					isghost_f2[0] = isghost_f2[1];
				}
				if (neigh_n.size()==0) NeighInterp=1;
			}
			else if (mini3 == deltaNode2) {
				grid.findNeighbours(OctPrevious, 2, 2, neigh_n, isghost_n);
				grid.findNeighbours(OctPrevious, 3, 1, neigh_f1, isghost_f1);
				grid.findNeighbours(OctPrevious, 0, 1, neigh_f2, isghost_f2);
				if (neigh_f2.size() == 2){
					neigh_f2[0] = neigh_f2[1];
					isghost_f2[0] = isghost_f2[1];
				}
				if (neigh_n.size()==0) NeighInterp=1;
			}
			else if (mini3 == deltaNode3){
				grid.findNeighbours(OctPrevious, 3, 2, neigh_n, isghost_n);
				grid.findNeighbours(OctPrevious, 3, 1, neigh_f1, isghost_f1);
				grid.findNeighbours(OctPrevious, 1, 1, neigh_f2, isghost_f2);
				if (neigh_f1.size() == 2){
					neigh_f1[0] = neigh_f1[1];
					isghost_f1[0] = isghost_f1[1];
				}
				if (neigh_f2.size() == 2){
					neigh_f2[0] = neigh_f2[1];
					isghost_f2[0] = isghost_f2[1];
				}
				if (neigh_n.size()==0) NeighInterp=1;
			}

		}

		vect[i] += 0;

		// Find global index
		int g_i=grid.getGlobalIdx(i);

		VecSetValues(RHS, 1,&g_i, &vect[i], INSERT_VALUES);
	}

	MPI_Barrier(PETSC_COMM_WORLD);

	VecAssemblyBegin(RHS);
	VecAssemblyEnd(RHS);
}


void Prediction::computeMatrixDiamond(PabloUniform& grid, double& dt, int& Re, Mat& Mat_poisson, int parameter){
	vector<uint32_t> owners(2,0);
	Intersection *inter;
	Octant *OctOwner, *Oct_n, *Oct_f2, *Oct_f1;
	PetscInt g_n0, g_n1, g_n3, g_n2, g_nn;
	PetscScalar value00, value11, value01, value10, value0N, value02, value03, value1N, value12, value13;
	double TaucenterX, TaucenterY;
	bool ifiner;
	darray3 center0, center1, node_n1, node_n2;
	darr3vector nodes;
	uint iface;
	vector<uint32_t> neigh_n, neigh_f2, neigh_f1;
	vector<bool> isghost_n, isghost_f2, isghost_f1;
	int TautbX, TautbY, NormalX, NormalY;
	vector<double> interpt, interpb;

	for (int i = 0; i < grid.getNumIntersections(); ++i) {
		inter=grid.getIntersection(i);
		owners=grid.getOwners(inter);
		iface=grid.getFace(inter);
		bool boundary=grid.getBound(inter);
		bool bound=0, nodeTop=0, face1=0;

		if (!boundary) {
			
			if (!grid.getIsGhost(inter)) {
				// Find Global index
				g_n0=grid.getGlobalIdx(owners[0]);
				g_n1=grid.getGlobalIdx(owners[1]);
				OctOwner=grid.getOctant(owners[1]);

				if (grid.getLevel(owners[0])==grid.getLevel(owners[1])) {
					value00 = 1./Re*dt/grid.getVolume(owners[0]);
					value11 = 1./Re*dt/grid.getVolume(owners[1]);
					value01 = -1./Re*dt/grid.getVolume(owners[0]);
					value10 = -1./Re*dt/grid.getVolume(owners[1]);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n1, &value01, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n0, &value10, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n1, &value11, ADD_VALUES);
				}
			}
			// ghost
			else {
				g_n0=grid.getGlobalIdx(owners[0]);
				g_n1=grid.getGhostGlobalIdx(owners[1]);
				OctOwner=grid.getGhostOctant(owners[1]);

				if (grid.getLevel(OctOwner)==grid.getLevel(owners[0])) {
					value00 = 1./Re*dt/grid.getVolume(owners[0]);
					value01 = -1./Re*dt/grid.getVolume(owners[0]);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n1, &value01, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
				}
			}

			if (grid.getLevel(OctOwner)!=grid.getLevel(owners[0])){
				ifiner=grid.getFiner(inter);
				center0=grid.getCenter(owners[0]);
				center1=grid.getCenter(OctOwner);
				double deltaCenter=sqrt( pow(center0[0]-center1[0],2) + pow(center0[1]-center1[1],2) );
				nodes=grid.getNodes(inter);
				TautbX=1/sqrt(pow(nodes[1][0]-nodes[0][0],2)+pow(nodes[1][1]-nodes[0][1],2))*(nodes[1][0]-nodes[0][0]);
				TautbY=1/sqrt(pow(nodes[1][0]-nodes[0][0],2)+pow(nodes[1][1]-nodes[0][1],2))*(nodes[1][1]-nodes[0][1]);
				TaucenterX=1/sqrt(pow(center1[0]-center0[0],2)+pow(center1[1]-center0[1],2))*(center1[0]-center0[0]);
				TaucenterY=1/sqrt(pow(center1[0]-center0[0],2)+pow(center1[1]-center0[1],2))*(center1[1]-center0[1]);
				
				// if owners[0] is the finest octant
				if (!ifiner) {
					NormalX=grid.getNormal(inter)[0];
					NormalY=grid.getNormal(inter)[1];

					if (iface==0){
						grid.findNeighbours(owners[0],2,1,neigh_f1,isghost_f1);
						grid.findNeighbours(owners[0],3,1,neigh_f2,isghost_f2);
						grid.findNeighbours(owners[0],0,2,neigh_n,isghost_n);
						if (neigh_n.size()==0) {
							grid.findNeighbours(owners[0],2,2,neigh_n,isghost_n);
							nodeTop=1;
							if (neigh_n.size()==0) bound=1;
						}
					}
					else if (iface==1){
						grid.findNeighbours(owners[0],2,1,neigh_f1,isghost_f1);
						grid.findNeighbours(owners[0],3,1,neigh_f2,isghost_f2);
						grid.findNeighbours(owners[0],1,2,neigh_n,isghost_n);
						if (neigh_n.size()==0) {
							grid.findNeighbours(owners[0],3,2,neigh_n,isghost_n);
							nodeTop=1;
							if (neigh_n.size()==0) bound=1;
						}
					}
					else if (iface==2){
						grid.findNeighbours(owners[0],0,1,neigh_f1,isghost_f1);
						grid.findNeighbours(owners[0],1,1,neigh_f2,isghost_f2);
						grid.findNeighbours(owners[0],0,2,neigh_n,isghost_n);
						if (neigh_n.size()==0) {
							grid.findNeighbours(owners[0],1,2,neigh_n,isghost_n);
							nodeTop=1;
							if (neigh_n.size()==0) bound=1;
						}
					}
					else if (iface==3){
						grid.findNeighbours(owners[0],0,1,neigh_f1,isghost_f1);
						grid.findNeighbours(owners[0],1,1,neigh_f2,isghost_f2);
						grid.findNeighbours(owners[0],2,2,neigh_n,isghost_n);
						if (neigh_n.size()==0) {
							grid.findNeighbours(owners[0],3,2,neigh_n,isghost_n);
							nodeTop=1;
							if (neigh_n.size()==0) bound=1;
						}
					}
				}
				// if owners[1] is the finest octant
				else {
					NormalX=-grid.getNormal(inter)[0];
					NormalY=-grid.getNormal(inter)[1];

					if (iface==0){
						node_n1=grid.getNode(owners[0],1);
						node_n2=grid.getNode(owners[0],3);
						if (node_n1[0]==nodes[0][0] && node_n1[1]==nodes[0][1]) {
							grid.findNeighbours(owners[0],1,2,neigh_n,isghost_n);
							grid.findNeighbours(owners[0],2,1,neigh_f2,isghost_f2);
							face1=1;
							if (neigh_f2.size()==2) {
								neigh_f2[0]=neigh_f2[1];
								isghost_f2[0]=isghost_f2[1];
							}
						}
						else if (node_n2[0]==nodes[1][0] && node_n2[1]==nodes[1][1]){
							grid.findNeighbours(owners[0],3,2,neigh_n,isghost_n);
							grid.findNeighbours(owners[0],3,1,neigh_f2,isghost_f2);
							nodeTop=1;
							if (neigh_f2.size()==2) {
								neigh_f2[0]=neigh_f2[1];
								isghost_f2[0]=isghost_f2[1];
							}
						}
						if (neigh_n.size()==0) bound=1;
						grid.findNeighbours(owners[0],1,1,neigh_f1,isghost_f1);
						if (isghost_f1[0]==1) g_n3=grid.getGhostGlobalIdx(neigh_f1[0]);
						else g_n3=grid.getGlobalIdx(neigh_f1[0]);
						if (g_n3==g_n1) {
							neigh_f1[0]=neigh_f1[1];
							isghost_f1[0]=isghost_f1[1];
						}
					}
					else if (iface==1){
						node_n1=grid.getNode(owners[0],0);
						node_n2=grid.getNode(owners[0],2);
						if (node_n1[0]==nodes[0][0] && node_n1[1]==nodes[0][1]) {
							grid.findNeighbours(owners[0],0,2,neigh_n,isghost_n);
							grid.findNeighbours(owners[0],2,1,neigh_f2,isghost_f2);
							face1=1;
						}
						else if (node_n2[0]==nodes[1][0] && node_n2[1]==nodes[1][1]){
							grid.findNeighbours(owners[0],2,2,neigh_n,isghost_n);
							grid.findNeighbours(owners[0],3,1,neigh_f2,isghost_f2);
							nodeTop=1;
						}
						if (neigh_n.size()==0) bound=1;
						grid.findNeighbours(owners[0],0,1,neigh_f1,isghost_f1);
						if (isghost_f1[0]==1) g_n3=grid.getGhostGlobalIdx(neigh_f1[0]);
						else g_n3=grid.getGlobalIdx(neigh_f1[0]);
						if (g_n3==g_n1) {
							neigh_f1[0]=neigh_f1[1];
							isghost_f1[0]=isghost_f1[1];
						}
					}
					else if (iface==2){
						node_n1=grid.getNode(owners[0],2);
						node_n2=grid.getNode(owners[0],3);
						if (node_n1[0]==nodes[0][0] && node_n1[1]==nodes[0][1]) {
							grid.findNeighbours(owners[0],2,2,neigh_n,isghost_n);
							grid.findNeighbours(owners[0],0,1,neigh_f2,isghost_f2);
							face1=1;
							if (neigh_f2.size()==2) {
								neigh_f2[0]=neigh_f2[1];
								isghost_f2[0]=isghost_f2[1];
							}
						}
						else if (node_n2[0]==nodes[1][0] && node_n2[1]==nodes[1][1]){
							grid.findNeighbours(owners[0],3,2,neigh_n,isghost_n);
							grid.findNeighbours(owners[0],1,1,neigh_f2,isghost_f2);
							nodeTop=1;
							if (neigh_f2.size()==2) {
								neigh_f2[0]=neigh_f2[1];
								isghost_f2[0]=isghost_f2[1];
							}
						}
						if (neigh_n.size()==0) bound=1;
						grid.findNeighbours(owners[0],3,1,neigh_f1,isghost_f1);
						if (isghost_f1[0]==1) g_n3=grid.getGhostGlobalIdx(neigh_f1[0]);
						else g_n3=grid.getGlobalIdx(neigh_f1[0]);
						if (g_n3==g_n1) {
							neigh_f1[0]=neigh_f1[1];
							isghost_f1[0]=isghost_f1[1];
						}
					}
					else if (iface==3){
						node_n1=grid.getNode(owners[0],0);
						node_n2=grid.getNode(owners[0],1);
						if (node_n1[0]==nodes[0][0] && node_n1[1]==nodes[0][1]) {
							grid.findNeighbours(owners[0],0,2,neigh_n,isghost_n);
							grid.findNeighbours(owners[0],0,1,neigh_f2,isghost_f2);
							face1=1;
						}
						else if (node_n2[0]==nodes[1][0] && node_n2[1]==nodes[1][1]){
							grid.findNeighbours(owners[0],1,2,neigh_n,isghost_n);
							grid.findNeighbours(owners[0],1,1,neigh_f2,isghost_f2);
							nodeTop=1;
						}
						if (neigh_n.size()==0) bound=1;
						grid.findNeighbours(owners[0],2,1,neigh_f1,isghost_f1);
						if (isghost_f1[0]==1) g_n3=grid.getGhostGlobalIdx(neigh_f1[0]);
						else g_n3=grid.getGlobalIdx(neigh_f1[0]);
						if (g_n3==g_n1) {
							neigh_f1[0]=neigh_f1[1];
							isghost_f1[0]=isghost_f1[1];
						}
					}
				} // end ifiner
				if (!bound) {
					if (isghost_n[0]==1) {
						g_nn=grid.getGhostGlobalIdx(neigh_n[0]);
						Oct_n=grid.getGhostOctant(neigh_n[0]);
					}	
					else {
						g_nn=grid.getGlobalIdx(neigh_n[0]);
						Oct_n=grid.getOctant(neigh_n[0]);
					}

					if (isghost_f2[0]==1) {
						g_n2=grid.getGhostGlobalIdx(neigh_f2[0]);
						Oct_f2=grid.getGhostOctant(neigh_f2[0]);
					}	
					else {
						g_n2=grid.getGlobalIdx(neigh_f2[0]);
						Oct_f2=grid.getOctant(neigh_f2[0]);
					}
					if (isghost_f1[0]==1) {
						g_n3=grid.getGhostGlobalIdx(neigh_f1[0]);
						Oct_f1=grid.getGhostOctant(neigh_f1[0]);
					}	
					else {
						g_n3=grid.getGlobalIdx(neigh_f1[0]);
						Oct_f1=grid.getOctant(neigh_f1[0]);
					}

					if (nodeTop){
						interpt = interp4p(grid, nodes[1], owners[0], OctOwner, Oct_f2, Oct_n);
						interpb = interp3p(grid, nodes[0], owners[0], OctOwner, Oct_f1);
						if (TautbX == 0) {
							value0N = 1./Re*dt/grid.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[3];
							value02 = 1./Re*dt/grid.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
							value03 = -1./Re*dt/grid.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];

							value1N = -1./Re*dt/grid.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[3];
							value12 = -1./Re*dt/grid.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
							value13 = 1./Re*dt/grid.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];
						}
						else {
							value0N = 1./Re*dt/grid.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[3];
							value02 = 1./Re*dt/grid.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
							value03 = -1./Re*dt/grid.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];

							value1N = -1./Re*dt/grid.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[3];
							value12 = -1./Re*dt/grid.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
							value13 = 1./Re*dt/grid.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];
						}
					}
					else if (face1) {
						interpt = interp3p(grid, nodes[1], owners[0], OctOwner, Oct_f1);
						interpb = interp4p(grid, nodes[0], owners[0], OctOwner, Oct_f2, Oct_n);
						if (TautbX == 0) {
							value0N = -1./Re*dt/grid.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[3];
							value02 = -1./Re*dt/grid.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];
							value03 = 1./Re*dt/grid.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];

							value1N = 1./Re*dt/grid.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[3];
							value12 = 1./Re*dt/grid.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];
							value13 = -1./Re*dt/grid.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
						}
						else {
							value0N = -1./Re*dt/grid.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[3];
							value02 = -1./Re*dt/grid.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];
							value03 = 1./Re*dt/grid.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];

							value1N = 1./Re*dt/grid.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[3];
							value12 = 1./Re*dt/grid.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];
							value13 = -1./Re*dt/grid.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
						}
					}
					else {
						interpt = interp3p(grid, nodes[1], owners[0], OctOwner, Oct_f2);
						interpb = interp4p(grid, nodes[0], owners[0], OctOwner, Oct_f1, Oct_n);
						if (TautbX == 0) {
							value0N = -1./Re*dt/grid.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[3];
							value02 = 1./Re*dt/grid.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
							value03 = -1./Re*dt/grid.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];

							value1N = 1./Re*dt/grid.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[3];
							value12 = -1./Re*dt/grid.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
							value13 = 1./Re*dt/grid.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];
						}
						else {
							value0N = -1./Re*dt/grid.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[3];
							value02 = 1./Re*dt/grid.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
							value03 = -1./Re*dt/grid.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];

							value1N = 1./Re*dt/grid.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[3];
							value12 = -1./Re*dt/grid.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
							value13 = 1./Re*dt/grid.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];
						}
					}

					if (grid.getIsGhost(inter)){
						value1N = 0;
						value12 = 0;
						value13 = 0;
					}

					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_nn, &value0N, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n3, &value03, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n2, &value02, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_nn, &value1N, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n3, &value13, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n2, &value12, ADD_VALUES);

					if (!grid.getIsGhost(inter)) {
						if (TautbX ==0 ){
							value00 = - 1./Re*dt/grid.getVolume(owners[0])*NormalX/TaucenterX*(-grid.getArea(inter)/deltaCenter + TaucenterY/TautbY*(interpb[0]-interpt[0]));
							value01 = - 1./Re*dt/grid.getVolume(owners[0])*NormalX/TaucenterX*(grid.getArea(inter)/deltaCenter + TaucenterY/TautbY*(interpb[1]-interpt[1]));
							value10 = - 1./Re*dt/grid.getVolume(OctOwner)*NormalX/TaucenterX*(grid.getArea(inter)/deltaCenter - TaucenterY/TautbY*(interpb[0]-interpt[0]));
							value11 = - 1./Re*dt/grid.getVolume(OctOwner)*NormalX/TaucenterX*(-grid.getArea(inter)/deltaCenter - TaucenterY/TautbY*(interpb[1]-interpt[1]));
						}
						else {
							value00 = - 1./Re*dt/grid.getVolume(owners[0])*NormalY/TaucenterY*(-grid.getArea(inter)/deltaCenter + TaucenterX/TautbX*(interpb[0]-interpt[0]));
							value01 = - 1./Re*dt/grid.getVolume(owners[0])*NormalY/TaucenterY*(grid.getArea(inter)/deltaCenter + TaucenterX/TautbX*(interpb[1]-interpt[1]));
							value10 = - 1./Re*dt/grid.getVolume(OctOwner)*NormalY/TaucenterY*(grid.getArea(inter)/deltaCenter - TaucenterX/TautbX*(interpb[0]-interpt[0]));
							value11 = - 1./Re*dt/grid.getVolume(OctOwner)*NormalY/TaucenterY*(-grid.getArea(inter)/deltaCenter - TaucenterX/TautbX*(interpb[1]-interpt[1]));
						}
					}
					else {
						if (TautbX ==0 ){
							value00 = - 1./Re*dt/grid.getVolume(owners[0])*NormalX/TaucenterX*(-grid.getArea(inter)/deltaCenter + TaucenterY/TautbY*(interpb[0]-interpt[0]));
							value10 = 0;
							value01 = - 1./Re*dt/grid.getVolume(owners[0])*NormalX/TaucenterX*(grid.getArea(inter)/deltaCenter + TaucenterY/TautbY*(interpb[1]-interpt[1]));
							value11= 0;
						}
						else {
							value00 = - 1./Re*dt/grid.getVolume(owners[0])*NormalY/TaucenterY*(-grid.getArea(inter)/deltaCenter + TaucenterX/TautbX*(interpb[0]-interpt[0]));
							value10 = 0;
							value01 = - 1./Re*dt/grid.getVolume(owners[0])*NormalY/TaucenterY*(grid.getArea(inter)/deltaCenter + TaucenterX/TautbX*(interpb[1]-interpt[1]));
							value11 = 0;
						}
					}
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n1, &value01, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n0, &value10, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n1, &value11, ADD_VALUES);

				} // end of !bound
				else{
					if (!grid.getIsGhost(inter)){
						value00 = 1./Re*dt/grid.getVolume(owners[0])*grid.getArea(inter)/deltaCenter;
						value11 = 1./Re*dt/grid.getVolume(OctOwner)*grid.getArea(inter)/deltaCenter;
						value01 = -1./Re*dt/grid.getVolume(owners[0])*grid.getArea(inter)/deltaCenter;
						value10 = -1./Re*dt/grid.getVolume(OctOwner)*grid.getArea(inter)/deltaCenter;
					}
					else {
						value00 = 1./Re*dt/grid.getVolume(owners[0])*grid.getArea(inter)/deltaCenter;
						value11 = 0;
						value01 = -1./Re*dt/grid.getVolume(owners[0])*grid.getArea(inter)/deltaCenter;
						value10 = 0;
					}
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n1, &value01, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n0, &value10, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n1, &value11, ADD_VALUES);
				}
			} // end of level difference
		} // end of !boundary
		else {
			if (parameter == 1) {
				if (iface == 0) value00 = 2*grid.getVolume(owners[0]);
			}
			else if (parameter == 2) {
				if (iface == 0 || iface == 2 || iface == 3) value00 = 2*grid.getVolume(owners[0]);
			}
			else value00 = 0;
			MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
		} // end of boundary
	} // end of intersection loop

	for (int j = 0; j < grid.getNumOctants(); ++j){
		if ( sqrt(pow(grid.getCenter(j)[0],2)+ pow(grid.getCenter(j)[1],2)) <= m_radius )
			m_pen=1;
		else
			m_pen =0;

		// Find Global index
		int g_j=grid.getGlobalIdx(j);
		value00= 1+m_pen*m_lambda*dt;

		MatSetValues(Mat_poisson, 1, &g_j, 1, &g_j, &value00, ADD_VALUES);
	}

	MatAssemblyBegin(Mat_poisson,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Mat_poisson,MAT_FINAL_ASSEMBLY);
	//MatView(Mat_poisson, PETSC_VIEWER_STDOUT_WORLD);
};

vector<double> Prediction::interp4p(PabloUniform& grid, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2, Octant*& Oct3){
	vector<double> value(4,0);
	double x0, y0, x1, y1, x2, y2, x3, y3;
	double intermed1, xdif, ydif, xydif, intermed2, intermed3, intermed4;
	double long1, long2, long3, long4, long5, long6, long7, long8, long9, long10, long11, long12, long13, long14;

	x0=grid.getCenter(owner)[0]; y0=grid.getCenter(owner)[1];
	x1=grid.getCenter(Oct1)[0]; y1=grid.getCenter(Oct1)[1];
	x2=grid.getCenter(Oct2)[0]; y2=grid.getCenter(Oct2)[1];
	x3=grid.getCenter(Oct3)[0]; y3=grid.getCenter(Oct3)[1];

	xdif = x0-x1; 	ydif = y0-y1; 	xydif = x0*y0-x1*y1;
	intermed1 = (x1-x2)/(x0-x1); 	intermed2 = intermed1*(y0-y1)-(y1-y2);
	intermed3 = x1*y1-x2*y2-(x1-x2)/(x0-x1)*xydif; 	intermed4 = (x2-x3)/intermed3*xydif/xdif-(x2*y2-x3*y3)/intermed3;

	long1 = ydif/xdif*(x2-x3) + (x2-x3)/intermed3*intermed2*xydif/xdif - (y2-y3) - (x2*y2-x3*y3)/intermed3*intermed2;
	long2 = (x2-x3)/xdif - (x2*y2-x3*y3)/intermed3*intermed1 + (x2-x3)/intermed3*intermed1*xydif/xdif;
	long3 = (x2*y2-x3*y3)/intermed3*(1+intermed1) - (x2-x3)/xdif - (x2-x3)/intermed3*(1+intermed1)*xydif/xdif;
	long4 = intermed2/long1*long2 - intermed1;	long5 = intermed2/long1*long3 + 1+intermed1;	long6 = intermed2/long1*(intermed4-1) - 1;
	long7 = 1 - ydif/long1*long2 - xydif/intermed3*long4;	long8 = -1 - ydif/long1*long3 - xydif/intermed3*long5;
	long9 = -ydif/long1*(intermed4-1) - xydif/intermed3*long6;	long10 = -ydif/long1 - xydif/intermed3*intermed2/long1;
	long11 = 1 - x0/xdif*long7 - y0/long1*long2 - x0*y0/intermed3*long4;	long12 = -x0/xdif*long8 - y0/long1*long3 - x0*y0/intermed3*long5;
	long13 = -x0/xdif*long9 - y0/long1*(intermed4-1) - x0*y0/intermed3*long6;	long14 = -x0/xdif*long10 - y0/long1 - x0*y0/intermed3*intermed2/long1;

	value[0] = long11 + nodes[0]/xdif*long7 + nodes[1]/long1*long2 + nodes[0]*nodes[1]/intermed3*long4;
	value[1] = long12 + nodes[0]/xdif*long8 + nodes[1]/long1*long3 + nodes[0]*nodes[1]/intermed3*long5;
	value[2] = long13 + nodes[0]/xdif*long9 + nodes[1]/long1*(intermed4-1) + nodes[0]*nodes[1]/intermed3*long6;
	value[3] = long14 + nodes[0]/xdif*long10 + nodes[1]/long1 + nodes[0]*nodes[1]/intermed3*intermed2/long1;

	return value;
};

vector<double> Prediction::interp3p(PabloUniform& grid, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2){
	vector<double> value(3,0);
	double x0, y0, x1, y1, x2, y2;
	double intermed1, xdif, ydif, intermed2;

	x0=grid.getCenter(owner)[0]; y0=grid.getCenter(owner)[1];
	x1=grid.getCenter(Oct1)[0]; y1=grid.getCenter(Oct1)[1];
	x2=grid.getCenter(Oct2)[0]; y2=grid.getCenter(Oct2)[1];

	xdif=x0-x1; ydif=y0-y1;
	intermed1=(x1-x2)/(x0-x1); intermed2=intermed1*(y0-y1)-(y1-y2);

	value[0] = 1 - x0/xdif + (x0*ydif)/(xdif*intermed2)*intermed1 - y0/intermed2*intermed1 + nodes[0]/xdif*(1-intermed1/intermed2*ydif) + nodes[1]/intermed2*intermed1;
	value[1] = x0/xdif + (x0*ydif)/(xdif*intermed2)*(-1-intermed1) + y0/intermed2*(1+intermed1) + nodes[0]/xdif*((1+intermed1)/intermed2*ydif-1) + nodes[1]/intermed2*(-1-intermed1);
	value[2] = (x0*ydif)/(xdif*intermed2) - y0/intermed2 - nodes[0]/xdif*ydif/intermed2 + nodes[1]/intermed2;

	return value;
};