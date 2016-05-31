// My Libraries
#include <libNavierStokes/HeatEq.hpp>

using namespace std;
using namespace bitpit;

vector<double> HeatEq::computeTNext(PabloUniform& gridCart, double& dt, double& temps, KSP& ksp, Mat& Mat_poisson, Vec& RHS, Vec& Solution){
	
	computeVectorT(gridCart,dt, RHS);

	if (temps == dt) {
		computeMatrixDiamond(gridCart, dt, Mat_poisson);

		KSPSetOperators(ksp, Mat_poisson, Mat_poisson);
  		KSPSetFromOptions(ksp);
  		KSPSetUp(ksp);
	}

	KSPSolve(ksp, RHS, Solution);

	PetscScalar *vecArray;
	VecGetArray(Solution, &vecArray);

	for (int i = 0; i < gridCart.getNumOctants(); ++i) {
		m_T[i] = (double)PetscRealPart(vecArray[i]);
	}

	VecRestoreArray(Solution, &vecArray);

	return m_T;
};


void HeatEq::computeVectorT(PabloUniform& gridCart, double& dt, Vec& RHS){
	
	for (int i = 0; i < gridCart.getNumOctants(); ++i){
		if ( sqrt(pow(gridCart.getCenter(i)[0],2)+pow(gridCart.getCenter(i)[1],2) )<= m_radius ) m_pen = 1;
		else m_pen =0;

		m_T[i] += m_pen*m_lambda*m_TBarre*dt;

		// Find global index
		int g_i=gridCart.getGlobalIdx(i);

		VecSetValues(RHS, 1,&g_i, &m_T[i], INSERT_VALUES);
	}

	MPI_Barrier(PETSC_COMM_WORLD);

	VecAssemblyBegin(RHS);
	VecAssemblyEnd(RHS);
};


void HeatEq::computeMatrix(PabloUniform& gridCart, double& dt, Mat& Mat_poisson){
	vector<uint32_t> owners(2,0);
	darray3 center0, center1;
	PetscScalar value00, value11, value01, value10;
	PetscInt g_n0, g_n1;
	double kappa=1.e0;
	Intersection *inter;
	Octant *ghostOct;

	for (int i = 0; i < gridCart.getNumIntersections(); ++i) {
		inter=gridCart.getIntersection(i);
		owners=gridCart.getOwners(inter);
		bool boundary=gridCart.getBound(inter);

		if (!boundary) {
			
			// no ghost
			if (!gridCart.getIsGhost(inter)){
				// Find global index
				g_n0=gridCart.getGlobalIdx(owners[0]);
				g_n1=gridCart.getGlobalIdx(owners[1]);

				center0=gridCart.getCenter(owners[0]);
				center1=gridCart.getCenter(owners[1]);
				double deltaCenter=sqrt( pow(center0[0]-center1[0],2) + pow(center0[1]-center1[1],2) );

				value00 = kappa*dt/gridCart.getArea(owners[0])*gridCart.getArea(inter)/deltaCenter;
				value11 = kappa*dt/gridCart.getArea(owners[1])*gridCart.getArea(inter)/deltaCenter;
				value01 = -kappa*dt/gridCart.getArea(owners[0])*gridCart.getArea(inter)/deltaCenter;
				value10 = -kappa*dt/gridCart.getArea(owners[1])*gridCart.getArea(inter)/deltaCenter;
			}
			else {
				g_n0=gridCart.getGlobalIdx(owners[0]);
				g_n1=gridCart.getGhostGlobalIdx(owners[1]);
				ghostOct = gridCart.getGhostOctant(owners[1]);

				center0=gridCart.getCenter(owners[0]);
				center1=gridCart.getCenter(ghostOct);
				double deltaCenter=sqrt( pow(center0[0]-center1[0],2) + pow(center0[1]-center1[1],2) );

				value00 = kappa*dt/gridCart.getArea(owners[0])*gridCart.getArea(inter)/deltaCenter;
				value11 = 0;
				value01 = -kappa*dt/gridCart.getArea(owners[0])*gridCart.getArea(inter)/deltaCenter;
				value10 = 0;
			}

			MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n1, &value01, ADD_VALUES);
			MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n0, &value10, ADD_VALUES);
			MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
			MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n1, &value11, ADD_VALUES);
		}
	}

	for (int j = 0; j < gridCart.getNumOctants(); ++j){
		if ( sqrt(pow(gridCart.getCenter(j)[0],2)+ pow(gridCart.getCenter(j)[1],2)) <= m_radius )
			m_pen=1;
		else
			m_pen =0;

		// Find Global index
		int g_j=gridCart.getGlobalIdx(j);
		value00= 1+m_pen*m_lambda*dt;

		MatSetValues(Mat_poisson, 1, &g_j, 1, &g_j, &value00, ADD_VALUES);
	}

	MatAssemblyBegin(Mat_poisson,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Mat_poisson,MAT_FINAL_ASSEMBLY);
	/*PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE);
	MatView(Mat_poisson, PETSC_VIEWER_STDOUT_WORLD);*/
}

void HeatEq::computeMatrixDiamond(PabloUniform& gridCart, double& dt, Mat& Mat_poisson){
	vector<uint32_t> owners(2,0);
	Intersection *inter;
	Octant *OctOwner, *Oct_n, *Oct_f2, *Oct_f1;
	PetscInt g_n0, g_n1, g_n3, g_n2, g_nn;
	PetscScalar value00, value11, value01, value10, value0N, value02, value03, value1N, value12, value13;
	double kappa=1.e0, TaucenterX, TaucenterY;
	bool ifiner;
	darray3 center0, center1, node_n1, node_n2;
	darr3vector nodes;
	uint iface;
	vector<uint32_t> neigh_n, neigh_f2, neigh_f1;
	vector<bool> isghost_n, isghost_f2, isghost_f1;
	int TautbX, TautbY, NormalX, NormalY;
	vector<double> interpt, interpb;

	for (int i = 0; i < gridCart.getNumIntersections(); ++i) {
		inter=gridCart.getIntersection(i);
		owners=gridCart.getOwners(inter);
		bool boundary=gridCart.getBound(inter);
		bool bound=0, nodeTop=0, face1=0;

		if (!boundary) {
			
			if (!gridCart.getIsGhost(inter)) {
				// Find Global index
				g_n0=gridCart.getGlobalIdx(owners[0]);
				g_n1=gridCart.getGlobalIdx(owners[1]);
				OctOwner=gridCart.getOctant(owners[1]);

				if (gridCart.getLevel(owners[0])==gridCart.getLevel(owners[1])) {
					value00 = kappa*dt/gridCart.getVolume(owners[0]);
					value11 = kappa*dt/gridCart.getVolume(owners[1]);
					value01 = -kappa*dt/gridCart.getVolume(owners[0]);
					value10 = -kappa*dt/gridCart.getVolume(owners[1]);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n1, &value01, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n0, &value10, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n1, &value11, ADD_VALUES);
				}
			}
			// ghost
			else {
				g_n0=gridCart.getGlobalIdx(owners[0]);
				g_n1=gridCart.getGhostGlobalIdx(owners[1]);
				OctOwner=gridCart.getGhostOctant(owners[1]);

				if (gridCart.getLevel(OctOwner)==gridCart.getLevel(owners[0])) {
					value00 = kappa*dt/gridCart.getVolume(owners[0]);
					value01 = -kappa*dt/gridCart.getVolume(owners[0]);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n1, &value01, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
				}
			}

			if (gridCart.getLevel(OctOwner)!=gridCart.getLevel(owners[0])){
				ifiner=gridCart.getFiner(inter);
				center0=gridCart.getCenter(owners[0]);
				center1=gridCart.getCenter(OctOwner);
				double deltaCenter=sqrt( pow(center0[0]-center1[0],2) + pow(center0[1]-center1[1],2) );
				nodes=gridCart.getNodes(inter);
				TautbX=1/sqrt(pow(nodes[1][0]-nodes[0][0],2)+pow(nodes[1][1]-nodes[0][1],2))*(nodes[1][0]-nodes[0][0]);
				TautbY=1/sqrt(pow(nodes[1][0]-nodes[0][0],2)+pow(nodes[1][1]-nodes[0][1],2))*(nodes[1][1]-nodes[0][1]);
				TaucenterX=1/sqrt(pow(center1[0]-center0[0],2)+pow(center1[1]-center0[1],2))*(center1[0]-center0[0]);
				TaucenterY=1/sqrt(pow(center1[0]-center0[0],2)+pow(center1[1]-center0[1],2))*(center1[1]-center0[1]);
				
				// if owners[0] is the finest octant
				if (!ifiner) {
					iface=gridCart.getFace(inter);
					NormalX=gridCart.getNormal(inter)[0];
					NormalY=gridCart.getNormal(inter)[1];

					if (iface==0){
						gridCart.findNeighbours(owners[0],2,1,neigh_f1,isghost_f1);
						gridCart.findNeighbours(owners[0],3,1,neigh_f2,isghost_f2);
						gridCart.findNeighbours(owners[0],0,2,neigh_n,isghost_n);
						if (neigh_n.size()==0) {
							gridCart.findNeighbours(owners[0],2,2,neigh_n,isghost_n);
							nodeTop=1;
							if (neigh_n.size()==0) bound=1;
						}
					}
					else if (iface==1){
						gridCart.findNeighbours(owners[0],2,1,neigh_f1,isghost_f1);
						gridCart.findNeighbours(owners[0],3,1,neigh_f2,isghost_f2);
						gridCart.findNeighbours(owners[0],1,2,neigh_n,isghost_n);
						if (neigh_n.size()==0) {
							gridCart.findNeighbours(owners[0],3,2,neigh_n,isghost_n);
							nodeTop=1;
							if (neigh_n.size()==0) bound=1;
						}
					}
					else if (iface==2){
						gridCart.findNeighbours(owners[0],0,1,neigh_f1,isghost_f1);
						gridCart.findNeighbours(owners[0],1,1,neigh_f2,isghost_f2);
						gridCart.findNeighbours(owners[0],0,2,neigh_n,isghost_n);
						if (neigh_n.size()==0) {
							gridCart.findNeighbours(owners[0],1,2,neigh_n,isghost_n);
							nodeTop=1;
							if (neigh_n.size()==0) bound=1;
						}
					}
					else if (iface==3){
						gridCart.findNeighbours(owners[0],0,1,neigh_f1,isghost_f1);
						gridCart.findNeighbours(owners[0],1,1,neigh_f2,isghost_f2);
						gridCart.findNeighbours(owners[0],2,2,neigh_n,isghost_n);
						if (neigh_n.size()==0) {
							gridCart.findNeighbours(owners[0],3,2,neigh_n,isghost_n);
							nodeTop=1;
							if (neigh_n.size()==0) bound=1;
						}
					}
				}
				// if owners[1] is the finest octant
				else {
					iface=gridCart.getFace(inter);
					NormalX=-gridCart.getNormal(inter)[0];
					NormalY=-gridCart.getNormal(inter)[1];

					if (iface==0){
						node_n1=gridCart.getNode(owners[0],1);
						node_n2=gridCart.getNode(owners[0],3);
						if (node_n1[0]==nodes[0][0] && node_n1[1]==nodes[0][1]) {
							gridCart.findNeighbours(owners[0],1,2,neigh_n,isghost_n);
							gridCart.findNeighbours(owners[0],2,1,neigh_f2,isghost_f2);
							face1=1;
							if (neigh_f2.size()==2) {
								neigh_f2[0]=neigh_f2[1];
								isghost_f2[0]=isghost_f2[1];
							}
						}
						else if (node_n2[0]==nodes[1][0] && node_n2[1]==nodes[1][1]){
							gridCart.findNeighbours(owners[0],3,2,neigh_n,isghost_n);
							gridCart.findNeighbours(owners[0],3,1,neigh_f2,isghost_f2);
							nodeTop=1;
							if (neigh_f2.size()==2) {
								neigh_f2[0]=neigh_f2[1];
								isghost_f2[0]=isghost_f2[1];
							}
						}
						if (neigh_n.size()==0) bound=1;
						gridCart.findNeighbours(owners[0],1,1,neigh_f1,isghost_f1);
						if (isghost_f1[0]==1) g_n3=gridCart.getGhostGlobalIdx(neigh_f1[0]);
						else g_n3=gridCart.getGlobalIdx(neigh_f1[0]);
						if (g_n3==g_n1) {
							neigh_f1[0]=neigh_f1[1];
							isghost_f1[0]=isghost_f1[1];
						}
					}
					else if (iface==1){
						node_n1=gridCart.getNode(owners[0],0);
						node_n2=gridCart.getNode(owners[0],2);
						if (node_n1[0]==nodes[0][0] && node_n1[1]==nodes[0][1]) {
							gridCart.findNeighbours(owners[0],0,2,neigh_n,isghost_n);
							gridCart.findNeighbours(owners[0],2,1,neigh_f2,isghost_f2);
							face1=1;
						}
						else if (node_n2[0]==nodes[1][0] && node_n2[1]==nodes[1][1]){
							gridCart.findNeighbours(owners[0],2,2,neigh_n,isghost_n);
							gridCart.findNeighbours(owners[0],3,1,neigh_f2,isghost_f2);
							nodeTop=1;
						}
						if (neigh_n.size()==0) bound=1;
						gridCart.findNeighbours(owners[0],0,1,neigh_f1,isghost_f1);
						if (isghost_f1[0]==1) g_n3=gridCart.getGhostGlobalIdx(neigh_f1[0]);
						else g_n3=gridCart.getGlobalIdx(neigh_f1[0]);
						if (g_n3==g_n1) {
							neigh_f1[0]=neigh_f1[1];
							isghost_f1[0]=isghost_f1[1];
						}
					}
					else if (iface==2){
						node_n1=gridCart.getNode(owners[0],2);
						node_n2=gridCart.getNode(owners[0],3);
						if (node_n1[0]==nodes[0][0] && node_n1[1]==nodes[0][1]) {
							gridCart.findNeighbours(owners[0],2,2,neigh_n,isghost_n);
							gridCart.findNeighbours(owners[0],0,1,neigh_f2,isghost_f2);
							face1=1;
							if (neigh_f2.size()==2) {
								neigh_f2[0]=neigh_f2[1];
								isghost_f2[0]=isghost_f2[1];
							}
						}
						else if (node_n2[0]==nodes[1][0] && node_n2[1]==nodes[1][1]){
							gridCart.findNeighbours(owners[0],3,2,neigh_n,isghost_n);
							gridCart.findNeighbours(owners[0],1,1,neigh_f2,isghost_f2);
							nodeTop=1;
							if (neigh_f2.size()==2) {
								neigh_f2[0]=neigh_f2[1];
								isghost_f2[0]=isghost_f2[1];
							}
						}
						if (neigh_n.size()==0) bound=1;
						gridCart.findNeighbours(owners[0],3,1,neigh_f1,isghost_f1);
						if (isghost_f1[0]==1) g_n3=gridCart.getGhostGlobalIdx(neigh_f1[0]);
						else g_n3=gridCart.getGlobalIdx(neigh_f1[0]);
						if (g_n3==g_n1) {
							neigh_f1[0]=neigh_f1[1];
							isghost_f1[0]=isghost_f1[1];
						}
					}
					else if (iface==3){
						node_n1=gridCart.getNode(owners[0],0);
						node_n2=gridCart.getNode(owners[0],1);
						if (node_n1[0]==nodes[0][0] && node_n1[1]==nodes[0][1]) {
							gridCart.findNeighbours(owners[0],0,2,neigh_n,isghost_n);
							gridCart.findNeighbours(owners[0],0,1,neigh_f2,isghost_f2);
							face1=1;
						}
						else if (node_n2[0]==nodes[1][0] && node_n2[1]==nodes[1][1]){
							gridCart.findNeighbours(owners[0],1,2,neigh_n,isghost_n);
							gridCart.findNeighbours(owners[0],1,1,neigh_f2,isghost_f2);
							nodeTop=1;
						}
						if (neigh_n.size()==0) bound=1;
						gridCart.findNeighbours(owners[0],2,1,neigh_f1,isghost_f1);
						if (isghost_f1[0]==1) g_n3=gridCart.getGhostGlobalIdx(neigh_f1[0]);
						else g_n3=gridCart.getGlobalIdx(neigh_f1[0]);
						if (g_n3==g_n1) {
							neigh_f1[0]=neigh_f1[1];
							isghost_f1[0]=isghost_f1[1];
						}
					}
				} // end ifiner
				if (!bound) {
					if (isghost_n[0]==1) {
						g_nn=gridCart.getGhostGlobalIdx(neigh_n[0]);
						Oct_n=gridCart.getGhostOctant(neigh_n[0]);
					}	
					else {
						g_nn=gridCart.getGlobalIdx(neigh_n[0]);
						Oct_n=gridCart.getOctant(neigh_n[0]);
					}

					if (isghost_f2[0]==1) {
						g_n2=gridCart.getGhostGlobalIdx(neigh_f2[0]);
						Oct_f2=gridCart.getGhostOctant(neigh_f2[0]);
					}	
					else {
						g_n2=gridCart.getGlobalIdx(neigh_f2[0]);
						Oct_f2=gridCart.getOctant(neigh_f2[0]);
					}
					if (isghost_f1[0]==1) {
						g_n3=gridCart.getGhostGlobalIdx(neigh_f1[0]);
						Oct_f1=gridCart.getGhostOctant(neigh_f1[0]);
					}	
					else {
						g_n3=gridCart.getGlobalIdx(neigh_f1[0]);
						Oct_f1=gridCart.getOctant(neigh_f1[0]);
					}

					if (nodeTop){
						interpt = interp4p(gridCart, nodes[1], owners[0], OctOwner, Oct_f2, Oct_n);
						interpb = interp3p(gridCart, nodes[0], owners[0], OctOwner, Oct_f1);
						if (TautbX == 0) {
							value0N = kappa*dt/gridCart.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[3];
							value02 = kappa*dt/gridCart.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
							value03 = -kappa*dt/gridCart.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];

							value1N = -kappa*dt/gridCart.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[3];
							value12 = -kappa*dt/gridCart.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
							value13 = kappa*dt/gridCart.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];
						}
						else {
							value0N = kappa*dt/gridCart.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[3];
							value02 = kappa*dt/gridCart.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
							value03 = -kappa*dt/gridCart.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];

							value1N = -kappa*dt/gridCart.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[3];
							value12 = -kappa*dt/gridCart.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
							value13 = kappa*dt/gridCart.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];
						}
					}
					else if (face1) {
						interpt = interp3p(gridCart, nodes[1], owners[0], OctOwner, Oct_f1);
						interpb = interp4p(gridCart, nodes[0], owners[0], OctOwner, Oct_f2, Oct_n);
						if (TautbX == 0) {
							value0N = -kappa*dt/gridCart.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[3];
							value02 = -kappa*dt/gridCart.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];
							value03 = kappa*dt/gridCart.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];

							value1N = kappa*dt/gridCart.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[3];
							value12 = kappa*dt/gridCart.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];
							value13 = -kappa*dt/gridCart.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
						}
						else {
							value0N = -kappa*dt/gridCart.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[3];
							value02 = -kappa*dt/gridCart.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];
							value03 = kappa*dt/gridCart.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];

							value1N = kappa*dt/gridCart.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[3];
							value12 = kappa*dt/gridCart.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];
							value13 = -kappa*dt/gridCart.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
						}
					}
					else {
						interpt = interp3p(gridCart, nodes[1], owners[0], OctOwner, Oct_f2);
						interpb = interp4p(gridCart, nodes[0], owners[0], OctOwner, Oct_f1, Oct_n);
						if (TautbX == 0) {
							value0N = -kappa*dt/gridCart.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[3];
							value02 = kappa*dt/gridCart.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
							value03 = -kappa*dt/gridCart.getVolume(owners[0])*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];

							value1N = kappa*dt/gridCart.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[3];
							value12 = -kappa*dt/gridCart.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpt[2];
							value13 = kappa*dt/gridCart.getVolume(OctOwner)*(NormalX*TaucenterY)/(TaucenterX*TautbY)*interpb[2];
						}
						else {
							value0N = -kappa*dt/gridCart.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[3];
							value02 = kappa*dt/gridCart.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
							value03 = -kappa*dt/gridCart.getVolume(owners[0])*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];

							value1N = kappa*dt/gridCart.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[3];
							value12 = -kappa*dt/gridCart.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpt[2];
							value13 = kappa*dt/gridCart.getVolume(OctOwner)*(NormalY*TaucenterX)/(TaucenterY*TautbX)*interpb[2];
						}
					}

					if (gridCart.getIsGhost(inter)){
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

					if (!gridCart.getIsGhost(inter)) {
						if (TautbX ==0 ){
							value00 = - kappa*dt/gridCart.getVolume(owners[0])*NormalX/TaucenterX*(-gridCart.getArea(inter)/deltaCenter + TaucenterY/TautbY*(interpb[0]-interpt[0]));
							value01 = - kappa*dt/gridCart.getVolume(owners[0])*NormalX/TaucenterX*(gridCart.getArea(inter)/deltaCenter + TaucenterY/TautbY*(interpb[1]-interpt[1]));
							value10 = - kappa*dt/gridCart.getVolume(OctOwner)*NormalX/TaucenterX*(gridCart.getArea(inter)/deltaCenter - TaucenterY/TautbY*(interpb[0]-interpt[0]));
							value11 = - kappa*dt/gridCart.getVolume(OctOwner)*NormalX/TaucenterX*(-gridCart.getArea(inter)/deltaCenter - TaucenterY/TautbY*(interpb[1]-interpt[1]));
						}
						else {
							value00 = - kappa*dt/gridCart.getVolume(owners[0])*NormalY/TaucenterY*(-gridCart.getArea(inter)/deltaCenter + TaucenterX/TautbX*(interpb[0]-interpt[0]));
							value01 = - kappa*dt/gridCart.getVolume(owners[0])*NormalY/TaucenterY*(gridCart.getArea(inter)/deltaCenter + TaucenterX/TautbX*(interpb[1]-interpt[1]));
							value10 = - kappa*dt/gridCart.getVolume(OctOwner)*NormalY/TaucenterY*(gridCart.getArea(inter)/deltaCenter - TaucenterX/TautbX*(interpb[0]-interpt[0]));
							value11 = - kappa*dt/gridCart.getVolume(OctOwner)*NormalY/TaucenterY*(-gridCart.getArea(inter)/deltaCenter - TaucenterX/TautbX*(interpb[1]-interpt[1]));
						}
					}
					else {
						if (TautbX ==0 ){
							value00 = - kappa*dt/gridCart.getVolume(owners[0])*NormalX/TaucenterX*(-gridCart.getArea(inter)/deltaCenter + TaucenterY/TautbY*(interpb[0]-interpt[0]));
							value10 = 0;
							value01 = - kappa*dt/gridCart.getVolume(owners[0])*NormalX/TaucenterX*(gridCart.getArea(inter)/deltaCenter + TaucenterY/TautbY*(interpb[1]-interpt[1]));
							value11= 0;
						}
						else {
							value00 = - kappa*dt/gridCart.getVolume(owners[0])*NormalY/TaucenterY*(-gridCart.getArea(inter)/deltaCenter + TaucenterX/TautbX*(interpb[0]-interpt[0]));
							value10 = 0;
							value01 = - kappa*dt/gridCart.getVolume(owners[0])*NormalY/TaucenterY*(gridCart.getArea(inter)/deltaCenter + TaucenterX/TautbX*(interpb[1]-interpt[1]));
							value11 = 0;
						}
					}
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n1, &value01, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n0, &value10, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n1, &value11, ADD_VALUES);

				} // end of !bound
				else{
					if (!gridCart.getIsGhost(inter)){
						value00 = kappa*dt/gridCart.getVolume(owners[0])*gridCart.getArea(inter)/deltaCenter;
						value11 = kappa*dt/gridCart.getVolume(OctOwner)*gridCart.getArea(inter)/deltaCenter;
						value01 = -kappa*dt/gridCart.getVolume(owners[0])*gridCart.getArea(inter)/deltaCenter;
						value10 = -kappa*dt/gridCart.getVolume(OctOwner)*gridCart.getArea(inter)/deltaCenter;
					}
					else {
						value00 = kappa*dt/gridCart.getVolume(owners[0])*gridCart.getArea(inter)/deltaCenter;
						value11 = 0;
						value01 = -kappa*dt/gridCart.getVolume(owners[0])*gridCart.getArea(inter)/deltaCenter;
						value10 = 0;
					}
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n1, &value01, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n0, &value10, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n0, 1, &g_n0, &value00, ADD_VALUES);
					MatSetValues(Mat_poisson, 1, &g_n1, 1, &g_n1, &value11, ADD_VALUES);
				}
			} // end of level difference
		} // end of !boundary
	} // end of intersection loop

	for (int j = 0; j < gridCart.getNumOctants(); ++j){
		if ( sqrt(pow(gridCart.getCenter(j)[0],2)+ pow(gridCart.getCenter(j)[1],2)) <= m_radius )
			m_pen=1;
		else
			m_pen =0;

		// Find Global index
		int g_j=gridCart.getGlobalIdx(j);
		value00= 1+m_pen*m_lambda*dt;

		MatSetValues(Mat_poisson, 1, &g_j, 1, &g_j, &value00, ADD_VALUES);
	}

	MatAssemblyBegin(Mat_poisson,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Mat_poisson,MAT_FINAL_ASSEMBLY);
	//MatView(Mat_poisson, PETSC_VIEWER_STDOUT_WORLD);
};


vector<int> HeatEq::getNeighbours(PabloUniform& gridCart,int i){
	vector<int> result(6,0);
	vector<uint32_t> owners(2,0);
	bool finer=gridCart.getFiner(gridCart.getIntersection(i));
	int int_owner=0;

	if (finer) int_owner=1;

	owners=gridCart.getOwners(gridCart.getIntersection(i));
	result[0]=owners[0+int_owner];
	result[1]=owners[1-int_owner];

	if (gridCart.getIsGhost(gridCart.getIntersection(i)))
	{
			if (finer){
				result[0]=-1;
				result[1]=owners[0];
			}
			else{
				result[1]=-1;
				result[0]=owners[0];

			}
			result[2]=owners[0+int_owner]*int_owner;
			result[3]=owners[1-int_owner]*(1-int_owner);
			result[4]=int_owner;
			result[5]=1-int_owner;
	}

	if (owners[0]==owners[1] && !gridCart.getIsGhost(gridCart.getIntersection(i))) 
	{
		if (finer)
		{
			result[0]=owners[1];
			result[1]=-1;
		}
		else{
			result[0]=owners[0];
			result[1]=-1;
		}
	}

	return result;
};


vector<double> HeatEq::interp4p(PabloUniform& gridCart, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2, Octant*& Oct3){
	vector<double> value(4,0);
	double x0, y0, x1, y1, x2, y2, x3, y3;
	double intermed1, xdif, ydif, xydif, intermed2, intermed3, intermed4;
	double long1, long2, long3, long4, long5, long6, long7, long8, long9, long10, long11, long12, long13, long14;

	x0=gridCart.getCenter(owner)[0]; y0=gridCart.getCenter(owner)[1];
	x1=gridCart.getCenter(Oct1)[0]; y1=gridCart.getCenter(Oct1)[1];
	x2=gridCart.getCenter(Oct2)[0]; y2=gridCart.getCenter(Oct2)[1];
	x3=gridCart.getCenter(Oct3)[0]; y3=gridCart.getCenter(Oct3)[1];

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

vector<double> HeatEq::interp3p(PabloUniform& gridCart, darray3 nodes, int owner, Octant*& Oct1, Octant*& Oct2){
	vector<double> value(3,0);
	double x0, y0, x1, y1, x2, y2;
	double intermed1, xdif, ydif, intermed2;

	x0=gridCart.getCenter(owner)[0]; y0=gridCart.getCenter(owner)[1];
	x1=gridCart.getCenter(Oct1)[0]; y1=gridCart.getCenter(Oct1)[1];
	x2=gridCart.getCenter(Oct2)[0]; y2=gridCart.getCenter(Oct2)[1];

	xdif=x0-x1; ydif=y0-y1;
	intermed1=(x1-x2)/(x0-x1); intermed2=intermed1*(y0-y1)-(y1-y2);

	value[0] = 1 - x0/xdif + (x0*ydif)/(xdif*intermed2)*intermed1 - y0/intermed2*intermed1 + nodes[0]/xdif*(1-intermed1/intermed2*ydif) + nodes[1]/intermed2*intermed1;
	value[1] = x0/xdif + (x0*ydif)/(xdif*intermed2)*(-1-intermed1) + y0/intermed2*(1+intermed1) + nodes[0]/xdif*((1+intermed1)/intermed2*ydif-1) + nodes[1]/intermed2*(-1-intermed1);
	value[2] = (x0*ydif)/(xdif*intermed2) - y0/intermed2 - nodes[0]/xdif*ydif/intermed2 + nodes[1]/intermed2;

	return value;
};