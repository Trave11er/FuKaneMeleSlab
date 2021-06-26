#include <iostream>
#include <armadillo>
#include <complex>
#include <cstdlib>
#include <ctime>
using namespace std;
using namespace arma;

int main() {
	// Input: sizeCells of lattice in z direction and # of k-points in the direction of kPointVec
	int sizeCells = 5;
	int k_points_max=1+19*1;
	rowvec kPointVec(3); // Defined below

	// Parameters for the Hamiltonian and data structures (Don't need to change those)
	double tHopping = 1;
	double deltatHopping = 0.4; 
	complex<double> spinOrbitCoupling(0,1); 

	// Misc
	double delta = 0.01; // small number for numerical comparisons
	int index =0;
	const complex<double> ii(0, 1);

	// Geometry
	double latVecNorm=1/sqrt(2); // cubic cell size a = 1
	double sublatVecNorm=sqrt(3)/4;
	int nSites=3*sizeCells;

	rowvec dummyVec(3), dummyVec2(3), dummyVec3(3), dummyVec4(3);
	cx_mat dummyMat(2,2);
	mat sublatOne(nSites, 4); // Sites at origin of a unit cell
	mat sublatTwo(nSites, 4); // Sites at sublatVec in a unit cell
	rowvec sublatVec(3); // 2nd lattice site within a unit cell
	rowvec latVec1(3), latVec2(3), latVec3(3), latVec4(3), latVec5(3); // lattice vectors
	rowvec bulklatVec1(3), bulklatVec2(3), bulklatVec3(3), GVec1(3), GVec2(3);
	sublatVec << 0.25 << 0.25 << 0.25; 
	bulklatVec1 << 0 << 0.5 << 0.5;
	bulklatVec2 << 0.5 << 0 << 0.5;
	bulklatVec3 << 0.5 << 0.5 << 0;

	// Variables below are needed for periodic b.c.s
	// Unit cell
	latVec1 << 0.5 << -0.5 << 0;
	latVec2 << 0 << 0.5 << -0.5;
	latVec3 << 1 << 1 << 1;
	// Two more atoms from sublattice One within the primary unit cell
	latVec4 << 0 << 0.5 << 0.5;
	latVec5 << 0.5 << 1 << 0.5;
	// Reciprocal lattice vectors for slab
	GVec1 << 4*M_PI/3 << -2*M_PI/3 << -2*M_PI/3;
	GVec2 << 2*M_PI/3 << 2*M_PI/3 << -4*M_PI/3;

	// Define Pauli matrices
	cx_mat unity(2,2), sigma_x(2,2), sigma_y(2,2), sigma_z(2,2); // Pauli matrices  
	mat dummy_unity(2,2), dummy_sigma_x(2,2), dummy_sigma_y(2,2), dummy_sigma_z(2,2);
	dummy_unity(0,0)=1; dummy_unity(1,1)=1; dummy_sigma_x(1,0)=1; dummy_sigma_x(0,1)=1;
	dummy_sigma_y(0,1)=-1; dummy_sigma_y(1,0)=1; dummy_sigma_z(0,0)=1; dummy_sigma_z(1,1)=-1;
	unity.set_real(dummy_unity);
	sigma_x.set_real(dummy_sigma_x);
	sigma_y.set_imag(dummy_sigma_y);
	sigma_z.set_real(dummy_sigma_z);

	// Set up data structures for the slab
	cx_mat Hamiltonian(4*nSites,4*nSites);
	Hamiltonian.zeros();
	vec eigvals(4*nSites);
	cx_mat eigvecs(4*nSites,4*nSites);
	mat bandstructure(k_points_max,4*nSites);

	// Set up primary unit cell of the slab
	for (int k=0; k<sizeCells; k++) {
		dummyVec=k*latVec3;
		sublatOne(index,0)=99; sublatOne(index+1,0)=99; sublatOne(index+2,0)=99;
		sublatOne(index,span(1,3))=dummyVec;
		sublatOne(index+1,span(1,3))=(dummyVec+latVec4);
		sublatOne(index+2,span(1,3))=(dummyVec+latVec5);
		sublatTwo(index,0)=99; sublatTwo(index+1,0)=99; sublatTwo(index+2,0)=99;
		sublatTwo(index,span(1,3))=(dummyVec+sublatVec);
		sublatTwo(index+1,span(1,3))=(dummyVec+latVec4+sublatVec);
		sublatTwo(index+2,span(1,3))=(dummyVec+latVec5+sublatVec);
		index=index+3;
	}

	// Loop over k-points
	for (int m=0; m<k_points_max; m++ ){
		kPointVec=GVec1*m/20;

		// Set up the Hamiltonian in a different way (results equivalent)
		for (int n=0; n<nSites; n++){
			// Go over unit cells in width (i,j; periodic), height (k) and across sites in each cell (l)
			for (int i=-2; i<3; i++) {
				for (int j=-2; j<3; j++) {
					for (int k=0; k<sizeCells; k++) {
						for (int l=0; l<3; l++) {
							switch(l) {
								case 0:
									dummyVec=i*latVec1+j*latVec2+k*latVec3-sublatOne(n,span(1,3));
									break;
								case 1:
									dummyVec=i*latVec1+j*latVec2+k*latVec3+latVec4-sublatOne(n,span(1,3));
									break;
								case 2:
									dummyVec=i*latVec1+j*latVec2+k*latVec3+latVec5-sublatOne(n,span(1,3));
									break;
								default:
									cout << "Error, l should be 0,1,2" << endl;
							}
							dummyVec2=dummyVec+sublatVec; // used to find hoppings between sublattices
							index=3*k+l; // used when have a site within primary cell 
							dummyVec3=i*latVec1+j*latVec2; // used if have an atom outside primary cell
							dummyVec4=dummyVec+sublatOne(n,span(1,3))-sublatTwo(n,span(1,3));

							// Set up nearest neighbour hoppings (between sublattices)
							// Set up hoppings into sublatOne
							if (dot(dummyVec2, dummyVec2) < sublatVecNorm*sublatVecNorm+delta) {
								// hopping to the n-th site from along [111] within the primary cell
								if (dot(dummyVec, dummyVec) < delta) {
									Hamiltonian(span(4*n+0,4*n+1),span(4*index+2,4*index+3))+=(tHopping+deltatHopping)*unity;
									if (n != index || i!=0 || j!=0) {cout << "Error for deltatHopping" << endl; }
									// hopping in other directions
								} else {
									if (i==0 && j==0) {
										// Hoppings withing the primary cell
										Hamiltonian(span(4*n+0,4*n+1),span(4*index+2,4*index+3))+=tHopping*unity;
										if (i!=0 || j!=0 || n==index) {cout << "Error for tHopping within primary cell" << endl; }
										// hopping to sublatOne from outside the primary cell
									} else {
										// hoppings from sublatTwo of index to sublatOne of n
										Hamiltonian(span(4*n+0,4*n+1),span(4*index+2,4*index+3))+=tHopping*exp(ii*dot(dummyVec3,kPointVec))*unity;
										if (n == index) {cout << "Error for tHopping" << endl; }
									}
								}
							}

							// Set up hoppings into sublatTwo
							if (dot(dummyVec4, dummyVec4) < sublatVecNorm*sublatVecNorm+delta) {
								// hopping to the n-th site from along [111] within the primary cell
								if (dot(dummyVec4+sublatVec,dummyVec4+sublatVec) < delta) {
									Hamiltonian(span(4*n+2,4*n+3),span(4*index+0,4*index+1))+=(tHopping+deltatHopping)*unity;
									if (n != index || i!=0 || j!=0) {cout << "Error for deltatHopping" << endl; }
									// hopping in other directions
								} else {
									if (i==0 && j==0) {
										// Hoppings withing the primary cell
										Hamiltonian(span(4*n+2,4*n+3),span(4*index+0,4*index+1))+=tHopping*unity;
										if (i!=0 || j!=0 || n==index) {cout << "Error for tHopping within primary cell" << endl; }
										// hopping to sublatTwo from outside the primary cell
									} else {
										// hoppings from sublatOne of index to sublatTwo of n
										Hamiltonian(span(4*n+2,4*n+3),span(4*index+0,4*index+1))+=tHopping*exp(ii*dot(dummyVec3,kPointVec))*unity;
										if (n == index) {cout << "Error for tHopping" << endl; }
									}
								}
							}
							if ( (dot(dummyVec,dummyVec) < latVecNorm*latVecNorm+delta)  && 
									(dot(dummyVec,dummyVec) > sublatVecNorm*sublatVecNorm+delta) ) { 
								// Matrix elements are direction-dependent
								if ( (dummyVec(0) > 0 && dummyVec(1) > 0) || (dummyVec(0) > 0 &&
											dummyVec(2) > 0) || (dummyVec(1) > 0 && dummyVec(2) > 0) ) {
									dummyVec = cross(sublatVec, dummyVec-sublatVec);
								} else if ( (dummyVec(1) < 0 && dummyVec(2) < 0) || (dummyVec(0) > 0 &&
											dummyVec(1) < 0) || (dummyVec(0) > 0 && dummyVec(2) < 0) ) {
									dummyVec = cross(-bulklatVec1+sublatVec, dummyVec+bulklatVec1-sublatVec);
								} else if ( (dummyVec(0) < 0 && dummyVec(1) > 0) || (dummyVec(0) < 0 &&
											dummyVec(2) < 0) || (dummyVec(1) > 0 && dummyVec(2) < 0) ) {
									dummyVec = cross(-bulklatVec2+sublatVec, dummyVec+bulklatVec2-sublatVec);
								} else if ( (dummyVec(0) < 0 && dummyVec(2) > 0) || (dummyVec(1) < 0 &&
											dummyVec(2) > 0) || (dummyVec(0) < 0 && dummyVec(1) < 0) ) {
									dummyVec = cross(-bulklatVec3+sublatVec, dummyVec+bulklatVec3-sublatVec);
								} else {cout << "Error in NNN hopping" << endl; }
								dummyMat=dummyVec(0)*sigma_x+dummyVec(1)*sigma_y+dummyVec(2)*sigma_z;
								if (i==0 && j==0) {
									Hamiltonian(span(4*n,4*n+1),span(4*index,4*index+1))+=spinOrbitCoupling*dummyMat;
									// Minus sign because for sublatTwo dummyVec is minus that of sublatOne
									Hamiltonian(span(4*n+2,4*n+3),span(4*index+2,4*index+3))+= -1.0*spinOrbitCoupling*dummyMat;
								} else {
									// Hoppings within sublatOne
									Hamiltonian(span(4*n,4*n+1),span(4*index,4*index+1))+=exp(ii*dot(dummyVec3, kPointVec))*spinOrbitCoupling*dummyMat;
									// Hoppings within sublatTwo
									Hamiltonian(span(4*n+2,4*n+3),span(4*index+2,4*index+3))+= -1.0*exp(ii*dot(dummyVec3,kPointVec))*spinOrbitCoupling*dummyMat;
								}

							}
						}
					}
				}
			}
		}  

		// abs returns a matrix, max returns a vector, 2nd max - largest value
		if (max(max(abs(Hamiltonian-Hamiltonian.t()))) > delta ) {
			cout << "Error, Hamiltonian is not Hermitian by at least " << delta << endl;
		}

		eig_sym(eigvals, eigvecs, Hamiltonian);
		bandstructure.row(m)=eigvals.t();
		cout << endl;
	//	Hamiltonian.save("Hamiltonian.txt", arma_ascii);
		cout << "Calculated k-point";
		kPointVec.print();
		Hamiltonian.zeros();
	}

	bandstructure.save("output_bs.txt", raw_ascii);
	//sublatOne.save("output.txt",raw_ascii);
	//sublatTwo.save("output2.txt",raw_ascii);
}

