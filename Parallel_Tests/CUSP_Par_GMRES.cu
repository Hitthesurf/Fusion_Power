#include <cusp/coo_matrix.h>
#include <cusp/monitor.h>
#include <cusp/io/matrix_market.h>
#include <cusp/krylov/gmres.h>
#include <iostream>
#include "My_StopWatch.hpp"


// where to perform the computation
typedef cusp::device_memory MemorySpace;
// which floating point type to use
typedef float ValueType;
int main(void)
{
	StopWatch My_Watch;
    // create an empty sparse matrix structure (coo format)
	cusp::coo_matrix<int, ValueType, MemorySpace> A;
	
    // load data from matrix file
	cusp::io::read_matrix_market_file(A, "thermal1.mtx");
	
    // allocate storage for solution (x) and right hand side (b)
    cusp::array1d<ValueType, MemorySpace> x(A.num_rows);
    cusp::array1d<ValueType, MemorySpace> b(A.num_rows);
	cusp::io::read_matrix_market_file(b, "thermal1_b.mtx");

	
    // set stopping criteria:
	std::cout<<"Start Solve"<<"\n";
	int its = 10000;
	float tols[4] = {1e-02, 1e-03, 1e-04, 1e-05};
	bool verbose = true;
	for (int i=0; i<4; ++i)
	{
		float tol = tols[i];
		std::cout<<"-----START-tol="<<tol<<"----------"<<std::endl;
		// set initial guess
		thrust::fill( x.begin(), x.end(), ValueType(0));
		
	    cusp::monitor<ValueType> monitor(b, its, tol, 0, verbose);
		int restart = 50;
		// solve the linear system A * x = b with the GMRES
		//Start Timer
		My_Watch.Start();
		cusp::krylov::gmres(A, x, b,restart, monitor);
		My_Watch.Stop();
		std::cout<<"Done Solve: " << My_Watch.Time() << std::endl;
	}

    return 0;
}