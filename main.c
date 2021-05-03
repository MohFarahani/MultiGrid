#include<iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include"solver.h"
#include"Grid.h"

void exact_solution(Grid& u , int type);
double l2norm_error(const Grid& u_exact , const Grid& u_approx);
//===================================================================================================================
int main( int argc, char** argv )
{
	   if(argc<4)
	    {
	      std::cerr << "mgbs <l> <n> <b>" << std::endl;
	     return EXIT_FAILURE;
	    }
	    int  L = atoi(argv[1]);
	    int  n = atoi(argv[2]);
	    int b = atoi(argv[3]);    // boundary type: "1" count as "Dirichlet" and other is Neumman
	    int Size ;
	    Size = (1<<L)+1 ;
	    Grid u(Size,Size) 	;      // "u" is our approximate solution (in agree with the notation of assignment sheet )
	    int Pre(2),Post(1) ;


	    Solver solve_mg(n,L,Pre,Post,b) ;

            timeval t;
	    //start of time measurement
    	    gettimeofday(&t, NULL);
    	    double start = t.tv_sec + t.tv_usec * 1.0e-6;
	    //solve
	    solve_mg.mg(u );
    	    gettimeofday(&t, NULL);
    	    double end = t.tv_sec + t.tv_usec * 1.0e-6;
	    //end of time measurement
	    solve_mg.write_plot(u,"solution");

	    Grid u_exact(Size,Size) ;
	    exact_solution(u_exact,b) ;



	    double l2n_error ;
	    l2n_error = l2norm_error(u_exact , u);
	    std::cout<<"The L2 norm of error is : "<<l2n_error<<std::endl;

	    double dt = end - start ;
	    std::cout<<"The time is  : "<<dt<<"  seconds" <<std::endl;


	    // In order to find error for level = 3:8 for both boundary ( just make "condition= true")
	    bool condition = false ;
	    if( condition == true )
	    {
	    	std::ofstream out ("Error.txt");
		    out<<"level : "<<"\t"<<"l2n_error_dirichlet"<<"\t"<<std::endl;
	    	for (int i=3 ;  i!= 9 ; ++i )
	    	{
			    int size ;
			    size = (1<<i)+1 ;
			    Grid u_dirichlet(size,size)  ;
			    int pre(2),post(1) ;
			    Solver mg_dirichlet(n,i,pre,post,1) ;
			    mg_dirichlet.mg(u_dirichlet );
			    Grid u_exact_dirichlet(size,size) ;

			    exact_solution(u_exact_dirichlet,1) ;
			    double l2n_error_dirichlet ;
			    l2n_error_dirichlet= l2norm_error(u_exact_dirichlet , u_dirichlet);
			    out<<"level : "<<i<<"l2n_error_dirichlet"<<l2n_error_dirichlet<<std::endl;


	    	}
		    out<<"******************************************************************************"<<std::endl;
		    out<< std::endl ;
		    out<<"level : "<<"\t"<<"l2n_error_neuman"<<"\t"<<std::endl;

	    	for (int i=3 ; i != 9 ; ++i )
	    	{
			    int size ;
			    size = (1<<i)+1 ;
			    Grid u_neuman(size,size)  ;
			    int pre(2),post(1) ;
			    Solver mg_neuman(n,i,pre,post,2) ;
			    mg_neuman.mg(u_neuman );
			    Grid u_exact_neuman(size,size) ;
			    exact_solution(u_exact_neuman,2) ;
			    double l2n_error_neuman  ;
			    l2n_error_neuman= l2norm_error(u_exact_neuman , u_neuman);
			    out<<"level : "<<i<<"l2n_error_neuman"<<l2n_error_neuman<<std::endl;

	    	}
	    	   out.close();

	    }
}

//===================================================================================================================
void exact_solution(Grid& u , int boundary_type){

	if (boundary_type==1)
	{
         for(int j=0; j<u.y();j++)
        	 for(int i=0; i<u.x();i++)
        	 {
        		 	 u(i,j) = sin(M_PI*i*u.dx())*sinh(M_PI*j*u.dy());
        	 }

	}
	else
	{
	  for( int j=0 ; j<u.y() ; j++ )
	     for( int i=0 ; i<u.x() ; i++ )
	       u(i,j) =i*u.dx()*(1-i*u.dx());

	}

}
//===================================================================================================================
double l2norm_error(const Grid& u_exact , const Grid& u_approx){

	Grid error(u_exact);
	error.fill(0.0);
	int M = error.x()-1 ;
	int N = error.y()-1 ;
	for (int j = 1 ; j != N ; ++j )
	{
		for (int i = 1 ; i != M ; ++i )
		{
			error(i,j) = u_exact(i,j) -u_approx(i,j) ;
		}
	}

	double l2n(0.0) ;
	  for(int j=1; j != N ; j++)
	   {
	     for(int i=1 ; i !=M ; i++)
	     {
	    	 l2n += error(i,j)*error(i,j);
	     }
	   }
	   l2n = sqrt(l2n/(M*N)) ;
	   return l2n ;

}
//===================================================================================================================

