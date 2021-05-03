#include"Grid.h"
#include<vector>
#include<iostream>
#include <fstream>
#include <string>
#include<cmath>
#define _USE_MATH_DEFINES


class Solver
{
  public:
	// Default Constructor
	   Solver(int cycle , int level , int pre , int post , int boundary);

	   void restriction( const Grid& res_fine , Grid& f_coarse);
	   void prolongation_correction(const Grid& coarse , Grid& fine );
	   void residual (const Grid& u, const Grid& f , Grid& r);
	   void l2_norm(double& l2n ,const  Grid& r);
	   void red_black_GS(Grid&,const Grid& ,const int& iteration);
	   void mg(Grid& u);
	   void mgm(int level);


	   // Boundary
	   void set_dirichlet_boundary(Grid& u) ;
	   void set_neuman_boundary(Grid& u) ;
	   // rhs
	   void set_rhs_for_DB(Grid& f) ;
	   void set_rhs_for_NB(Grid& f) ;
	   // write
	   void write_plot( Grid u ,  const std::string& filename ) ;

  private:
	   std::vector<Grid> u_ ; // contain whole "u_" for all level
	   std::vector<Grid> f_ ;
	   std::vector<Grid> res_ ;


	   int cycle_ ;
	   int level_;
	   int pre_ ;
	   int post_;
	   int boundary_type_ ;
};
