#include"solver.h"
//===================================================================================================================
Solver::Solver(int cycle , int level , int pre , int post , int boundary) :cycle_(cycle),level_(level), pre_(pre), post_(post),boundary_type_(boundary)
{
	u_.resize(level)   ;
	f_.resize(level)   ;
	res_.resize(level) ;

	for (int i = 0 ; i != level ; ++i)
	{
		int points ;
		points = ( 1<<(i+1) ) + 1 ;
		u_[i].Resize(points, points) ;
		f_[i].Resize(points, points) ;
		res_[i].Resize(points, points) ;
		u_[i].fill(0.0) ;
		f_[i].fill(0.0) ;
		res_[i].fill(0.0) ;

	}
}

//===================================================================================================================
void Solver::restriction(const Grid& res_fine , Grid& f_coarse)
{

	 for(int i=1 ; i<f_coarse.x()-1 ; i++)
    {
        for(int j=1 ; j<f_coarse.y()-1 ; j++)
        {
            int two_i = 2*i;
            int two_j = 2*j;

           f_coarse(i,j)=
        		   	     	 0.0625*res_fine(two_i-1,two_j +1)	 +  0.125*res_fine(two_i,two_j +1)	+	0.0625*res_fine(two_i+1,two_j+1)
           	   	   	   	   +
           	   	   	   	 0.125*res_fine(two_i-1,two_j )	 	 +  0.25*res_fine(two_i,two_j )		+	0.125*res_fine(two_i+1,two_j )
           	   	   	   	   +
           	   	   	   	 0.0625*res_fine (two_i-1,two_j -1)	 +  0.125*res_fine(two_i,two_j -1)	+	0.0625*res_fine(two_i+1,two_j -1)	;
        }
    }


}
//===================================================================================================================
void Solver::prolongation_correction(const Grid& coarse , Grid& fine )
{

    for(int i=1 ; i<coarse.x()-1 ; i++)
    {
        int two_i = 2*i;

        //calculate values at the lower nodes
        fine (two_i,1)+=0.5*(coarse(i,0)+coarse(i,1));
        fine (two_i+1,1)+=0.25*(coarse(i,0)+coarse(i+1,0)+coarse(i,1)+coarse(i+1,1));

        //calculate values at the inner nodes
        for(int j=1 ; j<coarse.y()-1 ; j++)
        {
            int two_j = 2*j;
            fine (two_i,two_j)+=coarse(i,j);
            fine (two_i+1,two_j)+=0.5*(coarse(i,j)+coarse(i+1,j));
            fine (two_i,two_j+1)+=0.5*(coarse(i,j)+coarse(i,j+1));
            fine (two_i+1,two_j+1)+=0.25*(coarse(i,j)+coarse(i+1,j)+coarse(i,j+1)+coarse(i+1,j+1));
        }
     }

    //calculate values at the upper nodes
    for(int j=1 ; j<coarse.y()-1 ; j++)
    {
        int two_j = 2*j;
        fine (1,two_j)	 +=0.5*(coarse(0,j)+coarse(1,j));
        fine (1,two_j+1) +=0.25*(coarse(0,j)+coarse(1,j)+coarse(0,j+1)+coarse(1,j+1));
    }

    //calculate values at the up-down nodes
    fine (1,1)+=0.25*(coarse(0,0)+coarse(1,0)+coarse(0,1)+coarse(1,1));

}

//===================================================================================================================
void Solver::residual(const Grid& u, const Grid& f , Grid& r)
{

	double dx2_inverse =(double) (1.0/r.dx())*(1.0/r.dx()) ;
	double dy2_inverse =(double) (1.0/r.dy())*(1.0/r.dy()) ;
	double coefficient = 2.0 * (dx2_inverse  + dy2_inverse ) ;


	for (int j = 1 ; j != u.y()-1 ; ++j )
	{
		for (int i = 1 ; i != u.x()-1 ; ++i )
		{
			r(i,j) = f(i,j) -
								(
								   - dx2_inverse*( u(i+1,j)+ u(i-1,j) )
								   - dy2_inverse*( u(i,j+1)+ u(i,j-1) )
								   + coefficient*u(i,j)
								);


		}

	}
	

}
//===================================================================================================================
void Solver::l2_norm(double& l2n , const Grid& res)
{
	
	  int M = res.x()-1;
	  int N = res.y()-1;
          double dx = 1.0 / M ;
	  l2n = 0.0 ;
	  for(int j=1; j != N ; j++)
	   {
	     for(int i=1 ; i != M ; i++)
	     {
	    	 l2n += res(i,j)*res(i,j);
	     }
	   }

	   l2n = dx*sqrt( l2n)  ;


}
//===================================================================================================================
void Solver::red_black_GS(Grid& u , const Grid& f , const int& iteration)
{

	double dx2_inverse =(double) (1.0/u.dx())*(1.0/u.dx()) ;
	double dy2_inverse =(double) (1.0/u.dy())*(1.0/u.dy()) ;
	double coefficient = 2.0 * (dx2_inverse  + dy2_inverse ) ;
        double coefficient_inverse = (double) 1.0 / coefficient ;
    	int M = u.x()-1;
    	int N = u.y()-1;

	

		//Loop
	  for(int itr=0; itr<iteration; itr++)
	    
  	{
   
           //red loop
	   for(int j=1; j<N;j++)
	   {
	     for(int i=1; i<M;i++)
	     {
	       if(((i+j)%2)==0)
	       {
		 double uE, uW, uN, uS;
		 uW=u(i-1,j);
		 uE=u(i+1,j);
		 uN=u(i,j+1);
		 uS=u(i,j-1);
		 u(i,j) = ( f(i,j) +
				 ( (dx2_inverse * (uW+uE) ) + (dy2_inverse * (uS+uN) ) )
				  )*  coefficient_inverse ;	 }
	     }
	   }
	   //black loop
	   for( int j=1 ; j<N ; j++)
	   {
	     for( int i=1 ; i<M ; i++)
	     {
	       if(((i+j)%2)!=0)
	       {
		 double uE, uW, uN, uS;
		 uW=u(i-1,j);
		 uE=u(i+1,j);
		 uN=u(i,j+1);
		 uS=u(i,j-1);
		 u(i,j) = ( f(i,j) +
				 	( (dx2_inverse * (uW+uE) ) + (dy2_inverse * (uS+uN) ) )
				  	)*  coefficient_inverse ;
	       }
	     }
	   }
	 }

}

//===================================================================================================================
void Solver::mg(Grid& u)
{
	int l = level_-1 ;

		if(boundary_type_ == 1)
		  {
			set_dirichlet_boundary(u_[l]);
			set_rhs_for_DB(f_[l]) ;
		  }
		else
		{
			set_neuman_boundary(u_[l]);
			set_rhs_for_NB(f_[l]) ;
		}
	

	double l2norm_new,l2norm_old ;
		Grid res( (1<<(level_) )+1 , (1<<(level_))+1 );
		res.fill(0.0) ;
		residual(u_[level_-1], f_[level_-1],res);
	double R(0.0)  ;
		l2_norm(R,res) ;
		l2norm_new = R ;
	for(int c = 0; c < cycle_ ; ++c)
	{
		mgm(level_ - 1);
		residual(u_[level_-1], f_[level_-1],res);


	    l2norm_old = l2norm_new ;
	    double r(0.0) ;
	    l2_norm(r,res) ;
	    l2norm_new = r ;
	    std::cout<<"L2 norm of residual after "<<c+1<<" v-cycle:"<< " in level :"<< level_ << "  ----> "<< r <<std::endl;
	    std::cout << "Convergence rate: " << l2norm_new/l2norm_old << std::endl;

	}
	u = u_[level_-1] ;
}
//===================================================================================================================
void Solver::mgm(int level)
{


	//pre-smoothing
	red_black_GS(u_[level], f_[level],pre_) ;

	
	//compute residual
	residual(u_[level], f_[level],res_[level]) ;

	//restriction on residual
	restriction(res_[level] , f_[level - 1] ) ;


	if(level == 1){
		red_black_GS(u_[0], f_[0],pre_) ;
	}
	
	else{
			// recusion		
			u_[level-1].fill(0.0) ;
			mgm(level- 1) ;
		}


	//prolongation and correction
	prolongation_correction(u_[level-1], u_[level]);

	//post-smoothing
	red_black_GS(u_[level], f_[level],post_);

}
//===================================================================================================================
//
void Solver::set_dirichlet_boundary( Grid& u)
{

	for(int i=0 ; i<u.x() ; i++)
    	{
			u(i,0)=0; //bottom boundary
			u(i,u.y()-1)=sin(M_PI*i*u.dx())*sinh(M_PI); // upper boundary
        }

        for(int j=1 ; j<u.y() ; j++)
        {
	       u(0,j)=0.0 ; 		// left boundary
	       u(u.x()-1,j)=0.0 ;   // right boundary
        }
}

//===================================================================================================================
//
void Solver::set_neuman_boundary(Grid& u)
{
	for(int i=0 ; i<u.x() ; i++)
    {
        u(i,0)=i*u.dx()*(1-(i*u.dx())); //bottom boundary
        u(i,u.y()-1)=i*u.dx()*(1-(i*u.dx())); // upper boundary
        }

        for(int j=1 ; j<u.y() ; j++)
        {

	       u(0,j)=u(1,j)-1*u.dx(); // left boundary
	       u(u.x()-1,j)=u(u.x()-2,j)-1*u.dx(); // right boundary
        }

}
//===================================================================================================================
void Solver::set_rhs_for_DB( Grid& f )
{
	f.fill(0.0) ;

}
//===================================================================================================================
void Solver::set_rhs_for_NB( Grid& f )
{
	f.fill(2.0) ;

}

//===================================================================================================================
void Solver::write_plot( Grid u ,  const std::string& filename )
{
   std::ofstream out(filename.c_str());

   for( int j=0; j<u.x(); ++j ) {
      for( int i=0; i<u.y(); ++i ) {
         out << i*u.dx() << "\t" << j*u.dy() << "\t" << u(i,j) << std::endl;
      }
      out << std::endl;
   }

   out.close();
}
