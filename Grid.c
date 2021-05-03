#include "Grid.h"




//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================
Grid::Grid()
{
   // Default construct
	xSize_ = 0 ;
	ySize_ = 0 ;
	vector_.clear();
	vector_.resize( 1,0.0) ;
}



Grid::Grid( const int& x_points, const int & y_points)

{
	// construct 2D array
	xSize_ = x_points;
	ySize_ = y_points ;

	vector_.clear();
	vector_.resize(xSize_*ySize_,0.0);
}


Grid::Grid(const Grid& s): xSize_(s.x()),ySize_ (s.y())
{
	vector_.clear();
	vector_= s.get_vector() ;
}


//===================================================================================================================
//
//  Overloaded operator
//
//===================================================================================================================
// assignment operator
Grid& Grid::operator= (const Grid& s)
{
	xSize_ = s.x() ;
	ySize_ = s.y() ;
	vector_.clear();
	vector_ = s.get_vector();
	return *this;
}
// Operator()
double& Grid::operator ()(int i,int j)
{
	assert (i <= xSize_ ) ;
	assert( i>= 0 ) ;
	assert (j <= ySize_ ) ;
	assert( j>= 0 ) ;

   return vector_[j*xSize_ + i ];
}


// Operator() const
const double& Grid::operator ()(int i,int j) const
{
	assert (i <= xSize_ ) ;
	assert( i>= 0 ) ;
	assert (j <= ySize_ ) ;
	assert( j>= 0 ) ;
   return vector_[j*xSize_ + i ];
}


//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================

//===================================================================================================================
//
double Grid::dx()
{
	double hx ;
	hx = (double) 1.0/(xSize_ -1);
	return hx ;

}
double Grid::dy()
{
	double hy ;
	hy = (double) 1.0/(ySize_-1) ;
	return hy ;
}
//===================================================================================================================
//initialize the whole array with a constant value
void Grid::fill( double value )
{
	  std::fill(vector_.begin(),vector_.end(),value);
}

//===================================================================================================================
// Print the whole array (for debugging purposes)
void Grid::print()
{

	for (int j=ySize_-1 ; j != -1 ; --j)
	{
		for (int i=0 ; i != xSize_ ; ++i)
		{
			std::cout<< vector_[j*xSize_+i] <<"\t";

		}
			std::cout<< std::endl;
	}

}

//===================================================================================================================
// return the vector_
std::vector<double> Grid::get_vector () const
{
	return vector_ ;
}

//===================================================================================================================
//

void Grid::Resize(const int& x_points, const int & y_points)
{
	xSize_ = x_points ;
	ySize_ = y_points ;
	vector_.clear();
	vector_.resize(xSize_*ySize_,0.0);
}

//===================================================================================================================


