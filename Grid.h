#ifndef GRID_H
#define GRID_H



#include <vector>
#include <algorithm>    // std::fill
#include <iostream>
#include <cassert>
//*******************************************************************************************************************
/*!  Grid class 2  dimensions
*
*    
*/
//*******************************************************************************************************************
class Grid
{
public:

   // Default Constructor
   Grid();
   ~Grid(){vector_.clear() ; }
   // Constructors
   Grid( const int& x_points, const int & y_points);


   // copy constructor
   Grid(const Grid& s);

   // assignment operator
   Grid& operator= (const Grid& s);


   // Access Operators
    double & operator () ( int i ,int j );
   // for const Arrays the following access operators are required
    const double&  operator () ( int i ,int j ) const;



   // initialize the whole array with a constant value
   void fill( double value );



   // return xSize &  ySize
   int x() {return xSize_ ; }
   int y() {return ySize_ ; }
   const int x() const {return xSize_ ; }
   const int y() const {return ySize_ ; }

   // return the "dx" and "dy" (one grid size)
   double dx() ;
   double dy() ;



   // resize the vector_
   std::vector<double>  get_vector () const;

   void Resize(const int& x_points, const int & y_points) ;


   // Print the whole content in grid ( for debugging purposes )
   void print();
   void write_plot( Grid u ,  const std::string& filename ) ;

private:

 std::vector<double>  vector_;
 int xSize_;
 int ySize_;

};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================






#endif //GRID_H

