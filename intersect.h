#ifndef intersect_H
#define intersect_H
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include "stuff.h"

using namespace std;
using namespace Eigen;



Matrix3d intersect_poly(RowVector3d raydirection, RowVector3d  origin,double t0 ,double t1,Matrix<double,9, 3> tria,HitRecord &hit)
{
	ret<< 0,0,0,0,0,0,0,0,0;	

	RowVector3d a = tria.row(0);
	RowVector3d b = tria.row(2);
	RowVector3d c = tria.row(4);


	double	xab = a[0] - b[0];
	double	xac = a[0] - c[0];
	double	x = raydirection[0];
	double	yab = a[1] - b[1];
	double	yac = a[1] - c[1];
	double	y = raydirection[1];
	double	zab = a[2] - b[2];
	double	zac = a[2] - c[2];
	double	z = raydirection[2];
	
		double m = (xab*((yac*z)-(y*zac))) + (yab*((x*zac)-(xac*z))) + (zab*((xac*y)-(yac*x)));

	double	xd = a[0] - origin[0];
	double	yd = a[1] - origin[1];
	double	zd = a[2] - origin[2];

		double gamma = ((z*((xab*yd)-(xd*yab))) + (y*((xd*zab)-(xab*zd))) + (x*((yab*zd)-(yd*zab)))) / m;
		    if(gamma < 0 || gamma > 1)
		    return ret;
		

		double beta = ((xd*((yac*z)-(y*zac))) + (yd*((x*zac)-(xac*z))) + (zd*((xac*y)-(yac*x))))/m;
		    if(beta < 0 || beta > (1 - gamma))
		    return ret;
		

		double t = -(((zac*((xab*yd)-(xd*yab))) + (yac*((xd*zab)-(xab*zd))) + (xac*((yab*zd)-(yd*zab)))))/m;
	
		 

	  	if(t > t0 && t < hit.t){
		hit.t = t;
		hit.intersected = true;
		hit.trio =tria;
		hit.intersection = origin + raydirection*t;
		hit.normal = ((b-a).cross(c-a)).normalized();
	
	}
      return temp;          
}

bool shade(RowVector3d raydirection, Matrix<double,9, 3> tria, RowVector3d  origin)
{
	

	RowVector3d a = tria.row(0);
	RowVector3d b = tria.row(2);
	RowVector3d c = tria.row(4);


	double	xab = a[0] - b[0];
	double	xac = a[0] - c[0];
	double	x = raydirection[0];
	double	yab = a[1] - b[1];
	double	yac = a[1] - c[1];
	double	y = raydirection[1];
	double	zab = a[2] - b[2];
	double	zac = a[2] - c[2];
	double	z = raydirection[2];
	
		double m = (xab*((yac*z)-(y*zac))) + (yab*((x*zac)-(xac*z))) + (zab*((xac*y)-(yac*x)));

	double	xd = a[0] - origin[0];
	double	yd = a[1] - origin[1];
	double	zd = a[2] - origin[2];

		double gamma = ((z*((xab*yd)-(xd*yab))) + (y*((xd*zab)-(xab*zd))) + (x*((yab*zd)-(yd*zab)))) / m;
		    if(gamma < 0 || gamma > 1)
		    return false;
		

		double beta = ((xd*((yac*z)-(y*zac))) + (yd*((x*zac)-(xac*z))) + (zd*((xac*y)-(yac*x))))/m;
		    if(beta < 0 || beta > (1 - gamma))
		    return false;
		

		double t = -(((zac*((xab*yd)-(xd*yab))) + (yac*((xd*zab)-(xab*zd))) + (xac*((yab*zd)-(yd*zab)))))/m;
		    
		if( t > 0.000000001)// comapring to a very small value to avoid acne
		{
		return true;
		}
		else
		{
		return false;
		}    

}


Matrix3d  intersect_sphere(RowVector3d raydirection, RowVector3d  origin,double t0 ,double t1,Matrix3d sphe)
{
    
			ret <<0,0,0,0,0,0,0,0,0;
			RowVector3d center= sphe.row(0);
			RowVector3d ec = origin-center;
			double A = raydirection.dot(raydirection);
			double B = 2.0*ec.dot(raydirection);
			double C = (ec.dot(ec))-(sphe(1,0)*sphe(1,0));
			
			if(B*B-4*A*C>=0)
			{
			

			
			double ts=(((-raydirection).dot(ec))+sqrt(((raydirection.dot(ec)*(raydirection.dot(ec))-((raydirection.dot(raydirection))*((ec.dot(ec))-sphe(1,0)*sphe(1,0)))))))/A;
			double ts2=(((-raydirection).dot(ec))-sqrt(((raydirection.dot(ec)*(raydirection.dot(ec))-((raydirection.dot(raydirection))*((ec.dot(ec))-sphe(1,0)*sphe(1,0)))))))/A;  
			

			if(ts < t0 || ts > t1	|| ts2 < t0 || ts2 > t1)
		    return ret;

			 if(ts2<ts)
			   {
				ts=ts2;
			   }


			   temp << sphe(2,0),sphe(2,1),sphe(2,2),ts,0,0,0,0,0;
      			return temp;      
			}
			
			else
			{ 
			return ret;	
			}
}		
#endif
