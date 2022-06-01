#ifndef stuff_H
#define stuff_H
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>


using namespace std;
using namespace Eigen;

class View{
	public:
		double b[3];
		double f[8];
		RowVector3d from;
		RowVector3d up;
		RowVector3d at;
		double angle;
		double hither;
		int resx,resy;
		

};
class HitRecord {
	public:
		HitRecord();
		HitRecord(double intialize);
		Matrix<double,9, 3> trio;
		double t;
		bool intersected;
		RowVector3d intersection;
		RowVector3d normal;
};
int raytracer();
RowVector3d findcolor(RowVector3d rayd, RowVector3d  origin, double t0 ,double t1, int depth);

vector<Matrix<double,9, 3>> tri;
vector<Matrix3d> spr;
vector<RowVector3d> light;
Matrix3d temp;
Matrix3d ret;

#endif
