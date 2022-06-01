#include <string>
#include <unistd.h>
#include <time.h>
#include <ctype.h>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <stdio.h>
#include <vector>
#include <math.h>       
#include <Eigen/Geometry>
#include "stuff.h"
#include "intersect.h"
using namespace Eigen;
using namespace std;

#define PI 3.14159265

View view;									     //View Object which stores our vectors and other important data
bool smooth=false;	
double n=1;
bool jitter=false;
double aperture = 0;									 

int main(int argc, char* argv[])
{	 

 int c;
while ((c = getopt (argc, argv, "apjs:")) != -1)  //Handling command line flags
    switch (c)
      {
      case 'a':
       aperture = atof(optarg);
        break;
      case 'p':
        smooth = true;
        break;
	  case 'j':
        jitter = true;
        break;	
      case 's':
        n = atof(optarg);
        break;
      case '?':
        if (optopt == 'c')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }

	
	if(smooth==true)
	{
	cout<<"Generating the smooth output.ppm";
	}

	else{
	if(argv[1])
	  cout<<"Generating the output.ppm";		//To take the Input file name
	else
	  cout << "Give proper arguments" << endl;
	}

    Matrix<double,9, 3> triangle; 		  		 // Temporary variable to store Traingle data while entering into vector
	Matrix3d sphere;							 // Temporary variable to store Sphere data while entering into vector
	RowVector3d l;								 // Temporary variable to store light data while entering into vector			
    fstream nff(argv[optind++]);           		 // Making an fstream for the nff file 
	string line;               			  		 // to get each line from nff file
	vector<string> tokens;     			  		 // vector to store the tokens
	string temp;                          		 // temporay string used in tokenization
	
	
	
	while (std::getline(nff, line)){        	 //to read the nff file line by line
	
		istringstream iss(line);
		stringstream li(line);          		 // storing the line in li of stringstream
	
	while(getline(li, temp, ' ')){       		 // tokenization with white spaces as delimter
        tokens.push_back(temp);           		 // storing the tokens in the token vector
    }

	}

	
	for(int i=0; i < tokens.size(); i++){  		//Iterating through tokens
	   
	   
		if(tokens[i]=="b")					    //Storing the background color
		{
			view.b[0]=stod(tokens[++i]);
			view.b[1]=stod(tokens[++i]);
			view.b[2]=stod(tokens[++i]);
		}
		
		if(tokens[i]=="v")					    //Storing the view
		{
			i++;			
			if(tokens[i]=="from")				//Storing the from coordinate
			{
			view.from<<stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]); 
			}
			i++;
			if(tokens[i]=="at")					//Storing the At coordinate
			{
			view.at<<stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]);  
			}
			i++;
			if(tokens[i]=="up")					//Storing the Up coordinate
			{
			view.up<<stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]);  
			}
			i++;
			if(tokens[i]=="angle")				//Storing the Angle
			{
			view.angle=stod(tokens[++i]);
			}
			i++;
			if(tokens[i]=="hither")				//Storing the Hither
			{
			tokens[++i]; 
			}
			i++;
			if(tokens[i]=="resolution")			//Storing the Resolution
			{
			view.resx=stoi(tokens[++i]); view.resy=stoi(tokens[++i]);
			}
		
		}
		
		if(tokens[i]=="l")						//Storing the lights
			{
			l<< stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]);
			light.push_back(l);
			}	



		if(tokens[i]=="f")					    //Storing the fill color
		{
			view.f[0]=stod(tokens[++i]);
			view.f[1]=stod(tokens[++i]);
			view.f[2]=stod(tokens[++i]);
			view.f[3]=stod(tokens[++i]);
			view.f[4]=stod(tokens[++i]);
			view.f[5]=stod(tokens[++i]);
			view.f[6]=stod(tokens[++i]);
			view.f[7]=stod(tokens[++i]);

		}
		
		
		if(tokens[i]=="p")					    //Storing the polygons
		{	
			double ver;							//temp vector
			int v = stoi(tokens[++i]) ;         //Converting token which is in string to integer for comparison
			vector<float> tr;

			if(v>3)								//To handle polygons with more than 3 vertices
			{ 
				for(int tx=0;tx<v*3;tx++)		//Storing all the vertices in temporary vector to split and store
				{ 
				ver = stod(tokens[++i]);
				tr.push_back(ver);
				}
								
				v=v-2;
				int t=0;
				for(double co=0;co < v;co++)
				{
				triangle << tr[0],tr[1],tr[2],0,0,0,tr[t+3],tr[t+4],tr[t+5],0,0,0,tr[t+6],tr[t+7],tr[t+8],0,0,0,view.f[0],view.f[1],view.f[2],view.f[3],view.f[4],view.f[5],view.f[6],view.f[7],1;
				tri.push_back(triangle);
				t=t+3;							//Here we saved our fill color in last row of the matrix
			    }
			}
			
			else
			{   								// To handle polygons with 3 vertices
				triangle << stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),0,0,0,stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),0,0,0,stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),0,0,0,view.f[0],view.f[1],view.f[2],view.f[3],view.f[4],view.f[5],view.f[6],view.f[7],0; 
				tri.push_back(triangle);		//Here we saved our fill color in last row of the matrix
				
			}
			
		}
	
		if(tokens[i]=="pp")					    //Storing the polygon patches
		{	
			double ver;							//temp vector
			int v = stoi(tokens[++i]) ;         //Converting token which is in string to integer for comparison
			vector<float> tr;

			if(v>3)								//To handle polygon patch with more than 3 vertices
			{ 
				cout<<"Currently not being handled";		//Storing all the vertices in temporary vector to split and store
				
			}
			
			else
			{   								// To handle polygons patches with 3 vertices
				triangle << stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),view.f[0],view.f[1],view.f[2],view.f[3],view.f[4],view.f[5],view.f[6],view.f[7],0; 
				tri.push_back(triangle);		//Here we saved our fill color in last row of the matrix
				
			}
			
		}	

		if(tokens[i]=="s")					 	//Storing the spheres
		{  
		 sphere<< stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),0,0,view.f[0],view.f[1],view.f[2];
         spr.push_back(sphere);					//Here we saved our fill color in last row of the matrix
		}
		
	}

	raytracer();      						   //Calling our raytracer as we have finished parsing

return 0;

}

HitRecord::HitRecord(double intitialT){
	t = intitialT;
	intersected = false;
}


int raytracer() {

unsigned char pixels[view.resx][view.resy][3];	  //Declalring pixels which is to store output
unsigned int hitrecord[view.resx][view.resy];	  //Declaring the hit record

RowVector3d distance =view.from-view.at;		  // Calculating Distance 
double	dist=sqrt((distance(0,0)*distance(0,0))+(distance(0,1)*distance(0,1))+(distance(0,2)*distance(0,2)));

RowVector3d w = (view.from-view.at).normalized(); // Calculating u,v,w
RowVector3d u = (view.up.cross(w)).normalized();
RowVector3d v = (w.cross(u)).normalized();


int p;
double l= (tan ((view.angle * PI / 180.0) /2))*dist; // Assuming that height is h ="1"
double wp= ((l*2)/view.resx);         			     // Calculating the pixel width
RowVector3d hit_intersection;
RowVector3d hit_normal;

srand( (unsigned)time( NULL ) );				 //For generating random values in 0,1

	for(int ix = 0; ix < view.resx; ix++)
	{  											 // calculating x,y of each pixel
		for(int iy = 0; iy <view.resy; iy++)
		{

		double x;
		double y;
		RowVector3d rgb={0,0,0};	
		RowVector3d rgbz={0,0,0};	
			

		if(jitter==true)								//if jitter
		{
			for(double p = 0;p< n;p++)
				{				 	
			for(double q = 0;q< n;q++)
				{
				float random=(float) rand()/RAND_MAX;	//Sending random rays(for stratified sampling ) where n is the number of samples
			
				 x= -l+wp/2+(ix + (p + random)/n)*wp;  
				 y= l-wp/2-(iy + (q + random)/n)*wp;  
				

				RowVector3d rayd = -w*dist+u*x+v*y;		 // Our Ray direction for the x,y

				rgbz = findcolor(rayd,view.from,0,0,0);
				rgb= rgb+rgbz;

				}
				}
			}
		else											//No jitter
		{

				 x= -l+wp/2+(ix)*wp;  
				 y= l-wp/2-(iy)*wp;  

				RowVector3d rayd = -w*dist+u*x+v*y;		 // Our Ray direction for the x,y

				rgbz = findcolor(rayd,view.from,0,0,0);
				rgb= rgb+rgbz;				
			}



		rgb=rgb/(n*n);
	
	

		if(rgb(0,0)>1) 									//Limit the colors to be in (0,1)
		{rgb(0,0)=1; }
		if(rgb(0,1)>1)
		{rgb(0,1)=1; }
		if(rgb(0,2)>1)
		{rgb(0,2)=1; }

		pixels[iy][ix][0]= rgb(0,0)*255;
		pixels[iy][ix][1]= rgb(0,1)*255;
		pixels[iy][ix][2]= rgb(0,2)*255;


		}
	}
	
		cout<<endl<<endl<<"Writing...";
		FILE *f = fopen("output.ppm","wb");		 //Saving the output to a ppm file
		fprintf(f, "P6\n%d %d\n%d\n", view.resx, view.resy, 255);
		fwrite(pixels, 1, view.resx*view.resy*3, f);
		fclose(f);	
		cout<<endl<<endl<<"done. ";
	
 return 0;
}		

RowVector3d findcolor(RowVector3d rayd, RowVector3d  origin, double t0 ,double t1, int depth)
{		double rgb[3]={0,0,0};
		RowVector3d colr;
		RowVector3d colr_temp;
		HitRecord hit = HitRecord(999);


		for(int i=0;i<tri.size(); i++){			 //Iterating through polygons
			
			Matrix3d tp = intersect_poly(rayd,origin,0.0000000001,999,tri[i],hit);
		}		
					
						
			 if(hit.intersected)
				{	rgb[0]=0;
					rgb[1]=0;
					rgb[2]=0;
					double rgbv1[3]={0,0,0};
					double rgbv2[3]={0,0,0};
					double rgbv3[3]={0,0,0};


					Matrix<double,9, 3> triang = hit.trio; 
					double intensity = 1/sqrt(light.size());
				
					
					for(int j = 0; j < light.size(); j++){  //looping through lights

					RowVector3d light_direction = light[j] - hit.intersection;
					light_direction= light_direction.normalized();
			
					RowVector3d view_direction = view.from - hit.intersection;
					view_direction= view_direction.normalized();

					RowVector3d vector_avg = view_direction + light_direction;
					vector_avg= vector_avg.normalized();	
					
				
					bool shadow = false;
					for(int s = 0; s < tri.size(); s++){	//Checking for Shadows
					if(shade(light_direction,tri[s],hit.intersection))	
					{	
						shadow= true;
					}

					}
					

					if(!shadow){							// Calculating color based on specular,shine,diffuse
					
					if(smooth==true && triang(8,2)==0 )		//smooth shading (here we are using last bit of triang to indicate if we have normal inforamtion or not)
					{
					RowVector3d a = triang.row(0);
					RowVector3d b = triang.row(2);
					RowVector3d c = triang.row(4);
					RowVector3d n1 = triang.row(1);
					RowVector3d n2 = triang.row(3);
					RowVector3d n3 = triang.row(5);
				

					double	xab = a[0] - b[0];
					double	xac = a[0] - c[0];
					double	x = rayd[0];
					double	yab = a[1] - b[1];
					double	yac = a[1] - c[1];
					double	y = rayd[1];
					double	zab = a[2] - b[2];
					double	zac = a[2] - c[2];
					double	z = rayd[2];					// Calculating alpha, beta, gamma for color interpolation
	
					double m = (xab*((yac*z)-(y*zac))) + (yab*((x*zac)-(xac*z))) + (zab*((xac*y)-(yac*x)));

					double	xd = a[0] - origin[0];
					double	yd = a[1] - origin[1];
					double	zd = a[2] - origin[2];

					double gamma = ((z*((xab*yd)-(xd*yab))) + (y*((xd*zab)-(xab*zd))) + (x*((yab*zd)-(yd*zab)))) / m;
		   
					double beta = ((xd*((yac*z)-(y*zac))) + (yd*((x*zac)-(xac*z))) + (zd*((xac*y)-(yac*x))))/m;
		    
		   
					double alpha=1-beta-gamma;
					RowVector3d bay(alpha,beta,gamma);
				
																//Calculating colors at vertices with their normal data

					double diffusev1 = triang(7,0)*max(0.0,(n1.dot(light_direction)));
					double specularv1 = max(0.0,(n1.dot(vector_avg)));
					specularv1 = triang(7,1)*pow(specularv1,triang(7,2));
					
					rgbv1[0] += triang(6,0)*diffusev1 * intensity;
					rgbv1[1] += triang(6,1)*diffusev1 * intensity;
					rgbv1[2] += triang(6,2)*diffusev1 * intensity;			

					rgbv1[0] += specularv1 * intensity;
					rgbv1[1] += specularv1 * intensity;
					rgbv1[2] += specularv1 * intensity;
					
					double diffusev2 = triang(7,0)*max(0.0,(n2.dot(light_direction)));
					double specularv2 = max(0.0,(n2.dot(vector_avg)));
					specularv2 = triang(7,1)*pow(specularv2,triang(7,2));
					
					rgbv2[0] += triang(6,0)*diffusev2 * intensity;
					rgbv2[1] += triang(6,1)*diffusev2 * intensity;
					rgbv2[2] += triang(6,2)*diffusev2 * intensity;			

					rgbv2[0] += specularv2 * intensity;
					rgbv2[1] += specularv2 * intensity;
					rgbv2[2] += specularv2 * intensity;
					
					double diffusev3 = triang(7,0)*max(0.0,(n3.dot(light_direction)));
					double specularv3 = max(0.0,(n3.dot(vector_avg)));
					specularv3 = triang(7,1)*pow(specularv3,triang(7,2));
					
					rgbv3[0] += triang(6,0)*diffusev3 * intensity;
					rgbv3[1] += triang(6,1)*diffusev3 * intensity;
					rgbv3[2] += triang(6,2)*diffusev3 * intensity;			

					rgbv3[0] += specularv3 * intensity;
					rgbv3[1] += specularv3 * intensity;
					rgbv3[2] += specularv3 * intensity;
					
					rgb[0]=bay(0,0)*rgbv1[0]+bay(0,1)*rgbv2[0]+bay(0,2)*rgbv3[0];
					rgb[1]=bay(0,0)*rgbv1[1]+bay(0,1)*rgbv2[1]+bay(0,2)*rgbv3[1];
					rgb[2]=bay(0,0)*rgbv1[2]+bay(0,1)*rgbv2[2]+bay(0,2)*rgbv3[2];	
												
					}



					else{
					double diffuse = triang(7,0)*max(0.0,(hit.normal.dot(light_direction)));
					double specular = max(0.0,(hit.normal.dot(vector_avg)));
					specular = triang(7,1)*pow(specular,triang(7,2));
					
					rgb[0] += triang(6,0)*diffuse * intensity;
					rgb[1] += triang(6,1)*diffuse * intensity;
					rgb[2] += triang(6,2)*diffuse * intensity;			

					rgb[0] += specular * intensity;
					rgb[1] += specular * intensity;
					rgb[2] += specular * intensity;

					
					}
					}

					}			
					
					if(depth < 5 && triang(7,1) > 0 ){      														//Calculating reflection 
							
					RowVector3d viewDirection = hit.intersection - origin;  										//Calculating View for reflection
					viewDirection=viewDirection.normalized();
					RowVector3d reflectionDirection = viewDirection - 2*(viewDirection.dot(hit.normal))*hit.normal; //Calculating reflection direction
					reflectionDirection =  reflectionDirection.normalized();

					colr_temp = findcolor(reflectionDirection, hit.intersection, 0.000000001, 999, depth + 1); 		//Recusrive call to calculate reflection 
				
															
					rgb[0] += triang(7,1)*colr_temp(0,0);
					rgb[1] += triang(7,1)*colr_temp(0,1);
					rgb[2] += triang(7,1)*colr_temp(0,2);
					
					}	

					colr<<rgb[0],rgb[1],rgb[2];
					return colr;
				}

				else{
					 																								//Returning Background color										
					rgb[0] += view.b[0];
					rgb[1] += view.b[1];
					rgb[2] += view.b[2];
					colr<<rgb[0],rgb[1],rgb[2];
					return colr;
				}	

    
					
}



		



