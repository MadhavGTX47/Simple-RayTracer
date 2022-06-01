# Simple-RayTracer
A Simple Ray Tracer using C++, Eigen with Reflections, shadows, anti-aliasing by stratified sampling.  

How to run:-

First run the make shade command to generate shade
then use any flags to generate desired image

Basic Assignment:-
 
Use:	```make shade ```
	
 to generate out file 
	the output is "shade"

Use:	```shade teapot-3.nff```
	for basic output

 the output is "output.ppm"

Flags:-

Use:	```shade -p teapot-3.nff```
	for Phong- Smooth shading     for extra credit smoothshading

 the output is "output.ppm"

Use:	```shade -j -s 3 teapot-3.nff ```
	for stratified sampling       for extra credit sampling
	
Here,-s 3 is the number of samples and can be used to specify the number of samples

the output is "output.ppm"

 
Resources used:-

Fundamentals of Computer Graphics, Fourth Edition
Pg 47
Pg 81-87
Pg 165
Pg 235-239
Pg 330-331
