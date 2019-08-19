/*

Numerical solution to the systems of ODEs by
fifth-order fixed step RK method.

Mai Quang Vinh ;)

*/

#include <cmath> //Header file for mathematical operations
#include <iomanip> //Header file for precession
#include <iostream> //Header file for cin & cout
#include <fstream> //Header file for writing the results 

#include "model.h" //Header file for model


using namespace std; //Calling the standard directory


//***************************
/* Main program */

int main ()
{
	/* Type variables */
  	/* ODE integration parameters */
  	
	int i;
	int nsteps = 10; //Number of integration steps
	double h; // Integration step size 
	
	double t0 = 0.0;
	double tf = 1.0;

	/* Initial conditions */
	double y0 [neqn]= {0.0, -1.0};

	/* Open file */
	ofstream myfile("results.txt");

	if(myfile.is_open())
	{

	myfile << "t" << "\t" << "y1" 
		   << "\t" << "y2" << endl;
		   
 	myfile << fixed << setprecision(4) << t0 << "\t"
 	 	   << y0[0] << "\t" << y0[1] << endl;	  	
	
	/* Step size */
	h = (tf - t0)/(float)nsteps;
		
	/* nsteps rkc4 steps */
	for(i = 1; i <= nsteps; i++)
	{
	 	/* Type variables */
    	double ytb[neqn], yb[neqn], tb;
    	double y4[neqn], y5[neqn], ye[neqn];
    	
    	double k1[neqn], k2[neqn], k3[neqn];
    	double k4[neqn], k5[neqn], k6[neqn];
    	int j;
			
		/* Derivative vector at initial (base) point */
		
		for(j = 0; j < neqn; j++)
		{
			yb[j] = y0[j];
		}
		//yb;
		tb = t0;
		
		f(ytb, yb, tb);
		/* k1; stepping for k2 */
		for(j = 0; j < neqn; j++)
		{
			k1[j] = h*ytb[j];
			
		}
		
		for(j = 0; j < neqn; j ++)
		{
			yb[j] = y0[j] + (0.25)*k1[j];
		}	
		
		f(ytb, yb, tb + (0.25)*h);
			   
		/* k2; stepping for k3 */
		for(j = 1; j <= neqn; j++)
		{
			k2[j] = h*ytb[j];
		}
		
		for(j = 0; j< neqn; j++)
		{
			yb[j] = y0[j] + (3.0/32.0)*k1[j]
    			  	      + (9.0/32.0)*k2[j]; 
		}
		
		f(ytb, yb, tb + (3.0/8.0)*h);

		/* k3; stepping for k4 */
		for(j = 0; j < neqn; j++)
		{
	    	k3[j] = h*ytb[j];
		}
		
		for(j = 0; j < neqn; j++)
		{
			yb[j] = y0[j] + (1932.0/2197.0)*k1[j]
    			  		  - (7200.0/2197.0)*k2[j]
    			  		  + (7296.0/2197.0)*k3[j];
		}
		
		f(ytb, yb, tb + (12.0/13.0)*h);

		/* k4; stepping for k5 */
	 
		for(j = 0; j < neqn; j++)
		{
		    k4[j] = h*ytb[j];
		}
		
		for(j = 0; j < neqn; j++)
		{
			yb[j] = y0[j] + ( 439.0/216.0)*k1[j]
		    	   		  - (         8.0)*k2[j]
		    	   		  + (3680.0/513.0)*k3[j]
		           		  - (845.0/4104.0)*k4[j];
		}
		
		f(ytb, yb, tb + h);
		
		/* k5; stepping for k6*/
		
		for(j = 0; j < neqn; j++)
		{
		    k5[j] = h*ytb[j];

		}
		
		for(j = 0; j < neqn; j++)
		{
			yb[j] = y0[j] - (     8.0/27.0)*k1[j]
		    	   		  + (          2.0)*k2[j]
		    	   		  - (3544.0/2565.0)*k3[j]
		           		  + (1859.0/4104.0)*k4[j]
		           		  - (    11.0/40.0)*k5[j];
		}
		
		f(ytb, yb, tb + (0.5)*h);
		
		/* k6 */
		for(j = 0; j < neqn; j++)
		{
		    k6[j] = h*ytb[j];
		}    
		
		    
	    /* Fourth-order RK method*/
	    for(j = 0; j < neqn; j++)
	    {
	    	y4[j] = y0[j] + (   25.0/216.0)*k1[j]
	    				  + (1408.0/2565.0)*k3[j]
	            		  + (2197.0/4101.0)*k4[j]
	            		  - (      1.0/5.0)*k5[j];
		}
		
	    /* Fifth-order RK method */
		for(j = 0; j < neqn; j++)
		{
		    y5[j] = y0[j] + (     16.0/135.0)*k1[j]
		    			  + ( 6656.0/12825.0)*k3[j]
		            	  + (28561.0/56430.0)*k4[j]
		            	  - (       9.0/50.0)*k5[j]
		            	  + (       2.0/55.0)*k6[j];
		}
		
		for(j = 0; j < neqn; j++)
		{
			y0[j] = y5[j];
		}
		    
		    //cout << y0 << endl;
		    //e = y0 - ye; 

		t0 = tb + h;
			
		//ye[0] = t0;
		//ye[1] = t0 - 1.0;
		
		//cout << setprecision(4) << fixed;
	    //cout << y0[0] << "\t" << y0[1] << "\t"
	    //	 << ye[0] << "\t" << ye[1] << endl; 	
		
		myfile << fixed << setprecision(4) << t0
			   << "\t"  << y0[0] << "\t" << y0[1] << endl;	
		
		/* End of if */
		}
 
	/* Close file */
	myfile.close();

	/* End of if(myfile.is_open()) */
	}

/* End of main */
return 0;
}

		
	
	
