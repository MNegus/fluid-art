#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "SGS/Vreman.h"
#include "tracer.h"
#include "diffusion.h"
#include "grid/adaptLES.h"


#define MAXLEVEL 5
#define temp 5

#define Re 1e3
#define sq(x)  x * x 

scalar T[];
scalar * tracers = {T};
mgstats mgT;
face vector av[];
double f = 1.39e-4;
int main() 
{
  	T[bottom]=dirichlet((265-(0.25*(t/3600))));
        T[top] = dirichlet (268);
        u.t[bottom]=dirichlet(0);
        u.t[top]=neumann(0);
        periodic (right);
        periodic (front);
        L0 = 400.;
        X0 = Y0 = Z0 = 0;
        DT = 10;
	init_grid (1 << MAXLEVEL);
        a=av;
        run();
}
face vector muv[];

event init(t=0)
{       
	adaptLESinit(MAXLEVEL); 
	mu=muv;
        foreach()
	{
                u.x[]= 8;
                T[]=((265+0.1*noise())*(y<50))+((y>50)*265) + (((0.01*(y-100)))*(y>100));
        }
        boundary(all);
}

event acceleration (i++)
{       
	face vector D = mu;
        mgT = diffusion (T, dt, D);
        foreach()
        {
          av.x[] = f*u.z[]-((y>300)*(u.x[]-8));
          av.y[] = 9.81*(T[]+T[0,-1,0])/(2*265)-((y>300)*u.y[]);
          av.z[] = f*(8-u.x[])-((y>300)*u.z[]);

        }
        boundary(all);
}

event SGS (i++)
{
        scalar muB[];
        eddyviscosity(0.17,u,muB);
        boundary({muB});
        foreach_face()
	{
	        muv.x[]=(muB[]+muB[-1,0,0])/2;
		
        }
        boundary((scalar *){muv});
}

event logfile (t <= temp; t += 1) 
{
	stats s = statsf (u.x);
	stats m = statsf(mu.x);
	fprintf (ferr, "t = %g\ti = %d\tdt = %g\tmin(Evis) = %g\nmax(Evis) = %g\tavg(Evis) = %g\tavg(ux) = %g \n",t, i, dt,m.min ,m.max,(float)m.sum/(L0*L0*L0), (float)s.sum/(L0*L0*L0));
}

event output (t += 1) 
{
	float Delt = L0/N;
	char name[80];
  	sprintf (name, "Tprof%d.dat",(int)roundf(t));
  	FILE * fp4 = fopen (name, "w");
  	sprintf (name, "uxprof%d.dat",(int)roundf(t));
  	FILE * fp5 = fopen (name, "w");
  	sprintf (name, "uzprof%d.dat",(int)roundf(t));
  	FILE * fp6 = fopen (name, "w");
	sprintf (name, "LEVELprof%d.dat",(int)roundf(t));
  	FILE * fp7 = fopen (name, "w");
	sprintf (name, "Yprof%d.dat",(int)roundf(t));
  	FILE * fp8 = fopen (name, "w");
	#ifdef _OPENMP
	#pragma omp parallel for reduction(+:F) default(shared)
	#endif
  	for (int k = 0; k < N; k++)
    	{
      		double uz=0 , ux=0 , F= 0, LEVEL=0, yz=0; ;
      		float yp = Delt * k + Y0 + Delt/2;
  		for (int i = 0; i < N; i++)
    		{
    			float xp = Delt*i + X0 + Delt/2.;
    			for (int j = 0; j < N; j++)
      			{
        			float zp = Delt*j + X0 + Delt/2.;
			        Point point = locate (xp, yp,zp);
        			F += T[];
        			uz+= u.z[];
        			ux+= u.x[];
				LEVEL += level;
				yz += y; 
			}
    		}
  		fprintf(fp4,"%g\t%g \n",yp,F/(N*N));
  		fprintf(fp5,"%g\t%g \n",yp,ux/(N*N));
  		fprintf(fp6,"%g\t%g \n",yp,uz/(N*N));
		fprintf(fp7,"%g\t%g \n",yp,LEVEL/(N*N));
		fprintf(fp8,"%g\t%g \n",yp,yz/(N*N));
	}
fclose(fp4);
fclose(fp5);
fclose(fp6);
fclose(fp7);
fclose(fp8);
}

void levels(double YY[],int LL[],int maxlevel)
{
	int DD = pow(2,maxlevel), l = 0;
		
	
	float xp = Z0+(L0/(DD*2)), zp = Z0 + (L0/(DD*2));  
	for (int k=0;k<DD;k++)
	{	
		float yp = Y0+(L0/(DD*2))+k*(L0/DD);	
		 Point point = locate (xp, yp,zp);
		if (k==0)
		{	
			YY[l]=y;
			LL[l]=level;
			l++;
		}	
		else if (y!=YY[l-1] )
		{
			YY[l]=y;
			LL[l]=level;
			l++; 
		}
	}
  
}

event adapt (t+=1)
{	
	double 	Yz[1<<MAXLEVEL] = {0};
	int	Lz[1<<MAXLEVEL] = {0};
	levels(Yz,Lz,MAXLEVEL);
	int 	m=0 , n=0; 
	

		
  
	while(Yz[m]>0.01) 
	{
	  	fprintf(ferr,"\t%d\t%g\t%d \n ",m+1,Yz[m] ,Lz[m]);	
	  	n+=pow(2,Lz[m]*2); 
	  	m++;
    	}
	fprintf(ferr,"\t#gridcells ~ %g ^ 3\n",round(pow(n,(0.33333))));
	adaptLESrun(MAXLEVEL,u,T,mu.y,Lz,Yz);
	boundary(all);
}
