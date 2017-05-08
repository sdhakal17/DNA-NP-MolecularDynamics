/* 
   Parallel Molecular Dynamics simulation of a DNA-Guided Nanoparticle Self-assembly

   Subas Dhakal

   Written for the Paper 
   Growth Dynamics for DNA-Guided Nanoparticle Crystallization[ACS Nano, 2013, 7 (12), pp 10948–10959]
   http://pubs.acs.org/doi/abs/10.1021/nn404476f

  Compile using 

  icc -xhost -openmp -O2 -o driver dna-np.cpp
  export OMP_NUM_THREADS=8
  ./driver

   Department of Material Science and Engineering, Northwestern University, Evanston, IL 60208
*/

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <vector>
# include <string>
# include <sstream>
# include <fstream>
# include <omp.h>
using std::vector;
using std::string;
using namespace std;
using std::ifstream;
using std::ofstream;

void  read_pot_force(double myr[],double p1[],double f1[],double p2[],double f2[],double p3[],double f3[])
{
  	int count = 0;
        string line;
	string filename="pot.txt";
	ifstream read_input(filename.c_str());	
	if(!read_input){
	}
	else{
		while (std::getline(read_input, line)){
		double x1,x2,x3,x4,x5,x6,x7;
		std::istringstream iss(line);		
		if(iss>>x1>>x2>>x3>>x4>>x5>>x6>>x7){
		myr[count]=x1;
		p1[count]=x2;
		f1[count]=x3;
		p2[count]=x4;
		f2[count]=x5;
		p3[count]=x6;
		f3[count]=x7;
		count = count + 1; }}
	
		read_input.close();		
       	}	

}
 	

void compute( int np, int nd, double pos[], double vel[],double mass,double f[],double box[],double *pot, double f1[],double p1[])
{
  double d;
  double d2;
  int i;
  int j;
  int k;
  double ke;
  double pe;
  double PI2 = 3.141592653589793 / 2.0;
  double rcut = 1.2000;
  double rij[3];

  pe = 0.0;
  ke = 0.0;
   for(k = 0; k < np; k++ ){
   f[k*nd] = 0.0;
   f[k*nd+1] = 0.0;
   f[k*nd+2] = 0.0;}

# pragma omp parallel shared(f,nd,np,pos,vel,box,f1,p1) private ( i, j, k, rij, d, d2 )
# pragma omp for 
  for(k = 0; k < np-1; k++ ){
    for ( j = k+1; j < np; j++ ){     
         double d,dx,dy,dz;
	dx = pos[k*nd]-pos[j*nd];
  	dy = pos[k*nd+1]-pos[j*nd+1];
  	dz = pos[k*nd+2]-pos[j*nd+2];

  	if (dx<0.0)  dx+=box[0]; 
  	if (dx>box[0])    dx-=box[0]; 
  	if (dy<0.0)  dy+=box[1]; 
  	if (dy>box[1])    dy-=box[1]; 
  	if (dz<0.0)  dz+=box[2]; 
  	if (dz>box[2])    dz-=box[2]; 

  	d = sqrt(dx*dx + dy*dy + dz*dz);
  	
	if(d<rcut){

        int myindex;
	myindex =  int ((d-0.1)/0.0002);
	//if(type[k]==1 && type[k]==1 || type[k]==1 && type[k]==1){
        pe = pe + 0.5 * p1[myindex];
        cout<<pe<<"	"<<1111111<<"\n";
	f[k*nd] = f[k*nd] + dx*f1[myindex]/d;
	f[k*nd+1] = f[k*nd+1] + dy*f1[myindex]/d;
	f[k*nd+2] = f[k*nd+2] + dz*f1[myindex]/d;
	f[j*nd] = f[j*nd] - dx*f1[myindex]/d;
	f[j*nd+1] = f[j*nd+1] - dy*f1[myindex]/d;
	f[j*nd+2] = f[j*nd+2] - dz*f1[myindex]/d;}
	//}     
      }
    }

   
  *pot = pe;

  return;
}

void initialize(int np,int nd,double box[],double pos[],double vel[],double force[],int particle_type[])
{
  	int count = 0,count1=0;
	string filename="myposition.xyzv";
	string line;
	ifstream read_input(filename.c_str());	
	if(!read_input){
	}
	else{
		while (std::getline(read_input, line)){
		double x1,x2,x3;
		string s1;
		std::istringstream iss(line);		
		if(iss>>s1>>x1>>x2>>x3){
		if(s1=="Na"){
		particle_type[count1] = 1;}
		else{
		particle_type[count1] = 2;}
		pos[count] = x1+0.5*box[0];
		pos[count+1] = x2+0.5*box[0];
		pos[count+2] = x3+0.5*box[0];
		count = count + 3;
		count1 = count1 + 1;
		}}	
		read_input.close();		
       	}

  	int i;
  	int j;

 	for(j = 0; j < np; j++ ){
    	for(i = 0; i < nd; i++ ){
      	vel[i+j*nd] = 0.2;    	
      	force[i+j*nd] = 1.0;
    	}}

  return;
}



void timestamp ( )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

void update_first( int np, int nd, double pos[], double vel[], double box[],double dt_2)
{
  int i;
  int j;
 # pragma omp parallel  shared (np,nd,pos,vel,box,dt_2)  private(i,j)
 # pragma omp for
  for(j = 0; j < np; j++ ){
       pos[j*nd] = pos[j*nd] + vel[j*nd] * dt_2;
       pos[j*nd+1] = pos[j*nd+1] + vel[j*nd+2] * dt_2;
       pos[j*nd+2] = pos[j*nd+2] + vel[j*nd+2] * dt_2;

      if (pos[j*nd]<0.0) { pos[j*nd]+=box[0]; }
      if (pos[j*nd]>box[0])   { pos[j*nd]-=box[0]; }
      if (pos[j*nd+1]<0.0) { pos[j*nd+1]+=box[1]; }
      if (pos[j*nd+1]>box[1])   { pos[j*nd+1]-=box[1]; }
      if (pos[j*nd+2]<0.0) { pos[j*nd+2]+=box[2]; }
      if (pos[j*nd+2]>box[2])   { pos[j*nd+2]-=box[2]; }   
  }
  return;
}

void update_second(int np,int nd,double pos[],double vel[],double force[],double dt_2,double *ke)
{
  double kinetic  = 0.0;
  int i;
  int j;
 # pragma omp parallel  shared (np,nd,pos,vel,force,dt_2,ke)  private(i,j)
 # pragma omp for
  for(j = 0; j < np; j++ ){

       vel[j*nd] = vel[j*nd]+dt_2*force[j*nd];
       vel[j*nd+1] = vel[j*nd+1]+dt_2*force[j*nd+1];
       vel[j*nd+2] = vel[j*nd+2]+dt_2*force[j*nd+2];

       pos[j*nd] = pos[j*nd] + vel[j*nd] * dt_2;
       pos[j*nd+1] = pos[j*nd+1] + vel[j*nd+2] * dt_2;
       pos[j*nd+2] = pos[j*nd+2] + vel[j*nd+2] * dt_2;  
       kinetic+=(vel[j*nd] * vel[j*nd]+ vel[j*nd+1] * vel[j*nd+1] + vel[j*nd+2] * vel[j*nd+2]);

	}

  *ke = kinetic;
  
  return;
}

void chain(double *KE,double dt,double dt_2,double dt_4,double dt_8,double *Q,double *xi,double *vxi,double *vel,int nl,int np,int nd,double T) 
{
  double G1, G2, s;
  int i;
  G2= (Q[0]*vxi[0]*vxi[0]-T);
  vxi[1]+=G2*dt_4;
  vxi[0]*=exp(-vxi[1]*dt_8);
  G1=(2*(*KE)-3*np*T)/Q[0];
  vxi[0]+=G1*dt_4;
  vxi[0]*=exp(-vxi[1]*dt_8);
  xi[0]+=vxi[0]*dt_2;
  xi[1]+=vxi[1]*dt_2;
  s=exp(-vxi[0]*dt_2);
  for (i=0;i<np;i++) {
    vel[i*nd]*=s; vel[1+i*nd]*=s; vel[2+i*nd]*=s;
  }
  (*KE)*=(s*s);
  vxi[0]*=exp(-vxi[1]*dt_8);
  G1=(2*(*KE)-3*np*T)/Q[0];
  vxi[0]+=G1*dt_4;
  vxi[0]*=exp(-vxi[1]*dt_8);
  G2=(Q[0]*vxi[0]*vxi[0]-T)/Q[1];
  vxi[1]+=G2*dt_4;
}

int main ( int argc, char *argv[] )
{

  double *force;
  double *box;
  double temp =0.95,tj_temp = 0.95;
  double dt = 0.001,dt2,dt_2,dt_4,dt_8;
  double e0;
  int i;
  int id;
  double PE,KE;
  double mass = 1.0;
  int nd = 3;
  int np = 8000;
  int NBIN =5532;
  double *pos;
  double potential;
  int seed = 123456789;
  int step;
  int step_num = 500;
  int step_print;
  int step_print_index;
  int step_print_num;
  double *vel;
  int *particle_type;
  double *myr,*p1,*p2,*p3,*f1,*f2,*f3;

  double wtime;

  timestamp ( );

  box = new double[nd];
  force = new double[nd*np];
  pos = new double[nd*np];
  vel = new double[nd*np];
  particle_type = new int[np];
  myr = new double[NBIN];
  p1 = new double[NBIN];
  f1 = new double[NBIN];
  p2 = new double[NBIN];
  f2 = new double[NBIN];
  p3 = new double[NBIN];
  f3 = new double[NBIN];
  read_pot_force(myr,p1,f1,p2,f2,p3,f3);

   /* Allocate arrays for the NH Chain */
  double *xi,*vxi,*Q;
  int nl = 2;
  Q = new double[nl];
  xi = new double[nl];
  vxi = new double[nl];
  Q[0] = Q[1] = 0.1;
  xi[0] = xi[1] = 0.0;
  vxi[0] = vxi[1] = 0.0;
  dt2  = dt*dt;
  dt_2 = 0.5*dt;
  dt_4 = 0.5*dt_2;
  dt_8 = 0.5*dt_4; 

  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  box[0] = box[1] = box[2] = 37.0;
 
  initialize(np,nd,box,pos,vel,force,particle_type);
  compute(np,nd,pos,vel,mass,force,box,&PE,f3,p3);

  cout <<"step"<< "  "<<"length"<<"	"<<"Temp"<<"	"<<"PE"<< "  " <<"KE"<< "  " <<"TotalE"<< "\n";
 
 
  step = 0;
  
  wtime = omp_get_wtime ( );

  for ( step = 1; step <= step_num; step++ )
  {
   
   chain(&KE,dt,dt_2,dt_4,dt_8,Q,xi,vxi,vel,nl,np,nd,temp);
    
    /* First integration half-step */
    KE = 0.0;

    update_first(np,nd,pos,vel,box,dt_2);

    compute(np,nd,pos,vel,mass,force,box,&PE,f3,p3);

    update_second(np,nd,pos,vel,force,dt_2,&KE);

   
    KE*=0.5;

    chain(&KE,dt,dt_2,dt_4,dt_8,Q,xi,vxi,vel,nl,np,nd,temp);
    if ( step %100==0)
    {
      /* do a temperature jump at the prescribed time */
    temp = tj_temp;	
    cout <<step<<"	"<<box[0]<<"	"<<KE*2/3./np<<"	"<<PE<<"	"<<KE<<"	"<<KE+PE<< "\n";
    ofstream subas("Mytrajectory.xyz",ios::app);
    subas<<np<<"\n"<<"Atoms. Timestep: "<<step<<"\n";
    for(int idx=0;idx<np;idx++){
    subas<<"HE"<<"	"<<pos[idx*nd]<<"	"<<pos[idx*nd+1]<<"	"<<pos[idx*nd+2]<<"\n";}}

  }

  	
  
  wtime = omp_get_wtime ( ) - wtime;
  cout << "\n";
  cout << "  Elapsed cpu time for main computation:\n";
  cout << "  " << wtime << " seconds.\n";

  delete [] box;
  delete [] force;
  delete [] pos;
  delete [] vel;

  timestamp ( );

  return 0;
}
