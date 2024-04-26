#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Simulation parameters
// Frej,PRB'2023

const double Ms=4*3.14159*572.957; //=7.2 Oe // Sat. Magn. [A/m] 
const double K=-1.0*250; //=-3.0e3 erg/cc Uniax. Magnetocryst. Anis. Const. [J/m3]
const double theta=0.0*M_PI/180.0; // Ext field polar angle [degrees]
const double phi=0.0*M_PI/180.0; //Ext field azimuthal angle [degrees]
const double hx=1.0/sqrt(2.0); //External field direction Ox
const double hy=1.0/sqrt(2.0); // Oy
const double hz=0.0; //Oz
const double Hk=2.0*K/Ms; // aprox. 700 Oe = 0.07 T . It already includes the (thin film) demag contribution; // Anis. field amplitude [Tesla]
const double K1=-1.0*840;  // Cubic Anis. K1 [J/m3]
const double K2=0.0; //Cubic Anis K2 [J/m3]
const double Hext=0.04; //External field amplitude [Tesla] 
const double gamma1 = 1.7595*1e11; //el. gyromagn. ratio [rad/sT]
const double alpha = 0.2; //damping
const double c1=gamma1/(1+alpha*alpha); //LLG prefactor
const double c2=c1*(alpha); //c1*(alpha/Ms); //LLG prefactor
const double ex=0.0; //Ox
const double ey=0.0; //Oy
const double ez=1.0; //Oz EA directions
const double dt=0.1e-12; // Integration time step
const double dto=1e-12; // Output time step
const double tmin=0; //Start simulation time
const double tmax=3e-9;// End simulation time
const double mx0=0.57735;//0.763089; // 111: 0.57735 everywhere Initial magn. conditions Ox
const double my0=0.57735;//0.598044; //Oy
const double mz0=0.57735;//0.245028; //Oz
const double tol=1e-12; //Tolerance for stop criterion
const double tau=20e-12; 
const double t0=0;
const double H2=0.04; 

// Func declaration
void derivs(const double t, std::vector<double> &y, std::vector<double> &dmdt, std::vector<double> &H);
void getH(double Hext, std::vector<double> &y, std::vector<double> &H, bool tfield, double t, double &Hl);
void rk4(std::vector<double> &m, std::vector<double> &dmdt, std::vector<double> &H, const double t,
const double h, std::vector<double> &mout, void derivs(const double, std::vector<double> &, std::vector<double> &, std::vector<double> &));
void heun(std::vector<double> &m, std::vector<double> &dmdt, std::vector<double> &H, const double t,
const double h, std::vector<double> &mout, void derivs(const double, std::vector<double> &, std::vector<double> &, std::vector<double> &));
void calctorque(std::vector<double> &y, std::vector<double> &H,std::vector<double> &torque);
void calcmfielddir(std::vector<double> &y, double &mpar);

int main(int argc, char const *argv[]) {

// Output file
  std::ofstream f1;
  f1.open("output.txt", std::ofstream::out);
  f1<<"time"<<" "<<"mx"<<" "<<"my"<<" "<<"mz"<<" "<<"Hx"<<" "<<"Hy"<<" "<<"Hz"<<" "<<"Hext"<<" "<<"torque modulus"<<" "<<"m_DOT_H"<<"\n";
// Vars
  std::vector<double>min,mout, H, torque, dmdt;
  double t, mpar;
  double Hl=0.0;
  int nsteps=(int)((tmax-tmin)/dt); // Calc. number of steps according to Start/End time and time-step
  int outstep=(int)(dto/dt); //Calc. file output step
// Set problem dims.
  const int dim=3;
  min.resize(dim,0.0);
  mout.resize(dim,0.0);
  H.resize(dim+1,0.0); // Total field Ox , Oy, Oz, modulus
  dmdt.resize(dim,0.0);
  torque.resize(dim+1,0.0); // torque Ox,Oy,Oz, modulus
  bool tfield = false;
// Set initial conditions, magnetisation:
  min[0]=mx0; //mx initial
  min[1]=my0; //my initial
  min[2]=mz0; //mz initial

// Equil LLG eq. : dm/dt
 for(int i=0; i<nsteps; i++)
  {
    t=i*dt; //Calc. real time step  
    getH(Hext,min,H, tfield, t,Hl); //Calc. total field
    rk4(min, dmdt, H, t, dt,mout,derivs); // Call Runge-Kutta 4th order
    for(int i=0; i<dim; i++)mout[i]/=sqrt(mout[0]*mout[0]+mout[1]*mout[1]+mout[2]*mout[2]); //Normalise m_{x,y,z}/|m| to maintain |m|=1
    for(int i=0; i<dim; i++)min[i]=mout[i]; //Set new initial cond.
//    calcmfielddir(mout, mpar);
    if(i%outstep==0)f1<<t<<" "<<mout[0]<<" "<<mout[1]<<" "<<mout[2]<<" "<<H[0]<<" "<<H[1]<<" "<<H[2]<<" "<<Hext<<" "<<torque[3]<<" "<<mpar<<" "<<Hl<<"\n"; //Print 
    calctorque(mout,H,torque); //Check torque condition and break if satisfied
    if(torque[3]<tol)break;  
  }

   double temp=t;
   tfield=true;
// Laser-induced dynamics LLG eq. : dm/dt
 for(int i=0; i<nsteps; i++)
  {
    t=i*dt; //Calc. real time step
    getH(Hext,min,H,tfield, t,Hl); //Calc. total field
    rk4(min, dmdt, H, t, dt,mout,derivs); // Call Runge-Kutta 4th order
    for(int i=0; i<dim; i++)mout[i]/=sqrt(mout[0]*mout[0]+mout[1]*mout[1]+mout[2]*mout[2]); //Normalise m_{x,y,z}/|m| to maintain |m|=1
    for(int i=0; i<dim; i++)min[i]=mout[i]; //Set new initial cond.
//    calcmfielddir(mout, mpar);
    if(i%outstep==0)f1<<temp+t<<" "<<mout[0]<<" "<<mout[1]<<" "<<mout[2]<<" "<<H[0]<<" "<<H[1]<<" "<<H[2]<<" "<<Hext<<" "<<torque[3]<<" "<<mpar<<" "<<Hl<<"\n"; //Print
    calctorque(mout,H,torque); //Check torque condition and break if satisfied
    if(torque[3]<tol)break;
  }

  return 0;
}
// Func definitions
void heun(std::vector<double> &m, std::vector<double> &dmdt, std::vector<double> &H, const double t,
const double h, std::vector<double> &mout, void derivs(const double, std::vector<double> &, std::vector<double> &, std::vector<double> &))
{//Simpler heun integrator


	int n=m.size(); //dimension of the problem. Here n=3 , for x,y,z.
	std::vector<double> dmt, dmm, mt;
	dmt.resize(n,0.0);
	mt.resize(n,0.0);

	for (int i=0; i<n; i++) mt[i]=m[i] + h*dmdt[i];
	derivs(t+h, mt, dmt, H); 
	for(int i=0; i<n; i++)mout[i]=m[i]+h*0.5*(dmdt[i]+dmt[i]);		                
}


void rk4(std::vector<double> &m, std::vector<double> &dmdt, std::vector<double> &H, const double t,
const double h, std::vector<double> &mout, void derivs(const double, std::vector<double> &, std::vector<double> &, std::vector<double> &))
{ //Runge-Kutta 4th order integration algorithm

	int n=m.size(); //dimension of the problem. Here n=3 , for x,y,z.
	std::vector<double> dmt, dmm, mt;
	dmt.resize(n,0.0);
	dmm.resize(n,0.0);
	mt.resize(n,0.0);

	double hh=h*0.5;
	double h6=h/6.0;
	double th=t+hh;

	for (int i=0; i<n; i++) mt[i]=m[i] + hh*dmdt[i]; //First step
	derivs(th, mt, dmt, H); //Second step - obtain k2
	for (int i=0; i<n; i++) mt[i]=m[i] + hh*dmt[i];
	derivs(th, mt, dmm, H); //Third step - obtain k3
	for (int i=0; i<n; i++) {
		mt[i]=m[i]+h*dmm[i];
		dmm[i]+=dmt[i];
	}

	derivs(t+h, mt, dmt, H); // Fourth step - obtain k4
	for(int i=0; i<n; i++)mout[i]=m[i]+h6*(dmdt[i]+dmt[i]+2.0*dmm[i]);
                
}


void getH(double Hext, std::vector<double> &y, std::vector<double> &H, bool tfield, double t, double &Hl)
{  // Func to calculate effective field
	H[0]= Hext*hx; //Hext*sin(theta)*cos(phi);
	H[1]= Hext*hy; //Hext*sin(theta)*sin(phi);
	H[2]= Hext*hz; //Hext*cos(theta);

// Uniax. anis.
	H[0]+= Hk*(y[0]*ex+y[1]*ey+y[2]*ez)*ex;
	H[1]+= Hk*(y[0]*ex+y[1]*ey+y[2]*ez)*ey;
	H[2]+= Hk*(y[0]*ex+y[1]*ey+y[2]*ez)*ez;

// Cubic anis.

	double D11=(-2.0/Ms)*(K1*(y[1]*y[1]+y[2]*y[2]) + K2*(y[1]*y[1]*y[2]*y[2]));
        double D22=(-2.0/Ms)*(K1*(y[0]*y[0]+y[2]*y[2]) + K2*(y[0]*y[0]*y[2]*y[2]));
	double D33=(-2.0/Ms)*(K1*(y[0]*y[0]+y[1]*y[1]) + K2*(y[0]*y[0]*y[1]*y[1]));

	H[0]+= D11*y[0];
	H[1]+= D22*y[1];
	H[2]+= D33*y[2];

// Laser field

	if(tfield && t<=tau)
//	if(tfield)
	{	Hl=H2; // To use for squared pulse
		//Hl=H2*exp(-((t-t0)/tau)*((t-t0)/tau));
		H[0]+= Hl;
		H[1]+= 0.0 ;
		H[2]+= 0.0 ;
	}
	else
	{
		Hl=0;
	}

// Calculate Effective field modulus	
        H[3]=sqrt(H[0]*H[0]+H[1]*H[1]+H[2]*H[2]);	




}

void derivs(const double t, std::vector<double> &y, std::vector<double> &dmdt,std::vector<double> &H)
{ // Func to calculate the rhs of the LLG equation dm/dt

  // LLG eq. along Ox, Oy and Oz
  dmdt[0] = c1*(H[1]*y[2] - H[2]*y[1]) + c2*( H[0]*(y[1]*y[1] + y[2]*y[2]) - y[0]*(H[1]*y[1] + H[2]*y[2]));
  dmdt[1] = c1*(H[2]*y[0] - H[0]*y[2]) + c2*( H[1]*(y[0]*y[0] + y[2]*y[2]) - y[1]*(H[0]*y[0] + H[2]*y[2]));
  dmdt[2] = c1*(H[0]*y[1] - H[1]*y[0]) + c2*( H[2]*(y[0]*y[0] + y[1]*y[1]) - y[2]*(H[0]*y[0] + H[1]*y[1]));

}
void calctorque(std::vector<double> &y, std::vector<double> &H,std::vector<double> &torque)
{
  // Calculate the torque
  torque[0]=(c1*(H[1]*y[2] - H[2]*y[1]) + c2*( H[0]*(y[1]*y[1] + y[2]*y[2]) - y[0]*(H[1]*y[1] + H[2]*y[2])));
  torque[1]=(c1*(H[2]*y[0] - H[0]*y[2]) + c2*( H[1]*(y[0]*y[0] + y[2]*y[2]) - y[1]*(H[0]*y[0] + H[2]*y[2])));
  torque[2]=(c1*(H[0]*y[1] - H[1]*y[0]) + c2*( H[2]*(y[0]*y[0] + y[1]*y[1]) - y[2]*(H[0]*y[0] + H[1]*y[1])));
  torque[3]=sqrt(torque[0]*torque[0]+torque[1]*torque[1]+torque[2]*torque[2]);
}
void calcmfielddir(std::vector<double> &y, double &mpar)
{
 //Calculate the projection of m on the dir. of the field
  mpar=(y[0]*sin(theta)*cos(phi)+y[1]*sin(theta)*sin(phi)+y[2]*cos(theta));
}

