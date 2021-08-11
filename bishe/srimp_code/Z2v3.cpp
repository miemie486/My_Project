/* 
	CALCULATION OF THE SURFACE IMPEDANCE OF SUPERCONDUCTORS
	THE ORIGINAL CODE DEVELOPED IN FORTRAN BY J. HALBRITTER, 
	 KERNFORSCHUNGSZENTRUM KARLSRUHE, F.R.D, 6.70
	MODIFIED IN C++ BY J. YAN, JLAB, 8/2002
*/
#include <stdio.h>
#include <iostream>
#include <complex>
#include <fstream>
#include <cstring>
 
#define NUM 100		// The number of temperatures

using namespace std;

void nonan();
void isum();
void haupt();
void wink();

// The data that will be used for program
class common_data {
public:
	int id[6], iq, ms, is, ip;
	double ak, a2, fl, o, gp;
	double pi, dq[6], dx[1000], dq1, ds[1000];
	complex<double> chs, chd, ci;
	complex<double> cds[1000], ct[1000], cs, cd;

} HL;

/*
		 The result data of calculation

	XS = PENETRATION DEPTH FOR SPECULAR REFLECTION
	XD = PENETRATION DEPTH FOR DIFFUSE REFLECTION
	RS = SURFACE RESISTANCE FOR SPECULAR REFLECTION
	RD = SURFACE RESISTANCE FOR DIFFUSE REFLECTION

	IT = # OF TEMPERATURE VALUES
	FO = FREQUENCY (GHz), TE(I) = TEMPERATURES (K)
	TC = CRITICAL TEMPERATURE
	FAK = OBSERVED DELTA(0) /( k * TC) 
	DLON = LONDON PENETRATION DEPTH AT T(0)
	XKOH = COHERENCE LENGTH
	FREI = MEAN FREE PATH OF NORMAL ELECTRONS
	MATERIAL PARAMETERS (UNITS : DEGREE KELVIN, ANGSTROM)
	"REST" SETS THE LEVEL OF ACCURACY OF THE COMPUTATION.	
*/

class result_data {
	int i;
public:
	int it;
	char ifname[64], ofname[64], matname[64];
	double g[6],ph,flo,al;
	double xs[NUM], xd[NUM], rs[NUM],rd[NUM];
	double og, fo, fak, tc, te[NUM], dlon, xkoh, frei, rest;
	complex<double> ochs[NUM], ochd[NUM];

public:
	void screen_input();
	int file_input();
	int save_data();
	void init_data();
	void screen_display();
	int file_save();
};

// Input data from terminal
void result_data::screen_input()
{
	cout << " PLEASE ENTER DATA INTERACTIVELY FROM KBD!" << endl;
	cout << " INPUT MATERIAL LABLE --> ";
	cin >> matname;
	cout << " EXCITATION FREQUENCY (GHz) IS .. ";
	cin >> fo;
	cout << " ENTER THE # OF TEMPERATURES AT WHICH CALCULATION IS PERFORMED .. ";
	cin >> it;
	for ( i=0; i< it; i++){
		cout << " TEMPERATURE (K) # ',I2, ' IS .. ";
		cin >> te[i];
	}
	cout << " CRITICAL TEMPERATURE (K) OF SAMPLE IS .. ";
	cin >> tc;
	cout << " OBSERVED VALUE OF -GAP(0) / k*Tc- IS .. ";
	cin >> fak;
	cout << " LONDON PENETRATION DEPTH (A) AT 0 K IS .. ";
	cin >> dlon;
	cout << " COHERENCE LENGTH (A) IS .. ";
	cin >> xkoh;
	cout << " NORMAL ELECTRON MEAN FREE PATH (A) IS .. ";
	cin >> frei;
	cout << " DESIRED % ACCURACY OF CALCULATION IS .. ";
	cin >> rest;
}

// Input data from a file
int result_data::file_input()
{
	cout << " ENTER FILE NAME --> ";
	cin >> ifname;
//	char ifname[20] = "Nb1500.txt";

//	cout<< ifname;

	ifstream ifin;
	ifin.open(ifname, ios::in);
	if (!ifin) {
		cout << "Cannot open input file.\n";
		return 0;
	}
	ifin >> matname;
	cout << matname << endl;
	ifin >> it;
	cout<< it << endl;
	ifin >> fo;
	cout << fo << "\t";
	for(int i=0; i<it;i++){
	 ifin >> te[i];
	 cout << te[i] <<"\t";
	}
	ifin >> tc >> fak >> dlon >> xkoh >>frei>>rest;
	cout <<tc<<"\t"<<fak<<"\t"<<dlon<<"\t"<<xkoh<<"\t"<<frei<<"\t"<<rest<<endl;

  return 1;
}

// Save the input data on a file
int result_data::save_data(){
	cout << " ENTER FILE NAME --> ";
	cin >> ofname;

	ofstream ifout;
	ifout.open(ofname, ios::out);
	if (!ifout) {
		cout << "Cannot open output file.\n";
		return 0;
	}
	ifout << matname << endl;
	ifout << it << endl;
	ifout << fo << "\t";
	for(int i=0; i<it;i++){
	 ifout << te[i]<<"\t";
	}
	ifout << endl;
	ifout<<tc<<"\t"<<fak<<"\t"<<dlon<<"\t"<<xkoh<<"\t"<<frei<<"\n"<<rest;
	return 1;
}

// Initialize the input data
void result_data::init_data() {

	fo = fo * 1.0e+9;
	rest = rest/10000.0;

//	cout << fo << "  " << te[0] << "  " <<te[1] <<"  " << te[2] << " ";
//	cout << te[3] << " " << te[4] <<" " << te[5] << "  ";
//	cout << tc << " " << fak <<" "<<endl;

// 	THE MOMENTUM INTEGRAL IS SPLIT IN 6 PARTS OVER THE INTERVAL 
//	g(I) TO g(I+1)
	double temp_g[6] = {1.0e-6, 1.0e-2, 0.2, 2.0, 5.0, 50.0};
	for (i = 0; i< 6; i++) g[i] = temp_g[i];
// 	THE INTERVALS ARE EVALUATED USING SIMPSON'S RULE.  
//	THE ABCISSA SPACING IS 1/ID(I).
	int id[6] = {5,7,15,20,20,9};

//  MODIFY INTERVAL SIZE FOR DESIRED COMPUTATIONAL ACCURACY
//	FOR "REST" = 0.0001 THE ACCURACY WILL BE ~ 1%.
//	THE LOWER LIMIT OF "REST" IS 1*10**(-1).
	HL.is = pow((6.0e-4 /rest),0.25) + 1.0;

	for (i=0; i< 6; i++)  { HL.id[i] = id[i] * HL.is;
//cout <<"ID are "<< HL.id[i] << "\t";
	}
// 	DEFINITION OF PARAMETERS & CONSTANTS
//	flo = COHERENCE LENGTH / MEAN FREE PATH
//	ph = h/k
//	al = GINZBERG-LANDAU PARAMETER
//	og = h*f/(2*PI*GAP)
	HL.pi = 3.1415926535;
	ph = 0.479*pow(10.0,-10);
	complex<double> temp_ci(0.0, 1.0);
	HL.ci = temp_ci;
	flo = xkoh / frei;
	al = dlon/xkoh;
	og = ph*fo / tc / fak;

//	SETUP STEP SIZES FOR CALCULATION
	HL.dx[0] = 0.0;
	HL.dx[1] = g[0];
	HL.iq = 1;
}

// DISPLAY RESULT ON SCREEN
void result_data::screen_display() {

	cout << endl;
	cout << "  SURFACE IMPEDANCE IN THE SUPERCONDUCTOR STATE"<< endl;
	cout << " -------------------------------------------------------"<<endl;
	cout << "  SAMPLE IDENTIFICATION : " << matname << endl;
	cout << "  MATERIAL PARAMETERS :" << endl;
	cout << "  TC = " << tc <<"K    GAP(T=0)/KTc = "<<fak  <<endl;
	cout << "  LONDON PEN. DEPTH(T=0, 1/L=0) = "<<dlon << " ANGSTROM"<<endl;
	cout << "  COHERENCE LENGTH(T=0, 1/L=0) = "<<xkoh << " ANGSTROM" <<endl;
	cout << "  MEAN FREE PATH = " << frei << " ANGSTROM" <<endl;
	cout << "  FOR F = " << fo<< "Hz" << endl;
	cout << "\t\t\t BCS-THEORY RESULTS" << endl;
	cout << "\t\t SPECULAR \t\t DIFFUSE REFLECTION" << endl<<endl;
	cout << "TEM \t PEN.\tSURFACE\tPEN.\tSURFACE\tSPECULAR\tDIFFUSE" <<endl;
	cout << "PERATURE DEPTH\tRESIST.\tDEPTH\tRESISTANCE\tSURFACE IMPEDANCE"<<endl;
	cout << " <K>\t   <A>\t <OHM> \t <A> \t <OHM>" <<endl;
	
//	cout.width(10);
	for (int i = 0; i<it; i++) 
	{
		cout.precision(4);
		cout.setf(ios::uppercase | ios::showpoint);
		cout << te[i]<<"\t";
		cout<<xs[i]<<"  "<<rs[i]<<"  "<<xd[i]<<"  "<<rd[i]<<"  "<<ochs[i]<<"  "<< ochd[i]<<endl;		
	}
}

// Save THE RESULT DATA
int result_data::file_save() {
	cout << " ENTER FILE NAME --> ";
	cin >> ofname;

	ofstream ofout;
	ofout.open(ofname, ios::out);
	if (!ofout) {
		cout << "Cannot open file.\n";
		return 0;
	}
	ofout << "SURFACE IMPEDANCE IN THE SUPERCONDUCTOR STATE"<< endl<<endl;
	ofout << "SAMPLE IDENTIFICATION : " << matname << endl;
	ofout << "MATERIAL PARAMETERS :" << endl;
	ofout << "TC = " << tc <<"K    GAP(T=0)/KTc = "<<fak  <<endl;
	ofout << "LONDON PEN. DEPTH(T=0, 1/L=0) = "<<dlon << " ANGSTROM"<<endl;
	ofout << "COHERENCE LENGTH(T=0, 1/L=0) = "<<xkoh << " ANGSTROM" <<endl;
	ofout << "MEAN FREE PATH = " << frei << " ANGSTROM" <<endl;
	ofout << "FOR F = " << fo<< "Hz" << endl;
	ofout << "\t\t\t BCS-THEORY RESULTS" << endl;
	ofout << "\t\t SPECULAR \t\t DIFFUSE \tREFLECTION" << endl<<endl;
	ofout << "TEM \t PEN.\tSURFACE\t\tPEN.\tSURFACE\tSPECULAR\t\tDIFFUSE" <<endl;
	ofout << "PERATURE DEPTH\tRESISTANCE\tDEPTH\tRESISTANCE\tSURFACE IMPEDANCE"<<endl;
	ofout << "<K>\t   <A>\t <OHM> \t\t <A> \t <OHM>" <<endl;
	for (int i = 0; i<it; i++) 
	{
		ofout << te[i]<<"\t"<<xs[i]<<"\t"<<rs[i]<<"\t"<<xd[i]<<"\t"<<rd[i]<<"\t"<<ochs[i]<<"  "<< ochd[i]<<endl;
	}
	return 1;
}

int main()
{
	result_data RD;
	int i, j;
	char iosw;

	cout << " INPUT MATERIAL DATA FROM DISK FILE? (Y) --> ";
	cin >> iosw;
//	iosw = 'Y';
	if (iosw != 'Y' ) {
// INPUT MATERIAL PARAMETERS INTERACTIVELY AND STORE FOR LATER USE.
		RD.screen_input();
	
// WRITE PARAMETER FILE TO DISK
		if (!RD.save_data()) {
			cout << "Cannot open save file.\n";
			exit(0) ;
		}
	}
	else {
		if (!RD.file_input()) {
			cout << "Cannot open input file.\n";
			exit(0);	
		}
	} 
	cout << "Working...\n";
    
	RD.init_data();

	int ll, d;
	double da, de, dz, delt, det;

//	LOOP OVER INTERVALS AND GET THE TOTAL STEP
	for(i=0; i<5; i++)
	{
		ll = HL.id[i] - 1;
		d = HL.id[i] *2;
		HL.id[i] = ll;
		da = RD.g[i];
		dz = RD.g[i+1];
		delt = (dz - da)/d;
		det = delt * 2.0;
		HL.dq[i] = delt/3.0;
		for (j = 0; j < ll; j++) 
		{
			da = da +det;
		    de = da - delt;
			HL.iq++;
			HL.dx[HL.iq] = de;
			HL.iq++;
			HL.dx[HL.iq] = da;
		}
		de = dz - delt;
		HL.iq++;
		HL.dx[HL.iq] = de;
		HL.iq++;
		HL.dx[HL.iq] = dz;
	}

	ll = HL.id[5] - 1;
	d = HL.id[5] * 2;
	HL.id[5] = ll;
	da = RD.g[0];
	dz = 1.0/RD.g[5];
	delt = (dz-da)/d;
	det  = delt * 2;
	HL.dq[5] = delt/3.0;
	HL.iq++;
	HL.dx[HL.iq] = 1.0/da;

	for (i=0;i<ll;i++){
		da = da + det;
		de = da - delt;
		HL.iq++;
		HL.dx[HL.iq] = 1.0/de;
		HL.iq++;
		HL.dx[HL.iq] = 1.0/da;
	}
	de = dz - delt;
	HL.iq++;
	HL.dx[HL.iq] = 1.0/de;
	HL.iq++;
	HL.dx[HL.iq] = 1.0/dz;
	HL.iq++;
	
	double t,b,ag;

// LOOP OVER TERMPERATURES TO BE INVESTIGATED
	for (i=0;i<RD.it;i++){
		t = RD.te[i]/RD.tc;
		// INCLUDE TEMPERATURE DEPENDENCE OF GAP
		b = HL.pi/2.0 * (1.0 - t*t);
		ag = sqrt(sin(b));
		HL.fl = RD.flo/ag;
		HL.o = RD.og/ag;
		if ( HL.o < 0.5){
			HL.gp = ag/t * RD.fak;
			HL.is = 1.0/sqrt(HL.o * RD.rest)/5.0;
			HL.ip = 1.0/sqrt((2.0 - HL.o) * RD.rest)/5.0;
			for (int j=0; j<HL.iq; j++){
				HL.ds[j] = 0.0;
				HL.cds[j] = (0.0,0.0);
			}

			nonan();

			HL.ms = HL.gp/(sqrt(RD.rest)*2.0*HL.pi);
			
			isum();

			HL.dq1 = abs(HL.cds[0]);
			double enp = 1.0/sqrt(HL.dq1);
			HL.ak = RD.al * enp * ag;
			HL.a2 = HL.ak * HL.ak;

			haupt();

			HL.chs = HL.chs *  enp;
			RD.ochs[i] = HL.chs;
			HL.chd = HL.chd * enp;
			RD.ochd[i] = HL.chd;
			RD.xs[i] = real(HL.chs) * RD.dlon;
			RD.xd[i] = real(HL.chd) * RD.dlon;
			b = (RD.fo/pow(10.0,17)*8.0*HL.pi*HL.pi*RD.dlon);
			RD.rs[i] = imag(HL.chs) * b;
			RD.rd[i] = imag(HL.chd) * b;
		}
	}
	
	RD.screen_display();
	if (!RD.file_save()) {
		cout << "Cannot open save file.\n";
		exit(0) ;
	}
}


// ROUTINE TO EVALUATE THE FUNCTION J(u)
void wink(complex<double> cl)
{
	int i;
	double d, u;
	complex<double> cc, c2;
	complex<double> ctmp[5];

	for (i=0; i<HL.iq; i++)
	{
		d = HL.dx[i];
		HL.cd = d/cl;
		c2  = HL.cd * HL.cd;
		u = abs(HL.cd);
		if(u > 0.2) {
			HL.cd = 1.0/HL.cd;
			ctmp[0] = log(HL.cd + HL.ci) - log(HL.cd - HL.ci);
		HL.ct[i] = 0.75/d*(-2.0*HL.cd + (1.0+HL.cd*HL.cd)/HL.ci*ctmp[0]);
		}
		else {
			ctmp[0] = 1.0/195.0 - c2/255.0;
			ctmp[1] = 1.0/143.0 - c2*ctmp[0];
			ctmp[2] = 1.0/99.0 - c2*ctmp[1];
			ctmp[3] = 1.0/63.0 - c2*ctmp[2];
			ctmp[4] = 1.0/35.0 - c2*ctmp[3];
			ctmp[5] = 1.0/15.0 - c2*ctmp[4];
			HL.ct[i] = 3.0/cl *(1.0/3.0 - c2*ctmp[5]);
		}
	}
}

// 	ROUTINE TO COMPUTE THE ANALYTICAL PORTION OF "Q" (Qa).
//	Qa IS REAL AND GIVES THE PENETRATION DEPTH FOR hw/2*PI*DELTA << 1.
//	Qa DESCRIBES THE CURRENT OF COOPER PAIRS

void isum()
{
	int i,j;
	double b;
	double g,dw;
	complex<double> cb,cw,cl,cp,ca;

	g = HL.pi /HL.gp;
	for (i=0; i<HL.iq; i++) HL.ds[i]=0.0;
	for (i=0;i<HL.ms; i++)
	{
		b = (2*(HL.ms - i-1)+1)*g;
		cb = b + HL.o * HL.ci;
		dw = sqrt(b * b + 1.0);
		cw = sqrt(cb * cb + 1.0);
		cl =(dw + cw)/2.0 + HL.fl;
		cp = dw * cw;
		ca = (1.0 - b * cb + cp)/cp;
		wink(cl);	
		for (j = 0; j<HL.iq; j++){
			HL.ds[j] = HL.ds[j] + (ca * HL.ct[j]).real();
		}
	}	
	for (i = 0; i<HL.iq; i++){
		HL.cds[i] = HL.ds[i] * g + HL.cds[i];
	}
}

// ROUTINE NONAN TO EVALUATE Q
void nonan()
{
	int i, j;
	double d, dv, dz, du, dg, da, dc ,dp, dm ,dw, dwo, ot, o2;
	double b, bs , dtmp[5];
	complex<double> cy;

	ot = HL.o/2.0 * HL.gp;

	if(ot <= 0.2){
		o2 = ot*ot;
		dv = ot*(1.0+o2/6.0*(1.0+o2/20.0 * (1.0+o2/42.0*(1.0+o2/72.0*(1.0+o2/110.0)))))*exp(-ot);
		ot = -ot * 2.0;
	}
	else {
		ot = -ot * 2.0;
		dv = (1.0 - exp(ot))/2.0;
	}

	b = 4 * HL.is;
	bs = HL.pi/b;
	for ( i = 0; i<(HL.is); i++)
	{
// cout << "*" ;
		b = 4*(HL.is - i-1) +2;
		d = b * bs;
		dtmp[0] = cos(d);
		dz = 0.5 * (1.0 + dtmp[0]);
		du = HL.gp / dz;
		da = exp(-du);
		dg = exp(-du + ot);
		dc = da/((1.0 + dg) * (1.0 + da)) * dv;
		dp = dz + 1.0;
		dm = 1.0 -dz;
		du = dp * dm;
		dw = sqrt(du);
		da  = dz * HL.o;
		dwo = sqrt(du + da*(2.0 + da));
		dp = dz * sqrt(dz*dp);
		du = da/du;
		if(du > 0.1) {
			dm = dwo - dw;
		}
		else {
			du = du * (2.0 + da) /2.0;
			dtmp[1] = (1.0 - 13.0 * du/9.0 * (1.0 - 1.4 * du));
			dtmp[2] = (1.0 - 11.0 * du/7.0 * (1.0 - 1.5 * du * dtmp[1]));
			dtmp[3] = (1.0 - 1.4 * du * (1.0 - 1.5 * du * dtmp[2]));
			dtmp[4] = (1.0 - 1.25 * du * dtmp[3]);
			dm = du * (1.0 - 0.5 * du * (1.0 - du * dtmp[4])) * dw;
		}
		cy = -HL.ci * dm/2.0/dz + HL.fl;
		wink(cy);

		du = dc/dp/dwo;
		double dwu = dwo * dw;
		dm = 1.0 + dz * (dz + HL.o);
		double dps = (dm + dwu) * du;
		double dms = (dm - dwu) * du;
		for (j=0; j<HL.iq; j++){
			HL.cds[j] = HL.cds[j] + dps * HL.ct[j];
		}
		dm = dwo + dw;
		cy = HL.ci * dm/2.0/dz + HL.fl;

		wink(cy);
		for (j =0; j<HL.iq; j++) HL.ds[j] = HL.ds[j] + dms * HL.ct[j].real();
	}
	cy = -HL.ci * bs * 4.0;
	for(j=0; j < HL.iq; j++) {
		HL.cds[j] = (HL.cds[j] + HL.ds[j]) * cy;
		HL.ds[j] = 0.0;
	}
	
	b = 4 * HL.ip;
	bs = HL.pi/b;
	for(i=0; i<HL.ip; i++)
	{
		b = 4*(HL.ip - i-1) +2;
		d = b * bs;
		dz = 1.0 + 0.5 * HL.o * (1.0 + cos(d));
		du = HL.gp * dz;
		da = exp(- du);
		dg = exp(- du - ot);
		dc = dg/((1.0 + dg) * (1.0 + da)) * dv;
		dp = dz + 1.0;
		dm = dz - 1.0;
		dw = sqrt(dp*dm);
		da  = HL.o - dm;
		du = dp - HL.o;
		dwo = sqrt(da * du);	

		cy = dwo + HL.ci * dw + HL.fl;

		wink(cy);
		cy = (1.0 + dz*(dz - HL.o) + HL.ci * dw * dwo)/sqrt(dp*du) * dc;	
		for(j=0; j<HL.iq; j++) HL.ds[j] = HL.ds[j] + (cy * HL.ct[j]).real();
	}
	b = +bs * 4.0;

	for (j=0; j<HL.iq; j++) {
		HL.cds[j] = HL.cds[j] + HL.ds[j]*b;
	}
}

//     ROUTINE TO CALCULATE THE INTEGRAND IN THE MOMENTUM INTEGRATION
//     FOR THE SURFACE IMPEDANCE.   THE MOMENTA HAVE BEEN TRANSFORMED 
//     TO q-SPACE (i.e. MOMENTUM MEASURED IN UNITS OF COHERENCE LENGTH).
void intq(int l, int m)
{
	double dx2, da;
	complex<double> cdg;
	
	cdg = HL.cds[l-1] / HL.dq1;
	dx2 = HL.dx[l-1] * HL.dx[l-1];
	da = dx2 * HL.a2;

	if (m == 2) {
		HL.cs = 1.0 / (HL.a2 + cdg /dx2);
		HL.cd = log(1.0 + cdg/da) * dx2;
	}
	else {
		HL.cs = 1.0 /(da + cdg);
		HL.cd = log(1.0 + cdg/da);
	}
}


//	SUBROUTINE HAUPT
//
//	ROUTINE TO PERFORM THE INTEGRATION IN THE SPECULAR AND 
//	DIFFUSE REFLECTION TERMS FOR THE SURFACE IMPEDANCE.

void haupt()
{
	int i,j,k, ll, m;
	double a;
	complex<double> cis, cid, ces, ced, cas, cad;

	m = 1;
	k = 2;
	 
	intq(k,m);
	k++;
	cis = HL.dx[1]* HL.cs;
	cid = HL.dx[1] * HL.cd + 2.0 * HL.dx[1];

	for(j = 0 ; j<6; j++)
	{
		if (j >= 5) {
			m=2;
			intq(k,m);
			k++;
			cis = cis + HL.cs*HL.dx[1];
			cid = cid + HL.cd*HL.dx[1];
		}
		ces = (0.0, 0.0);
		ced = (0.0, 0.0);
		cas = 0.5 * HL.cs;
		cad = 0.5 * HL.cd;
		ll = HL.id[j];
		for ( i = 0; i<ll; i++){
			intq(k,m);
			k++;
			ces = ces + HL.cs;
			ced = ced + HL.cd;
			intq(k,m);
			k++;
			cas = cas + HL.cs;
			cad = cad + HL.cd;
		}
		
		intq(k,m);
		k++;
		ces = ces + HL.cs;
		ced = ced + HL.cd;
		intq(k,m);
		k++;
		cis= cis + HL.dq[j] * (4.0 *ces + (2.0 * cas + HL.cs));
		cid = cid + HL.dq[j]*(4.0*ced + (2.0* cad + HL.cd));
	}	
		a = 2.0/HL.pi* HL.ak;
		HL.chs = a * cis;
		HL.chd = 2.0/a/cid;	
}
