#include <stdio.h>
#include <math.h>

int STEPNO = 30000;

double x,z,t,U,T,D,S;
double x_a,U_a,T_a,D_a,S_a,z_a;
double x_b,U_b,T_b,D_b,S_b,z_b;
double x_c,U_c,T_c,D_c,S_c,z_c;
double x_d,U_d,T_d,D_d,S_d,z_d;
double dt = 0.01;

//Parameters
double k1=0.5, k2=1, k4=0.5;
double d1=0.0385, d3=0.5;
double k0TU=0.21,k0UT=0,kaTU=0.0798,kaUT=0.479;
double k0DT=0,   k0TD=0,kaDT=0.173, kaTD=0.213;
double k0DS=0.31*0.91,k0SD=0,kaDS=-0.32, kaSD=0.5057;
double k0SU=0.11,k0US=0,kaSU=-0.133,kaUS=0.0532;
double K=0.1, V=0.0078;
double vd=0.022;
double ks=2.7;

double KaiA=1.3, khalf=0.43;
int n = 10;

double xx[100000];


void peakDetect(FILE *f, double i_starttime, double i_endtime){
	char str[1000];
	sprintf(str,"%f\t%f\t",i_starttime,i_endtime);

	for(int i=1;i<STEPNO;i++){
		if (xx[i-1] < xx[i] && xx[i] > xx[i+1]){
			printf("peak: %f\n",i*dt);
			sprintf(str,"%s\t%f",str,i*dt);
		}
	}
	sprintf(str,"%s\n",str);
	fprintf(f,"%s\n",str);
}


void init(){

	t = 0;
	x = 8;
	U = 2.23;
	T = 0.83;
	D = 0.31;
	S = 0.48;
	z = 0;
	for(int i=0;i<sizeof(xx)/sizeof(xx[0]);i++){
		xx[i]=U;
	}

}


double A(double S){
	return fmax(0,KaiA-2*S);
}

double kxy(double k0, double ka, double S){
	return k0 + ka*A(S)/(khalf + A(S));
}

double f_x(double x, double U, double T, double D, double S, double z){
	double ans = 0;
	ans = k1/(k2+pow(z,n)) - d1 * x;
	return ans;
}

double f_U(double x, double U, double T, double D, double S, double z){
	double ans = 0;
	ans = kxy(k0TU,kaTU,S)*T + kxy(k0SU,kaSU,S)*S - kxy(k0UT,kaUT,S)*U - kxy(k0US,kaUS,S)*U - V*U/(K+U) - vd*U + ks*x;
	return ans;
}

double f_T(double x, double U, double T, double D, double S, double z){
	double ans = 0;
	ans = kxy(k0UT,kaUT,S)*U + kxy(k0DT,kaDT,S)*D - kxy(k0TU,kaTU,S)*T - kxy(k0TD,kaTD,S)*T - V*T/(K+T) - vd*T;
	return ans;
}

double f_D(double x, double U, double T, double D, double S, double z){
	double ans = 0;
	ans = kxy(k0TD,kaTD,S)*T + kxy(k0SD,kaSD,S)*S - kxy(k0DT,kaDT,S)*D - kxy(k0DS,kaDS,S)*D - V*D/(K+D) - vd*D;
	return ans;
}

double f_S(double x, double U, double T, double D, double S, double z){
	double ans = 0;
	ans = kxy(k0US,kaUS,S)*U + kxy(k0DS,kaDS,S)*D - kxy(k0SU,kaSU,S)*S - kxy(k0SD,kaSD,S)*S - V*S/(K+S) - vd*S;
	return ans;
}

double f_z(double x, double U, double T, double D, double S, double z){
	double ans = 0;
	ans = k4*(D) -d3*z;
	return ans;
}

void calA(){
	x_a = f_x(x, U, T, D, S, z) * dt;
	U_a = f_U(x, U, T, D, S, z) * dt;
	T_a = f_T(x, U, T, D, S, z) * dt;
	D_a = f_D(x, U, T, D, S, z) * dt;
	S_a = f_S(x, U, T, D, S, z) * dt;
	z_a = f_z(x, U, T, D, S, z) * dt;
}

void calB(){
	x_b = f_x(x+x_a/2.0, U+U_a/2.0, T+T_a/2.0, D+D_a/2.0, S+S_a/2.0, z+z_a/2.0) * dt;
	U_b = f_U(x+x_a/2.0, U+U_a/2.0, T+T_a/2.0, D+D_a/2.0, S+S_a/2.0, z+z_a/2.0) * dt;
	T_b = f_T(x+x_a/2.0, U+U_a/2.0, T+T_a/2.0, D+D_a/2.0, S+S_a/2.0, z+z_a/2.0) * dt;
	D_b = f_D(x+x_a/2.0, U+U_a/2.0, T+T_a/2.0, D+D_a/2.0, S+S_a/2.0, z+z_a/2.0) * dt;
	S_b = f_S(x+x_a/2.0, U+U_a/2.0, T+T_a/2.0, D+D_a/2.0, S+S_a/2.0, z+z_a/2.0) * dt;
	z_b = f_z(x+x_a/2.0, U+U_a/2.0, T+T_a/2.0, D+D_a/2.0, S+S_a/2.0, z+z_a/2.0) * dt;
}


void calC(){
	x_c = f_x(x+x_b/2.0, U+U_b/2.0, T+T_b/2.0, D+D_b/2.0, S+S_b/2.0, z+z_b/2.0) * dt;
	U_c = f_U(x+x_b/2.0, U+U_b/2.0, T+T_b/2.0, D+D_b/2.0, S+S_b/2.0, z+z_b/2.0) * dt;
	T_c = f_T(x+x_b/2.0, U+U_b/2.0, T+T_b/2.0, D+D_b/2.0, S+S_b/2.0, z+z_b/2.0) * dt;
	D_c = f_D(x+x_b/2.0, U+U_b/2.0, T+T_b/2.0, D+D_b/2.0, S+S_b/2.0, z+z_b/2.0) * dt;
	S_c = f_S(x+x_b/2.0, U+U_b/2.0, T+T_b/2.0, D+D_b/2.0, S+S_b/2.0, z+z_b/2.0) * dt;
	z_c = f_z(x+x_b/2.0, U+U_b/2.0, T+T_b/2.0, D+D_b/2.0, S+S_b/2.0, z+z_b/2.0) * dt;
}

void calD(){
	x_d = f_x(x+x_c, U+U_c, T+T_c, D+D_c, S+S_c, z+z_c) * dt;
	U_d = f_U(x+x_c, U+U_c, T+T_c, D+D_c, S+S_c, z+z_c) * dt;
	T_d = f_T(x+x_c, U+U_c, T+T_c, D+D_c, S+S_c, z+z_c) * dt;
	D_d = f_D(x+x_c, U+U_c, T+T_c, D+D_c, S+S_c, z+z_c) * dt;
	S_d = f_S(x+x_c, U+U_c, T+T_c, D+D_c, S+S_c, z+z_c) * dt;
	z_d = f_z(x+x_c, U+U_c, T+T_c, D+D_c, S+S_c, z+z_c) * dt;
}

void newVal(){
	x = x + (x_a + 2 * x_b + 2 * x_c + x_d)/6.0;
	U = U + (U_a + 2 * U_b + 2 * U_c + U_d)/6.0;
	T = T + (T_a + 2 * T_b + 2 * T_c + T_d)/6.0;
	D = D + (D_a + 2 * D_b + 2 * D_c + D_d)/6.0;
	S = S + (S_a + 2 * S_b + 2 * S_c + S_d)/6.0;
	z = z + (z_a + 2 * z_b + 2 * z_c + z_d)/6.0;
}

void oneStep(){
	calA();
	calB();
	calC();
	calD();
	newVal();
	t = t + dt;

}



void disp(){
	printf("%f,%.20f,%.20f,%.20f%.20f\n",t,U,T,D,S);
}

void fw(FILE *f){

	 fprintf(f, "%f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",t,x,U,T,D,S,z);
}


void runRK(int STEP, FILE *f){
	xx[0] = U;
	init();

	for(int i=0;i<STEP;i++){
		oneStep();
		xx[i+1] = U;
		fw(f);
	}
}




int main(void){
	V=0.001;
	vd=0.001;
	ks=0.001;
	k4=1;


	FILE *f = fopen("./output/wt_k0DSr.txt", "w");
	runRK(STEPNO, f);
	fclose(f);
	return 0;
}
