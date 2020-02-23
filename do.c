#include <stdio.h>
#include <math.h>

// 刺激なしのときのperiodは36.13-25.57=10.79とする。(3rdpeak-2ndpeak);
double period = 10.79;
double first_peak = 14.78;

int STEPNO = 10000;

double x,y,z,t;
double x_a,y_a,z_a;
double x_b,y_b,z_b;
double x_c,y_c,z_c;
double x_d,y_d,z_d;
double dt = 0.01;

double dx = 0.1;
double dy = 0.5;
double k = 1;
int n = 10;

double i_starttime = 0;
double i_endtime = 0;

double input_strength = 0;

double xx[100000];


void peakDetect(FILE *f){
	char str[1000];
	sprintf(str,"%f\t%f\t",i_starttime,i_endtime);

	for(int i=1;i<STEPNO;i++){
		if (xx[i-1] < xx[i] && xx[i] > xx[i+1]){
			printf("peak: %f\n",i*dt);
			sprintf(str,"%s\t%f",str,i*dt);
		}
	}
	sprintf(str,"%s\n",str);
	fprintf(f,str);
}


void init(){

	t = 0;
	x = 0.5;
	y = 1;
	z = 1;

	for(int i=0;i<sizeof(xx)/sizeof(xx[0]);i++){
		xx[i]=0;
	}

}

double i(t){
	double ans=0;
	if(i_starttime <= t && t <= i_endtime){
		ans = 1;
	}

	return ans;
}


double f(double x, double y, double z){
	double ans = 0;
	ans = 1.0/(k+pow(z,n)) - dx * x + input_strength * i(t);
	return ans;
}

double g(double x, double y, double z){
	double ans = 0;
	ans = x - dy * y;
	return ans;
}

double h(double x, double y, double z){
	double ans = 0;
	ans = y - z;
	return ans;
}


void calA(){
	x_a = f(x,y,z) * dt;
	y_a = g(x,y,z) * dt;
	z_a = h(x,y,z) * dt;
}

void calB(){
	x_b = f(x+x_a/2.0,y+y_a/2.0, z+z_a/2.0) * dt;
	y_b = g(x+x_a/2.0,y+y_a/2.0, z+z_a/2.0) * dt;
	z_b = h(x+x_a/2.0,y+y_a/2.0, z+z_a/2.0) * dt;
}

void calC(){
	x_c = f(x+x_b/2.0,y+y_b/2.0, z+z_b/2.0) * dt;
	y_c = g(x+x_b/2.0,y+y_b/2.0, z+z_b/2.0) * dt;
	z_c = h(x+x_b/2.0,y+y_b/2.0, z+z_b/2.0) * dt;
}

void calD(){
	x_d = f(x+x_c, y+ y_c, z+z_c) * dt;
	y_d = g(x+x_c, y+ y_c, z+z_c) * dt;
	z_d = h(x+x_c, y+ y_c, z+z_c) * dt;
}

void newVal(){
	x = x + (x_a + 2 * x_b + 2 * x_c + x_d)/6.0;
	y = y + (y_a + 2 * y_b + 2 * y_c + y_d)/6.0;
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
	printf("%f,%.20f,%.20f,%.20f\n",t,x,y,z);
}

void fw(FILE *f){

	 fprintf(f, "%f\t%.10f\t%.10f\t%.10f\n",t,x,y,z);
}


void runRK(int STEP, FILE *f){
	xx[0] = x;
	init();

	for(int i=0;i<STEP;i++){
		oneStep();
		xx[i+1] = x;
		fw(f);
	}
}

int main(void){
	FILE *f = fopen("./output/do.txt", "w");
	double input_hours = 2;
	FILE *g = fopen("./output/peak0.02_2.txt", "w");

	runRK(STEPNO, f);
	peakDetect(g);
	input_strength = 0.02;

	for(int i = 0;i<30;i++){
		i_starttime = first_peak+ i*(period / 24.0);
		i_endtime = i_starttime + input_hours*(period / 24.0);
		runRK(STEPNO, f);
		peakDetect(g);
	}
	fclose(f);
	fclose(g);

	return 0;
}
