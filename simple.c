#include <stdio.h>
#include <math.h>

int STEPNO = 10000;

double x,y,z,t;
double x_a,y_a,z_a;
double x_b,y_b,z_b;
double x_c,y_c,z_c;
double x_d,y_d,z_d;
double dt = 0.01;

double lambda = 0.1, w0=2*M_PI/1.2, I=1, w=2*M_PI/1, e=1;

int n = 10;

double i_starttime = 0;
double i_endtime = 0;

double input_strength = 0;

double xx[100000];


double amplitude(int starttime, int endtime){
	double min=100000, max=0;
		for(int i=starttime;i<=endtime;i++){
			if(xx[i] < min) min = xx[i];
			if(xx[i] > max) max = xx[i];
		}
	return max-min;
}

void peakDetect(FILE *f){
	char str[1000];
	double first =0;
	double second=0;
	for(int i=1;i<STEPNO;i++){
		if (xx[i-1] < xx[i] && xx[i] > xx[i+1]){
			printf("%f\n",i*dt);
			sprintf(str,"%s\t%f\t%f\n",str,i*dt,xx[i]);
		}

	}
	fprintf(f,"%s",str);
}

double periodDetect(){
	double first =0;
	double second=0;
	for(int i=1;i<STEPNO;i++){
		if (xx[i-1] < xx[i] && xx[i] > xx[i+1]){

			first = second;
			second = i*dt;
		}

	}

	return second-first;

}

double ampDetect(){
	return amplitude((int)(STEPNO*0.9), STEPNO);
}

void ampPeriodDetect(){
	printf("%f\t%f\n",periodDetect(),ampDetect());
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



double f(double x, double y, double z){
	double ans = 0;
	ans = y;
	return ans;
}

double g(double x, double y, double z){
	double ans = 0;
	ans = -w0*w0 * x - 2 * lambda * y + I * cos (z);
	return ans;
}

double h(double x, double y, double z){
	double ans = 0;
	ans = w + e * sin(z) * x;
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
	z = fmod(z, 2*M_PI);

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

	char filename[1000];

	for(int i=0;i<=10;i++){
		I = i*0.1;
		e=1-I;
		sprintf(filename, "./output/simple_I_%.1f.txt",I);
		FILE *f = fopen(filename, "w");
		runRK(STEPNO, f);
		printf("%.1f\t",I);
		ampPeriodDetect();
		fclose(f);
	}

	return 0;
}
