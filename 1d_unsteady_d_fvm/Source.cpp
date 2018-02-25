#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string.h>
using namespace std;
typedef struct temparature_data {
	int boundary_condition;
	double x,y,thermal_conductivity,temparature;
}temp_data;
int conv_stat;
temp_data *temp = new temp_data[101];
double delta_t,delta_x,delta_y,aW,aE,aP,aP0,a[101],b[101],d[101],source[101],var[101],error[101],Conductivity,Density,SpecificHeat;
void data_initialization();
void creat_grid();
void check_convergence();
void tdma();
void boundary_condition_assingment();
int const n = 101;
int main() {
	ofstream myfile2;
	myfile2.open ("result.txt");
	int i,j,k,done = 0;
	data_initialization();
	creat_grid();
	boundary_condition_assingment();
	for(int iter=0;(iter<100000)&&(done != 1);iter++) {
		for(int xx=1;xx<(n-1);xx++) {
			source[xx] = aP0*temp[xx].temparature;
		}
		tdma();
		for(int xx=1;xx<(n-1);xx++) {
			error[xx] = abs(var[xx] - temp[xx].temparature);
			temp[xx].temparature = var[xx];
		}
		if(iter%100 == 0) {
			for(int xx=0;xx<n;xx++) {
				myfile2 << temp[xx].temparature<<endl;
			}
		}
		conv_stat = 0;
		check_convergence();
		if(conv_stat == 0) {
			cout << "Solution converged at timestep"<<" "<<iter+1<<endl;
			for(int xx=0;xx<n;xx++) {
				myfile2 << temp[xx].temparature<<endl;
			}
			done = 1;
		}
		cout <<"colution completed at timestep "<< iter+1<<endl;;
	}
	myfile2.close();
	if(done == 1) {
		cout << "Solution converged !!!\n";
	}else
		cout << "Solution Not Converged More iterations required !!\n";

}

void creat_grid() {
	ofstream myfile1;
	myfile1.open ("grid.txt");
	for(int i = 1;i<n;i++) {
		temp[i].x = delta_x*(i);
		temp[i].boundary_condition = 0;
	}
	temp[0].boundary_condition = 1; temp[n-1].boundary_condition = 1;	
	for(int i=0;i<n;i++) {
		myfile1 << temp[i].x<<"  "<<temp[i].boundary_condition<<"  "<<temp[i].temparature<<endl;
	}
	myfile1.close();
}
void data_initialization() {
	int i,j;
	delta_x = 0.01; Conductivity = 10; delta_t = 0.12,delta_y = 0.01,Density = 100,SpecificHeat = 200;
	aW = Conductivity/delta_x;
	aE = Conductivity/delta_x;
	aP0 = Density*SpecificHeat*(delta_x/delta_t);
	aP = aW + aE + aP0;
	for(i=1;i<(n-1);i++) {
		d[i] = aP;
	}
	d[0] = 1; d[n-1] = 1;
	for(i=1;i<(n-1);i++) {
		a[i] = aW;
	}
	a[0] = 0.0; a[n-1] = 0.0;
	for(i=1;i<(n-1);i++) {
		b[i] = aE;
	}
	b[0] = 0.0; b[n-1] = 0.0;
	source[0] = 100; source[n-1] = 0;
	for(int xx=0;xx<n;xx++) {
			temp[xx].temparature = 0;
	}
}
void boundary_condition_assingment() {
	temp[0].temparature = 100; temp[n-1].temparature = 0;
}
void tdma() {
	int i,j;
	double aa[n],cc[n];
	aa[0] = 0; aa[n-1] = 0;
	for(i=1;i<(n-1);i++) {
		aa[i] = a[i]/(d[i] - (b[i]*aa[i-1]));
	}
	cc[0] = source[0]; cc[n-1] = source[n-1];
	for(i=1;i<(n-1);i++) {
		cc[i] = (b[i]*cc[i-1] + source[i])/(d[i] - (b[i]*aa[i-1]));
	}
	var[n-1] = source[n-1];
	for(i=n-1;i>0;i--) {
		var[i] = (aa[i]*var[i+1]) + cc[i];
	}
}
void check_convergence() {
	for(int xx=0;xx<n;xx++) {
			if(error[xx] >= 0.00001) {
				conv_stat += 1;
			}
	}
}
//void write_data() {
	//myfile2 << "Writing this to a file.\n";
//	for(int i=0;i<n;i++) {
//		myfile2 << temp[i].temparature<<endl;
//	}
//	myfile2.close();
//}
