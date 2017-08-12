#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>

#define MAX 20000
#define EPSILON 1.0E-6

using namespace std;
bool is_equal(double x, double y)
{
	return fabs(x - y) < EPSILON;
}

int main() {
	
	int Low, Up, Num, K, M = 0, D;
	double x[MAX], t[MAX], y[MAX], e[MAX];
	double  a, b, Vmin, Vmax, A, B, alpha, Ermsd ;

	//get();
	stringstream sPoint[MAX];
	string line[12], space, Point[MAX];
	//get value for a, b , Num, K, alpha
	for (int i = 1; i <= 11; i++) {
		getline(cin, line[i]);
	}
	stringstream snum(line[4]), salpha(line[5]), sA(line[6]), sB(line[7]), sK(line[8]);
	snum >> space  >> Num;
	salpha >> space >> alpha;
	sA >> space >> A;
	sB >> space >> B;
	sK >> space >> K;
	//get value for Point
	for (int i = 0; getline(cin, Point[i]); i++) {
		sPoint[i] << Point[i];
		sPoint[i] >> x[i] >> t[i];
		M++;
	}
	D = M/K;

	for (int k = 1; k <= K; k++) {
		Low = (k - 1)*D;
		if (k == K)
			Up = M;
		else
			Up = k*D;
		//find_a_b(Low, Up);
		double La, Lb, Lenght;
		a = A;
		b = B;
		for (int j = 0; j < Num; j++) {
			La = 0;
			Lb = 0;
			for (int i = 0; i < M; i++) {
				if (i < Low || i >= Up) {
					La += ((a*x[i] + b - t[i])*x[i]);
					Lb += (a*x[i] + b - t[i]);
				}
			}
			Lenght = sqrt(La*La + Lb*Lb);
			//ga = 
			//gb = 
			a = a - alpha*La / Lenght;
			b = b - alpha*Lb / Lenght;
		}

		//find_Ermsd(Low, Up)

		double Sum = 0;
		for (int i = Low; i < Up; i++) {
			Sum += (a*x[i] + b - t[i])*(a*x[i] + b - t[i]);
		}
		Ermsd = sqrt(Sum / (Up - Low));

		//find_y_e()

		for (int i = 0; i < M; i++) {
			e[i] = y[i] - t[i];
			y[i] = a*x[i] + b;	
		}

		//find_histogram(Low, Up)

		double Arg=0, Sumvariance = 0, variance;
		// Find the error e
		
		for (int i = Low; i < Up; i++) {
			e[i] = a*x[i] + b - t[i];
			Arg += e[i];
		}
		Arg = Arg / int(Up - Low);
		// Find the variance

		for (int i = Low; i < Up; i++) {
			Sumvariance += (e[i] - Arg)*(e[i] - Arg);
		}
		variance = sqrt(Sumvariance / int(Up - Low));
		// Find upper bound, lower bound

		Vmin = Arg - 3.0*variance;
		Vmax = Arg + 3.0*variance;
		
		// Split and find 9 points on the interval

		double Inter = (Vmax - Vmin) / 10.0;//length
		double Lpoint = Vmin;
		double Npoint = Lpoint + Inter;
		double hisarr[10] = {};
		int hissum = 0;
		for (int j = 0; j < 10; j++)
		{
			for (int k = Low; k < Up; k++){
				switch (j){
				case 9:{
					if ((e[k] > Lpoint && e[k] < Vmax) || is_equal(Vmax, e[k]) || is_equal(Lpoint, e[k])){
						hissum += 1;
						hisarr[j] += 1.0;
					}
					break;
				}
				default:{
					if ((e[k] > Lpoint && e[k] < Npoint) || is_equal(Lpoint, e[k])){
						hissum += 1;
						hisarr[j] += 1.0;
					}
					break;
				}
				}
			}
			Lpoint = Npoint;
			Npoint += Inter;
		}

		for (int j = 0; j < 10; j++){
			hisarr[j] = hisarr[j] / hissum;
		}


		cout << fixed << setw(10)  << setprecision(5) << right << a<< fixed << setw(10)  << setprecision(5) << right << b
			<< fixed << setw(10)  << setprecision(5) << right << Ermsd;
		for (int j = 0; j < 10; j++) {
			cout << fixed << setw(10)  << setprecision(5) << right << hisarr[j];
		}
		cout << endl;
	}
	return 0;
}