#include<iostream>
#include<math.h>
#include<stdio.h>  //printf

#define pi 3.14156
#define alpha 1/137
#define L 30000000 //[nm]
#define z 1.0
#define beta 1.0
#define n 1.0113
#define f(x) 2*pi*alpha*L*z*z*(1-1/beta/beta/n/n)/x/x //function

using namespace std;

class Calc
{
	static const int m = 100; //number of partitions

	double h;
	double x;
	double s;
	int k;

	public:
		void calcIntegral(double a, double b);
};

void Calc::calcIntegral(double a, double b)
{
	h = (b-a)/m;

	x = a;
	s = 0;

	for ( k=1; k <= m-1; k++){
		x = x + h;
		s = s + f(x);
	}

	s = h * ((f(a) + f(b))/2 + s);

	cout<<" result = "<< s <<endl;
}

int main()
{

	double a, b;

	try
	{
		cout<<"Start:"; cin>>a;
		cout<<"End:"; cin>>b;
	//	a = 300; //nm
	//	b = 650; //nm

		Calc objCalc;

		objCalc.calcIntegral(a, b);
	}
	catch(...){
		cout << "Error!" << endl;
		return -1;
		}

	return 0;
}
