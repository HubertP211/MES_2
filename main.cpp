#include <iostream>
#include <fstream>
using namespace std;


int liczba_elementow;
int liczba_wezlow;

double global_k;

double global_temp_otoczenia;
double global_alfa;

double global_r_max;
double global_r_delta;
double global_tau_delta;
double global_c;
double global_ro;
double global_temp_pocz;
char wybor;
int wybor2;
double tau;
double czas = 0;


ifstream plik;
ofstream plik2;
ofstream plik3;


struct wezel
{
	long double temp;
	int warunek_brzegowy;
	double r;
};

struct element
{
	int nod1;
	int nod2;
	
	double k_element;
	double r_delta;
	double h_lokalne[2][2];
	double p_lokalne[2];

};


bool loadData(string nazwaPliku)
{
	
	plik.open(nazwaPliku.c_str());
	if(!plik.good())
	{
		return false;
	}

	while(true)
	{
		plik >> liczba_elementow >>   global_k >>  global_temp_otoczenia >> global_alfa >>  global_r_max >>
			global_r_delta >> global_tau_delta >> global_c >> global_ro >> global_temp_pocz;
		break;
	}
	
	liczba_wezlow = liczba_elementow + 1;

}

bool loadData_element(string nazwaPliku, element *tablica_elementow, wezel *tablica_wezlow)
{
	
	if(!plik.good())
	{
		return false;
	}

	for(int i = 0; i<liczba_elementow; i++)
	{
		while(true)
		{
			plik >> tablica_elementow[i].nod1;
			plik >> tablica_elementow[i].nod2;
			plik >> tablica_elementow[i].k_element;
			tablica_elementow[i].r_delta = tablica_wezlow[i+1].r - tablica_wezlow[i].r;
			
			break;
		}
	}
}

bool loadData_wezel(string nazwaPliku, wezel *tablica_wezlow)
{
	
	if(!plik.good())
	{
		return false;
	}

	for(int i = 0; i<liczba_wezlow; i++)
	{
		while(true)
		{
			plik >> tablica_wezlow[i].temp;
			plik >> tablica_wezlow[i].warunek_brzegowy;
			plik >> tablica_wezlow[i].r;
			break;
		}
	}
}


void generateLocalMatrixForElement(element *tablica_elementow, wezel *tablica_wezlow, long double *warunki)
{

	for(int i = 0; i<liczba_wezlow; i++)
	{
		warunki[i] = 0;
	}

	for(int i = 0; i<liczba_elementow; i++)
	{
		double r1 = 0;
		double r2 = 0;
		double w1 = 1;
		double w2 = 1;
		double ni1 = 0.7885;;
		double ni2 = 0.2115;
		double nj1 = 0.2115;
		double nj2 = 0.7885;

		r1 = (0.5 * (1 + 0.577 ) * tablica_wezlow[tablica_elementow[i].nod1].r ) + ( 0.5 * (1 - 0.577 ) * tablica_wezlow[tablica_elementow[i].nod2].r );
		r2 = (0.5 * (1 -  0.577 ) * tablica_wezlow[tablica_elementow[i].nod1].r ) + ( 0.5 * (1 + 0.577 ) * tablica_wezlow[tablica_elementow[i].nod2].r );
		

		tablica_elementow[i].h_lokalne[0][0] = ((global_k / tablica_elementow[i].r_delta) * (r1*w1 + r2*w2)) + 
			(((global_c*global_ro*tablica_elementow[i].r_delta)/global_tau_delta) * (ni1*ni1*r1*w1 + ni2*ni2*r2*w2)); 
		
		tablica_elementow[i].h_lokalne[0][1] = ((-1)*(global_k / tablica_elementow[i].r_delta) * (r1*w1 + r2*w2)) + 
			(((global_c*global_ro*tablica_elementow[i].r_delta)/global_tau_delta) * (ni1*nj1*r1*w1 + ni2*nj2*r2*w2));
		
		tablica_elementow[i].h_lokalne[1][0] = ((-1)*(global_k / tablica_elementow[i].r_delta) * (r1*w1 + r2*w2)) + 
			(((global_c*global_ro*tablica_elementow[i].r_delta)/global_tau_delta) * (ni1*nj1*r1*w1 + ni2*nj2*r2*w2));  

		tablica_elementow[i].h_lokalne[1][1] = ((global_k / tablica_elementow[i].r_delta) * (r1*w1 + r2*w2)) + 
			(((global_c*global_ro*tablica_elementow[i].r_delta)/global_tau_delta) * (nj1*nj1*r1*w1 + nj2*nj2*r2*w2)); 

		
		
		
		warunki[i] += (-1)*((global_c*global_ro*tablica_elementow[i].r_delta)/global_tau_delta)*( (ni1*tablica_wezlow[i].temp + 
			nj1*tablica_wezlow[i+1].temp)*(ni1*r1*w1) + (ni2*tablica_wezlow[i].temp + 
			nj2*tablica_wezlow[i+1].temp)*(ni2*r2*w2) ); //f1

		warunki[i+1] += (-1)*((global_c*global_ro*tablica_elementow[i].r_delta)/global_tau_delta)*( (ni1*tablica_wezlow[i].temp + 
			nj1*tablica_wezlow[i+1].temp)*(nj1*r1*w1) + (ni2*tablica_wezlow[i].temp + 
			nj2*tablica_wezlow[i+1].temp)*(nj2*r2*w2) ); //f2

		
		
		
		

	}

	warunki[liczba_wezlow-1] -= 2*global_alfa*global_r_max*global_temp_otoczenia;	
	

}


void checkElement(element *tablica_elementow)
{
	cout<<"//////"<<endl;
	cout<<"sprawdzam elementy "<<endl;
	cout<<"//////"<<endl;

	for(int i = 0; i<liczba_elementow; i++)
	{
		cout<<"nod1: "<<tablica_elementow[i].nod1<<endl;
		cout<<"nod2: "<<tablica_elementow[i].nod2<<endl;
		cout<<"k element: "<<tablica_elementow[i].k_element<<endl;
		cout<<"r delta: "<<tablica_elementow[i].r_delta<<endl;
		cout<<endl;
	}
}

void checkWezel(wezel *tablica_wezlow)
{
	cout<<"//////"<<endl;
	cout<<"sprawdzam wezly"<<endl;
	cout<<"//////"<<endl;

	for(int i = 0; i<liczba_wezlow; i++)
	{
		cout<<"temp: "<<tablica_wezlow[i].temp<<endl;
		cout<<"warunek brzegowy: "<<tablica_wezlow[i].warunek_brzegowy<<endl;
		cout<<"r: "<<tablica_wezlow[i].r<<endl;
		cout<<endl;

	}
}

void checkData()
{
	cout<<"//////"<<endl;
	cout<<"liczba elementow: "<<liczba_elementow<<endl;
	cout<<"liczba wezlow: "<<liczba_wezlow<<endl;
	cout<<"global_k: "<<global_k<<endl;
	cout<<"temperatura otoczenia: "<<global_temp_otoczenia<<endl;
	cout<<"global alfa: "<<global_alfa<<endl;
	cout<<"global r_max: "<<global_r_max<<endl;
	cout<<"global r_delta: "<<global_r_delta<<endl;
	cout<<"global tau_delta: "<<global_tau_delta<<endl;
	cout<<"global c: "<<global_c<<endl;
	cout<<"global ro: "<<global_ro<<endl;
	cout<<"global temp_pocz: "<<global_temp_pocz<<endl;
	cout<<"//////"<<endl;
}

void generateFEM_grid(element *tablica_elementow, wezel *tablica_wezlow, long double *warunki)
{
	generateLocalMatrixForElement(tablica_elementow, tablica_wezlow, warunki);
}

void calculateMatrixH(long double *matrix_h[], element *tablica_elementow, long double *warunki)
{
	
	for(int i = 0; i<liczba_elementow+1; i++)
	{
		for(int j = 0; j<liczba_elementow+1; j++)
		{
			matrix_h[i][j] = 0;
		}
	}

	for(int k = 0; k<liczba_elementow; k++)
	{
		for(int i = 0; i<2; i++)
		{
			for(int j = 0; j<2; j++)
			{
				matrix_h[i+k][j+k] += tablica_elementow[k].h_lokalne[i][j];
								
			}			
		}
		
		
	}
	
	matrix_h[liczba_elementow][liczba_elementow] += 2*global_alfa*global_r_max;   

		
	

}

void calculate(long double *matrix_h[], element *tablica_elementow, wezel *tablica_wezlow, long double *warunki)
{
	long double **gauss = new long double *[liczba_elementow+1];
	for(int i = 0; i<liczba_elementow+1; i++)
	{
		gauss[i] = new long double[liczba_elementow+2];
	}
	

	for(int i = 0; i<liczba_elementow+1; i++)
	{
		for(int j = 0; j<liczba_elementow+1; j++)
		{
			gauss[i][j] = matrix_h[i][j];
		}
	}
	
	for(int i = 0; i<liczba_elementow+1; i++)
	{
		 gauss[i][liczba_elementow+1] = (-1)*warunki[i];
		
	}

	//trojkat
	for(int i = 1; i<liczba_wezlow; i++)
	{
		for(int k=i; k<liczba_wezlow; k++)
		{
			for(int j=i; j<liczba_wezlow+1; j++)
			{
				gauss[k][j] = gauss[k][j] - (gauss[i-1][j]*(gauss[k][i-1]/gauss[i-1][i-1]));
			}
			gauss[k][i-1]=0;
		}
	}


	tablica_wezlow[liczba_wezlow-1].temp = gauss[liczba_wezlow-1][liczba_wezlow] / gauss[liczba_wezlow-1][liczba_wezlow-1];
	
	for(int i = liczba_wezlow-2; i>=0; i--)
	{
		tablica_wezlow[i].temp = (gauss[i][liczba_wezlow] - (gauss[i][i+1] * tablica_wezlow[i+1].temp))/gauss[i][i];


	}

	if(wybor=='k')
	{
		cout<<"Czas: "<<czas<<" s Temperatura wewnatrz: "<<tablica_wezlow[0].temp<<" K Temperatura powierzchni: "<<
			tablica_wezlow[liczba_wezlow-1].temp<<" K"<<endl;
	}
	if(wybor=='c')
	{
		cout<<"Czas: "<<czas<<" s Temperatura wewnatrz: "<<tablica_wezlow[0].temp -273.15<<" C Temperatura powierzchni: "<<
			tablica_wezlow[liczba_wezlow-1].temp -273.15<<" C"<<endl;
	}
		
	

	
}



int main()
{
	
	

	if(!loadData("dane.txt")) cout<<"Blad przy otwieraniu pliku"<<endl;
	checkData();
	

	element *tablica_elementow = new element[liczba_elementow];
	wezel *tablica_wezlow = new wezel[liczba_elementow+1];

	
	long double *warunki = new long double[liczba_wezlow];

	long double **matrix_h = new long double * [liczba_elementow+1];
	for(int i = 0; i<liczba_elementow+1; i++)
	{
		matrix_h[i] =new long double[liczba_elementow+1];
	}

	if(!loadData_wezel("dane.txt", tablica_wezlow)) cout<<"Blad przy otwieraniu pliku"<<endl;
	checkWezel(tablica_wezlow);

	if(!loadData_element("dane.txt", tablica_elementow, tablica_wezlow)) cout<<"Blad przy otwieraniu pliku"<<endl;
	checkElement(tablica_elementow);

	
	double tau = 0;
	//global_tau_delta = 0;

	cout<<"Stopnie: Kelvin(k) / Celsjusz(c): "<<endl;
	cin>>wybor;

	cout<<"Podaj czas testu (w sekundach): "<<endl;
	cin>>tau;

	if(wybor=='k')
	{
		cout<<"Czas: "<<0<<" s Temperatura wewnatrz: "<<tablica_wezlow[0].temp<<" K Temperatura powierzchni: "<<
			tablica_wezlow[liczba_wezlow-1].temp<<" K"<<endl;
	}
	if(wybor=='c')
	{
		cout<<"Czas: "<<0<<" s Temperatura wewnatrz: "<<tablica_wezlow[0].temp -273.15<<" C Temperatura powierzchni: "<<
			tablica_wezlow[liczba_wezlow-1].temp -273.15<<" C"<<endl;
	}

	plik2.open("wynik_wnetrze.txt");
	plik3.open("wynik_powierzchnia.txt");

	int tmp1, tmp2;

	
	for(int i = 1; i<=tau; i++)
	{
		czas = i;
		generateFEM_grid(tablica_elementow, tablica_wezlow, warunki);	
		calculateMatrixH(matrix_h, tablica_elementow, warunki);
		calculate(matrix_h, tablica_elementow, tablica_wezlow, warunki);

		tmp1 = tablica_wezlow[0].temp -273.15;
		tmp2 = tablica_wezlow[liczba_wezlow-1].temp - 273.15;
		plik2 << tmp1 << " " << i << endl;
		plik3 << tmp2 << " " << i << endl;

		//plik2 << tmp1 << endl;
		//plik3 << tmp2 << endl;
				
		
	}

	

	system("pause");


}