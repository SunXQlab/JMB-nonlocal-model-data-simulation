#include<iostream>
#include<math.h>
#include<stdlib.h>
#include <string.h>
#include<stdio.h>
#include<fstream>
using namespace std;


int main()
{
	double alpha, beta, nu, gamma, ksi, kappa, eta, lambda, phi, V_u;
	double  tau, p, q, G, h1, h3, m1, m3, m4, * x, * y;
	double *** C1, *** C2, *** C11, *** C12, ** U1, ** U2, ** U11, ** D1, ** D2, ** X1, ** X2;
	double* a1, * b1, * c1, * d1, * a2, * b2, * c2, * d2, * a3, * b3, * c3, * d3,  * ans;
	int P, Q,  N, i, j, k, n;
	double temp;
	int pos[2][4000]={0};    
	double initial(double x);
	double sum(double a[], int l, int r);
	double product(double a1[], double a2[], int l, int r);
	double max(double** a, double** b, int m, int n);
	double* chase_algorithm(double* a, double* b, double* c, double* d, int n);


	//parameter values
	alpha = 1e-13;
	beta = 1e-13;
	nu = 2e-5;
	gamma = 1e-4;
	ksi = 1.5e-7;
	kappa = 2e-13;
	eta = 2e-11;
	lambda = 2.78e-6;
	phi = 1e-5;
	V_u = 6.3996e-7;

	N = 3e7;  //number of iterations 
	tau = 0.5;   // time step
	P = 100; 
	Q = 50;  
	p = 0.5 / P;   // spatial step size
	q = 1.0 / Q;  //phenotypic step size
	h1 = beta * tau / (2.0 * pow(p, 2));
	h3 = alpha * tau / (2.0 * pow(q, 2));
	m1 = nu * tau / (2.0 * pow(p, 2));
	m3 = lambda * tau;

	double mal[P+1][P+1];  
	double phe[P+1][P+1];  
	int vas[P+1][P+1]={0}; 
	
	//tridiagonal matrix and RHS to solve C^{l+1/3} and C^{l+2/3}
	a1 = (double*)malloc(sizeof(double) * (P));
	b1 = (double*)malloc(sizeof(double) * (P + 1));
	c1 = (double*)malloc(sizeof(double) * (P));
	d1 = (double*)malloc(sizeof(double) * (P + 1));
	for (i = 1; i <= P - 1; i++)
	{
		a1[i - 1] = -h1;  
		b1[i] = 1 + 2*h1;  
		c1[i] = -h1;  
	}
	a1[P - 1] = -2*h1;
	b1[0] = 1+2*h1;
	b1[P] = 1+2*h1;
	c1[0] = -2*h1;
    
	//tridiagonal matrix and RHS to solve C^{l+1}
	a2 = (double*)malloc(sizeof(double) * (Q));
	b2 = (double*)malloc(sizeof(double) * (Q + 1));
	c2 = (double*)malloc(sizeof(double) * (Q));
	d2 = (double*)malloc(sizeof(double) * (Q + 1));
	for (i = 1; i <= Q - 1; i++)
	{
		a2[i - 1] = -h3;
		b2[i] = 1+2*h3;
		c2[i] = -h3;
	}
	a2[Q - 1] = -2*h3;
	b2[0] = 1+2*h3;;
	b2[Q] = 1+2*h3;
	c2[0] = -2*h3;


	//tridiagonal matrix and RHS to solve U^{l+1/2} and U^{l+1}
	a3 = (double*)malloc(sizeof(double) * (P));
	b3 = (double*)malloc(sizeof(double) * (P + 1));
	c3 = (double*)malloc(sizeof(double) * (P));
	d3 = (double*)malloc(sizeof(double) * (P + 1));
	for (i = 1; i <= P - 1; i++)
	{
		a3[i - 1] = -m1;
		b3[i] = 1 + 2 * m1;
		c3[i] = -m1;
	}
	a3[P - 1] = -2 * m1;
	b3[0] = 1 + 2 * m1;
	b3[P] = 1 + 2 * m1;
	c3[0] = -2 * m1;

	x = (double*)malloc(sizeof(double) * (Q + 1));  //discretized phenotypic state x
	y = (double*)malloc(sizeof(double) * (Q + 1));  //1-x.^2
	for (i = 0; i <= Q; i++)
	{
		x[i] = i * q;
		y[i] = 1 - pow(x[i], 2);
	}

	//phenotypic distribution of tumor
	C1 = (double***)malloc(sizeof(double*) * ((P + 1) * (P + 1) * (Q + 1)));  //C^{l}
	C2 = (double***)malloc(sizeof(double*) * ((P + 1) * (P + 1) * (Q + 1)));  //C^{l+1}
	C11 = (double***)malloc(sizeof(double*) * ((P + 1) * (P + 1) * (Q + 1))); //C^{l+1/3}
	C12 = (double***)malloc(sizeof(double*) * ((P + 1) * (P + 1) * (Q + 1)));  //C^{l+2/3}
	for (i = 0; i <= P; i++)
	{
		C1[i] = (double**)malloc(sizeof(double*) * ((P + 1) * (Q + 1)));
		C2[i] = (double**)malloc(sizeof(double*) * ((P + 1) * (Q + 1)));
		C11[i] = (double**)malloc(sizeof(double*) * ((P + 1) * (Q + 1)));
		C12[i] = (double**)malloc(sizeof(double*) * ((P + 1) * (Q + 1)));
	}
	for (i = 0; i <= P; i++)
	{
		for (j = 0; j <= P; j++)
		{
			C1[i][j] = (double*)malloc(sizeof(double) * (Q + 1));
			C2[i][j] = (double*)malloc(sizeof(double) * (Q + 1));
			C11[i][j] = (double*)malloc(sizeof(double) * (Q + 1));
			C12[i][j] = (double*)malloc(sizeof(double) * (Q + 1));
		}
	}


	//oxygen concentration
	U1 = (double**)malloc(sizeof(double*) * (P + 1));
	U2 = (double**)malloc(sizeof(double*) * (P + 1));
	U11 = (double**)malloc(sizeof(double*) * (P + 1));
	for (i = 0; i <= P; i++)
	{
		U1[i] = (double*)malloc(sizeof(double) * (P + 1));
		U2[i] = (double*)malloc(sizeof(double) * (P + 1));
		U11[i] = (double*)malloc(sizeof(double) * (P + 1));
	}

	//tumor cell density and mean phenotypic state
	D1 = (double**)malloc(sizeof(double*) * (P + 1));
	D2 = (double**)malloc(sizeof(double*) * (P + 1));
	X1 = (double**)malloc(sizeof(double*) * (P + 1));
	X2 = (double**)malloc(sizeof(double*) * (P + 1));
	for (i = 0; i <= P; i++)
	{
		D1[i] = (double*)malloc(sizeof(double) * (P + 1));
		D2[i] = (double*)malloc(sizeof(double) * (P + 1));
		X1[i] = (double*)malloc(sizeof(double) * (P + 1));
		X2[i] = (double*)malloc(sizeof(double) * (P + 1));
	}

	//initial conditions
	FILE* fp = fopen("vascular_pos.csv", "r");
	i = 0;
    char line1[10000];
	char *tok;
	while (fgets(line1, 10000, fp)!=NULL)
    {
    	j=0;
        char* tmp = strdup(line1);
        for (tok = strtok(line1, ","); tok && *tok; j++, tok = strtok(NULL, ",\n")){
        	pos[i][j] = atoi(tok);
		}
        i++;
        free(tmp);
    }
    fclose(fp);

	int col = j-1;

	for (i = 0; i <= P; i++)    
		for (j = 0; j <= P; j++)
		{
			U1[i][j]=0;
			vas[i][j]=0;
			for(k=0;k<=col;k++)
			{
				if (j == pos[0][k]&& i ==pos[1][k])
				{
					U1[i][j] = V_u;
					vas[i][j]=1;
					break;
				}
			}
		}


	fp = fopen("malignant.csv", "r");
	i=0;
	char line2[10000];
	while (fgets(line2, 10000, fp)!=NULL)
    {
    	j=0;
        char* tmp = strdup(line2);
        for (tok = strtok(line2, ","); tok && *tok; j++, tok = strtok(NULL, ",\n")){
        	mal[i][j] = atof(tok);
		}
        i++;
        free(tmp);
	}

	fp = fopen("phenotype.csv", "r");
	i=0;
	char line3[10000];
	while (fgets(line3, 10000, fp)!=NULL)
    {
    	j=0;
        char* tmp = strdup(line3);
        for (tok = strtok(line3, ","); tok && *tok; j++, tok = strtok(NULL, ",\n")){
        	phe[i][j] = atof(tok);//转换成整数
		}//字符串拆分操作 
        i++;
        free(tmp);
	}
	for(i=0;i<=P;i++)
	{
		for(j=0;j<=P;j++)
		{
                if(phe[i][j]==0)
                {
                    for(k=1;k<=Q;k++)
                        C1[i][j][k]=0;
                    C1[i][j][0] = mal[i][j];
                }
				else if(phe[i][j]==1)
				{
					for(k=0;k<Q;k++)
                        C1[i][j][k] = 0;
					C1[i][j][Q] = mal[i][j];
				}

                else if(phe[i][j]<0.5)
                {
                    for(k=0;k<=Q;k++)
                        C1[i][j][k]=(x[k]<= 2 * phe[i][j]) * 1e8 * exp(-pow(x[k] - phe[i][j], 2) / 0.1);
                }
                else
                {
                    for(k=0;k<=Q;k++)
                        C1[i][j][k]=(x[k]>= 2 * phe[i][j]-1)*1e8 * exp(-pow(x[k] - phe[i][j], 2) / 0.1);
                }
		}
	}

    for(i=0;i<=P;i++)
    {
        for(j=0;j<=P;j++)
        {
            if(mal[i][j]>0)
            {
                temp = sum(C1[i][j],0,Q);
                for(k=0;k<=Q;k++)
                {
                    C1[i][j][k]=C1[i][j][k] / temp * mal[i][j];
                }
            }
        }
    }
	

	//start of iteration	
	for (n = 0; n <= N; n++)
	{
		//calculate for C
		for (j = 0; j <= P; j++)
		{
			for (k = 0; k <= Q; k++)
			{
				for (i = 0; i <= P; i++)
				{
					G = phi * (1 - pow(1 - x[k], 2)) + gamma * U1[i][j] / (ksi + U1[i][j]) * (1 - pow(x[k], 2)) - kappa * q/2 * (2*sum(C1[i][j], 1, Q-1) + C1[i][j][0] + C1[i][j][Q]);
					d1[i] = h1 * h1 * h3 * C1[(i!=0)*(i - 1)+(i==0)][(j!=0)*(j - 1)+(j==0)][(k!=0)*(k - 1)+(k==0)] + (h1 * h1 - 2 * h1 * h1 * h3) * C1[(i!=0)*(i - 1)+(i==0)][(j!=0)*(j - 1)+(j==0)][k] + h1 * h1 * h3 * C1[(i!=0)*(i - 1)+(i==0)][(j!=0)*(j - 1)+(j==0)][(k!=Q)*(k + 1)+(k==Q)*(Q-1)]
							+ (h1 * h3 - 2 * h1 * h1 * h3) * C1[(i != 0) * (i - 1)+(i==0)][j][(k != 0) * (k - 1)+(k==0)] + (h1 - 2 * h1 * h1 - 2 * h1 * h3 + 4 * h1 * h1 * h3) * C1[(i != 0) * (i - 1)+(i==0)][j][k] + (h1 * h3 - 2 * h1 * h1 * h3) * C1[(i != 0) * (i - 1)+(i==0)][j][(k != Q) * (k + 1) + (k == Q) * (Q - 1)]
							+h1 * h1 * h3 * C1[(i != 0) * (i - 1)+(i==0)][(j!=P)*(j + 1)+(j==P)*(P-1)][(k != 0) * (k - 1)+(k==0)] + (h1 * h1 - 2 * h1 * h1 * h3) * C1[(i != 0) * (i - 1)+(i==0)][(j != P) * (j + 1) + (j == P) * (P - 1)][k] + h1 * h1 * h3 * C1[(i != 0) * (i - 1)+(i==0)][(j != P) * (j + 1) + (j == P) * (P - 1)][(k != Q) * (k + 1) + (k == Q) * (Q - 1)]
							+ (h1 * h3 - 2 * h1 * h1 * h3) * C1[i][(j != 0) * (j - 1)+(j==0)][(k != 0) * (k - 1)+(k==0)] + (h1 - 2 * h1 * h1 - 2 * h1 * h3 + 4 * h1 * h1 * h3) * C1[i][(j != 0) * (j - 1)+(j==0)][k] + (h1 * h3 - 2 * h1 * h1 * h3) * C1[i][(j != 0) * (j - 1)+(j==0)][(k != Q) * (k + 1) + (k == Q) * (Q - 1)]
							+ (h3 - 2 * h1 * h3 - 2 * h1 * h3 + 4 * h1 * h1 * h3) * C1[i][j][(k != 0) * (k - 1)+(k==0)] + (1 - 4 * h1 - 2 * h3 + 4 * h1 * h1 + 8 * h1 * h3 - 8 * h1 * h1 * h3) * C1[i][j][k] + (h3 - 4 * h1 * h3 + 4 * h1 * h1 * h3) * C1[i][j][(k != Q) * (k + 1) + (k == Q) * (Q - 1)]
							+ (h1 * h3 - 2 * h1 * h1 * h3) * C1[i][(j != P) * (j + 1) + (j == P) * (P - 1)][(k != 0) * (k - 1)+(k==0)] + (h1 - 2 * h1 * h1 - 2 * h1 * h3 + 4 * h1 * h1 * h3) * C1[i][(j != P) * (j + 1) + (j == P) * (P - 1)][k] + (h1 * h3 - 2 * h1 * h1 * h3) * C1[i][(j != P) * (j + 1) + (j == P) * (P - 1)][(k != Q) * (k + 1) + (k == Q) * (Q - 1)]
							+ h1 * h1 * h3 * C1[(i!=P)*(i + 1)+(i==P)*(P-1)][(j != 0) * (j - 1)+(j==0)][(k != 0) * (k - 1)+(k==0)] + (h1 * h1 - 2 * h1 * h1 * h3) * C1[(i != P) * (i + 1) + (i == P) * (P - 1)][(j != 0) * (j - 1)+(j==0)][k] + h1 * h1 * h3 * C1[(i != P) * (i + 1) + (i == P) * (P - 1)][(j != 0) * (j - 1)+(j==0)][(k != Q) * (k + 1) + (k == Q) * (Q - 1)]
							+ (h1 * h3 - 2 * h1 * h1 * h3) * C1[(i != P) * (i + 1) + (i == P) * (P - 1)][j][(k != 0) * (k - 1)+(k==0)] + (h1 - 2 * h1 * h1 - 2 * h1 * h3 + 4 * h1 * h1 * h3) * C1[(i != P) * (i + 1) + (i == P) * (P - 1)][j][k] + (h1 * h3 - 2 * h1 * h1 * h3) * C1[(i != P) * (i + 1) + (i == P) * (P - 1)][j][(k != Q) * (k + 1) + (k == Q) * (Q - 1)]
							+ h1 * h1 * h3 * C1[(i != P) * (i + 1) + (i == P) * (P - 1)][(j != P) * (j + 1) + (j == P) * (P - 1)][(k != 0) * (k - 1)+(k==0)] + (h1 * h1 - 2 * h1 * h1 * h3) * C1[(i != P) * (i + 1) + (i == P) * (P - 1)][(j != P) * (j + 1) + (j == P) * (P - 1)][k] + h1 * h1 * h3 * C1[(i != P) * (i + 1) + (i == P) * (P - 1)][(j != P) * (j + 1) + (j == P) * (P - 1)][(k != Q) * (k + 1) + (k == Q) * (Q - 1)]
							+ tau * G * C1[i][j][k];
				}
				ans = chase_algorithm(a1, b1, c1, d1, P + 1);
				for (i = 0; i <= P; i++)
				{
					C11[i][j][k] = ans[i];    
				}
				free(ans);
			}
		}  
	
		for (k = 0; k <= Q; k++)
		{
			for (i = 0; i <= P; i++)
			{
				for (j = 0; j <= P; j++)
					d1[j] = C11[i][j][k];
				ans = chase_algorithm(a1, b1, c1, d1, P + 1);
				for (j = 0; j <= P; j++)
				{
					C12[i][j][k] = ans[j];   
				}
				free(ans);
			}
		}

		for (j = 0; j <= P; j++)
		{
			for (i = 0; i <= P; i++)
			{
				for (k = 0; k <= Q; k++)
					d2[k] = C12[i][j][k];
				ans = chase_algorithm(a2, b2, c2, d2, Q + 1);
				for (k = 0; k <= Q; k++)
				{
					C2[i][j][k] = ans[k];    
				}
				free(ans);
			}
		}
		
		//calculate for U
		for (j = 0; j <= P; j++)
		{
			for (i =0; i <= P; i++)
			{
				m4= tau * eta * gamma * q/2 *(2 * product(y, C1[i][j], 1, Q-1) + y[0]*C1[i][j][0] + y[Q]*C1[i][j][Q]);
				d3[i] = m1 * m1 * U1[(i!=0)*(i - 1)+(i==0)][(j!=0)*(j - 1)+(j==0)] + m1 * (1 - 2 * m1) * U1[(i != 0) * (i - 1)+(i==0)][j] + m1 * m1 * U1[(i != 0) * (i - 1)+(i==0)][(j!=P)*(j + 1)+(j==P)*(P-1)]
					+ (m1 - 2 * m1 * m1) * U1[i][(j != 0) * (j - 1)+(j==0)] + (1 - 4 * m1 - m3 + 4 * m1 * m1 ) * U1[i][j] + (m1 - 2 * m1 * m1) * U1[i][(j != P) * (j + 1) + (j == P) * (P - 1)]
					+ m1 * m1 * U1[(i!=P)*(i + 1)+(i==P)*(P-1)][(j != 0) * (j - 1)+(j==0)] + m1 * (1 - 2 * m1 ) * U1[(i != P) * (i + 1) + (i == P) * (P - 1)][j] + m1 * m1 * U1[(i != P) * (i + 1) + (i == P) * (P - 1)][(j != P) * (j + 1) + (j == P) * (P - 1)]
					-m4 * U1[i][j] / (ksi + U1[i][j])
					+ V_u * tau * vas[i][j];
			}
			ans = chase_algorithm(a3, b3, c3, d3, P + 1);
			for (i = 0; i <= P; i++)
			{
				U11[i][j] = ans[i];    
			}
			free(ans);
		}

		for (i = 0; i <= P; i++)
		{
			for (j = 0; j <= P; j++)
			{
				d3[j] = U11[i][j];
			}
			ans = chase_algorithm(a3, b3, c3, d3, P + 1);
			for (j = 0; j <= P; j++)
			{
				U2[i][j] = ans[j];   
			}
			free(ans);
		}

		//save data
		if (n % 100 == 0)
		{
			ofstream outFile;
			outFile.open("tumor.csv", ios::out);
			outFile << "iteration number" << ',' << n  << endl;
			for (int i = 0; i <= P; i++)
			{
				for (int j = 0; j <= P; j++)
				{
					for (int k = 0; k <= Q; k++)
						outFile << C1[i][j][k] << ',';
				}
				outFile << endl;
			}
			outFile.close();
			outFile.open("oxygen.csv", ios::out);
			outFile << "iteration number" << ',' << n  << endl;
			for (int i = 0; i <= P; i++)
			{
				for (int j = 0; j <= P; j++)

					outFile << U1[i][j] << ',';
				outFile << endl;
			}
			outFile.close();
		}

		//threshold check
		float max_diff1 = 0, max_diff2 = 0, max_diff3 = 0;
		for (i = 0; i <= P; i++)
			for (j = 0; j <= P; j++)
			{
				D1[i][j] = q/2 * (2*sum(C1[i][j], 1, Q-1) + C1[i][j][0] + C1[i][j][Q]);
				D2[i][j] = q/2 * (2*sum(C2[i][j], 1, Q-1) + C2[i][j][0] + C2[i][j][Q]);
				X1[i][j] = q/2 * (2*product(x, C1[i][j], 1, Q-1) + x[0] * C1[i][j][0] + x[Q] * C1[i][j][Q]) /D1[i][j];
				X2[i][j] = q/2 * (2*product(x, C2[i][j], 1, Q-1) + x[0] * C2[i][j][0] + x[Q] * C2[i][j][Q]) /D2[i][j];
				max_diff1 = max_diff1 > fabs ( U1[i][j]-U2[i][j] ) ? max_diff1 : fabs ( U1[i][j] - U2[i][j] );   
				max_diff2 = max_diff2 > fabs ( D1[i][j]-D2[i][j] ) ? max_diff2 : fabs ( D1[i][j] - D2[i][j] );   
				max_diff3 = max_diff3 > fabs ( X1[i][j]-X2[i][j] ) ? max_diff3 : fabs ( X1[i][j] - X2[i][j] );  
			}
		if (max_diff1 < 1.5e-18 && max_diff2 < 7.5e-4  && max_diff3 < 1.5e-10 ) {
        printf("Early stop");
        break;
    	}


		for (i = 0; i <= P; i++)
			for (j = 0; j <= P; j++)
			{
				U1[i][j] = U2[i][j];
			}
		for (i = 0; i <= P; i++)
			for (j = 0; j <= P; j++)
				for (k = 0; k <= Q; k++)
				{
					C1[i][j][k] = C2[i][j][k];
				}
        
	}

	free(x); free(y);
	free(a1); free(b1); free(c1); free(d1);
	free(a2); free(b2); free(c2); free(d2);
	free(a3); free(b3); free(c3); free(d3);
	for (i = 0; i <= P; i++)
	{
		for (j = 0; j <= P; j++)
		{
			free(C1[i][j]); free(C2[i][j]);
			free(C11[i][j]); free(C12[i][j]);
		}
	}
	for (i = 0; i <= P; i++)
	{
		free(U1[i]); free(U2[i]);  free(U11[i]);
		free(C1[i]); free(C2[i]);  free(C11[i]); free(C12[i]);
	}
	return 0;
}




double initial(double x)
{
	return 1e8 * exp(-pow(x - 0.5, 2) / 0.1);
}


double sum(double a[], int l, int r)
{
	double s = 0;
	int i;
	for (i = l; i <= r; i++)
		s += a[i];
	return s;
}


double product(double a1[], double a2[], int l, int r)
{
	double p = 0;
	for (int i = l; i <= r; i++)
	{
		p += a1[i] * a2[i];
	}
	return p;
}


double* chase_algorithm(double* a, double* b, double* c, double* d, int n)
{
	double* x = NULL, * y = NULL, * alpha = NULL, * beta = NULL;
	if (n <= 3)
		exit(-1);
	int i;
	y = (double*)malloc(sizeof(double) * n);
	x = (double*)malloc(sizeof(double) * n);
	alpha = (double*)malloc(sizeof(double) * n);
	beta = (double*)malloc(sizeof(double) * (n - 1));
	if (x == NULL || y == NULL || alpha == NULL || beta == NULL)
		exit(-1);
	alpha[0] = b[0]; beta[0] = c[0] / b[0];
	for (i = 1; i <= n - 2; i++)
	{
		alpha[i] = b[i] - a[i - 1] * beta[i - 1];
		beta[i] = c[i] / alpha[i];
	}
	alpha[n - 1] = b[n - 1] - a[n - 2] * beta[n - 2];
	y[0] = d[0] / b[0];
	for (i = 1; i <= n - 1; i++)
	{
		y[i] = (d[i] - a[i - 1] * y[i - 1]) / alpha[i];
	}
	x[n - 1] = y[n - 1];
	for (i = n - 2; i >= 0; i--)
	{
		x[i] = y[i] - beta[i] * x[i + 1];
	}

	free(alpha);
	free(beta);
	free(y);
	return x;
}

