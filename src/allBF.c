#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>


// Gonzalo Garcia-Donato and Anabel Forte, 2010
// The Bayes factor obtained with the the g-prior (in favour of 
// a model and against the null with only the intercept)

double gBF21fun(int n, int k2, double Q)
{
	
		double BF21=0.0;
    BF21 = exp(((n-k2)/2.0)*log(1.0+n)-((n-1.0)/2.0)*log(1.0+n*Q));
		return BF21;
	
}



/*-------------Robust--------------------*/
//Sirve para calcular con o sin tess, cuando es sin tess el ultimo objeto tiene que ser n.


/*----------funcion auxiliar para la hypergeo cuando z esta muy cerca de 1------*/

/*Structure for Parameters*/
struct parrob {
	double a;
	double b;
	double c;
	double z;
};


/*auxiliar functions for integration */

double robint_aux (double x, void *p){
	struct parrob * params=(struct parrob *)p;/*Defino un puntero a una estructura del tipo par*/
	/*Inicializo el puntero en la dirección de memoria en la que están los parametros que le estoy
	 pasando, p obligando a que sean de tipo struct par*/
	
	/*Defino los parametros que son los que estaran en la estructura que le pasamos*/
	double a=(params->a);
	double b=(params->b);
	double c=(params->c);
	double z=(params->z);
	
	/*Calculo el valor de la función y lo devuelvo*/
	double l=pow(x,b-1.0)*pow((1.0-x),c-b-1.0)*pow((1.0-x*z),-a);
	//printf("valor de la funcion %f\n",l);
	return l;
}


/*Integrated functions the arguments will be n,k,Qi0*/
double robint (double a,double b, double c,double z){
	/*guardamos espacio de memoria para realizar la integracion, este 10000 es el que luego va en la función
	 de integracion*/
	gsl_integration_workspace * w=gsl_integration_workspace_alloc(10000);
	
	double result=0.0;
	double error=0.0;
	
	/*Ponemos los parametros en la forma que necesitamos para la funcion*/
	struct parrob  params={a,b,c,z};
	
	/*Definimos cual es la funcion que vamos a usar y le pasamos los parametros*/
	gsl_function F;
	F.function = &robint_aux;
	F.params = &params;
	
	/*integramos y guardamos el resultado en result y el error en error*/
	gsl_integration_qags(&F, 0,1, 0, 1e-9,10000,w,&result,&error);
	
	
	/*Liberamos el espacio de trabajo*/
	gsl_integration_workspace_free (w);
	
	/*devolvemos el resultado*/
	return result*gsl_sf_gamma(c)/(gsl_sf_gamma(b)*gsl_sf_gamma(c-b));
}


/*----Caculo de los factores bayes robustos------*/

double RobustBF21fun(int n, int k2, double Q)
{

	double  T1=0.0, T2=0.0, T3=0.0;
	double arg=0.0;
	double R1=0.0;
	double z=0.0;
	double rho=0.0;
  rho=pow((k2+1.0),-1.0);
  
  double k2aux=0.0;
  k2aux=k2+1.0;
  
  double Qaux=0.0;
  Qaux=pow(Q,-1.0);
	
	// Qaux aquí representa Q_0i
	T1=pow(Qaux,(n-1.0)/2.0);
	T2=pow(rho*(n+1),-(k2/2.0))/k2aux;
	//for the hypergeometric factor we disdtinguish whether the argument is
	//<-1 or not
	
	
	arg=(1.0-Qaux)/(rho*(1.0+n));
	gsl_sf_result result;
	int STATUS=0;
	
	
	if (arg>=-1.0)
	{
		T3=gsl_sf_hyperg_2F1(k2aux/2.0, (n-1.0)/2.0, (k2aux/2.0)+1.0, arg);
	}
	else 
	{
		z=arg/(arg-1.0);
		STATUS=gsl_sf_hyperg_2F1_e((n-1.0)/2.0, 1.0, (k2aux/2.0)+1.0, z,&result);
		
		
		if (STATUS==0) T3=pow((1.0-arg),(1.0-n)/2.0)*result.val; //succed
		else //gsl_hyper failed, then numerical approx of the log(2F1)
		{
							
			T3=pow((1.0-arg),(1.0-n)/2.0)*robint((n-1.0)/2.0,1.0, (k2aux/2.0)+1.0, z);
					}
	}
	
	
	R1=T1*T2*T3;
	
	return(R1);
}

/* FUNCION QUE USAREMOS EN EL main.c*/
/* Llamada al factor bayes robusto con el Rho del paper 1/(ki+k0) */
/* NO sera necesaria, la defino arriba más arregladita 
double RobustBF21fun(int n, int k2, double Q)
{
    double RobustBF21=0.0; 
    double rho=0.0; 
    rho = pow((k2+1.0),-1.0);
    RobustBF21 = Robustfun((double)k2, (double)n, Q, rho, n);
    return (RobustBF21);
    
}
*/


/*---------------------LIANG----------------*/
/*Structure for Parameters, for Liang and ZS*/
struct par {
	double n;
	double k_i;
	double Q_i0;
};


double liang_aux (double x, void *p){
	struct par * params=(struct par *)p;/*Defino un puntero a una estructura del tipo par*/
	/*Inicializo el puntero en la dirección de memoria en la que están los parametros que le estoy
	 pasando, p obligando a que sean de tipo struct par*/
	
	/*Defino los parametros que son los que estaran en la estructura que le pasamos*/
	double n=(params->n);
	double k=(params->k_i);
	double Q=(params->Q_i0);
	
	/*Calculo el valor de la función y lo devuelvo*/
	double l=pow((1.0+x), (n-1.0-k)/2.0)*pow((1.0+Q*x), (1.0-n)/2)*(1.0/(2.0*n))*pow(1.0+x/n, -1.5);
	
	return l;
}

/*Esta es igual que la anterior pero para el FB de Liang*/
double liang (double n,double k, double Q){
	gsl_integration_workspace * w=gsl_integration_workspace_alloc(10000);
	
	double result=0.0, error=0.0;
	
	struct par  params={n,k,Q};
	
	gsl_function F;
	F.function = &liang_aux;
	F.params = &params;
	
	gsl_integration_qagiu(&F, 0, 0, 1e-9,10000,w,&result,&error);
	
	gsl_integration_workspace_free (w);
	
	return result;
}



/* FUNCION QUE USAREMOS EN EL main.c*/
double LiangBF21fun(int n, int k2, double Q)
{
    double LiangBF21=0.0;
    LiangBF21 = liang((double) n, (double) k2, Q);
    return LiangBF21;
    
}


/* ---------------JZS------------------*/


/*auxiliar functions for integration */

double zell_aux (double x, void *p){
	struct par * params=(struct par *)p;/*Defino un puntero a una estructura del tipo par*/
	/*Inicializo el puntero en la dirección de memoria en la que están los parametros que le estoy
	 pasando, p obligando a que sean de tipo struct par*/
	
	/*Defino los parametros que son los que estaran en la estructura que le pasamos*/
	double n=(params->n);
	double k=(params->k_i);
	double Q=(params->Q_i0);
	
	/*Calculo el valor de la función y lo devuelvo*/
	double l=pow((1.0+x), (n-1.0-k)/2.0)*pow((1.0+Q*x), (1.0-n)/2)*pow((n/(2.0*M_PI)),0.5)*pow(x, -1.5)*exp(-n/(2.0*x));
	
	return l;
}


/*Integrated functions the arguments will be n,k,Qi0*/
double zell (double n,double k, double Q){
	/*guardamos espacio de memoria para realizar la integracion, este 10000 es el que luego va en la función
	 de integracion*/
	gsl_integration_workspace * w=gsl_integration_workspace_alloc(10000);
	
	double result=0.0, error=0.0;
	
	/*Ponemos los parametros en la forma que necesitamos para la funcion*/
	struct par  params={n,k,Q};
	
	/*Definimos cual es la funcion que vamos a usar y le pasamos los parametros*/
	gsl_function F;
	F.function = &zell_aux;
	F.params = &params;
	
	/*integramos y guardamos el resultado en result y el error en error*/
	gsl_integration_qagiu(&F, 0, 0, 1e-9,10000,w,&result,&error);
	
	//printf("errorJZS=%f \n", error);
	/*Liberamos el espacio de trabajo*/
	gsl_integration_workspace_free (w);
	
	/*devolvemos el resultado*/
	return result;
}



/* FUNCION QUE USAREMOS EN EL main.c*/
double ZSBF21fun(int n, int k2, double Q)
{
    double ZSBF21=0.0;
    ZSBF21 = zell ((double) n,(double) k2, Q);
    return ZSBF21;
    
}
