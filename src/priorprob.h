double Constpriorprob(int p, int dimensioni);
double SBpriorprob(int p, int dimensioni);
double SBConstpriorprob(gsl_vector * indexfr, gsl_vector * positionsx, gsl_matrix * positions, int nofvars, gsl_vector * levels,
	                   int p, gsl_vector * isfactor);
double SBSBpriorprob(gsl_vector * indexfr, gsl_vector * positionsx, gsl_matrix * positions, int nofvars, gsl_vector * levels,
	                   int p, gsl_vector * isfactor);
double ConstConstpriorprob(gsl_vector * indexfr, gsl_vector * positionsx, gsl_matrix * positions, int nofvars, gsl_vector * levels,
										 int p, gsl_vector * isfactor);

