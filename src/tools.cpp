
#include "tools.h"


/**@brief ascending on first VERTEX_TYPE of pair<VERTEX_TYPE,VERTEX_TYPE*> comparator*/
bool compKeepIndex(pair<VERTEX_TYPE,VERTEX_TYPE*> a, pair<VERTEX_TYPE,VERTEX_TYPE*> b)
{
	return a.first < b.first;
}

///////////////////////////////////////////////////////////////////

bool index_double_comparator(const indexed_double &l, const indexed_double &r)
{ 
	return l.first < r.first; 
}

///////////////////// MISCELLANEAOUS ////////////////////////////////

/** @brief test if point proj is in cube whose center is point and side is side variable.
 If a toll is provided, the test is performed up to a tollerance value.
 By default the tollerance is 0*/
bool testInCube(const double proj0, const double proj1, const double proj2,
				const double point0, const double point1, const double point2,
				const double side, const double toll)
{
	double ss = side*0.5;
	double minx = point0-ss;
	double miny = point1-ss;
	double minz = point2-ss;
	double maxx = point0+ss;
	double maxy = point1+ss;
	double maxz = point2+ss;

	if (proj0>minx-toll && proj0<maxx+toll &&
        proj1>miny-toll && proj1<maxy+toll &&
        proj2>minz-toll && proj2<maxz+toll)
		return true;
	else
		return false;
}

double rintp(const double x)
{
  return floor(x + 0.5);
}


string toLowerCase(string str) 
{
    for (unsigned int i=0; i<strlen(str.c_str()); i++)
		if (str[i] >= 0x41 && str[i] <= 0x5A)
			str[i] = str[i] + 0x20;
    return str;
}

void cleanLine()
{
	printf("\r                                                                          ");
}


// randnumber between 0 and 1
double randnum()
{
	long int aa, mm, qq, rr, hh, lo, test;
	double reslt;
	aa= 16807;
	mm= 2147483647;
	qq= 127773;
	rr= 2836;
	hh= (long int)(SEED/qq);
	lo= SEED - hh*qq;
	test= aa*lo-rr*hh;
	if (test >= 0)
		SEED = test;
	else
		SEED = test + mm;
	reslt = SEED/(double)mm;
	return( reslt );
}


/** get the real roots by computing the companion matrix and then extracting the eigenvalues. A root is real
if its imaginary part in absolute value is less than a given threshold. Usually this threshold is rather conservative
such that possibly unprecise real roots are not lost*/
void getRealRootsCompanion(double *const poly, const int degree, double *const roots, int &numroots)
{	
	int zero_roots = 0, ind=0, ind2=0;
	int numelem = degree;
	int truedegree = degree;
	//int* notNull = allocateVector<int>(degree+1);
	int notNull[MAX_POLY_DEGREE+1];

	// find non zeros elements and record indexes
	for (int i=0; i<degree+1; i++)
	{
		if (poly[i] != 0)
		{
			notNull[ind2] = i;
			ind2++;
		}
	}

	// remove all zeros and record the zeros roots; in place rewrite the true polynomial
	ind = 0;
	for (int i=notNull[0]; i <= notNull[ind2-1]; i++)
	{
		poly[ind] = poly[i];
		ind++;
	}
	numelem = ind;
	truedegree = notNull[ind2-1];
	// number of zero roots
    zero_roots = notNull[0];
	
	//double *localpoly = allocateVector<double>(numelem-1);
	double localpoly[MAX_POLY_DEGREE+1];

	// get the real degree of the polynomial by stripping almost zero leading coefficients
	// almost-zero leading coefficient are those which produce 1/c = INF; this enhances
	// numerical stability in case bad coefficients are supplied
	bool oneinf = true;
	while (oneinf)
	{
		oneinf = false;
		// try to normalize coefficients
		for (int i=0; i<numelem-1; i++)
		{
			localpoly[i] = poly[i]/poly[numelem-1];
			// normalization failed due to 1/c = INF
			if (numeric_limits<double>::infinity() == localpoly[i])
				oneinf = true;
		}		
		if (oneinf)
		{
			numelem--;
			truedegree--;
		}
	}

	numroots = truedegree; // true degree is the final number of valid roots

	//Array2D<double> M(6,6);
	TNT::Array2D<double> M(numelem-1, numelem-1);
	
	// build the companion matrix of the filtered polynomial
	for (int i=0; i<numelem-1; i++)
		for (int j=0; j<numelem-1; j++)
			M[i][j] = 0;

	for (int i=0; i<numelem-1; i++)
		M[i][numelem-2] = -localpoly[i];

	for (int i=1; i<numelem-1; i++)
		for (int j=0; j<numelem-2; j++)
			if (i-1 == j)
				M[i][j] = 1;
	
	JAMA::Eigenvalue<double> eig(M);
	
	// get the roots by the companion matrix, skip eigenvectors computation	
	int i=0;
	int final = numelem-1;
	numroots = 0;

	//Array1D<double> real(6);
	//Array1D<double> img(6);
	TNT::Array1D<double> real(numelem-1);
	TNT::Array1D<double> img(numelem-1);

	eig.getRealEigenvalues(real);
	eig.getImagEigenvalues(img);

	for (i=0; i<final; i++)
	{
		// identify real roots by using a prescribed tollerance
		//printf("\n re %lf img %lf",real[i],img[i]);
		if (fabs(img[i]) < 1e-1)
		{
			roots[numroots]=real[i];
			numroots++;
		}
		else
			numelem--;
	}	

	// add the zero roots
	for (int j=numroots; j < zero_roots + numroots; j++)
	{
		roots[j] = 0;
		numroots++;
	}
}


/** get real roots by using Sturm method. Directly search real roots. Much faster, often less accurate
than companion matrix*/
void getRealRootsSturm(const double *const polyy, const int degree, double *const roots, int &numrootss)
{
	int nchanges, np, atmin, atmax,nroots,i;
	double min,max;

	poly sseq[MAX_POLY_DEGREE+1];

	for (int i=0;i<degree+1;i++)
		sseq[0].coef[i] = polyy[i];

	np = buildsturm(degree, sseq);

	// get the number of real roots	 
	nroots = numroots(np, sseq, &atmin, &atmax);

	if (nroots == 0) 
	{
		numrootss = 0;
		return;
	}
		
	//calculate the bracket that the roots live in
	min = -1.0;
	nchanges = numchanges(np, sseq, min);
	for (i = 0; nchanges != atmin && i != MAXPOW; i++) 
	{ 
		min *= 10.0;
		nchanges = numchanges(np, sseq, min);
	}

	if (nchanges != atmin) 
	{
		printf("solve: unable to bracket all negative roots\n");
		atmin = nchanges;
	}

	max = 1.0;
	nchanges = numchanges(np, sseq, max);
	for (i = 0; nchanges != atmax && i != MAXPOW; i++) 
	{ 
		max *= 10.0;
		nchanges = numchanges(np, sseq, max);
	}

	if (nchanges != atmax) 
	{
		printf("solve: unable to bracket all positive roots\n");
		atmax = nchanges;
	}
	nroots = atmin - atmax;

	// perform the bisection.
	sbisect(np, sseq, min, max, atmin, atmax, roots);

	numrootss = nroots;
}


/** plane by 3 points routine*/
void plane3points(const double p1[3], const double p2[3], const double p3[3], double w[4], const bool normalize)
{
	w[0] = p1[1]*(p2[2] - p3[2]) + p2[1]*(p3[2] - p1[2]) + p3[1]*(p1[2] - p2[2]);
	w[1] = p1[2]*(p2[0] - p3[0]) + p2[2]*(p3[0] - p1[0]) + p3[2]*(p1[0] - p2[0]);
	w[2] = p1[0]*(p2[1] - p3[1]) + p2[0]*(p3[1] - p1[1]) + p3[0]*(p1[1] - p2[1]);
	w[3] = -DOT(w,p1);

	if (normalize)
	{
		double norm = sqrt(DOT(w,w));
		w[0] /= norm;
		w[1] /= norm;
		w[2] /= norm;
		w[3] /= norm;
	}
	return;
}

/** point to plane projection*/
void point2plane(const double p[3], double w[4], double *const dist, double proj[3])
{
	double den = (DOT(w,w));
	double d = sqrt(den);
	double val = DOT(w,p)+w[3];
	double c = (val/(den));
	proj[0] = p[0]-w[0]*c;
	proj[1] = p[1]-w[1]*c;
	proj[2] = p[2]-w[2]*c;
	(*dist) = fabs(val/d);
}

/** in place inversion of 4x4 matrix*/
void inplace_invert4x4(double M[4][4])
{
	double a00 = M[0][0]; double a01 = M[0][1]; double a02 = M[0][2]; double a03 = M[0][3];
	double a10 = M[1][0]; double a11 = M[1][1]; double a12 = M[1][2]; double a13 = M[1][3];
	double a20 = M[2][0]; double a21 = M[2][1]; double a22 = M[2][2]; double a23 = M[2][3];
	double a30 = M[3][0]; double a31 = M[3][1]; double a32 = M[3][2]; double a33 = M[3][3];
	
	M[0][0] =      a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22;
	M[0][1] =    - a01*a22*a33 + a01*a23*a32 + a21*a02*a33 - a21*a03*a32 - a31*a02*a23 + a31*a03*a22;
	M[0][2] =      a01*a12*a33 - a01*a13*a32 - a11*a02*a33 + a11*a03*a32 + a31*a02*a13 - a31*a03*a12;
	M[0][3] =    - a01*a12*a23 + a01*a13*a22 + a11*a02*a23 - a11*a03*a22 - a21*a02*a13 + a21*a03*a12;
	M[1][0] =    - a10*a22*a33 + a10*a23*a32 + a20*a12*a33 - a20*a13*a32 - a30*a12*a23 + a30*a13*a22;
	M[1][1] =      a00*a22*a33 - a00*a23*a32 - a20*a02*a33 + a20*a03*a32 + a30*a02*a23 - a30*a03*a22;
	M[1][2] =    - a00*a12*a33 + a00*a13*a32 + a10*a02*a33 - a10*a03*a32 - a30*a02*a13 + a30*a03*a12;
	M[1][3] =      a00*a12*a23 - a00*a13*a22 - a10*a02*a23 + a10*a03*a22 + a20*a02*a13 - a20*a03*a12;
	M[2][0] =      a10*a21*a33 - a10*a23*a31 - a20*a11*a33 + a20*a13*a31 + a30*a11*a23 - a30*a13*a21;
	M[2][1] =    - a00*a21*a33 + a00*a23*a31 + a20*a01*a33 - a20*a03*a31 - a30*a01*a23 + a30*a03*a21;
	M[2][2] =      a00*a11*a33 - a00*a13*a31 - a10*a01*a33 + a10*a03*a31 + a30*a01*a13 - a30*a03*a11;
	M[2][3] =    - a00*a11*a23 + a00*a13*a21 + a10*a01*a23 - a10*a03*a21 - a20*a01*a13 + a20*a03*a11;
	M[3][0] =    - a10*a21*a32 + a10*a22*a31 + a20*a11*a32 - a20*a12*a31 - a30*a11*a22 + a30*a12*a21;
	M[3][1] =      a00*a21*a32 - a00*a22*a31 - a20*a01*a32 + a20*a02*a31 + a30*a01*a22 - a30*a02*a21;
	M[3][2] =    - a00*a11*a32 + a00*a12*a31 + a10*a01*a32 - a10*a02*a31 - a30*a01*a12 + a30*a02*a11;
	M[3][3] =      a00*a11*a22 - a00*a12*a21 - a10*a01*a22 + a10*a02*a21 + a20*a01*a12 - a20*a02*a11;
	
	double D = a00*M[0][0] + a10*M[0][1] + a20*M[0][2] + a30*M[0][3];
      
	if (D)
	{
		M[0][0] /= D; M[0][1] /= D; M[0][2] /= D; M[0][3] /= D;
		M[1][0] /= D; M[1][1] /= D; M[1][2] /= D; M[1][3] /= D;
		M[2][0] /= D; M[2][1] /= D; M[2][2] /= D; M[2][3] /= D;
		M[3][0] /= D; M[3][1] /= D; M[3][2] /= D; M[3][3] /= D;
	}
	else
		cout << endl << "Singular 4x4 matrix inversion!";
}

/** unrolled 4x4 matrix multiply*/
void Matrix4x4MultiplyBy4x4 (const double src1[4][4], const double src2[4][4], double dest[4][4])
{
	dest[0][0] = src1[0][0] * src2[0][0] + src1[0][1] * src2[1][0] + src1[0][2] * src2[2][0] + src1[0][3] * src2[3][0]; 
	dest[0][1] = src1[0][0] * src2[0][1] + src1[0][1] * src2[1][1] + src1[0][2] * src2[2][1] + src1[0][3] * src2[3][1]; 
	dest[0][2] = src1[0][0] * src2[0][2] + src1[0][1] * src2[1][2] + src1[0][2] * src2[2][2] + src1[0][3] * src2[3][2]; 
	dest[0][3] = src1[0][0] * src2[0][3] + src1[0][1] * src2[1][3] + src1[0][2] * src2[2][3] + src1[0][3] * src2[3][3]; 
	dest[1][0] = src1[1][0] * src2[0][0] + src1[1][1] * src2[1][0] + src1[1][2] * src2[2][0] + src1[1][3] * src2[3][0]; 
	dest[1][1] = src1[1][0] * src2[0][1] + src1[1][1] * src2[1][1] + src1[1][2] * src2[2][1] + src1[1][3] * src2[3][1]; 
	dest[1][2] = src1[1][0] * src2[0][2] + src1[1][1] * src2[1][2] + src1[1][2] * src2[2][2] + src1[1][3] * src2[3][2]; 
	dest[1][3] = src1[1][0] * src2[0][3] + src1[1][1] * src2[1][3] + src1[1][2] * src2[2][3] + src1[1][3] * src2[3][3]; 
	dest[2][0] = src1[2][0] * src2[0][0] + src1[2][1] * src2[1][0] + src1[2][2] * src2[2][0] + src1[2][3] * src2[3][0]; 
	dest[2][1] = src1[2][0] * src2[0][1] + src1[2][1] * src2[1][1] + src1[2][2] * src2[2][1] + src1[2][3] * src2[3][1]; 
	dest[2][2] = src1[2][0] * src2[0][2] + src1[2][1] * src2[1][2] + src1[2][2] * src2[2][2] + src1[2][3] * src2[3][2]; 
	dest[2][3] = src1[2][0] * src2[0][3] + src1[2][1] * src2[1][3] + src1[2][2] * src2[2][3] + src1[2][3] * src2[3][3]; 
	dest[3][0] = src1[3][0] * src2[0][0] + src1[3][1] * src2[1][0] + src1[3][2] * src2[2][0] + src1[3][3] * src2[3][0]; 
	dest[3][1] = src1[3][0] * src2[0][1] + src1[3][1] * src2[1][1] + src1[3][2] * src2[2][1] + src1[3][3] * src2[3][1]; 
	dest[3][2] = src1[3][0] * src2[0][2] + src1[3][1] * src2[1][2] + src1[3][2] * src2[2][2] + src1[3][3] * src2[3][2]; 
	dest[3][3] = src1[3][0] * src2[0][3] + src1[3][1] * src2[1][3] + src1[3][2] * src2[2][3] + src1[3][3] * src2[3][3]; 
}


/** Function useful to analytically compute torus vs ray intersections */
inline double realCubicSolution(double b, double c, double d)
{
	double p = -(1.0/3.0)*b*b + c;
	double q = (2./27.)*b*b*b - (1./3.)*b*c + d;

	double D3p = (1./27.)*p*p*p;
	double D3 = D3p + 0.25*q*q;

	if (fabs(c) < b*b*1.E-6 && fabs(d) < fabs(b)*1.E-6)
	{
		return c/b - b - d/(b*b);
	}
	else if (abs(d/c) < 4.E-4)
	{
		double ratio = d/c;

		if (D3 <= 0.)
		{
			return 0.5*(ratio - b + sqrt((b-ratio)*(b-ratio) - 4.*(c - ratio*(b-ratio))));
		}
		else
		{
			double eps = -ratio*ratio*(b-ratio) / (3.*ratio*ratio - 2.*b*ratio + c);
			return eps - ratio;
		}
	}
	else
	{
		if (D3 >= 0.)
		{
			double radD3 = sqrt(D3);
			double A = cbrt(radD3 - 0.5*q);
			double psu3A = cbrt(radD3 + 0.5*q);
			return A - psu3A - (1./3.)*b;
		}
		else
		{
			double Aarg = acos(-0.5*q / sqrt(-D3p));
			return 2.*sqrt(-(1./3.)*p) * cos((1./3.)*Aarg) - (1./3.)*b;
		}
	}
}


/** Function useful to analytically compute torus vs ray intersections */
/*
inline double realCubicSolution(double b, double c, double d)
{
	double reduced_delta = c*c - 4.*b*d;

	if (fabs(c) < b*b*1.E-6 && fabs(d) < fabs(b)*1.E-6)
	{
		return c/b - b - d/(b*b);
	}
	else if (reduced_delta >= 0. && fabs(d/c) < 1.E-4)
	{
		return -(d+d) / (c + sqrt(reduced_delta));

		// double r = sqrt(b*b - 4.*c);
		// double x1 = 0.5*(r - b);
		// return x1 - d/(x1*r);
	}
	else
	{
		double p = (-1./3.)*b*b + c;
		double q = (2./27.)*b*b*b - (1./3.)*b*c + d;

		double D3p = (1./27.)*p*p*p;
		double D3 = D3p + 0.25*q*q;

		if (D3 >= 0.)
		{
			double radD3 = sqrt(D3);
			double A = cbrt(radD3 - 0.5*q);
			double psu3A = cbrt(radD3 + 0.5*q);
			return A - psu3A - (1./3.)*b;
		}
		else
		{
			double Aarg = acos(-0.5*q / sqrt(-D3p));
			return 2.*sqrt(-(1./3.)*p) * cos((1./3.)*Aarg) - (1./3.)*b;
		}
	}
}
*/

/** Function useful to analytically compute torus vs ray intersections */
/*
inline double realCubicSolution(double b, double c, double d)
{
	double reduced_delta = c*c - 4*b*d;

	if (fabs(c) < b*b*1.E-6 && fabs(d) < fabs(b)*1.E-6)
	{
		return c/b - b - d/(b*b);
	}
	else if (reduced_delta >= 0. && fabs(d/c) < 1.E-4)
	{
		return -(d+d) / (c + sqrt(reduced_delta));
	}
	else
	{
		double p = (-1./3.)*b*b + c;
		double q = (2./27.)*b*b*b - (1./3.)*b*c + d;

		double D3 = (1./27.)*p*p*p + 0.25*q*q;

		if (D3 >= 0)
		{
			double radD3 = sqrt(D3);
			double A = cbrt(radD3 - 0.5*q);
			double psu3A = cbrt(radD3 + 0.5*q);
			return A - psu3A - (1./3.)*b;
		}
		else
		{
			complex<double> temp = 1i * std::sqrt(-D3) - 0.5*q;
			complex<double> A = std::pow(temp, 1./3.);
			return 2.*A.real() - (1./3.)*b;
		}
	}
}
*/

/*
inline double realCubicSolution(double b, double c, double d)
{
	double p = (-1./3.)*b*b + c;
	double q = (2./27.)*b*b*b - (1./3.)*b*c + d;

	double D3 = (1./27.)*p*p*p + 0.25*q*q;

	if (D3 >= 0)
	{
		double A = cbrt(sqrt(D3) - 0.5*q);
		return A - (1./3.)*(p/A + b);
	}
	else
	{
		// A = (-q/2+sqrt(D3))^(1/3);
		// complex<double> temp = std::sqrt(std::complex<double>(D3)) - 0.5*q;
		complex<double> temp = 1i * std::sqrt(-D3) - 0.5*q;
		complex<double> A = std::pow(temp, 1./3.);
		return 2.*A.real() - (1./3.)*b;
	}
}
*/


/*
inline double realCubicSolution(double b, double c, double d)
{
	double p = (-1./3.)*b*b + c;
	double q = (2./27.)*b*b*b - (1./3.)*b*c + d;

	double D3 = (1./27.)*p*p*p + 0.25*q*q;

	if (D3 >= 0)
	{
		double radD3 = sqrt(D3);
		double A = cbrt(radD3 - 0.5*q);
		double psu3A = cbrt(radD3 + 0.5*q);
		return A - psu3A - (1./3.)*b;
	}
	else
	{
		complex<double> temp = 1i * std::sqrt(-D3) - 0.5*q;
		complex<double> A = std::pow(temp, 1./3.);
		return 2.*A.real() - (1./3.)*b;
	}
}
*/


/** Function useful to analytically compute torus vs ray intersections */
inline void quarticEqSolutions(double roots[4], double b, double c, double d, double e, int *num_sol)
{
	// The following 3 lines are rewritten below ...
	// double p = c - (3.0/8.0)*(b*b);
	// double q = 0.125*(b*b*b) - 0.5*b*c + d;
	// double r = e - (3./(4.*4.*4.*4.))*((b*b*b)*b) - 0.25*b*d + 0.0625*c*(b*b);

	double b_2 = 0.5*b;
	double b_4 = 0.25*b;
	double p = c - 1.5*(b_2*b_2);
	// double q = accurateSum3(0.125*(b*b*b), -0.5*b*c, d);
	double q = d + b_2 * (b_2*b_2 - c);
	// double r = accurateSum4(e, -(3./(4.*4.*4.*4.))*((b*b*b)*b), -0.25*b*d, 0.0625*c*(b*b));
	double r = e + b_4 * (b_4 * (c - 3.0*b_4*b_4) - d);

	double f = p*p - 4.*r;

	// double delta = 4.*f * accurateSum3(4.*r*(p*p), -p*(q*q), -16.*(r*r)) + (q*q) * (128.*p*r - 27.*(q*q));
	double delta = 4.*f * (4.*r*(p*p) - p*(q*q) - 16.*(r*r)) + (q*q) * (128.*p*r - 27.*(q*q));

	bool change_variant = false;

	if (e != 0. && (fabs(p) > 3.E+3 || fabs(r) > 3.E+3 || fabs(q) > 3.E+3 || fabs(delta) > 1.E+14))
	{
		double dsu2e = 0.5*d / e;

		change_variant = true;

		p = -1.5* dsu2e*dsu2e + c/e;
		// q = accurateSum3(dsu2e*dsu2e*dsu2e, -dsu2e*c/e, b/e);
		q = dsu2e * (dsu2e*dsu2e - c/e) + b/e;
		// r = accurateSum4(-3./16.*dsu2e*dsu2e*dsu2e*dsu2e, 0.25*dsu2e*dsu2e*(c/e), -0.5*dsu2e*b/e, 1./e);
		r = dsu2e * (-3./16.*(dsu2e*dsu2e*dsu2e) + 0.25*dsu2e*(c/e) - 0.5*b/e) + 1./e;
		f = p*p - 4.*r;
		// delta = 4.*f * accurateSum3(4.*r*(p*p), -p*(q*q), -16.*(r*r)) + (q*q) * (128.*p*r - 27.*(q*q));
		delta = 4.*f * (4.*r*(p*p) - p*(q*q) - 16.*(r*r)) + (q*q) * (128.*p*r - 27.*(q*q));
	}

	if (delta < 0. || (delta == 0. && (p >= 0. || f < 0.)))
		*num_sol = 2;
	else if ((p < 0. && f > 0.) || (delta == 0. && p == 0. && r == 0.))
		*num_sol = 4;
	else {
		*num_sol = 0;
		return;
	}


	double m = realCubicSolution(p, 0.25*(p*p) - r, (-0.125)*(q*q));

	if (m <= 0)
	{
		*num_sol = 0;
		return;
	}

	double coeff;

	if (*num_sol >= 2)
	{
		double sol_delta_min, sol_delta_max;
		double p1, p2, p3;

		if (m < 1.E-14)
		{
			p1 = p+p;
			p2 = (q*q)/f;
			p3 = 2.*sqrt(f);

			// sol_delta_max = accurateSum3(p3, -p1, -p2);
			sol_delta_max = p3 - p1 - p2;
			coeff = 2.*fabs(q)/f;
		}
		else
		{
			p1 = -2.*(p+m);
			p2 = sqrt(2.*(q*q)/m);

			sol_delta_max = p1 + p2;
			coeff = sqrt(2.*m);
		}
		roots[0] = 0.5 * (((q>=0.) ? -coeff : coeff) - sqrt(sol_delta_max));
		roots[1] = 0.5 * (((q>=0.) ? -coeff : coeff) + sqrt(sol_delta_max));

		if (*num_sol == 4)
		{
			if (m < 1.E-14)
			{
				// sol_delta_min = -accurateSum3(p1, p2, p3);
				sol_delta_min = -p1 - p2 - p3;
			}
			else
			{
				sol_delta_min = p1 - p2;
			}
			roots[2] = 0.5 * (((q>=0.) ? coeff : -coeff) - sqrt(sol_delta_min));
			roots[3] = 0.5 * (((q>=0.) ? coeff : -coeff) + sqrt(sol_delta_min));
		}
	}

	if (change_variant)
	{
		for (int i = 0; i < *num_sol; i++)
		{
			roots[i] = 1. / (roots[i] - 0.25*d/e);
		}
	}
	else
	{
		for (int i = 0; i < *num_sol; i++)
		{
			roots[i] -= 0.25*b;
		}
	}
}


/** Function useful to analytically compute torus vs ray intersections */
/*
inline void quarticEqSolutions(double roots[4], double b, double c, double d, double e, int *num_sol)
{
	double p = c - (3.0/8.0)*(b*b);
	double q = 0.125*(b*b*b) - 0.5*b*c + d;
	double r = e - (3./(4.*4.*4.*4.))*((b*b*b)*b) - 0.25*b*d + 0.0625*c*(b*b);

	double delta = 16.*(p*p*p*p)*r - 4.*(p*p*p)*(q*q) - 128.*(p*p)*(r*r) + 144.*p*(q*q)*r - 27.*(q*q)*(q*q) + 256.*(r*r*r);

	double f = p*p - 4.*r;

	if (delta < 0. || (delta == 0. && (p >= 0. || f < 0.)))
		*num_sol = 2;
	else if ((p < 0. && f > 0.) || (delta == 0. && p == 0. && r == 0.))
		*num_sol = 4;
	else {
		*num_sol = 0;
		return;
	}

	double m = realCubicSolution(p, 0.25*(p*p) - r, (-0.125)*(q*q));

	if (m <= 0)
	{
		*num_sol = 0;
		return;
	}

	double coeff;

	if (*num_sol >= 2)
	{
		double sol_delta_min, sol_delta_max;
		double p1, p2, p3;

		if (m < 1.E-14)
		{
			p1 = p+p;
			p2 = (q*q)/f;
			p3 = 2.*sqrt(f);

			sol_delta_max = p3 - p1 - p2;
			coeff = 2.*fabs(q)/f;
		}
		else
		{
			p1 = -2.*(p+m);
			p2 = sqrt(2.*(q*q)/m);

			sol_delta_max = p1 + p2;
			coeff = sqrt(2.*m);
		}
		roots[0] = 0.5 * (((q>=0.) ? -coeff : coeff) - sqrt(sol_delta_max)) - 0.25*b;
		roots[1] = 0.5 * (((q>=0.) ? -coeff : coeff) + sqrt(sol_delta_max)) - 0.25*b;

		if (*num_sol == 4)
		{
			if (m < 1.E-14)
			{
				sol_delta_min = -p1 - p2 - p3;
			}
			else
			{
				sol_delta_min = p1 - p2;
			}
			roots[2] = 0.5 * (((q>=0.) ? -coeff : coeff) - sqrt(sol_delta_min)) - 0.25*b;
			roots[3] = 0.5 * (((q>=0.) ? -coeff : coeff) + sqrt(sol_delta_min)) - 0.25*b;
		}
	}
}
*/

/** Function useful to analytically compute torus vs ray intersections */
/*
inline void quarticEqSolutions(double roots[4], double b, double c, double d, double e, int *num_sol)
{
	double p = c - (3.0/8.0)*(b*b);
	double q = 0.125*(b*b*b) - 0.5*b*c + d;
	double r = e - (3./(4.*4.*4.*4.))*((b*b*b)*b) - 0.25*b*d + 0.0625*c*(b*b);

	double delta = 16.*(p*p*p*p)*r - 4.*(p*p*p)*(q*q) - 128.*(p*p)*(r*r) + 144.*p*(q*q)*r - 27.*(q*q)*(q*q) + 256.*(r*r*r);

	if (delta < 0.)
		*num_sol = 2;
	else if (p < 0. && r < 0.25*(p*p))
		*num_sol = 4;
	else {
		*num_sol = 0;
		return;
	}
	// double tol = 1.E-7;
	// if (fabs(delta) < tol)
	//     cout << "possible coincident solutions" << endl;

	double m = realCubicSolution(p, 0.25*(p*p) - r, (-0.125)*(q*q));

	if (m <= 0)
	{
		*num_sol = 0;
		return;
	}

	double coeff = sqrt(m+m);
	double sol_delta;

	if (*num_sol == 2)
	{
		sol_delta = sqrt(2.*q*q/m) - (p+p) - (m+m);

		if (sol_delta < 0.)
			*num_sol = 0;
		else
		{
			roots[0] = 0.5 * (((q>=0.) ? -coeff : coeff) - sqrt(sol_delta)) - 0.25*b;
			roots[1] = 0.5 * (((q>=0.) ? -coeff : coeff) + sqrt(sol_delta)) - 0.25*b;
		}
	}
	else
	{
		sol_delta = 2.0 * (q/coeff - p - m);
		roots[0] = 0.5 * (-sqrt(sol_delta) - coeff) - 0.25*b;
		roots[1] = 0.5 * ( sqrt(sol_delta) - coeff) - 0.25*b;
		sol_delta = -2.0 * (p + m + q/coeff);
		roots[2] = 0.5 * (coeff - sqrt(sol_delta)) - 0.25*b;
		roots[3] = 0.5 * (coeff + sqrt(sol_delta)) - 0.25*b;
	}
}
*/

/** Function useful to analytically compute torus vs ray intersections */
/*
inline void quarticEqSolutions(double roots[4], double b, double c, double d, double e, int *num_sol)
{
	double p = c - (3.0/8.0)*(b*b);
	double q = 0.125*(b*b*b) - 0.5*b*c + d;
	double r = e - (3./(4.*4.*4.*4.))*((b*b*b)*b) - 0.25*b*d + 0.0625*c*(b*b);

	double delta = 16.*(p*p*p*p)*r - 4.0*(p*p*p)*(q*q) - 128.*(p*p)*(r*r) + 144.*p*(q*q)*r - 27.*(q*q)*(q*q) + 256.*(r*r*r);

	if (delta < 0.)
		*num_sol = 2;
	else if (p < 0. && r < 0.25*(p*p))
		*num_sol = 4;
	else {
		*num_sol = 0;
		return;
	}
	// double tol = 1.E-7;
	// if (fabs(delta) < tol)
	//     cout << "possible coincident solutions" << endl;

	double m = realCubicSolution(p, 0.25*(p*p) - r, (-0.125)*(q*q));

	if (m <= 0)
	{
		*num_sol = 0;
		return;
	}

	double coeff = sqrt(m+m);
	double sol_delta;

	if (*num_sol == 2)
	{
		sol_delta = sqrt(2.*q*q/m) - (p+p) - (m+m);

		if (sol_delta < 0.)
			*num_sol = 0;
		else
		{
			roots[0] = 0.5 * (((q>=0.) ? -coeff : coeff) - sqrt(sol_delta)) - 0.25*b;
			roots[1] = 0.5 * (((q>=0.) ? -coeff : coeff) + sqrt(sol_delta)) - 0.25*b;
		}
	}
	else
	{
		sol_delta = 2.0 * (q/coeff - p - m);
		roots[0] = 0.5 * (-sqrt(sol_delta) - coeff) - 0.25*b;
		roots[1] = 0.5 * ( sqrt(sol_delta) - coeff) - 0.25*b;
		sol_delta = -2.0 * (p + m + q/coeff);
		roots[2] = 0.5 * (coeff - sqrt(sol_delta)) - 0.25*b;
		roots[3] = 0.5 * (coeff + sqrt(sol_delta)) - 0.25*b;
	}
}
*/


#if !defined(USE_NEW_RAY_VS_SPHERE_ALGORITHM)

/** ray/sphere intersection. ray is o+t*dir */
bool raySphere (const double *const orig, const double *const dir, const double *const sphere_center,
				const double sphere_radius, double *const t1, double *const t2)
{
	// perform sphere intersection test
	double temp[3], A, B, C;

	A = DOT(dir,dir);
	SUB(temp,orig,sphere_center)
	B = DOT(temp,dir);
	C = DOT(temp,temp) - sphere_radius*sphere_radius;
	double det = B*B - A*C;

	// no intersection
	if (det < 0.)
		return false;

	det = sqrt(det);

	*t1 = (-B-det) / A;
	*t2 = (-B+det) / A;

	if (*t2 < *t1)
	{
		double tt = *t1;
		*t1 = *t2;
		*t2 = tt;
	}
	return true;
}

#else // USE_NEW_RAY_VS_SPHERE_ALGORITHM

bool raySphere (const double *const orig, const double *const dir, const double *const sphere_center,
				const double sphere_radius, double *const t1, double *const t2)
{
	// perform sphere intersection test.
	// This version is more accurate than the previous one. It is an optimised version
	// of the algorithm provided in:
	// Hearn, D. D., and Baker, M. P. Computer Graphics with OpenGL, third ed. Pearson, 2004
	// (reviewed in https://link.springer.com/content/pdf/10.1007/978-1-4842-4427-2_7.pdf)
	double f[3], f_prime[3], A, B, B_prime;

	A = DOT(dir,dir);
	SUB(f,orig,sphere_center)
	B = DOT(f,dir);
	B_prime = B / A;

	f_prime[0] = f[0] - B_prime*dir[0];
	f_prime[1] = f[1] - B_prime*dir[1];
	f_prime[2] = f[2] - B_prime*dir[2];

	double det = sphere_radius*sphere_radius - DOT(f_prime,f_prime);

	// no intersection
	if (det < 0.)
		return false;

	det = sqrt(A * det);

	*t1 = (-B-det) / A;
	*t2 = (-B+det) / A;

	if (*t2 < *t1)
	{
		double tt = *t1;
		*t1 = *t2;
		*t2 = tt;
	}
	return true;
}

#endif // USE_NEW_RAY_VS_SPHERE_ALGORITHM


bool rayTorus (int analytical, double invrot[3][3], double torus_center[3], double sphere_center[3],
			   double probe_radius, double major_radius, double radius,
			   int panel, double orig[3], double dir[3], double t[4], int &numInt)
{
	// Attention: t[0] and t[1] have to be provided following the intersection of the ray with
	// the clipping sphere and are rewritten below if there are more than two intersections
	double roots[4], pp[3], temp[3];
	double orig1[3];
	// double dir_norm;

	// start the ray from the clipping sphere intersection
	// ADD_MUL(orig1,orig,dir,t[0])
	// ADD_MUL(orig2,orig,dir,t[1])
	// SUB(tdir,orig2,orig1)

	// some optimised code ...
	int varying_coord = (panel == 0) ? 0 : ((panel == 1) ? 2 : 1);
	double ray_dir = dir[varying_coord];

	ASSIGN(orig1,orig);
	orig1[varying_coord] += t[0] * ray_dir;
	// dir_norm = (t[1] - t[0]) * ray_dir;

	// build torus roots equation
	SUB(temp,orig1,torus_center)

	for (int i=0;i<3;i++)
	{
		// roto-translate origin point
		pp[i] = invrot[i][0]*temp[0] + invrot[i][1]*temp[1] + invrot[i][2]*temp[2];
		// rotate ray direction
		// ddir[i] = invrot[i][varying_coord]*dir_norm;
	}
	double r2 = probe_radius * probe_radius;
	double R2 = major_radius * major_radius;

	double pd,p2,coeff;

	// pd = DOT(pp,ddir);
	pd = pp[0]*invrot[0][varying_coord] + pp[1]*invrot[1][varying_coord] + pp[2]*invrot[2][varying_coord];
	p2 = DOT(pp,pp);
	coeff = p2-r2-R2;

	double B = 4*pd;
	double C = 2*coeff + 4*pd*pd + 4*R2*invrot[2][varying_coord]*invrot[2][varying_coord];
	double D = 4*pd*coeff + 8*R2*pp[2]*invrot[2][varying_coord];
	double E = coeff*coeff - 4*R2*(r2-pp[2]*pp[2]);

	int numroots;

	#if !defined(CHECK_ACCURACY_DIFF)
	if (analytical)
	{
		quarticEqSolutions(roots, B, C, D, E, &numroots);
	}
	else
	{
		double polyy[5];

		numroots = 4;

		polyy[0] = E;
		polyy[1] = D;
		polyy[2] = C;
		polyy[3] = B;
		polyy[4] = 1.0;

		int order = 4;

		getRealRootsSturm(polyy,order,roots,numroots);
		// getRealRootsCompanion(polyy,order,roots,numroots);
	}
	#else // CHECK_ACCURACY_DIFF

	double temp_roots[4];
	numroots = 0;
	int temp_numroots;

	quarticEqSolutions(temp_roots, B, C, D, E, &temp_numroots);

	for (int i=0; i<temp_numroots; i++)
	{
		if (temp_roots[i] < 0.) continue;

		double point[3], dd;

		ASSIGN(point,orig1)
		point[varying_coord] += temp_roots[i]; // roots[i]*dir_norm
		DIST2(dd,sphere_center,point)

		if (dd >= radius*radius) continue;

		roots[numroots++] = temp_roots[i];
	}

	double polyy[5];

	polyy[0] = E;
	polyy[1] = D;
	polyy[2] = C;
	polyy[3] = B;
	polyy[4] = 1.0;

	double temp_roots2[4];
	int temp_numroots2 = 4;

	getRealRootsSturm(polyy,4,temp_roots2,temp_numroots2);

	if (numroots == 0 && temp_numroots2 == 0)
		return false;

	if (numroots != temp_numroots2)
	{
		double roots2[4];
		int numroots2 = 0;

		for (int i=0; i<temp_numroots2; i++)
		{
			if (temp_roots2[i] < 0.) continue;
			// if (temp_roots2[i] < 0.) { numroots2 = 0; break; }

			double point[3], dd;

			ASSIGN(point,orig1)
			point[varying_coord] += temp_roots2[i]; // roots[i]*dir_norm
			DIST2(dd,sphere_center,point)

			if (dd >= radius*radius) continue;

			roots2[numroots2++] = temp_roots2[i];
		}
		if (numroots2 != numroots)
		{
			double ddir[3], dir_norm;

			dir_norm = (t[1] - t[0]) * ray_dir;

			for (int i=0;i<3;i++)
				ddir[i] = invrot[i][varying_coord]*dir_norm;

			printf ("\n");
			printf ("%.1f %.3f %.3f ", probe_radius, major_radius, radius);
			printf ("%.14e %.14e %.14e ", pp[0], pp[1], pp[2]);
			printf ("%.14e %.14e %.14e ", ddir[0], ddir[1], ddir[2]);
			printf ("%.14f %.14f %.14f %.14f %i ", B, C, D, E, numroots);

			for (int i=0; i<numroots; i++)
				printf ("%.14f ", roots[i]);

			printf ("%i ", numroots2);

			for (int i=0; i<numroots2; i++)
				printf ("%.14f ", roots2[i]);

			printf ("\n");
		}
	}
	#endif // CHECK_ACCURACY_DIFF

	if (numroots == 0)
		return false;

	double t0 = t[0];
	numInt = 0;

	for (int i=0; i<numroots; i++)
	{
		if (roots[i] < 0)
			continue;

		double point[3], dd;

		// roots[i] /= dir_norm;
		// ADD_MUL(point,orig1,tdir,roots[i])
		ASSIGN(point,orig1)
		point[varying_coord] += roots[i]; // roots[i]*dir_norm
		DIST2(dd,sphere_center,point)
		// acceptable
		if (dd < radius*radius)
		{
			// t[ numInt++ ] = (point[v_coord]-orig[v_coord]) / ray_dir;
			//               = (orig1[v_coord] (=pa[v_coord]+t[0]*ray_dir) + roots[i] - pa[v_coord]) / ray_dir =
			//               = t0 + roots[i]/ray_dir
			t[ numInt++ ] = t0 + roots[i]/ray_dir;
		}
	}

	if (numInt == 0)
		return false;

	return true;
}


/** get the normal to a sphere*/
void getNormalToSphere(const double *const y, const double *const center,
					   const double radius, double *const normal)
{
	SUB(normal,y,center);
	normal[0] /= radius;
	normal[1] /= radius;
	normal[2] /= radius;
}

/** project y to sphere and returns the normal*/
void projectToSphere(const double *const y, const double *const center,
					 const double radius, double *const proj, double &dist)
{
	DIST(dist,y,center)
	dist = radius / dist;
	SUB(proj,y,center)
	ADD_MUL(proj,center,proj,dist)
	DIST(dist,proj,y)
}

void planeByplane(const double w1[4], const double w2[4], double Q[3][3], double a[3], double &c)
{
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			Q[i][j] = w1[i]*w2[j];

	for (int i=0; i<3; i++)
		a[i] = w1[3]*w2[i] + w2[3]*w1[i];

	c = w1[3]*w2[3];
}

void getMemSpace (double &current_mem_in_MB, double &peak_mem_in_MB)
{
	#if !defined(AVOID_MEM_CHECKS)
	/*
	// old operating system-dependent method to get mem. consumption, see below new method
	system("pgrep -x NanoShaper >> NS_job_id.txt");

	FILE *fp = fopen("NS_job_id.txt", "r");

	char job_id_string[256];
	fscanf (fp, "%s", job_id_string);

	char s1[1024] = "cat /proc/";

	strcat (s1, job_id_string);
	strcat (s1, "/status | grep -e VmHWM -e VmRSS >> memory.txt");
	system (s1);

	fclose(fp);

	system ("rm NS_job_id.txt");

	char memory_string[1024];

	fp = fopen("memory.txt", "r");

	int current_memory_in_kB, peak_memory_in_kB;

	fscanf (fp, "%s", memory_string);
	fscanf (fp, "%i", &peak_memory_in_kB);
	fscanf (fp, "%s\n", memory_string);
	fscanf (fp, "%s", memory_string);
	fscanf (fp, "%i", &current_memory_in_kB);

	fclose(fp);

	system ("rm memory.txt");

	current_mem_in_MB = (double)current_memory_in_kB / 1024.;
	peak_mem_in_MB    = (double)peak_memory_in_kB    / 1024.;
	*/

	struct rusage current_mem_usage;

	static double peak_mem = 0.;


	getrusage(RUSAGE_SELF, &current_mem_usage);

	current_mem_in_MB = current_mem_usage.ru_maxrss / 1024.;

	peak_mem = fmax(peak_mem, current_mem_in_MB);

	peak_mem_in_MB = peak_mem;
	#endif
}


void initTBB (int num_threads)
{
	// auto mp = tbb::global_control::max_allowed_parallelism;
	// tbb::global_control gc(mp, num_threads);
	static tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);
}
