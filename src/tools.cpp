
#include <tools.h>
#include <logging.h>

/**@brief ascending on first double of pair<double,double*> comparator*/
bool compKeepIndex(pair<double,double*> a, pair<double,double*> b)
{
	return a.first < b.first;
}

///////////////////////////////////////////////////////////////////

bool index_double_comparator( const indexed_double& l, const indexed_double& r)
{ 
	return l.first < r.first; 
}

///////////////////// MISCELLANEAOUS ////////////////////////////////

/** @brief test if point proj is in cube whose center is point and side is side variable.
 If a toll is provided, the test is performed up to a tollerance value.
 By default the tollerance is 0*/
bool testInCube(const double proj0,const double proj1,const double proj2,const double point0,const double point1,const double point2,const double side,const double toll)
{
	double ss = side/2;
	double minx = point0-ss;
	double miny = point1-ss;
	double minz = point2-ss;
	double maxx = point0+ss;
	double maxy = point1+ss;
	double maxz = point2+ss;

	if (proj0>minx-toll && proj0<maxx+toll 
        && proj1>miny-toll && proj1<maxy+toll
        && proj2>minz-toll && proj2<maxz+toll
        )
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
    for (unsigned int i=0;i<strlen(str.c_str());i++)
		if (str[i] >= 0x41 && str[i] <= 0x5A)
			str[i] = str[i] + 0x20;
    return str;
}

void cleanLine()
{
	logging::log<logging::level::info>("--------------------------------------------------------------------------");
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
	if (test>=0)
		SEED = test;
	else
		SEED = test +mm;
	reslt = SEED/(double)mm;
	return( reslt );
}


/** get the real roots by computing the companion matrix and then extracting the eigenvalues. A root is real
if its imaginary part in absolute value is less than a given threshold. Usually this threshold is rather conservative
such that possibly unprecise real roots are not lost*/
void getRealRootsCompanion(double *const poly,const int degree,double *const roots,int& numroots)
{	
	int zero_roots = 0,ind=0,ind2=0;
	bool first = true;
	int numelem = degree;	
	int truedegree = degree;
	//int* notNull = allocateVector<int>(degree+1);
	int notNull[MAX_POLY_DEGREE+1];

	// find non zeros elements and record indexes
	for (int i=0;i<degree+1;i++)
	{
		if (poly[i]!=0)
		{
			notNull[ind2]=i;	
			ind2++;
		}
	}

	// remove all zeros and record the zeros roots; in place rewrite the true polynomial
	ind = 0;
	for (int i=notNull[0];i<=notNull[ind2-1];i++)
	{
		poly[ind]=poly[i];
		ind++;
	}
	numelem = ind;
	truedegree = notNull[ind2-1];
	// number of zero roots	
    zero_roots = notNull[0];	
	
	//double* localpoly=allocateVector<double>(numelem-1);
	double localpoly[MAX_POLY_DEGREE+1];

	// get the real degree of the polynomial by stripping almost zero leading coefficients
	// almost-zero leading coefficient are those which produce 1/c = INF; this enhances
	// numerical stability in case bad coefficients are supplied
	bool oneinf = true;
	while(oneinf)
	{
		oneinf = false;
		// try to normalize coefficients
		for (int i=0;i<numelem-1;i++)
		{
			localpoly[i]=poly[i]/poly[numelem-1];
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
	TNT::Array2D<double> M(numelem-1,numelem-1);
	
	// build the companion matrix of the filtered polynomial
	for (int i=0;i<numelem-1;i++)
		for (int j=0;j<numelem-1;j++)			
			M[i][j]=0;

	for (int i=0;i<numelem-1;i++)
		M[i][numelem-2]=-localpoly[i];

	for (int i=1;i<numelem-1;i++)
		for (int j=0;j<numelem-2;j++)
			if ((i-1)==j)
				M[i][j]=1;
	
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

	for (i=0;i<final;i++)
	{
		// identify real roots by using a prescribed tollerance
		if (fabs(img[i])<1e-1)
		{
			roots[numroots]=real[i];
			numroots++;			
		}
		else
			numelem--;
	}	

	// add the zero roots
	for (int j=numroots;j<zero_roots+numroots;j++)
	{
		roots[j]=0;
		numroots++;
	}
}


/** get real roots by using Sturm method. Directly search real roots. Much faster, often less accurate
than companion matrix*/
void getRealRootsSturm(const double *const polyy,const int degree,double *const roots,int& numrootss)
{
	int nchanges, np, atmin, atmax,nroots,i;	
	double min,max;

	poly	sseq[MAX_POLY_DEGREE+1];

	for (int i=0;i<degree+1;i++)
		sseq[0].coef[i]=polyy[i];

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
		logging::log<logging::level::info>("solve: unable to bracket all negative roots");
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
		logging::log<logging::level::info>("solve: unable to bracket all positive roots");
		atmax = nchanges;
	}

	nroots = atmin - atmax;

	// perform the bisection.
	sbisect(np, sseq, min, max, atmin, atmax, roots);

	numrootss = nroots;
}


/** plane by 3 points routine*/
void plane3points(const double p1[3],const double p2[3],const double p3[3],double w[4],const bool normalize)
{
	w[0] = p1[1]*(p2[2] - p3[2]) + p2[1]*(p3[2] - p1[2]) + p3[1]*(p1[2] - p2[2]);
	w[1] = p1[2]*(p2[0] - p3[0]) + p2[2]*(p3[0] - p1[0]) + p3[2]*(p1[0] - p2[0]);
	w[2] = p1[0]*(p2[1] - p3[1]) + p2[0]*(p3[1] - p1[1]) + p3[0]*(p1[1] - p2[1]);
	w[3] = -DOT(w,p1);

	if (normalize)
	{
		double norm = sqrt(DOT(w,w));
		w[0]/=norm;
		w[1]/=norm;
		w[2]/=norm;		
		w[3]/=norm;
	}
	return;
}

/** point to plane projection*/
void point2plane(const double p[3],double w[4],double* const dist, double proj[3])
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
	double a00 = M[0][0]; double a01 = M[0][1];double a02 = M[0][2];double a03 = M[0][3];
    double a10 = M[1][0];double a11 = M[1][1];double a12 = M[1][2];double a13 = M[1][3];
	double a20 = M[2][0];double a21 = M[2][1];double a22 = M[2][2];double a23 = M[2][3];
	double a30 = M[3][0];double a31 = M[3][1];double a32 = M[3][2];double a33 = M[3][3];

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
    
	double D = a00*M[0][0] + a10*M[0][1] +  a20*M[0][2] + a30*M[0][3];
      
	if(D)
    {
        M[0][0] /=D; M[0][1] /=D; M[0][2] /=D; M[0][3] /=D;
        M[1][0] /=D; M[1][1] /=D; M[1][2] /=D; M[1][3] /=D;
        M[2][0] /=D; M[2][1] /=D; M[2][2] /=D; M[2][3] /=D;
        M[3][0] /=D; M[3][1] /=D; M[3][2] /=D; M[3][3] /=D;
    }
	else
    logging::log<logging::level::info>("Singular 4x4 matrix inversion!");
}

/** unrolled 4x4 matrix multiply*/
void Matrix4x4MultiplyBy4x4 (const double src1[4][4],const double src2[4][4], double dest[4][4])
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

/** ray/sphere intersection. ray is o+t*dir */	
bool raySphere(const double *const orig,const double *const dir,const double *const sphere_center,const double sphere_radius,double *const t1,double *const t2)
{
	// perform sphere intersection test
	double A,B,C,temp[3],/*ndir[3],*/tt;
	A = DOT(dir,dir);
	SUB(temp,orig,sphere_center)
	B = DOT(temp,dir);
	B*=2;
	C= DOT(temp,temp)-sphere_radius*sphere_radius;
	double det = B*B-4*A*C;

	// no intersection
	if (det<0)
		return false;

	det = sqrt(det);

	(*t1)= (-B-det)/(2*A);
	(*t2)= (-B+det)/(2*A);

	if ((*t2)<(*t1))
	{
		tt = (*t1);
		(*t1)=(*t2);
		(*t2)=tt;
	}

	return true;
}

/** get the normal to a sphere*/
void getNormalToSphere(const double *const y,const double *const center,const double radius,double *const normal)
{
	SUB(normal,y,center);
	normal[0]/=radius;
	normal[1]/=radius;
	normal[2]/=radius;
}		

/** project y to sphere and returns the normal*/
void projectToSphere(const double* const y,const double *const center,const double radius,double *const proj,double& dist)
{
	DIST(dist,y,center)
	dist = radius/dist;
	SUB(proj,y,center)
	ADD_MUL(proj,center,proj,dist)
	DIST(dist,proj,y)
}
