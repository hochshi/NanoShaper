
#include "ExternalSurface.h"

void ExternalSurface::init()
{
	surfType = GENERIC_SURFACE;
}

void ExternalSurface::init(ConfigFile* cf)
{
	
}

ExternalSurface::~ExternalSurface()
{
	clear();
}

void ExternalSurface::clear()
{	
}

ExternalSurface::ExternalSurface(DelPhiShared* ds):Surface()
{
	// set environment
	delphi = ds;	
	init();
}

ExternalSurface::ExternalSurface(ConfigFile* cf,DelPhiShared* ds)
{
	// set environment
	delphi = ds;	
	init();
	init(cf);
}

ExternalSurface::ExternalSurface():Surface()
{
	init();
}

bool ExternalSurface::build()
{
	return true;
}


bool ExternalSurface::load(char* fileName) 
{
	return true;
}


bool ExternalSurface::save(char* fileName)
{
	cout << endl << INFO << "Saving externally loaded surface by saving DelphiShared Object";
	delphi->saveEpsMaps(fileName);
	delphi->saveStatus(fileName);
	delphi->saveBGP(fileName);
	return true;
}

void ExternalSurface::printSummary()
{
	return;
}

void ExternalSurface::getLocalArea(double gridPoint[3],double* area)
{
	cout << endl << WARN << "Cannot compute area in an external surface!";
	(*area)=0;
}
	
double ExternalSurface::getSurfaceArea()
{
	cout << endl << WARN << "Cannot compute area in an external surface!";
	return 0.0;
}


double ExternalSurface::getVolume()
{
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	int ix,iy,iz;
	
	double cubeVol = delphi->side*delphi->side*delphi->side;
	int count = 0;	
	int value[6];
	for (ix=0;ix<NX;ix++)
		for (iy=0;iy<NY;iy++)
			for (iz=0;iz<NZ;iz++)
			{
				/*
				value[0] = delphi->EPSMAP(ix,iy,iz,0,NX,NY,NZ);
				value[1] = delphi->EPSMAP(ix,iy,iz,1,NX,NY,NZ);
				value[2] = delphi->EPSMAP(ix,iy,iz,2,NX,NY,NZ);
				value[3] = delphi->EPSMAP(ix-1,iy,iz,0,NX,NY,NZ);
				value[4] = delphi->EPSMAP(ix,iy-1,iz,1,NX,NY,NZ);
				value[5] = delphi->EPSMAP(ix,iy,iz-1,2,NX,NY,NZ);
				*/
				
				value[0] = read4DVector<int>(delphi->epsmap,ix,iy,iz,0,NX,NY,NZ,3);
				value[1] = read4DVector<int>(delphi->epsmap,ix,iy,iz,1,NX,NY,NZ,3);
				value[2] = read4DVector<int>(delphi->epsmap,ix,iy,iz,2,NX,NY,NZ,3);
				value[3] = read4DVector<int>(delphi->epsmap,ix-1,iy,iz,0,NX,NY,NZ,3);
				value[4] = read4DVector<int>(delphi->epsmap,ix,iy-1,iz,1,NX,NY,NZ,3);
				value[5] = read4DVector<int>(delphi->epsmap,ix,iy,iz-1,2,NX,NY,NZ,3);

				count = 0;
				for (int j=0;j<6;j++)
				{
					if (value[j]!=0)
						count++;
				}
				totalVolume+=cubeVol*(count/6);
			}
	return totalVolume;
}


bool ExternalSurface::getSurf(bool fillCav,double vol)
{
	cout << endl << INFO << "Loading external surface...";

	FILE *fepsx,*fepsy,*fepsz,*fproj,*fstatus;
	fepsx = fopen("epsmapx.txt","r");
	fepsy = fopen("epsmapy.txt","r");
	fepsz = fopen("epsmapz.txt","r");
	fproj = fopen("projections.txt","r");
	fstatus = fopen("status.txt","r");

	if (fepsx==NULL || fepsy==NULL || fepsz==NULL || fproj==NULL || fstatus==NULL)
	{	
		cout << endl << WARN << "Cannot open one of the external surface file";
		return false;
	}

	int tempint;
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	int ix,iy,iz;
	int nbgp;
	
	for (ix=0;ix<NX;ix++)
		for (iy=0;iy<NY;iy++)
			for (iz=0;iz<NZ;iz++)
			{
				int check = fscanf(fepsx,"%d",&tempint);
				if (check<=0)
				{
					cout << endl << ERR << "Error in reading from external surface";
					exit(-1);
				}
				if (tempint)
					inside = tempint;
				
				//delphi->EPSMAP(ix,iy,iz,0,NX,NY,NZ)=tempint;
				write4DVector<int>(delphi->epsmap,tempint,ix,iy,iz,0,NX,NY,NZ,3);
			}
	
	if (ix!=NX)
	{
		cout << endl << WARN << "Error reading epsmapx!";
		return false;
	}
	fclose(fepsx);
	
	for (ix=0;ix<NX;ix++)
		for (iy=0;iy<NY;iy++)
			for (iz=0;iz<NZ;iz++)
			{
				int check = fscanf(fepsy,"%d",&tempint);
				if (check<=0)
				{
					cout << endl << ERR << "Error in reading from external surface";
					exit(-1);
				}
				//delphi->EPSMAP(ix,iy,iz,1,NX,NY,NZ)=tempint;
				write4DVector<int>(delphi->epsmap,tempint,ix,iy,iz,1,NX,NY,NZ,3);
			}
	
	if (ix!=NX)
	{
		cout << endl << WARN << "Error reading epsmapy!";
		return false;
	}
	fclose(fepsy);


	for (ix=0;ix<NX;ix++)
		for (iy=0;iy<NY;iy++)
			for (iz=0;iz<NZ;iz++)
			{
				int check = fscanf(fepsz,"%d",&tempint);
				if (check<=0)
				{
					cout << endl << ERR << "Error in reading from external surface";
					exit(-1);
				}
				//delphi->EPSMAP(ix,iy,iz,2,NX,NY,NZ)=tempint;
				write4DVector<int>(delphi->epsmap,tempint,ix,iy,iz,2,NX,NY,NZ,3);
			}
	
	if (ix!=NX)
	{
		cout << endl << WARN << "Error reading epsmapz!";
		return false;
	}
	fclose(fepsz);
	
	int check = fscanf(fproj,"%d",&nbgp);
	if (check<=0)
	{
		cout << endl << ERR << "Error in reading from external surface";
		exit(-1);
	}

	if (nbgp<=0)
	{
		cout << endl << WARN << "Null number of bgps!";
		return false;
	}

	delphi->nbgp = nbgp;

	delphi->scspos	= allocateVector<double>(3*delphi->nbgp);
	delphi->scsnor	= allocateVector<double>(3*delphi->nbgp);
	delphi->scsarea = allocateVector<double>(delphi->nbgp);
	delphi->ibgp	= allocateVector<int>(3*delphi->nbgp);

	int ind[3];
	double bgpdata[6];
	// load bgp data
	for (int i=0;i<nbgp;i++)
	{
		int check = fscanf(fproj,"%d %d %d %lf %lf %lf %lf %lf %lf",&(ind[0]),&(ind[1]),
			&(ind[2]),&(bgpdata[0]),&(bgpdata[1]),&(bgpdata[2]),&(bgpdata[3]),
			&(bgpdata[4]),&(bgpdata[5]));

		if (check<=0)
		{
			cout << endl << ERR << "Error in reading from external surface";
			exit(-1);
		}

		ind[0]-=1;
		ind[1]-=1;
		ind[2]-=1;

		if (ind[0]<0 || ind[1]<0 || ind[2]<0)
		{
			cout << endl << WARN << "Wrong bgp indexing!";
			return false;
		}

		delphi->ibgp[3*i]=ind[0];
		delphi->ibgp[3*i+1]=ind[1];
		delphi->ibgp[3*i+2]=ind[2];

		delphi->scspos[3*i]=bgpdata[0];
		delphi->scspos[3*i+1]=bgpdata[1];
		delphi->scspos[3*i+2]=bgpdata[2];

		delphi->scsnor[3*i]=bgpdata[3];
		delphi->scsnor[3*i+1]=bgpdata[4];
		delphi->scsnor[3*i+2]=bgpdata[5];
	}
	fclose(fproj);

	// load status
	for (ix=0;ix<NX;ix++)
	{
		for (iy=0;iy<NY;iy++)
		{
			for (iz=0;iz<NZ;iz++)
			{
				int check = fscanf(fstatus,"%d",&tempint);
				if (check<=0)
				{
					cout << endl << ERR << "Error in reading from external surface";
					exit(-1);
				}
				//delphi->status[ix][iy][iz]=((char)tempint);
				//delphi->STATUSMAP(ix,iy,iz,NX,NY)=((short)tempint);
				write3DVector<short>(delphi->status,(short)tempint,ix,iy,iz,NX,NY,NZ);
			}				
		}
	}
	
	if (ix!=NX)
	{
		cout << endl << WARN << "Error reading status map!";
		return false;
	}

	fclose(fstatus);
	cout << endl << INFO << "Number bgps..." << delphi->nbgp;		

	if (fillCav)
	{
		int cav = getCavities();
		cout << endl << INFO << "Detected " << cav << " cavitiy[ies]";		
		fillCavities(vol);
		cout << endl << INFO << "Filtering bgps...";		
	
		// remove false bgp from bgp list; this filter is due to conditional filling
		// of the cavity detector
		int index = 0;
		for (int i=0;i<nbgp;i++)
		{
			ix=delphi->ibgp[3*i];
			iy=delphi->ibgp[3*i+1];
			iz=delphi->ibgp[3*i+2];

			//int kost = delphi->EPSMAP(ix,iy,iz,0,NX,NY,NZ);
			const int konst = read4DVector<int>(delphi->epsmap,ix,iy,iz,0,NX,NY,NZ,3);
			if ((ix-1<0) || (iy-1<0) || (iz-1<0))
				continue;

			// true bgp is saved
			/*
			if (kost != delphi->EPSMAP(ix,iy,iz,1,NX,NY,NZ) ||
				kost != delphi->EPSMAP(ix,iy,iz,2,NX,NY,NZ) ||
				kost != delphi->EPSMAP((ix-1),iy,iz,0,NX,NY,NZ) ||
				kost != delphi->EPSMAP(ix,(iy-1),iz,1,NX,NY,NZ) ||
				kost != delphi->EPSMAP(ix,iy,(iz-1),2,NX,NY,NZ) )
			*/

			if (konst != read4DVector<int>(delphi->epsmap,ix,iy,iz,1,NX,NY,NZ,3) ||
				konst != read4DVector<int>(delphi->epsmap,ix,iy,iz,2,NX,NY,NZ,3) ||
				konst != read4DVector<int>(delphi->epsmap,ix-1,iy,iz,0,NX,NY,NZ,3) ||
				konst != read4DVector<int>(delphi->epsmap,ix,iy-1,iz,1,NX,NY,NZ,3) ||
				konst != read4DVector<int>(delphi->epsmap,ix,iy,iz-1,2,NX,NY,NZ,3) )
			{

				delphi->ibgp[3*index]=delphi->ibgp[3*i];
				delphi->ibgp[3*index+1]=delphi->ibgp[3*i+1];
				delphi->ibgp[3*index+2]=delphi->ibgp[3*i+2];

				delphi->scspos[3*index]=delphi->scspos[3*i];
				delphi->scspos[3*index+1]=delphi->scspos[3*i+1];
				delphi->scspos[3*index+2]=delphi->scspos[3*i+2];

				delphi->scsnor[3*index]=delphi->scsnor[3*i];
				delphi->scsnor[3*index+1]=delphi->scsnor[3*i+1];
				delphi->scsnor[3*index+2]=delphi->scsnor[3*i+2];			

				index++;

				/*printf("\n %d %d %d %d %d %d (%d %d %d)",delphi->EPSMAP(ix,iy,iz,0,NX,NY,NZ),delphi->EPSMAP(ix,iy,iz,1,NX,NY,NZ),
					delphi->EPSMAP(ix,iy,iz,2,NX,NY,NZ),delphi->EPSMAP((ix-1),iy,iz,0,NX,NY,NZ),delphi->EPSMAP(ix,(iy-1),iz,1,NX,NY,NZ),
					delphi->EPSMAP(ix,iy,(iz-1),2,NX,NY,NZ),ix,iy,iz);
				getchar();*/
			}
			else
			{
				/*
				printf("\n Discarding bgp..");
				printf("\n %d %d %d %d %d %d",delphi->EPSMAP(ix,iy,iz,0,NX,NY,NZ),delphi->EPSMAP(ix,iy,iz,1,NX,NY,NZ),
					delphi->EPSMAP(ix,iy,iz,2,NX,NY,NZ),delphi->EPSMAP((ix-1),iy,iz,0,NX,NY,NZ),delphi->EPSMAP(ix,(iy-1),iz,1,NX,NY,NZ),
					delphi->EPSMAP(ix,iy,(iz-1),2,NX,NY,NZ));
				*/
			}
		}
		delphi->ibgp = (int*)realloc(delphi->ibgp,3*sizeof(int)*index);
		delphi->scspos = (double*)realloc(delphi->scspos,3*sizeof(double)*index);
		delphi->scsnor = (double*)realloc(delphi->scsnor,3*sizeof(double)*index);

		delphi->nbgp = index;
		cout << endl << INFO << "Number bgps after conditional filling..." << delphi->nbgp;		
	}

	/** TODO salt for now not supported, always false*/
	for (int i=0;i<NX;i++)
	  for (int j=0;j<NY;j++)
	    for (int k=0;k<NZ;k++)
	      //delphi->IDEBMAP(i,j,k,NX,NY)=false;
		  write3DVector<bool>(delphi->idebmap,false,i,j,k,NX,NY,NZ);

	return true;
}

	
inline void ExternalSurface::getRayIntersection(double pa[3],double pb[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals)
{
	intersections.clear();
	cout << endl << WARN << "Cannot perform ray intersection with an externally loaded surface";
}					

bool ExternalSurface::getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3)
{
	cout << endl << WARN << "Cannot perform projection with an externally loaded surface";
	return false;
}
