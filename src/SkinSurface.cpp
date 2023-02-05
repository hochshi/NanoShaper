#include "SkinSurface.h"

void SkinSurface::clear()
{
	if (gridMixedCellMap!=NULL)	
		deleteVector<int>(gridMixedCellMap);
	if (x!=NULL)
		deleteVector<double>(x);
	if (y!=NULL)
		deleteVector<double>(y);
	if (z!=NULL)
		deleteVector<double>(z);	
	if (ind!=NULL)
		deleteMatrix3D<unsigned short>(nx,ny,ind);

	if (mixedComplex!=NULL)
	{
		for (int i=0;i<numMixedCells;i++)
			delete mixedComplex[i];

		free(mixedComplex);
	}

	// deallocate acceleration grid 2D
	if (gridMixedCellMap2D!=NULL)
		deleteVector<int>(gridMixedCellMap2D);

	if (ind_2d!=NULL)
		deleteMatrix2D<unsigned int>(last_cols_ind,ind_2d);

	for (unsigned int i=0;i<pendingCells.size();i++)
		delete pendingCells[i];
}

SkinSurface::~SkinSurface()
{
	clear();
}

void SkinSurface::init()
{
	gridMixedCellMap = NULL;
	gridMixedCellMap2D = NULL;
	mixedComplex = NULL;
	numMixedCells = 0;
	ind = NULL;
	ind_2d = NULL;
	x = NULL;
	y = x;
	z = y;
	s =0.5;
	savePovRay = false;
	// aggressive settings can be used in a 64 bit machine with 6 GB of memory
	AUX_GRID_DIM_SKIN = 100;
	MAX_MIXEDCELLS = 400;
	AUX_GRID_DIM_SKIN_2D = 50; // aggressive setting is 150
	MAX_MIXEDCELLS_2D = (400*AUX_GRID_DIM_SKIN_2D); // aggressive setting is (200*AUX_GRID_DIM_SKIN_2D)
	surfType = MOLECULAR_SURFACE;
	fastProjection = false;
	providesAnalyticalNormals = true;
}

SkinSurface::SkinSurface(DelPhiShared* ds):Surface()
{
	init();
	// set environment
	delphi = ds;	
}

SkinSurface::SkinSurface():Surface()
{
	delphi = NULL;
	init();
}

void SkinSurface::init(ConfigFile* cf)
{
	double skin_s = cf->read<double>( "Skin_Surface_Parameter", 0.45 );
	unsigned int maxSkinDim = cf->read<unsigned int>( "Max_skin_patches_auxiliary_grid_size", 100 );
	unsigned int maxSkinPatches = cf->read<unsigned int>( "Max_skin_patches_per_auxiliary_grid_cell", 400 );
	unsigned int maxSkinDim2D = cf->read<unsigned int>( "Max_skin_patches_auxiliary_grid_2d_size", 50 );
	unsigned int maxSkinPatches2D = cf->read<unsigned int>( "Max_skin_patches_per_auxiliary_grid_2d_cell", 400 );
	bool useFastProjection = cf->read<bool>("Skin_Fast_Projection",false);
	bool savePovRay = cf->read<bool>("Save_PovRay",false);

	setAuxGrid(maxSkinDim,maxSkinPatches);
	setAuxGrid2D(maxSkinDim2D,maxSkinPatches2D);
	setFastProjection(useFastProjection);
	inside = 5;	
	setShrinking(skin_s);	
	setSavePovRay(savePovRay);
}

SkinSurface::SkinSurface(ConfigFile* cf,DelPhiShared* ds):Surface(cf)
{
	init();
	init(cf);
	// set environment
	delphi = ds;	
}


void SkinSurface::postRayCasting()
{
	// remove 2d grid for ray casting 
	if (gridMixedCellMap2D!=NULL)
		deleteVector<int>(gridMixedCellMap2D);	
	
	if (ind_2d!=NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,ind_2d);

	gridMixedCellMap2D = NULL;
	ind_2d = NULL;
}

bool SkinSurface::preBoundaryProjection()
{
	// 3d grid is necessary only for boundary grid points projection
	if (projBGP)
		return buildAuxiliaryGrid();
	return false;
}

bool SkinSurface::build()
{
	bool f=false;
	#ifdef ENABLE_CGAL
		f = buildSkinCGAL();
	#else
		cout << endl << ERR << "Skin surface cannot be used, please install CGAL and rebuild\n";
		exit(-1);
	#endif
	if (!f)
	{
		cout << endl << ERR << "Error during skin build-up";
	}
	return f;
}

#ifdef ENABLE_CGAL
bool SkinSurface::buildSkinCGAL()
{
	vector<Weighted_point> l;
	float x,y,z,r;
	float max_x=-1e6,min_x=1e6,max_y=-1e6,min_y=1e6,max_z=-1e6,min_z=1e6;
	// plane orientation test
	double sign;
	double RM_[4][4];
	double T_[4][4];
	double TEMP_[4][4];
	double C_[4][4];
	double Cinv_[4][4];
	double M_[4][4];
	double Q_[4][4];

	// clean mixed complex
	if (mixedComplex!=NULL)
	{
		for (int i=0;i<numMixedCells;i++)
			delete mixedComplex[i];

		free(mixedComplex);
	}

	l.reserve(delphi->numAtoms+6);

	time_t start, end;
	time(&start);	

	cout << endl << INFO << "Building skin surface..";
			
	for (int i=0;i<delphi->numAtoms;i++)
	{
		delphi->atoms[i]->pos[0] = delphi->atoms[i]->pos[0]+randDisplacement*(randnum()-0.5);
		delphi->atoms[i]->pos[1] = delphi->atoms[i]->pos[1]+randDisplacement*(randnum()-0.5);
		delphi->atoms[i]->pos[2] = delphi->atoms[i]->pos[2]+randDisplacement*(randnum()-0.5);

		x = (float)delphi->atoms[i]->pos[0];
		y = (float)delphi->atoms[i]->pos[1];
		z = (float)delphi->atoms[i]->pos[2];
		r = (float)delphi->atoms[i]->radius;
		
		l.push_back(Weighted_point(Point3(x,y,z),(r*r)/s,i)); 
		max_x = max(max_x,x+r);
		max_y = max(max_y,y+r);
		max_z = max(max_z,z+r);

		min_x = std::min(min_x,x-r);
		min_y = std::min(min_y,y-r);
		min_z = std::min(min_z,z-r);	
	}
	
	// Regular Triangulation object 
	Rt rT;

	float mid_x = (max_x+min_x)/2.f;
	float mid_y = (max_y+min_y)/2.f;
	float mid_z = (max_z+min_z)/2.f;

	min_x -= fabs(mid_x-min_x)*2;
	min_y -= fabs(mid_y-min_y)*2;
	min_z -= fabs(mid_z-min_z)*2;

	max_x += fabs(mid_x-max_x)*2;
	max_y += fabs(mid_y-max_y)*2;
	max_z += fabs(mid_z-max_z)*2;


	vector<MixedCell*> tempMixedComplex;
	tempMixedComplex.reserve(4*delphi->numAtoms);

	// add bounding box using -1 weights to let easy detection of these virtual atoms
	l.push_back(Weighted_point(Point3(min_x,mid_y,mid_z),-1,delphi->numAtoms));
	l.push_back(Weighted_point(Point3(max_x,mid_y,mid_z),-1,delphi->numAtoms+1));
	l.push_back(Weighted_point(Point3(mid_x,min_y,mid_z),-1,delphi->numAtoms+2));
	l.push_back(Weighted_point(Point3(mid_x,max_y,mid_z),-1,delphi->numAtoms+3));
	l.push_back(Weighted_point(Point3(mid_x,mid_y,min_z),-1,delphi->numAtoms+4));
	l.push_back(Weighted_point(Point3(mid_x,mid_y,max_z),-1,delphi->numAtoms+5));

	cout << endl << INFO << "Regular triangulation....";

	rT.insert (l.begin(), l.end());

	assert( rT.is_valid() );
	assert( rT.dimension() == 3 );

	cout << "ok!";

	cout << endl << INFO << "Computing 3-Delaunay patches...";
  
	int currentCell=0;
	int num_d3_v0_patches = 0;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// Compute Voronoi points and reduced Delaunay solids (Skin) ///////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Finite_Cells_Iterator fcit = rT.finite_cells_begin(); 

	Del0Cell** atomPatches = NULL;
	atomPatches = allocateVector<Del0Cell*>(delphi->numAtoms);

	for (;fcit!=rT.finite_cells_end();fcit++,currentCell++)    
	{
		const Point3& p = rT.geom_traits().construct_weighted_circumcenter_3_object()(fcit->vertex(0)->point(),fcit->vertex(1)->point(),fcit->vertex(2)->point(),fcit->vertex(3)->point());
		Del3Cell*& mc = (Del3Cell*&)fcit->info();
		mc = new Del3Cell();
		mc->surface_type = DELAUNAY_TETRA_CELL;	
		mc->vor[0] = p.x();
		mc->vor[1] = p.y();
		mc->vor[2] = p.z();
		double* vor = mc->vor;

		// minkowsky sum of tethraedron vertices and voronoi point. Record the atom from which
		// it was computed
		double a,b,c;
		double radius2=INFINITY;

		double cellPoints[4][3];

		for (int i=0;i<4;i++)
		{
			Weighted_point& wp = fcit->vertex(i)->point();

			mc->ids[i] = wp.index();
			mc->reduced[i][0] = wp.x();
			mc->reduced[i][1] = wp.y();
			mc->reduced[i][2] = wp.z();

			mc->points[i][0] = wp.x();
			mc->points[i][1] = wp.y();
			mc->points[i][2] = wp.z();

			cellPoints[i][0]  = wp.x();
			cellPoints[i][1]  = wp.y();
			cellPoints[i][2]  = wp.z();

			a = vor[0]-mc->reduced[i][0];				
			mc->reduced[i][0] = mc->reduced[i][0]+s*a;
				
			b = vor[1]-mc->reduced[i][1];
			mc->reduced[i][1] = mc->reduced[i][1]+s*b;

			c = vor[2]-mc->reduced[i][2];
			mc->reduced[i][2] = mc->reduced[i][2]+s*c;
				
			if (i==0) 
			{	
				double norm2 = a*a+b*b+c*c;
				radius2 = -(wp.weight()-norm2);
			}
		}

		// avoid storing non sense patch, but keep its voronoi point
		if (radius2<0)
		{
			pendingCells.push_back(mc);
			// avoid deletion
			mc->quadric = NULL;
			continue;
		}

		tempMixedComplex.push_back(mc);

		////////////////////////// store equation and simplex ////////////////////////////////////////////
					
		plane3points(mc->reduced[2],mc->reduced[1],mc->reduced[0],mc->planes[0],false);
		plane3points(mc->reduced[3],mc->reduced[2],mc->reduced[0],mc->planes[1],false);
		plane3points(mc->reduced[0],mc->reduced[1],mc->reduced[3],mc->planes[3],false);
		plane3points(mc->reduced[3],mc->reduced[1],mc->reduced[2],mc->planes[2],false);

		plane3points(mc->points[2],mc->points[1],mc->points[0],mc->planes_or[0],false);
		plane3points(mc->points[3],mc->points[2],mc->points[0],mc->planes_or[1],false);
		plane3points(mc->points[0],mc->points[1],mc->points[3],mc->planes_or[3],false);
		plane3points(mc->points[3],mc->points[1],mc->points[2],mc->planes_or[2],false);

		mc->mc_points.push_back(mc->reduced[0]);
		mc->mc_points.push_back(mc->reduced[1]);
		mc->mc_points.push_back(mc->reduced[2]);
		mc->mc_points.push_back(mc->reduced[3]);

		// NOTE:
		// these planes are already correctly oriented because the choosen
		// order of points in producing the planes leverages the combinatoric properties CGAL internal data structures.
		// In a few words, using points on a tethraedron in that order one is assured to get correctly oriented planes
		// In a non CGAL implementation these planes must be oriented

		double* QQ = allocateVector<double>(10);
		mc->quadric=QQ;
	
		QQ[0]=1;
		QQ[1]=1;		
		QQ[2]=1;		
		QQ[3]= -radius2*(1-s)+vor[0]*vor[0]+vor[1]*vor[1]+vor[2]*vor[2];		
		QQ[4]=0;		
		QQ[7]=0;		
		QQ[9]=-vor[0];		
		QQ[5]=0;
		QQ[8]=-vor[1];
		QQ[6]=-vor[2];

		////////////////////////// end store equation and simplex ////////////////////////////////////////////	
	}

	
	cout << "ok!";
    currentCell--;
	type[DELAUNAY_TETRA_CELL] = currentCell;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// Compute Voronoi (VDW) or reduced Voronoi solids (Skin) //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << endl << INFO << "Computing 0-Delaunay patches...";

	vector<Cell_handle> cells;
	cells.reserve(1000);

	Finite_Vertex_Iterator fvit = rT.finite_vertices_begin();
	int cellCount = 0;
	for(;fvit!=rT.finite_vertices_end();fvit++)
	{	  
		const Weight& w0 = fvit->point().weight();
		
		// skip virtual atoms
		if ((w0)==-1)
			continue;
	  
		int numPoints = 0;

		if (rT.is_infinite(fvit->cell()))
			continue;

		// current atom id around which we are moving
		const int& currentId = fvit->point().index();
		
		cells.clear();
		rT.incident_cells(fvit,std::back_inserter(cells));

		bool infiniteCell = false;

		// start check if all thetraedra are feasible
		std::vector<Cell_handle>::iterator it;
		for (it=cells.begin();it!=cells.end();it++)
		{
			if (rT.is_infinite((*it)))
			{
				infiniteCell = true;
				break;
			}
		}
	 
		if (infiniteCell)
			continue;

		// this cell is finite, it is acceptable
		cellCount++;
	  
		Del0Cell* mc = new Del0Cell();
		tempMixedComplex.push_back(mc);
		mc->id = currentId;
		mc->surface_type = DELAUNAY_POINT_CELL;

		atomPatches[currentId] = mc;

		double* QQ = allocateVector<double>(10);
		mc->quadric=QQ;
	
		QQ[0]=1;
		QQ[1]=1;		
		QQ[2]=1;		
		QQ[3]= -s*(fvit->point().weight())+fvit->point().x()*fvit->point().x()+fvit->point().y()*fvit->point().y()+fvit->point().z()*fvit->point().z();
		QQ[4]=0;		
		QQ[7]=0;		
		QQ[9]=-fvit->point().x();		
		QQ[5]=0;
		QQ[8]=-fvit->point().y();
		QQ[6]=-fvit->point().z();

		////////////////////////// end store equation and simplex ////////////////////////////////////////////
	  
		numPoints = 0;
	}
  
	cout << "ok!";
	type[DELAUNAY_POINT_CELL] = cellCount;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// Compute Delaunay Edges + Voronoi facets /////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << endl << INFO << "Computing 1-Delaunay patches...";

	// Skin Surface needs the remaining solids
	int num_d1_v2_patches = 0;

	vector<double*> upperPoints;
	vector<double*> lowerPoints;
	vector<double*> vorFace;

	vorFace.reserve(50);
	upperPoints.reserve(50);
	lowerPoints.reserve(50);

	// Delaunay Edge vs Voronoi facets
	Finite_Edges_Iterator fei = rT.finite_edges_begin();
	for (;fei!=rT.finite_edges_end();fei++)
	{
		vector<Cell_handle> cells;
		Cell_circulator cit = rT.incident_cells(*fei);
		int num = (int)circulator_size(cit);

		// up and down edge vertices
		const Vertex_handle& v1 = fei->first->vertex(fei->second);
		const Vertex_handle& v2 = fei->first->vertex(fei->third);

		// connecting one or two non atoms
		if (v1->point().weight()==-1 || v2->point().weight()==-1)
			continue;

		const int refAtom1 = v1->point().index();
		const int refAtom2 = v2->point().index();

		double v2p[3],v1p[3];
		
		v2p[0] = v2->point().x();
		v2p[1] = v2->point().y();
		v2p[2] = v2->point().z();

		v1p[0] = v1->point().x();
		v1p[1] = v1->point().y();
		v1p[2] = v1->point().z();

		vorFace.clear();
		upperPoints.clear();
		lowerPoints.clear();

		Del1Cell* mc = new Del1Cell();
		mc->ids[0] = refAtom1;
		mc->ids[1] = refAtom2;

		// collect all the points
		for (int i=0;i<num;i++)
		{
			const Cell_handle& ch = cit;
			const Del3Cell* ci = ch->info();

			double* p1=NULL;
			
			// search the point with the same id of the upper atom
			for (int j=0;j<4;j++)
			{
				if (ci->ids[j]==refAtom1)
				{
					p1 = (double*)ci->reduced[j];
					break;
				}
			}

			if (p1==NULL)
			{
				for (int j=0;j<4;j++)
					cout << endl << ERR << ci->ids[j];

				cout << endl << ERR << "Index not found";
				cout << endl << ERR << refAtom1;

				for (unsigned int i=0; i<tempMixedComplex.size();i++)
					delete tempMixedComplex[i];
				return false;
			}

			// store the pointer to the reduced point just identified
			mc->mc_points.push_back(p1);
			upperPoints.push_back(p1);
			// we have the three points to compute the upper facet
			// remember that the first point is up, the second down, the third up and so on
			if (i==2)
			{
				plane3points(upperPoints[0],upperPoints[1],upperPoints[2],mc->upper);
				// now the corresponding reduced voronoi cell must link to this plane
				if (refAtom1<delphi->numAtoms)
					atomPatches[refAtom1]->planes.push_back(mc->upper);
			}
			
			double * p2=NULL;
			// the same for the inferior point
			// search the point with the same id
			for (int j=0;j<4;j++)
			{
				if (ci->ids[j]==refAtom2)
				{
					p2 = (double*)ci->reduced[j];
					break;
				}
			}
			
			if (p2==NULL)
			{
				for (int j=0;j<4;j++)
					cout << endl << ERR << ci->ids[j];

				cout << endl << ERR << "Index not found";
				cout << endl << ERR << refAtom2;
				for (unsigned int i=0; i<tempMixedComplex.size();i++)
					delete tempMixedComplex[i];
				return false;
			}
			
			mc->mc_points.push_back(p2);
			lowerPoints.push_back(p2);
			// we have the three points to compute the lower facet
			// remember that the first point is up, the second down, the third up and so on
			if (i==2)
			{
				plane3points(lowerPoints[0],lowerPoints[1],lowerPoints[2],mc->lower);
				// now the corresponding reduced voronoi cell must link to this plane
				if (refAtom2<delphi->numAtoms)
					atomPatches[refAtom2]->planes.push_back(mc->lower);
			}

			//printf("\n vor %x",ci->vor);
			// store the vornoi points to infer voronoi original plane
			// this is later used to compute the focus of the cell
			vorFace.push_back((double*)ci->vor);

			cit++;
		}

		// plane line intersection point to get focus		
		double mat[3][3];
		double invmat[3][3];

		mat[0][0] = v2p[0]-v1p[0];
		mat[1][0] = v2p[1]-v1p[1];
		mat[2][0] = v2p[2]-v1p[2];

		mat[0][1] =	vorFace[1][0]-vorFace[0][0];
		mat[1][1] = vorFace[1][1]-vorFace[0][1];
		mat[2][1] = vorFace[1][2]-vorFace[0][2];

		mat[0][2] =	vorFace[2][0]-vorFace[0][0];
		mat[1][2] = vorFace[2][1]-vorFace[0][1];
		mat[2][2] = vorFace[2][2]-vorFace[0][2];

		double det;

		INVERT_3X3(invmat,det,mat)

		double v[3];
		
		v[0] = v2p[0]-vorFace[0][0];
		v[1] = v2p[1]-vorFace[0][1];
		v[2] = v2p[2]-vorFace[0][2];

		double res=0;
		for (int j=0;j<3;j++)
			res+=invmat[0][j]*v[j];

		double focus[3];

		focus[0] = v2p[0] + res*(v1p[0]-v2p[0]);
		focus[1] = v2p[1] + res*(v1p[1]-v2p[1]);
		focus[2] = v2p[2] + res*(v1p[2]-v2p[2]);

		// generate a reference inner point
		double inner[3],t1[3],t2[3];
		int n = 2;

		t1[0] = (upperPoints[0][0]+upperPoints[1][0]+upperPoints[n][0])/3.;
		t1[1] = (upperPoints[0][1]+upperPoints[1][1]+upperPoints[n][1])/3.;
		t1[2] = (upperPoints[0][2]+upperPoints[1][2]+upperPoints[n][2])/3.;

		t2[0] = (lowerPoints[0][0]+lowerPoints[1][0]+lowerPoints[n][0])/3.;
		t2[1] = (lowerPoints[0][1]+lowerPoints[1][1]+lowerPoints[n][1])/3.;
		t2[2] = (lowerPoints[0][2]+lowerPoints[1][2]+lowerPoints[n][2])/3.;

		inner[0] = (t1[0]+t2[0])/2.;
		inner[1] = (t1[1]+t2[1])/2.;
		inner[2] = (t1[2]+t2[2])/2.;

		// decide orientation based on the position of the focus of the computed planes
		sign = DOT(mc->upper,inner)+mc->upper[3];
		if (sign>0) ASSIGN4(mc->upper,-mc->upper)
		
		sign = DOT(mc->lower,inner)+mc->lower[3];
		if (sign>0) ASSIGN4(mc->lower,-mc->lower)
		
		// build lateral planes
		for (int i=0;i<num-1;i++)
		{
			double* plane = allocateVector<double>(4);
			plane3points(upperPoints[i],lowerPoints[i],upperPoints[i+1],plane);
			// orient
			sign = DOT(plane,inner)+plane[3];
			if (sign>0) ASSIGN4(plane,-plane)
			mc->planes.push_back(plane);
		}
		
		// add the last one		
		double* plane = allocateVector<double>(4);
		plane3points(upperPoints[num-1],lowerPoints[num-1],upperPoints[0],plane);
		// orient
		
		sign = DOT(plane,inner)+plane[3];
		if (sign>0) ASSIGN4(plane,-plane)
		
		mc->planes.push_back(plane);

		// compute the reduced cell radius
		double a = v1p[0]-focus[0];
		double b = v1p[1]-focus[1];
		double c = v1p[2]-focus[2];
		double radius2 = v1->point().weight()-(a*a+b*b+c*c);
		
		num_d1_v2_patches++;
		
		double o2[3],o1[3],o3[3],ttt;

		// vector parallel to voronoi face
		SUB(o2,vorFace[0],focus)
		NORMALIZE(o2,ttt)

		// vector parallel to delaunay edge
		SUB(o1,v2p,v1p)
		NORMALIZE(o1,ttt)

		CROSS(o3,o1,o2)		
		NORMALIZE(o3,ttt)

		// o1,o2,o3 is an orthonormal set

		for (int ii=0;ii<4;ii++)
			for (int jj=0;jj<4;jj++)
			{
				C_[ii][jj]		=	0;
				Cinv_[ii][jj]	=	0;
				M_[ii][jj]		=	0;
				Q_[ii][jj]		=	0;
				TEMP_[ii][jj]	=	0;
				RM_[ii][jj]		=	0;
				T_[ii][jj]		=	0;
			}


		RM_[0][0] = o1[0];
		RM_[1][0] = o1[1];
		RM_[2][0] = o1[2];

		RM_[0][1] = o3[0];
		RM_[1][1] = o3[1];
		RM_[2][1] = o3[2];

		RM_[0][2] = o2[0];
		RM_[1][2] = o2[1];
		RM_[2][2] = o2[2];
	 
		RM_[3][3] = 1.;

		T_[0][0] = 1.;
		T_[1][1] = 1.;
		T_[2][2] = 1.;
		T_[0][3] = focus[0];
		T_[1][3] = focus[1];
		T_[2][3] = focus[2];
		T_[3][3] = 1.;
	 		
		Matrix4x4MultiplyBy4x4(T_,RM_,Cinv_);

		for (int ii=0;ii<4;ii++)
			for (int jj=0;jj<4;jj++)
				C_[ii][jj]=Cinv_[ii][jj];

		inplace_invert4x4(Cinv_);
		
		M_[0][0] = -1./(1.-s);
		M_[1][1] = 1./s;
		M_[2][2] = 1./s;
		M_[3][3] = -radius2;

		Matrix4x4MultiplyBy4x4(M_,Cinv_,TEMP_);

		double temp;
		for (int ii = 0; ii<4; ii++)
			for (int jj = ii+1; jj<4; jj++)
			{
				temp=Cinv_[ii][jj];
				Cinv_[ii][jj]=Cinv_[jj][ii];
				Cinv_[jj][ii]=temp;
			}
  

		Matrix4x4MultiplyBy4x4(Cinv_,TEMP_,Q_);

		////////////////////////// store equation and simplex ////////////////////////////////////////////
			
		tempMixedComplex.push_back(mc);
		mc->surface_type = DELAUNAY_EDGE_CELL;
		
		double* QQ = allocateVector<double>(10);
		mc->quadric=QQ;
		
		// Symmetric matrix coding
		//	Q0	Q4	Q7	Q9
		//		Q1	Q5	Q8
		//			Q2	Q6
		//				Q3
		
		QQ[0]=Q_[0][0];
		QQ[1]=Q_[1][1];		
		QQ[2]=Q_[2][2];		
		QQ[3]=Q_[3][3];
		QQ[4]=Q_[0][1];		
		QQ[7]=Q_[0][2];		
		QQ[9]=Q_[0][3];		
		QQ[5]=Q_[1][2];
		QQ[8]=Q_[1][3];
		QQ[6]=Q_[2][3];

				
		////////////////////////// end store equation and simplex ////////////////////////////////////////////

	}

	type[DELAUNAY_EDGE_CELL] = num_d1_v2_patches;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// Compute Delaunay Facets + Voronoi Edges /////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout << "ok!";
	int num_d2_v1_patches = 0;

	cout << endl << INFO << "Computing 2-Delaunay patches...";

	upperPoints.reserve(3);
	lowerPoints.reserve(3);

	Finite_Facets_Iterator ffi = rT.finite_facets_begin();
	for (;ffi!=rT.finite_facets_end();ffi++)
	{
		const Cell_handle& ch = ffi->first;
		// get the mirror facet/cell go on if infinite
		const Facet& mirrorFace = rT.mirror_facet(*ffi);		
		const Cell_handle& ch2 = mirrorFace.first;

		Del3Cell* ci = ch->info();
		Del3Cell* ci2 = ch2->info();

		// it is the local index of the tethraedron vertex
		int f_index1 = ffi->second;
		int f_index2 = mirrorFace.second;

		double* vor1 = ci->vor;
		double* vor2 = ci2->vor;

		if (rT.is_infinite(ch2) || rT.is_infinite(ch))
			continue;

		Del2Cell* mc = new Del2Cell();		
		mc->surface_type = DELAUNAY_FACET_CELL;
		
		int ind = 0;

		upperPoints.clear();
		lowerPoints.clear();

		double weight;
		double delFace[3][3];

		// recover the three reduced common vertices
		for (int i=0;i<4;i++)
		{		
			// this vertex is opposite and cannot be shared, so skip
			if (i==f_index1)
				continue;

			const int index1 = ch->vertex(i)->point().index();

			mc->ids[ind]=index1;
			delFace[ind][0] = ch->vertex(i)->point().x();
			delFace[ind][1] = ch->vertex(i)->point().y();
			delFace[ind][2] = ch->vertex(i)->point().z();

			// store the weight of the del point that will be used to compute the focus
			if (ind==0)
				weight = ch->vertex(i)->point().weight();				
			
			ind++;

			double* p1 = NULL;
			double* p2 = NULL;

			for (int j=0;j<4;j++)
			{
				// search the corresponding reduced point
				if (ci->ids[j]==index1)
				{
					p1 = ci->reduced[j];
					break;
				}
			}

			// search the corresponding reduced point into the other tethraedron
			for (int j=0;j<4;j++)
			{
				if (ci2->ids[j]==index1)
				{
					p2 = ci2->reduced[j];
					break;
				}
			}

			// weird..
			if (p2==NULL || p1==NULL)			
			{
				cout << endl << ERR << "Two adjacent tethraedra does not share the same 3 atoms!";
				for (unsigned int i=0; i<tempMixedComplex.size();i++)
					delete tempMixedComplex[i];
				return false;
			}	

			mc->mc_points.push_back(p1);
			upperPoints.push_back(p1);
			mc->mc_points.push_back(p2);
			lowerPoints.push_back(p2);
		}

		// plane line intersection point to get focus
		double mat[3][3];
		double invmat[3][3];
		
		mat[0][0] = vor2[0]-vor1[0];
		mat[1][0] = vor2[1]-vor1[1];
		mat[2][0] = vor2[2]-vor1[2];

		mat[0][1] =	delFace[1][0]-delFace[0][0];
		mat[1][1] = delFace[1][1]-delFace[0][1];
		mat[2][1] = delFace[1][2]-delFace[0][2];

		mat[0][2] =	delFace[2][0]-delFace[0][0];
		mat[1][2] = delFace[2][1]-delFace[0][1];
		mat[2][2] = delFace[2][2]-delFace[0][2];

		double det;
		INVERT_3X3(invmat,det,mat)

		double v[3];
		
		v[0] = vor2[0]-delFace[0][0];
		v[1] = vor2[1]-delFace[0][1];
		v[2] = vor2[2]-delFace[0][2];

		double res=0;
		for (int j=0;j<3;j++)
			res+=invmat[0][j]*v[j];

		double focus[3];

		// ok this mc really exists and it is acceptable
		tempMixedComplex.push_back(mc);
		focus[0] = vor2[0] + res*(vor1[0]-vor2[0]);
		focus[1] = vor2[1] + res*(vor1[1]-vor2[1]);
		focus[2] = vor2[2] + res*(vor1[2]-vor2[2]);

		num_d2_v1_patches++;
		
		// compute the reduced cell radius
		double a = delFace[0][0]-focus[0];
		double b = delFace[0][1]-focus[1];
		double c = delFace[0][2]-focus[2];
		double radius2 = weight-(a*a+b*b+c*c);

		// generate a reference inner point
		double inner[3],t1[3],t2[3];

		t1[0] = (upperPoints[0][0]+upperPoints[1][0]+upperPoints[2][0])/3.;
		t1[1] = (upperPoints[0][1]+upperPoints[1][1]+upperPoints[2][1])/3.;
		t1[2] = (upperPoints[0][2]+upperPoints[1][2]+upperPoints[2][2])/3.;

		t2[0] = (lowerPoints[0][0]+lowerPoints[1][0]+lowerPoints[2][0])/3.;
		t2[1] = (lowerPoints[0][1]+lowerPoints[1][1]+lowerPoints[2][1])/3.;
		t2[2] = (lowerPoints[0][2]+lowerPoints[1][2]+lowerPoints[2][2])/3.;

		inner[0] = (t1[0]+t2[0])/2.;
		inner[1] = (t1[1]+t2[1])/2.;
		inner[2] = (t1[2]+t2[2])/2.;

		// compute the del facets planes and orient
		plane3points(mc->mc_points[0],mc->mc_points[2],mc->mc_points[4],mc->upper);
		// orient
		sign = DOT(mc->upper,inner)+mc->upper[3];
		if (sign>0) ASSIGN4(mc->upper,-mc->upper)
		
		plane3points(mc->mc_points[1],mc->mc_points[3],mc->mc_points[5],mc->lower);
		// orient
		sign = DOT(mc->lower,inner)+mc->lower[3];
		if (sign>0) ASSIGN4(mc->lower,-mc->lower)

		// compute the lateral planes and orient
		double* plane;
		plane = mc->planes[0];
		plane3points(mc->mc_points[0],mc->mc_points[1],mc->mc_points[2],plane);
		// orient
		sign = DOT(plane,inner)+plane[3];
		if (sign>0) ASSIGN4(plane,-plane)

		plane = mc->planes[1];
		plane3points(mc->mc_points[2],mc->mc_points[3],mc->mc_points[4],plane);
		// orient
		sign = DOT(plane,inner)+plane[3];
		if (sign>0) ASSIGN4(plane,-plane)

		plane = mc->planes[2];
		plane3points(mc->mc_points[4],mc->mc_points[5],mc->mc_points[0],plane);
		// orient
		sign = DOT(plane,inner)+plane[3];
		if (sign>0) ASSIGN4(plane,-plane)


		double o1[3],o2[3],o3[3],ttt;
		
		// vector parallel to voronoi edge
		SUB(o1,vor1,vor2)
		NORMALIZE(o1,ttt)

		// vector parallel to delaunay face
		SUB(o2,delFace[0],delFace[1])
		NORMALIZE(o2,ttt)

		// cross product
		CROSS(o3,o1,o2)
		NORMALIZE(o3,ttt)

		// o1,o2,o3 is an orthonormal set
		// compute equation

		for (int ii=0;ii<4;ii++)
			for (int jj=0;jj<4;jj++)
			{
				C_[ii][jj]		=	0;
				Cinv_[ii][jj]	=	0;
				M_[ii][jj]		=	0;
				Q_[ii][jj]		=	0;
				TEMP_[ii][jj]	=	0;
				RM_[ii][jj]		=	0;
				T_[ii][jj]		=	0;
			}


		RM_[0][0] = o2[0];
		RM_[1][0] = o2[1];
		RM_[2][0] = o2[2];

		RM_[0][1] = o3[0];
		RM_[1][1] = o3[1];
		RM_[2][1] = o3[2];

		RM_[0][2] = o1[0];
		RM_[1][2] = o1[1];
		RM_[2][2] = o1[2];
	 
		RM_[3][3] = 1.;

		T_[0][0] = 1.;
		T_[1][1] = 1.;
		T_[2][2] = 1.;
		T_[0][3] = focus[0];
		T_[1][3] = focus[1];
		T_[2][3] = focus[2];
		T_[3][3] = 1.;
	 
		Matrix4x4MultiplyBy4x4(T_,RM_,Cinv_);

		for (int ii=0;ii<4;ii++)
			for (int jj=0;jj<4;jj++)
				C_[ii][jj]=Cinv_[ii][jj];

		inplace_invert4x4(Cinv_);

		M_[0][0] = -1./(1.-s);
		M_[1][1] = -1./(1.-s);
		M_[2][2] = 1./s;
		M_[3][3] = -radius2;

		Matrix4x4MultiplyBy4x4(M_,Cinv_,TEMP_);

		double temp;
		for (int ii = 0; ii<4; ii++)
			for (int jj = ii+1; jj<4; jj++)
			{
				temp=Cinv_[ii][jj];
				Cinv_[ii][jj]=Cinv_[jj][ii];
				Cinv_[jj][ii]=temp;
			}
  

		Matrix4x4MultiplyBy4x4(Cinv_,TEMP_,Q_);

		////////////////////////// store equation and simplex ////////////////////////////////////////////
					
		double* QQ = allocateVector<double>(10);
		mc->quadric=QQ;
		
		/* Symmetric matrix coding
			Q0	Q4	Q7	Q9
				Q1	Q5	Q8
					Q2	Q6
						Q3
		*/
		QQ[0]=Q_[0][0];
		QQ[1]=Q_[1][1];		
		QQ[2]=Q_[2][2];		
		QQ[3]=Q_[3][3];
		QQ[4]=Q_[0][1];		
		QQ[7]=Q_[0][2];		
		QQ[9]=Q_[0][3];		
		QQ[5]=Q_[1][2];
		QQ[8]=Q_[1][3];
		QQ[6]=Q_[2][3];

		////////////////////////// end store equation and simplex ////////////////////////////////////////////

	}

	type[DELAUNAY_FACET_CELL] = num_d2_v1_patches;

	// assign pointers, consolidate mixed complex.
	numMixedCells = (int)tempMixedComplex.size();
	
	mixedComplex = allocateVector<MixedCell*>(numMixedCells);	
	int k=0,empty=0;
	for(vector<MixedCell*>::iterator it = tempMixedComplex.begin();tempMixedComplex.end()!=it;it++,k++)	
		mixedComplex[k] = (*it);

	cout << "ok!";
	double duration;
	time(&end);	
	duration = difftime(end, start);
	cout << endl << INFO << "Surface build-up time.. " << duration << " [s]"; 	

	// remove references to atom patches
	free(atomPatches);
	
	if (savePovRay)
	{
		ofstream of;
		cout << endl << INFO << "Saving surface in Pov-Ray in skin.pov...";
		cout.flush();
		of.open("skin.pov");	
		of << "#include \"shapes.inc\" ";
		of << "\n#include \"colors.inc\" ";
		of << "\nglobal_settings {max_trace_level 3}";
		of << "\nbackground { color Black }"; 
		of << "\ncamera {";
		of << "\n\tlocation  <" << mid_x << ","  << max_y <<"," << mid_z <<">";
		of << "\nlook_at  <" << mid_x << "," << mid_y << "," << mid_z << "> translate-<" << mid_x << "," << mid_y << "," << mid_z << ">" << " rotate<0,0,clock> "<< "translate<" << mid_x << "," << mid_y << "," << mid_z << ">}";		
		of << "\nlight_source {<" << mid_x << ","  << max_y <<"," << mid_z <<">" << " color White "<< "translate-<" << mid_x << "," << mid_y << "," << mid_z << ">" << " rotate<0,0,clock> "<< "translate<" << mid_x << "," << mid_y << "," << mid_z << ">}";		
		for(int k=0;k<numMixedCells;k++)	
		{			
			MixedCell* mc = mixedComplex[k];
			saveSkinPatch(of,mc,k,l);
		}
		cout << "ok!";
		of.close();
	}
		
	printSummary();
	return true;
}
#endif

void SkinSurface::saveSkinPatch(ofstream& of,MixedCell* mc,int index,vector<Weighted_point>& la)
{
	char clipplanes[1000],buff2[1000];

	if (mc->surface_type == DELAUNAY_EDGE_CELL)
	{	
		Del1Cell* ec = (Del1Cell*)mc;
		of << "\n// ---------------------------------------------//";
		of << "\n// ------------- Delunay edge cell -------------//";
		of << "\n// ---------------------------------------------//";
		of << "\n// atoms " << ec->ids[0] << "," << ec->ids[1];

		sprintf(clipplanes,"Clipping_Mesh%d",index);		
		of << "\n#declare " << clipplanes << "=\n mesh{";

		// this is slow but allows Pov-Ray to work well
		Polyhedron poly;	  
		// little noise to avoid problems and to let povray work
		vector<Point3> local;
		for (unsigned int i=0;i<ec->mc_points.size();i++)
		{
			double* p = ec->mc_points[i];			
			local.push_back(Point3(p[0]+1e-3*(randnum()-0.5),p[1]+1e-3*(randnum()-0.5),p[2]+1e-3*(randnum()-0.5)));
		}

		// compute convex hull to recover mesh triangles
		CGAL::convex_hull_3(local.begin(), local.end(), poly);	

		Polyhedron::Facet_iterator fit;
		HF_circulator h;		
		// save the mesh 
		for(fit=poly.facets_begin();fit!=poly.facets_end();fit++)
		{	
			of << "\n triangle{\n";			
			int g=0;
			for (h = fit->facet_begin();g<3;h++)
			{
				of << "<" << h->vertex()->point().x() << "," << h->vertex()->point().y() << "," << h->vertex()->point().z() << ">";
				if (g<2) of << ",";
				g++;
			}
			of << "}";
		}

		// generate a reference inner point
		double inner[3],t1[3],t2[3];
		
		t1[0] = (ec->mc_points[0][0]+ec->mc_points[2][0]+ec->mc_points[4][0])/3.;
		t1[1] = (ec->mc_points[0][1]+ec->mc_points[2][1]+ec->mc_points[4][1])/3.;
		t1[2] = (ec->mc_points[0][2]+ec->mc_points[2][2]+ec->mc_points[4][2])/3.;

		t2[0] = (ec->mc_points[1][0]+ec->mc_points[3][0]+ec->mc_points[5][0])/3.;
		t2[1] = (ec->mc_points[1][1]+ec->mc_points[3][1]+ec->mc_points[5][1])/3.;
		t2[2] = (ec->mc_points[1][2]+ec->mc_points[3][2]+ec->mc_points[5][2])/3.;

		inner[0] = (t1[0]+t2[0])/2.;
		inner[1] = (t1[1]+t2[1])/2.;
		inner[2] = (t1[2]+t2[2])/2.;

		of  << "\n inside_vector<" << inner[0] << "," << inner[1] << "," << inner[2] << ">";
		of << "}";

			
		// this is formally correct but Pov-Ray seems not to be able to manage this code for 
		// more than very small molecules
		/*
		// add clipping planes
		of << "\n\n #declare " << clipplanes << " =" ;		
		sprintf(buff3,"\n\n // clipping planes \n\n intersection { \n ");
		of << buff3;
		
		for (unsigned int l=0;l<ec->planes.size();l++)
		{
			sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",ec->planes[l][0],ec->planes[l][1],ec->planes[l][2],-(ec->planes[l][3])/(sqrt((DOT(ec->planes[l],ec->planes[l])))));				
			of << buff3;
		}

		sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",ec->lower[0],ec->lower[1],ec->lower[2],-(ec->lower[3])/(sqrt((DOT(ec->lower,ec->lower)))));				
		of << buff3;

		sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",ec->upper[0],ec->upper[1],ec->upper[2],-(ec->upper[3])/(sqrt((DOT(ec->upper,ec->upper)))));				
		of << buff3;
		of << "\n}";
		*/
		
		double A,B,C,D,E,F,G,H,I,J;
		A = mc->quadric[0];		
		B = mc->quadric[1];		
		C = mc->quadric[2];		
		J = mc->quadric[3];
		D = 2*mc->quadric[4];		
		E = 2*mc->quadric[7];		
		G = 2*mc->quadric[9];		
		F = 2*mc->quadric[5];
		H = 2*mc->quadric[8];
		I = 2*mc->quadric[6];

		sprintf(buff2,"\n\n quadric { \n <%f,%f,%f>,<%f,%f,%f>,<%f,%f,%f>,%f",A,B,C,D,E,F,G,H,I,J);
		of << buff2;
		of << "\n pigment{color Yellow}";
		of << "\n bounded_by { " << clipplanes << " } \n clipped_by{ bounded_by}}"; 
		of << "\n";
	}
	else if (mc->surface_type == DELAUNAY_POINT_CELL)
	{		
		Del0Cell* pc = (Del0Cell*)mc;
		of << "\n// ---------------------------------------------//";
		of << "\n// ------------- Delunay point cell ------------//";
		of << "\n// ---------------------------------------------//";
		of << "\n// atom " << pc->id;

		sprintf(buff2,"\n\n sphere { \n <%f,%f,%f>,%f",la[pc->id].x(),la[pc->id].y(),la[pc->id].z(),sqrt(la[pc->id].weight()*s));	  						
		of << buff2;
		of << "\n pigment{color Green}}";
	}
	else if (mc->surface_type == DELAUNAY_FACET_CELL)
	{
		Del2Cell* fc = (Del2Cell*)mc;
		of << "\n// ----------------------------------------------//";
		of << "\n// ------------- Delunay facet cell -------------//";
		of << "\n// ----------------------------------------------//";
		of << "\n// atoms " << fc->ids[0] << "," << fc->ids[1] << "," << fc->ids[2];

		sprintf(clipplanes,"Clipping_Mesh%d",index);		
		of << "\n#declare " << clipplanes << "=\n mesh{";

		// this is slow but allows Pov-Ray to work well
		Polyhedron poly;	  
		// little noise to avoid problems
		vector<Point3> local;
		for (unsigned int i=0;i<fc->mc_points.size();i++)
		{
			double* p = fc->mc_points[i];			
			local.push_back(Point3(p[0]+1e-3*(randnum()-0.5),p[1]+1e-3*(randnum()-0.5),p[2]+1e-3*(randnum()-0.5)));
		}

		// compute convex hull to recover mesh triangles
		CGAL::convex_hull_3(local.begin(), local.end(), poly);	

		Polyhedron::Facet_iterator fit;
		HF_circulator h;		
		// save the mesh 
		for(fit=poly.facets_begin();fit!=poly.facets_end();fit++)
		{	
			of << "\n triangle{\n";			
			int g=0;
			for (h = fit->facet_begin();g<3;h++)
			{
				of << "<" << h->vertex()->point().x() << "," << h->vertex()->point().y() << "," << h->vertex()->point().z() << ">";
				if (g<2) of << ",";
				g++;
			}
			of << "}";
		}

		// generate a reference inner point
		double inner[3],t1[3],t2[3];
		
		t1[0] = (fc->mc_points[0][0]+fc->mc_points[2][0]+fc->mc_points[4][0])/3.;
		t1[1] = (fc->mc_points[0][1]+fc->mc_points[2][1]+fc->mc_points[4][1])/3.;
		t1[2] = (fc->mc_points[0][2]+fc->mc_points[2][2]+fc->mc_points[4][2])/3.;

		t2[0] = (fc->mc_points[1][0]+fc->mc_points[3][0]+fc->mc_points[5][0])/3.;
		t2[1] = (fc->mc_points[1][1]+fc->mc_points[3][1]+fc->mc_points[5][1])/3.;
		t2[2] = (fc->mc_points[1][2]+fc->mc_points[3][2]+fc->mc_points[5][2])/3.;

		inner[0] = (t1[0]+t2[0])/2.;
		inner[1] = (t1[1]+t2[1])/2.;
		inner[2] = (t1[2]+t2[2])/2.;

		of  << "\n inside_vector<" << inner[0] << "," << inner[1] << "," << inner[2] << ">";
		of << "}";

		// this is formally correct but Pov-Ray seems not to be able to manage this code for 
		// more than very small molecules
		/*
		// add clipping planes
		of << "\n\n #declare " << clipplanes << " =" ;		
		sprintf(buff3,"\n\n // clipping planes \n\n intersection { \n ");
		of << buff3;
		
		for (int l=0;l<3;l++)
		{
			sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",fc->planes[l][0],fc->planes[l][1],fc->planes[l][2],-(fc->planes[l][3])/(sqrt((DOT(fc->planes[l],fc->planes[l])))));				
			of << buff3;
		}
		
		sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",fc->lower[0],fc->lower[1],fc->lower[2],-(fc->lower[3])/(sqrt((DOT(fc->lower,fc->lower)))));				
		of << buff3;

		sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",fc->upper[0],fc->upper[1],fc->upper[2],-(fc->upper[3])/(sqrt((DOT(fc->upper,fc->upper)))));				
		of << buff3;
		of << "\n}";
		*/

		double A,B,C,D,E,F,G,H,I,J;
		A = mc->quadric[0];		
		B = mc->quadric[1];		
		C = mc->quadric[2];		
		J = mc->quadric[3];
		D = 2*mc->quadric[4];		
		E = 2*mc->quadric[7];		
		G = 2*mc->quadric[9];		
		F = 2*mc->quadric[5];
		H = 2*mc->quadric[8];
		I = 2*mc->quadric[6];

		sprintf(buff2,"\n\n quadric { \n <%f,%f,%f>,<%f,%f,%f>,<%f,%f,%f>,%f",A,B,C,D,E,F,G,H,I,J);
		of << buff2;
		of << "\n pigment{color Magenta}";
		of << "\n bounded_by { " << clipplanes << " } \n clipped_by{ bounded_by}}"; 
		of << "\n";

	}
	else if (mc->surface_type == DELAUNAY_TETRA_CELL)
	{		
		Del3Cell* tc = (Del3Cell*)mc;
		of << "\n// ----------------------------------------------------//";
		of << "\n// ------------- Delunay tethraedron cell -------------//";
		of << "\n// ----------------------------------------------------//";
		of << "\n// atoms " << tc->ids[0] << "," << tc->ids[1] << "," << tc->ids[2] << "," << tc->ids[3];
		
		sprintf(clipplanes,"Clipping_Mesh%d",index);		
		of << "\n#declare " << clipplanes << "=\n mesh{";

		// this is slow but allows Pov-Ray to work well
		Polyhedron poly;	  
		// little noise to avoid problems
		vector<Point3> local;
		double pp[3];
		for (unsigned int i=0;i<4;i++)
		{
			pp[0] = tc->reduced[i][0];
			pp[1] = tc->reduced[i][1];
			pp[2] = tc->reduced[i][2];
			local.push_back(Point3(pp[0]+1e-3*(randnum()-0.5),pp[1]+1e-3*(randnum()-0.5),pp[2]+1e-3*(randnum()-0.5)));
		}

		// compute convex hull to recover mesh triangles
		CGAL::convex_hull_3(local.begin(), local.end(), poly);	

		Polyhedron::Facet_iterator fit;
		HF_circulator h;		
		// save the mesh 
		for(fit=poly.facets_begin();fit!=poly.facets_end();fit++)
		{	
			of << "\n triangle{\n";			
			int g=0;
			for (h = fit->facet_begin();g<3;h++)
			{
				of << "<" << h->vertex()->point().x() << "," << h->vertex()->point().y() << "," << h->vertex()->point().z() << ">";
				if (g<2) of << ",";
				g++;
			}
			of << "}";
		}

		of  << "\n inside_vector<" << tc->vor[0] << "," << tc->vor[1] << "," << tc->vor[2] << ">";
		of << "}";

		// this is formally correct but Pov-Ray seems not to be able to manage this code for 
		// more than very small molecules
		/*
		sprintf(clipplanes,"Clipping_Planes%d",index);
		of << "\n\n #declare " << clipplanes << " =" ;	
		of << "\n intersection {";
		for (int l=0;l<4;l++)
		{
			sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",tc->planes[l][0],tc->planes[l][1],tc->planes[l][2],-(tc->planes[l][3])/(sqrt((DOT(tc->planes[l],tc->planes[l])))));				
			of << buff3;
		}
		of << "\n}";
		*/
		double radius2,v[3],norm2,w,radius;
		v[0] = la[tc->ids[0]].x()-tc->vor[0];
		v[1] = la[tc->ids[0]].y()-tc->vor[1];
		v[2] = la[tc->ids[0]].z()-tc->vor[2];
		norm2 = DOT(v,v);
		w = la[tc->ids[0]].weight();
		radius2 = -(w-norm2);
		radius = sqrt(radius2*(1-s));
		sprintf(buff2,"\n\n sphere { \n <%f,%f,%f>,%f",tc->vor[0],tc->vor[1],tc->vor[2],radius);
		of << buff2;
		of << "\n pigment{color Red}";
		of << "\n bounded_by { " << clipplanes << " } \n clipped_by{bounded_by}}";
	}		
}

void SkinSurface::preProcessPanel()
{
	unsigned int igrid=delphi->nx;	
	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption

	// auxiliary grid is consistent with delphi grid in order to speed-up tracing
	int gridMul=1;
	while(igrid>AUX_GRID_DIM_SKIN_2D)
	{
		gridMul+=2;
		int digrid = delphi->nx;
		while(1)
		{
			// get nearest odd multiple
			int fixedPoint = (digrid+(gridMul-1))/gridMul;
			igrid = fixedPoint*gridMul;
			if (igrid%2==0)
				digrid=igrid+1;
			else
			{
				igrid = fixedPoint;
				break;
			}
		}
	}
	// auxiliary scale
	scale_2d = delphi->scale/((double)gridMul);

	//cout << endl << INFO << "Auxiliary grid is " << igrid;
	
	xmin_2d = delphi->baricenter[0]-(igrid-1)/(2*scale_2d);
	ymin_2d = delphi->baricenter[1]-(igrid-1)/(2*scale_2d);
	zmin_2d = delphi->baricenter[2]-(igrid-1)/(2*scale_2d);

	xmax_2d = delphi->baricenter[0]+(igrid-1)/(2*scale_2d);
	ymax_2d = delphi->baricenter[1]+(igrid-1)/(2*scale_2d);
	zmax_2d = delphi->baricenter[2]+(igrid-1)/(2*scale_2d);

	side_2d = 1/scale_2d;
	nx_2d = igrid;
	ny_2d = igrid;
	nz_2d = igrid;

	if (gridMixedCellMap2D!=NULL)
		deleteVector<int>(gridMixedCellMap2D);	
	
	if (ind_2d!=NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,ind_2d);

	if (panel==0)
	{
		last_rows_ind = ny_2d;
		last_cols_ind = nz_2d;
	}
	else if (panel==1)
	{
		last_rows_ind = nx_2d;
		last_cols_ind = ny_2d;
	}
	else
	{
		last_rows_ind = nx_2d;
		last_cols_ind = nz_2d;	
	}
	
	gridMixedCellMap2D = allocateVector<int>(last_rows_ind*last_cols_ind*MAX_MIXEDCELLS_2D);

	ind_2d = allocateMatrix2D<unsigned int>(last_rows_ind,last_cols_ind);
	for (int i=0;i<last_rows_ind;i++)
		for (int j=0;j<last_cols_ind;j++)			
			ind_2d[i][j]=0;

	
	// build a bounding box for cell and map it to
	// the auxiliary grid
	unsigned int max_t = 0;

	if (mixedComplex == NULL)
	{	
		cout << endl << WARN << "Cannot get surface without a computed mixed complex!";
		return;
	}	

	for (int it=0;it<numMixedCells;it++)
	{
		// mixed cell points
		vector<double*>& ll = mixedComplex[it]->mc_points;		 

		/*
		if (ll.size()==0)
		{
			cout << endl << ERR << "Empty mixed cell!";
			exit(-1);
		}
		*/

		// compute bounding box object
		double downx = INFINITY;
		double downy = -INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = INFINITY;
		double upz = -INFINITY;

		// use the cube of the atom as bounding box object for voronoi cells
		if (mixedComplex[it]->surface_type == DELAUNAY_POINT_CELL)
		{
			int index = ((Del0Cell*)mixedComplex[it])->id;
			
			// avoid mapping a non atom (bounding box virtual atom)
			if (index>=delphi->numAtoms)
			{				
				continue;
			}

			double radius = delphi->atoms[index]->radius;
			double* sphere_center = delphi->atoms[index]->pos;

			downx = sphere_center[0]-radius;
			downy = sphere_center[1]+radius;
			downz = sphere_center[2]-radius;

			upx = sphere_center[0]+radius;
			upy = sphere_center[1]-radius;
			upz = sphere_center[2]+radius;		
		}
		else
		{
			for (unsigned int pind=0;pind<ll.size();pind++)
			{
				double* pp = ll[pind];
				downx = MIN(downx,pp[0]);
				downy = MAX(downy,pp[1]);
				downz = MIN(downz,pp[2]);

				upx = MAX(upx,pp[0]);
				upy = MIN(upy,pp[1]);
				upz = MAX(upz,pp[2]);
			}
		}     
		// resolve which are the grid cubes that
		// are cut by the bounding box object. These
		// cubes see the mixed cell bounding box    
		int ix_start = (int)rintp((downx-xmin_2d)/side_2d);
		int iy_start = (int)rintp((downy-ymin_2d)/side_2d);
		int iz_start = (int)rintp((downz-zmin_2d)/side_2d);
	     	     
		int ix_end = (int)rintp((upx-xmin_2d)/side_2d);
		int iy_end = (int)rintp((upy-ymin_2d)/side_2d);
		int iz_end = (int)rintp((upz-zmin_2d)/side_2d);

		if (ix_start<0)
		 ix_start=0;
		if (iy_start<0)
		 iy_start=0;
		if (iz_start<0)
		 iz_start=0;

		if (ix_end<0)
		 ix_end=0;
		if (iy_end<0)
		 iy_end=0;
		if (iz_end<0)
		 iz_end=0;

		if (ix_end>=nx_2d)
		 ix_end=nx_2d-1;
		if (iy_end>=ny_2d)
		 iy_end=ny_2d-1;	          
		if (iz_end>=nz_2d)
		 iz_end=nz_2d-1;	          

		if (ix_start>=nx_2d)
		 ix_start=nx_2d-1;
		if (iy_start>=ny_2d)
		 iy_start=ny_2d-1;	          
		if (iz_start>=nz_2d)
		 iz_start=nz_2d-1;	          

		// map the points into the 2D matrices of each plane
		
		if (panel==0)
		{
			// plane yz
			for (int iz=iz_start;iz<=iz_end;iz++)
			{
				for (int iy=iy_start;iy>=iy_end;iy--)			 
				{
					if (ind_2d[iy][iz]>max_t)
						max_t =ind_2d[iy][iz];			
				
					if (ind_2d[iy][iz]>=MAX_MIXEDCELLS_2D)
					{
						cout << endl << ERR << "Number of mixed cells is superior to maximum allowed, please increase Max_skin_patches_per_auxiliary_grid_2d_cell";						
						exit(-1);
					}
					GRID_MIXEDCELLMAP_2D(iy,iz,(ind_2d[iy][iz]),ny_2d,nz_2d) = it;
					ind_2d[iy][iz]++;
				}
			}
		}

		else if (panel==1)
		{
			// plane xy
			for (int ix=ix_start;ix<=ix_end;ix++)
			{
				for (int iy=iy_start;iy>=iy_end;iy--)			 
				{
					if (ind_2d[ix][iy]>max_t)
						max_t =ind_2d[ix][iy];			
				
					if (ind_2d[ix][iy]>=MAX_MIXEDCELLS_2D)
					{
						cout << endl << ERR << "Number of mixed cells is superior to maximum allowed, please increase Max_skin_patches_per_auxiliary_grid_2d_cell";						
						exit(-1);
					}
					GRID_MIXEDCELLMAP_2D(ix,iy,(ind_2d[ix][iy]),nx_2d,ny_2d) = it;
					ind_2d[ix][iy]++;
				}
			}
		}
		else
		{
			// plane xz
			for (int ix=ix_start;ix<=ix_end;ix++)
			{
				for (int iz=iz_start;iz<=iz_end;iz++)
				{
					if (ind_2d[ix][iz]>max_t)
						max_t =ind_2d[ix][iz];			
				
					if (ind_2d[ix][iz]>=MAX_MIXEDCELLS_2D)
					{
						cout << endl << ERR << "Number of mixed cells is superior to maximum allowed, please increase Max_skin_patches_per_auxiliary_grid_2d_cell";						
						exit(-1);
					}
					GRID_MIXEDCELLMAP_2D(ix,iz,(ind_2d[ix][iz]),nx_2d,nz_2d) = it;
					ind_2d[ix][iz]++;
				}
			}
		}
	}    	
}

bool SkinSurface::buildAuxiliaryGrid()
{
	// Perform pre-processing to speed-up intersections and projections.
	// Compute bounding box for each mixed cell in the surface, and map each
	// bounding box to the proper grid point.
	// This routine uses an auxiliary grid. Delphi grid is not directly used
	// because it may be too memory consuming. To speed up computations
	// a maximal number of mixed cells is allowed in each auxiliary grid cell
	// the macro is MAX_MIXEDCELLS in SkinSurface.h.
	// The macro AUX_GRID_DIM_SKIN sets the maximally allowed grid size.
	unsigned int igrid=delphi->nx;	
	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption

	/*
	double layer = 0.0;

	if (igrid>AUX_GRID_DIM_SKIN)
		igrid = AUX_GRID_DIM_SKIN;					

	scale = ((double)(igrid-1))*delphi->perfill/(100.*(delphi->rmaxdim+layer));*/

	// auxiliary grid is consistent with delphi grid in order to speed-up tracing
	int gridMul=1;
	while(igrid>AUX_GRID_DIM_SKIN)
	{
		gridMul+=2;
		int digrid = delphi->nx;
		while(1)
		{
			// get nearest odd multiple
			int fixedPoint = (digrid+(gridMul-1))/gridMul;
			igrid = fixedPoint*gridMul;
			if (igrid%2==0)
				digrid=igrid+1;
			else
			{
				igrid = fixedPoint;
				break;
			}
		}
	}
	// auxiliary scale
	scale = delphi->scale/((double)gridMul);

	cout << endl << INFO << "Auxiliary grid is " << igrid;
	
	xmin = delphi->baricenter[0]-(igrid-1)/(2*scale);
	ymin = delphi->baricenter[1]-(igrid-1)/(2*scale);
	zmin = delphi->baricenter[2]-(igrid-1)/(2*scale);

	xmax = delphi->baricenter[0]+(igrid-1)/(2*scale);
	ymax = delphi->baricenter[1]+(igrid-1)/(2*scale);
	zmax = delphi->baricenter[2]+(igrid-1)/(2*scale);

	////////////// Allocate memory for the grid, planes and maps //////////////////////

	if (x!=NULL)
		delete x;
	if (y!=NULL)
		delete y;
	if (z!=NULL)
		delete z;

	if (ind!=NULL)
		deleteMatrix3D<unsigned short>(nx,ny,ind);

	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);

	side = 1/scale;
	nx = igrid;
	ny = igrid;
	nz = igrid;

	if (gridMixedCellMap!=NULL)
		deleteVector<int>(gridMixedCellMap);	

	cout << endl << INFO << "Allocating " << (nx*ny*nz*MAX_MIXEDCELLS)*sizeof(int)/1024.0/1024.0 << " MB" << " for the auxiliary grid...";
	gridMixedCellMap = allocateVector<int>(nx*ny*nz*MAX_MIXEDCELLS);

	ind = allocateMatrix3D<unsigned short>(nx,ny,nz);
	for (int i=0;i<nx;i++)
		for (int j=0;j<ny;j++)
			for (int k=0;k<nz;k++)
				ind[i][j][k]=0;

	for (int i=0;i<nx;i++)
		x[i]=xmin+i*side;
	
	for (int i=0;i<ny;i++)
		y[i]=ymin+i*side;

	for (int i=0;i<nz;i++)
		z[i]=zmin+i*side;
	
	cout << "ok!";
	//////////////////////////////////////////////////////////////////////////

	if (mixedComplex == NULL)
	{	
		cout << endl << WARN << "Cannot get surface without a computed mixed complex!";
		return false;
	}	
	// build a bounding box for each mixed cell and map it to
	// the auxiliary grid
	int max_t = 0;

	cout << endl << INFO << "Mapping auxiliary grid...";
	
	for (int it=0;it<numMixedCells;it++)
	{
		// mixed cell points
		vector<double*>& ll = mixedComplex[it]->mc_points;		 

		/*
		if (ll.size()==0)
		{
			cout << endl << ERR << "Empty mixed cell!";
			exit(-1);
		}
		*/
		 
		// compute bounding box object
		double downx = INFINITY;
		double downy = -INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = INFINITY;
		double upz = -INFINITY;

		// use the cube of the atom as bounding box object for voronoi cells
		if (mixedComplex[it]->surface_type == DELAUNAY_POINT_CELL)
		{
			int index = ((Del0Cell*)mixedComplex[it])->id;
			
			// avoid mapping a non atom (bounding box virtual atom)
			if (index>=delphi->numAtoms)
				continue;

			double radius = delphi->atoms[index]->radius;
			double* sphere_center = delphi->atoms[index]->pos;

			downx = sphere_center[0]-radius;
			downy = sphere_center[1]+radius;
			downz = sphere_center[2]-radius;

			upx = sphere_center[0]+radius;
			upy = sphere_center[1]-radius;
			upz = sphere_center[2]+radius;		
		}
		else
		{
			for (unsigned int pind=0;pind<ll.size();pind++)
			{
				double* pp = ll[pind];
				downx = MIN(downx,pp[0]);
				downy = MAX(downy,pp[1]);
				downz = MIN(downz,pp[2]);

				upx = MAX(upx,pp[0]);
				upy = MIN(upy,pp[1]);
				upz = MAX(upz,pp[2]);
			}
		}     
		// resolve which are the grid cubes that
		// are cut by the bounding box object. These
		// cubes see the mixed cell bounding box    
		int ix_start = (int)rintp((downx-xmin)/side);
		int iy_start = (int)rintp((downy-ymin)/side);
		int iz_start = (int)rintp((downz-zmin)/side);
	     	     
		int ix_end = (int)rintp((upx-xmin)/side);
		int iy_end = (int)rintp((upy-ymin)/side);
		int iz_end = (int)rintp((upz-zmin)/side);

		if (ix_start<0)
		 ix_start=0;
		if (iy_start<0)
		 iy_start=0;
		if (iz_start<0)
		 iz_start=0;

		if (ix_end<0)
		 ix_end=0;
		if (iy_end<0)
		 iy_end=0;
		if (iz_end<0)
		 iz_end=0;

		if (ix_end>=nx)
		 ix_end=nx-1;
		if (iy_end>=ny)
		 iy_end=ny-1;	          
		if (iz_end>=nz)
		 iz_end=nz-1;	          

		if (ix_start>=nx)
		 ix_start=nx-1;
		if (iy_start>=ny)
		 iy_start=ny-1;	          
		if (iz_start>=nz)
		 iz_start=nz-1;	          

		//printf("\n in ");
		for (int iz=iz_start;iz<=iz_end;iz++)
		 for (int iy=iy_start;iy>=iy_end;iy--)
			 for (int ix=ix_start;ix<=ix_end;ix++)
			 {
				 //printf("-");
				 if (ind[ix][iy][iz]>max_t)
					 max_t =ind[ix][iy][iz];

				 if (ind[ix][iy][iz]>=MAX_MIXEDCELLS)
				 {
					cout << endl << ERR << "Number of mixed cells is superior to maximum allowed, please increase Max_skin_patches_per_auxiliary_grid_cell";						
					exit(-1);
				 }
				 GRIDMIXEDCELLMAP(ix,iy,iz,(ind[ix][iy][iz]),nx,ny,nz) = it;
				 ind[ix][iy][iz]++;
			 }				 

	 //printf(" out");
	}    	

	cout << "ok!";
	cout << endl << INFO << "Max mixed cells per auxiliary cell -> " << max_t;

	return true;
}

/** save in .skin format*/
bool SkinSurface::save(char* fileName)
{
	ofstream fout;
    fout.open(fileName,ios::out);

	cout << endl << INFO << "Writing skin in .skin file format in " << fileName << "...";

	if (fout.fail())
	{
		cout << endl << WARN << "Cannot write file " << fileName;
		return false;
	}

	if (mixedComplex==NULL)
	{
		cout << endl << WARN << "Cannot write null mesh!";
		return false;
	}

	// save TODO
	return true;
}


bool SkinSurface::load(char* fileName) 
{
	int len = (int)strlen(fileName);
	if (len==0)
	{
		cout << WARN << "Cannot load with empty file name!";
		return false;
	}
			
	//cout << endl << INFO << "Loading Skin Surface in file " << fileName << "...";

	// TODO LOAD
	return true;
}

void SkinSurface::printSummary()
{
	cout << endl << INFO << "Shrinking value " << getShrinking();
	if (mixedComplex == NULL)
	{
		cout << endl << WARN << "Skin surface not loaded!";
	}
	else
	{
		cout << endl << INFO << "Number of mixed cells -> " << numMixedCells;
		cout << endl << INFO << "Number of del_point/vor_cell -> " << type[0];
		cout << endl << INFO << "Number of del_edge/vor_facet -> " << type[1];
		cout << endl << INFO << "Number of del_facet/vor_edge -> " << type[2];
		cout << endl << INFO << "Number of del_cell/vor_point -> " << type[3];		

		if (internals!=NULL)
		{
			*internals << endl << "mixedcells " << numMixedCells;
			*internals << endl << "del_point " << type[0];
			*internals << endl << "del_edge " << type[1];
			*internals << endl << "del_facet " << type[2];	
			*internals << endl << "del_cell " << type[3];	
		}
	}
}

bool SkinSurface::isFeasible(MixedCell* mc,double* point)
{	
	if (mc->surface_type == DELAUNAY_POINT_CELL)
	{
		Del0Cell* mc0 = (Del0Cell*)mc;
		int num = (int)mc0->planes.size();		
		for (int it=0;it<num;it++)
		{
			double* plane = mc0->planes[it];
			double val = DOT(plane,point)+plane[3];
			if (val < 0)
				return false;
		}
		return true;
	}
	else if (mc->surface_type == DELAUNAY_EDGE_CELL)
	{
		Del1Cell* mc1 = (Del1Cell*)mc;
		int num = (int)mc1->planes.size();
		for (int it = 0;it<num;it++)
		{
			double* plane = mc1->planes[it];
			double val = DOT(plane,point)+plane[3];
			if (val > 0)
				return false;
		}
		
		double val = DOT(mc1->upper,point)+mc1->upper[3];
		
		if (val > 0)
			return false;

		val = DOT(mc1->lower,point)+mc1->lower[3];
		if (val > 0)
			return false;

		
		return true;
	}
	else if (mc->surface_type == DELAUNAY_FACET_CELL)
	{
		Del2Cell* mc2 = (Del2Cell*)mc;
		for (int it = 0;it<3;it++)
		{
			double* plane = mc2->planes[it];
			double val = DOT(plane,point)+plane[3];
			if (val > 0)
				return false;
		}
		
		double val = DOT(mc2->upper,point)+mc2->upper[3];
		
		if (val > 0)
			return false;

		val = DOT(mc2->lower,point)+mc2->lower[3];

		if (val > 0)
			return false;

		
		return true;
	}
	else if (mc->surface_type == DELAUNAY_TETRA_CELL)
	{
		Del3Cell* mc3 = (Del3Cell*)mc;
		for (int it = 0;it<4;it++)
		{
			double* plane = mc3->planes[it];
			double val = DOT(plane,point)+plane[3];
			if (val > 0)
				return false;
		}

		return true;
	}
	else
	{
		cout << endl << ERR << "Unrecognized patch type!";
		return false;
	}
}

inline bool SkinSurface::rayQuadricIntersection(double* orig,double* dir,double*Q,double* i1,double* i2,int id,double* cache)
{
	double a=0,b=0,c=0;
	
	//for (int i=0;i<4;i++)		
	//	for (int j=0;j<4;j++)
	//		a+=dir[i]*Q[i][j]*dir[j];
	
	//inline using symmetry and caching
	//double q44=Q[3][3],q33=Q[2][2],q22=Q[1][1],q11=Q[0][0],q34=Q[2][3],q24=Q[1][3],q23=Q[1][2],q14=Q[0][3],q13=Q[0][2],q12=Q[0][1];	
	double q44=Q[3],q33=Q[2],q22=Q[1],q11=Q[0],q34=Q[6],q24=Q[8],q23=Q[5],q14=Q[9],q13=Q[7],q12=Q[4];	

	a = cache[0]*q44+cache[1]*q34+cache[2]*q33+cache[3]*q24+cache[4]*q23+cache[5]*q22+cache[6]*q14+cache[7]*q13+cache[8]*q12+cache[9]*q11;

	// bug fix (extremely rare)
	if (a==0)
		return false;

	//for (int i=0;i<4;i++)		
	//	for (int j=0;j<4;j++)
	//		b+=orig[i]*Q[i][j]*dir[j];

	b = cache[20]*q44+cache[21]*q34+cache[22]*q33+cache[23]*q24+cache[24]*q23+cache[25]*q22+cache[26]*q14+cache[27]*q13+cache[28]*q12+cache[29]*q11;

	//for (int i=0;i<4;i++)		
	//	for (int j=0;j<4;j++)
	//		c+=orig[i]*Q[i][j]*orig[j];	

	c = cache[10]*q44+cache[11]*q34+cache[12]*q33+cache[13]*q24+cache[14]*q23+cache[15]*q22+cache[16]*q14+cache[17]*q13+cache[18]*q12+cache[19]*q11;
			
	double delta = b*b-a*c;

	if (delta<0)
		return false;
	else
	{
		delta = sqrt(delta);
		(*i1) = (-b-delta)/a;
		(*i2) = (-b+delta)/a;
	}
	return true;	
}

inline void SkinSurface::getRayIntersection(double pa[3],double pb[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals)
{
	double dir[3];
	double aposx = pa[0],aposy = pa[1],aposz = pa[2];
	dir[0]=pb[0]-aposx;
	dir[1]=pb[1]-aposy;
	dir[2]=pb[2]-aposz;
	double apos[3];
	apos[0]=aposx;
	apos[1]=aposy;
	apos[2]=aposz;
	
	intersections.reserve(2000);

	int iyaux,izaux,ixaux;
	
	int numInt = 0;

	// this cache exploits the particular structure of skin ray-tracing in order to improve speed.
	// all intersections are with quadrics so the ray invariant part is precomputed once and for all.
	double cache[30];

	if (panel==0)
	{
		int t1 = (int)rintp((pa[1]-ymin_2d)/side_2d);
		int t2 = (int)rintp((pa[2]-zmin_2d)/side_2d);

		// direction and origin is constant
		double pat[4],dirt[4];				
		ASSIGN(pat,pa);
		ASSIGN(dirt,dir);
		pat[3]=1;
		dirt[3]=0;

		// load cache with constant-ray info
		// this is the quadric-ray intersection that is constant
		// given the current ray
		double d1=dirt[0],d2=dirt[1],d3=dirt[2],d4=dirt[3];
		cache[0]=d4*d4;
		cache[1]=2*d3*d4;
		cache[2]=d3*d3;
		cache[3]=2*d2*d4;
		cache[4]=2*d2*d3;
		cache[5]=d2*d2;
		cache[6]=2*d1*d4;
		cache[7]=2*d1*d3;
		cache[8]=2*d1*d2;
		cache[9]=d1*d1;

		double o1=pat[0],o2=pat[1],o3=pat[2],o4=pat[3];
		cache[10]=o4*o4;
		cache[11]=2*o3*o4;
		cache[12]=o3*o3;
		cache[13]=2*o2*o4;
		cache[14]=2*o2*o3;
		cache[15]=o2*o2;
		cache[16]=2*o1*o4;
		cache[17]=2*o1*o3;
		cache[18]=2*o1*o2;
		cache[19]=o1*o1;

		cache[20]=d4*o4;
		cache[21]=d3*o4+d4*o3;
		cache[22]=d3*o3;
		cache[23]=d2*o4+d4*o2;
		cache[24]=d2*o3+d3*o2;
		cache[25]=d2*o2;
		cache[26]=d1*o4+d4*o1;
		cache[27]=d1*o3+d3*o1;
		cache[28]=d1*o2+d2*o1;
		cache[29]=d1*o1;

		iyaux = t1;
		izaux = t2;
		
		int numcells = ind_2d[iyaux][izaux];
		for (int iter=0;iter<numcells;iter++)
		{				
			int it = GRID_MIXEDCELLMAP_2D(iyaux,izaux,iter,ny_2d,nz_2d);
			double t1,t2;
			double intPoint[3];

			// ray quadric intersection
			bool ff = rayQuadricIntersection(pat,dirt,mixedComplex[it]->quadric,&t1,&t2,thdID,cache);

			numInt++;
				
			// no intersection
			if (!ff)
				continue;

			intPoint[0]=apos[0]+dir[0]*t1;
			intPoint[1]=apos[1];
			intPoint[2]=apos[2];

			// feasibility test: check the intersection is inside the mixed complex cell			
			bool testFeasibility = isFeasible(mixedComplex[it],intPoint);

			// accept the intersection if the intersection is in the cube
			// and the intersection exists 				
			if(testFeasibility)		
			{
				// normal
				if (computeNormals)
				{
					double* normal = allocateVector<double>(3);
					MixedCell* mc = mixedComplex[it];
					double* Q = mc->quadric;
					normal[0] = 2*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
					normal[1] = 2*(Q[1]*intPoint[1]+Q[4]*intPoint[0]+Q[5]*intPoint[2]+Q[8]);
					normal[2] = 2*(Q[2]*intPoint[2]+Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[6]);
					if (mc->surface_type==DELAUNAY_TETRA_CELL)
					{
						CHANGE_SIGN(normal)
					}
					intersections.push_back(pair<double,double*>(t1,normal));
				}
				else
					intersections.push_back(pair<double,double*>(t1,(double*)NULL));				
			}

			intPoint[0]=apos[0]+dir[0]*t2;
			intPoint[1]=apos[1];
			intPoint[2]=apos[2];

			// feasibility test: check the intersection is inside the mixed complex cell			
			testFeasibility = isFeasible(mixedComplex[it],intPoint);
				
			// accept the intersection if the intersection is in the cube
			// and the intersection exists 				
			if(testFeasibility)
			{
				// normal
				if (computeNormals)
				{
					double* normal = allocateVector<double>(3);
					MixedCell* mc = mixedComplex[it];
					double* Q = mc->quadric;
					normal[0] = 2*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
					normal[1] = 2*(Q[1]*intPoint[1]+Q[4]*intPoint[0]+Q[5]*intPoint[2]+Q[8]);
					normal[2] = 2*(Q[2]*intPoint[2]+Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[6]);					
					if (mc->surface_type==DELAUNAY_TETRA_CELL)
					{
						CHANGE_SIGN(normal)
					}
					intersections.push_back(pair<double,double*>(t2,normal));
				}
				else
					intersections.push_back(pair<double,double*>(t2,(double*)NULL));	
			}
		}
	}
	else if(panel==1)
	{
		int t1 = (int)rintp((pa[1]-ymin_2d)/side_2d);
		int t2 = (int)rintp((pa[0]-xmin_2d)/side_2d);

		// direction and origin is constant
		double pat[4],dirt[4];				
		ASSIGN(pat,pa);
		ASSIGN(dirt,dir);
		pat[3]=1;
		dirt[3]=0;

		// load cache with constant-ray info
		double d1=dirt[0],d2=dirt[1],d3=dirt[2],d4=dirt[3];
		cache[0]=d4*d4;
		cache[1]=2*d3*d4;
		cache[2]=d3*d3;
		cache[3]=2*d2*d4;
		cache[4]=2*d2*d3;
		cache[5]=d2*d2;
		cache[6]=2*d1*d4;
		cache[7]=2*d1*d3;
		cache[8]=2*d1*d2;
		cache[9]=d1*d1;

		double o1=pat[0],o2=pat[1],o3=pat[2],o4=pat[3];
		cache[10]=o4*o4;
		cache[11]=2*o3*o4;
		cache[12]=o3*o3;
		cache[13]=2*o2*o4;
		cache[14]=2*o2*o3;
		cache[15]=o2*o2;
		cache[16]=2*o1*o4;
		cache[17]=2*o1*o3;
		cache[18]=2*o1*o2;
		cache[19]=o1*o1;

		cache[20]=d4*o4;
		cache[21]=d3*o4+d4*o3;
		cache[22]=d3*o3;
		cache[23]=d2*o4+d4*o2;
		cache[24]=d2*o3+d3*o2;
		cache[25]=d2*o2;
		cache[26]=d1*o4+d4*o1;
		cache[27]=d1*o3+d3*o1;
		cache[28]=d1*o2+d2*o1;
		cache[29]=d1*o1;

		// get auxiliary grid index
		iyaux = t1;
		ixaux = t2;

		int numcells = ind_2d[ixaux][iyaux];
		for (int iter = 0;iter<numcells;iter++)
		{				
			int it = GRID_MIXEDCELLMAP_2D(ixaux,iyaux,iter,nx_2d,ny_2d);
			double t1,t2;
			double intPoint[3];

			// ray quadric intersection
			bool ff = rayQuadricIntersection(pat,dirt,mixedComplex[it]->quadric,&t1,&t2,thdID,cache);
				
			// no intersection
			if (!ff)
				continue;

			intPoint[0]=apos[0];
			intPoint[1]=apos[1];
			intPoint[2]=apos[2]+dir[2]*t1;

			// feasibility test: check the intersection is inside the mixed complex cell
			bool testFeasibility = isFeasible(mixedComplex[it],intPoint);

			if(testFeasibility)
			{
				// normal
				if (computeNormals)
				{
					double* normal = allocateVector<double>(3);
					MixedCell* mc = mixedComplex[it];
					double* Q = mc->quadric;
					normal[0] = 2*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
					normal[1] = 2*(Q[1]*intPoint[1]+Q[4]*intPoint[0]+Q[5]*intPoint[2]+Q[8]);
					normal[2] = 2*(Q[2]*intPoint[2]+Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[6]);
					if (mc->surface_type==DELAUNAY_TETRA_CELL)
					{
						CHANGE_SIGN(normal)
					}
					intersections.push_back(pair<double,double*>(t1,normal));
				}
				else
					intersections.push_back(pair<double,double*>(t1,(double*)NULL));
			}

			intPoint[0]=apos[0];
			intPoint[1]=apos[1];
			intPoint[2]=apos[2]+dir[2]*t2;

			// feasibility test: check the intersection is inside the mixed complex cell
			testFeasibility = isFeasible(mixedComplex[it],intPoint);

			if(testFeasibility)
			{
				// normal
				if (computeNormals)
				{
					double* normal = allocateVector<double>(3);
					MixedCell* mc = mixedComplex[it];
					double* Q = mc->quadric;
					normal[0] = 2*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
					normal[1] = 2*(Q[1]*intPoint[1]+Q[4]*intPoint[0]+Q[5]*intPoint[2]+Q[8]);
					normal[2] = 2*(Q[2]*intPoint[2]+Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[6]);
					if (mc->surface_type==DELAUNAY_TETRA_CELL)
					{
						CHANGE_SIGN(normal)
					}
					intersections.push_back(pair<double,double*>(t2,normal));
				}
				else
					intersections.push_back(pair<double,double*>(t2,(double*)NULL));
			}
		}	
	}
	else
	{

		int t1 = (int)rintp((pa[0]-xmin_2d)/side_2d);
		int t2 = (int)rintp((pa[2]-zmin_2d)/side_2d);

		// direction and origin is constant
		double pat[4],dirt[4];				
		ASSIGN(pat,pa);
		ASSIGN(dirt,dir);
		pat[3]=1;
		dirt[3]=0;

		// load cache with constant-ray info
		double d1=dirt[0],d2=dirt[1],d3=dirt[2],d4=dirt[3];
		cache[0]=d4*d4;
		cache[1]=2*d3*d4;
		cache[2]=d3*d3;
		cache[3]=2*d2*d4;
		cache[4]=2*d2*d3;
		cache[5]=d2*d2;
		cache[6]=2*d1*d4;
		cache[7]=2*d1*d3;
		cache[8]=2*d1*d2;
		cache[9]=d1*d1;

		double o1=pat[0],o2=pat[1],o3=pat[2],o4=pat[3];
		cache[10]=o4*o4;
		cache[11]=2*o3*o4;
		cache[12]=o3*o3;
		cache[13]=2*o2*o4;
		cache[14]=2*o2*o3;
		cache[15]=o2*o2;
		cache[16]=2*o1*o4;
		cache[17]=2*o1*o3;
		cache[18]=2*o1*o2;
		cache[19]=o1*o1;

		cache[20]=d4*o4;
		cache[21]=d3*o4+d4*o3;
		cache[22]=d3*o3;
		cache[23]=d2*o4+d4*o2;
		cache[24]=d2*o3+d3*o2;
		cache[25]=d2*o2;
		cache[26]=d1*o4+d4*o1;
		cache[27]=d1*o3+d3*o1;
		cache[28]=d1*o2+d2*o1;
		cache[29]=d1*o1;

		// get auxiliary grid index
		ixaux = t1;
		izaux = t2;

		int numcells = ind_2d[ixaux][izaux];
		for (int iter = 0;iter<numcells;iter++)
		{				
			int it = GRID_MIXEDCELLMAP_2D(ixaux,izaux,iter,nx_2d,nz_2d);
			double t1,t2;
			double intPoint[3];

			// ray quadric intersection
			bool ff = rayQuadricIntersection(pat,dirt,mixedComplex[it]->quadric,&t1,&t2,thdID,cache);
				
			// no intersection
			if (!ff)
				continue;

			intPoint[0]=apos[0];
			intPoint[1]=apos[1]+dir[1]*t1;
			intPoint[2]=apos[2];

			// feasibility test: check the intersection is inside the mixed complex cell
			bool testFeasibility = isFeasible(mixedComplex[it],intPoint);

			if(testFeasibility)
			{
				// normal
				if (computeNormals)
				{
					double* normal = allocateVector<double>(3);
					MixedCell* mc = mixedComplex[it];
					double* Q = mc->quadric;
					normal[0] = 2*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
					normal[1] = 2*(Q[1]*intPoint[1]+Q[4]*intPoint[0]+Q[5]*intPoint[2]+Q[8]);
					normal[2] = 2*(Q[2]*intPoint[2]+Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[6]);
					if (mc->surface_type==DELAUNAY_TETRA_CELL)
					{
						CHANGE_SIGN(normal)
					}
					intersections.push_back(pair<double,double*>(t1,normal));
				}
				else
					intersections.push_back(pair<double,double*>(t1,(double*)NULL));
			}

			intPoint[0]=apos[0];
			intPoint[1]=apos[1]+dir[1]*t2;
			intPoint[2]=apos[2];

			// feasibility test: check the intersection is inside the mixed complex cell
			testFeasibility = isFeasible(mixedComplex[it],intPoint);

			if(testFeasibility)
			{
				// normal
				if (computeNormals)
				{
					double* normal = allocateVector<double>(3);
					MixedCell* mc = mixedComplex[it];
					double* Q = mc->quadric;
					normal[0] = 2*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
					normal[1] = 2*(Q[1]*intPoint[1]+Q[4]*intPoint[0]+Q[5]*intPoint[2]+Q[8]);
					normal[2] = 2*(Q[2]*intPoint[2]+Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[6]);
					if (mc->surface_type==DELAUNAY_TETRA_CELL)
					{
						CHANGE_SIGN(normal)
					}
					intersections.push_back(pair<double,double*>(t2,normal));
				}
				else
					intersections.push_back(pair<double,double*>(t2,(double*)NULL));
			}
		}
	}
	if (intersections.size()>0)
		sort(intersections.begin(), intersections.end(),compKeepIndex);	
}					

bool SkinSurface::getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3)
{
	// get the mixed cells that are associated to this grid point
	// by querying the auxiliary grid
	double dist;
	set<int>cells;
	double hside = delphi->hside;

	double aposx = p[0];
	double aposy = p[1];
	double aposz = p[2];

	// move from delphi grid to auxiliary grid
	int irefx = (int)rintp((aposx-xmin)/side);
	int irefy = (int)rintp((aposy-ymin)/side);
	int irefz = (int)rintp((aposz-zmin)/side);

	int ixaux = irefx;
	int iyaux = irefy;
	int izaux = irefz;

	for (int i=0;i<ind[ixaux][iyaux][izaux];i++)
		cells.insert((GRIDMIXEDCELLMAP(ixaux,iyaux,izaux,i,nx,ny,nz)));		

	// keep the nearest patch
	set<int>::iterator it;
	double locProj[3], minDist=INFINITY, locNorm[3];
	bool ff=false;

	dist = 0;
	locProj[0]=0;
	locProj[1]=0;
	locProj[2]=0;

	locNorm[0]=0;
	locNorm[1]=0;
	locNorm[2]=0;

	for (it=cells.begin();it!=cells.end();it++)
	{						
		double* Q = mixedComplex[(*it)]->quadric;
		int patchType = mixedComplex[(*it)]->surface_type;
		projectToQuadric(p,Q,patchType,locProj,locNorm,dist);
		
		if (dist==INFINITY)
			continue;

		bool flag = isFeasible(mixedComplex[(*it)],locProj);

		if (!flag)
			continue;

		if (dist<minDist)
		{
			minDist = dist;
			(*proj1)=locProj[0];
			(*proj2)=locProj[1];
			(*proj3)=locProj[2];

			(*normal1)=locNorm[0];
			(*normal2)=locNorm[1];
			(*normal3)=locNorm[2];
		}
	}
	if (cells.size()==0 || minDist == INFINITY)
	{
		{  
			#ifdef ENABLE_BOOST_THREADS
				boost::mutex::scoped_lock scopedLock(mutex);
			#endif
			(*errorStream) << endl << WARN << "Approximating bgp with grid point";
		}
		(*proj1)=p[0];
		(*proj2)=p[1];
		(*proj3)=p[2];
		return true;
	}
	
	// further check that the projection is in the cube	
	//if (!testInCube((*proj1),(*proj2),(*proj3),p[0],p[1],p[2],delphi->side))
	//	cout << endl << WARN << "Out of cube projection in SkinSurface::getProjection!";
	
	return true;
}

void SkinSurface::projectToQuadric(double* y,double* Q,int type,double* proj,double* normal,double& dist)
{	
	// sphere case performs fast projection
	if (type==0 || type==3)
	{
		double cc[3];
		cc[0] = -Q[9]/Q[0];
		cc[1] = -Q[8]/Q[0];
		cc[2] = -Q[6]/Q[0];		
		double k = Q[3]/Q[0];
		double r2 = -k+cc[0]*cc[0]+cc[1]*cc[1]+cc[2]*cc[2];

		// immaginary sphere, skip
		if (r2<0)
		{
			dist = INFINITY;		
			return;
		}
		double r = sqrt(r2);
		
		// existing sphere is projected
		double pp[3];
		double ttt;
		pp[0]=0;
		pp[1]=0;
		pp[2]=0;
		SUB(pp,y,cc);
		NORMALIZE_S(pp,ttt);
		ADD_MUL(proj,cc,pp,r);

		// get the normal depending on the insidness value
		double insideness=0;
		double yi,yj;
		for (int i=0;i<4;i++)
		{
			if (i==3) 
				yi = 1;
			else
				yi = y[i];

			for (int j=0;j<4;j++)
			{						
				if (j==3)
					yj = 1;
				else
					yj = y[j];

				double qij = 0;
				
				if (i==0 && j==0)
					qij = Q[0];
				else if (i==1 && j==1)
					qij = Q[1];
				else if (i==2 && j==2)
					qij = Q[2];
				else if (i==3 && j==3)
					qij = Q[3];
				if ((i==1 && j==0) || (i==0 && j==1))
					qij = Q[4];
				else if ((i==2 && j==0) || (i==0 && j==2))
					qij = Q[7];
				else if ((i==3 && j==0) || (i==0 && j==3))
					qij = Q[9];
				else if ((i==1 && j==2) || (i==2 && j==1))
					qij = Q[5];
				else if ((i==2 && j==3) || (i==3 && j==2))
					qij = Q[6];
				else if ((i==1 && j==3) || (i==3 && j==1))
					qij = Q[8];
				
				insideness+=yi*yj*qij;
			}
		}
     
		if (insideness>0)
			insideness = 1;
		else
			insideness = -1;
	     
		if (type== 3)
			insideness = -insideness;		      
			
		DIST(dist,proj,y);
		SUB(normal,y,proj);
		if (insideness<0)
		{
			normal[0]=-normal[0];
			normal[1]=-normal[1];
			normal[2]=-normal[2];
		}
		double t;
		NORMALIZE(normal,t);	

/*
		printf("\n r2 = %f\n",r2);				
		printf("\ny %f %f %f\n",y[0],y[1],y[2]);
		printf("\nc %f %f %f r2 %f\n",cc[0],cc[1],cc[2],r2);
		printf("\nproj %f %f %f \n",proj[0],proj[1],proj[2]);
		printf("\ndist %f r %f\n",dist,r);
	    getchar();
*/		
		return;
	}

	// generic quadric projection

	double b[3];
	double poly[7];

	/*double c = Q[3][3];
	b[0] = 2*Q[3][0];
	b[1] = 2*Q[3][1];
	b[2] = 2*Q[3][2];

	// compute sixth degree polynomial
	double a11 = Q[0][0];
	double a12 = Q[0][1];
	double a13 = Q[0][2];
	double a22 = Q[1][1];
	double a33 = Q[2][2];
	double a23 = Q[1][2];*/

	double c = Q[3];
	b[0] = 2*Q[9];
	b[1] = 2*Q[8];
	b[2] = 2*Q[6];

	// compute sixth degree polynomial
	double a11 = Q[0];
	double a12 = Q[4];
	double a13 = Q[7];
	double a22 = Q[1];
	double a33 = Q[2];
	double a23 = Q[5];

	double a112 = a11*a11;
	double a222 = a22*a22;
	double a332 = a33*a33;

	double a232 = a23*a23;
	double a132 = a13*a13;
	double a122 = a12*a12;

	double b02 = b[0]*b[0];
	double b12 = b[1]*b[1];
	double b22 = b[2]*b[2];

	double y02 = y[0]*y[0];
	double y12 = y[1]*y[1];
	double y22 = y[2]*y[2];

	double a234 = a232*a232;
	double a134 = a132*a132;
	double a123 = a12*a12*a12;
	double a124 = a123*a12;

	double a233 = a232*a23;
	double a133 = a132*a13;

	// t^6 coefficient
	poly[6] = 16* (a11*a22*a33-a122*a33-a11*a232+2*a12*a13*a23-a132*a22)*
					(4*a11*a22*a33*c-4*a122*a33*c-4*a11*a232*c+8*a12*a13*a23*c
					-4*a132*a22*c-b02*a22*a33+2*b[0]*
					b[1]*a12*a33-b12*a11*a33+b02*a232-2*b[0]*b[1]*a13*a23-
					2*b[0]*b[2]*a12*a23+2*b[1]*b[2]*a11*a23+2*b[0]*b[2]*a13*a22
					-b22*a11*a22+b12*
					a132-2*b[1]*b[2]*a12*a13+b22*a122);   

	// t^5 coefficient
	poly[5] = 16*(a22*a33+a11*a33-a232+a11*a22-a132-a122)*
					(4*a11*a22*a33*c-4*a122*a33*c-4*a11*a232*c+8*a12*a13*a23*c
					-4*a132*a22*c-b02*a22*a33+2*b[0]*b[1]*a12*a33-
					b12*a11*a33+b02*a232-2*b[0]*b[1]*a13*a23-
					2*b[0]*b[2]*a12*a23+2*b[1]*b[2]*a11*a23+2*b[0]*b[2]*a13*a22
					-b22*a11*a22+b12*a132-2*b[1]*
					b[2]*a12*a13+b22*a122);

	// t^4 coefficient
	poly[4] = +(((16*a222+64*a11*a22-32*a122+16*a112)*a332+
				((-32*a22-64*a11)*a232+64*a12*a13*a23+64*a11*a222+
				(-64*a132-64*a122+64*a112)*a22-32*a11*a132-64*a11*a122)*a33+16*a234
				+(-64*a11*a22+32*a132+32*a122-32*a112)*a232+
				(64*a12*a13*a22+64*a11*a12*a13)*a23+(16*a112-32*a132)*a222
				+(-64*a11*a132-32*a11*a122)*a22+16*a134+32*a122*a132+16*a124)*c+((16*y02*a11+16*b[0]*y[0])*a222+
				(-16*y02*a122+(-32*y[0]*y[1]*a11-16*b[0]*y[1]-16*y[0]*b[1])*
				a12+16*y12*a112-4*b12-16*b02)*a22+32*y[0]*y[1]*a123+(-16*y12*a11+16*b[1]*y[1]+16*b[0]*y[0])*a122+
				((-16*b[0]*y[1]-16*y[0]*b[1])*a11+24*b[0]*b[1])*a12+16*b[1]*y[1]*a112+(-16*b12-4*b02)*a11)*a332+(
				((-32*y02*a11-32*b[0]*y[0])*a22+16*y02*a122+(32*y[0]*y[1]*a11+16*b[0]*y[1]+
				16*y[0]*b[1])*a12-16*y12*a112+4*b12+16*b02)*a232+(
				((32*y02*a12+32*y[0]*y[1]*a11+16*b[0]*y[1]+16*y[0]*b[1])*a13+
				(32*y[0]*y[2]*a11+16*b[0]*y[2]+16*y[0]*b[2])*a12-32*y[1]*y[2]*a112+8*b[1]*b[2])*a22+
				(-96*y[0]*y[1]*a122+(32*y12*a11-32*b[1]*y[1]-32*b[0]*y[0])*a12+(16*b[0]*y[1]+
				16*y[0]*b[1])*a11-24*b[0]*b[1])*a13-32*y[0]*y[2]*a123+
				(32*y[1]*y[2]*a11-16*b[1]*y[2]-16*y[1]*b[2])*a122+((16*b[0]*y[2]+16*y[0]*b[2])*a11-
				24*b[0]*b[2])*a12+(-16*b[1]*y[2]-16*y[1]*b[2])*a112+32*b[1]*b[2]*a11)*a23+
				(-16*y02*a132+(-32*y[0]*y[2]*a11-16*b[0]*y[2]-16*y[0]*b[2])*a13+16*y22*a112
				-4*b22-16*b02)*a222+((32*y[0]*y[1]*a12-32*y12*a11)*a132+
				(32*y[0]*y[2]*a122+(32*y[1]*y[2]*a11+16*b[1]*y[2]+16*y[1]*b[2])*a12+32*b[0]*b[2])*a13-
				32*y22*a11*a122+32*b[0]*b[1]*a12+(-16*b22-16*b12-16*b02)*a11)*a22+
				(16*y12*a122+(16*b[0]*y[1]+16*y[0]*b[1])*a12-32*b[1]*y[1]*a11+16*b12+4*b02)*a132+
				(-32*y[1]*y[2]*a123+(-16*b[0]*y[2]-16*y[0]*b[2])*a122+((16*b[1]*y[2]+
				16*y[1]*b[2])*a11-24*b[1]*b[2])*a12+8*b[0]*b[2]*a11)*a13+16*y22*a124+8*b22*a122+32*b[0]*b[1]*a11* 
				a12+(-4*b22-16*b12)*a112)*a33+(16*y02*a11+16*b[0]*y[0])*a234+
				((-32*y02*a12-32*y[0]*y[1]*a11-16*b[0]*y[1]-16*y[0]*b[1])*a13+(-32*y[0]*y[2]*a11-
				16*b[0]*y[2]-16*y[0]*b[2])*a12+32*y[1]*y[2]*a112-8*b[1]*b[2])*a233+(
				(16*y02*a132+(32*y[0]*y[2]*a11+16*b[0]*y[2]+16*y[0]*b[2])*a13-16*y22*a112+4*b22+
				16*b02)*a22+(64*y[0]*y[1]*a12+16*y12*a11+16*b[1]*y[1]+16*b[0]*y[0])*a132+
				(64*y[0]*y[2]*a122+(-96*y[1]*y[2]*a11+16*b[1]*y[2]+16*y[1]*b[2])*a12+(-16*b[0]*y[2]-
				16*y[0]*b[2])*a11-8*b[0]*b[2])*a13+(16*y22*a11+16*b[2]*y[2]+16*b[0]*y[0])*a122+
				((-16*b[0]*y[1]-16*y[0]*b[1])*a11-8*b[0]*b[1])*a12+(16*b[2]*y[2]+16*b[1]*y[1])*a112+
				8*b02*a11)*a232+((-32*y[0]*y[1]*a133+
				(-96*y[0]*y[2]*a12+32*y[1]*y[2]*a11-16*b[1]*y[2]-16*y[1]*b[2])*a132+(
				(32*y22*a11-32*b[2]*y[2]-32*b[0]*y[0])*a12+(16*b[0]*y[1]+16*y[0]*b[1])*a11-24*b[0]*b[1])*a13+
				((16*b[0]*y[2]+16*y[0]*b[2])*a11-24*b[0]*b[2])*a12+(-16*b[1]*y[2]-16*y[1]*b[2])*a112+32*b[1]*b[2]*a11)*a22
				+(-32*y12*a12-16*b[0]*y[1]-16*y[0]*b[1])*a133+
				(64*y[1]*y[2]*a122+(16*b[0]*y[2]+16*y[0]*b[2])*a12+(16*b[1]*y[2]+16*y[1]*b[2])*a11-8*b[1]*b[2])*a132+
				(-32*y22*a123+(16*b[0]*y[1]+16*y[0]*b[1])*a122+((-32*b[2]*y[2]-32*b[1]*y[1])*a11+8*b22+8*b12+
				8*b02)*a12-24*b[0]*b[1]*a11)*a13+(-16*b[0]*y[2]-16*y[0]*b[2])*a123+
				((16*b[1]*y[2]+16*y[1]*b[2])*a11-8*b[1]*b[2])*a122-24*b[0]*b[2]*a11*a12+24*b[1]*b[2]*a112)*a23+
				(32*y[0]*y[2]*a133+(-16*y22*a11+16*b[2]*y[2]+16*b[0]*y[0])*a132+((-16*b[0]*y[2]-
				16*y[0]*b[2])*a11+24*b[0]*b[2])*a13+16*b[2]*y[2]*a112+(-16*b22-4*b02)*a11)*a222+(16*
				y12*a134-32*y[1]*y[2]*a12*a133+(16*y22*a122+(-16*b[0]*y[1]-16*y[0]*b[1])*a12+8*b12)*a132+
				((16*b[0]*y[2]+16*y[0]*b[2])*a122+((16*b[1]*y[2]+16*y[1]*b[2])*a11-24*b[1]*b[2])*a12+32*b[0]*b[2]*a11)*a13
				+(-32*b[2]*y[2]*a11+16*b22+4*b02)*a122+8*b[0]*b[1]*a11*a12+
				(-16*b22-4*b12)*a112)*a22+16*b[1]*y[1]*a134+((-16*b[1]*y[2]-16*y[1]*b[2])*a12-8*b[0]*b[2])*a133
				+((16*b[2]*y[2]+16*b[1]*y[1])*a122-8*b[0]*b[1]*a12+(4*b22+16*b12)*a11)*
				a132+((-16*b[1]*y[2]-16*y[1]*b[2])*a123-8*b[0]*b[2]*a122-24*b[1]*b[2]*a11*a12)*a13+16*b[2]*y[2]*a124
				-8*b[0]*b[1]*a123+(16*b22+4*b12)*a11*a122);

	// t^3 coefficient
	poly[3] = +(((16*a22+16*a11)*
				a332+(-16*a232+16*a222+64*a11*a22-16*a132-32*a122+16*a112)*a33
				+(-16*a22-32*a11)*a232+32*a12*a13*a23+16*a11*a222+(-32*a132-16*a122+16*a112)*a22-16*a11*a132-16*a11*
				a122)*c+(((16*y12+16*y02)*a11+16*b[0]*y[0])*a22+(-16*y12-16*y02)*a122+(-16*b[0]*y[1]
				-16*y[0]*b[1])*a12+16*b[1]*y[1]*a11-4*b12-4*b02)*a332+(
				((-16*y12-16*y02)*a11-16*b[0]*y[0])*a232+(-32*y[1]*y[2]*a11*a22+
				((32*y12+32*y02)*a12+16*b[0]*y[1]+16*y[0]*b[1])*a13+32*y[1]*y[2]*a122+(16*b[0]*y[2]+16*y[0]*b[2])*a12
				+(-16*b[1]*y[2]-16*y[1]*b[2])*a11+8*b[1]*b[2])*a23+((16*y22+16*y02)*a11+16*b[0]*y[0])*a222+
				((-16*y12-16*y02)*a132+(-32*y[0]*y[2]*a11-16*b[0]*y[2]-16*y[0]*b[2])*a13+
				(-16*y22-16*y02)*a122+(-32*y[0]*y[1]*a11-16*b[0]*y[1]-16*y[0]*b[1])*a12+(16*y22+16*y12)*a112
				-4*b22-4*b12-16*b02)*a22-16*b[1]*y[1]*a132+
				(32*y[0]*y[2]*a122+(16*b[1]*y[2]+16*y[1]*b[2])*a12+8*b[0]*b[2])*a13+32*y[0]*y[1]*a123+
				((-16*y22-16*y12)*a11+16*b[1]*y[1]+16*b[0]*y[0])*a122+
				((-16*b[0]*y[1]-16*y[0]*b[1])*a11+24*b[0]*b[1])*a12+16*b[1]*y[1]*a112+(-4*b22-
				16*b12-4*b02)*a11)*a33+32*y[1]*y[2]*a11*a233+
				(((-16*y22-16*y02)*a11-16*b[0]*y[0])*a22+(32*y[0]*y[2]*a11-64*y[1]*y[2]*a12)*a13+32*y[0]*y[1]*a11*a12
				+(-16*y22-16*y12)*a112+(16*b[2]*y[2]+16*b[1]*y[1])*a11+8*b02)* 
				a232+((32*y[1]*y[2]*a132+((32*y22+32*y02)*a12+16*b[0]*y[1]+16*y[0]*b[1])*a13+(16*b[0]*y[2]
				+16*y[0]*b[2])*a12+(-16*b[1]*y[2]-16*y[1]*b[2])*a11+8*b[1]*b[2])*a22-64*y[0]* 
				y[2]*a12*a132+(-64*y[0]*y[1]*a122+((32*y22+32*y12)*a11-32*b[2]*y[2]-32*b[1]*y[1]
				-32*b[0]*y[0])*a12+(16*b[0]*y[1]+16*y[0]*b[1])*a11-16*b[0]*b[1])*a13+
				((16*b[0]*y[2]+16*y[0]*b[2])*a11-16*b[0]*b[2])*a12+(-16*b[1]*y[2]-16*y[1]*b[2])*a112+24*b[1]*b[2]*a11)*a23+
				((-16*y22-16*y02)*a132+(-16*b[0]*y[2]-16*y[0]*b[2])*a13+16*b[2]*y[2]*a11-4*b22-4*b02)*a222+(32*y[0]*y[2]*a133+
				(32*y[0]*y[1]*a12+(-16*y22-16*y12)*a11+16*b[2]*y[2]+16*b[0]*y[0])*a132+((16*b[1]*y[2]+16*y[1]*b[2])*a12
				+(-16*b[0]*y[2]-16*y[0]*b[2])*a11+24*b[0]*b[2])*a13-16*b[2]*y[2]*
				a122+8*b[0]*b[1]*a12+16*b[2]*y[2]*a112+(-16*b22-4*b12-4*b02)*a11)*a22+(8*b12-16*b[1]*y[1]*a11)*a132+
				(((16*b[1]*y[2]+16*y[1]*b[2])*a11-16*b[1]*b[2])*a12+8*b[0]*b[2]*a11)*a13+(8*b22-16*b[2]*y[2]*a11)*a122
				+8*b[0]*b[1]*a11*a12+(-4*b22-4*b12)*a112);

	// t^2 coefficient
	poly[2] =	+((4*a332+(16*a22+16*a11)*a33-8*a232+4*a222+16*a11*a22-8*a132-8*a122+4*a112)*c
				+(4*y12*a22+8*y[0]*y[1]*a12+4*y02*a11+4*b[1]*y[1]+4*b[0]*y[0])*a332+(-4*y12*a232+
				(-8*y[1]*y[2]*a22-8*y[0]*y[1]*a13-8*y[0]*y[2]*a12-4*b[1]*y[2]-4*y[1]*b[2])*a23+4*y22*a222+
				((16*y22+16*y12+16*y02)*a11+16*b[0]*y[0])*a22-4*y02*a132+
				(-8*y[1]*y[2]*a12-8*y[0]*y[2]*a11-4*b[0]*y[2]-4*y[0]*b[2])*a13+(-8*y22-16*y12-16*y02)*a122
				+(-16*b[0]*y[1]-16*y[0]*b[1])*a12+4*y22*a112+16*b[1]*y[1]*a11-b22-4*b12 
				-4*b02)*a33+8*y[1]*y[2]*a233+(-4*y22*a22+8*y[0]*y[2]*a13+8*y[0]*y[1]*a12+(-16*y22-16*y12-
				8*y02)*a11+4*b[2]*y[2]+4*b[1]*y[1]-8*b[0]*y[0])*a232+(
				(-8*y[0]*y[1]*a13-8*y[0]*y[2]*a12-4*b[1]*y[2]-4*y[1]*b[2])*a22+8*y[1]*y[2]*a132+
				((24*y22+24*y12+24*y02)*a12-8*y[0]*y[1]*a11+12*b[0]*y[1]+12*y[0]*b[1])*a13+8*y[1]*
				y[2]*a122+(-8*y[0]*y[2]*a11+12*b[0]*y[2]+12*y[0]*b[2])*a12+8*y[1]*y[2]*a112
				+(-16*b[1]*y[2]-16*y[1]*b[2])*a11+6*b[1]*b[2])*a23+
				(8*y[0]*y[2]*a13+4*y02*a11+4*b[2]*y[2]+4*b[0]*y[0])*a222+((-16*y22-8*y12-16*y02)*a132
				+(-8*y[1]*y[2]*a12-16*b[0]*y[2]-16*y[0]*b[2])*a13-4*y02*a122+
				(-8*y[0]*y[1]*a11-4*b[0]*y[1]-4*y[0]*b[1])*a12+4*y12*a112+16*b[2]*y[2]*a11-4*b22-b12-4*b02)*a22+8*y[0]*y[2]*a133+
				(8*y[0]*y[1]*a12-4*y22*a11+4*b[2]*y[2]-8*b[1]*y[1]+4*b[0]*y[0])*a132+
				(8*y[0]*y[2]*a122+(-8*y[1]*y[2]*a11+12*b[1]*y[2]+12*y[1]*b[2])*a12+(-4*b[0]*y[2]
				-4*y[0]*b[2])*a11+6*b[0]*b[2])*a13+8*y[0]*y[1]*a123+
				(-4*y12*a11-8*b[2]*y[2]+4*b[1]*y[1]+4*b[0]*y[0])*a122+((-4*b[0]*y[1]-4*y[0]*b[1])*a11+6*b[0]*b[1])*a12+
				(4*b[2]*y[2]+4*b[1]*y[1])*a112+(-4*b22-4*b12-b02)*a11);

	// t coefficient
	poly[1] =((4*a33+4*a22+4*a11)*c+((4*y22+4*y12)*a22+8*y[0]*y[1]*a12+(4*y22+4*y02)*a11+
				4*b[1]*y[1]+4*b[0]*y[0])*a33+(-4*y22-4*y12)*a232+
				(-8*y[0]*y[1]*a13-8*y[0]*y[2]*a12+8*y[1]*y[2]*a11-4*b[1]*y[2]-4*y[1]*b[2])*a23
				+(8*y[0]*y[2]*a13+(4*y12+4*y02)*a11+4*b[2]*y[2]+4*b[0]*y[0])*a22+(-4*y22-4*y02)*a132+
				(-8*y[1]*y[2]*a12-4*b[0]*y[2]-4*y[0]*b[2])*a13+(-4*y12-4*y02)*a122+(-4*b[0]*y[1]
				-4*y[0]*b[1])*a12+(4*b[2]*y[2]+4*b[1]*y[1])*a11-b22-b12-b02);

	// t^0 coefficient
	poly[0] = +c+y22*a33+2*y[1]*y[2]*a23+y12*a22+2*y[0]*y[2]*a13+2*y[0]*y[1]*a12+y02*a11+b[2]*y[2]+b[1]*y[1]+b[0]*y[0];

	double roots[6];
	int num=6;
	
	if (fastProjection)
		getRealRootsSturm(poly,6,roots,num);
	else
		getRealRootsCompanion(poly,6,roots,num);

	// skip non properly projected patch (typical case of distant patch of no interest)
	if (num==0)
	{
		//printf("\n %lf %lf %lf %lf %lf %lf %lf",poly[6],poly[5],poly[4],poly[3],poly[2],poly[1],poly[0]);
		dist = INFINITY;
		//getchar();
		return;
	}

	double projt[6][3];
	int opt=7;
	double minDist = INFINITY;

	// get solution for each root and choose the nearest point
	for (int i=0;i<num;i++)
	{				
		double jj = y[0]-roots[i]*b[0];
		double kk = y[1]-roots[i]*b[1];
		double ll = y[2]-roots[i]*b[2];		

		/*double aa=Q[0][0]*2*roots[i]+1;
		double ee=Q[1][1]*2*roots[i]+1;
		double ii=Q[2][2]*2*roots[i]+1;
		
		double dd = Q[0][1]*2*roots[i];
		double gg = Q[0][2]*2*roots[i];
		double ff = Q[1][2]*2*roots[i];		*/

		// solve a 3x3 linear system by Cramer rule
		// perform all inline and exploit symmetry		
		double aa=Q[0]*2*roots[i]+1;
		double ee=Q[1]*2*roots[i]+1;
		double ii=Q[2]*2*roots[i]+1;
		
		double dd = Q[4]*2*roots[i];
		double gg = Q[7]*2*roots[i];
		double ff = Q[5]*2*roots[i];

		double ff2 = ff*ff;
		double fg = ff*gg;

		double den = aa*ee*ii+2*dd*fg-gg*gg*ee-dd*dd*ii-aa*ff2;
		double num1 = jj*ee*ii+dd*ff*ll+kk*fg-gg*ee*ll-dd*kk*ii-jj*ff2;
		double num2 = aa*kk*ii+jj*fg+gg*dd*ll-gg*gg*kk-jj*dd*ii-aa*ff*ll;
		double num3 = aa*ee*ll+dd*kk*gg+jj*dd*ff-jj*ee*gg-dd*dd*ll-aa*kk*ff;
		
		projt[i][0]=num1/den;
		projt[i][1]=num2/den;
		projt[i][2]=num3/den;
		
		if (numeric_limits<double>::infinity() == projt[i][0])
			projt[i][0]=num1/(den+1e-20);
		
		if (numeric_limits<double>::infinity() == projt[i][1])
			projt[i][1]=num2/(den+1e-20);

		if (numeric_limits<double>::infinity() == projt[i][2])
			projt[i][2]=num3/(den+1e-20);
		
		//printf("\n %f %f %f\n",projt[i][0],projt[i][1],projt[i][2]);
		double a = (projt[i][0]-y[0]);
		double b = (projt[i][1]-y[1]);
		double c = (projt[i][2]-y[2]);
		double distl = a*a+b*b+c*c;
		if (distl<minDist)
		{
			minDist = distl;
			opt = i;
		}
	}
	
	proj[0]=projt[opt][0];
	proj[1]=projt[opt][1];
	proj[2]=projt[opt][2];
	dist = sqrt(minDist);
	//printf("\nopt %f %f %f\n",proj[0],proj[1],proj[2]);

	// get the normal depending on the insidness value
	double insideness=0;
	double yi,yj;
	for (int i=0;i<4;i++)
	{
		if (i==3) 
			yi = 1;
		else
			yi = y[i];

		for (int j=0;j<4;j++)
		{						
			if (j==3)
				yj = 1;
			else
				yj = y[j];

			double qij = 0;
			
			if (i==0 && j==0)
				qij = Q[0];
			else if (i==1 && j==1)
				qij = Q[1];
			else if (i==2 && j==2)
				qij = Q[2];
			else if (i==3 && j==3)
				qij = Q[3];
			if ((i==1 && j==0) || (i==0 && j==1))
				qij = Q[4];
			else if ((i==2 && j==0) || (i==0 && j==2))
				qij = Q[7];
			else if ((i==3 && j==0) || (i==0 && j==3))
				qij = Q[9];
			else if ((i==1 && j==2) || (i==2 && j==1))
				qij = Q[5];
			else if ((i==2 && j==3) || (i==3 && j==2))
				qij = Q[6];
			else if ((i==1 && j==3) || (i==3 && j==1))
				qij = Q[8];
			
			//insideness+=yi*yj*Q[i][j];
			insideness+=yi*yj*qij;
		}
	}

	// check the correctness of the projection
	double test = 0;
	for (int i=0;i<4;i++)
	{
		if (i==3) 
			yi = 1;
		else
			yi = proj[i];

		for (int j=0;j<4;j++)
		{						
			if (j==3)
				yj = 1;
			else
				yj = proj[j];

			double qij = 0;
			
			if (i==0 && j==0)
				qij = Q[0];
			else if (i==1 && j==1)
				qij = Q[1];
			else if (i==2 && j==2)
				qij = Q[2];
			else if (i==3 && j==3)
				qij = Q[3];
			if ((i==1 && j==0) || (i==0 && j==1))
				qij = Q[4];
			else if ((i==2 && j==0) || (i==0 && j==2))
				qij = Q[7];
			else if ((i==3 && j==0) || (i==0 && j==3))
				qij = Q[9];
			else if ((i==1 && j==2) || (i==2 && j==1))
				qij = Q[5];
			else if ((i==2 && j==3) || (i==3 && j==2))
				qij = Q[6];
			else if ((i==1 && j==3) || (i==3 && j==1))
				qij = Q[8];

			//test+=yi*yj*Q[i][j];
			test+=yi*yj*qij;
		}
	}	
	if (fabs(test)>1e-3)
	{
		//Projected point does not belong to quadric
		//Unreliable projection detected is discarded (probably a distant patch)
		dist = INFINITY;		
		return;
	}

	// concave sphere
	if (type==3)
		insideness = -insideness;

	SUB(normal,y,proj);
	if (insideness<0)
	{
		normal[0]=-normal[0];
		normal[1]=-normal[1];
		normal[2]=-normal[2];
	}
	double t;
	NORMALIZE(normal,t);	
}

