
#include "DelphiShared.h"
#include <cstdlib>
#include <sstream>

void DelPhiShared::init()
{
	epsmap = NULL;
	idebmap = NULL;
	buildEpsMap = false;
	buildStatus = false;
	x = NULL;
	y = NULL;
	z = NULL;
	// atoms = NULL;
	scspos = NULL;
	scsnor = NULL;
	scsarea = NULL;
	ibgp = NULL;
	status = NULL;
	bilevel_status = NULL;
	cavitiesVec = NULL;
	delphiBinding = false;
	multi_diel = false;
	tempStatus = NULL;
	bilevel_tempStatus = NULL;
	isAvailableAtomInfo = false;
	atsurf = NULL;
	maxNumAtoms = -1;
	domainShrinkage = 0.;
	optimizeGrids = true;
}


DelPhiShared::DelPhiShared()
{
	init();
}


// SD: PB_NEW
// new constructor for API
DelPhiShared::DelPhiShared(double scale, double *coord_min, double *coord_max,
						   const std::vector<Atom> &ats, int max_num_atoms, double domain_shrinkage, bool optimize_grids, bool eps_flag, bool stat_flag, bool multi,double tol)
{
	init();
	init(scale,coord_min,coord_max,ats,max_num_atoms,domain_shrinkage,optimize_grids,eps_flag,stat_flag,multi,tol);
}


// SD: PB_NEW
// new constructor for API
void DelPhiShared::init(double scale, double *coord_min, double *coord_max,
						const std::vector<Atom> &ats, int max_num_atoms, double domain_shrinkage, bool optimize_grids,
						bool eps_flag, bool stat_flag, bool multi,double tol)
{
	cout << endl << INFO << "Starting grid new PB initialization";
	atoms = ats;
	buildEpsMap = eps_flag;
	buildStatus = stat_flag;
	multi_diel = multi;
	isAvailableAtomInfo = false;
	maxNumAtoms = max_num_atoms;
	domainShrinkage = domain_shrinkage;
	optimizeGrids = optimize_grids;

	bool flag = buildGrid(scale,coord_min,coord_max,tol);
	if (!flag)
	{
		cout << endl << ERR << "DelPhiShared object initialization failed";
		exit(-1);
	}
	cout << endl << INFO << "Initialization completed";
}


void DelPhiShared::init(int max_num_atoms, double domain_shrinkage, bool optimize_grids,
						double scale, double perfill, string fn, bool eps_flag,
						bool stat_flag, bool multi, bool atinfo)
{
	buildEpsMap = eps_flag;
	buildStatus = stat_flag;
	multi_diel = multi;
	isAvailableAtomInfo = atinfo;
	maxNumAtoms = max_num_atoms;
	domainShrinkage = domain_shrinkage;
	optimizeGrids = optimize_grids;
	
	cout << endl << INFO << "Loading atoms....";

	bool flag = loadAtoms(fn);

	if (!flag)
	{
		cout << endl << ERR << "Missing or corrupt atoms file. Initialization failed";
		exit(-1);
	}
	flag = buildGrid(scale, perfill);

	if (!flag)
	{
		cout << endl << ERR << "Initialization failed";
		exit(-1);
	}
	cout << endl << INFO << "Initialization completed";
}


void DelPhiShared::init(double scale, double perfill, const std::vector<Atom> &ats,
						int max_num_atoms, double domain_shrinkage, bool optimize_grids,
						bool eps_flag, bool stat_flag, bool multi, bool atinfo)
{
	atoms = ats;
	buildEpsMap = eps_flag;
	buildStatus = stat_flag;
	multi_diel = multi;
	isAvailableAtomInfo = atinfo;
	maxNumAtoms = max_num_atoms;
	domainShrinkage = domain_shrinkage;
	optimizeGrids = optimize_grids;

	bool flag = buildGrid(scale,perfill);

	if (!flag)
	{
		cout << endl << ERR << "DelPhiShared object initialization failed";
		exit(-1);
	}
}


DelPhiShared::DelPhiShared(int max_num_atoms, double domain_shrinkage, bool optimize_grids,
						   double scale, double perfill, string fn,
						   bool eps_flag, bool stat_flag, bool multi, bool atinfo)
{
	init();
	init(max_num_atoms,domain_shrinkage,optimize_grids,scale,perfill,fn,eps_flag,stat_flag,multi,atinfo);
}


DelPhiShared::DelPhiShared(double scale, double perfill, const std::vector<Atom> &ats,
						   int max_num_atoms, double domain_shrinkage, bool optimize_grids,
						   bool eps_flag, bool stat_flag, bool multi)
{
	init();

    bool atInfo = false;

    init(scale,perfill,ats,max_num_atoms,domain_shrinkage,optimize_grids,eps_flag,stat_flag,multi,atInfo);
}


void DelPhiShared::DelPhiBinding(	double xmin,double ymin,double zmin,
									double xmax,double ymax,double zmax,
									double c1,double c2,double c3,double rmax,double perf,
									int *local_i_epsmap,int igrid,double scale,
									double *local_i_scspos,double *local_i_scsnor,
									bool *local_i_idebmap,int *local_i_ibgp,int maxbgp,
									bool bstatus,int *_atsurf)
{
	cout << endl << INFO << "DelPhi Binding...";
	delphiBinding = true;
	buildStatus = bstatus;
	atsurf = _atsurf;
	
	if (x != NULL)
		deleteVector<double>(x);
	if (y != NULL)
		deleteVector<double>(y);
	if (z != NULL)
		deleteVector<double>(z);

	if (buildStatus)
	{
		if (!optimizeGrids)
		{
			if (status != NULL)
				deleteVector<int>(status);
		}
		else
		{
			if (bilevel_status != NULL)
			{
				deleteBilevelGridCells<int>(bilevel_status, nx, ny, nz);
				deleteVector<int *>(bilevel_status);
			}
		}
	}

	this->xmin = xmin;
	this->ymin = ymin;
	this->zmin = zmin;

	this->rmaxdim = rmax;
	this->perfill = perf;
	this->xmax = xmax;
	this->ymax = ymax;
	this->zmax = zmax;

	this->baricenter[0] = c1;
	this->baricenter[1] = c2;
	this->baricenter[2] = c3;

	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);

	this->scale = scale;
	side = 1./scale;
	hside = 0.5/scale;
	A = side*side;
	nx = igrid;
	ny = igrid;
	nz = igrid;

	/*
	if (atoms != NULL)
	{
		for (int i=0; i<numAtoms; i++)
			delete atoms[i];
		delete[] atoms;
	}
	atoms = NULL;
	*/
    if (atoms.size() != 0)
        atoms.clear();
	
	// assuming all maps are already allocated before binding 
	this->epsmap  = local_i_epsmap;
	this->idebmap = local_i_idebmap;
	this->scsnor = local_i_scsnor;
	this->scsarea = NULL;
	this->scspos = local_i_scspos;
	this->ibgp = local_i_ibgp;

	for (int i=0; i<nx; i++)
		x[i] = xmin + i*side;

	for (int i=0; i<ny; i++)
		y[i] = ymin + i*side;

	for (int i=0; i<nz; i++)
		z[i] = zmin + i*side;

	// allocate grid point status memory
	// and initialize all to temporary outside status
	if (buildStatus)
	{
		if (!optimizeGrids)
		{
			/*
			status = allocateMatrix3D<int>(nx,ny,nz);
			for (int k=0; k<nz; k++)
				for (int j=0; j<ny; j++)
					for (int i=0; i<nx; i++)
						status[k][j][i] = 2;
			*/
			int64_t tot = nx*ny*nz;;
			status = allocateVector<int>(tot);

			if (status == NULL)
			{
				cout << endl << ERR << "Not enough memory to allocate status map";
				cout << endl;
				exit(-1);
			}
			for (int64_t i=0; i<tot; i++)
				status[i] = STATUS_POINT_TEMPORARY_OUT;
		}
		else
		{
			bilevel_status = allocateBilevelGridCells<int>(nx, ny, nz);
		}
	}
	// set max bgp
	this->maxbgp = maxbgp;
	
	cout << "ok!";
	cout << endl << INFO << "Max number of bgps was set from DelPhi to " << maxbgp;
}


void DelPhiShared::finalizeBinding(int *ibnum)
{
	// update indexes to be conformant to Fortran
	for (int i=0; i<nbgp; i++)
	{
	  ibgp[3*i]++;
	  ibgp[3*i+1]++;
	  ibgp[3*i+2]++;
	}

	// update indexes to be conformant to Fortran
	if (atsurf != NULL)
		for (int i=0; i<nbgp; i++)
			atsurf[i]++;

	epsmap = NULL; //avoid deletion
	scspos = NULL; //avoid deletion
	scsnor = NULL; //avoid deletion
	idebmap = NULL; //avoid deletion
	//scsarea = NULL; //for now not used, so it must be destroied
	ibgp = NULL;   //avoid deletion
	// save ibnum
	(*ibnum) = nbgp;
}


bool DelPhiShared::loadAtoms(string fn)
{
	// load atoms file
	ifstream fin;
	fin.open(fn.c_str(), ios::in);

	char buffer[BUFLEN];
	double max_rad = 0;


	if (multi_diel)
		cout << endl << INFO << "Atom dielectric info is available..";

	if (isAvailableAtomInfo)
		cout << endl << INFO << "Atom Info is available..";

	if (fin.fail())
	{
		cout << endl << WARN << "Cannot read file " << fn;
		return false;
	}

	// vector<Atom*> list;
	// Atom *atm;
	if (atoms.size() != 0)
		atoms.clear();

	bool is_pqr_file = false;

	string molFileName = fn.c_str();

	if (molFileName.find(".pqr") != string::npos)
	{
		is_pqr_file = true;
		isAvailableAtomInfo = true;
	}

	auto parse_int_field = [](const string &field, int &value) -> bool
	{
		char *end_ptr = NULL;
		long parsed = strtol(field.c_str(), &end_ptr, 10);
		if (end_ptr == field.c_str() || *end_ptr != '\0')
			return false;
		value = (int)parsed;
		return true;
	};

	auto parse_double_field = [](const string &field, double &value) -> bool
	{
		char *end_ptr = NULL;
		value = strtod(field.c_str(), &end_ptr);
		if (end_ptr == field.c_str() || *end_ptr != '\0')
			return false;
		return true;
	};

	auto parse_pqr_line = [&](const char *line, bool &is_atom_record, int &serial,
							 string &atom_name, string &res_name, string &chain_name, int &res_id,
							 double &x_coord, double &y_coord, double &z_coord,
							 double &charge, double &radius) -> bool
	{
		is_atom_record = false;
		serial = -1;
		vector<string> fields;
		{
			istringstream line_stream(line);
			string field;
			while (line_stream >> field)
				fields.push_back(field);
		}

		if (fields.size() == 0)
			return true;

		if (fields[0] != "ATOM" && fields[0] != "HETATM")
			return true;

		is_atom_record = true;

		if (fields.size() == 11)
		{
			atom_name = fields[2];
			res_name = fields[3];
			chain_name = fields[4];
			return parse_int_field(fields[1], serial) &&
				   parse_int_field(fields[5], res_id) &&
				   parse_double_field(fields[6], x_coord) &&
				   parse_double_field(fields[7], y_coord) &&
				   parse_double_field(fields[8], z_coord) &&
				   parse_double_field(fields[9], charge) &&
				   parse_double_field(fields[10], radius);
		}

		if (fields.size() == 10)
		{
			atom_name = fields[2];
			res_name = fields[3];
			chain_name = "X";
			return parse_int_field(fields[1], serial) &&
				   parse_int_field(fields[4], res_id) &&
				   parse_double_field(fields[5], x_coord) &&
				   parse_double_field(fields[6], y_coord) &&
				   parse_double_field(fields[7], z_coord) &&
				   parse_double_field(fields[8], charge) &&
				   parse_double_field(fields[9], radius);
		}

		return false;
	};

	double max_x=-1e6, min_x=1e6, max_y=-1e6, min_y=1e6, max_z=-1e6, min_z=1e6;

	if (domainShrinkage > 0. && domainShrinkage < 1.)
	{
		double x, y, z, radius;

		// Read the file just to find the atomic system extrema that are then used to avoid
		// storing the atoms which are not positioned within the chosen shrinked domain

		while (!fin.eof( ))
		{
			fin.getline(buffer,BUFLEN);

			if (strlen(buffer) <= 1)
				// skip empty line
				continue;

				if (is_pqr_file)
				{
					string atomName, resName, chain;
					int serial = -1;
					int resId = 0;
					double charge = 0.0;
					bool is_atom_record = false;
					if (!parse_pqr_line(buffer, is_atom_record, serial, atomName, resName, chain,
										resId, x, y, z, charge, radius))
					{
						cout << endl << ERR << "Cannot parse pqr line: " << buffer;
						cout << endl;
						exit(-1);
					}
					if (!is_atom_record)
						continue;

					if (radius < 1e-20)
						continue;
            }
            else
            {
                if (!multi_diel && !isAvailableAtomInfo)
				{
					sscanf(buffer, "%lf %lf %lf %lf", &x, &y, &z, &radius);
				}
				else if (multi_diel)
				{
					int d = -1;

					sscanf(buffer, "%lf %lf %lf %lf %d", &x, &y, &z, &radius, &d);

					if (d == -1)
					{
						cout << endl << ERR << "Cannot get the dielectric value of the atom number " << atoms.size() << "; please add it after radius entry to the xyzr input file";
						cout << endl;
						exit(-1);
					}
				}
				else if (isAvailableAtomInfo)
				{
					int d = -1;
					char name[10],resName[10],chain[2];
					int resid;

					sscanf(buffer, "%lf %lf %lf %lf %d %s %s %d %s",
						   &x, &y, &z, &radius, &d, name, resName, &resid, chain);

					if (d == -1)
					{
						cout << endl << ERR << "Cannot get the dielectric value of the atom number " << atoms.size() << "; please add it after radius entry to the xyzr input file";
						cout << endl;
						exit(-1);
					}
				}
			}
			max_x = max(max_x, x);
			max_y = max(max_y, y);
			max_z = max(max_z, z);

			min_x = min(min_x, x);
			min_y = min(min_y, y);
			min_z = min(min_z, z);
		}
		fin.close();

		fin.open(fn.c_str(), ios::in);

		double mid_x = 0.5 * (max_x + min_x);
		double mid_y = 0.5 * (max_y + min_y);
		double mid_z = 0.5 * (max_z + min_z);
	}

	while (!fin.eof( ))
	{
		fin.getline(buffer,BUFLEN);

		if (strlen(buffer) <= 1)
			// skip empty line
			continue;

		// atm = new Atom();
		Atom atm;

			if (is_pqr_file)
			{
				string atomName, resName, chain;
				int serial = -1;
				int resId = 0;
				double charge = 0.0;
				bool is_atom_record = false;
				if (!parse_pqr_line(buffer, is_atom_record, serial, atomName, resName, chain,
									resId, atm.pos[0], atm.pos[1], atm.pos[2], charge, atm.radius))
				{
					cout << endl << ERR << "Cannot parse pqr line: " << buffer;
					cout << endl;
					exit(-1);
				}
				if (!is_atom_record)
					continue;

				if (atm.radius < 1e-20)
					continue;
        
                atm.radius2 = atm.radius*atm.radius;
                atm.charge = charge;
				atm.ai.setName(atomName);
				atm.ai.setResName(resName);
				atm.ai.setResNum(resId);
				atm.ai.setChain(chain);
				atm.ai.setSerial(serial);
			}
		else
		{
			if (!multi_diel && !isAvailableAtomInfo)
			{
				// sscanf(buffer, "%lf %lf %lf %lf",
				// &atm->pos[0], &atm->pos[1], &atm->pos[2], &atm->radius);
				// atm->radius2 = atm->radius*atm->radius;
				sscanf(buffer, "%lf %lf %lf %lf",
					   &atm.pos[0], &atm.pos[1], &atm.pos[2], &atm.radius);
			}
			else if (multi_diel)
			{
				int d = -1;

				// sscanf(buffer, "%lf %lf %lf %lf %d", &atm->pos[0], &atm->pos[1], &atm->pos[2], &atm->radius, &d);
				// atm->radius2 = atm->radius*atm->radius;
				sscanf(buffer, "%lf %lf %lf %lf %d",
					   &atm.pos[0], &atm.pos[1], &atm.pos[2], &atm.radius, &d);
				atm.radius2 = atm.radius*atm.radius;

				if (d == -1)
				{
					cout << endl << ERR << "Cannot get the dielectric value of the atom number " << atoms.size() << "; please add it after radius entry to the xyzr input file";
					cout << endl;
					exit(-1);
				}
				else
					// atm->dielectric = d;
					atm.dielectric = d;
			}
			else if (isAvailableAtomInfo)
			{
				int d = -1;
				char name[10],resName[10],chain[2];
				int resid;

				// sscanf(buffer, "%lf %lf %lf %lf %d %s %s %d %s", &(atm->pos[0]),&(atm->pos[1]),&(atm->pos[2]),&(atm->radius),&d,name,resName,&resid,chain);
				sscanf(buffer, "%lf %lf %lf %lf %d %s %s %d %s",
						&atm.pos[0],&atm.pos[1],&atm.pos[2],&atm.radius,&d,name,resName,&resid,chain);

				// atm->radius2 = atm->radius*atm->radius;
				atm.radius2 = atm.radius*atm.radius;

				if (d == -1)
				{
					// cout << endl << ERR << "Cannot get the dielectric value of the atom number " << list.size() << "; please add it after radius entry to the xyzr input file";
					cout << endl << ERR << "Cannot get the dielectric value of the atom number " << atoms.size() << "; please add it after radius entry to the xyzr input file";
					cout << endl;
					exit(-1);
				}
				else
					// atm->dielectric = d;
					atm.dielectric = d;

				// atm->ai.setName(string(name));
				// atm->ai.setResName(string(resName));
				// atm->ai.setResNum(resid);
				// atm->ai.setChain(string(chain));

				atm.ai.setName(string(name));
				atm.ai.setResName(string(resName));
				atm.ai.setResNum(resid);
				atm.ai.setChain(string(chain));
			}
		}

		if (domainShrinkage > 0. && domainShrinkage < 1.)
		{
			double reduced_min_x, reduced_min_y, reduced_min_z;
			double reduced_max_x, reduced_max_y, reduced_max_z;

			double mid_x = 0.5 * (max_x + min_x);
			double mid_y = 0.5 * (max_y + min_y);
			double mid_z = 0.5 * (max_z + min_z);

			reduced_min_x = mid_x + domainShrinkage * (min_x - mid_x);
			reduced_min_y = mid_y + domainShrinkage * (min_y - mid_y);
			reduced_min_z = mid_z + domainShrinkage * (min_z - mid_z);

			reduced_max_x = mid_x + domainShrinkage * (max_x - mid_x);
			reduced_max_y = mid_y + domainShrinkage * (max_y - mid_y);
			reduced_max_z = mid_z + domainShrinkage * (max_z - mid_z);

			if (atm.pos[0] < reduced_min_x) continue;
			if (atm.pos[1] < reduced_min_y) continue;
			if (atm.pos[2] < reduced_min_z) continue;

			if (atm.pos[0] > reduced_max_x) continue;
			if (atm.pos[1] > reduced_max_y) continue;
			if (atm.pos[2] > reduced_max_z) continue;
		}

		// max_rad = MAX(atm->radius,max_rad);
		max_rad = MAX(atm.radius,max_rad);
		// atoms.push_back(*atm);
		// slower but cleaner
		const Atom &tt = atm;
		atoms.push_back(tt);

		if (maxNumAtoms > 0 && atoms.size() >= maxNumAtoms)
			break;

		// list.push_back(atm);
	}
	fin.close();
	
	// this is not anymore needed

	// // freeze list
	// vector<Atom*>::iterator it = list.begin();
	// if (atoms != NULL)
	// {
	// 	for (int i=0; i<numAtoms; i++)
	// 		delete atoms[i];
	// 	delete[] atoms;
	// }

	// atoms = new Atom *[list.size()];
	
	// for (int i=0; it != list.end(); it++, i++)
	// 	atoms[i] = (*it);

	// numAtoms = (int)list.size();
	// list.clear();


    unsigned numAtoms = atoms.size();
    cout << endl << INFO << "Read " << numAtoms << " atoms";

	if (numAtoms < 4)
	{
		cout << endl << ERR << "NanoShaper needs at least 4 atoms to work.";
		cout << endl << REMARK << "To emulate 4 atoms you can place dummy atoms with null radius";
		cout << endl << REMARK << "at the same centers of the real atoms";
		exit(-1);
	}

	if (max_rad == 0)
	{
		cout << endl << ERR << "All null radii? If you are using place holders atoms please set at least one radius > 0";
		exit(-1);
	}

	#if !defined(AVOID_MEM_CHECKS)
	if (!conf.parallelPocketLoop)
	{
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after atoms' loading " << current_mem_in_MB << " MB";
	}
	#endif

	return true;
}


bool DelPhiShared::loadAtoms(int na, double *x, double *r, double *q, int *d, char *atinf)
{
	cout << endl << INFO << "Read " << na << " atoms";
	/*
	if (atoms != NULL)
	{
		for (int i=0; i<numAtoms; i++)
			delete atoms[i];
		delete[] atoms;
	}
	*/
	if (atoms.size() != 0)
		atoms.clear();

	if (r == NULL || x == NULL)
	{
		cout << endl << WARN << "Atoms info not available!";
		return false;
	}

	// atoms = new Atom *[na];
	// numAtoms = na;

	int lastChar = 0;

	for (int i=0; i<na; i++)
	{
		double qq;
		int dd;
		if (q == NULL)
			qq = 0;
		else
			qq = q[i];

		if (d == NULL)
			dd = 0;
		else
			dd = d[i];

		// parse atinf variable and store if available
		if (atinf != NULL)
		{
			isAvailableAtomInfo = true;
			char buff[100];
			for (int kk=0; kk<15; kk++)
				buff[kk] = atinf[lastChar++];
			buff[15] = '\0';
			//cout << endl << buff;
			string str(buff);
			string name = str.substr(0,5);
			// remove white spaces
			name.erase( remove_if( name.begin(), name.end(), static_cast<int(*)(int)>(isspace) ), name.end() );
			string resName = str.substr(6,3);
			// remove white spaces
			resName.erase( remove_if( resName.begin(), resName.end(), static_cast<int(*)(int)>(isspace) ), resName.end() );
			// atoms[i] = new Atom(x[i*3],x[i*3+1],x[i*3+2],r[i],qq,dd,name,resName,atoi(str.substr(11,4).c_str()),str.substr(10,1));

            const Atom at(x[i*3],x[i*3+1],x[i*3+2],r[i],qq,dd,name,resName,atoi(str.substr(11,4).c_str()),str.substr(10,1));
            atoms.push_back(at);
		}
		else
		{
			// atoms[i] = new Atom(x[i*3],x[i*3+1],x[i*3+2],r[i],qq,dd);
			const Atom at(x[i*3],x[i*3+1],x[i*3+2],r[i],qq,dd);
			atoms.push_back(at);
		}
	}
	return true;
}


void DelPhiShared::getBounds(double *cmin,double *cmax)
{
	cmin[0] = INFINITY;
	cmin[1] = INFINITY;
	cmin[2] = INFINITY;
	
	cmax[0] = -INFINITY;
	cmax[1] = -INFINITY;
	cmax[2] = -INFINITY;

	unsigned numAtoms = atoms.size();

	for (int ia=0; ia<numAtoms; ia++)
	{
		cmin[0] = min(cmin[0],atoms[ia].pos[0] - atoms[ia].radius);
		cmin[1] = min(cmin[1],atoms[ia].pos[1] - atoms[ia].radius);
		cmin[2] = min(cmin[2],atoms[ia].pos[2] - atoms[ia].radius);

		cmax[0] = max(cmax[0],atoms[ia].pos[0] + atoms[ia].radius);
		cmax[1] = max(cmax[1],atoms[ia].pos[1] + atoms[ia].radius);
		cmax[2] = max(cmax[2],atoms[ia].pos[2] + atoms[ia].radius);
	}
}


// SD: PB_NEW
bool DelPhiShared::buildGrid(double scale, double *cmin, double *cmax,double TOLL_ON_INTEGER_GRID)
{
	
	if (atoms.size() == 0)
	{
		cout << endl << ERR << "Cannot build grid with no atoms!";
		return false;
	}

	if (scale < 0)
	{
		cout << endl << WARN << "Cannot use a <0 scale: setting " << DEFAULT_SCALE;
		scale = DEFAULT_SCALE;
	}

	this->scale = scale;
	this->perfill = 100;
	
	side = 1./scale;

	double oldmid[3];

	cmin[0] -= side/2.0;
	cmin[1] -= side/2.0;
	cmin[2] -= side/2.0;
	
	cmax[0] += side/2.0;
	cmax[1] += side/2.0;
	cmax[2] += side/2.0;

	double dx,dy,dz;

	dx = cmax[0] - cmin[0];
	dy = cmax[1] - cmin[1];
	dz = cmax[2] - cmin[2];

	double tnx = dx*scale;
	double tny = dy*scale;
	double tnz = dz*scale;

	nx = round(dx*scale);
	ny = round(dy*scale);
	nz = round(dz*scale);

	if (fabs(tnx-nx) > TOLL_ON_INTEGER_GRID || fabs(tny-ny) > TOLL_ON_INTEGER_GRID || fabs(tnz-nz) > TOLL_ON_INTEGER_GRID)
	{
		cout << endl << ERR << "Non conformant grid!";
		cout << endl << ERR << "Expected integer nx grid size " << nx << " obtained real result " << tnx ;
		cout << endl << ERR << "Expected integer ny grid size " << ny << " obtained real result " << tny ;
		cout << endl << ERR << "Expected integer nz grid size " << nz << " obtained real result " << tnz ;
		return false;
	}
	
	cout << endl << INFO << "Input nx = " << nx-1;
	cout << endl << INFO << "Input ny = " << ny-1;
	cout << endl << INFO << "Input nz = " << nz-1;
	
	unsigned igrid = 0;

	int winnerIndex = 0;

	// get max dimension
	if (nx > ny)
	{
		igrid = nx;
		winnerIndex = 0;
	}
	else
	{
		igrid = ny;
		winnerIndex = 1;
	}
	
	if (nz > igrid)
	{
		igrid = nz;
		winnerIndex = 2;
	}

	// make it odd if necessary
	if ((igrid%2) == 0)
	{
		cout << endl << INFO << "Making the grid size odd" ;

		++igrid;

		dx += side;
		dy += side;
		dz += side;
	}

	// adjust the max so as to make a cube
	if (winnerIndex == 0)
	{
		cmax[1] = cmin[1] + dx;
		cmax[2] = cmin[2] + dx;
	}
	if (winnerIndex == 1)
	{
		cmax[0] = cmin[0] + dy;
		cmax[2] = cmin[2] + dy;
	}
	if (winnerIndex == 2)
	{
		cmax[0] = cmin[0] + dz;
		cmax[1] = cmin[1] + dz;
	}
	oldmid[0] = (cmax[0] + cmin[0])/2.;
	oldmid[1] = (cmax[1] + cmin[1])/2.;
	oldmid[2] = (cmax[2] + cmin[2])/2.;

	cout << endl << INFO << "Adjusted geometric baricenter ->  " << oldmid[0] <<
		" " << oldmid[1] << " " << oldmid[2];
				
	cout << endl << INFO << "Grid is " << igrid;

	double v[6];
	v[0] = fabs(cmax[0] - oldmid[0]);
	v[1] = fabs(cmin[0] - oldmid[0]);
	v[2] = fabs(cmax[1] - oldmid[1]);
	v[3] = fabs(cmin[1] - oldmid[1]);
	v[4] = fabs(cmax[2] - oldmid[2]);
	v[5] = fabs(cmin[2] - oldmid[2]);

	rmaxdim = -INFINITY;

	for (int i=0; i<6; i++)
		if (rmaxdim < v[i])
			rmaxdim = v[i];

	rmaxdim = 2*rmaxdim;


	baricenter[0] = oldmid[0];
	baricenter[1] = oldmid[1];
	baricenter[2] = oldmid[2];

	xmin = oldmid[0] - (igrid-1)/(2*scale);
	ymin = oldmid[1] - (igrid-1)/(2*scale);
	zmin = oldmid[2] - (igrid-1)/(2*scale);

	xmax = oldmid[0] + (igrid-1)/(2*scale);
	ymax = oldmid[1] + (igrid-1)/(2*scale);
	zmax = oldmid[2] + (igrid-1)/(2*scale);

	cout << endl << INFO << "MAX " << xmax << " " << ymax << " " << zmax;
	cout << endl << INFO << "MIN " << xmin << " " << ymin << " "  << zmin;

	cout << endl << INFO << "Allocating memory...";
	cout.flush();

	if (x!=NULL)
		deleteVector<double>(x);
	if (y!=NULL)
		deleteVector<double>(y);
	if (z!=NULL)
		deleteVector<double>(z);

	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);

	side = 1./scale;
	hside = 0.5/scale;
	A = side*side;

	nx = igrid;
	ny = igrid;
	nz = igrid;

	for (int i=0;i<nx;i++)
		x[i] = xmin + i*side;

	for (int i=0;i<ny;i++)
		y[i] = ymin + i*side;

	for (int i=0;i<nz;i++)
		z[i] = zmin + i*side;

	// allocate epsmap memory
	bool flag = true;
	
	if (buildEpsMap)
	{
		flag = clearAndAllocEpsMaps();
		
		if (!flag)
			return false;

		if (idebmap != NULL)
			deleteVector<bool>(idebmap);

		int64_t tot = nx*ny*nz;
		idebmap = allocateVector<bool>(tot);
		
		for (int64_t i=0; i<tot; i++)
			idebmap[i] = true;
	}

	// allocate grid point status memory
	// and initialize all to temporary outside status
	if (buildStatus)
	{
		if (!optimizeGrids)
		{
			if (status != NULL)
				deleteVector<int>(status);

			int64_t tot = nx*ny*nz;
			status = allocateVector<int>(tot);
			// int alignementBits = 64;
			// status = allocateAlignedVector<int>(tot,alignementBits);

			if (status == NULL)
			{
				cout << endl << ERR << "Not enough memory to allocate status map";
				cout << endl;
				exit(-1);
			}
			for (int64_t i=0; i<tot; i++)
				status[i] = STATUS_POINT_TEMPORARY_OUT;
		}
		else
		{
			if (bilevel_status != NULL)
			{
				deleteBilevelGridCells<int>(bilevel_status, nx, ny, nz);
				deleteVector<int *>(bilevel_status);
			}
			bilevel_status = allocateBilevelGridCells<int>(nx, ny, nz);
		}
	}
	cout << "ok!";
	return true;
}


bool DelPhiShared::buildGrid(double scale, double perfill)
{
	// if (atoms == NULL)
	if (atoms.size() == 0)
	{
		cout << endl << ERR << "Cannot build grid with no atoms!";
		return false;
	}

	if (scale < 0)
	{
		cout << endl << WARN << "Cannot use a <0 scale: setting " << DEFAULT_SCALE;
		scale = DEFAULT_SCALE;
	}

	if (perfill < 0 || perfill > 100)
	{
		cout << endl << WARN << "Perfil is in (0,100). Setting " << DEFAULT_PERFIL;
		scale = DEFAULT_PERFIL;
	}

	this->scale = scale;
	this->perfill = perfill;

	double cmin[3],cmax[3],oldmid[3];

	getBounds(cmin,cmax);

	oldmid[0] = (cmax[0] + cmin[0])/2.;
	oldmid[1] = (cmax[1] + cmin[1])/2.;
	oldmid[2] = (cmax[2] + cmin[2])/2.;

	cout << endl << INFO << "Geometric baricenter ->  " << oldmid[0] << 
		" " << oldmid[1] << " " << oldmid[2];

	double v[6];
	v[0] = fabs(cmax[0] - oldmid[0]);
	v[1] = fabs(cmin[0] - oldmid[0]);
	v[2] = fabs(cmax[1] - oldmid[1]);
	v[3] = fabs(cmin[1] - oldmid[1]);
	v[4] = fabs(cmax[2] - oldmid[2]);
	v[5] = fabs(cmin[2] - oldmid[2]);

	rmaxdim = -INFINITY;

	for (int i=0; i<6; i++)
		if (rmaxdim < v[i])
			rmaxdim = v[i];

	rmaxdim = 2*rmaxdim;

	       
	unsigned int igrid = (unsigned int)floor(scale*100./perfill*rmaxdim);

	if ((igrid%2)==0)
		igrid++;
	
	cout << endl << INFO << "Grid is " << igrid;

	baricenter[0] = oldmid[0];
	baricenter[1] = oldmid[1];
	baricenter[2] = oldmid[2];

	xmin = oldmid[0] - (igrid-1)/(2*scale);
	ymin = oldmid[1] - (igrid-1)/(2*scale);
	zmin = oldmid[2] - (igrid-1)/(2*scale);

	xmax = oldmid[0] + (igrid-1)/(2*scale);
	ymax = oldmid[1] + (igrid-1)/(2*scale);
	zmax = oldmid[2] + (igrid-1)/(2*scale);

	cout << endl << INFO << "MAX " << xmax << " " << ymax << " " << zmax;
	cout << endl << INFO << "MIN " << xmin << " " << ymin << " "  << zmin;
	cout << endl << INFO << "Perfil " << perfill << " %";
	cout << endl << INFO << "Rmaxdim " << rmaxdim;

	cout << endl << INFO << "Allocating memory...";
	cout.flush();

	if (x != NULL)
		deleteVector<double>(x);
	if (y != NULL)
		deleteVector<double>(y);
	if (z != NULL)
		deleteVector<double>(z);

	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);

	side = 1./scale;
	hside = 0.5/scale;
	A = side*side;

	nx = igrid;
	ny = igrid;
	nz = igrid;

	for (int i=0; i<nx; i++)
		x[i] = xmin + i*side;

	for (int i=0; i<ny; i++)
		y[i] = ymin + i*side;

	for (int i=0; i<nz; i++)
		z[i] = zmin + i*side;

	// allocate epsmap memory
	bool flag = true;

	if (buildEpsMap)
	{
		flag = clearAndAllocEpsMaps();
		
		if (!flag)
			return false;

		if (idebmap != NULL)
			deleteVector<bool>(idebmap);

		int64_t tot = nx*ny*nz;
		idebmap = allocateVector<bool>(tot);

		for (int64_t i=0; i<tot; i++)
			idebmap[i] = true;
	}

	// allocate grid point status memory
	// and initialize all to temporary outside status
	if (buildStatus)
	{
		if (!optimizeGrids)
		{
			if (status != NULL)
				deleteVector<int>(status);

			int64_t tot = nx*ny*nz;
			status = allocateVector<int>(tot);
			// int alignementBits = 64;
			// status = allocateAlignedVector<int>(tot,alignementBits);

			if (status == NULL)
			{
				cout << endl << ERR << "Not enough memory to allocate status map";
				cout << endl;
				exit(-1);
			}
			for (int64_t i=0; i<tot; i++)
				status[i] = STATUS_POINT_TEMPORARY_OUT;
		}
		else
		{
			if (bilevel_status != NULL)
			{
				deleteBilevelGridCells<int>(bilevel_status, nx, ny, nz);
				deleteVector<int *>(bilevel_status);
			}
			bilevel_status = allocateBilevelGridCells<int>(nx, ny, nz);
		}
	}

	cout << "ok!";
	return true;
}


bool DelPhiShared::clearAndAllocEpsMaps()
{
	if (epsmap != NULL)
		deleteVector<int>(epsmap);

	int64_t tot = nx*ny*nz*3;
	epsmap = allocateVector<int>(tot);
	
	if (epsmap == NULL)
	{
		cout << endl << ERR << "Error allocating epsmap!";
		return false;
	}
	
	// init epsmap
	for (int64_t i=0; i<tot; i++)
		epsmap[i] = 0;
	
	return true;

}


void DelPhiShared::clearEpsMaps()
{
	if (epsmap == NULL)
	{
		cout << endl << ERR << "Cannot clear null map!";
		return;
	}
	
	int64_t tot = nx*ny*nz*3;

	for (int64_t i = 0; i<tot; i++)
		epsmap[i] = 0;

	/*
	for (int k=0; k<nz; k++)
		for (int j=0; j<ny; j++)
			for (int i=0; i<nx; i++)
			{
				EPSMAP(i,j,k,0,nx,ny,nz) = 0;
				EPSMAP(i,j,k,1,nx,ny,nz) = 0;
				EPSMAP(i,j,k,2,nx,ny,nz) = 0;
			}
	*/

	return;
}


void DelPhiShared::saveIdebMap(char *fname)
{
	if  (idebmap == NULL)
	{
		cout << endl << WARN << "Cannot save null idebmap!";
		return;
	}

	char f1[BUFLEN];
	// if (debug_stern == -1)
	sprintf(f1,"%s%s.idebmap.txt",conf.rootFile.c_str(),fname);
	//else
	//	sprintf(f1,"%s.idebmap.stern.txt",fname);
	FILE *fp = fopen(f1,"w");

	for (int i=0; i<nx; i++)
	{
		for (int j=0; j<ny; j++)
		{
			for (int k=0; k<nz; k++)
			{
				// fprintf(fp, "%d ", (int)IDEBMAP(i,j,k,nx,ny));
				fprintf(fp, "%d ", read3DVector<bool>(idebmap,i,j,k,nx,ny,nz));
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}


void DelPhiShared::saveEpsMaps(char *fname)
{
	if  (epsmap == NULL)
	{
		cout << endl << WARN << "Cannot save null epsmap!";
		return;
	}
	
	char f1[BUFLEN],f2[BUFLEN],f3[BUFLEN];
	sprintf(f1, "%s%s.epsmapx.txt", conf.rootFile.c_str(),fname);
	sprintf(f2, "%s%s.epsmapy.txt", conf.rootFile.c_str(),fname);
	sprintf(f3, "%s%s.epsmapz.txt", conf.rootFile.c_str(),fname);

	FILE *fp = fopen(f1,"w");
	for (int i=0; i<nx; i++)
	{
		for (int j=0; j<ny; j++)
		{
			for (int k=0; k<nz; k++)
			{
				// fprintf(fp, "%d ", EPSMAP(i,j,k,0,nx,ny,nz));
				fprintf(fp, "%d ", read4DVector<int>(epsmap,i,j,k,0,nx,ny,nz,3));
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen(f2,"w");
	for (int i=0; i<nx; i++)
	{
		for (int j=0; j<ny; j++)
		{
			for (int k=0; k<nz; k++)
			{
				//fprintf(fp, "%d ", EPSMAP(i,j,k,1,nx,ny,nz));
				fprintf(fp, "%d ", read4DVector<int>(epsmap,i,j,k,1,nx,ny,nz,3));
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen(f3,"w");
	for (int i=0; i<nx; i++)
	{
		for (int j=0; j<ny; j++)
		{
			for (int k=0; k<nz; k++)
			{
				//fprintf(fp, "%d ", EPSMAP(i,j,k,2,nx,ny,nz));
				fprintf(fp, "%d ", read4DVector<int>(epsmap,i,j,k,2,nx,ny,nz,3));
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}


DelPhiShared::~DelPhiShared()
{
	clear();
}


void DelPhiShared::clear()
{
	if (epsmap != NULL)
		deleteVector<int>(epsmap);

	if  (ibgp != NULL)
		deleteVector<int>(ibgp);
	if  (scspos != NULL)
		deleteVector<double>(scspos);
	if  (scsnor != NULL)
		deleteVector<double>(scsnor);
	if  (scsarea != NULL)
		deleteVector<double>(scsarea);
	if (x != NULL)
		deleteVector<double>(x);
	if (y != NULL)
		deleteVector<double>(y);
	if (z != NULL)
		deleteVector<double>(z);

	if (!optimizeGrids)
	{
		if (status != NULL)
			// deleteMatrix3D<int>(nx,ny,nz,status);
			// deleteAlignedVector<int>(status);
			deleteVector<int>(status);
	}
	else
	{
		if (bilevel_status != NULL)
		{
			deleteBilevelGridCells<int>(bilevel_status, nx, ny, nz);
			deleteVector<int*>(bilevel_status);
		}
	}

	if (idebmap != NULL)
		deleteVector<bool>(idebmap);

    if (atoms.size() != 0)
        atoms.clear();

	// free memory 
	if (cavitiesVec != NULL)
	{
		for (vector<vector<int*>*>::iterator it = cavitiesVec->begin(); it != cavitiesVec->end(); it++)
		{
			vector<int*> *inner = (*it);
			vector<int*>::iterator it2;
			for (it2 = inner->begin(); it2 != inner->end(); it2++)
				deleteVector<int>((*it2));
			delete inner;
		}
		delete cavitiesVec;
	}

	if (cav2atoms.size() != 0)
	{
		for (unsigned int ic = 0; ic < cav2atoms.size(); ic++)
			delete cav2atoms[ic];
	}
	cav2atoms.clear();
}


void DelPhiShared::saveBGP(char *fname)
{
	if  (scspos == NULL || scsnor == NULL)
	{
		cout << endl << WARN << "Cannot save null scspos or scsnor!";
		return;
	}

	char ff[BUFLEN];
	sprintf(ff, "%s%s.projections.txt", conf.rootFile.c_str(),fname);
	FILE *fp = fopen(ff, "w");
	sprintf(ff, "%ssurf.txt", conf.rootFile.c_str());
	FILE *fp2 = fopen(ff, "w");
	fprintf(fp, "%d\n", nbgp);

	for (int iibgp=0; iibgp<nbgp; iibgp++)
	{
		fprintf(fp,"%d %d %d %lf %lf %lf %lf %lf %lf",
				ibgp[iibgp*3]+1, ibgp[iibgp*3+1]+1, ibgp[iibgp*3+2]+1,
				scspos[iibgp*3], scspos[iibgp*3+1], scspos[iibgp*3+2],
				scsnor[iibgp*3], scsnor[iibgp*3+1], scsnor[iibgp*3+2]);

		fprintf(fp2,"%lf %lf %lf %lf %lf %lf",
				scspos[iibgp*3], scspos[iibgp*3+1], scspos[iibgp*3+2],
				scsnor[iibgp*3], scsnor[iibgp*3+1], scsnor[iibgp*3+2]);

		if (iibgp < nbgp-1)
		{
			fprintf(fp, "\n");
			fprintf(fp2, "\n");
		}
	}
	fclose(fp);
	fclose(fp2);
}


void DelPhiShared::saveStatus(char *fname)
{
	if (!optimizeGrids)
	{
		if  (status == NULL)
		{
			cout << endl << WARN << "Cannot save null status!";
			return;
		}
	}
	else
	{
		if  (bilevel_status == NULL)
		{
			cout << endl << WARN << "Cannot save null bilevel status!";
			return;
		}
	}

	char ff[BUFLEN];
	sprintf(ff, "%s%s.status.txt", conf.rootFile.c_str(),fname);
	FILE *fp = fopen(ff, "w");

	
	for (int i=0; i<nx; i++)
	{		
		for (int j=0; j<ny; j++)
		{
			for (int k=0; k<nz; k++)
			{
				// fprintf(fp, "%d ", status[k][j][i]);
				// fprintf(fp, "%d ", STATUSMAP(i,j,k,nx,ny));
				if (!optimizeGrids)
					fprintf(fp, "%d ", read3DVector<int>(status,i,j,k,nx,ny,nz));
				else
					fprintf(fp, "%d ", readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,nx,ny,nz));
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	saveCavities(false);
}


void DelPhiShared::initCav2Atoms()
{
	int nc = (int)cavitiesVec->size();

	if (cav2atoms.size() != 0)
		clearCav2Atoms();

	cav2atoms.reserve(nc);
	
	for (int ic=0; ic<nc; ic++)
		cav2atoms.push_back(new set<int>());
}


void DelPhiShared::clearCav2Atoms()
{
	if (cav2atoms.size() != 0)
	{
		for (unsigned int ic = 0; ic < cav2atoms.size(); ic++)
			delete cav2atoms[ic];
	}
	cav2atoms.clear();
}


void DelPhiShared::markPockets(int *external_status, vector<bool> &isPocket)
{
	int num = (int)cavitiesVec->size();
	
	if (external_status == NULL)
	{
		cout << endl << WARN << "Cannot mark pockets with a reference status map!";
		return;
	}

	if (external_status == status)
	{
		cout << endl << WARN << "In marking pockets the passed status map must be from another reference object; here the same was passed";
		return;
	}

	int *status = external_status;

	if (num == 0)
		return;

	int sync_i = 0;
	
	for (vector<vector<int*>*>::iterator it = cavitiesVec->begin(); it != cavitiesVec->end(); it++, sync_i++)
	{
		if (cavitiesFlag[sync_i] == true)
		{
			continue;
		}
		
		vector<int*>::iterator it2;
		vector<int*> *vec = (*it);

		bool isPocketFlag = false;
		
		for (it2 = vec->begin(); it2 != vec->end(); it2++)
		{
			int *vec = *it2;
			int value = read3DVector<int>(status, vec[0], vec[1], vec[2], nx, ny, nz);
			// if at least one was out, then that's a pocket
			// if at least one was out, then that's a pocket
			// if (STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny) == STATUS_POINT_TEMPORARY_OUT ||
			// 	STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny) == STATUS_POINT_OUT)
			if (value == STATUS_POINT_TEMPORARY_OUT || value == STATUS_POINT_OUT)
			{
				isPocketFlag = true;
				break;
			}
		}
		isPocket.push_back(isPocketFlag);
	}	
	return;
}


// This version is identical to the above one except for the use of bilevel grids for memory space reduction
void DelPhiShared::markPockets(int **external_bilevel_status, vector<bool> &isPocket)
{
	int num = (int)cavitiesVec->size();
	
	if (external_bilevel_status == NULL)
	{
		cout << endl << WARN << "Cannot mark pockets with a reference bilevel status map!";
		return;
	}

	if (external_bilevel_status == bilevel_status)
	{
		cout << endl << WARN << "In marking pockets the passed bilevel status map must be from another reference object; here the same was passed";
		return;
	}

	int **bilevel_status = external_bilevel_status;

	if (num == 0)
		return;

	int sync_i = 0;
	
	for (vector<vector<int*>*>::iterator it = cavitiesVec->begin(); it != cavitiesVec->end(); it++, sync_i++)
	{
		if (cavitiesFlag[sync_i] == true)
		{
			continue;
		}
		
		vector<int*>::iterator it2;
		vector<int*> *vec = (*it);

		bool isPocketFlag = false;
		
		for (it2 = vec->begin(); it2 != vec->end(); it2++)
		{
			int *vec = *it2;
			
			int value = readBilevelGrid<int>(bilevel_status,
											 STATUS_POINT_TEMPORARY_OUT, vec[0], vec[1], vec[2], nx, ny, nz);
			
			if (value == STATUS_POINT_TEMPORARY_OUT || value == STATUS_POINT_OUT)
			{
				isPocketFlag = true;
				break;
			}
		}
		isPocket.push_back(isPocketFlag);
	}
	return;
}


void DelPhiShared::parseAtomInfo(char *atinf, string mol, int natom, double *xn1, double *rad)
{
	// get atoms and info from delphi
	FILE *ff = fopen(mol.c_str(),"w");
	int lastChar = 0;
		
	// parse atinfo into a file
	for (int i=0; i<natom; i++)
	{
		char buff[100];
		for (int kk=0; kk<15; kk++)
			buff[kk] = atinf[lastChar++];

		buff[15] = '\0';
		string str(buff);
		string name = str.substr(0,5);
		// avoid disambiguation
		name.erase( remove_if(name.begin(), name.end(), static_cast<int(*)(int)>(isspace)), name.end() );
		string resName = str.substr(6,3);
		resName.erase( remove_if(resName.begin(), resName.end(), static_cast<int(*)(int)>(isspace)), resName.end() );

		fprintf(ff, "%f %f %f %f 0 %s %s %d %s\n",
				xn1[3*i], xn1[3*i+1], xn1[3*i+2], rad[i], name.c_str(), resName.c_str(),
				atoi(str.substr(11,4).c_str()), str.substr(10,1).c_str());
	}
	fclose(ff);
}


int DelPhiShared::cavitiesToAtoms(double rad)
{
	int num = (int)cavitiesVec->size();

	if (num == 0)
		return 0;

	char buff[100];
	int i = 0, sync_i = 0;
	int nonActive = 0;
	FILE *fp;

	for (vector<vector<int*>*>::iterator it = cavitiesVec->begin(); it != cavitiesVec->end(); it++, i++, sync_i++)
	{
		if (cavitiesFlag[sync_i] == true)
		{
			nonActive++;
			i--;
			continue;
		}
		sprintf(buff, "%scav%d.txt", conf.rootFile.c_str(), i);
		fp = fopen(buff, "w");
		sprintf(buff, "%sall_cav%d.txt", conf.rootFile.c_str(), i);
		FILE *fp2 = fopen(buff,"w");

		vector<int*> *vec = (*it);
		int cavityId = sync_i + STATUS_FIRST_CAV;

		int count = 0;
		double lastx,lasty,lastz;

		for (vector<int*>::iterator it2 = vec->begin(); it2 != vec->end(); it2++)
		{
			int *vec = *it2;
			fprintf(fp2, "%f %f %f %f\n", x[vec[0]], y[vec[1]], z[vec[2]], rad);

			// save only support cavity points
			// if (STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny) == -cavityId)
			if (read3DVector<int>(status,vec[0],vec[1],vec[2],nx,ny,nz) == -cavityId)
			{
				// printf("\n %d",STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny));
				fprintf(fp, "%f %f %f %f\n", x[vec[0]],y[vec[1]],z[vec[2]],rad);
				// printf("\n%d %d %d %f %f %f",vec[0],vec[1],vec[2],x[vec[0]],y[vec[1]],z[vec[2]]);
				lastx = x[vec[0]];
				lasty = y[vec[1]];
				lastz = z[vec[2]];
				count++;
			}
		}

		if (count == 0)
		{
			cout << endl << ERR << "Zero support atoms to save";
			exit(-1);
		}

		// clone atoms to let NanoShaper work on this set of points
		// NanoShaper needs at least 4 atoms for triangulation
		if (count < 4)
		{
			count = 4 - count;
			for (int i=0; i<count; i++)
				fprintf(fp, "%f %f %f %f\n", lastx,lasty,lastz,rad);
		}
		fclose(fp);
		fclose(fp2);
	}
	char temp[1024];
	sprintf(temp, "%snumcav.txt", conf.rootFile.c_str());
	fp = fopen(temp, "w");
	fprintf(fp, "%d", num-nonActive);

	fclose(fp);

	return num - nonActive;
}


int DelPhiShared::cavitiesToAtomsWithBilevelStatusMap(double rad)
{
	int num = (int)cavitiesVec->size();

	if (num == 0)
		return 0;

	char buff[100];
	int i = 0, sync_i = 0;
	int nonActive = 0;
	FILE *fp;

	for (vector<vector<int*>*>::iterator it = cavitiesVec->begin(); it != cavitiesVec->end(); it++, i++, sync_i++)
	{
		if (cavitiesFlag[sync_i] == true)
		{
			nonActive++;
			i--;
			continue;
		}
		sprintf(buff, "%scav%d.txt", conf.rootFile.c_str(), i);
		fp = fopen(buff, "w");
		sprintf(buff, "%sall_cav%d.txt", conf.rootFile.c_str(), i);
		FILE *fp2 = fopen(buff,"w");
		
		vector<int*> *vec = (*it);
		int cavityId = sync_i + STATUS_FIRST_CAV;

		int count = 0;
		double lastx,lasty,lastz;

		for (vector<int*>::iterator it2 = vec->begin(); it2 != vec->end(); it2++)
		{
			int *vec = *it2;
			fprintf(fp2, "%f %f %f %f\n", x[vec[0]], y[vec[1]], z[vec[2]], rad);
			
			// save only support cavity points
			// if (STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny) == -cavityId)
			if (readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,vec[0],vec[1],vec[2],nx,ny,nz) == -cavityId)
			{
				fprintf(fp, "%f %f %f %f\n", x[vec[0]],y[vec[1]],z[vec[2]],rad);
				// printf("\n%d %d %d %f %f %f",vec[0],vec[1],vec[2],x[vec[0]],y[vec[1]],z[vec[2]]);
				lastx = x[vec[0]];
				lasty = y[vec[1]];
				lastz = z[vec[2]];
				count++;
			}
		}

		if (count == 0)
		{
			cout << endl << ERR << "Zero support atoms to save";
			exit(-1);
		}
		
		// clone atoms to let NanoShaper work on this set of points
		// NanoShaper needs at least 4 atoms for triangulation
		if (count < 4)
		{
			count = 4 - count;
			for (int i=0; i<count; i++)
				fprintf(fp, "%f %f %f %f\n", lastx,lasty,lastz,rad);
		}
		fclose(fp);
		fclose(fp2);
	}	
	char temp[1024];
	sprintf(temp, "%snumcav.txt", conf.rootFile.c_str());
	fp = fopen(temp, "w");
	fprintf(fp, "%d", num-nonActive);

	fclose(fp);

	return num - nonActive;
}


void DelPhiShared::saveCavities2(bool onlySaveNonFilled,string sysName)
{
	if (isAvailableAtomInfo)
	{
		FILE *fp;
		char buff[100];
		sprintf(buff, "%s%s.pocket", conf.rootFile.c_str(), sysName.c_str());
		fp = fopen(buff, "w");
		if (fp == NULL)
		{
			cout << endl << WARN << "Cannot open output file " << buff;
			return;
		}
		int savedIndex = 1;
		int savedAtoms = 0;
		int savedCavities = 0;
		int skippedInvalidIndices = 0;

		for (unsigned int i=0; i<cav2atoms.size(); i++)
		{
			set<int>::iterator setIt;
			set<int> *setPt = cav2atoms[i];

			// save only non filled if required, or save all of them
			if ((onlySaveNonFilled && !cavitiesFlag[i]) || (!onlySaveNonFilled))
			{
				bool savedCurrentCavity = false;
				for (setIt = setPt->begin(); setIt!=setPt->end();setIt++)
				{
					// printf("\nIndex %d",i);
					fflush(stdout);
					int ii = (*setIt);
					if (ii < 0 || ii >= (int)atoms.size())
					{
						skippedInvalidIndices++;
						continue;
					}
					// Atom *at = atoms[ii];
					// AtomInfo &ai = at->ai;
					Atom at = atoms[ii];
					AtomInfo ai = at.ai;

					const int atom_serial = (ai.getSerial() > 0) ? ai.getSerial() : (ii + 1);

					if ((ai.getName()).size() == 4)
						// fprintf(fp, "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f                 %6d POC\n", ii+1, ai.getName().c_str(),ai.getResName().c_str(), ai.getChain().c_str(), ai.getResNum(), at->pos[0], at->pos[1], at->pos[2],savedIndex);
						fprintf(fp, "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f                 %6d POC\n", atom_serial, ai.getName().c_str(), ai.getResName().c_str(), ai.getChain().c_str(), ai.getResNum(), at.pos[0],at.pos[1], at.pos[2],savedIndex);
					else
						// fprintf(fp, "ATOM  %5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f                 %6d POC\n", ii+1, ai.getName().c_str(), ai.getResName().c_str(), ai.getChain().c_str(), ai.getResNum(), at->pos[0], at->pos[1], at->pos[2],savedIndex);
						fprintf(fp, "ATOM  %5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f                 %6d POC\n", atom_serial, ai.getName().c_str(), ai.getResName().c_str(), ai.getChain().c_str(), ai.getResNum(), at.pos[0], at.pos[1], at.pos[2], savedIndex);
					savedAtoms++;
					savedCurrentCavity = true;
				}
				if (savedCurrentCavity)
					savedCavities++;
				savedIndex++;
			}
		}
		fclose(fp);
		cout << endl << INFO << "Cavity atom records saved in " << buff << " (" << savedAtoms << " atoms in " << savedCavities << " cavities)";
		if (skippedInvalidIndices > 0)
			cout << endl << WARN << "Skipped " << skippedInvalidIndices << " invalid cavity atom indices while saving " << buff;
	}
}

void DelPhiShared::saveCavities(bool onlySaveNonFilled)
{
	if (!optimizeGrids)
	{
		if  (status == NULL)
		{
			cout << endl << WARN << "Cannot save cavities with null status!";
			return;
		}
	}
	else
	{
		if  (bilevel_status == NULL)
		{
			cout << endl << WARN << "Cannot save cavities with null status!";
			return;
		}
	}
	
	char temp[1024];
	sprintf(temp, "%scavities.txt", conf.rootFile.c_str());
	FILE *fp2 = fopen(temp,"w");
	
	for (int k=0; k<nz; k++)
	{
		for (int j=0; j<ny; j++)
		{
			for (int i=0; i<nx; i++)
			{
				// if (status[k][j][i] > 3)
				// if ((STATUSMAP(i,j,k,nx,ny))>STATUS_POINT_OUT)
				if (!optimizeGrids)
				{
					if (read3DVector<int>(status,i,j,k,nx,ny,nz) > STATUS_POINT_OUT)
						fprintf(fp2, "\n %f %f %f", x[i],y[j],z[k]);
				}
				else
				{
					if (readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,nx,ny,nz) > STATUS_POINT_OUT)
						fprintf(fp2, "\n %f %f %f", x[i],y[j],z[k]);
				}
			}
		}
	}
	fclose(fp2);
	
	sprintf(temp, "%scavitiesSize.txt", conf.rootFile.c_str());
	fp2 = fopen(temp,"w");
	
	int i = 0;
	
	fprintf(fp2, "%f\n", (float)cavitiesSize.size());
	
	// for each cavity save the baricenter and the estimated volume
	// save all cavities; add a flag for the filled one
	for (vector<double>::iterator itCav = cavitiesSize.begin(); itCav != cavitiesSize.end(); itCav++, i++)
	{
		double barx,bary,barz;
		barx = 0;
		bary = 0;
		barz = 0;
		for (vector<int*>::iterator itCavVec = cavitiesVec->at(i)->begin(); itCavVec != cavitiesVec->at(i)->end(); itCavVec++)
		{
			barx += x[(*itCavVec)[0]];
			bary += y[(*itCavVec)[1]];
			barz += z[(*itCavVec)[2]];
		}
		double volSize = cavitiesSize.at(i);
		int size = (int)cavitiesVec->at(i)->size();
		barx /= size;
		bary /= size;
		barz /= size;
		int fillFlag;
		
		if (cavitiesFlag[i])
			fillFlag = 1;
		else
			fillFlag = 0;
		// save baricenter,volume, and filling flag
		fprintf(fp2, "%lf %lf %lf %lf %d\n", barx, bary, barz, volSize, fillFlag);
	}
	fclose(fp2);
	
	sprintf(temp, "%scavAtomsSerials.txt", conf.rootFile.c_str());
	fp2 = fopen(temp,"w");
	
	for (unsigned int i=0; i<cav2atoms.size(); i++)
	{
		// save only non filled if required, or save all of them
		if ((onlySaveNonFilled && !cavitiesFlag[i]) || (!onlySaveNonFilled))
		{
			set<int> *setPt = cav2atoms[i];
				for (set<int>::iterator setIt = setPt->begin(); setIt != setPt->end(); setIt++)
				{
					const int atom_index = (*setIt);
					int atom_serial = atom_index + 1;
					if (atom_index >= 0 && atom_index < (int)atoms.size() && atoms[atom_index].ai.getSerial() > 0)
						atom_serial = atoms[atom_index].ai.getSerial();
					fprintf(fp2, "%d ", atom_serial);
				}
			
			fprintf(fp2,"\n");
		}
	}
	fclose(fp2);
	
	if (isAvailableAtomInfo)
	{
		FILE *fp3;
		//cout << endl << INFO << "Atoms info is available I am saving also residues";
		char temp[1024];
		sprintf(temp, "%sresidues.txt", conf.rootFile.c_str());
		fp2 = fopen(temp,"w");
		sprintf(temp, "%sresidues_vmd.txt", conf.rootFile.c_str());
		fp3 = fopen(temp, "w");
		
		for (unsigned int i=0; i<cav2atoms.size(); i++)
		{
			// only save the non filled cavities
			if (onlySaveNonFilled && !cavitiesFlag[i])
			{
				// detected the involved residues in the cavity
				set<int> involvedResidues;
				set<int> *setPt = cav2atoms[i];
				set<int>::iterator setIt;
				
				for (setIt = setPt->begin(); setIt != setPt->end(); setIt++)
				{
					int atomIndex = (*setIt);
					// involvedResidues.insert(atoms[atomIndex]->ai.getResNum());
					involvedResidues.insert(atoms[atomIndex].ai.getResNum());
				}
				
				int ll, tot = (int)involvedResidues.size();
				
				for (ll=0, setIt = involvedResidues.begin(); setIt != involvedResidues.end(); setIt++, ll++)
				{
					int resIndex = (*setIt);
					fprintf(fp2, "%d ",resIndex);
					fprintf(fp3, "(resid %d and not water) ", resIndex);
					
					if (ll != tot-1)
						fprintf(fp3, " or ");
				}
				fprintf(fp2, "\n");
				fprintf(fp3, "\n");
			}
		}
		fclose(fp2);
		fclose(fp3);
	}
}
