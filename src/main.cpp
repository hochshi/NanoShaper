/*
  * Copyright (C) 2015 Istituto Italiano di Tecnologia - IIT Genova (Italy)
  * Author: Sergio Decherchi
  * Email: sergio.decherchi@iit.it
  * Permission is granted to copy, distribute, and/or modify this program
  * under the terms of the GNU General Public License, version 3 or any
  * later version published by the Free Software Foundation.
  *
  * This program is distributed in the hope that it will be useful, but
  * WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
  * Public License for more details
  * 
  * You should have received a copy of the GNU General Public License
  * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#include "globals.h"

// #define NS_TEST

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif

#include "SurfaceFactory.h"
#include "Surface.h"
#include "tools.h"
#include "DelphiShared.h"
#include "ConfigFile.h"

#ifdef NS_TEST
	#include "SkinSurface.h"
#endif

// have to use a pointer because Swig will try to make a copy and fstream does
// not have a copy semantic
fstream *errorStream;
fstream *internals;

// global conf file
Configuration conf;

#if defined(USE_VIS_TOOLS)
Vis vis;
#endif

// initial path where NanoShaper is run
#ifdef BOOST_ENABLED
	boost::filesystem::path initialPath;
#endif

ConfigFile *init(string argv);
void dispose(ConfigFile *cf);
void stopDebug();
void restartDebug();
void cite();
void normalMode(Surface *surf,DelPhiShared *dg);
void membfitMode(Surface *surf,DelPhiShared *dg);
void pocketMode(bool hasAtomInfo,ConfigFile *cf);


class pocketWrapper
{
public:

	double *volumes, *areas, *areas2;

private:
	struct AtomAreaStats
	{
		int noWinnerVertices;
		int invalidAtomIndex;
		int degenerateTriangles;

		AtomAreaStats() : noWinnerVertices(0), invalidAtomIndex(0), degenerateTriangles(0) {}
	};

	static double dot3(const double *a, const double *b)
	{
		return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	}

	static void sub3(const double *a, const double *b, double *out)
	{
		out[0] = a[0] - b[0];
		out[1] = a[1] - b[1];
		out[2] = a[2] - b[2];
	}

	static double norm3(const double *a)
	{
		return sqrt(dot3(a, a));
	}

	static double triArea(const double *p0, const double *p1, const double *p2)
	{
		double e0[3], e1[3], c[3];
		sub3(p1, p0, e0);
		sub3(p2, p0, e1);
		c[0] = e0[1]*e1[2] - e0[2]*e1[1];
		c[1] = e0[2]*e1[0] - e0[0]*e1[2];
		c[2] = e0[0]*e1[1] - e0[1]*e1[0];
		return 0.5 * norm3(c);
	}

	static bool computeVoronoiContrib(const double *p0, const double *p1, const double *p2,
									  double contrib[3], double &areaOut)
	{
		const double epsArea = 1e-16;
		const double epsCot = 1e-16;

		double e01[3], e02[3], e10[3], e12[3], e20[3], e21[3];
		sub3(p1, p0, e01);
		sub3(p2, p0, e02);
		sub3(p0, p1, e10);
		sub3(p2, p1, e12);
		sub3(p0, p2, e20);
		sub3(p1, p2, e21);

		double c[3] = {
			e01[1]*e02[2] - e01[2]*e02[1],
			e01[2]*e02[0] - e01[0]*e02[2],
			e01[0]*e02[1] - e01[1]*e02[0]
		};
		double crossNorm = norm3(c);
		areaOut = 0.5 * crossNorm;

		contrib[0] = 0.0;
		contrib[1] = 0.0;
		contrib[2] = 0.0;

		if (areaOut <= epsArea)
			return false;

		double dot0 = dot3(e01, e02);
		double dot1 = dot3(e10, e12);
		double dot2 = dot3(e20, e21);

		bool obt0 = (dot0 < 0.0);
		bool obt1 = (dot1 < 0.0);
		bool obt2 = (dot2 < 0.0);

		if (obt0 || obt1 || obt2)
		{
			if (obt0)
			{
				contrib[0] = 0.5 * areaOut;
				contrib[1] = 0.25 * areaOut;
				contrib[2] = 0.25 * areaOut;
			}
			else if (obt1)
			{
				contrib[0] = 0.25 * areaOut;
				contrib[1] = 0.5 * areaOut;
				contrib[2] = 0.25 * areaOut;
			}
			else
			{
				contrib[0] = 0.25 * areaOut;
				contrib[1] = 0.25 * areaOut;
				contrib[2] = 0.5 * areaOut;
			}
			return true;
		}

		if (crossNorm <= epsCot)
			return false;

		double cot0 = dot0 / crossNorm;
		double cot1 = dot1 / crossNorm;
		double cot2 = dot2 / crossNorm;

		double l01 = dot3(e01, e01);
		double l02 = dot3(e02, e02);
		double l12 = dot3(e12, e12);

		contrib[0] = (l01*cot2 + l02*cot1) / 8.0;
		contrib[1] = (l12*cot0 + l01*cot2) / 8.0;
		contrib[2] = (l02*cot1 + l12*cot0) / 8.0;

		for (int k = 0; k < 3; k++)
		{
			if (!std::isfinite(contrib[k]))
				return false;
			if (contrib[k] < 0.0 && contrib[k] > -1e-12)
				contrib[k] = 0.0;
		}

		double sumC = contrib[0] + contrib[1] + contrib[2];
		if (sumC <= epsArea)
			return false;

		double scale = areaOut / sumC;
		contrib[0] *= scale;
		contrib[1] *= scale;
		contrib[2] *= scale;

		return true;
	}

	static void addContributionToAtom(int atomIdx, double area, vector<double> &atomArea,
									  double &failedAssignmentArea, AtomAreaStats &stats)
	{
		if (area == 0.0)
			return;

		if (atomIdx >= 0 && atomIdx < (int)atomArea.size())
			atomArea[atomIdx] += area;
		else
		{
			failedAssignmentArea += area;
			stats.invalidAtomIndex++;
		}
	}

	static void computeAtomAreasForMesh(Surface *mesh, Surface *atomRef, const vector<bool> &externalMask,
										vector<double> &atomAreaFull, vector<double> &atomAreaBody,
										double &failedAssignmentAreaFull, double &failedAssignmentAreaBody,
										double &filteredBodyArea, AtomAreaStats &statsFull, AtomAreaStats &statsBody,
										vector<char> &atomTouched)
	{
		int nv = mesh->getTriangulatedNumVertices();
		int nt = mesh->getTriangulatedNumTriangles();
		int numAtoms = atomRef->getNumLoadedAtoms();
		bool validExternalMask = ((int)externalMask.size() == nv);

		atomAreaFull.assign(numAtoms, 0.0);
		atomAreaBody.assign(numAtoms, 0.0);
		atomTouched.assign(numAtoms, 0);

		vector<double> vertices(3 * nv, 0.0);
		vector<int> vertexAtomFull(nv, -1);
		vector<int> vertexAtomBody(nv, -1);

		for (int i = 0; i < nv; i++)
		{
			double *p = &vertices[i * 3];
			mesh->getTriangulatedVertex(i, p);
			int winner = -1;
			atomRef->nearestAtomForPoint(p, winner);
			vertexAtomFull[i] = winner;
			vertexAtomBody[i] = winner;
			if (winner >= 0 && winner < numAtoms)
				atomTouched[winner] = 1;
			if (validExternalMask && externalMask[i])
				vertexAtomBody[i] = -1;
			if (winner < 0)
				statsFull.noWinnerVertices++;
		}

		for (int t = 0; t < nt; t++)
		{
			int v0, v1, v2;
			mesh->getTriangulatedTriangle(t, v0, v1, v2);

			double *p0 = &vertices[v0 * 3];
			double *p1 = &vertices[v1 * 3];
			double *p2 = &vertices[v2 * 3];

			double contrib[3];
			double area = 0.0;
			bool valid = computeVoronoiContrib(p0, p1, p2, contrib, area);
			bool keepBody = true;
			if (validExternalMask)
				keepBody = (!externalMask[v0] && !externalMask[v1] && !externalMask[v2]);

			if (!valid)
			{
				double tri_area = triArea(p0, p1, p2);
				failedAssignmentAreaFull += tri_area;
				failedAssignmentAreaBody += tri_area;
				statsFull.degenerateTriangles++;
				statsBody.degenerateTriangles++;
				continue;
			}

			addContributionToAtom(vertexAtomFull[v0], contrib[0], atomAreaFull, failedAssignmentAreaFull, statsFull);
			addContributionToAtom(vertexAtomFull[v1], contrib[1], atomAreaFull, failedAssignmentAreaFull, statsFull);
			addContributionToAtom(vertexAtomFull[v2], contrib[2], atomAreaFull, failedAssignmentAreaFull, statsFull);

			if (keepBody)
			{
				addContributionToAtom(vertexAtomBody[v0], contrib[0], atomAreaBody, failedAssignmentAreaBody, statsBody);
				addContributionToAtom(vertexAtomBody[v1], contrib[1], atomAreaBody, failedAssignmentAreaBody, statsBody);
				addContributionToAtom(vertexAtomBody[v2], contrib[2], atomAreaBody, failedAssignmentAreaBody, statsBody);
			}
			else
			{
				filteredBodyArea += (contrib[0] + contrib[1] + contrib[2]);
			}
		}
	}

	static bool writeCombinedAtomAreasFile(int cavityIndex,
										   const vector<double> &atomAreaFull,
										   const vector<double> &atomAreaBody,
										   const vector<int> &atomOutputIds,
										   double fullMeshArea, double bodyMeshArea,
										   double filteredBodyArea,
										   double failedAssignmentAreaFull,
										   double failedAssignmentAreaBody,
										   const vector<char> &atomTouched,
										   bool includeZeroAtoms = true)
	{
		char fullName[1024];
		sprintf(fullName, "%scav_tri%d.atomSurfaceAreas.txt", conf.rootFile.c_str(), cavityIndex);
		FILE *fp = fopen(fullName, "w");
		if (fp == NULL)
		{
			cout << endl << WARN << "Cannot write file " << fullName;
			return false;
		}

		double assignedFull = 0.0;
		double assignedBody = 0.0;
		for (int i = 0; i < (int)atomAreaFull.size(); i++)
			assignedFull += atomAreaFull[i];
		for (int i = 0; i < (int)atomAreaBody.size(); i++)
			assignedBody += atomAreaBody[i];

		fprintf(fp, "# cavity_index %d\n", cavityIndex);
		fprintf(fp, "# summary full_mesh_area_A2 %.10f\n", fullMeshArea);
		fprintf(fp, "# summary body_mesh_area_A2 %.10f\n", bodyMeshArea);
		fprintf(fp, "# summary full_assigned_atoms_A2 %.10f\n", assignedFull);
		fprintf(fp, "# summary body_assigned_atoms_A2 %.10f\n", assignedBody);
		fprintf(fp, "# summary body_filtered_external_A2 %.10f\n", filteredBodyArea);
		fprintf(fp, "# summary full_failed_assignment_A2 %.10f\n", failedAssignmentAreaFull);
		fprintf(fp, "# summary body_failed_assignment_A2 %.10f\n", failedAssignmentAreaBody);
		fprintf(fp, "# section atom_index area_A2\n");

		int writtenFull = 0;
		int writtenBody = 0;
		for (int i = 0; i < (int)atomAreaFull.size(); i++)
		{
			bool inCavityAtomSet = (i < (int)atomTouched.size() && atomTouched[i] != 0);
			if ((includeZeroAtoms && inCavityAtomSet) || (!includeZeroAtoms && atomAreaFull[i] != 0.0))
			{
				const int atomOutputId = (i < (int)atomOutputIds.size() && atomOutputIds[i] > 0) ? atomOutputIds[i] : (i + 1);
				fprintf(fp, "full %d %.10f\n", atomOutputId, atomAreaFull[i]);
				writtenFull++;
			}
		}
		fprintf(fp, "full -2 %.10f\n", failedAssignmentAreaFull);

		for (int i = 0; i < (int)atomAreaBody.size(); i++)
		{
			bool inCavityAtomSet = (i < (int)atomTouched.size() && atomTouched[i] != 0);
			if ((includeZeroAtoms && inCavityAtomSet) || (!includeZeroAtoms && atomAreaBody[i] != 0.0))
			{
				const int atomOutputId = (i < (int)atomOutputIds.size() && atomOutputIds[i] > 0) ? atomOutputIds[i] : (i + 1);
				fprintf(fp, "body %d %.10f\n", atomOutputId, atomAreaBody[i]);
				writtenBody++;
			}
		}
		fprintf(fp, "body -1 %.10f\n", filteredBodyArea);
		fprintf(fp, "body -2 %.10f\n", failedAssignmentAreaBody);

		fclose(fp);

		cout << endl << INFO << "Per-atom surface areas saved in " << fullName
			 << " (full atoms " << writtenFull << ", body atoms " << writtenBody
			 << ", full area " << setprecision(10) << assignedFull
			 << " [A^2], body area " << assignedBody << " [A^2], -1 " << filteredBodyArea
			 << ", full -2 " << failedAssignmentAreaFull << ", body -2 " << failedAssignmentAreaBody << ")";
		return true;
	}

	static void cleanupCavityAreaSidecars(int cavityIndex)
	{
		char fileName[160];
		sprintf(fileName, "%scav_tri%dtriangleAreas.txt", conf.rootFile.c_str(), cavityIndex);
		::remove(fileName);
		sprintf(fileName, "%scav_tri%dtriangleAreas.bin", conf.rootFile.c_str(), cavityIndex);
		::remove(fileName);
		sprintf(fileName, "%scav_tri_body%dtriangleAreas.txt", conf.rootFile.c_str(), cavityIndex);
		::remove(fileName);
		sprintf(fileName, "%scav_tri_body%dtriangleAreas.bin", conf.rootFile.c_str(), cavityIndex);
		::remove(fileName);
		sprintf(fileName, "%scav_tri_body%d.atomSurfaceAreas.txt", conf.rootFile.c_str(), cavityIndex);
		::remove(fileName);
	}

	static void cleanupTemporaryCavityFiles(int cavityIndex)
	{
		char cavFile[100];
		char allCavFile[100];
		sprintf(cavFile, "%scav%d.txt", conf.rootFile.c_str(), cavityIndex);
		sprintf(allCavFile, "%sall_cav%d.txt", conf.rootFile.c_str(), cavityIndex);
		::remove(cavFile);
		::remove(allCavFile);
	}

public:
	void buildPocketCavityEnvelope(Surface *surf2, vector<bool> &isPocket,
								   int cav_start, int cav_end, int num_threads, ConfigFile *cf, bool saveEntranceInfo)
{
	for (int i = cav_start; i < cav_end; i += num_threads)
	{
		char mol[100], tri_[100];
		sprintf(mol, "%scav%d.txt", conf.rootFile.c_str(), i);
		sprintf(tri_, "cav_tri%d", i);

		DelPhiShared *dg_temp = new DelPhiShared(conf.maxNumAtoms, conf.domainShrinkage, conf.optimizeGrids,
												 conf.scale, 50, mol, false, false, false);

		// reset seed
		// srand(conf.currentSeed);

		Surface *cs_temp = surfaceFactory().create(cf, dg_temp);
		// relatively big such that non degenerate configurations are avoided
		// due to the high packing of atoms. this has a minimal impact
		// on the estimation of the pocket/surface volume
		cs_temp->setRandDisplacement(0.1);
		cs_temp->setAccurateTriangulationFlag(true);
		cs_temp->setCheckDuplicatedVertices(false);
		cs_temp->build();

		cs_temp->getSurf(&volumes[i], conf.optimizeGrids, false);

		if (conf.tri)
		{
			cout << endl << INFO << "Triangulating enveloped cavity/pockets " << i;
			areas[i] = cs_temp->triangulateSurface(true, true, 0.0, tri_);
			areas2[i] = areas[i];

			vector<bool> externalMask;
			cs_temp->triangulationPointsAreCompletelyOut(surf2, externalMask);

				vector<double> atomAreaFull, atomAreaBody;
				vector<char> atomTouched;
				vector<int> atomOutputIds;
				double failedAssignmentAreaFull = 0.0;
				double failedAssignmentAreaBody = 0.0;
				double filteredBodyArea = 0.0;
				AtomAreaStats statsFull, statsBody;

				computeAtomAreasForMesh(cs_temp, surf2, externalMask,
										atomAreaFull, atomAreaBody,
										failedAssignmentAreaFull, failedAssignmentAreaBody,
										filteredBodyArea, statsFull, statsBody, atomTouched);

				const int numLoadedAtoms = surf2->getNumLoadedAtoms();
				atomOutputIds.assign(numLoadedAtoms, -1);
				for (int ai = 0; ai < numLoadedAtoms; ai++)
					atomOutputIds[ai] = surf2->getLoadedAtomOutputSerial(ai);

			// NanoShaper identifies the set of grid points that are nearest
			// to the given triangulation points and save those points that are completely
			// solvent exposed (i.e. external) in surf2 (the slim surface).
			// Those points represent a point-wise representation of the entrance of the pocket
			if (saveEntranceInfo)
			{
				if (!isPocket[i])
				{
					// cout << endl << ERR << "Cannot find the entrance of a cavity; it must be a pocket.";
					if (!conf.cavAndPockets)
					{
						cout << endl << REMARK << "You have to enable Cavities and pockets flag to do a distinction between a pocket and a cavity";
						exit(-1);
					}
				}
				else
				{
					// it's a pocket, we can identify the entrance;
					// true means that it is part of the entrance
					vector<bool> &results = externalMask;

					// the points which are completely out are those that are represent the entrance
					// together with them we must save the normals;
					// surf2 selectively save those normals and those vertices of this surface
					char oldEntrance[160];
					sprintf(oldEntrance, "%sentrance%d.xyz", conf.rootFile.c_str(), i);
					::remove(oldEntrance);

					char buff[100];
					sprintf(buff, "entrance%d.xyzn", i);
					cs_temp->saveSelectedPoints(results, buff, NULL);

					char triSubset[100];
					sprintf(triSubset, "cav_tri_body%d", i);

					// save triangulation and estimate area removing entrance triangulation points
					vector<bool> results2;
					for (unsigned int i=0; i<results.size(); i++)
						results2.push_back(!results[i]);

					bool revert = false;
					double reducedArea = cs_temp->saveTriSubSet(triSubset, results2, revert);

					areas2[i] = reducedArea;
				}
			}

				writeCombinedAtomAreasFile(i, atomAreaFull, atomAreaBody,
										   atomOutputIds,
										   areas[i], areas2[i],
										   filteredBodyArea,
										   failedAssignmentAreaFull,
										   failedAssignmentAreaBody,
										   atomTouched);
			cleanupCavityAreaSidecars(i);
			}
			delete cs_temp;
			delete dg_temp;
			cleanupTemporaryCavityFiles(i);
		}
	}
	};



#ifndef DELPHI_BIND
/** @mainpage
	@brief
	NanoShaper is a framework to analyze an arbitrary surface; a particular attention is given to molecular surfaces.
	NanoShaper presents several features:
	<ul>
	<li> Loading of a triangulated mesh in OFF, PLY or MSMS format
	<li> Analytical build and triangulation of the Skin surface, the Blobby surface and the Connolly surface
	<li> Colour with in/out info of a 3D grid
	<li> Pockets and cavities detection.
	<li> Easy to plug a new surface: it is sufficient to expand Surface class, no other modifications are required,
	just writing the .h and .cpp. The surface registration is automatically managed by a template based mechanism.
	<li> It can be compiled as a Standalone, python lib or as a DelPhi plug-in.
	</ul>

	History:

	<ul>
	<li> Version 0.3.1:
	<ul>
	<li> Much faster Skin Surface build-up. Acceleration grid 2d also for the mesh ray casting.
	</ul>

	<li> Version 0.3.2:
	<ul>
	<li> bug fix on Connolly Surface.
	</ul>

	<li> Version 0.4:
	<ul>
	<li> Using octrees to reduce memory requirements for triangulation and using an improved multi-threading strategy.
	<li> Atom based Multi-dielectric is supported (DelPhi).
	<li> Stern Layer (DelPhi).
	<li> An efficient, ray-casting based algorithm to convert a set of surface samples to a set of balls whose envelope approximatess the surface.
	<li> Atoms in the neighboroud of a detected cavity are returned.
	</ul>

	<li> Version 0.5:
	<ul>
	<li> "-=" Filtered Operator to support difference map to detect pockets not detectable by cavity detector
	<li> Operative_Mode keyword. Values = normal, difference, membrane
	<li> Introduced experimental Fast Van Der Waals surface. Analytical for FD PDE, not analytical triangulation
	<li> Introduced experimental Coulombic surface
	<li> From DelPhi atinf and atsurf are imported. atsurf is filled  (DelPhi)
	<li> Fast projection flag to switch between fast (approximate) and slow (accurate) SkinSurface bgp projections
	<li> Output for pocket detection in ProShape style added when used with DelPhi
	<li> Improved threads load balancing
	<li> Every vertex has its own reference atoms (OFF+A format defined)
	<li> Faster flood fill algorithm
	<li> Analytical normals computation or approximation by triangulation
	<li> OFF+N, OFF+N+A, MSMS formats
	<li> Vertex normals are computed, analytically if possible.
	</ul>

	<li> Version 0.6:
	<ul>
	<li> Refactoring using metaprogramming to instantiate surfaces.
		To introduce a new surface, just define that class and register it using the registration macro.
		Now surface constructors access directly the configuration file structure. Any new surface thus
		can define each own keywords.
	</ul>

	<li> Version 0.7
	<ul>
	<li> Parallel Marching Cubes: 2 steps MC, cache friendly and trivially parallel
	<li> Dangling vertices now are removed a priori. Dangling vertices are those that were sampled, gave
		an in/out map that is in contrast with another ray for numerical reasons. If this happen these
		vertices, after ray tracing, are identified and removed before marching cubes. In this way
		the exact number of MC vertices is known a priori.
		Another pathological situation that may happen on previous versions is that two vertices are false
		vertices on the same edge (let's say an high frequency) now these vertices are removed.
		The previous version of NanoShaper had this added spurious vertices in memory; however they are not
		present in the mesh, so the previous meshes are still correct.
		Moreover now a slightly faster way to color the grid in ray tracing is defined.
	<li> If a normal cannot be computed analytically due to any kind of problem, than this is approximated
		in a subsequent step such that any file with normals will have all the normals computed. Most
		of them will be analytical (usually less than 1 over 10000 are not analytical).
	<li> High accuracy timing boost::chrono is optionally available
	<li> Bug fix in blobby surface due to surface factory refactoring
	<li> Introduced the ExampleSurface to let the user better understand how to add a surface
	<li> Using templates instead of macros for memory access/write, check bounds and added p=NULL semantic on deletion.
		from now on vectors of primitive types should be allocated, accessed, written, deleted with the given template functions
	<li> Refactoring of main method.
	<li> Introduced pure virtuals init/clear methods to force the Developer to write clean C++ code for the Surfaces
	<li> Optimization on the floodfill. Ping pong to get as much as possible cache friendly.
	<li> Bug fix of ConnollySurface, now Probe_Radius is not ignored.
	</ul>

	<li> Version 0.7.4
	<ul>
	<li> Introduced Max_Atoms_Multi_Grid keyword to control atomsmap max size.
	<li> Introduced an algorithm to detected the entrance of a pocket. The entrance is saved as a set of points in entranceX.xyz
	and entranceX.xyzn where the first file can be visualized in vmd.
	<li> If the Connolly surface is computed a file name called exposed.xyz saves as Carbons in xyz format (VMD) the list of exposed
	atoms. Additionally a file named exposedIndices.txt is the list of 0-based indices of exposed atoms.
	<li> The area of pockets, excluding the mouth is returned. Triangulated files are named cav_triX_body.TRI where is triangulation
	extension and X is the index of the pocket
	<li> Bugfix for MSMS save for blobby

	</ul>

	<li> Version 0.7.5
	<ul>
	<li> Bug fix: Save_Status_map flag was not read. Fixed
	</ul>

	<li> Version 0.7.6
	<ul>
	<li> Bug fix: Idebmap is corrected twice due to cavity detection. Now idebmap
		 is directly written when epsmap is created. It is consistent with fillCavities
		 and no more corrected based on statusmap.
	</ul>

	<li> Version 0.7.7
	<ul>
	<li> Added "membfit" operative mode. This mode allows to compute a score for the protein
	     current orientation that gives an indication on how much is reasonable that the current
		 orientation is feasible with a membrane. The membrane hight is given in the new keyword
		 "Memb_Height" and the score moves the membrane along the z-axis thus as a result also the
		 z value of the lower plane of the membrane is given. The membrane is shifted with a discrete
		 size, the lower the value, the higher the accuracy. This shift value along the z-axis is
		 given by "Memb_Shift". It is assumed that in the input xyzr file the epsilon
		 column is present. The value epsilon = 1 means hydrophobic residue, the value epsilon = 2
		 means charged residue. Based on these values on the solvent exposed atoms a score is returned.

	</ul>

	<li> Version 0.7.8
	<ul>
	<li> Added generation of triangleAreas.txt file here
	<li> HINT: to statically link on Linux:
		 CGAL ->  cmake . -DBUILD_SHARED_LIBS=OFF -DWITH_examples=false -DWITH_CGAL_Qt4=false -DWITH_CGAL_Qt3=false -DWITH_CGAL_ImageIO=false
		 NS -> compile with make verbose, then at link time replace the original with adding "-static", ".so" -> ".a", add "-lpthread"
	<li> Added Output_Folder keyword: if present this keyword, everything will be saved in this folder
	<li> Added Root_FileName: every file will be saved with this root in front. Modified all surface classes
	<li> Bug fix in MeshSurface; computeNormals duplicated entry removed
	<li> By default in normal mode Memory clean up at exit is not done, this increases speed. New keyword "Skip_Mem_CleanUp"
		 that allows to skip final memory clean up.
	</ul>

	</ul>

	<li> Version 0.8
	<ul>
	<li> Updated cgal version, link with TBB for parallel Regular Triangulation Support.
	<li> Atom is a vector of objects not pointers.
	<li> Added API support for guest codes.
	<li> Added support from API to match arbitrary external grids to allow PB solution.
	<li> Added from API retrieval of intersections and normals in the grid
	</ul>

	</ul>

	<li> Version 0.8.1
	<ul>
	<li> Now it supports CGAL 5.6.2
	<li> Removed the need for CGAL patching, now it uses the _with_info feature of CGAL Vertex class.
	<li> Dockerfile updated accordingly, also using Rockylinux
	</ul>

	</ul>

	@author Sergio Decherchi
	@date 04-10-2024
	@version 0.8.1
	*/

#ifndef NS_TEST

int main(int argc, char *argv[])
{
	#ifdef  DBGMEM_CRT
	// _crtBreakAlloc = 16222;
	_CrtSetDbgFlag (_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	#ifdef BOOST_ENABLED
	initialPath = boost::filesystem::initial_path();
	#endif

	int numargs = argc-1;
	char confFile[BUFLEN];

	if (!numargs)
		strcpy(confFile, "surfaceConfiguration.prm");
	else
		strcpy(confFile,argv[1]);

	if (numargs > 1)
		cout << endl << WARN << "Ignoring additional non required parameters";

	// check configuration consistency, init error stream, get configuration
	ConfigFile *cf = init(string(confFile));

	// operative modes are the NanoShaper lib utilization modes.

	int num_threads = conf.numThreads;

	#if defined(DEBUG_SURFACE) || !defined(MULTITHREADING)
	num_threads = 1;
	#endif
	#if defined(ENABLE_CGAL) && defined(CGAL_LINKED_WITH_TBB)
	if (num_threads <= 0)
	{
		num_threads = (int)std::thread::hardware_concurrency();
		cout << endl << INFO << "Detected " << num_threads << " logical cores";
		cout << endl << INFO << "Setting " << num_threads << " threads";
	}
	else
	{
		cout << endl << INFO << "User selected num threads " << num_threads;
	}
	initTBB (num_threads);
	#endif

	conf.numThreads = num_threads;


	if (conf.buildStatus)
		cout << endl << INFO << "Status map building is enabled... ";
	else
		cout << endl << INFO << "Status map building is disabled... ";

	if (!conf.optimizeGrids)
		cout << endl << INFO << "Optimized grids are disabled... ";
	else
		cout << endl << INFO << "Optimized grids are enabled... ";

	// just build the surface
	if (!conf.operativeMode.compare("normal"))
	{
		// Set up DelPhi-like environment
		DelPhiShared *dg = new DelPhiShared(conf.maxNumAtoms, conf.domainShrinkage, conf.optimizeGrids,
											conf.scale, conf.perfill, conf.molFile, conf.buildEpsmaps,
											conf.buildStatus, conf.multi_diel);

		// Get surface
		Surface *surf = surfaceFactory().create(cf, dg);

		normalMode(surf, dg);

		#if defined(USE_VIS_TOOLS)

		surf->initVisualization (argc, argv, &vis);

		#ifndef NO_OPENGL
		vis.editing.pause = 0;
		/*
		glutReshapeFunc (surf->Reshape);
		glutIdleFunc (surf->Display);
		glutDisplayFunc (surf->Display);
		glutKeyboardFunc (surf->keybordFunction);
		glutMouseFunc (surf->mouseFunction);
		glutMotionFunc (surf->motionFunction);
		glutMainLoop ();
		*/
		#endif // NO_OPENGL

		#endif // USE_VIS_TOOLS

		bool skipCleanUp = cf->read<bool>("Skip_Mem_CleanUp", true);

		if (!skipCleanUp)
		{
			auto chrono_start = chrono::high_resolution_clock::now();

			cout << endl << INFO << "Cleaning memory...";
			cout.flush();

			delete surf;
			delete dg;
			cout << "ok!";

			auto chrono_end = chrono::high_resolution_clock::now();

			chrono::duration<double> cleaning_time = chrono_end - chrono_start;
			cout << endl << INFO << "Cleaning time... ";
			printf ("%.4e [s]", cleaning_time.count());
		}
		else
		{
			cout << endl << INFO << "Not cleaning memory...";
		}
	}
	// detect pockets
	else if (!conf.operativeMode.compare("pockets"))
	{
		pocketMode(false, cf);
	}
	// just build the surface
	else if (!conf.operativeMode.compare("membfit"))
	{
		// this must be always true
		conf.multi_diel = true;
		conf.buildEpsmaps = true;
		// Set up DelPhi-like environment
		DelPhiShared *dg = new DelPhiShared(conf.maxNumAtoms, conf.domainShrinkage, conf.optimizeGrids,
											conf.scale, conf.perfill, conf.molFile, conf.buildEpsmaps, conf.buildStatus, conf.multi_diel);
		// Get surface
		Surface *surf = surfaceFactory().create(cf, dg);
		membfitMode(surf, dg);

		auto chrono_start = chrono::high_resolution_clock::now();

		cout << endl << INFO << "Cleaning memory...";
		cout.flush();

		delete surf;
		delete dg;
		cout << "ok!";

		auto chrono_end = chrono::high_resolution_clock::now();

		chrono::duration<double> cleaning_time = chrono_end - chrono_start;
		cout << endl << INFO << "Cleaning time... ";
		printf ("%.4e [s]", cleaning_time.count());
	}
	else
	{
		cout << endl << INFO << "Unknown operative mode";
		cout << endl;
		return -1;
	}

	#if !defined(AVOID_MEM_CHECKS)
	double current_mem_in_MB, peak_mem_in_MB;
	getMemSpace (current_mem_in_MB, peak_mem_in_MB);
	cout << endl << INFO << "Peak memory is " << peak_mem_in_MB << " MB";
	#endif

	cite();

	dispose(cf);

	cout << endl << endl;

	// in any case reset path to the original
    #ifdef BOOST_ENABLED
	boost::filesystem::current_path(initialPath);
	#endif

	// Memory leak detection 
	#ifdef DBGMEM_CRT
		_CrtDumpMemoryLeaks();
	#endif

	return 0;
}

#else // NS_TEST

int main(int argc, char *argv[])
{
	SkinSurface *sf = new SkinSurface();
	double votes[8],isolevel = 0.5;
	int cube;
	
	votes[0] = 0;
	votes[1] = 0;
	votes[2] = 0;
	votes[3] = 1;
	votes[4] = 0;
	votes[5] = 0;
	votes[6] = 1;
	votes[7] = 0;
	cube = sf->classifyCube(votes,isolevel);
	cout << endl << "case 3 " << cube << endl;

	for (int i=0; i<15; i+=3)
		printf("\t %d %d %d\n", triTable[cube][i],triTable[cube][i+1],triTable[cube][i+2]);
	
	delete sf;
	getchar();
}

#endif // NS_TEST
#endif // DELPHI_BIND


#ifdef DELPHI_BIND
	
	/** This is the DelPhi Fortran-C binding. Returns the number of bgp */
	// check visual, borland c++, xlc, watcom, for windows
	#if  (defined _WIN32) || (defined __WIN32__) || (defined __TOS_WIN__) || (defined __WINDOWS__)
			extern "C" __declspec(dllexport)
	// linux,unix,apple,android
	#elif (defined __linux) || (defined __unix) || (defined macintosh) || (defined Macintosh) || (defined __APPLE__ && defined __MACH__) || (defined __ANDROID__)
			extern "C"
	#endif


void epsmakemodule(double xmin,double ymin,double zmin,
				   double xmax,double ymax,double zmax,
				   double c1,double c2,double c3, double rmax,double perf,
				   int *i_epsmap,int igrid,double scale,double *i_scspos,
				   double *i_scsnor,bool *i_idebmap,int *i_ibgp,int inside,int *ibnum,int natom,
				   double *xn1,double *rad, int *d,int maxbgp,double probe,double exrad,char *atinf,int *atsurf)
// epsmap: per wise dielectric value
// i_scspos: projected points on the n cells which are cut from the surface (unset, in this version)
// i_scsnor: surface normals
// i_idebmap: inside/outside map on inflated surface (screening out ions)
// xn1: atoms' coords
// probe: size of ball-like probe (water molecule)
{
	// check configuration consistency, init error stream and get configuration
	ConfigFile *cf = init(string("surfaceConfiguration.prm"));

	// override settings coming from conf with that coming from DelPhi
	if (cf->keyExists("Grid_scale"))
		cf->remove("Grid_scale");
	cf->add<double>("Grid_scale",scale);

	if (cf->keyExists("Grid_perfil"))
		cf->remove("Grid_perfil");
	cf->add<double>("Grid_perfil",perf);

	conf.scale = scale;
	conf.perfill = perf;

	ConfigFile *cf2;
	// get data from custom configuration file
	try
	{
		cf2 = new ConfigFile("custom.prm");
	}
	catch (...)
	{
		cout << endl << ERR << "Cannot read custom.prm";
		cout.flush();
		exit(-1);
	}
	string mode = cf2->read<string>("Surface", "ses");
	delete cf2;

	if (cf->keyExists("Surface"))
		cf->remove("Surface");

	cf->add<string>("Surface",mode);

	if (conf.operativeMode == "normal")
	{
		cout << endl << INFO << "Binding with DelPhi...";

		// Set up delphi environment
		DelPhiShared *dg = new DelPhiShared();
		dg->DelPhiBinding(	xmin,ymin,zmin,xmax,ymax,zmax,
							c1,c2,c3,rmax,perf,
							i_epsmap,igrid,scale,
							i_scspos,i_scsnor,
							i_idebmap,i_ibgp,
							maxbgp,conf.buildStatus,atsurf);

		// load atoms info if present
		dg->loadAtoms(natom,xn1,rad,NULL,NULL,atinf);

		// here only means populate epsilon map
		dg->buildEpsmap(true);

		cout << endl << INFO << "DelPhi grid is " << igrid;

		Surface *surf = surfaceFactory().create(cf, dg);
		// projection of boundary grid points
		surf->setProjBGP(true);
		surf->setInsideCode(inside);

		if (exrad > 0)
		{
			cout << endl << INFO << "Setting Stern Layer -> " << exrad << " Angstrom";
			surf->setSternLayer(exrad);
		}

		// normal protocol
		normalMode(surf, dg);

		// avoid to destroy things that now belong to DelPhi
		dg->finalizeBinding(ibnum);

		cout << endl << INFO << "Cleaning memory...";
		cout.flush();

		delete surf;
		delete dg;
		
		dispose(cf);

		cout << "ok!";
		cout << endl << endl;
		return;
	}
	// In this mode use the difference on grid operator to infer pockets and cavities (if requested)
	else if (!conf.operativeMode.compare("pockets"))
	{
		DelPhiShared::parseAtomInfo(conf.maxNumAtoms,atinf,conf.molFile,natom,xn1,rad);

		pocketMode(true, cf);

		dispose(cf);

		cout << endl << INFO << "Brute force exiting to avoid DelPhi solver";
		cout << endl << INFO << "Closing " << PROGNAME << "\n";
		exit(-1);
	}
	else
	{
		cout << endl << INFO << "Unknown operative mode: " << conf.operativeMode;
		exit(-1);
	}

	cite();
	dispose(cf);
	cout << endl << endl;
}

#endif // DELPHI_BIND

// init streams, check configuration file for errors and read variables
ConfigFile *init(string confFile)
{
	cout << endl << endl << INFO << "Starting " << PROGNAME << " " << VERSION;

	ConfigFile *cf;

	// get data from configuration file
	try
	{
		cf = new ConfigFile(confFile.c_str());
	}
	catch (...)
	{
		cout << endl << ERR << "Cannot read " << confFile << endl;
		exit(-1);
	}

	// if specific path is given, it creates the output folder if absent
	// and change the current directory
	string outputFolder = cf->read<string>("Output_Folder", "");
	if (outputFolder.size() != 0)
	{
		#ifdef BOOST_ENABLED
		boost::filesystem::path fs(outputFolder);
		if (!boost::filesystem::exists(fs))
			boost::filesystem::create_directory(fs);
		boost::filesystem::current_path(fs);
		boost::filesystem::current_path(fs);
		boost::filesystem::current_path(fs);
		#endif
	}

	string rootFile = cf->read<string>("Root_FileName", "");
	string tt = rootFile + "stderror.txt";
	errorStream = new fstream(tt.c_str(),fstream::out);

	conf.buildEpsmaps = cf->read<bool>("Build_epsilon_maps", false);
	conf.saveEpsmaps = cf->read<bool>("Save_eps_maps", false);
		
	conf.saveStatusMap = cf->read<bool>("Save_Status_map", false);

	conf.buildStatus = cf->read<bool>("Build_status_map", false);
	conf.saveIdebmap = cf->read<bool>("Save_ideb_map", false);
	conf.fillCavities = cf->read<bool>("Cavity_Detection_Filling", false);

	conf.projBGP = cf->read<bool>("Project_boundary_grid_points", false);

	conf.accTri = cf->read<bool>("Accurate_Triangulation", false);
	conf.tri = cf->read<bool>("Triangulation", false);

	conf.operativeMode = cf->read<string>("Operative_Mode", "normal");

	bool dbg = cf->read<bool>("Debug_Internals", false);

	conf.patchBased = cf->read<bool>("Patch_Based_Algorithm", true);
	conf.analyticalRayVsTorusIntersection = cf->read<bool>("Analytical_Ray_Vs_Torus_Intersection", true);
	conf.collectFaceIntersections = cf->read<bool>("Collect_face_intersections", false);
	conf.forceSerialBuild = cf->read<bool>("Force_Serial_Build", false);
	conf.parallelBuildupHaloThickness = cf->read<double>("Parallel_Buildup_Halo_Thickness", 0.0);
	conf.maxNumAtoms = cf->read<int>("Max_Num_Atoms", -1);
	conf.domainShrinkage = cf->read<double>("Domain_Shrinkage", 0);
	conf.optimizeGrids = cf->read<bool>("Optimize_Grids", true);
	
	if (dbg)
		internals = new fstream("internals.txt", fstream::out);
	else
		internals = NULL;

	if (cf != NULL)
	{
		if (!conf.buildEpsmaps && conf.saveEpsmaps)
		{
			cout << endl << ERR << "Asked to save epsmap without builiding it";
			cout << endl << REMARK << "Please set Build_epsilon_maps = true";
			cout << endl;
			exit(-1);
		}
		if (!conf.buildEpsmaps && conf.projBGP)
		{
			cout << endl << ERR << "Cannot project boundary grid points without an epsilon map.";
			cout << endl << REMARK << "Please set Build_epsilon_maps = true";
			cout << endl;
			exit(-1);
		}
		if (!conf.accTri && !conf.buildStatus && conf.tri)
		{
			// status map is needed to deduced in/out vertices
			cout << endl << ERR << "If non analytical triangulation is enabled status map is needed.";
			cout << endl << REMARK << "Please set Build_status_map = true";
			cout << endl;
			exit(-1);
		}
		if (conf.fillCavities && !conf.buildStatus)
		{
			// status map is needed to search cavities
			cout << endl << ERR << "If cavity detection is enabled status map is needed.";
			cout << endl << REMARK << "Please set Build_status_map = true";
			cout << endl;
			exit(-1);
		}
		if (conf.saveIdebmap && !conf.buildEpsmaps)
		{
			cout << endl << ERR << "Idebmap is computed only if epsilon map is enabled";
			cout << endl << REMARK << "Please set Build_epsilon_maps = true";
			cout << endl;
			exit(-1);
		}
		if (!conf.operativeMode.compare("pockets") && !conf.buildStatus)
		{
			cout << endl << WARN << "Cannot do pocket detection without status map";
			cout << endl << REMARK << "Please set Build_status_map = true";
			cout << endl;
			exit(-1);
		}
		if (!conf.operativeMode.compare("pockets"))
		{
			cout << endl << WARN << "Status map space is not optimised in pocket mode because of slower runs";
			cout << endl;

			conf.optimizeGrids = false;
		}
		string surf_build_type = cf->read<string>("Surface", "ses");

		if (conf.patchBased &&
			(surf_build_type.compare("ses") && surf_build_type.compare("skin") && surf_build_type.compare("mesh")))
		{
			cout << endl << ERR << "This surface type does not support patch based ray-tracing";
			cout << endl << REMARK << "Please set Patch_Based_Algorithm = false";
			cout << endl;
			exit(-1);
		}
	}

	conf.cavVol = cf->read<double>("Conditional_Volume_Filling_Value", 11.4);
	conf.numMol = cf->read<int>("Num_Wat_Pocket", 2);
	
	// grid (DelPhi) params
	conf.scale = cf->read<double>("Grid_scale", 2.0);
	conf.perfill = cf->read<double>("Grid_perfil", 80.0);
	conf.molFile = cf->read<string>("XYZR_FileName", "temp.xyzr");
	conf.multi_diel = cf->read<bool>("Multi_Dielectric", false);

	// tri
	conf.smoothing = cf->read<bool>("Smooth_Mesh", false);
	conf.tri2balls = cf->read<bool>("Tri2Balls", false);
	
	// save data
	conf.saveEpsmaps = cf->read<bool>("Save_eps_maps", false);
	conf.saveBgps = cf->read<bool>("Save_bgps", false);
	conf.saveCavities = cf->read<bool>("Save_Cavities", false);

	// globals
	conf.sysName = cf->read<string>("Sys_Name","mol");
	conf.numThreads = cf->read<int>("Number_thread", -1);
	conf.printAvailSurf = cf->read<bool>("Print_Available_Surfaces", false);
	conf.currentSeed = cf->read<int>("Seed", 1);

	// pocket detection
	conf.cavAndPockets = cf->read<bool>("Pockets_And_Cavities", true);
	conf.linkPockets = cf->read<bool>("Link_Pockets", false);
	conf.pocketRadiusBig = cf->read<double>("Pocket_Radius_Big", 3.0);
	conf.pocketRadiusSmall = cf->read<double>("Pocket_Radius_Small", 1.4);
	conf.pocketRadiusLink = cf->read<double>("Pocket_Radius_Link", 1.0);

	// membfit parameters
	conf.membHeight = cf->read<double>("Memb_Height", 25.0);
	conf.membShift = cf->read<double>("Memb_Shift", 1.0);
			
	// root file name parseing
	conf.rootFile = cf->read<string>("Root_FileName","");

	return cf;
}


void dispose(ConfigFile *cf)
{
	delete cf;
	errorStream->close();
	delete errorStream;

	if (internals != NULL)
	{
		internals->close();
		delete internals;
	}
}


void stopDebug()
{
	if (internals != NULL)
	{
		internals->close();
		internals = NULL;
	}
}


void restartDebug()
{
	if (internals == NULL)
		internals = new	fstream("internals.txt",fstream::app);	
}


void cite()
{
	cout << endl;
	cout << endl << INFO << "If you use NanoShaper please cite these works:";
	cout << endl << CITE << "\tS. Decherchi, W. Rocchia, \"A general and Robust Ray-Casting-Based Algorithm for Triangulating Surfaces at the Nanoscale\"; PlosOne";
	cout << endl << CITE << "\tlink: http://www.plosone.org/article/metrics/info%3Adoi%2F10.1371%2Fjournal.pone.0059744";
	cout << endl << CITE << "\tS. Decherchi, A. Spitaleri, J. Stone, W. Rocchia, \"NanoShaper-VMD interface: computing and visualizing surfaces, pockets and channels in molecular systems\"; Bioinformatics";
	cout << endl << CITE << "\t2019 Apr 1;35(7):1241-1243. doi: 10.1093/bioinformatics/bty761. PMID: 30169777; PMCID: PMC6449750.";
	cout << endl;
}


/** the set of operations in the usual mode of usage. This function is not responsible for
Surface or grid memory. The caller is the responsible.*/
// surface computation (?)
void normalMode(Surface *surf, DelPhiShared *dg)
{
	if (conf.printAvailSurf)
		surfaceFactory().print();

	auto chrono_total_time_start = chrono::high_resolution_clock::now();

	auto chrono_start = chrono::high_resolution_clock::now();

	double surf_area, surf_volume;

	char refName[100];
	strcpy(refName,conf.sysName.c_str());

	// cout << endl << "Compute normals surf " << surf->getComputeNormals();

	// Pre-process surface
	bool outsurf = surf->build();

	if (!outsurf)
	{	
		cout << endl << ERR << "Surface construction failed!" << endl;
		exit(-1);
	}

	// Build DelPhi stuff
	surf->getSurf(&surf_volume, conf.optimizeGrids, conf.fillCavities, conf.cavVol);

	auto chrono_end = chrono::high_resolution_clock::now();

	chrono::duration<double> surface_computation_time = chrono_end - chrono_start;
	cout << endl << INFO << "Surface computation time... ";
	printf ("%.4e [s]", surface_computation_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	double current_mem_in_MB, peak_mem_in_MB;
	getMemSpace (current_mem_in_MB, peak_mem_in_MB);
	cout << endl << INFO << "Memory required after call of getSurf() in normalMode() is " << current_mem_in_MB << " MB";
	#endif

	cout << endl << INFO << "Estimated volume " << setprecision(10) << surf_volume << " [A^3]";

	if (internals != NULL)
		(*internals) << endl << "volume " << surf_volume;

	if (conf.tri)
	{
		cout << endl << INFO << "Triangulating Surface...";

		auto chrono_start = chrono::high_resolution_clock::now();

		if (conf.smoothing)
			surf_area = surf->triangulateSurface(false, false);
		else
			surf_area = surf->triangulateSurface(true, true);

		if (internals != NULL)
		{
			(*internals) << endl << "area " << surf_area;
			(*internals) << endl << "nv " << surf->getNumVertices();
			(*internals) << endl << "nt " << surf->getNumTriangles();
		}
		if (conf.smoothing)
		{
			cout << endl << INFO << "Smoothing surface...";
			surf->smoothSurface(true, true);
		}
		cout << " ...ok!";

		auto chrono_end = chrono::high_resolution_clock::now();

		chrono::duration<double> triangulation_time = chrono_end - chrono_start;

		cout << endl << INFO << "Total triangulation time ";
		printf ("%.4e [s]", triangulation_time.count());

		#if !defined(AVOID_MEM_CHECKS)
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after triangulation is " << current_mem_in_MB << " MB";
		#endif
	}

	if (conf.tri2balls)
	{
		auto chrono_start = chrono::high_resolution_clock::now();

		cout << endl << INFO << "Converting triangulation to balls...";
		surf->tri2Balls();
		cout << "ok!";

		auto chrono_end = chrono::high_resolution_clock::now();

		chrono::duration<double> tri2balls_time = chrono_end - chrono_start;

		cout << endl << INFO << "tri2balls() time ";
		printf ("%.4e [s]", tri2balls_time.count());

		#if !defined(AVOID_MEM_CHECKS)
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after tri2balls() is " << current_mem_in_MB << " MB";
		#endif
	}

	chrono_start = chrono::high_resolution_clock::now();

	if (conf.saveEpsmaps)
	{
		cout << endl << INFO << "Saving epsmaps...";
		// Save epsmap
		dg->saveEpsMaps(refName);
		cout << "ok!";
	}

	if (conf.saveBgps)
	{
		cout << endl << INFO << "Saving bgpmap...";
		dg->saveBGP(refName);
		cout << "ok!";
	}

	if (conf.saveStatusMap)
	{
		cout << endl << INFO << "Saving statusmap and cavities...";
		dg->saveStatus(refName);
		cout << "ok!";
	}

	if (conf.saveCavities)
	{
		cout << endl << INFO << "Saving cavities...";
		dg->saveCavities(false);
		cout << "ok!";
	}
		
	if (conf.saveIdebmap)
	{
		cout << endl << INFO << "Saving idebmap...";
		dg->saveIdebMap(refName);
		cout << "ok!";
	}

	chrono_end = chrono::high_resolution_clock::now();

	chrono::duration<double> file_saving_time = chrono_end - chrono_start;

	cout << endl << INFO << "Files' saving time (in normalMode()): ";
	printf ("%.4e [s]", file_saving_time.count());

	auto chrono_total_time_end = chrono::high_resolution_clock::now();

	chrono::duration<double> total_computation_time = chrono_total_time_end - chrono_total_time_start;
	cout << endl << INFO << "Surface + triangulation (+ smoothing + files outputting) time: ";
	printf ("%.4e [s]", total_computation_time.count());
}


void pocketMode(bool hasAtomInfo, ConfigFile *cf)
{
	bool localEpsMap = false;
	bool localStatusMap = true;
	bool localMulti = false;

	// only debug pockets/cavities
	if (conf.debug) stopDebug();

	auto chrono_start = chrono::high_resolution_clock::now();

	Surface *surf1, *surf2, *surf3;

	double surf_area[] = {0, 0, 0};
	double surf_volume[] = {0, 0, 0};

	char mol[100] = "temp.txt";
	// char refName[100] = "mol";


	// For now, the following forces the cofde to use the full flat status-map grids until conf.optimizeGrids is re-updated
	bool old_optimizeGrids = conf.optimizeGrids;
	// conf.optimizeGrids = true;

	// Set up Surface 1 (fat probe)
	cout << endl;
	cout << endl << INFO << "Step 1 -> fat probe";

	DelPhiShared *dg1 = new DelPhiShared(conf.maxNumAtoms, conf.domainShrinkage, conf.optimizeGrids,
										 conf.scale, conf.perfill, conf.molFile, localEpsMap,
										 localStatusMap, localMulti, hasAtomInfo);

    hasAtomInfo = dg1->isAvailableAtomInfo;
	int natom = dg1->getNumAtoms();
	cf->remove(string("Surface"));
	cf->add<string>("Surface","ses");
	surf1 = surfaceFactory().create(cf, dg1);
	surf1->setProjBGP(false);
	surf1->setProbeRadius(conf.pocketRadiusBig);
	surf1->setKeepWellShapedCavities(false);
	surf1->setInsideCode(5);

	// Pre-process surface
	bool outsurf = surf1->build();

	if (!outsurf)
	{
		cout << endl << ERR << "Surface 1 construction failed!";
		exit(-1);
	}

	// fat connolly cancel each cavity (anyway they are smaller than they should be at the end)
	surf1->getSurf(&surf_volume[0], conf.optimizeGrids, true, INFINITY);

	//////////////////////////////////////////////////////////////////////////////////////////
	// the subsequent surfaces don't need to read atom info

	// Set up Surface 2 (regular probe 1.4)
	// Set up DelPhi grid  2 
	cout << endl;
	cout << endl << INFO << "Step 2 -> small probe";

	DelPhiShared *dg2 = new DelPhiShared(conf.maxNumAtoms, conf.domainShrinkage, conf.optimizeGrids,
										 conf.scale, conf.perfill, conf.molFile, localEpsMap, localStatusMap, localMulti, false);

	surf2 = surfaceFactory().create(cf, dg2);
	surf2->setProjBGP(false);
	surf2->setProbeRadius(conf.pocketRadiusSmall);
	surf2->setKeepWellShapedCavities(false);
	surf2->setInsideCode(10);
	
	// Pre-process surface
	outsurf = surf2->build();

	if (!outsurf)
	{
		cout << endl << ERR << "Surface 2 construction failed!";
		exit(-1);
	}

	// if cav and pockets together -> do not perform cavity detection (keep all)
	// if only pockets then perform cavity detection and remove all cavities
	// by default keep both cavities and pockets
	surf2->getSurf(&surf_volume[1], conf.optimizeGrids, !conf.cavAndPockets, INFINITY);

	// if triangulation is enabled we triangulate the small probe surface, get surface and area
	if (conf.tri)
	{
		surf_area[1] = surf2->triangulateSurface(false, false);
		surf2->smoothSurface(true, true);
	}

	// to check percolation
	/*
	surf1->difference(surf2);
	surf_area[0] = surf1->triangulateSurface(true, true, 0.0, "diffmap.off", true);
	exit(-1);
	*/

	DelPhiShared *dg3 = NULL;

	if (conf.linkPockets)
	{
		// Set up Surface 3 (accessibility probe)
		dg3 = new DelPhiShared(conf.maxNumAtoms, conf.domainShrinkage, conf.optimizeGrids,
							   conf.scale, conf.perfill, conf.molFile, localEpsMap,
							   localStatusMap, localMulti, false);

		surf3 = surfaceFactory().create(cf, dg3);
		surf3->setProjBGP(false);
		surf3->setProbeRadius(conf.pocketRadiusSmall);
		surf3->setKeepWellShapedCavities(false);
		surf3->setInsideCode(15);

		// Pre-process surface
		outsurf = surf3->build();

		if (!outsurf)
		{
			cout << endl << ERR << "Surface 3 construction failed!";
			exit(-1);
		}
		// keep original surface
		surf3->getSurf(&surf_volume[2], conf.optimizeGrids, false);
	}

	cout << endl;
	cout << endl << INFO << "Step 3 -> differential map";

	cout << endl << INFO << "Building pockets by difference map...";

	(*surf1) -= (*surf2);

	/////////////// recover split cavities links ///////////////
	if (conf.linkPockets)
	{
		int nr = 1;

		auto chrono2_start = chrono::high_resolution_clock::now();

		cout << endl;
		cout << endl << INFO << "Step 3.1 -> linking";

		cout << endl << INFO << "Linking cavities/pockets...";
		cout.flush();

		while (nr != 0)
		{
			// check cavities links, use the link status map as reference map
			// to check accessibility
			if (!conf.optimizeGrids)
				nr = surf1->linkCavities(dg2->status, dg3->status);
			else
				nr = surf1->linkCavities(dg2->bilevel_status, dg3->bilevel_status);

			cout << endl << INFO << "Merged " << nr << " cavities";
			cout.flush();
		}
		auto chrono2_end = chrono::high_resolution_clock::now();

		chrono::duration<double> linking_cavities_time = chrono2_end - chrono2_start;
		cout << endl << INFO << "Diff. Step 4 ";
		printf ("%.4e [s]", linking_cavities_time.count());

		#if !defined(AVOID_MEM_CHECKS)
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after diff. step 4 is " << current_mem_in_MB << " MB";
		#endif

		delete surf3;
		delete dg3;
	}
	///////////////////////////////////////////////////////////////////

	cout << endl;
	cout << endl << INFO << "Step 4 -> filtering, envelope building";

	///////////////////// volume filter ///////////////////////////////
	surf1->fillCavities(11.4*conf.numMol,false);

	if (!conf.optimizeGrids)
		surf1->getCavitiesAtoms();
	else
		surf1->getCavitiesAtomsWithBilevelStatusMap();

	cout << endl << INFO << "Saving cavities info..";

	if (hasAtomInfo)
	{
		// save cavities in ProShape format
		dg1->saveCavities2(true, conf.sysName);
		// Keep atom serials/residue helper files synchronized with .pocket output.
		dg1->saveCavities(true);
	}
	else
	{
		// save cavities in NanoShaper format
		bool saveOnlyNonFilled = true;
		dg1->saveCavities(saveOnlyNonFilled);
	}

	///////////////////// mark support atoms //////////////////////////
	// perform again filtering only to mark the support atoms of each cavity
	if (!conf.optimizeGrids)
		surf1->filterCavities(true);
	else
		surf1->filterCavitiesWithBilevelStatusMap(true);

	// write cavities as atoms on files to subsequent triangulation
	int nc;
	if (!conf.optimizeGrids)
		nc = dg1->cavitiesToAtoms(1.4);
	else
		nc = dg1->cavitiesToAtomsWithBilevelStatusMap(1.4);

	vector<bool> isPocket;

	if (conf.cavAndPockets)
	{
		cout << endl;
		// do cavity detection to recover which is cavity in the reference
		cout << endl << INFO << "Recovering cavity/pocket distinction..";
		cout.flush();

		if (!conf.optimizeGrids)
			surf2->getCavities();
		else
			surf2->getCavitiesWithBilevelStatusMap();

		// mark which are pockets and which are cavities
		if (!dg2->optimizeGrids)
			dg1->markPockets(dg2->status, isPocket);
		else
			dg1->markPockets(dg2->bilevel_status, isPocket);
	}

	cout << endl << INFO << "Build the envelope of each cavity/pocket";
	cout.flush();

	cf->remove("Surface");
	cf->add<string>("Surface","skin");

	// TODO make it an option
	bool saveEntranceInfo = true;

	pocketWrapper *pw = new pocketWrapper();

	// areas of pockets including the pocket closure
	pw->areas   = new double [nc];
	// areas of pockets upon removal of the pocket closure surface
	pw->areas2  = new double [nc];
	pw->volumes = new double [nc];

	if (conf.tri)
		surf2->prepareAtomsMap();

	#if !defined(MULTITHREADED_POCKET_LOOP)
	pw->buildPocketCavityEnvelope(surf2, isPocket, 0, nc, 1, cf, saveEntranceInfo);
	#else // MULTITHREADED_POCKET_LOOP

	int num_old_threads = conf.numThreads;
	#if defined(AVOID_INNER_POCKET_LOOP_MULTITHREADING)
	conf.numThreads = 1;
	#endif

	conf.parallelPocketLoop = true;

	#ifdef ENABLE_BOOST_THREADS
	boost::thread_group thdGroup;
	#endif
	for (int thread_id=0; thread_id < num_old_threads; thread_id++)
	{
		#ifdef ENABLE_BOOST_THREADSoptimizeGrids
		thdGroup.create_thread(boost::bind(&pocketWrapper::buildPocketCavityEnvelope,pw,
										   surf2,isPocket,thread_id,nc,num_old_threads,cf,saveEntranceInfo));
		#else
		pw->buildPocketCavityEnvelope(surf2,isPocket,thread_id,nc,num_old_threads,cf,saveEntranceInfo);
		#endif
	}
	#ifdef ENABLE_BOOST_THREADS
	thdGroup.join_all();
	#endif

	conf.parallelPocketLoop = false;

	conf.numThreads = num_old_threads;

	conf.optimizeGrids = old_optimizeGrids;

	#endif // MULTITHREADED_POCKET_LOOP

	if (conf.tri)
		surf2->clearAtomsMap();

	auto chrono_end = chrono::high_resolution_clock::now();

	if (conf.debug) restartDebug();

	if (hasAtomInfo)
	{
		char name2[100];
		sprintf(name2, "%s.info", conf.sysName.c_str());
		FILE *fp = fopen(name2,"w");
		fprintf(fp, "\n Protein             : ");
		fprintf(fp, "\n Probe radius        : %.2f", conf.pocketRadiusSmall);
		fprintf(fp, "\n Number of atoms     : %d", natom);
		fprintf(fp, "\n Total Surface Area  : %.9e", surf_area[1]);
		fprintf(fp, "\n Total Volume        : %.9e", surf_volume[1]);
		fprintf(fp, "\n\n");
		fprintf(fp, "\n Pockets :");
		fprintf(fp, "\n");
		fprintf(fp, "\n Id\t\tN_mth\tSurface\t\t\tVolume\t\tBody\n");
		
		double sumA = 0, sumV = 0;

		for (int i=0; i<nc; i++)
		{
			if (conf.tri)
			{
				if (conf.cavAndPockets)
				{
					if (isPocket[i] == true)
						fprintf(fp, " %d\t\t>0\t\t%.4f\t\t%.4f\t\t%.4f\n", i+1,pw->areas[i],pw->volumes[i], pw->areas2[i]);
					else
						fprintf(fp, " %d\t\t=0\t\t%.4f\t\t%.4f\t\t%.4f\n", i+1,pw->areas[i],pw->volumes[i], pw->areas2[i]);
				}
				else
					fprintf(fp, " %d\t\t>0\t\t%.4f\t\t%.4f\t\t%.4f\n", i+1,pw->areas[i],pw->volumes[i], pw->areas2[i]);

				sumA += pw->areas[i];
				sumV += pw->volumes[i];
			}
			else
			{
				if (conf.cavAndPockets)
				{
					if (isPocket[i] == true)
						fprintf(fp, "\t%d\t>0\t\t%.4f\t\t%.4f\n", i+1,0.0,pw->volumes[i]);
					else
						fprintf(fp, "\t%d\t0\t\t%.4f\t\t%.4f\n", i+1,0.0,pw->volumes[i]);
				}
				else
					fprintf(fp, " %d\t\t>0\t\t%.4f\t\t%.4f\n", i+1,pw->areas[i],pw->volumes[i]);

				sumV += pw->volumes[i];
			}
		}
		fprintf(fp, "\nTot\t\t%.4f\t\t%.4f\n", sumA, sumV);
		fclose(fp);
	}
	else
	{
		cout << endl << INFO;
		cout << endl << INFO << "------------------------------------";
		cout << endl << INFO << "     Pocket Detection Summary       ";
		cout << endl << INFO << "------------------------------------";
		cout << endl << INFO;
		cout << endl << INFO << "Detected a total of " << nc << " pockets/cavities having at least the volume of " << conf.numMol << " water molecules" ;
			
		for (int i=0; i<nc; i++)
		{
			if (conf.cavAndPockets)
			{
				if (conf.tri)
				{
					if (isPocket[i])
					{
						#if !defined(MULTITHREADED_POCKET_LOOP)
						cout << endl << INFO << "Pocket " << i << " vol " << pw->volumes[i] << " area " << pw->areas[i] << " body area "<< pw->areas2[i];
						if (internals != NULL)
						{
							*internals << endl << "pocket_vol " << pw->volumes[i];
							*internals << endl << "pocket_area " << pw->areas[i] << " pocket_body_area " << pw->areas2[i];
						}
						#else
						cout << endl << INFO << "Pocket " << i << " vol " << pw->volumes[i] << " area " << pw->areas[i] << " body area "<< pw->areas2[i];
						if (internals != NULL)
						{
							*internals << endl << "pocket_vol " << pw->volumes[i];
							*internals << endl << "pocket_area " << pw->areas[i] << " pocket_body_area " << pw->areas2[i];
						}
						#endif
					}
					else
					{
						cout << endl << INFO << "Cavity " << i << " vol " << pw->volumes[i] << " area " << pw->areas[i];
						if (internals != NULL)
						{
							*internals << endl << "cav_vol " << pw->volumes[i];
							*internals << endl << "cav_area " << pw->areas[i];
						}
					}
				}
				else
				{
					if (isPocket[i])
					{
						cout << endl << INFO << "Pocket " << i << " vol " << pw->volumes[i];
						if (internals != NULL)
							*internals << endl << "pocket_vol " << pw->volumes[i];
					}
					else
					{
						cout << endl << INFO << "Cavity " << i << " vol " << pw->volumes[i];
						if (internals != NULL)
							*internals << endl << "cav_vol " << pw->volumes[i];
					}
				}
			}
			else
			{
				if (conf.tri)
				{
					#if !defined(MULTITHREADED_POCKET_LOOP)
					cout << endl << INFO << "Pocket " << i << " vol " << pw->volumes[i] << " area " << pw->areas[i] << " body area " << pw->areas2[i];
					#else
					cout << endl << INFO << "Pocket " << i << " vol " << pw->volumes[i] << " area " << pw->areas[i] << " body area " << pw->areas2[i];
					#endif
					if (internals != NULL)
					{
						*internals << endl << "pocket_vol " << pw->volumes[i];
						*internals << endl << "pocket_area " << pw->areas[i];
					}
				}
				else
				{
					cout << endl << INFO << "Pocket " << i << " vol " << pw->volumes[i];
					if (internals != NULL)
						*internals << endl << "pocket_vol " << pw->volumes[i];
				}
			}
		}
	}

	chrono::duration<double> pocket_detection_time = chrono_end - chrono_start;
	cout << endl << INFO << "Pocket detection time.. ";
	printf ("%.4e [s]", pocket_detection_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	double current_mem_in_MB, peak_mem_in_MB;
	getMemSpace (current_mem_in_MB, peak_mem_in_MB);
	cout << endl << INFO << "Memory required after pocket detection is " << current_mem_in_MB << " MB";
	#endif

	cout << endl << INFO << "Cleaning memory...";
	cout.flush();

	delete[] pw->areas;
	delete[] pw->areas2;
	delete[] pw->volumes;

	delete pw;

	delete surf1;
	delete dg1;

	delete surf2;
	delete dg2;
	cout << "ok!";
}


/** it computes and returns the membFit score*/
void membfitMode(Surface *surf, DelPhiShared *dg)
{
	auto chrono_start = chrono::high_resolution_clock::now();

	char refName[100];
	strcpy(refName,conf.sysName.c_str());

	double surf_volume;

	// Pre-process surface
	bool outsurf = surf->build();

	if (!outsurf)
	{	
		cout << endl << ERR << "Surface construction failed!" << endl;
		exit(-1);
	}

	if (conf.membHeight <= 0)
	{
		cout << endl << ERR << "Cannot have a negative or null membrane height\n";
		exit(-1);
	}

	if (conf.membHeight <= 0)
	{
		cout << endl << ERR << "Cannot have a negative or null membrane height\n";
		exit(-1);
	}

	// Build DelPhi stuff
	surf->getSurf(&surf_volume, conf.optimizeGrids, conf.fillCavities, conf.cavVol);

	auto chrono_end = chrono::high_resolution_clock::now();

	chrono::duration<double> surface_time = chrono_end - chrono_start;
	cout << endl << INFO << "Surface computation time.. ";
	printf ("%.4e [s]", surface_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	double current_mem_in_MB, peak_mem_in_MB;
	getMemSpace (current_mem_in_MB, peak_mem_in_MB);
	cout << endl << INFO << "Memory required after call of getSurf() in membfitMode() is " << current_mem_in_MB << " MB";
	#endif

	cout << endl << INFO << "Estimated volume " << setprecision(10) << surf_volume << " [A^3]";
	
	// get the minimum Z from which to start the membrane shift
	double lowerZ = 0, upperZ = 0;
	double minBox[3], maxBox[3];
	double maxScore = -INFINITY;

	dg->getBounds(minBox, maxBox);
	
	lowerZ = minBox[2];
	upperZ = lowerZ + conf.membHeight;
	
	cout << endl << INFO << "Scanning along Z direction...";
	int nshift = 0;
	double bestLower = +INFINITY;
	double bestUpper = -INFINITY;

	while (upperZ < maxBox[2])
	{
		// get the starting z value
		int z_start = (int)rintp((lowerZ-dg->zmin)/dg->side);
		int z_end   = (int)rintp((upperZ-dg->zmin)/dg->side);

		cout << endl << INFO << "Testing Planes z=" << lowerZ << ",z=" << upperZ;
		unsigned int iscore = 0;
		unsigned int tot = 0;
		
		for (int iz=z_start; iz<z_end; iz++)
			for (int iy=0; iy<dg->ny; iy++)
				for (int ix=0; ix<dg->nx; ix++)
				{
					int value = read4DVector<int>(dg->epsmap,ix,iy,iz,0,dg->nx,dg->ny,dg->nz,3);
					
					// skip solvent
					if (!value)
						continue;

					int vote = 0;

					// consider it only if it is solvent exposed, thus check the nearby cubes
					for (unsigned int k=0; k<SHIFT_MAP; k++)
					{
						int cx = ix + shift_map[k][0];
						int cy = iy + shift_map[k][1];
						int cz = iz + shift_map[k][2];

						if (cx>(dg->nx-1) || cy>(dg->ny-1) || cz>(dg->nz-1) || cx<0 || cy<0 || cz<0)
							continue;

						// if solvent exposed at least one nearby cube, then the grid cube is judged as
						// solvent exposed
						if (!read4DVector<int>(dg->epsmap,cx,cy,cz,0,dg->nx,dg->ny,dg->nz,3))
							vote++;
					}

					//cout << endl << vote  << " " << SHIFT_MAP/2.;
					if (vote > (SHIFT_MAP>>1))
					{
						tot++;
						if (value == HYDROPHOBIC)
							iscore++;
					}
				}
		
		double score = (double)(100*iscore)/tot;
		
		if (score > maxScore)
		{
			maxScore = score;
			bestLower = lowerZ;
			bestUpper = upperZ;
		}
		
		// get the final z value
		lowerZ += conf.membShift;
		upperZ += conf.membShift;

		nshift++;
	}
	cout << endl << INFO << "Best Upper " << bestUpper;
	cout << endl << INFO << "Best Lower " << bestLower;
	cout << endl << INFO << "Max Score " << maxScore;
}
