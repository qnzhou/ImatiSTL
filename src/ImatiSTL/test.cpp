#include "imatistl.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

using namespace IMATI_STL;

void usage()
{
 printf("\n%s V%s - by %s\n------\n", ImatiSTL::app_name, ImatiSTL::app_version, ImatiSTL::app_authors);
 printf("Usage: %s inmeshfile [outmeshfile] [-e] [-E val] [-f] [-g grid_size] [-j] [-o val] [-q grid_size] [-r] [-s] [-t secs] [-v] [-x]\n", ImatiSTL::app_name);
 printf("  Processes 'inmeshfile' and saves the result to 'outmeshfile'\n");
 printf("  If 'outmeshfile' is not specified 'inmeshfile_fixed.stl' will be produced\n");
 printf("  Option '-e' computes the exact outer hull (slow, use with CAD-like models). This is the default switch.\n");
 printf("  Option '-E val' as '-e' but thickens visible sheets (val-offset) and re-iterates (max 3 times).\n");
 printf("  Option '-S val' as '-e' but thickens all the triangles in advance.\n");
 printf("  Option '-r' rotates the mesh before fixing.\n");
 printf("  Option '-o val' produces an offset of the input at radius 'val'.\n");
 printf("  Option '-f' fills surface holes before starting.\n");
 printf("  Option '-q' quantizes coordinates within a grid_size^3 grid before starting.\n");
 printf("  Option '-s' just saves the output file without modifications.\n");
 printf("  Option '-g' makes the algorithm remesh the input within a uniform grid whose\n");
 printf("  size is 'grid_size'.\n");
 printf("  Option '-j' = output files in OFF format insted of STL\n");
 printf("  Option '-F' = output files in EFF format (exact) insted of STL\n");
 printf("  Option '-t secs' interrupts computation after 'secs' seconds\n");
 printf("  Option '-v' saves additional files for debugging purposes.\n");
 printf("  Option '-x' exits if output file already exists.\n");
 printf("  Accepted input formats are STL, PLY and OFF.\n");
 printf("\nHIT ENTER TO EXIT.\n");
 getchar();
 exit(0);
}

char *createFilename(const char *iname, const char *subext, char *oname, const char *newextension)
{
	static char tname[2048];
	strcpy(tname, iname);
	for (int n = strlen(tname) - 1; n>0; n--) if (tname[n] == '.') { tname[n] = '\0'; break; }
	sprintf(oname, "%s%s%s", tname, subext, newextension);
	return oname;
}


int main(int argc, char *argv[])
{
 ImatiSTL::init(); // This is mandatory
 ImatiSTL::app_name = "STLFix";
 ImatiSTL::app_version = "1.3";
 ImatiSTL::app_year = "2016";
 ImatiSTL::app_authors = "Marco Attene";
 ImatiSTL::app_maillist = "attene@ge.imati.cnr.it";

 clock_t beginning = clock();

 // Uncomment the following to prevent message reporting
 // ImatiSTL::quiet = true;

 TriMesh tin;
 int grid_size = 0;
 int quantization = 0;
 int timeout_secs = 0;
 bool justsave = false;
 bool rotate = false;
 bool holefill = false;
 bool exact = true;
 bool off_output = false;
 bool eff_output = false;
 bool Exact = false;
 bool Shelled_Exact = false;
 bool skip_if_fixed = false;
 bool debug = false;
 double offset = 0.0;
 coord thickness = 0.12;

 char infilename[2048], outfilename[2048], dglfilename[2048], extension[] = ".stl";

 if (argc < 2) usage();

 float par;
 int i = 2;
 if (argc > 2 && argv[2][0] == '-') i--;

 for (; i<argc; i++)
 {
  if (i<argc-1) par = (float)atof(argv[i+1]); else par = 0;
  if (!strcmp(argv[i], "-g"))
  {
	  if (par < 2) ImatiSTL::error("Grid_size must be > 1.\n");
	  if (par > 10000) ImatiSTL::error("Grid_size must be < 10000.\n");
	  grid_size = (int)par;
  }
  else if (!strcmp(argv[i], "-q"))
  {
	  if (par < 2) ImatiSTL::error("Grid_size must be > 1.\n");
	  if (par > 1000000) ImatiSTL::error("Grid_size must be < 1000000.\n");
	  quantization = (int)par;
  }
  else if (!strcmp(argv[i], "-E"))
  {
	  if (par < 0) ImatiSTL::error("thickening value must be >= 0.\n");
	  thickness = par;
	  Exact = true;
  } else if (!strcmp(argv[i], "-S"))
  {
	  if (par <= 0) ImatiSTL::error("thickening value must be > 0.\n");
	  thickness = par;
	  Shelled_Exact = true;
  } else if (!strcmp(argv[i], "-t"))
  {
	  if (par <= 0) ImatiSTL::error("timeout must be positive.\n");
	  timeout_secs = (int)par;
  } 
  else if (!strcmp(argv[i], "-s")) justsave = true;
  else if (!strcmp(argv[i], "-r")) rotate = true;
  else if (!strcmp(argv[i], "-x")) skip_if_fixed = true;
  else if (!strcmp(argv[i], "-o")) offset = par;
  else if (!strcmp(argv[i], "-f")) holefill = true;
  else if (!strcmp(argv[i], "-e")) exact = true;
  else if (!strcmp(argv[i], "-j")) off_output = true;
  else if (!strcmp(argv[i], "-F")) eff_output = true;
  else if (!strcmp(argv[i], "-v")) debug = true;
  else if (argv[i][0] == '-') ImatiSTL::warning("%s - Unknown operation.\n", argv[i]);

  if (par) i++;
 }

 sprintf(infilename, "%s", argv[1]);
 if (off_output) strcpy(extension, ".off");
 if (eff_output) strcpy(extension, ".eff");
 if (argc>2 && argv[2][0] != '-') sprintf(outfilename, "%s", argv[2]);
 else createFilename(infilename, "_fixed", outfilename, extension);
 createFilename(infilename, "_dangling", dglfilename, extension);

 if (skip_if_fixed && fopen(outfilename, "r")) ImatiSTL::error("Output file already exists (-x option specified).");

 // The loader automatically reconstructs a manifold triangle connectivity
 if (tin.load(infilename) != 0) ImatiSTL::error("Can't open file.\n");

 tin.rebuildConnectivity();
 tin.printReport();

 //Node *n;
 //Triangle *t;
 //tin.deselectTriangles();
 //FOREACHVTTRIANGLE((&tin.T), t, n) if (t->getColor()==DEFAULT_TRIANGLE_COLOR)
 //{
	// tin.selectConnectedComponent(t);
	// Node *m;
	// uint32_t col = Triangle::random_color();
	// FOREACHVTTRIANGLE((&tin.T), t, m) if (IS_VISITED(t)) {
	//	 t->setColor(col); UNMARK_VISIT(t);
	// }
 //}

 if (!justsave)
 {
   if (rotate)
   {
		 ImatiSTL::info("Rotating mesh...\n");
		 Matrix4x4 rot;
		 rot.setRotation(sqrt(3.0) / 2.0, sqrt(3.0) / 4.0, 0.25, 0.0);
		 tin.transform(rot);
   }

   if (offset)
   {
	   ImatiSTL::info("Offsetting mesh...\n");
	   tin.toRawOffset(coord(offset));
   }

   if (holefill)
   {
	   if (tin.fillSmallBoundaries(tin.E.numels()))
	   {
		   ImatiSTL::warning("Surface holes have been patched.\n");
		   tin.deselectTriangles();
	   }
   }

   if (quantization)
   {
	   ImatiSTL::info("Quantizing mesh...\n");
	   tin.quantize(quantization);
	   ImatiSTL::useRationals(false);
   }

   if (grid_size)
   {
	if (Shelled_Exact)
	{
		ImatiSTL::info("Thickening ...\n");
		tin.toThinShell(thickness);
	}
    ImatiSTL::info("Remeshing ...\n");
	tin.latticeRemesh(grid_size, true);
   }
   else if (Exact)
   {
	   int n=0, iteration = 1;
	   if (thickness==0) thickness = tin.bboxLongestDiagonal() / 1000.0; // This is a parameter to be set depending on application
	   double epsilon_area = tin.area()*1.0e-6;				 // This also should be considered to be app-dependent
	   double epsilon_squelen = TMESH_TO_DOUBLE(thickness*thickness)*1.0e-6;				 // This also should be considered to be app-dependent

	   do
	   {
		printf("ITERATION %d\n", iteration++);
		tin.safeCoordBackApproximation();	// When iterating, this ensures that the process does not slow down too much
		TriMesh *dangling = tin.computeOuterHull(debug, timeout_secs);
		n = dangling->removeSmallestComponents(epsilon_area); if (n) printf("%d small dangling shells removed\n", n);
		n = tin.removeSmallestComponents(epsilon_area); if (n) printf("%d small solid shells removed\n", n);
		n = dangling->T.numels();

		tin.removeRedundantVertices(); tin.delaunizeFlatAreas();
		dangling->removeRedundantVertices(); dangling->delaunizeFlatAreas();

		if (n != 0 && iteration<4) { dangling->toThinShell(thickness); tin.moveMeshElements(dangling); }
		tin.deselectTriangles();
	   } while (n!=0 && iteration<4);
   }
   else if (exact || Shelled_Exact)
   {
	   if (Shelled_Exact) tin.toThinShell(thickness);
	   TriMesh *dangling = tin.computeOuterHull(debug, timeout_secs);

	   if (dangling->T.numels() && dangling->area() > tin.area()*1.0e-9)
	   {
		   ImatiSTL::warning("Dangling open surfaces remained that do not bound any solid. These are saved to '%s'.\n", dglfilename);
		   if (!eff_output) dangling->safeCoordBackApproximation();
		   dangling->save(dglfilename);
	   }
   }
   else
   {
	   // Keep only the largest component (i.e. with most triangles)
	   int sc = tin.removeSmallestComponents();
	   if (sc) ImatiSTL::warning("Removed %d small components\n",sc);

	   // Fill holes
	   if (tin.boundaries())
	   {
		ImatiSTL::warning("Patching holes\n");
		tin.fillSmallBoundaries(tin.E.numels());
	   }

	   // Run geometry correction
	   if (!tin.boundaries()) ImatiSTL::warning("Fixing degeneracies and intersections...\n");
	   if (tin.boundaries() || !tin.meshclean()) ImatiSTL::warning("STLFix could not fix everything.\n", sc);
   }
 }

 //if (!eff_output)
 //{
	// ImatiSTL::useRationals(true);
	// tin.delaunizeFlatAreas();
	// tin.safeCoordBackApproximation();
	// ImatiSTL::useRationals(false);
	// ImatiSTL::info("Checking output mesh ...\n");
	// if (tin.selectIntersectingTriangles(50U, true) != 0) ImatiSTL::warning("\n *********** ASCII conversion introduced new intersections !!! *********************** \n\n");
 //}

 ImatiSTL::info("Saving output mesh ...\n");
 tin.save(outfilename);

 printf("Elapsed time: %d ms\n", clock() - beginning);

 tin.printReport();

 return 0;
}
