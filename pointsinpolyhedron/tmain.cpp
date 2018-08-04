#include "pinpolyhedron.h"
#include <string>
//#include "mex.h"
using namespace std;

PointInPolyhedron* ptpoly = 0;
extern int absolute;
extern void randompoint(double p[3],double bd[6]);
extern void boxOfPoints( double (*p)[3] , int num ,double box[6] );

void boxOfPoints(double(*p)[3], int num, double box[6]) {
  int i, j;

  if (num <= 0)
    jf_error("boxofP");

  for (i = 0; i < 3; ++i)
    box[i] = box[i+3] = p[0][i];

  for (j = 1; j < num; ++j) {
    for (i = 0; i < 3; ++i) {
     if (p[j][i] < box[i])
       box[i] = p[j][i];
     if (p[j][i] > box[i+3])
       box[i + 3] = p[j][i];
    }
  }

  double a = max(box[5] - box[2], max(box[3] - box[0], box[4] - box[1]));
  for (i = 0; i < 3; ++i) {
    box[i]   -= 0.01 * a;
    box[i+2] += 0.01 * a;
  }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  int i, j; 
  string errmsg;

  int numvert, numtri, nummaterial;
  double (*vert)[3];
  int (*tris)[3], (*twoma)[2];

  vert  = NULL;
  tris  = NULL;
  twoma = NULL;
  
  absolute=0;
  ////////////////////
  // Check Params
  //
  if (nrhs != 3)
    mexErrMsgTxt("Three input arguement required.");
  if (nlhs > 1)
    mexErrMsgTxt("One output arguement required.");


  ////////////////////
  // Ectract Params
  //

  // Extract param 0£º list of  perps
  int perpm = mxGetM(prhs[0]);
  int perpn = mxGetN(prhs[0]);
  if (!mxIsDouble(prhs[0]))
     mexErrMsgTxt("The classID of first input must be double.");
  if (perpn != 3)
    mexErrMsgTxt("The columns of first input must be 3.");
  double* perp = mxGetPr(prhs[0]);
  
  // Extract param 1: list of vertexs
  int vertm = mxGetM(prhs[1]);
  int vertn = mxGetN(prhs[1]);
  if (!mxIsDouble(prhs[1]))
    mexErrMsgTxt("The classID of second input must be double.");
  if (vertn != 3)
    mexErrMsgTxt("The columns of second input must be 3.");
  numvert = vertm;

  double* vert1 = mxGetPr(prhs[1]);
  vert = (double(*)[3]) new double[3 * numvert];
  for (j = 0; j < 3; ++j)
    for (i = 0; i < numvert; ++i)
      vert[i][j] = vert1[j * numvert + i];

  // Extract param 1: list of triangles
  int trism = mxGetM(prhs[2]);
  int trisn = mxGetN(prhs[2]);
  if (!mxIsDouble(prhs[2]) && !mxIsInt32(prhs[2])) {
    delete[] vert;  
    mexErrMsgTxt("The classID of third input must be double or int32.");
  }
  if (trisn != 5 && trisn != 3) {
    delete[] vert;
    mexErrMsgTxt("The columns of third input must be 3 or 5.");
  } 
  numtri = trism;
  
  tris  = (int(*)[3]) new int[3 * numtri];
  twoma = (int(*)[2]) new int[2 * numtri];
  for (i = 0; i < numtri; ++i)
    twoma[i][0] = -1, twoma[i][1] = 0;

  if (!mxIsDouble(prhs[2]) ) {
    int* pertris_int = (int*) mxGetData(prhs[2]);
	int min_int=numeric_limits<int>::max();
	for ( i = 0; i < trism * 3; ++i )
		if(min_int > pertris_int[i] )
			min_int = pertris_int[i];
    for (j = 0; j < 3; ++j)
	   for (i = 0; i < numtri; ++i)
		   tris[i][j] = pertris_int[j * numtri + i] - min_int;
    if (trisn == 5)
      for (j = 3; j < 5; ++j)
        for (i = 0; i < numtri; ++i)
          twoma[i][j - 3] = pertris_int[j * numtri + i];
  } else {
    double* pertris_double = (double*) mxGetPr(prhs[2]);
	double min_double=numeric_limits<double>::max();
	for (i = 0; i < trism * 3; ++i)
		if (min_double > pertris_double[i])
			min_double = pertris_double[i];
    for (j = 0; j < 3; ++j)
      for (i = 0; i < numtri; ++i)
		tris[i][j] = pertris_double[j * numtri + i] - min_double;
    if (trisn == 5)
      for (j = 3; j < 5; ++j)
        for (i = 0; i < numtri; ++i)
          twoma[i][j - 3] = pertris_double[j * numtri + i];
  }
  
  // Count number of materials
  int min = 0, max = 0;
  for (j = 0; j < 2; ++j) {
    for (i = 0; i < numtri; ++i){
      if (twoma[i][j] < min) 
        min = twoma[i][j];
      if (twoma[i][j] > max) 
        max = twoma[i][j];
         }
    }
  nummaterial = max - min;
  mexPrintf("Your model is composed of %d different materials.\n",nummaterial);

 
    ptpoly = new PointInPolyhedron(vert, numvert, tris, twoma, numtri, nummaterial);

  
  double bd[6];
  boxOfPoints(vert, numvert, bd);  


  ////////////////////
  // Output results
  //
  plhs[0] = mxCreateDoubleMatrix(perpm, 4, mxREAL);
  double* out = mxGetPr(plhs[0]);
  for (i = 0; i < perpm * perpn; ++i)
     out[i] = perp[i];

  double p[3];
  for (i = 0; i < perpm; ++i) {
    p[0] = perp[i];
    p[1] = perp[i + perpm];
    p[2] = perp[i + perpm * 2];
    out[i + perpm * 3] = (double) ptpoly->isPinPolyhedron(p);
  }
  
  delete[] vert;
  delete[] tris;
  delete[] twoma;
  delete ptpoly;
  return;
}
