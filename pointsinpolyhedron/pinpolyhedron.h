#include <vector>
#include "kodtree.h"
//#include "mex.h"
using namespace std;

class PointInPolyhedron{

public:
	int isPinPolyhedron(double p[3]); //-2,-1,0,1,...nummaterial-1.-2=>on boundary.
	bool isPointInGreyOfPolyhedron(double p[3]);
	PointInPolyhedron(double vert[][3], int numvi,int (*tris)[3],int (*twoma)[2],int numtris,int nummat=1); //,double epsi);
	double absoluteClosestSqDistance(Point p){
		CellNode *cnode=polytree->findaLeafCellContainingPoint(polytree->getRoot(),p);
		int ntri; double dist;
		getAbsoluteClosestTriForPointInGCell(p,cnode,ntri,dist);
		return dist;
	}
	void resetAttrib(void);
	~PointInPolyhedron();
	static double (*vertcoord)[3];
	static int numvert;
	static int (*trips)[3],(*twomaterialsoftri)[2],numtri;
private:
	int testPinPolyhedronForPinGcell(double p[3],CellNode *cnode);
	void getVisibleEntityForPointInGCell( double p[3],CellNode *cnode,int &id,int &nentity, int &ntri,double &dist);
	void getRelativeClosestTriForPointInGCell(double p[3],CellNode *cnode, int &tri, double &dist);
	void getAbsoluteClosestTriForPointInGCell(double p[3],CellNode *cnode, int &tri, double &dist);
	void getTheClosestTriNonLeaf(double p[3],double dist0,CellNode *pcell,int &tri,double &dist);
	void getTheClosestTriAmongCell(double p[3],CellNode *pcell, double &dist,int &ntri);
	void getEndPointOfTri(int tri, double p0[3],double p1[3],double p2[3]);
	void setGCellAttribOfSubTree(CellNode *pcell);
	void getCellSeqWithUnknownAttribFromaCell(CellNode *cnode,vector<CellNode *> *&pcellseq,
					CellNode * &pcellm,int &ia,double pm[3]);
	int getDetermingTriFromVertOfaVisibleTri(double p[3], int tri,int vert);
	int getDetermingTriFromEdgeOfaVisibleTri(double p[3],int va,int vb);
	void getVertsAroundaVertMulti(int v, int * &nbverts,int &numnbv);
	int nextTriOfVert_g(int v,int &k);
	int  nextVertOfTri(int tri,int v);
	void  getThePointFormingLeastAngleWith2Points(double p[3],int v, int *nbverts,
		int numnbv,double &maxcosa,double &mincosa,int &vridge);
	void  get2TriCom2Vert(int va, int vb, int &ta, int &tb);
	int indexOfNeighbTriToTri(int tria,int trib);
	int indexOfVertAtTri(int v, int tri[3]);
	void wrapPointsUpasVerts(void  ** &vti);
	void  getEdgeOfTri(int np[3], int index, int &a, int &b);
	void PointInPolyhedron::sortTrianglesForMultiMaterial(int nummat);
	void sortTrianglesOuterNormAndRecNeighb(double (*vertb)[3],int numvertb,int (*tripsb)[3],int numtrib);
	int getAndSortaLowestTri(double (*vertb)[3],int numvertb,int (*tripsb)[3],bool firstshell);
	void sort1ShellFromaTri(int (*tripsb)[3],int tri);
	int getNeighbTriWithoutTopology(int (*tripsb)[3],int tri,int ind);
	bool triSortAs2Nodes(int tri3p[3],int va, int vb);
	void get2TriCom2NodesWithoutTopology(int (*tripsb)[3],int va, int vb, int &ta, int &tb);
	static void pofvforcoordnodes3(double p[3],void *pv);
	static bool ifexinfooverlapbox(void *info,int infotype,const Box &bd,double eps);
	static bool ifexinfoshouldbeincell(void *info,int infotype,CellNode *cnode);
	void getTheTrisAmongCellsRayCrossing(double p[3],CellNode *pcell0,std::vector<int> &tris);
	int PointInPolyhedron::countRayCrossingFromGcell(double p[3],CellNode *pcell);
	CellNode* getNextCell(CellNode *cnode,double ps[3],double pe[3]);
//	double signed_distance(double p[3],int tri);
private:
	static const double epsilonon;
	static const double epsoverlap;
	static const double epscoplanar;
	//double epscell;
	Kodtree *polytree;
//	double (*tri_plane_edge_norm)[3][3],(*edge_len)[3],(*face_norm)[3],(*edge)[3][3];
	int (*tneighb)[3];
	int *triofnode;
	int *vertattrib;
	bool *vertused;
	int *numtriofnode;
	int *tripositionofnode;
	int *trilist;
	int *trisort;
	int *triused;
	vector<int> testedtris;
};