#ifndef _kodtree_
#define _kodtree_


#include<list>
#include<vector>
#include <limits>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include "mex.h"
using namespace std;

typedef double Point[3];//自定义point为3个double数组
typedef double Box[6];
struct WpInfo{
	void *info;//定义通用型指针
	int infotype;
	bool get;
	int rcount;
	WpInfo( void *inf,int infot):info(inf),infotype(infot),get(false),rcount(0){}
};

struct WpVert{
	void *vt;
	bool vget;
	int rcount;
	WpVert(void *vin):vt(vin),vget(false),rcount(0){}
	WpVert():vt(0),vget(false),rcount(0){}
};
typedef WpVert * PtWpVert;
 struct CellNode;
template<typename   T>   
T* renew(T* &p,size_t oldsize,size_t newsize){   
   T* tmp=new T[newsize];   
   memcpy(tmp,p,oldsize*sizeof(T));   
   delete [] p;   
   return p=tmp;   
}   

class Kodtree{
public:
	typedef void (*Funcpointofvert)(Point p,void *v);
	typedef bool (*Funcexinfoshouldbeincell)(void *info, int infotype, CellNode *cnode);
	typedef bool (*Funcexinfooverlapbox) (void *info, int infotype, const Box &bd,double epsi);
	Kodtree(const Box &bd,Funcpointofvert pofv, const int capacity=10,double epsi=0);
	Kodtree(void **vti, int numvi,Funcpointofvert pofv,const int capacity=10,double epsi=0);
	Kodtree(void **vti, int numvi,const Box &bd,Funcpointofvert pofv, const int capacity=10,double epsi=0);
	~Kodtree();
	CellNode *getRoot(void){ return root;}
	double getEpsCell(void){return epscell;}
	WpVert *insertVert(void *v){
		WpVert *nvert= new WpVert(v);
		Point p;
		pofv(p,v);
		insertWpVertInSubTree(p,nvert,root);
		if(nvert->rcount==0){
			delete nvert;
			return 0;
		}
		return nvert;
	}
	bool deleteVert(void *v){
		Point p;
		pofv(p,v);
		if(!isVertRecordInSubTree(p,v,root)) return false;
		deleteVertInSubTree(p,v,root);
		checkAndMergeSubTreeAfterDelete(p,root);
		return true;
	}
	WpInfo * insertExinfo(void *info,int infotype){
		WpInfo *nwinf=new WpInfo(info,infotype);
		insertWpInfoInSubTree(nwinf,root);
		if(nwinf->rcount==0){
			delete nwinf;
			return 0;
		}
		return nwinf;
	}
	void deleteExinfo(void *info,int infotype){
		deleteExinfoInSubTree(info,infotype,root);
	}
	void collectVertsWithBox(const Box &bd, std::list<void *> &lvert);
	void collectVertsWithCell(CellNode *cnode, std::vector<void *> &vecvert);
	void collectExinfoWithBox(const Box &bd, int infotype,std::list<void *> &linfo);
	void collectExinfoWithCell(CellNode *cnode, int infotype,std::list<void *> &lexinfo);
	void setFuncExinfoShouldbeInCell(Funcexinfoshouldbeincell infunc) {ifExinfoShouldbeInCell=infunc;}
	void setFuncExinfoOverlapBox(Funcexinfooverlapbox infunc){ifExinfoOverlapBox=infunc;}
	CellNode *findaLeafCellContainingPoint(CellNode *pcell,Point p);
	CellNode *findTheNearestAncestorContainingPoint(CellNode *pcell,Point pcha);
	void freeSubTree(CellNode *pcell);
	Kodtree();
private:
	bool isVertRecordInSubTree( const Point &p, void *v,CellNode *cnode);
	void insertWpVertInSubTree( const Point &p, WpVert *nv,CellNode *cnode);
	void insertWpInfoInSubTree(WpInfo *pwinfo, CellNode *cnode);
	void collectWpVertsWithBoxInSubTree(CellNode *cnode,const Box &bd,std::list<WpVert *> &lvert);
	void collectWpinfoWithBoxInSubTree(CellNode *cnode,const Box &bd,int infotype,std::list<WpInfo *> &lwpinfo);
	void deleteVertInSubTree(const Point &p,void *v,CellNode *cnode);
	void deleteExinfoInSubTree(void *info,int infotype, CellNode *cnode);
	void checkAndRemoveSurplusWpInfoAfterMerge(CellNode *cnode);
	void checkAndMergeSubTreeAfterDelete(const Point &p,CellNode *cnode);
	void mergeSubTree(CellNode *cnode);
	void merge2SubCellWpVert(CellNode *cnode);
	void merge2SubCellWpInfo(CellNode *cnode);
	bool if2CellNeighb(CellNode *pcell0, CellNode *pcell1);
	bool canNodeSplit(CellNode *cnode);
	void splitNode(CellNode *cnode);

private:
//	static const double epsilonon;
//	static const double epscoplanar;
	double epsoverlap;
	int cellcapacity;
	Funcpointofvert pofv;
	Funcexinfoshouldbeincell ifExinfoShouldbeInCell;
	Funcexinfooverlapbox ifExinfoOverlapBox;
	double epscell;
	CellNode *root;
};

struct CellNode{
	WpVert **vert;
	int numvert;
	int nodecapacity;
	std::list<WpInfo *> *lpwpinfo;
	Box bound;
	CellNode *child[2];
	CellNode *parent;
	int inoutattrib;
	CellNode(const Box &bd);
	~CellNode();
	CellNode *anotherChild(CellNode *pcell){
		if(pcell==0) return this;
		if(child[0]==pcell) return child[1];
		else return child[0];
	}
	bool isEmpty(){ return vert==0;}
	bool isLeaf(){ return !child[0];}
};
 
extern void jf_error(char *ch);
extern bool ifBoxContainPoint( Point p,const Box &bound,const Box &rootbound);
extern bool if2BoxOverlap(const Box &a,const Box &b);
extern bool if2BoxNeighb(const Box &a,const Box &b);
extern bool ifPointOverlapWithBox(const Point &p,const Box &bd,const Box &rootbound,double eps);
extern void boxOfVerts( void **v , int num ,Box b ,void (*Funcpointofv)(Point p,void *v));
extern double sqdistPointToBox(Point p,const Box &bd);
extern void getTheLongestDistOfBox(const Box &b,int &di, double *pdist=0);
extern void copy3DPoint(const Point &pfr,Point &pto); //&?
extern int IsTriangleBoxInt(Point p1 ,Point p2 ,Point p3 , double box[6] );
extern double sqdistInnerPointToBoxBound(Point p,const Box &bd);

extern void vec_2p(double * , double * , double *) ;
extern int vec_uni(double *) ;
void vec_neg(double *vector);
extern double vec_val( double * ) ;
extern double vec_sqval( double *vec );
extern double vec_dotp(double * ,double *) ;
extern void vec_crop(double * ,double * ,double * ) ;
extern double Distance3D( double * ,double *) ;
extern double SqDistance3D( double * ,double * ) ;
extern double vec_blep( double * , double * , double * ) ;
extern void norm_3p( double * , double * ,double * ,double * ) ;
extern double VolumOf4p(Point ,Point ,Point ,Point ) ;
extern bool isTriangleBoxOver(Point p1 ,Point p2 ,Point p3 ,const double bd[6],double eps );
extern bool isTriRayIntersect(double p0[3],double p1[3],double p2[3],double orig[3],double dir[3]);//, double &t){

#endif //_kodtree_