
#include <vector>
#include <limits>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <time.h>
#include "pinpolyhedron.h"
//#include "mex.h"
extern  void jf_error(char *);
using namespace std;

const double PointInPolyhedron::epsilonon=0.00000000000001;
const double PointInPolyhedron::epsoverlap=0.000001;
const double PointInPolyhedron::epscoplanar=0.000001;

double (*PointInPolyhedron::vertcoord)[3];
int PointInPolyhedron:: numvert;
int (*PointInPolyhedron::trips)[3];
int (*PointInPolyhedron::twomaterialsoftri)[2]=0;
int PointInPolyhedron::numtri;
int absolute;
int *startaddress=(int *)1;

extern int positionOfPointProjectToTri(double p[3],double p0[3],double p1[3],double p2[3]);
extern void getClosestPointOnTriangle(double p[3],double p0[3],double p1[3],double p2[3],double pnear[3]);
extern double sqDistPointToTri(double p[3],double p0[3],double p1[3],double p2[3]);//,double sqdist0);
extern double sqDistPointToSeg3D(double p[3],double p0[3],double p1[3]);
int triIndexFromPt(void *ptri){
	int *ptr=(int *)ptri;
	return ptr-startaddress;
}
int vertIndexFromPt(void *pv){
	int *pt=(int *)pv;
	return pt-startaddress;
}

void  PointInPolyhedron::pofvforcoordnodes3(double p[3],void *pv){

//	extern static int *startaddress;

	int nd=vertIndexFromPt(pv);
	p[0]=vertcoord[nd][0];
	p[1]=vertcoord[nd][1];
	p[2]=vertcoord[nd][2];
}
void PointInPolyhedron::wrapPointsUpasVerts(void  ** &vti){

	vti=new void *[numvert];
	for(int i=0; i<numvert; i++)
		vti[i]=startaddress+i;
}

bool   PointInPolyhedron::ifexinfooverlapbox(void *info,int infotype,const Box &bd,double eps){
	  if(infotype==1){
		  int tri =triIndexFromPt(info);
		  return isTriangleBoxOver(vertcoord[trips[tri][0]],vertcoord[trips[tri][1]],vertcoord[trips[tri][2]],bd,eps);
	  }
	  return false;
  }

bool   PointInPolyhedron::ifexinfoshouldbeincell(void *info,int infotype,CellNode *cnode){
	  if(infotype==1){
		  int tri=triIndexFromPt(info);
		  for(int i=0; i<cnode->numvert; i++){
			  int v=vertIndexFromPt(cnode->vert[i]->vt);
			  if(v==trips[tri][0]||v==trips[tri][1]||v==trips[tri][2])
				  return false;
		  }
	  }
	  return true;
}
bool PointInPolyhedron::isPointInGreyOfPolyhedron(double p[3]){
	
	int rt;
	CellNode *pcell=polytree->findaLeafCellContainingPoint(polytree->getRoot(),p);
	if(pcell->inoutattrib ==-2) 
		return true;
	else 
		return false;
}
int PointInPolyhedron::isPinPolyhedron( double p[3]){

	int rt;
	CellNode *pcell;
	vector<CellNode *> *pcellseq;

	pcell=polytree->findaLeafCellContainingPoint(polytree->getRoot(),p);
	if(pcell==0) 
		return -1;
	if(pcell->inoutattrib!=-3&&pcell->inoutattrib!=-2) //==-1||pcell->inoutattrib==1)-3=>initial(non-characterized empty cell),-2=>gray cell.-1=>out
		return pcell->inoutattrib;
	else if(pcell->inoutattrib ==-2)
		return testPinPolyhedronForPinGcell(p,pcell);
	//else if(pcell->inoutattrib==-3)
	//if((pcell->inoutattrib=testPinPolyhedronForPinGcell(p,pcell))==-2)
	//	jf_error("err ispointin");
	//else return pcell->inoutattrib;
	double pm[3];
	CellNode *pcellm=0;
	getCellSeqWithUnknownAttribFromaCell(pcell,pcellseq,pcellm,rt,pm);
	if(rt==-2)
		rt=testPinPolyhedronForPinGcell(pm,pcellm);
	if(rt==-2)
		jf_error("ispinoPolyhedron");
	if(pcellseq!=0) // what's the meaning, will it be colored when only pcellm alone?
		for(unsigned i=0; i<pcellseq->size(); i++)
			(*pcellseq)[i]->inoutattrib=rt;
	delete pcellseq;
	return rt;

}

int PointInPolyhedron::testPinPolyhedronForPinGcell(double p[3],CellNode *cnode){

	if(absolute==2){ 
		int count=countRayCrossingFromGcell(p,cnode);	
		return count%2==0?-1:0;
	}
	int id,nentity,tri,rt;
	double dist,p0[3],p1[3],p2[3];

	getVisibleEntityForPointInGCell(p,cnode,id,nentity,tri,dist);
	if(dist<=epsilonon)
		return -2;
	//for multiRCT
	if(id==0){
		if(vertattrib[nentity]!=-3&&vertattrib[nentity]!=-2)
			return vertattrib[nentity];
		else if(vertattrib[nentity]==-3)
			tri=getDetermingTriFromVertOfaVisibleTri(p,tri,nentity); //vertattrib[] maybe changed inside
	}else if(id==1)
		tri=getDetermingTriFromEdgeOfaVisibleTri(p,trips[tri][(nentity+1)%3],trips[tri][(nentity+2)%3]);
	if(id!=1&&id!=0&&id!=2)
		jf_error("err ispoinPolyhedron");
	getEndPointOfTri(tri,p0,p1,p2); //id==2||id==1&&coplanar at edge||id==0&&coplanar at vertex
	if(VolumOf4p(p0,p1,p2,p)<0)
		rt= twomaterialsoftri[tri][1];
	else 
		rt= twomaterialsoftri[tri][0];
	if(id==0&&vertattrib[nentity]==-3) vertattrib[nentity]=rt;
	return rt;
}


void PointInPolyhedron::getVisibleEntityForPointInGCell( double p[3],CellNode *cnode,int &id,
														   int &nentity, int &ntri,double &dist){

	int ip;
	double p0[3],p1[3],p2[3];

	if(absolute==0)
		getRelativeClosestTriForPointInGCell(p,cnode,ntri,dist);
	else
		getAbsoluteClosestTriForPointInGCell(p,cnode,ntri,dist);
	//return  ; // for multiRCT
	if(dist==numeric_limits<double>::max())
		jf_error("err getrelativeclosetentityforpingcell");
	getEndPointOfTri(ntri,p0,p1,p2);
	if((ip=positionOfPointProjectToTri(p,p0,p1,p2))==6){
		nentity=ntri;
		id=2;
	}else if(ip<3){
		nentity=trips[ntri][ip];
		id=0;
	}else{
//		nentity=neighbOfTri(ntri,ip-3);
//		if(nentity<0) jf_error("getrealvie");
		nentity=ip-3;
		id=1;
	}
}
/*
void PointInPolyhedron::getRelativeClosestTriForPointInGCell(double p[3],CellNode *cnode, int &tri, double &dist){

	CellNode *pcell0=0;
	CellNode *pcell=cnode;
//	if(!pcell||pcell->isEmpty())
//		jf_error("getrelativeclosesttri");
	dist=numeric_limits<double>::max();
	tri=-1;
	double pcha[3]={0,0,dist},p0[3],p1[3],p2[3];
	
	while(pcell){
		int trin;
		double distn;
		bool triupdate;
		triupdate=false;
		getTheClosestTriNonLeaf(p,dist,pcell->anotherChild(pcell0),trin,distn);
		if(distn<dist){
			dist=distn; tri=trin; triupdate=true;
		}
		if(sqdistInnerPointToBoxBound(p,pcell->bound)>=dist) break;
		if(triupdate){
			getEndPointOfTri(tri,p0,p1,p2);
			getClosestPointOnTriangle(p,p0,p1,p2,pcha);
		}
		if(ifBoxContainPoint(pcha,pcell->bound,pcell->bound)) break;
		pcell0=pcell;
		pcell=pcell->parent;
	}
	for(unsigned i=0; i<testedtris.size(); i++)
		triused[testedtris[i]]=0;
	testedtris.clear();
}*/
void findOutPointofBox(double ps[3],double pe[3],double *pl,double *ph,double eps,double px[3]){

	for(int i=0; i<3; i++){
		if(ps[i]<=pe[i]){
			if(pe[i]<=ph[i]) px[i]=pe[i];
			else px[i]=ph[i]+eps;
		}else{
			if(pe[i]<pl[i]) px[i]=pl[i]-eps;
			else px[i]=pe[i];
		}
	}
}
CellNode* PointInPolyhedron::getNextCell(CellNode *cnode,double ps[3],double pe[3]){

	double px[3];
	CellNode *pancestor,*pcell;

	findOutPointofBox(ps,pe,cnode->bound,cnode->bound+3,polytree->getEpsCell(),px);
	if((pancestor=polytree->findTheNearestAncestorContainingPoint(cnode,px))==0)
		return 0;
	else pcell=polytree->findaLeafCellContainingPoint(pancestor,px);
	if(pcell==cnode)
		jf_error("err epscell, contact the developer please,liujianfei@pku.edu.cn");
	return pcell;
}
void PointInPolyhedron::getRelativeClosestTriForPointInGCell(double p[3],CellNode *cnode, int &tri, double &dist){

//	CellNode *pcell0=0;
	CellNode *pcell=cnode;
//	if(!pcell||pcell->isEmpty())
//		jf_error("getrelativeclosesttri");
	dist=numeric_limits<double>::max();
	tri=-1;
	double pcha[3]={0,0,dist},p0[3],p1[3],p2[3];
	
	for(;;){
		int trin;
		double distn;
		bool triupdate;
		triupdate=false;
		getTheClosestTriAmongCell(p,pcell,distn,trin);
		if(distn<dist){
			dist=distn; tri=trin; triupdate=true;
			getEndPointOfTri(tri,p0,p1,p2);
			getClosestPointOnTriangle(p,p0,p1,p2,pcha);
		}
		if(ifBoxContainPoint(pcha,pcell->bound,pcell->bound)) break;
		if(triupdate)
			pcell=getNextCell(cnode,p,pcha);
		else
			pcell=getNextCell(pcell,p,pcha);
		if(pcell==0) jf_error(" err getrelative");
	}
	for(unsigned i=0; i<testedtris.size(); i++)
		triused[testedtris[i]]=0;
	testedtris.clear();
}

void PointInPolyhedron::getAbsoluteClosestTriForPointInGCell(double p[3],CellNode *cnode, int &tri, double &dist){

	CellNode *pcell0=0;
	CellNode *pcell=cnode;
//	if(!pcell||pcell->isEmpty())
//		jf_error("getrelativeclosesttri");
	dist=numeric_limits<double>::max();
	tri=-1;
	
	while(pcell){
		int trin;
		double distn;
		getTheClosestTriNonLeaf(p,dist,pcell->anotherChild(pcell0),trin,distn);
		if(distn<dist){
			dist=distn; tri=trin;
		}
		if(sqdistInnerPointToBoxBound(p,pcell->bound)>=dist)
			break;
		pcell0=pcell;
		pcell=pcell->parent;
	}
	for(unsigned i=0; i<testedtris.size(); i++)
		triused[testedtris[i]]=0;
	testedtris.clear();
}


void PointInPolyhedron::getTheClosestTriNonLeaf(double p[3],double dist0,
										   CellNode *pcell,int &tri,double &dist){
	
	double distn;
	int trin;

	dist=dist0, tri=-1;
	if(sqdistPointToBox(p,pcell->bound)>=dist0) return;
	if(pcell->isLeaf()){
		getTheClosestTriAmongCell(p,pcell,distn,trin);
		if(distn<dist){	dist=distn;	tri=trin; return;}
	}else{
		CellNode *sortsub[2]={pcell->child[0],pcell->child[1]};
		if(sqdistPointToBox(p,pcell->child[0]->bound)>sqdistPointToBox(p,pcell->child[1]->bound)){
			sortsub[0]=pcell->child[1];
			sortsub[1]=pcell->child[0];
		}
		for(int i=0; i<2; i++){
			getTheClosestTriNonLeaf(p,dist,sortsub[i],trin,distn);
			if(distn<dist){ dist=distn; tri=trin; }
		}
	}
}

int PointInPolyhedron::nextTriOfVert_g(int v,int &k){

	k++;
	int currentposi;
	if(k>numtriofnode[v]){cout<<k<<"=k "<<numtriofnode[v]<<"=num "<<v<<endl;jf_error("err nextriofvertg");}
	if(k==numtriofnode[v]) 
		currentposi=tripositionofnode[v];
	else
		currentposi=tripositionofnode[v]+k;
	return trilist[currentposi];

}
void PointInPolyhedron::getTheClosestTriAmongCell(double p[3],CellNode *pcell, double &dist,int &ntri){

	int tri;
	double distemp,p0[3],p1[3],p2[3];

	dist=numeric_limits<double>::max();
	if(!pcell||!pcell->isLeaf())
		jf_error("error gettheclosettriamongcell");
	if(pcell->lpwpinfo!=0)
		for(std::list<WpInfo *>::iterator ite=pcell->lpwpinfo->begin();ite!=pcell->lpwpinfo->end(); ite++){
			if((*ite)->infotype!=1)	continue;
			tri=triIndexFromPt((*ite)->info);
			if(triused[tri]) continue;
			else{ triused[tri]=1;testedtris.push_back(tri);}
			getEndPointOfTri(tri,p0,p1,p2);
	//		sqDistPointToTri(p,p0,p1,p2);
			if((distemp=sqDistPointToTri(p,p0,p1,p2))<dist){
//			if((distemp=signed_distance(p,tri))<dist){
				dist=distemp;
				ntri=tri;
			}
		}
	for(int i=0; i<pcell->numvert; i++){
		int v=vertIndexFromPt(pcell->vert[i]->vt);
		int tri0,tri;
		tri0=tri=triofnode[v];
		int k=0; //multiRCT
		do{
			if(triused[tri]) continue;
			else{ triused[tri]=1;testedtris.push_back(tri);}
			getEndPointOfTri(tri,p0,p1,p2);
		//	sqDistPointToTri(p,p0,p1,p2);
			if((distemp=sqDistPointToTri(p,p0,p1,p2))<dist){
//			if((distemp=signed_distance(p,tri))<dist){
				dist=distemp;
				ntri=tri;
			}
		}while((tri=nextTriOfVert_g(v,k))!=tri0); //multiRCT
	}
}
int PointInPolyhedron::countRayCrossingFromGcell(double p[3],CellNode *pcell){

	std::vector<int> tris;
	getTheTrisAmongCellsRayCrossing(p,pcell,tris);
	int count=0;
	double dir[3]={1,0,0},p0[3],p1[3],p2[3];
	for(unsigned i=0; i<tris.size(); i++){
		getEndPointOfTri(tris[i],p0,p1,p2);
		if(isTriRayIntersect(p0,p1,p2,p,dir)) count++;
	}
	return count;
}

void PointInPolyhedron::getTheTrisAmongCellsRayCrossing(double p[3],CellNode *pcell0,std::vector<int> &tris){
	int tri;
	double pm[3]={p[0],p[1],p[2]};
	CellNode *pcell=pcell0,*pancestor;
	if(!pcell||!pcell->isLeaf()) jf_error("error gettheclosettriamongcell");
	do{
		if(pcell->lpwpinfo!=0)
			for(std::list<WpInfo *>::iterator ite=pcell->lpwpinfo->begin();ite!=pcell->lpwpinfo->end(); ite++){
				if((*ite)->infotype!=1)	continue;
				tri=triIndexFromPt((*ite)->info);
				if(triused[tri]==0){tris.push_back(tri);triused[tri]=1;}
			}
		for(int i=0; i<pcell->numvert; i++){
			int v=vertIndexFromPt(pcell->vert[i]->vt);
			int tri0,tri;
			tri0=tri=triofnode[v];
			int k=0; //multiRCT
			do{
				if(triused[tri]==0){ tris.push_back(tri);triused[tri]=1;}
			}while((tri=nextTriOfVert_g(v,k))!=tri0); //multiRCT
		}
		pm[0]=pcell->bound[3]+polytree->getEpsCell();
//		pm[2]=pcell->bound[5]+polytree->getEpsCell();
		pancestor=polytree->findTheNearestAncestorContainingPoint(pcell,pm);
		if(pancestor==0) break;
	}while((pcell=polytree->findaLeafCellContainingPoint(pancestor,pm))!=0);
	for(unsigned i=0; i<tris.size(); i++)
		triused[tris[i]]=0;
}
void PointInPolyhedron::getEndPointOfTri(int tri, double p0[3],double p1[3],double p2[3]){

	p0[0]=vertcoord[trips[tri][0]][0];
	p0[1]=vertcoord[trips[tri][0]][1];
	p0[2]=vertcoord[trips[tri][0]][2];
	p1[0]=vertcoord[trips[tri][1]][0];
	p1[1]=vertcoord[trips[tri][1]][1];
	p1[2]=vertcoord[trips[tri][1]][2];
	p2[0]=vertcoord[trips[tri][2]][0];
	p2[1]=vertcoord[trips[tri][2]][1];
	p2[2]=vertcoord[trips[tri][2]][2];
}
int positionOfPointProjectToSeg3D(double p[3],double p0[3],double p1[3]){

	double vp0p[3],vp0p1[3],vp1p[3];

	vec_2p(p0,p,vp0p);

	vec_2p(p0,p1,vp0p1);
	if(vec_dotp(vp0p,vp0p1)<=0)
		return -1;
	vec_2p(p1,p,vp1p);
	if(vec_dotp(vp1p,vp0p1)>=0)
		return 1;
	return 0;  
}
int positionOfPointProjectToTri(double p[3],double p0[3],double p1[3],double p2[3]){

	double v0p[3],v20[3],v01[3];
	vec_2p(p0,p,v0p);
	vec_2p(p2,p0,v20);
	vec_2p(p0,p1,v01);

	double d0p20=vec_dotp(v0p,v20);
	double d0p01=vec_dotp(v0p,v01);
	if(d0p20>=0&&d0p01<=0) return 0;
	
	double v1p[3],v12[3];
	vec_2p(p1,p,v1p);
	vec_2p(p1,p2,v12);
	double d1p01=vec_dotp(v1p,v01);
	double d1p12=vec_dotp(v1p,v12);
	if(d1p01>=0&&d1p12<=0) return 1;

	double v2p[3];
	vec_2p(p2,p,v2p);
	double d2p12=vec_dotp(v2p,v12);
	double d2p20=vec_dotp(v2p,v20);
	if(d2p12>=0&&d2p20<=0) return 2;

	double nm012[3],nm01p[3],nm12p[3],nm20p[3];
	vec_crop(v20,v01,nm012);

	vec_crop(v01,v0p,nm01p);
	double dt01=vec_dotp(nm012,nm01p);
	if(dt01<=0&&d0p01>=0&&d1p01<=0)
		return 5; // rt==0;
		
	vec_crop(v12,v1p,nm12p);
	double dt12=vec_dotp(nm012,nm12p);
	if(dt12<=0&&d1p12>=0&&d2p12<=0)
		return 3; // rt==0;
	
	vec_crop(v20,v2p,nm20p);
	double dt20=vec_dotp(nm012,nm20p);
	if(dt20<=0&&d2p20>=0&&d0p20<=0)
		return 4; // rt==0;

	if(dt01>0&&dt12>0&&dt20>0) return 6;
	else jf_error("asdf posiotin");
}
void PointProjectToPlane(double pa[3],double d ,double norm[3] ,double pj[3] ){

  int i ;
  double dp ; // vec_uni( norm ) ;

  dp=norm[0]*pa[0]+norm[1]*pa[1]+norm[2]*pa[2] +d ;
  for(i=0 ; i<3 ; i++ )
    pj[i]=pa[i]-dp*norm[i] ;
}

void getTheClosestPointAt3DSeg(double p[3],double p0[3],double p1[3],double pcha[3]){

	double ratio,dp1,dp2,vp0p[3],vp0p1[3],vp1p[3];

	vec_2p(p0,p,vp0p);
	vec_2p(p0,p1,vp0p1);
	if((dp1=vec_dotp(vp0p,vp0p1))<=0){
		pcha[0]=p0[0];
		pcha[1]=p0[1];
		pcha[2]=p0[2];
		return;
	}
	vec_2p(p1,p,vp1p);
	if((dp2=vec_dotp(vp1p,vp0p1))>=0){
		pcha[0]=p1[0];
		pcha[1]=p1[1];
		pcha[2]=p1[2];
		return;
	}
	ratio=dp1/(dp1-dp2);
	pcha[0]=p0[0]+ratio*vp0p1[0];
	pcha[1]=p0[1]+ratio*vp0p1[1];
	pcha[2]=p0[2]+ratio*vp0p1[2];
}

void getClosestPointOnTriangle(double p[3],double p0[3],double p1[3],double p2[3],double pnear[3]){

	int ip=positionOfPointProjectToTri(p,p0,p1,p2);
	switch(ip){
		case 0: 
			memcpy(pnear,p0,3*sizeof(double));
			return;
		case 1: 
			memcpy(pnear,p1,3*sizeof(double));
			return;
		case 2: 
			memcpy(pnear,p2,3*sizeof(double));
			return;
		case 3:
			getTheClosestPointAt3DSeg(p,p1,p2,pnear);
			return;
		case 4:
			getTheClosestPointAt3DSeg(p,p2,p0,pnear);
			return;
		case 5:
			getTheClosestPointAt3DSeg(p,p0,p1,pnear);
			return;
		case 6:
			double norm[3];
			norm_3p(p0,p1,p2,norm);
		    if(vec_uni(norm)==0) 
				jf_error( " failed since rounderror change sizes try again") ;
			double d0=-vec_dotp(norm,p0);
			PointProjectToPlane(p,d0,norm,pnear);
	}
}


void BoxOf3p( double p1[3] ,double p2[3] ,double p3[3] , double box[6] ){

  //int i ;
  //double a ;

  box[0]=min( p1[0] , min(p2[0] , p3[0])) ;
  box[1]=min( p1[1] , min(p2[1] , p3[1])) ;
  box[2]=min( p1[2] , min(p2[2] , p3[2])) ;
  box[3]=max( p1[0] , max(p2[0] , p3[0])) ;
  box[4]=max( p1[1] , max(p2[1] , p3[1])) ;
  box[5]=max( p1[2] , max(p2[2] , p3[2])) ;
  //a=Max3( box[1]-box[0] ,box[3]-box[2] ,box[5]-box[4] ) ;
  //for( i=0 ; i<3 ; i++ ){
    //box[2*i] -= 0.01*a ;
    //box[2*i+1] += 0.01*a ;
  //}
}
//	OutPutOff((vertcoord,numvert,trips,numtri);
PointInPolyhedron::PointInPolyhedron(double (*vti)[3], int numvi,int (*tris)[3],int (*twoma)[2],int numti,int nummat){//,double epsion=0){

//	epsilonon=epsion;
	clock_t ts=clock();
	if(nummat==0) nummat=1;
	numvert=numvi;
	vertcoord=(double (*)[3]) new double[3*numvert];
	memcpy(vertcoord,vti,sizeof(double)*3*numvert);
	numtri=numti;
	trips=(int (*)[3])new int[3*numtri];
	memcpy(trips, tris,sizeof(int)*3*numtri);
//	vertcoord=vti; trips=tris;
	twomaterialsoftri=(int (*)[2]) new int[2*numtri];
	if( twoma){ memcpy(twomaterialsoftri,twoma,sizeof(int)*2*numtri);}
	else{
		for(int i=0; i<numtri; i++){
			twomaterialsoftri[i][0]=-1;
			twomaterialsoftri[i][1]=0; //1;
		}
	}	
	tneighb=(int (*)[3]) new int [3*numtri];//* for multiRCT
	trisort=new int[numtri];
	numtriofnode=new int[numvert];
	tripositionofnode=new int[numvert];
	trilist=new int[3*numtri];
	triofnode=(int *) new int[numvert];
	sortTrianglesForMultiMaterial(nummat);
//	MeshtxtOut(vertcoord,numvert,trips,numtri);
	delete [] tneighb;
	delete [] trisort;
	vertattrib= new int [numvert];
	for(int i=0; i<numvert; i++)
		vertattrib[i]=-3;
	void **wvti;
	wrapPointsUpasVerts(wvti);
	polytree=new Kodtree(wvti,numvert,pofvforcoordnodes3,2,epsoverlap);
	delete [] wvti;
    polytree->setFuncExinfoShouldbeInCell(ifexinfoshouldbeincell);
	polytree->setFuncExinfoOverlapBox(ifexinfooverlapbox);
	for(int i=0; i<numtri; i++)
		polytree->insertExinfo(i+startaddress,1);
	setGCellAttribOfSubTree(polytree->getRoot());
	vertused=new bool[numvert];
	for(int i=0; i<numvert; i++)
		vertused[i]=false;
	triused=new  int[numtri];
	for(int i=0; i<numtri; i++)
		triused[i]=0;
	if(nummat==1){	
		for(int i=0; i<numvert; i++)
			if(numtriofnode[i]<3){ 
				cout<<i<<" a vertex with less than 3 triangles,check your model "<<numtriofnode[i]<<endl; exit(0);
			}
		clock_t te=clock();	cout<<"preprocess :"<<(te-ts)/1000.<<" sec."<<endl;return;
	}
	int i;
	for(i=0; i<numvert; i++)  //record numbers of triangles around each node
		numtriofnode[i]=0;
	for(i=0; i<numtri; i++){
		for(int j=0; j<3; j++)
			numtriofnode[trips[i][j]]++;
	}
	tripositionofnode[0]=0;
	for( i=1; i<numvert; i++)  //positions of first triangles in the trilist for each node
		tripositionofnode[i]=tripositionofnode[i-1]+numtriofnode[i-1];
	for( i=0; i<numtri; i++){
		for(int j=0; j<3; j++){
			trilist[tripositionofnode[ trips[i][j] ]]=i;
			tripositionofnode[ trips[i][j] ]++;
		}
	}
	tripositionofnode[0]=0;
	for( i=1; i<numvert; i++) //recover the first positions
		tripositionofnode[i]=tripositionofnode[i-1]+numtriofnode[i-1];
	for(int i=0; i<numvert; i++)
		triofnode[i]=trilist[tripositionofnode[i]];
	for(int i=0; i<numvert; i++)
		if(numtriofnode[i]<3){ 
			cout<<i<<" a vertex with less than 3 triangles,check your model "<<numtriofnode[i]<<endl; exit(0);
		}
	clock_t te=clock();
	cout<<"preprocess :"<<(te-ts)/1000.<<" sec."<<endl;
}
extern void getBodyOf1Material(int mt,double (*vertall)[3],int numvertall,int (*trisall)[3],int numtriall,int (*twomoftri)[2],
					   double (*&vert)[3],int &numvert,int (*&tris)[3],int &numtri,int *&gtriindof1mat,int *&ltogofvert);
void getBodyOf1Material(int mt,double (*vertall)[3],int numvertall,int (*trisall)[3],int numtriall,int (*twomoftri)[2],
					   double (*&vert)[3],int &numvert,int (*&tris)[3],int &numtri,int *&gtriindof1mat,int *&ltogofvert){

   int *gtolofvert=new int[numvertall];
   gtriindof1mat=new int[numtriall];
   ltogofvert=new int[numvertall];
   for(int i=0; i<numvertall; i++) gtolofvert[i]=-1;
   int tricount=0,vertcount=0;
   for(int i=0; i<numtriall; i++){
	   if(twomoftri[i][0]!=mt&&twomoftri[i][1]!=mt) continue;
	   gtriindof1mat[tricount++]=i;
	   for(int j=0; j<3; j++){
		   int nv=trisall[i][j];
		   if(gtolofvert[nv]!=-1) continue;
		   gtolofvert[nv]=vertcount;
		   ltogofvert[vertcount++]=nv;
	   }
   }
   vert=(double (*)[3]) new double[3*vertcount];
   tris=(int (*)[3]) new int[3*tricount];
   numvert=vertcount;
   numtri=tricount;
   for(int i=0; i<numvertall; i++){
	   int nv=gtolofvert[i];
	   if(nv!=-1)
		  memcpy( vert[nv],vertall[i],sizeof(double)*3);
   }
   for(int i=0; i<numtri; i++){
	   int nt=gtriindof1mat[i];
	   for(int j=0; j<3; j++)
		   tris[i][j]=gtolofvert[trisall[nt][j]];
//	   if(twomoftri[nt][1]==mt) swap(tris[i][0],tris[i][1]);
   }
   delete [] gtolofvert;
}


void PointInPolyhedron::sortTrianglesForMultiMaterial(int nummat){

	double (*vertb)[3];
	int numvertb,numtrib,(*trisb)[3],*gtriindof1mat,*ltogofvert;
	if(nummat==1){
		sortTrianglesOuterNormAndRecNeighb(vertcoord,numvert,trips,numtri);
		return;
	}
	for(int i=0; i<nummat; i++){
		getBodyOf1Material(i,vertcoord,numvert,trips,numtri,twomaterialsoftri,vertb,numvertb,trisb,numtrib,gtriindof1mat,ltogofvert);
		sortTrianglesOuterNormAndRecNeighb(vertb,numvertb,trisb,numtrib);
		for(int j=0; j<numtrib; j++){
			int nt=gtriindof1mat[j];
			for(int k=0; k<3; k++)
				trips[nt][k]=ltogofvert[trisb[j][k]];
			twomaterialsoftri[nt][0]=twomaterialsoftri[nt][0]+twomaterialsoftri[nt][1]-i;
			twomaterialsoftri[nt][1]=i;
		}
		delete [] vertb; delete [] trisb;
		delete [] gtriindof1mat; delete [] ltogofvert;
	}
}

void PointInPolyhedron::setGCellAttribOfSubTree(CellNode *pcell){

	if(!pcell) return;
	if(!pcell->isLeaf())
		for(int i=0; i<2; i++)
			setGCellAttribOfSubTree(pcell->child[i]);
	else if(pcell->lpwpinfo!=0||pcell->numvert!=0)
		pcell->inoutattrib=-2;
	else
		pcell->inoutattrib=-3;
}
void PointInPolyhedron::resetAttrib(void){
	for(int i=0; i<numvert; i++)
		vertattrib[i]=-3;
	setGCellAttribOfSubTree(polytree->getRoot());
	for(int i=0; i<numvert; i++)
		vertused[i]=false;
	for(int i=0; i<numtri; i++)
		triused[i]=0;
}
PointInPolyhedron::~PointInPolyhedron(){

	delete [] vertcoord;
	delete [] trips;
	delete [] twomaterialsoftri;
	delete [] vertattrib;
	delete [] triofnode;
	delete [] numtriofnode;
	delete [] tripositionofnode;
	delete [] trilist;
	delete [] vertused;
//	delete [] tneighb;
	delete [] triused;
	delete polytree; //->freeSubTree(polytree->getRoot());
}


void PointInPolyhedron::getCellSeqWithUnknownAttribFromaCell(CellNode *cnode,vector<CellNode *>  * &pcellseq,
														CellNode * &pcellm,int &ia,double pm[3]){

	CellNode *pcellt,*pcell,*pancestor;

	if(cnode==0) return;
	pcellseq=new vector<CellNode *>;
	pcellseq->push_back(cnode);
	pcellt=cnode;
	pm[0]=pcellt->bound[3]; pm[1]=pcellt->bound[4];pm[2]=pcellt->bound[5];
	for(;;){
////		pcell=getTheNeighbOfCellAtSpeciDirectWithRefPoint(pcellt,-1,0,pm);
		pm[0]=pcellt->bound[3]+polytree->getEpsCell();
		if((pancestor=polytree->findTheNearestAncestorContainingPoint(pcellt,pm))==0)
			pcell=0;
		else pcell=polytree->findaLeafCellContainingPoint(pancestor,pm);

		if(pcell==0){
			ia=-1; pcellm=0;
			return;
		}else if(pcell->inoutattrib!=-3){
			ia=pcell->inoutattrib; pcellm=pcell;
			return;
		}
		pcellseq->push_back(pcell);
		pcellt=pcell;
	}
}
int PointInPolyhedron::getDetermingTriFromVertOfaVisibleTri(double p[3], int tri,int vert){
	int numbv, *neighbverts;
	int vertridge;
	double maxcosa,mincosa;
//	int vert=trips[tri][nentity];

	getVertsAroundaVertMulti(vert,neighbverts,numbv);
	getThePointFormingLeastAngleWith2Points(p,vert,neighbverts,numbv,maxcosa,mincosa,vertridge);
	delete neighbverts;
	if(maxcosa>epscoplanar)
		jf_error("getdetermingtrifromvert");
	if(mincosa>-epscoplanar){vertattrib[vert]=-2;return tri;}
	int rt=getDetermingTriFromEdgeOfaVisibleTri(p,vert,vertridge);
	return rt;
}
int PointInPolyhedron::getDetermingTriFromEdgeOfaVisibleTri(double p[3],int va,int vb){

	double vavb[3],vap[3],normapb[3],normaxb[3],vax[3];
	vec_2p(vertcoord[va],vertcoord[vb],vavb);
//	vec_uni(vavb);
	vec_2p(vertcoord[va],p,vap);
//	vec_uni(vap);
	vec_crop(vap,vavb,normapb);
	vec_uni(normapb);
	double maxcos2fangle=-1.,cos2fa;
	int tri=-1;
	for(int i=0; i<numtriofnode[va]; i++){
		int ip=tripositionofnode[va]+i;
		int ctri=trilist[ip];
		if(trips[ctri][0]!=vb&&trips[ctri][1]!=vb&&trips[ctri][2]!=vb) continue;
		int vx=trips[ctri][0]+trips[ctri][1]+trips[ctri][2]-va-vb;
		vec_2p(vertcoord[va],vertcoord[vx],vax);
//		vec_uni(vax);
		vec_crop(vax,vavb,normaxb);
		vec_uni(normaxb);
		if((cos2fa=vec_dotp(normaxb,normapb))>maxcos2fangle){maxcos2fangle=cos2fa;tri=ctri;}
	}
	return tri;
}
void CompRatioPoint( double p0[3] ,double p1[3] ,double dlt,double p[3] ){

  int i ;

  for( i=0 ; i<3 ; i++ )
    p[i]=p0[i]+dlt*(p1[i]-p0[i]) ;

}
double sqDistPointToRatioPoint(double p[3],double p0[3],double p1[3],double dlt){
	double prt[3];
	CompRatioPoint(p0,p1,dlt,prt);
	return SqDistance3D(p,prt);
}
extern double sqDistPointToTritemp(double p[3],double p0[3],double p1[3],double p2[3]);
//double sqDistPointToTri(double p[3],double p0[3],double p1[3],double p2[3]){
//	sqDistPointToTritemp(p,p0,p1,p2);
//	return sqDistPointToTritemp(p,p0,p1,p2);
//}
double sqDistPointToTri(double p[3],double p0[3],double p1[3],double p2[3]){//,double sqdist0){
	extern double signed_distance(double p[3],double p0[3],double p1[3],double p2[3]);
	//return signed_distance(p,p0,p1,p2);

	double v0p[3],v20[3],v01[3];
	vec_2p(p0,p,v0p);
	vec_2p(p2,p0,v20);
	vec_2p(p0,p1,v01);
//double nm012[3];
//vec_crop(v20,v01,nm012);
//double sqdptoinner=vec_dotp(nm012,v0p);
//sqdptoinner*=sqdptoinner/vec_sqval(nm012);
//if(sqdptoinner>=sqdist0) return sqdptoinner;

	double d0p20=vec_dotp(v0p,v20);
	double d0p01=vec_dotp(v0p,v01);
	if(d0p20>=0&&d0p01<=0) return SqDistance3D(p,p0);
	
	double v1p[3],v12[3];
	vec_2p(p1,p,v1p);
	vec_2p(p1,p2,v12);
	double d1p01=vec_dotp(v1p,v01);
	double d1p12=vec_dotp(v1p,v12);
	if(d1p01>=0&&d1p12<=0) return SqDistance3D(p,p1);

	double v2p[3];
	vec_2p(p2,p,v2p);
	double d2p12=vec_dotp(v2p,v12);
	double d2p20=vec_dotp(v2p,v20);
	if(d2p12>=0&&d2p20<=0) return SqDistance3D(p,p2);

	double 
		nm012[3],
		nm01p[3],nm12p[3],nm20p[3];
	vec_crop(v20,v01,nm012);

	vec_crop(v01,v0p,nm01p);
	double dt01=vec_dotp(nm012,nm01p);
	if(dt01<=0&&d0p01>=0&&d1p01<=0)
		return sqDistPointToRatioPoint(p,p0,p1,d0p01/(d0p01-d1p01)); 
//		return sqDistPointToSeg3D(p,p0,p1); // rt==0;
		
	vec_crop(v12,v1p,nm12p);
	double dt12=vec_dotp(nm012,nm12p);
	if(dt12<=0&&d1p12>=0&&d2p12<=0)
		return sqDistPointToRatioPoint(p,p1,p2,d1p12/(d1p12-d2p12)); 
//		return sqDistPointToSeg3D(p,p1,p2);
	
	vec_crop(v20,v2p,nm20p);
	double dt20=vec_dotp(nm012,nm20p);
	if(dt20<=0&&d2p20>=0&&d0p20<=0)
		return sqDistPointToRatioPoint(p,p2,p0,d2p20/(d2p20-d0p20)); 
//		return sqDistPointToSeg3D(p,p2,p0);

	if(dt01>=0&&dt12>=0&&dt20>=0){
		double a=vec_dotp(nm012,v0p);
		return a*a/vec_sqval(nm012);
//		return sqdptoinner;
	}else
		jf_error("asdf posiotin");
}

void PointInPolyhedron::getVertsAroundaVertMulti(int v, int * &nbverts,int &numnbv){

	int count=0;
	nbverts=(int *) new int[numtriofnode[v]];
	for(int i=0; i<numtriofnode[v]; i++){
		int ip=tripositionofnode[v]+i;
		int ctri=trilist[ip];
		for(int j=0; j<3; j++){
			int cv=trips[ctri][j];
			if(cv==v||vertused[cv]) continue;
			nbverts[count++]=cv;
			vertused[cv]=true;
		}
	}
	numnbv=count;
	for(int i=0; i<numnbv; i++)
		vertused[nbverts[i]]=false;
}
int PointInPolyhedron::indexOfVertAtTri(int v, int tri[3]){

	if(tri[0]==v) return 0;
	else if(tri[1]==v) return 1;
	else if(tri[2]==v) return 2;
	else 
		jf_error("indexoftri");
}

int PointInPolyhedron::nextVertOfTri(int tri,int v){
	if(v==trips[tri][0]) return trips[tri][1];
	else if(v==trips[tri][1]) return trips[tri][2];
	else if(v==trips[tri][2]) return trips[tri][0];
	else jf_error("nextvoftri");
}
void PointInPolyhedron::getThePointFormingLeastAngleWith2Points(double p[3],int v, int *nbverts,int numnbv,double &maxcosa,double &mincosa,int &vridge){

	double pv[3],pvp[3],pvpi[3],dp;
	maxcosa=-1.;
	mincosa=1.;
	copy3DPoint(vertcoord[v],pv);
	vec_2p(pv,p,pvp);
	vec_uni(pvp);
	for(int i=0; i<numnbv; i++){
		vec_2p(pv,vertcoord[nbverts[i]],pvpi);
		vec_uni(pvpi);
		dp=vec_dotp(pvpi,pvp);
		if(dp>maxcosa){
			if(dp>epscoplanar)
				dp=dp;
			maxcosa=dp;
			vridge=nbverts[i];
		}
		if(dp<mincosa) mincosa=dp;
	}
}
double sqDistPointToSeg3D(double p[3],double p0[3],double p1[3]){

	double vp0p[3],vp0p1[3],vp1p[3];

	vec_2p(p0,p,vp0p);
	vec_2p(p0,p1,vp0p1);
	if(vec_dotp(vp0p,vp0p1)<=0)
		return SqDistance3D(p0,p);
	vec_2p(p1,p,vp1p);
	double prjp1p=vec_dotp(vp1p,vp0p1);
	double sqdp1p=SqDistance3D(p1,p);
	if(prjp1p>=0)
		return sqdp1p;
	double sqd=sqdp1p-prjp1p*prjp1p/vec_sqval(vp0p1);
	if(sqd<0){
		cout<<sqd<<" less than 0"<<endl;
		sqd=0;
	}
	return sqd;  //?
}

void  PointInPolyhedron::getEdgeOfTri(int np[3], int index, int &a, int &b){
	if(index==0){
		a=np[1]; b=np[2];
	}else	if(index==1){
		a=np[2]; b=np[0];
	}else	if(index==2){
		a=np[0]; b=np[1];
	}else
		jf_error("error getedgeoftri");
}
void  PointInPolyhedron::sortTrianglesOuterNormAndRecNeighb(double (*vertb)[3],int numvertb,int (*tripsb)[3],int numtrib){

	int i;
	for(i=0; i<numvertb; i++)  //record numbers of triangles around each node
		numtriofnode[i]=0;
	for(i=0; i<numtrib; i++){
		for(int j=0; j<3; j++)
			numtriofnode[tripsb[i][j]]++;
	}
	tripositionofnode[0]=0;
	for( i=1; i<numvertb; i++)  //positions of first triangles in the trilist for each node
		tripositionofnode[i]=tripositionofnode[i-1]+numtriofnode[i-1];
	for( i=0; i<numtrib; i++){
		for(int j=0; j<3; j++){
			trilist[tripositionofnode[ tripsb[i][j] ]]=i;
			tripositionofnode[ tripsb[i][j] ]++;
		}
	}
	tripositionofnode[0]=0;
	for( i=1; i<numvertb; i++) //recover the first positions
		tripositionofnode[i]=tripositionofnode[i-1]+numtriofnode[i-1];
	for(int i=0; i<numvertb; i++)
		triofnode[i]=trilist[tripositionofnode[i]];
	for(i=0; i<numtrib; i++){
		trisort[i]=0;
		for(int j=0; j<3; j++)
			tneighb[i][j]=-1;
	}
	int tri,shellcount=0;
	while((tri=getAndSortaLowestTri(vertb,numvertb,tripsb,++shellcount==1))!=-1)
		sort1ShellFromaTri(tripsb,tri);

}
	
int  PointInPolyhedron::getAndSortaLowestTri(double (*vertb)[3],int numvertb,int (*tripsb)[3],bool firstshell){

	double z0=std::numeric_limits<double>::max();
	int vert=-1;

	for(int i=0; i<numvertb; i++){
		if(trisort[triofnode[i]]==1) continue;
		if(vertb[i][2]<z0){ z0=vertb[i][2]; vert=i;}
	}
	if(vert==-1) return -1;

	int vertridge=-1;
	int tria,trib;

	z0=1.;
	for(int i=0; i<numtriofnode[vert]; i++){
		int ip=tripositionofnode[vert]+i;
		int tri=trilist[ip];
		for(int j=0; j<3; j++){
			int nd=tripsb[tri][j];
			if(nd==vert) continue;
			double z,p[3];
			vec_2p(vertb[vert],vertb[nd],p);
			if((z=p[2]*p[2]/(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]))<z0){
				z0=z; vertridge=nd;
			}
		}
	}
	if(vertridge==-1) jf_error("err getand sortalow");
	get2TriCom2NodesWithoutTopology(tripsb,vert,vertridge,tria,trib);
	int nd=tripsb[trib][0]+tripsb[trib][1]+tripsb[trib][2]-vertridge-vert;
	double vol=VolumOf4p(vertb[tripsb[tria][0]],vertb[tripsb[tria][1]],vertb[tripsb[tria][2]],vertb[nd]);
	bool orientflag=false;
	if(fabs(vol)<=0.000000000001){// exact computation should used in the future at here.
		double norm[3];
		{ double nm1[3],nm2[3];
		  int nc=tripsb[tria][0]+tripsb[tria][1]+tripsb[tria][2]-vertridge-vert;
		  norm_3p(vertb[vert],vertb[vertridge],vertb[nc],nm1);
		  norm_3p(vertb[vert],vertb[vertridge],vertb[nd],nm2);
		  if(vec_dotp(nm1,nm2)>=0) jf_error("err getandsortalow");
		}
		norm_3p(vertb[tripsb[tria][0]],vertb[tripsb[tria][1]],vertb[tripsb[tria][2]],norm);
		if(norm[2]<0) orientflag=true;
	}else
		if(vol<0) orientflag=true;
	if(!firstshell) orientflag=!orientflag;
	if(!orientflag)
		swap(tripsb[tria][0],tripsb[tria][1]);
	return tria;
}

void  PointInPolyhedron::sort1ShellFromaTri(int (*tripsb)[3],int tri){

	std::queue<int> quet;
	quet.push(tri);
	trisort[tri]=1;
	while(!quet.empty()){
		int ctri=quet.front();
		quet.pop();
		for(int i=0;i<3;i++){
			if(tneighb[ctri][i]>=0) continue;
			int tnb=getNeighbTriWithoutTopology(tripsb,ctri,i);
			if(trisort[tnb]==0){
				if(!triSortAs2Nodes(tripsb[tnb],tripsb[ctri][(i+2)%3],tripsb[ctri][(i+1)%3]))
					swap(tripsb[tnb][0],tripsb[tnb][1]);
				trisort[tnb]=1;
				quet.push(tnb);
			}
			tneighb[ctri][i]=tnb;
			int ind=indexOfVertAtTri(tripsb[ctri][(i+1)%3],tripsb[tnb]);
			tneighb[tnb][(ind+1)%3]=ctri;
		}
	}
}
int  PointInPolyhedron::getNeighbTriWithoutTopology(int (*tripsb)[3],int tri,int ind){

	int a=tripsb[tri][(ind+1)%3];
	int b=tripsb[tri][(ind+2)%3];
	for(int i=0; i<numtriofnode[a]; i++){
		int ip=tripositionofnode[a]+i;
		int ctri=trilist[ip];
		if(ctri==tri)continue;
		if(tripsb[ctri][0]==b||tripsb[ctri][1]==b||tripsb[ctri][2]==b)
			return ctri;
	}
	return -1;
}
bool  PointInPolyhedron::triSortAs2Nodes(int tri3p[3],int va, int vb){

	if(tri3p[0]==va&&tri3p[1]==vb||
		tri3p[1]==va&&tri3p[2]==vb||
		tri3p[2]==va&&tri3p[0]==vb)
		return true;
	else
		return false;
}
void  PointInPolyhedron::get2TriCom2NodesWithoutTopology(int (*tripsb)[3],int va, int vb, int &ta, int &tb){

	ta=tb=-1;
	for(int i=0; i<numtriofnode[va]; i++){
		int ip=tripositionofnode[va]+i;
		int tri=trilist[ip];
		if(tripsb[tri][0]==vb||tripsb[tri][1]==vb||tripsb[tri][2]==vb){
			if(ta==-1)
				ta=tri;
			else{
				tb=tri;
				return;
			}
		}
	}
	if(tb==-1)	jf_error("get2triwith");
}

