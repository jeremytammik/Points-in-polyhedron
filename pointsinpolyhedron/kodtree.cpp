#include "kodtree.h"
#include <queue>
#include "mex.h"
typedef double Point[3];
extern void jf_error(char *ch);
extern bool ifBoxContainPoint( Point p,const Box &bound,const Box &rootbound);
extern bool if2BoxOverlap(const Box &a,const Box &b);
extern bool if2BoxNeighb(const Box &a,const Box &b);
extern bool ifPointOverlapWithBox(const Point &p,const Box &bd,const Box &rootbound,double eps);
extern void boxOfVerts( void **v , int num ,Box b ,void (*Funcpointofv)(Point p,void *v));
extern void copy3DPoint(const Point &pfr,Point &pto); //&?
extern int IsTriangleBoxInt(Point p1 ,Point p2 ,Point p3 , double box[6] );
//const double Kodtree::epsilonon=0.00000000000001;
//double Kodtree:: epsoverlap;
//const double Kodtree:: epscoplanar=0.000001;
//int Kodtree::cellcapacity;

 CellNode * Kodtree:: findTheNearestAncestorContainingPoint(CellNode *pcell0,Point pcha){
	
	CellNode *pcell=pcell0;
	for(;;){
		if(pcell==0) return 0;
		if(ifBoxContainPoint(pcha,pcell->bound,root->bound)) return pcell;//need or not to follow the convention?
		else pcell=pcell->parent;
	}
}

CellNode *  Kodtree::findaLeafCellContainingPoint(CellNode *pcell,Point p){

	CellNode *rtpcell;

	if(!pcell||!ifBoxContainPoint(p,pcell->bound,root->bound))//need or not to follow the convention?
		return 0;
	if(pcell->isLeaf())
		return pcell;
	for( int i=0; i<2; i++)
		if((rtpcell=findaLeafCellContainingPoint(pcell->child[i],p))!=0)
			return rtpcell;
	jf_error("err findaleafcellcontainp");
}

bool  Kodtree::if2CellNeighb(CellNode *pcell0, CellNode *pcell1){

	if(!pcell0||!pcell1) 
		jf_error("err is2cellneigh");
	if(if2BoxNeighb(pcell0->bound,pcell1->bound))
		return true;
	else
		return false;
}


Kodtree::Kodtree(const Box &bd,Funcpointofvert pofvin,int capacity,double epsi){

	double lcube=max(bd[3]-bd[0],max(bd[4]-bd[1],bd[5]-bd[2]));
	epscell=numeric_limits<double>::epsilon()*(1+lcube*10.);
	root=new CellNode(bd);
	pofv=pofvin;
	cellcapacity=capacity;
	epsoverlap=epsi;
}

Kodtree::Kodtree(void **vert, int numvert,Funcpointofvert pofvin,int capacity,double epsi){

	Box bd;
	boxOfVerts(vert,numvert,bd,pofvin);
	double lcube=max(bd[3]-bd[0],max(bd[4]-bd[1],bd[5]-bd[2]));
	epscell=numeric_limits<double>::epsilon()*(1+lcube*10.);
	root=new CellNode(bd);
	pofv=pofvin;
	cellcapacity=capacity;
	epsoverlap=epsi;
	for( int i=0; i<numvert; i++)
		insertVert(vert[i]);
}

Kodtree::Kodtree(void **vert, int numvert,const Box &bd,Funcpointofvert pofvin,int capacity,double epsi){

//	Box bd;
//	boxOfPoints(vert,numvert,bd);
	double lcube=max(bd[3]-bd[0],max(bd[4]-bd[1],bd[5]-bd[2]));
	epscell=numeric_limits<double>::epsilon()*(1+lcube*10.);
	root=new CellNode(bd);
	pofv=pofvin;
	cellcapacity=capacity;
	epsoverlap=epsi;
	for( int i=0; i<numvert; i++)
		insertVert(vert[i]);
}



Kodtree::~Kodtree(){

	freeSubTree(root);
}



void Kodtree::insertWpVertInSubTree(const Point &p, WpVert *v, CellNode *cnode){

	if(!cnode)
		jf_error("err insvinst");
	if(!ifPointOverlapWithBox(p,cnode->bound,root->bound ,epsoverlap ))
		return ;
	if(!cnode->isLeaf()){
		for(int i=0; i<2; i++)
			insertWpVertInSubTree(p,v,cnode->child[i]);
		cnode->numvert++;
		return ;
	}
	if(!cnode->vert){
		cnode->vert =(PtWpVert *) new PtWpVert[Kodtree::cellcapacity];
		cnode->nodecapacity=Kodtree::cellcapacity;
	}else if(cnode->numvert==cnode->nodecapacity&&!canNodeSplit(cnode)){
		cnode->vert=renew(cnode->vert,cnode->nodecapacity,2*cnode->nodecapacity);
		cnode->nodecapacity*=2;
	}
	if(cnode->numvert<cnode->nodecapacity){
		cnode->vert[cnode->numvert++]=v;
		v->rcount ++;
	}else{
		splitNode(cnode);
		for(int i=0; i<2; i++)
			insertWpVertInSubTree(p,v,cnode->child[i]);
		cnode->numvert++;
		return;
	}
}


bool Kodtree::isVertRecordInSubTree(const Point &p,void *v, CellNode *cnode){

	if(!cnode)
		jf_error("err insvinst");
	if(cnode->numvert<=0||!ifPointOverlapWithBox(p,cnode->bound,root->bound ,epsoverlap ))
		return false;
	if(!cnode->isLeaf()){
		for(int i=0; i<2; i++)
			if(isVertRecordInSubTree(p,v,cnode->child[i])) return true;
		return false;
	}
	if(!cnode->vert)
		jf_error("err insvinst");
	for(int i=0; i<cnode->numvert; i++)
		if(cnode->vert[i]->vt==v) return true;
	return false;
}

int comWpVertNum(CellNode *cnode, CellNode *cnsib){

	int num=0;
	for(int i=0; i<cnsib->numvert; i++){
		for(int j=0; j<cnode->numvert; j++)
			if(cnsib->vert[i]==cnode->vert[j]){num++; break;}
	}
	return num;
}

void Kodtree::deleteVertInSubTree(const Point &p,void *v, CellNode *cnode){

	if(!cnode)
		jf_error("err insvinst");
	if(!ifPointOverlapWithBox(p,cnode->bound,root->bound ,epsoverlap ))
		return ;
	cnode->numvert--;
	if(!cnode->isLeaf()){
		for(int i=0; i<2; i++)
			deleteVertInSubTree(p,v,cnode->child[i]);
		return ;
	}
	if(!cnode->vert)
		jf_error("err deletevertinsubtree");
	int i;
	for(i=0; i<cnode->numvert; i++)
		if(cnode->vert[i]->vt==v)
			break;
	if(--(cnode->vert[i]->rcount)<=0) delete cnode->vert[i];
	if(i!=cnode->numvert)
		cnode->vert[i]=cnode->vert[cnode->numvert];
	if(cnode->numvert==0){
		delete cnode->vert;
		cnode->vert=0;
	}
	
}

void Kodtree::insertWpInfoInSubTree(WpInfo *pwinfo, CellNode *cnode){

	if(!cnode)
		jf_error("err insvinst");
	if(!ifExinfoOverlapBox(pwinfo->info,pwinfo->infotype ,cnode->bound,epsoverlap ))
		return ;
	if(!cnode->isLeaf()){
		for(int i=0; i<2; i++)
			insertWpInfoInSubTree(pwinfo,cnode->child[i]);
		return ;
	}
	if(!ifExinfoShouldbeInCell(pwinfo->info,pwinfo->infotype ,cnode ))
		return ;
	if(!cnode->lpwpinfo){
		cnode->lpwpinfo=new std::list<WpInfo *>;
	}
	cnode->lpwpinfo->push_back(pwinfo);
	pwinfo->rcount++;
}

void Kodtree::deleteExinfoInSubTree(void *info,int infotype, CellNode *cnode){

	if(!cnode)
		jf_error("err insvinst");
	if(!ifExinfoOverlapBox(info,infotype ,cnode->bound,epsoverlap ))
		return ;
	if(!cnode->isLeaf()){
		for(int i=0; i<2; i++)
			deleteExinfoInSubTree(info,infotype,cnode->child[i]);
		return ;
	}
	if(!ifExinfoShouldbeInCell(info,infotype ,cnode ))
		return ;
	if(!cnode->lpwpinfo) return;
	std::list<WpInfo *>::iterator ite,iten;
	for( ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite=iten){
		iten=ite; iten++;
		if((*ite)->info==info&&(*ite)->infotype==infotype){
			if(--(*ite)->rcount<=0) delete *ite;
			cnode->lpwpinfo->erase(ite);
		}
	}
	if((cnode->lpwpinfo)->empty()){
		delete cnode->lpwpinfo;
		cnode->lpwpinfo=0;
	}
}



void Kodtree::checkAndMergeSubTreeAfterDelete(const Point &p,CellNode *cnode){

	if(!cnode||cnode->isLeaf()|| ! ifPointOverlapWithBox(p,cnode->bound,root->bound ,epsoverlap ))
		return;
	else if(cnode->numvert<=Kodtree::cellcapacity){
		mergeSubTree(cnode);
		checkAndRemoveSurplusWpInfoAfterMerge(cnode);
	}else
		for(int i=0; i<2; i++)
			checkAndMergeSubTreeAfterDelete(p,cnode->child[i]);
}


void Kodtree::checkAndRemoveSurplusWpInfoAfterMerge(CellNode *cnode){

	if(!cnode->lpwpinfo) return;
	std::list<WpInfo *>::iterator ite,iten;
	for( ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite=iten){
		iten=ite; iten++;
		if(!ifExinfoShouldbeInCell((*ite)->info,(*ite)->infotype,cnode)){
			if(--(*ite)->rcount<=0) delete *ite;
			cnode->lpwpinfo->erase(ite);
		}
	}
	if((cnode->lpwpinfo)->empty()){
		delete cnode->lpwpinfo;
		cnode->lpwpinfo=0;
	}
}


void Kodtree::mergeSubTree(CellNode *cnode){

	if(cnode==0) jf_error("err mergecellup");
	if(cnode->isLeaf()) return;
	for(int i=0; i<2; i++)
		mergeSubTree(cnode->child[i]);
	merge2SubCellWpVert(cnode);
	merge2SubCellWpInfo(cnode);
	for(int i=0; i<2; i++){
		delete cnode->child[i];
		cnode->child[i]=0;
	}
}

void Kodtree::merge2SubCellWpVert(CellNode *cnode){

	cnode->vert =(PtWpVert *) new PtWpVert[Kodtree::cellcapacity];
	if(cnode->isLeaf()) jf_error("err merge2subcellvert");
	for(int i=0; i<cnode->child[0]->numvert; i++){
		cnode->vert[i]=cnode->child[0]->vert[i];
		cnode->vert[i]->vget=true;
		cnode->vert[i]->rcount++;
	}
	int count=cnode->child[0]->numvert;
	for(int i=0; i<cnode->child[1]->numvert; i++){
		WpVert *v=cnode->child[1]->vert[i];
		if(v->vget ==false){cnode->vert[count++]=v; v->rcount++;}
//		else v->rcount--;
	}
	for(int i=0; i<count; i++)
		cnode->vert[i]->vget=false;
	if(cnode->numvert!=count) jf_error("err merge2subcellvert1");
}

void Kodtree::merge2SubCellWpInfo(CellNode *cnode){

	if(cnode->isLeaf()) jf_error("err merge2subcellwpinfo");
	CellNode *left=cnode->child[0],*right=cnode->child[1];
	if(left->lpwpinfo==0&&right->lpwpinfo==0){
		cnode->lpwpinfo=0;
		return;
	}
	if(left->lpwpinfo!=0){
		if(right->lpwpinfo!=0){
			for(std::list<WpInfo *>::iterator ite=left->lpwpinfo->begin();ite!=left->lpwpinfo->end(); ite++){
				(*ite)->get=true;
				//(*ite)->rcount++;
			}
			for(std::list<WpInfo *>::iterator iten, ite=right->lpwpinfo->begin();ite!=right->lpwpinfo->end(); ite=iten){
			  iten=ite, iten++;
			  if(!(*ite)->get){
				  //(*ite)->rcount++;
				  left->lpwpinfo->splice(left->lpwpinfo->end(),*(right->lpwpinfo),ite);
			  }
			}
		    for(std::list<WpInfo *>::iterator ite=left->lpwpinfo->begin();ite!=left->lpwpinfo->end(); ite++)
			  (*ite)->get=false;
		}
		cnode->lpwpinfo=left->lpwpinfo;
		left->lpwpinfo=0;
		//delete right->lpwpinfo;
		//right->lpwpinfo=0;
	}else{
		cnode->lpwpinfo=right->lpwpinfo;
		right->lpwpinfo=0;
	}
}



void Kodtree::collectVertsWithBox(const Box &bd, list<void *> &lvert){

	std::list<WpVert *> lwpvert;
	collectWpVertsWithBoxInSubTree(root,bd,lwpvert); // may be not unique.
	for(std::list<WpVert *>::iterator ite=lwpvert.begin();ite!=lwpvert.end(); ite++){
		lvert.push_back((*ite)->vt);
		(*ite)->vget=false;
	}
}

void Kodtree::collectWpVertsWithBoxInSubTree(CellNode *cnode,const Box &bd,list<WpVert *> &lvert){
	if(!cnode) return;
	if(!if2BoxOverlap(bd,cnode->bound)) return;
	if(!cnode->isLeaf()){
		collectWpVertsWithBoxInSubTree(cnode->child[0],bd,lvert);
		collectWpVertsWithBoxInSubTree(cnode->child[1],bd,lvert);
	}else{
		for(int i=0; i<cnode->numvert; i++){
			Point p;
			if(cnode->vert[i]->vget==true) continue;
			pofv(p,cnode->vert[i]->vt);
			if((p,bd,bd)){
				lvert.push_back(cnode->vert[i]); //maybe not unique.
				cnode->vert[i]->vget =true;
			}
		}
	}
}

void Kodtree::collectVertsWithCell(CellNode *cnode, std::vector<void *> &vecvert){

	for(int i=0; i<cnode->numvert; i++)
		vecvert.push_back(cnode->vert[i]->vt);
}


void Kodtree::collectExinfoWithBox(const Box &bd, int infotype,list<void *> &lexinfo){

	std::list<WpInfo *> lwpinfo;
	collectWpinfoWithBoxInSubTree(root,bd,infotype,lwpinfo);
	for(std::list<WpInfo *>::iterator ite=lwpinfo.begin();ite!=lwpinfo.end(); ite++){
		lexinfo.push_back((*ite)->info);
		(*ite)->get=false;
	}
}

void Kodtree::collectWpinfoWithBoxInSubTree(CellNode *cnode,const Box &bd,int infotype,list<WpInfo *> &lwpinfo){
	if(!cnode) return;
	if(!if2BoxOverlap(bd,cnode->bound)) return;
	if(!cnode->isLeaf()){
		collectWpinfoWithBoxInSubTree(cnode->child[0],bd,infotype,lwpinfo);
		collectWpinfoWithBoxInSubTree(cnode->child[1],bd,infotype,lwpinfo);
	}else{
		if(cnode->lpwpinfo==0) return;
		for(std::list<WpInfo *>::iterator ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite++){
			if((*ite)->infotype!=infotype||(*ite)->get==true)
				continue;
			if(ifExinfoOverlapBox((*ite)->info,infotype,bd,epsoverlap)){
				lwpinfo.push_back((*ite));
				(*ite)->get=true;
			}

		}
	}
}

void Kodtree::collectExinfoWithCell(CellNode *cnode, int infotype,list<void *> &lexinfo){

	if(cnode->lpwpinfo==0) return;
	for(std::list<WpInfo *>::iterator ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite++)
		if((*ite)->infotype==infotype)
			lexinfo.push_back((*ite)->info);
}
bool Kodtree::canNodeSplit(CellNode *cnode){
	return (cnode->bound[3]-cnode->bound[0])*100000.> root->bound[3]-root->bound[0];
}
void Kodtree::splitNode(CellNode *cnode){

	for(int i=0; i<2; i++){
		cnode->child[i]=new CellNode(cnode->bound);
		cnode->child[i]->parent=cnode;
	}
	int di;
	getTheLongestDistOfBox(cnode->bound,di);
	cnode->child[1]->bound[di]=cnode->child[0]->bound[di+3]=(cnode->bound[di]+cnode->bound[di+3])/2.;
//	if(cnode->vert==0)
//		return;
	for(int i=0; i<cnode->numvert; i++){
		Point p;
		pofv(p,cnode->vert[i]->vt);
		for(int j=0; j<2; j++)
			insertWpVertInSubTree(p,cnode->vert[i],cnode->child[j]);
	}
	for(int i=0; i<cnode->numvert ; i++)
		cnode->vert[i]->rcount--;
	delete [] cnode->vert;
	cnode->vert=0;
	if(cnode->lpwpinfo==0) return;
	for(std::list<WpInfo *>::iterator ite=cnode->lpwpinfo->begin();ite!=cnode->lpwpinfo->end(); ite++){
		(*ite)->rcount--;
		for(int i=0; i<2; i++)
			insertWpInfoInSubTree(*ite,cnode->child[i]);
	}
	delete cnode->lpwpinfo;
	cnode->lpwpinfo=0;
//	cnode->numvert=0;
}


void Kodtree::freeSubTree(CellNode *pcell){

	if(pcell ==0) return;
	for(int i=0; i<2; i++)
		freeSubTree(pcell->child[i]);
	delete pcell;
}


CellNode ::CellNode(const Box &bd){

//	psegar=0;
	vert=0;
	numvert=0;
	nodecapacity=0;
	lpwpinfo=0;
	inoutattrib=-3;
	for(int i=0; i<6; i++)
		bound[i]=bd[i];
	child[0]=child[1]=0;
	parent=0;
}

CellNode ::~CellNode (){
//	delete psegar;
	if(vert!=0)
	    for(int i=0; i<numvert; i++)
			if(--(vert[i]->rcount)<=0) delete vert[i];
	if(lpwpinfo!=0)
		for(std::list<WpInfo *>::iterator ite=lpwpinfo->begin();ite!=lpwpinfo->end(); ite++)
			if(--(*ite)->rcount<=0)	delete *ite;
	delete [] vert;
	delete lpwpinfo;
}

double sqdistPointToBox( Point p,const Box &bd){

	double a[3];
	double q=0;
	for(int i=0; i<3; i++){
		if(p[i]>bd[i+3]) a[i]=p[i]-bd[i+3];
		else if(p[i]<bd[i]) a[i]=bd[i]-p[i];
		else a[i]=0;
		q+=a[i]*a[i];
	}
	return q;
}
double sqdistInnerPointToBoxBound( Point p,const Box &bd){

	double a=min(p[0]-bd[0],bd[3]-p[0]);
	double b=min(p[1]-bd[1],bd[4]-p[1]);
	double c=min(p[2]-bd[2],bd[5]-p[2]);
	double d= min(c,min(a,b));
	return d*d;
}
void getTheLongestDistOfBox(const Box &b,int &di, double *pdist){

	di=0;
	double dist=0.;
	for(int i=0; i<3; i++)
		if(b[i+3]-b[i]>dist){
			dist=b[i+3]-b[i];
			di=i;
		}
	if(pdist!=0) *pdist=dist;
}
bool if2BoxOverlap(const Box &a,const Box &b){

	if(a[0]>b[3]||a[1]>b[4]||a[2]>b[5]||a[3]<b[0]||a[4]<b[1]||a[5]<b[2])
		return false;
	return true;
}


void copy3DPoint(const Point &pfr,Point &pto){ //const or not?

	pto[0]=pfr[0];
	pto[1]=pfr[1];
	pto[2]=pfr[2];
}

void jf_error(char *ch){

		 mexErrMsgTxt("Maybe input model is invalid or illegal.");
}


bool ifBoxContainPoint( Point p,const Box &bound,const Box &rootbound){

	if(p[0]<bound[0]||p[1]<bound[1]||p[2]<bound[2]||p[0]>bound[3]||p[1]>bound[4]||p[2]>bound[5])
		return false;
	else if(bound[0]!=rootbound[0]&&p[0]==bound[0]||
		    bound[1]!=rootbound[1]&&p[1]==bound[1]||
			bound[2]!=rootbound[2]&&p[2]==bound[2] ) //need or not to follow the convention?
		return false;
	else 
		return true;
}

bool ifPointOverlapWithBox(const Point &p,const Box &bd,const Box &rootbound,double eps){

	Box bound;

	double a[3];
	for(int i=0; i<3; i++)
		a[i]=bd[i+3]-bd[i];
	for(int i=0; i<3; i++){
		bound[i]=bd[i]-eps*a[i];
		bound[i+3]=bd[i+3]+eps*a[i];
	}//convention
	if(p[0]<bound[0]||p[1]<bound[1]||p[2]<bound[2]||p[0]>bound[3]||p[1]>bound[4]||p[2]>bound[5])
		return false;
	else if(bound[0]!=rootbound[0]&&p[0]==bound[0]||
		    bound[1]!=rootbound[1]&&p[1]==bound[1]||
			bound[2]!=rootbound[2]&&p[2]==bound[2] )
		return false;
	else 
		return true;
}

bool if2BoxNeighb(const Box &a,const Box &b){

	if(a[0]>b[3]||a[1]>b[4]||a[2]>b[5]||a[3]<b[0]||a[4]<b[1]||a[5]<b[2])
		return false;
	return true;
}

void boxOfVerts( void **v, int num ,Box box,void (*pofv)(Point p,void *v) ){

  int i , j ;
  double a ;

  Point p;
 
  pofv(p,v[0]);
  for(i=0 ; i<3 ; i++ ){
	  box[i]=box[i+3]=p[i] ;
  }
  for(j=1 ; j<num ; j++ ){
	  pofv(p,v[j]);
    for( i=0 ; i<3 ; i++ ){
	 if( p[i]<box[i] ) box[i]=p[i] ;
	 if( p[i]>box[i+3] ) box[i+3]=p[i] ;
    }
  }
  a=max( box[3]-box[0] ,max(box[4]-box[1],box[5]-box[2]) ) ;
  for( i=0 ; i<3 ; i++ ){
    box[i] -= 0.01*a ;
    box[i+3] += 0.01*a ; //keep unchanged for a undegenerate 3D box or use a,b and c?
  }
}


void vec_2p(double *pointa , double *pointb , double *vector)
{
  int i;
  for(i=0; i<3; i++)
		vector[i]=pointb[i]-pointa[i];
}
void
vec_neg(double *vector)
{
  vector[0]= -vector[0];
  vector[1]= -vector[1];
  vector[2]= -vector[2];
}
int
vec_uni(double *vector)
{
   double len;  int i;
   len=(double)sqrt(vector[0]*vector[0]+vector[1]*vector[1]
		+vector[2]*vector[2]);
   if(len<=0.00000000001) return(0);        /* eps_2 zero vector bound :min_siz/2. */
   for(i=0; i<3; i++)   vector[i]/=len;
   return(1);
}
double vec_val( double *vec )
{
  return( sqrt(vec[0] * vec[0] +vec[1] * vec[1] + vec[2] *vec[2] )) ;
}
double vec_sqval( double *vec )
{
  return vec[0] * vec[0] +vec[1] * vec[1] + vec[2] *vec[2]  ;
}
double
vec_dotp(double *vector1,double *vector2)
{
   return(vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2]);
}
void
vec_crop(double *vector1,double *vector2,double *vector3)
{
  vector3[0]=vector1[1]*vector2[2]-vector2[1]*vector1[2];
  vector3[1]=vector2[0]*vector1[2]-vector1[0]*vector2[2];
  vector3[2]=vector1[0]*vector2[1]-vector2[0]*vector1[1];
}
	
double
Distance3D( double *v1,double *v2){
  return( sqrt( (v1[0]-v2[0])*(v1[0]-v2[0]) +
			 (v1[1]-v2[1])*(v1[1]-v2[1]) +
			 (v1[2]-v2[2])*(v1[2]-v2[2]) )
	   ) ;
}
double
SqDistance3D( double *v1,double *v2){
  return( (v1[0]-v2[0])*(v1[0]-v2[0]) +
			 (v1[1]-v2[1])*(v1[1]-v2[1]) +
			 (v1[2]-v2[2])*(v1[2]-v2[2])
	   ) ;
}


double vec_blep( double *vec1 , double *vec2 , double *vec3 )
{
   double vecp[3] ;
   vec_crop( vec2 , vec3 , vecp ) ;
   return( vec_dotp(vec1 , vecp ) ) ;
}
void norm_3p( double *p1 , double *p2 ,double *p3 ,double *normal )
{
 double v12[3] , v13[3] ;
 vec_2p( p1 , p2 , v12 ) ;
 vec_2p( p1 , p3 , v13 ) ; // the direction of crop has changed ?!
 vec_crop( v12 , v13 , normal ) ;
}

double VolumOf4p(Point p0,Point p1,Point p2,Point p3)
{
  Point vec ,p01,p02,p03 ;
  p01[0] = p1[0] -p0[0] ;   p01[1] = p1[1] -p0[1] ;   p01[2] = p1[2] -p0[2] ;
  p02[0] = p2[0] -p0[0] ;   p02[1] = p2[1] -p0[1] ;   p02[2] = p2[2] -p0[2] ;
  p03[0] = p3[0] -p0[0] ;   p03[1] = p3[1] -p0[1] ;   p03[2] = p3[2] -p0[2] ;
  vec_crop(p01,p02,vec) ;                            /* valu  or parem ? */
  return vec_dotp(vec,p03) ;
}

int triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3]);

bool isTriangleBoxOver(Point p1 ,Point p2 ,Point p3 ,const double bd[6],double eps ){

//  int i ;
 
	double bound[6];
	double a[3];
	for(int i=0; i<3; i++)
		a[i]=bd[i+3]-bd[i];
	for(int i=0; i<3; i++){
		bound[i]=bd[i]-eps*a[i];
		bound[i+3]=bd[i+3]+eps*a[i];
	}
  double boxcenter[3],boxhalfsize[3],triverts[3][3];
  for(int i=0; i<3; i++){
	  boxcenter[i]=(bound[i]+bound[i+3])/2.;
	  boxhalfsize[i]=(bound[i+3]-bound[i])/2;
	  triverts[0][i]=p1[i];
 	  triverts[1][i]=p2[i];
	  triverts[2][i]=p3[i];
 }
  if(triBoxOverlap(boxcenter,boxhalfsize,triverts)==0) return false;
  else return true;
//  for(i=0 ; i<3 ; i++ )
  //  if( (p1[i]<box[i]&&p2[i]<box[i]&&p3[i]<box[i])||
	//   (p1[i]>box[i+3]&&p2[i]>box[i+3]&&p3[i]>box[i+3]) )return false ;
 // return true ;
}
bool isTriRayIntersect(double p0[3],double p1[3],double p2[3],double orig[3],double dir[3]){//, double &t){
	double tvec[3], pvec[3], qvec[3]; //revised from gel
	double det,inv_det;
	double edge[3][3];
//	return true;
	vec_2p(p0,p1,edge[0]);
	vec_2p(p0,p2,edge[2]);
   /* begin calculating determinant - also used to calculate U parameter */
    vec_crop(dir,edge[2],pvec);
   /* if determinant is near zero, ray lies in plane of triangle */
   det = vec_dotp(edge[0], pvec);

   if (det > -0.00000000001 && det < 0.00000000001)
     return false;
   inv_det = 1.0 / det;

   /* calculate distance from v0 to ray origin */
   vec_2p(p0,orig,tvec);

   /* calculate U parameter and test bounds */
   double u = vec_dotp(tvec, pvec) * inv_det;
   if (u < 0.0 || u > 1.0)
     return false;

   /* prepare to test V parameter */
    vec_crop(tvec, edge[0],qvec);

   /* calculate V parameter and test bounds */
   double v = vec_dotp(dir, qvec) * inv_det;
   if (v < 0.0 || u + v > 1.0)
     return false;

   /* calculate t, ray intersects triangle */
   return  vec_dotp(edge[2], qvec) * inv_det>=0;

//   return true;
}


/********************************************************/

/* AABB-triangle overlap test code                      */

/* by Tomas Akenine-Moller                              */

/* Function: int triBoxOverlap(double boxcenter[3],      */

/*          double boxhalfsize[3],double triverts[3][3]); */

/* History:                                             */

/*   2001-03-05: released the code in its first version */

/*   2001-06-18: changed the order of the tests, faster */

/*                                                      */

/* Acknowledgement: Many thanks to Pierre Terdiman for  */

/* suggestions and discussions on how to optimize code. */

/* Thanks to David Hunt for finding a ">="-bug!         */

/********************************************************/

#include <math.h>

#include <stdio.h>



#define X 0

#define Y 1

#define Z 2



#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0]; 



#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])



#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2]; 



#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;



int planeBoxOverlap(double normal[3], double vert[3], double maxbox[3])	// -NJMP-

{

  int q;

  double vmin[3],vmax[3],v;

  for(q=X;q<=Z;q++)

  {

    v=vert[q];					// -NJMP-

    if(normal[q]>0.0f)

    {

      vmin[q]=-maxbox[q] - v;	// -NJMP-

      vmax[q]= maxbox[q] - v;	// -NJMP-

    }

    else

    {

      vmin[q]= maxbox[q] - v;	// -NJMP-

      vmax[q]=-maxbox[q] - v;	// -NJMP-

    }

  }

  if(DOT(normal,vmin)>0.0f) return 0;	// -NJMP-

  if(DOT(normal,vmax)>=0.0f) return 1;	// -NJMP-

  

  return 0;

}





/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			       	   \
	p2 = a*v2[Y] - b*v2[Z];			       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



#define AXISTEST_X2(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			           \
	p1 = a*v1[Y] - b*v1[Z];			       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p2 = -a*v2[X] + b*v2[Z];	       	       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



#define AXISTEST_Y1(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p1 = -a*v1[X] + b*v1[Z];	     	       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



/*======================== Z-tests ========================*/



#define AXISTEST_Z12(a, b, fa, fb)			   \
	p1 = a*v1[X] - b*v1[Y];			           \
	p2 = a*v2[X] - b*v2[Y];			       	   \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;



#define AXISTEST_Z0(a, b, fa, fb)			   \
	p0 = a*v0[X] - b*v0[Y];				   \
	p1 = a*v1[X] - b*v1[Y];			           \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;



int triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3])

{



  /*    use separating axis theorem to test overlap between triangle and box */

  /*    need to test for overlap in these directions: */

  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */

  /*       we do not even need to test these) */

  /*    2) normal of the triangle */

  /*    3) crossproduct(edge from tri, {x,y,z}-directin) */

  /*       this gives 3x3=9 more tests */

   double v0[3],v1[3],v2[3];

//   double axis[3];

   double min,max,p0,p1,p2,rad,fex,fey,fez;		// -NJMP- "d" local variable removed

   double normal[3],e0[3],e1[3],e2[3];



   /* This is the fastest branch on Sun */

   /* move everything so that the boxcenter is in (0,0,0) */

   SUB(v0,triverts[0],boxcenter);

   SUB(v1,triverts[1],boxcenter);

   SUB(v2,triverts[2],boxcenter);



   /* compute triangle edges */

   SUB(e0,v1,v0);      /* tri edge 0 */

   SUB(e1,v2,v1);      /* tri edge 1 */

   SUB(e2,v0,v2);      /* tri edge 2 */



   /* Bullet 3:  */

   /*  test the 9 tests first (this was faster) */

   fex = fabs(e0[X]);

   fey = fabs(e0[Y]);

   fez = fabs(e0[Z]);

   AXISTEST_X01(e0[Z], e0[Y], fez, fey);

   AXISTEST_Y02(e0[Z], e0[X], fez, fex);

   AXISTEST_Z12(e0[Y], e0[X], fey, fex);



   fex = fabs(e1[X]);

   fey = fabs(e1[Y]);

   fez = fabs(e1[Z]);

   AXISTEST_X01(e1[Z], e1[Y], fez, fey);

   AXISTEST_Y02(e1[Z], e1[X], fez, fex);

   AXISTEST_Z0(e1[Y], e1[X], fey, fex);



   fex = fabs(e2[X]);

   fey = fabs(e2[Y]);

   fez = fabs(e2[Z]);

   AXISTEST_X2(e2[Z], e2[Y], fez, fey);

   AXISTEST_Y1(e2[Z], e2[X], fez, fex);

   AXISTEST_Z12(e2[Y], e2[X], fey, fex);



   /* Bullet 1: */

   /*  first test overlap in the {x,y,z}-directions */

   /*  find min, max of the triangle each direction, and test for overlap in */

   /*  that direction -- this is equivalent to testing a minimal AABB around */

   /*  the triangle against the AABB */



   /* test in X-direction */

   FINDMINMAX(v0[X],v1[X],v2[X],min,max);

   if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;



   /* test in Y-direction */

   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);

   if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;



   /* test in Z-direction */

   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);

   if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;



   /* Bullet 2: */

   /*  test if the box intersects the plane of the triangle */

   /*  compute plane equation of triangle: normal*x+d=0 */

   CROSS(normal,e0,e1);

   // -NJMP- (line removed here)

   if(!planeBoxOverlap(normal,v0,boxhalfsize)) return 0;	// -NJMP-



   return 1;   /* box and triangle overlaps */

}




