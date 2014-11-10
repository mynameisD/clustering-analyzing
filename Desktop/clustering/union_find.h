void searchedge(SHEET *,SHEET *,POLY *);
void searchnode(POLY *,SHEET *,POLY *);
int findunion(SHEET *, POLY *);
void group(SHEET *,POLY *);
int spanbeta(SHEET *,POLY *);

SHEET * chainroot(POLY *);
SHEET * noderoot(SHEET *);

void groupanalyze(int,SHEET *,POLY *);
int midsplitunion(int,int,SHEET *,int,POLY *);
void midsplitedge(int,SHEET *,int,SHEET *,int,POLY *);
void midsplitnode(int,POLY *,int,SHEET *,int,POLY *);
int pbcsplitunion(int,int,SHEET *,int,POLY *);
void pbcsplitedge(int,SHEET *,int,SHEET *,int,POLY *);
void pbcsplitnode(int,POLY *,int,SHEET *,int,POLY *);
int partitiontest(int,int,SHEET *,int,POLY *);
void clearunion(int,SHEET *,int,POLY *);

int findunion(SHEET *s,POLY *p)
{
  int i,j,k;
  int countsheet=0,countchain=0;
  for(i=0;i<N_CLAY;i++){
    if(s[i].parent==NULL){
      if(s[i].group_id!=0) printf("ERROR,SHEET ROOT WITH GROUP ID NOT 0\n");
      countsheet++;
      s[i].group_id = countsheet;//mark itself as parent,it become root.
      searchedge(&s[i],s,p);//depth first search, root node always
    }
  }
  //printf("number of root sheet found: %d\n",countsheet);
  for(i=0;i<N_CHAIN;i++){
    if(p[i].parent==NULL){
      if(p[i].group_id!=0) printf("ERROR,CHAIN ROOT WITH GROUP ID NOT 0\n");
      countchain++;
      p[i].group_id = countchain+countsheet;
      searchnode(&p[i],s,p);
    }
  }
  //printf("number of root chain found: %d\n",countchain);
  return countchain+countsheet;
}
void searchedge(SHEET *this,SHEET *s,POLY *p)
{
  int i,j;
  int count=0;
  for(i=0;i<N_CHAIN;i++){
    if(p[i].parent==NULL){//only test if the the chain has no parent and not root
      if(p[i].group_id!=0) printf("ERROR,SEARCHING EDGE WITH GROUP NOT 0\n");
      if(dis(this,&p[i])<=CUT_OFF){
	p[i].parent=this;
	p[i].group_id=this->group_id;
	searchnode(&p[i],s,p);
      }
    }
  }
  //printf("EDGE FOUND %d \n",count);
}
void searchnode(POLY *this,SHEET *s,POLY *p)
{
  int i,j;
  int count=0;
  for(i=0;i<N_CLAY;i++){
    if(s[i].parent==NULL&&s[i].group_id==0){//only if not root and has no parent
      if(dis(&s[i],this)<=CUT_OFF){
	s[i].parent=this;
	s[i].group_id=this->group_id;
	searchedge(&s[i],s,p);
      }
    }
  }
}

SHEET * chainroot(POLY *this)
{
  if(this->parent == NULL) return NULL;
  else return noderoot(this->parent);
}
SHEET * noderoot(SHEET *this)
{
  if(this->parent==NULL) return this;
  else return chainroot(this->parent);
}

void groupanalyze(int step,SHEET *s,POLY *p)
{
  int i,j,k,l;
  int countsheet,countchain,group_id=1;
  do{
    countsheet = 0;
    countchain = 0;
    for(i=0;i<N_CLAY;i++){
      if(s[i].group_id==group_id) countsheet++;
    }
    for(i=0;i<N_CHAIN;i++){
      if(p[i].group_id==group_id) countchain++;
    }
    if(countsheet>0||countchain>0){
      //printf("Group: %d, %d sheets, %d chains, %d in all\n",group_id,countsheet,countchain,countsheet*SIZE_SHEET+countchain*SIZE_CHAIN);
      group_id++;
    }
  } while(countsheet>0||countchain>0);

  //new array of each group
  group_id--;//this is the number of groups
  int **groupsize=(int **)malloc(sizeof(int*)*group_id);
  for(i=0;i<group_id;i++){
    groupsize[i]=(int *)malloc(sizeof(int)*2);
    groupsize[i][0]=0;groupsize[i][1]=0;
  }
  int *spanlog=(int *)malloc(sizeof(int)*group_id); for(i=0;i<group_id;i++) spanlog[i]=0;
  int spandim[3];
  for(i=0;i<N_CLAY;i++) groupsize[s[i].group_id-1][0]++;
  for(i=0;i<N_CHAIN;i++) groupsize[p[i].group_id-1][1]++;


  //analyze group sequentially
  int left,right;
  for(i=0;i<group_id;i++){
    countsheet=0;
    countchain=0;
    if(groupsize[i][0]>1&&groupsize[i][1]>1){
      SHEET *temps=(SHEET *)malloc(groupsize[i][0]*sizeof(SHEET));
      POLY *tempc=(POLY *)malloc(groupsize[i][1]*sizeof(POLY));
      for(j=0;j<N_CLAY;j++){
	if(s[j].group_id==i+1){
	  sheetcopy(&temps[countsheet],&s[j]);
	  countsheet++;
	}
      }
      for(j=0;j<N_CHAIN;j++){
	if(p[j].group_id==i+1){
	  chaincopy(&tempc[countchain],&p[j]);
	  countchain++;
	}
      }
      vmdgroup(step,i+1,groupsize[i][0],temps,groupsize[i][1],tempc);
      for(k=0;k<3;k++){
	spandim[k]=0;
	clearunion(groupsize[i][0],temps,groupsize[i][1],tempc);
	if(partitiontest(k,groupsize[i][0],temps,groupsize[i][1],tempc)==1&&pbcsplitunion(k,groupsize[i][0],temps,groupsize[i][1],tempc)==1){
	  for(j=0;j<countsheet;j++){
	    temps[j].parent=NULL;
	    temps[j].group_id=0;
	  }
	  for(j=0;j<countchain;j++){
	    tempc[j].parent=NULL;
	    tempc[j].group_id=0;
	  }
	  if(midsplitunion(k,groupsize[i][0],temps,groupsize[i][1],tempc)==1){
	    spandim[k]=1;
	  }
	}
	else spandim[k]=0;	
      }
      if(spandim[0]+spandim[1]+spandim[2]>2){
	spanlog[i]=1;
	printf("SPANNING GROUP %d, size %d\n",i+1,groupsize[i][0]*SIZE_SHEET+groupsize[i][1]*SIZE_CHAIN);
      }
      else spanlog[i]=0;
      free(temps);
      free(tempc);
    }
    else spanlog[i]=0;
  }
  free(groupsize);
}

int midsplitunion(int dim,int ns,SHEET *s,int nc,POLY *p)
{
  int i,j,k;
  int countsheet=0,countchain=0;
  for(i=0;i<ns;i++){
    if(s[i].parent==NULL){
      if(s[i].group_id!=0) printf("ERROR,SHEET ROOT WITH GROUP ID NOT 0\n");
      countsheet++;
      s[i].group_id = countsheet;//mark itself as parent,it become root.
      midsplitedge(dim,&s[i],ns,s,nc,p);//depth first search, root node always
    }
  }
  printf("midsplit %d lead to %d group(origin size %d)\n",dim,countsheet,ns*SIZE_SHEET+nc*SIZE_CHAIN);
  //getchar();
  //printf("number of root chain found: %d\n",countchain);
  return countsheet;
}
void midsplitedge(int dim,SHEET *this,int ns,SHEET *s,int nc,POLY *p)
{
  int i,j;
  int count=0;
  int sidechain=1,sidesheet=0;
  if(this->com[dim]<0.5*PBC&&this->com[dim]>0.25*PBC) sidesheet=-1;//this sheet is in mid left range
  else if(this->com[dim]>0.5*PBC&&this->com[dim]<0.75*PBC) sidesheet=1;//this sheet is in mid right range 

  for(i=0;i<nc;i++){
    if(p[i].parent==NULL){//only test if the the chain has no parent and not root
      if(p[i].group_id!=0) printf("ERROR,SEARCHING EDGE WITH GROUP NOT 0\n");
      sidechain=1;
      for(j=0;j<SIZE_CHAIN;j++){
	if(p[i].elm[j][dim]<0.5*PBC){sidechain=-1;break;}//this chain is in left box
      }
      if(sidechain*sidesheet!=-1&&dis(this,&p[i])<=CUT_OFF){
	p[i].parent=this;
	p[i].group_id=this->group_id;
	midsplitnode(dim,&p[i],ns,s,nc,p);
      }
    }
  }
  //printf("EDGE FOUND %d \n",count);
}
void midsplitnode(int dim,POLY *this,int ns,SHEET *s,int nc,POLY *p)
{
  int i,j;
  int count=0;
  int sidechain=1,sidesheet=0;
  for(i=0;i<SIZE_CHAIN;i++){
    if(this->elm[i][dim]<0.5*PBC){sidesheet=-1;break;}//this chain is in left half
  }
  for(i=0;i<ns;i++){
    if(s[i].parent==NULL&&s[i].group_id==0){//only if not root and has no parent
      sidesheet=0;
      if(s[i].com[dim]<0.5*PBC&&s[i].com[dim]>0.25*PBC) sidechain=-1;//this sheet is in mid right range
      else if(s[i].com[dim]>0.5*PBC&&s[i].com[dim]<0.75*PBC) sidechain=1;

      if(sidechain*sidesheet!=-1&&dis(&s[i],this)<=CUT_OFF){
	s[i].parent=this;
	s[i].group_id=this->group_id;
	midsplitedge(dim,&s[i],ns,s,nc,p);
      }
    }
  }
}
int pbcsplitunion(int dim,int ns,SHEET *s,int nc,POLY *p)
{
  int i,j,k;
  int countsheet=0,countchain=0;
  for(i=0;i<ns;i++){
    if(s[i].parent==NULL){
      if(s[i].group_id!=0) printf("ERROR,SHEET ROOT WITH GROUP ID NOT 0\n");
      countsheet++;
      s[i].group_id = countsheet;//mark itself as parent,it become root.
      pbcsplitedge(dim,&s[i],ns,s,nc,p);//depth first search, root node always
    }
  }
  printf("pbcsplit %d lead to %d group(origin size %d)\n",dim,countsheet,ns*SIZE_SHEET+nc*SIZE_CHAIN);
  //getchar();
  //printf("number of root chain found: %d\n",countchain);
  return countsheet;
}
void pbcsplitedge(int dim,SHEET *this,int ns,SHEET *s,int nc,POLY *p)
{
  int i,j;
  int count=0;
  int sidechain=1,sidesheet=0;
  if(this->com[dim]<0.25*PBC) sidesheet=-1;//this sheet is in left edge
  else if(this->com[dim]>0.75*PBC) sidesheet=1;//this sheet is in right edge 

  for(i=0;i<nc;i++){
    if(p[i].parent==NULL){//only test if the the chain has no parent and not root
      if(p[i].group_id!=0) printf("ERROR,SEARCHING EDGE WITH GROUP NOT 0\n");
      sidechain=1;
      for(j=0;j<SIZE_CHAIN;j++){
	if(p[i].elm[j][dim]<0.5*PBC){sidechain=-1;break;}//this chain is in left edge
      }
      if(sidechain*sidesheet!=-1&&dis(this,&p[i])<=CUT_OFF){
	p[i].parent=this;
	p[i].group_id=this->group_id;
	pbcsplitnode(dim,&p[i],ns,s,nc,p);
      }
    }
  }
  //printf("EDGE FOUND %d \n",count);
}

void pbcsplitnode(int dim,POLY *this,int ns,SHEET *s,int nc,POLY *p)
{
  int i,j;
  int count=0;
  int sidechain=1,sidesheet=0;
  for(i=0;i<SIZE_CHAIN;i++){
    if(this->elm[i][dim]<0.5*PBC){sidesheet=-1;break;}//this chain is in left half
  }
  for(i=0;i<ns;i++){
    if(s[i].parent==NULL&&s[i].group_id==0){//only if not root and has no parent
      sidesheet=0;
      if(s[i].com[dim]<0.25*PBC) sidechain=-1;//this sheet is in mid right range
      else if(s[i].com[dim]>0.75*PBC) sidechain=1;

      if(sidechain*sidesheet!=-1&&dis(&s[i],this)<=CUT_OFF){
	s[i].parent=this;
	s[i].group_id=this->group_id;
	pbcsplitedge(dim,&s[i],ns,s,nc,p);
      }
    }
  }
}
int partitiontest(int dim,int ns,SHEET *s,int nc,POLY *p)
{
  int i,j,k;
  int *bin=(int *)malloc(sizeof(int)*PBC/2);
  for(i=0;i<PBC/2;i++){
    bin[i]=0;
  }
  for(i=0;i<ns;i++){
    for(j=0;j<SIZE_SHEET;j++){
      bin[(int)(s[i].elm[j][dim]/2)]++;
    }
  }
  for(i=0;i<nc;i++){
    for(j=0;j<SIZE_CHAIN;j++){
      bin[(int)(p[i].elm[j][dim]/2)]++;
    }
  }
  for(i=0;i<PBC/2;i++){
    if(bin[i]==0){free(bin); return 0;}
  }
  free(bin);
  return 1;
}
void clearunion(int ns,SHEET *s,int nc,POLY *p)
{
  int i;
  for(i=0;i<ns;i++){
    s[i].parent=NULL;
    s[i].group_id=0;
  }
  for(i=0;i<nc;i++){
    p[i].parent=NULL;
    p[i].group_id=0;
  }
}

