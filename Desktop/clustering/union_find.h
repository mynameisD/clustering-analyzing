void searchedge(SHEET *,SHEET *,POLY *);
void searchnode(POLY *,SHEET *,POLY *);
void findunion(SHEET *, POLY *);
void groupanalyze(SHEET *,POLY *);

SHEET * chainroot(POLY *);
SHEET * noderoot(SHEET *);

void findunion(SHEET *s,POLY *p)
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

