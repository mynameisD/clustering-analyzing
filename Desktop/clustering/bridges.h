#include <stdlib.h>
#include <stdio.h>

void bridges(SHEET *, POLY *);
void clearbridges();

void clearbridges()
{
  FILE *chaininfo=fopen("chain_info.txt","w");
  fprintf(chaininfo,"N_BRIDGES  N_DANGLE  N_FREE\n");
  fclose(chaininfo);
}
void bridges(SHEET *s,POLY *p)
{
  int i,j;
  int count;
  int nbridge=0,ndangling=0,nfree=0;
  for(i=0;i<N_CHAIN;i++){
    count=0;
    for(j=0;j<N_CLAY;j++){
      if(dis(&s[j],&p[i])<=CUT_OFF) count++;
      if(count==2) break;
    }
    if(count==2) nbridge++;
    else if(count==1) ndangling++;
    else if(count==0) nfree++;
    else printf("ERROR:SOMETHING WRONG IN BRIDGES\n");
  }
  if(ndangling+nfree+nbridge!=N_CHAIN) printf("ERROR:TOTAL NUMBER OF CHAIN NOT MATCHING\n");
  FILE *chaininfo=fopen("chain_info.txt","a");
  fprintf(chaininfo,"%d %d %d\n",nbridge,ndangling,nfree);
  fclose(chaininfo);
}
void vmdbackbone(int step,int gid,int gssize,SHEET *s,int gcsize,POLY *p)
{
  int i,j;
  int count;
  int nbridge=0,ndangling=0,nfree=0;
  int *bridge=(int *)malloc(sizeof(int)*gcsize);
  for(i=0;i<gcsize;i++){bridge[i]=0;}
  for(i=0;i<gcsize;i++){
    count=0;
    for(j=0;j<gssize;j++){
      if(dis(&s[j],&p[i])<=CUT_OFF) count++;
      if(count==2){
	bridge[i]=1;//mark as bridge
	break;
      }
    }
    if(count==2) nbridge++;
    else if(count==1) ndangling++;
    else if(count==0) printf("WARNING:FREE CHAIN FOUND IN UNION GROUP\n");
    else printf("ERROR:SOMETHING WRONG IN BRIDGES\n");
  }
  if(ndangling+nfree+nbridge!=gcsize) printf("ERROR:TOTAL NUMBER OF CHAIN NOT MATCHING\n");

  //start dump vmd
  int countbond=0;
  char str[35],str2[35];
  sprintf(str,"bb_s%dg%d_sheet.data",step,gid);
  FILE *fsheet=fopen(str,"w");
  sprintf(str2,"bb_s%dg%d_poly.data",step,gid);
  FILE *fpoly=fopen(str2,"w");

  fprintf(fsheet,"\n%d atoms\n",gssize*SIZE_SHEET);
  fprintf(fsheet,"%d bonds\n",0);
  fprintf(fsheet,"0 angles\n");
  fprintf(fsheet,"1 atom types\n");
  fprintf(fsheet,"0 bond types\n\n");
  fprintf(fsheet,"%lf %lf xlo xhi\n",0.00,PBC);
  fprintf(fsheet,"%lf %lf ylo yhi\n",0.00,PBC);
  fprintf(fsheet,"%lf %lf zlo zhi\n",0.00,PBC);
  fprintf(fsheet,"\n Atoms\n\n");
  for(i=0;i<gssize;i++){
    for(j=0;j<SIZE_SHEET;j++){
      fprintf(fsheet,"%d %d %d %lf %lf %lf\n",s[i].index[j],s[i].mol_id,s[i].type[j],s[i].elm[j][0],s[i].elm[j][1],s[i].elm[j][2]);
    }
  }

  int nbonds=0;//calculate number of bonds
  for(i=0;i<gcsize;i++){      
    if(bridge[i]==1){
      for(j=0;j<SIZE_CHAIN-1;j++){
	if(fabs(p[i].elm[j][0]-p[i].elm[j+1][0])<5.0&&fabs(p[i].elm[j][1]-p[i].elm[j+1][1])<5.0&&fabs(p[i].elm[j][2]-p[i].elm[j+1][2])<5.0)
	  nbonds++;
      }
    }
  }
  
  fprintf(fpoly,"\n%d atoms\n",nbridge*SIZE_CHAIN);
  fprintf(fpoly,"%d bonds\n",nbonds);
  fprintf(fpoly,"0 angles\n");
  fprintf(fpoly,"2 atom types\n");
  fprintf(fpoly,"1 bond types\n\n");
  fprintf(fpoly,"%lf %lf xlo xhi\n",0.00,PBC);
  fprintf(fpoly,"%lf %lf ylo yhi\n",0.00,PBC);
  fprintf(fpoly,"%lf %lf zlo zhi\n",0.00,PBC);
  fprintf(fpoly,"\n Atoms\n\n");
  for(i=0;i<gcsize;i++){
    if(bridge[i]==1){
      for(j=0;j<SIZE_CHAIN;j++){
	fprintf(fpoly,"%d %d %d %lf %lf %lf\n",p[i].index[j],p[i].mol_id,p[i].type[j],p[i].elm[j][0],p[i].elm[j][1],p[i].elm[j][2]);
      }
    }
  }
  fprintf(fpoly,"\n Bonds\n\n");
  nbonds=0;
  for(i=0;i<gcsize;i++){
    if(bridge[i]==1){
      for(j=0;j<SIZE_CHAIN-1;j++){
	if(fabs(p[i].elm[j][0]-p[i].elm[j+1][0])<5.0&&fabs(p[i].elm[j][1]-p[i].elm[j+1][1])<5.0&&fabs(p[i].elm[j][2]-p[i].elm[j+1][2])<5.0)
	  {
	    nbonds++;
	    fprintf(fpoly,"%d %d %d %d\n",nbonds,1,p[i].index[j],p[i].index[j+1]);
	  }
      }
    }
  }
  
  free(bridge);
  fclose(fsheet);
  fclose(fpoly);
}
