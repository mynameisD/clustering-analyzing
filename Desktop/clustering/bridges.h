#include <stdlib.h>
#include <stdio.h>

void bridges(SHEET *, POLY *);
void clearbridges();
void vmdbackbone(int,int,int,temps,groupsize[i][1],tempc);

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
  for(i=0;i<N_CHAIN;i++){
    count=0;
    for(j=0;j<N_CLAY;j++){
      if(dis(&s[j],&p[i])<=CUT_OFF) count++;
      if(count==2) break;
    }
    if(count==2) nbridge++;
    else if(count==1) ndangling++;
    else if(count==0) printf("WARNING:FREE CHAIN FOUND IN GROUP");
    else printf("ERROR:SOMETHING WRONG IN BRIDGES\n");
  }
  if(ndangling+nfree+nbridge!=N_CHAIN) printf("ERROR:TOTAL NUMBER OF CHAIN NOT MATCHING\n");
  FILE *chaininfo=fopen("chain_info.txt","a");
  fprintf(chaininfo,"%d %d %d\n",nbridge,ndangling,nfree);
  fclose(chaininfo);
}
