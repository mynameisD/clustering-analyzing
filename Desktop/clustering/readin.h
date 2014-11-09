#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define N_SHEETS 3800
#define ENTRIES_BEFORESHEAR 1000
#define PBC 71.538792
#define N_POLYMER 80000
#define TYPESHEET 3
#define TYPESTICKER 2
#define TYPECHAIN 1
#define CUT_OFF 1.12
#define SIZE_CHAIN 20
#define SIZE_SHEET 19
#define N_CLAY (N_SHEETS/SIZE_SHEET)
#define N_CHAIN (N_POLYMER/SIZE_CHAIN)

struct sheet{
  double elm[SIZE_SHEET][3];
  int type[SIZE_SHEET];
  int mol_id;
  int count;
  int group_id;
  double com[3];
  double norm[3];
  struct poly *parent;//the id of parent sheet node
};

struct poly{
  double elm[SIZE_CHAIN][3];
  int type[SIZE_CHAIN];
  int mol_id;
  int count;
  int group_id;
  struct sheet *parent;//the id of parent sheet node
};

typedef struct poly POLY;
typedef struct sheet SHEET;

void clearinfo(struct sheet *,struct poly *);
void readfiller(FILE *,struct sheet *);
void readchain(FILE *,struct poly *);
double dis(SHEET *,POLY *);
void systeminfo();

void systeminfo()
{
  printf("THE PERIODIC BOUNDARY: %lf\n",PBC);
  printf("THE CUTOFF VALUE: %lf\n",CUT_OFF);
  printf("THE NUMBER OF POLYMERS: %d\n",N_POLYMER);
  printf("THE NUMBER OF CHAINS: %d\n",N_CHAIN);
  printf("THE NUMBER OF SHEETS: %d\n",N_SHEETS);
  printf("THE NUMBER OF SHEET PALATLET: %d\n",N_CLAY);
  printf("PRESS TO CONTINUE\n");
  getchar();
}
void clearinfo(struct sheet *s,struct poly *c)
{
  int i,j;
  for(i=0;i<N_CLAY;i++){
    s[i].group_id=0;
    s[i].mol_id=0;
    s[i].parent=NULL;
    s[i].count=0;
    for(j=0;j<SIZE_SHEET;j++){
      s[i].elm[j][0]=-10000.0;
      s[i].elm[j][1]=-10000.0;
      s[i].elm[j][2]=-10000.0;
      s[i].type[j]=0;
    }
    for(j=0;j<3;j++){
      s[i].com[j]=-10000.0;
      s[i].norm[j]=-10000.0;
    }
  }
  for(i=0;i<N_CHAIN;i++){
    c[i].group_id=0;
    c[i].mol_id=0;
    c[i].count=0;
    c[i].parent=NULL;
    for(j=0;j<SIZE_SHEET;j++){
      c[i].elm[j][0]=-10000.0;
      c[i].elm[j][1]=-10000.0;
      c[i].elm[j][2]=-10000.0;
      c[i].type[j]=0;
    }
  }
}
void readfiller(FILE *f,struct sheet *s)
{
  /*read in cordinates in wraped form*/
  int i,j;
  int index,molid,type;
  int seq;
  double cordx,cordy,cordz;
  size_t size_buffer = 200;
  char *skip = (char *)malloc(size_buffer);
  for(i=0;i<9;i++){
    getline(&skip,&size_buffer,f);
    if(i==3&&atoi(skip)!=N_SHEETS){printf("ERROR: The number of sheets not match record");exit(55);}
    if(i==1) printf("TIMESTEP:%s",skip);
  }
  for(i=0;i<N_SHEETS;i++){
    fscanf(f,"%d %d %d %lf %lf %lf\n",&index,&molid,&type,&cordx,&cordy,&cordz);
    molid=molid-N_CHAIN;//the order of sheets
    if(molid>N_CLAY){printf("ERROR:index of clay exceed limite\n");exit(67);}
    seq=s[molid-1].count;
    s[molid-1].mol_id=molid+N_CHAIN;
    s[molid-1].type[seq]=type;
    s[molid-1].elm[seq][0]=cordx;
    s[molid-1].elm[seq][1]=cordy;
    s[molid-1].elm[seq][2]=cordz;
    s[molid-1].count++;
  }
  for(i=0;i<N_CLAY;i++){if(s[i].count!=SIZE_SHEET){printf("ERROR, THE NUMBER IN CLAY %d NOT MATCH SIZE OF SHEETS",i+1+N_CHAIN);exit(66);}}
  free(skip);
}
void readchain(FILE *f,struct poly *c)
{
  int i,j;
  int index,molid,type;
  int seq;
  double cordx,cordy,cordz;
  size_t size_buffer = 200;
  char *skip=(char *)malloc(size_buffer);
  for(i=0;i<9;i++){
    getline(&skip,&size_buffer,f);
    if(i==3&&atoi(skip)!=N_POLYMER){printf("ERROR:THE NUMBER OF POLYMERS NOT MATCH FILE");exit(87);}
  }
  for(i=0;i<N_POLYMER;i++){
    fscanf(f,"%d%d%d%lf%lf%lf\n",&index,&molid,&type,&cordx,&cordy,&cordz);
    molid=molid;if(molid>N_CHAIN){printf("ERROR:MOLID OF POLYMER EXCEED LIMIT\n");exit(93);}
    seq=c[molid-1].count;
    c[molid-1].type[seq]=type;
    c[molid-1].mol_id=molid;
    c[molid-1].elm[seq][0]=cordx;
    c[molid-1].elm[seq][1]=cordy;
    c[molid-1].elm[seq][2]=cordz;
    c[molid-1].count++;
  }
  for(i=0;i<N_CHAIN;i++){if(c[i].count!=SIZE_CHAIN){printf("ERROR,THE SIZE OF CHAIN NOT MATCH\n");exit(102);}}
  free(skip);
}

double dis(SHEET *s,POLY *p)
{
  double dist;
  double tempx,tempy,tempz;
  double distance=PBC;
  int i,j,k;
  for(j=0;j<SIZE_CHAIN;j++){
    if(p->type[j]==TYPESTICKER){//only test sticker
      for(i=0;i<SIZE_SHEET;i++){
	tempx=fabs(s->elm[i][0] - p->elm[j][0]);
	tempy=fabs(s->elm[i][1] - p->elm[j][1]);
	tempz=fabs(s->elm[i][2] - p->elm[j][2]);
	tempx=tempx>(PBC/2.0)? (PBC-tempx):tempx;
	tempy=tempy>(PBC/2.0)? (PBC-tempy):tempy;
	tempz=tempz>(PBC/2.0)? (PBC-tempz):tempz;
	dist=tempx*tempx+tempy*tempy+tempz*tempz;
	dist=sqrt(dist);
	distance=distance>dist?dist:distance;
      }
    }
  }
  return distance;
}
