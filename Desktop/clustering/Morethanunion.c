/*this code use union find method to */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "readin.h"
#include "union_find.h"
#include "bridges.h"
void main()
{
  int i,j;
  clock_t time;
  SHEET *clay=(SHEET*)malloc(sizeof(SHEET)*N_CLAY);
  POLY *chain=(POLY*)malloc(sizeof(POLY)*N_CHAIN);
  FILE *f_filler=fopen("..//filler.lammpstrj","r");
  FILE *f_chain=fopen("..//chains.lammpstrj","r");  
  FILE *gel=fopen("percent_gel.txt","w");
  fclose(gel);
  int n_entries=0;
  clearbridges();
  systeminfo();
  while(!feof(f_filler)&&!feof(f_chain)&&n_entries<ENTRIES_BEFORESHEAR){
    time=clock();
    if(feof(f_filler)||feof(f_chain)){printf("ERROR: Two file not of same length, or reading speed not sync\n");exit(17);}
    clearinfo(clay,chain);
    readfiller(f_filler,clay);
    readchain(f_chain,chain);
    findcom(clay);
    bridges(clay,chain);
    findunion(clay,chain);
    groupanalyze(n_entries,clay,chain);
    n_entries++;
    printf("TIME COST: %f \n",((float)(clock()-time))/CLOCKS_PER_SEC);
    //getchar();
  }
}
