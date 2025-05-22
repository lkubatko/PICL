/* Functions to read and store gene trees */

#include <stdio.h>
#include <math.h>
#include <string.h>


int FindParentI(int target) {

  int i, j, parent=0;
  
  for (j=0; j<2; j++) {
    for (i=0; i<ntaxa; i++) {
      if (ppTwoRow_temp[j][i]==target) {
	parent = i+ntaxa+1;
	break;
      }
    }
  }

  /*if (parent==0) { 
     printf("Could not find parent for node %d\n",target); 
     }*/
  return(parent);

}


int NextNode(int next) {

  int stop=0;

  while (stop!=1) {
    if (filled_ind[next+1-(ntaxa+1)] != 1) {
      next+=1;
      stop=1;
    }
    else next+=1;
  }
  return next;

}

double ReadLength(FILE *fp) {

  float length=0.0;

  fscanf(fp,"%18f",&length);
  return (double)length;

}

int AddToNode(int current,int next,int last_tip,FILE *fp,double read_length) {

  int sss, junk, i, j, k, ss2, ss3, toadd;
  double temp;
  char tname[15];

  sss = fgetc(fp);

  if (ppTwoRow_temp[0][current-(ntaxa+1)] == 0) i = 0;
  else i=1;

  if (i==0) {

    ppTwoRow_temp[i][current-(ntaxa+1)] = last_tip;
    ppLengthMat[current][last_tip] = ppLengthMat[last_tip][current] = read_length;

    if (sss!='(') {

      ss2 = sss;
      tname[0] = ss2;
      tname[1] = 00;
      //printf("tname is %s\n",tname);
      //printf("It really is\n");
      //fflush(0);
      j=0;
      while (ss2 != ':') {
	ss2=fgetc(fp);
	j++;
	tname[j] = ss2;
      }
      tname[j] = 00;
      
      k=0;
      while (strcmp(tname,taxname[k])!=0 && k<=ntaxa) k++;
      if (k>=ntaxa) {
	printf("\n\t ERROR: Taxon name in gene tree (%s) does not match any name given in settings file.\n \t Exiting ...\n\n",tname);
	exit(1);
      }
      toadd = k+1;
      temp = ReadLength(fp);

      ppTwoRow_temp[1][current-(ntaxa+1)] = toadd;
      ppLengthMat[current][toadd] = ppLengthMat[toadd][current] = temp;
      sss = fgetc(fp);
      filled_ind[current-(ntaxa+1)] = 1;
      
    }

  }

  if (i==1 && sss!='(') {

    ss2=sss;
    tname[0] = ss2;
    j=0;
    while (ss2 != ':') {
      ss2=fgetc(fp);
      j++;
      tname[j] = ss2;
    }
    tname[j] = 00;

    k=0;
    while (strcmp(tname,taxname[k])!=0 && k<=ntaxa) k++;
    if (k>=ntaxa) {
      printf("\n\t ERROR: Taxon name in gene tree (%s) does not match any name given in settings file.\n \t Exiting ...\n\n",tname);
      exit(1);
    }
    toadd = k+1;
    temp = ReadLength(fp);
    
    ppTwoRow_temp[1][current-(ntaxa+1)] = toadd;
    ppLengthMat[current][toadd] = ppLengthMat[toadd][current] = temp;
    filled_ind[current-(ntaxa+1)] = 1;
    sss = fgetc(fp);

  }

  return sss;

}

int FindLastOpen(int current) {

  int node, stop=0;

  node = current-1;
  while(node>=ntaxa+1 && stop==0) {
    if (filled_ind[node-(ntaxa+1)]==1) node = node-1;
    else stop=1;
  }
  return node;

}

int CloseBack(int open, int current, int last_char, int last_tipp, FILE *fp, double read_length) {

  int sss;
  int j;

  sss=getc(fp);

  if (sss!=';') {
  
    if (ppTwoRow_temp[0][open-(ntaxa+1)] == 0) j=0;
    else {
      
      j=1;
      filled_ind[open-(ntaxa+1)] = 1;
      
    }
    
    if (sss == ':') {

      read_length = ReadLength(fp);
      sss = fgetc(fp);

    }
      
    ppTwoRow_temp[j][open-(ntaxa+1)] = current;  
    ppLengthMat[open][current] = ppLengthMat[current][open] = read_length;
   
  }
    
  return sss;
 
}

double ReadTree(FILE *fp) {

  int ss=0, ss2, ss3, last_ss, lastlast_ss, last_tip=0, j, k;
  int next_avail = ntaxa+1, last_open=ntaxa+1, curr_node=ntaxa+1;
  int parent;
  double read_length=0.0;
  float rate = 2.0;
  char tname[15];

  for (j=0; j<ntaxa+1; j++) {
    for (k=0; k<ntaxa+1; k++) {
      ppLengthMat[j][k] = ppLengthMat[k][j] = 1000;
    }
  }

//  while (ss!='[') ss=fgetc(fp);
//  fscanf(fp,"%14f",&rate);
//  ss=fgetc(fp);

  while (ss!='(') ss=fgetc(fp);
  last_ss = ss;
  ss = fgetc(fp);

  while (!feof(fp) && ss!=';') {

    if (ss == '(') {

      next_avail = NextNode(next_avail);
      curr_node = next_avail;
      last_ss = ss;
      ss = fgetc(fp);

    }

    else {

      if (ss == ',') {

	last_ss = ss;
	ss = AddToNode(curr_node,next_avail,last_tip,fp,read_length);
    
      }

      else {

	if (ss == ')') {

	  last_open = FindLastOpen(curr_node);
	  lastlast_ss = last_ss;
	  last_ss = ss;
	  ss = CloseBack(last_open,curr_node,lastlast_ss,last_tip,fp,read_length);
	  if (ss!=';' && filled_ind[last_open-(ntaxa+1)]==1) curr_node = FindParentI(curr_node);
	  else curr_node = last_open;
	 	  
	}
	
	else {

	  ss2=ss;
	  tname[0] = ss;
	  j=0;
	  while (ss2 != ':') {
	    ss2=fgetc(fp);
	    j++;
	    tname[j] = ss2;
	  }
	  tname[j] = 00;

	  k=0;
	  while (strcmp(tname,taxname[k])!=0 && k<=ntaxa) k++;
	  if (k>=ntaxa) {
	    printf("\n\t ERROR: Taxon name in gene tree (%s) does not match any name given in settings file.\n \t Exiting ...\n\n",tname);
	    exit(1);
	  }
	  last_tip = k+1;
    
	  last_ss = last_tip;
	  ss = ss2;
	  if (ss == ':') read_length = ReadLength(fp);
	  parent = FindParentI(last_tip);
	  if (parent != 0) ppLengthMat[parent][last_tip] = ppLengthMat[last_tip][parent] = read_length;
	  ss = fgetc(fp);
	  
	}

      }

    }
      
  }

  return rate;

}

double FindTotalTime() {

  int start=1, last_parent, parent;
  float sum_time = 0.0;

  /*  tree check */
  /*printf("tree check in find total time\n");
  for (i=0; i<2; i++) {
    for (j=0; j<ntaxa; j++) printf("%d ",ppTwoRow_temp[i][j]);
    printf("\n");
    }*/

  parent = start;
   while (parent != ntaxa+1) {
    last_parent = parent;
    parent = FindParentI(parent);
    sum_time += ppLengthMat[last_parent][parent];
  }

   //printf("total time is %f\n",sum_time);

   return (double)sum_time;

}

void CalcTimeVec(double total, double rate) {

  int last_parent, parent, i;

  for (i=1; i<ntaxa+1; i++) {
    parent = i;
    TimeVec_temp[parent] = total;
    while (parent != ntaxa+1) {
      last_parent = parent;
      parent = FindParentI(parent);
      if (ppLengthMat[parent][last_parent] < 0.00000001) ppLengthMat[parent][last_parent] = ppLengthMat[last_parent][parent] = 0.00000001;
      TimeVec_temp[parent] = TimeVec_temp[last_parent] - ppLengthMat[parent][last_parent];
      if (TimeVec_temp[parent] >= TimeVec_temp[last_parent]) TimeVec_temp[parent] = TimeVec_temp[last_parent] - 0.00000001;
    }
  }

  for (i=1; i<2*ntaxa+1; i++) TimeVec_temp[i] = (total - TimeVec_temp[i]);
  for (i=ntaxa+1; i<2*ntaxa+1; i++) {
    if (TimeVec_temp[FindParentI(i)] == TimeVec_temp[i]) {
      if (TimeVec_temp[i] - TimeVec_temp[ppTwoRow_temp[0][i-(ntaxa+1)]] > 0.00000001 && TimeVec_temp[i]-TimeVec_temp[ppTwoRow_temp[1][i-(ntaxa+1)]] > 0.00000001) TimeVec_temp[i] = TimeVec_temp[i] - 0.00000001;
      else {
	if (TimeVec_temp[ppTwoRow_temp[0][i-(ntaxa+1)]]>TimeVec_temp[ppTwoRow_temp[1][i-(ntaxa+1)]]) TimeVec_temp[i] = TimeVec_temp[i] - 0.5*(TimeVec_temp[i] - TimeVec_temp[ppTwoRow_temp[0][i-(ntaxa+1)]]);
	else TimeVec_temp[i] = TimeVec_temp[i] - 0.5*(TimeVec_temp[i] - TimeVec_temp[ppTwoRow_temp[1][i-(ntaxa+1)]]);
      }
    }
  }
  //for (i=1; i<2*ntaxa+1; i++) printf("%d %1.20f\n",i,TimeVec_temp[i]); 
}

