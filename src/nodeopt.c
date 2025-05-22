/**** Functions for optimization of nodes times ****/

/* Populates ppNodeChildren row 0 with left children of node mynode; [0][0] entry is number of children */
void AllMyChildrenLeft(int mynode) {

  int i, j=0, currentpar;
  int leftpar = ppTwoRow[0][mynode-(ntaxa+1)];

  for (i=1; i<ntaxa+1; i++)  {

	currentpar = i;
	while (currentpar != leftpar && currentpar != ntaxa+1) currentpar = find_parent(currentpar);
	if (currentpar!=ntaxa+1) {
		j++;
		ppNodeChildren[0][j] = i;
	}
  }
  ppNodeChildren[0][0] = j;

}

/* Populates ppNodeChildren row 1 with right children of node mynode; [1][0] entry is number of children */
void AllMyChildrenRight(int mynode) {

  int i, j=0, currentpar;
  int rightpar = ppTwoRow[1][mynode-(ntaxa+1)];

  for (i=1; i<ntaxa+1; i++)  {

        currentpar = i;
        while (currentpar != rightpar && currentpar != ntaxa+1) currentpar = find_parent(currentpar);
        if (currentpar != ntaxa+1) {
		j++;
		ppNodeChildren[1][j] = i;
	}
  }
  ppNodeChildren[1][0] = j;

}


/* Note that this function finds the children of the parent on mynode on the "other" side */
/* Populates row 2 of ppNodeChildren; [2][0] entry is number of children                  */

void AllMyChildren(int mynode) {

  int i, j=0, currentpar;
  int par, startnode;

  par = find_parent(mynode);
  if (ppTwoRow[0][par-(ntaxa+1)] == mynode) startnode = ppTwoRow[1][par-(ntaxa+1)];
  else startnode = ppTwoRow[0][par-(ntaxa+1)];
   
  for (i=1; i<ntaxa+1; i++) {

	currentpar = i;
	while(currentpar != startnode && currentpar != ntaxa+1) currentpar = find_parent(currentpar);
	if (currentpar != ntaxa+1) {
		j++;
		ppNodeChildren[2][j] = i;
	}

  }
  ppNodeChildren[2][0] = j;
}


double GetCompMOM(int ltax1, int ltax2, int ltax3, int ltax4, int issymmetric, int whichtau) {
  
  int count_array[16][16], p[15];
  int i, j, count_noambigs=0, sumps=0;
  double m = 4.0/3.0, t=0.01; // t is 2*theta
  double q[15], tau;

  for (i=0; i<16; i++) for(j=0; j<16; j++) count_array[i][j]=0;

    for (i=0; i<num_unique; i++) {
      if (ppBase_unique[ltax1][i]<4 && ppBase_unique[ltax2][i]<4 && ppBase_unique[ltax3][i]<4 && ppBase_unique[ltax4][i]<4) {        
            count_noambigs+=site_counter[i];
            count_array[ppBase_unique[ltax1][i]*4+ppBase_unique[ltax2][i]][ppBase_unique[ltax3][i]*4+ppBase_unique[ltax4][i]] += site_counter[i];
      }
    }

  //printf("\n\n count_noambigs is %d\n\n",count_noambigs);
    // order is:
    // 0 = xxxx
    // 1 = xxxy
    // 2 = xxyx
    // 3 = xyxx
    // 4 = yxxx
    // 5 = xyxy
    // 6 = yxxy
    // 7 = xxyy
    // 8 = xyxz
    // 9 = xyzx
    // 10 = yxxz
    // 11 = yxzx
    // 12 = xxyz
    // 13 = yzxx
    // 14 = xyzw

    p[0] = count_array[0][0] + count_array[5][5] + count_array[10][10] + count_array[15][15];
    p[1] = count_array[0][1] + count_array[0][2] + count_array[0][3] + count_array[5][4] + count_array[5][6] + count_array[5][7] + count_array[10][8] + count_array[10][9] + count_array[10][11] + count_array[15][12] + count_array[15][13] + count_array[15][14];
    p[2] = count_array[0][4] + count_array[0][8] + count_array[0][12] + count_array[5][1] + count_array[5][9] + count_array[5][13] + count_array[10][2] + count_array[10][6] + count_array[10][14] + count_array[15][3] + count_array[15][7] + count_array[15][11];
    p[7] = count_array[0][5] + count_array[0][10] + count_array[0][15] + count_array[5][0] + count_array[5][10] + count_array[5][15] + count_array[10][0] + count_array[10][5] + count_array[10][15] + count_array[15][0] + count_array[15][5] + count_array[15][10];
    p[12] = count_array[0][6] + count_array[0][7] + count_array[0][9] + count_array[0][11] + count_array[0][13] + count_array[0][14]
    + count_array[5][2] + count_array[5][3] + count_array[5][8] + count_array[5][11] + count_array[5][12] + count_array[5][14]
    + count_array[10][1] + count_array[10][3] + count_array[10][4] + count_array[10][7] + count_array[10][12] + count_array[10][13]
    + count_array[15][1] + count_array[15][2] + count_array[15][4] + count_array[15][6] + count_array[15][8] + count_array[15][9];
    p[3] = count_array[1][0] + count_array[2][0] + count_array[3][0] + count_array[4][5] + count_array[6][5] + count_array[7][5] + count_array[8][10] + count_array[9][10] + count_array[11][10] + count_array[12][15] + count_array[13][15] + count_array[14][15];
    p[5] = count_array[1][1] + count_array[2][2] + count_array[3][3] + count_array[4][4] + count_array[6][6] + count_array[7][7] + count_array[8][8] + count_array[9][9] + count_array[11][11] + count_array[12][12] + count_array[13][13] + count_array[14][14];
    p[8] = count_array[1][2] + count_array[1][3] + count_array[2][1] + count_array[2][3] + count_array[3][1] + count_array[3][2] + count_array[4][6] + count_array[4][7] + count_array[6][4] + count_array[6][7] + count_array[7][4] + count_array[7][6] + count_array[8][9] + count_array[8][11] + count_array[9][8] + count_array[9][11] + count_array[11][8] + count_array[11][9] + count_array[12][13] + count_array[12][14] + count_array[13][12] + count_array[13][14] + count_array[14][12] + count_array[14][13];
    p[6] = count_array[1][4] + count_array[2][8] + count_array[3][12] + count_array[4][1] + count_array[6][9] + count_array[7][13] + count_array[8][2] + count_array[9][6] + count_array[11][14] + count_array[12][3] + count_array[13][7] + count_array[14][11];
    p[4] = count_array[4][0] + count_array[8][0] + count_array[12][0] + count_array[1][5] + count_array[9][5] + count_array[13][5] + count_array[2][10] + count_array[6][10] + count_array[14][10] + count_array[3][15] + count_array[7][15] + count_array[11][15];
    p[10] = count_array[1][6] + count_array[1][7] + count_array[2][9] + count_array[2][11] + count_array[3][13] + count_array[3][14] + count_array[4][2] + count_array[4][3] + count_array[6][8] + count_array[6][11] + count_array[7][12] + count_array[7][14] + count_array[8][1] + count_array[8][3] + count_array[9][4] + count_array[9][7] + count_array[11][12] + count_array[11][13] + count_array[12][1] + count_array[12][2] + count_array[13][4] + count_array[13][6] + count_array[14][8] + count_array[14][9];
    p[9] = count_array[8][6] + count_array[12][7] + count_array[4][9] + count_array[12][11] + count_array[4][13] + count_array[8][14] + count_array[9][2] + count_array[13][3] + count_array[1][8] + count_array[13][11] + count_array[1][12] + count_array[9][14] + count_array[6][1] + count_array[14][3] + count_array[2][4] + count_array[14][7] + count_array[2][12] + count_array[6][13] + count_array[7][1] + count_array[11][2] + count_array[3][4] + count_array[11][6] + count_array[3][8] + count_array[7][9];
    p[11] = count_array[4][8] + count_array[4][12] + count_array[8][4] + count_array[8][12] + count_array[12][4] + count_array[12][8] + count_array[1][9] + count_array[1][13] + count_array[9][1] + count_array[9][13] + count_array[13][1] + count_array[13][9] + count_array[2][6] + count_array[2][14] + count_array[6][2] + count_array[6][14] + count_array[14][2] + count_array[14][6] + count_array[3][7] + count_array[3][11] + count_array[7][3] + count_array[7][11] + count_array[11][3] + count_array[11][7];
    p[13] = count_array[6][0] + count_array[7][0] + count_array[9][0] + count_array[11][0] + count_array[13][0] + count_array[14][0] + count_array[2][5] + count_array[3][5] + count_array[8][5] + count_array[11][5] + count_array[12][5] + count_array[14][5] + count_array[1][10] + count_array[3][10] + count_array[4][10] + count_array[7][10] + count_array[12][10] + count_array[13][10] + count_array[1][15] + count_array[2][15] + count_array[4][15] + count_array[6][15] + count_array[8][15] + count_array[9][15];
    p[14] = count_array[1][11] + count_array[1][14] + count_array[2][7] + count_array[2][13] + count_array[3][6] + count_array[3][9] + count_array[4][11] + count_array[4][14] + count_array[6][3] + count_array[6][12] + count_array[7][2] + count_array[7][8] + count_array[8][7] + count_array[8][13] + count_array[9][3] + count_array[9][12] + count_array[11][1] + count_array[11][4] + count_array[12][6] + count_array[12][9] + count_array[13][2] + count_array[13][8] + count_array[14][1] + count_array[14][4];


    sumps = 0;
    for (i=0; i<15; i++) sumps += p[i];
    if (fabs((1.0/count_noambigs)*sumps - 1.0) > 0.005) {
        printf("There was a problem counting site patterns ... exiting.");
        exit(1);
    }

    printf("Site pattern counts are:\n");
    for (i=0; i<15; i++) printf("%d ",p[i]);
    printf("\n\n");

    for (i=0; i<15; i++) q[i] = (double)p[i]/count_noambigs;
    printf("%d\n",count_noambigs);

	printf("Probs are:\n");
    for (i=0; i<15; i++) printf("%f ",q[i]); 
    printf("\n\n");


    if (issymmetric==1) {
	
	if (whichtau==1) tau = -1.0*log(pow(4*(1+m*t)*(q[0]/4-2*(q[1]/12+q[2]/12)+3*q[7]/12-2*q[12]/24+6*(q[3]/12+q[4]/12)-2*(q[5]/12+q[6]/12)-8*(q[8]/24+q[9]/24+q[10]/24+q[11]/24)-2*q[14]/24+6*q[13]/24),1.0/(2*m)));
	else if (whichtau==2) tau = -1.0*log(pow(4.0*(1.0+m*t)*(q[0]/4+6*(q[1]/12+q[2]/12)+3*q[7]/12+6*q[12]/24-2*(q[3]/12+q[4]/12)-2*(q[5]/12+q[6]/12)-8*(q[8]/24+q[9]/24+q[10]/24+q[11]/24)-2*q[14]/24-2*q[13]/24),1.0/(2*m)));
	else tau = -1.0*log(pow(4*(1+m*t)*(q[0]/4+2*(q[1]/12+q[2]/12)-q[7]/12-2*q[12]/24+2*(q[3]/12+q[4]/12)+2*(q[5]/12+q[6]/12)-2*q[14]/24-2*q[13]/24),1.0/(2*m)));
    }



  return tau;

}



/* Now identify quartets that we have selected that contain the node time of interest   */
/* Species to pick from are stored in ppNodeChildren                                    */
/* The sets we pick from determines the shape and location of the node time of interest */
/* Refer to Figure 1 of Chifman and Kubatko (2015) for shapes and taus                  */

double GetMomentEstimate(int mynode) {

  int i, j, k, l, i2, j2, k2, l2, tau_counter=0;
  int tax1, tax2, tax3, tax4;
  double tau=0.0;

  AllMyChildrenLeft(mynode);
  AllMyChildrenRight(mynode);
  if (mynode != ntaxa+1) AllMyChildren(mynode);
  else ppNodeChildren[2][0]=0;

  /* one left, one right, two parent */
  /* symmetric, tau1                 */

  if (ppNodeChildren[2][0]>1) {

  for (i=1; i<ppNodeChildren[0][0]+1; i++){
    for (j=1; j<ppNodeChildren[1][0]+1; j++) {
      for (k=1; k<ppNodeChildren[2][0]+1; k++) {
        for (l=k+1; l<ppNodeChildren[2][0]+1; l++) {
          printf("Quartet %d %d %d %d\n",ppNodeChildren[0][i],ppNodeChildren[1][j],ppNodeChildren[2][k],ppNodeChildren[2][l]);
          // now loop over choices within that quartet
          for (i2=0; i2<seq_counter[ppNodeChildren[0][i]-1]; i2++) {
            for (j2=0; j2<seq_counter[ppNodeChildren[1][j]-1]; j2++) {
              for (k2=0; k2<seq_counter[ppNodeChildren[2][k]-1]; k2++) {
                for (l2=0; l2<seq_counter[ppNodeChildren[2][l]-1]; l2++) {
		   printf("    Sequences are %d %d %d %d\n",ppSp_assign[ppNodeChildren[0][i]-1][i2],ppSp_assign[ppNodeChildren[1][j]-1][j2],ppSp_assign[ppNodeChildren[2][k]-1][k2],ppSp_assign[ppNodeChildren[2][l]-1][l2]);
		   if (TimeVec[mynode] < TimeVec[ppNodeChildren[2][k]] && TimeVec[mynode] < TimeVec[ppNodeChildren[2][l]]) tau += GetCompMOM(ppSp_assign[ppNodeChildren[0][i]-1][i2],ppSp_assign[ppNodeChildren[1][j]-1][j2],ppSp_assign[ppNodeChildren[2][k]-1][k2],ppSp_assign[ppNodeChildren[2][l]-1][l2],1,1);
		   else tau += GetCompMOM(ppSp_assign[ppNodeChildren[0][i]-1][i2],ppSp_assign[ppNodeChildren[1][j]-1][j2],ppSp_assign[ppNodeChildren[2][k]-1][k2],ppSp_assign[ppNodeChildren[2][l]-1][l2],1,2);
		   tau_counter++;
		   printf("tau estimate is %f\n\n\n",tau);
	        }
	      }
	    }
	  }
        }
      }
    }
  }

  }

  /* two left, one right, one parent */
  /* asymmetric, tau2                */

  if (ppNodeChildren[0][0]>1) {
  for (i=1; i<ppNodeChildren[0][0]+1; i++){
    for (j=i+1; j<ppNodeChildren[0][0]+1; j++) {
      for (k=1; k<ppNodeChildren[1][0]+1; k++) {
        for (l=1; l<ppNodeChildren[2][0]+1; l++) {
	  //printf("i=%d, j=%d, k=%d, l=%d\n",i,j,k,l);
          printf("Quartet %d %d %d %d\n",ppNodeChildren[0][i],ppNodeChildren[0][j],ppNodeChildren[1][k],ppNodeChildren[2][l]);  
          // now loop over choices within that quartet
          for (i2=0; i2<seq_counter[ppNodeChildren[0][i]-1]; i2++) {
            for (j2=0; j2<seq_counter[ppNodeChildren[0][j]-1]; j2++) {
              for (k2=0; k2<seq_counter[ppNodeChildren[1][k]-1]; k2++) {
                for (l2=0; l2<seq_counter[ppNodeChildren[2][l]-1]; l2++) {
                   printf("    Sequences are %d %d %d %d\n",ppSp_assign[ppNodeChildren[0][i]-1][i2],ppSp_assign[ppNodeChildren[0][j]-1][j2],ppSp_assign[ppNodeChildren[1][k]-1][k2],ppSp_assign[ppNodeChildren[2][l]-1][l2]);
                }
              }
            }
          }
        }                                                                                                                             
      }
    }
  }

  }
	

  /* one left, two right, one parent */
  /* asymmetric, tau2                */

  if (ppNodeChildren[1][0]>1) {
  for (i=1; i<ppNodeChildren[0][0]+1; i++){
    for (j=1; j<ppNodeChildren[1][0]+1; j++) {
      for (k=j+1; k<ppNodeChildren[1][0]+1; k++) {
        for (l=1; l<ppNodeChildren[2][0]+1; l++) {  
          //printf("i=%d, j=%d, k=%d, l=%d\n",i,j,k,l);
          printf("Quartet %d %d %d %d\n",ppNodeChildren[0][i],ppNodeChildren[1][j],ppNodeChildren[1][k],ppNodeChildren[2][l]);
          // now loop over choices within that quartet
          for (i2=0; i2<seq_counter[ppNodeChildren[0][i]-1]; i2++) {
            for (j2=0; j2<seq_counter[ppNodeChildren[1][j]-1]; j2++) {
              for (k2=0; k2<seq_counter[ppNodeChildren[1][k]-1]; k2++) {
                for (l2=0; l2<seq_counter[ppNodeChildren[2][l]-1]; l2++) {
                   printf("    Sequences are %d %d %d %d\n",ppSp_assign[ppNodeChildren[0][i]-1][i2],ppSp_assign[ppNodeChildren[1][j]-1][j2],ppSp_assign[ppNodeChildren[1][k]-1][k2],ppSp_assign[ppNodeChildren[2][l]-1][l2]);
                }
              }
            }
          }
        }
      }
    }
  }
       
  }



  /* two left, two right, no parent */
  /* symmetric, tau3                */

  if (ppNodeChildren[0][0]>1 && ppNodeChildren[1][0]>1) {
  for (i=1; i<ppNodeChildren[0][0]+1; i++){
    for (j=i+1; j<ppNodeChildren[0][0]+1; j++) {  
      for (k=1; k<ppNodeChildren[1][0]+1; k++) {
        for (l=k+1; l<ppNodeChildren[1][0]+1; l++) {
          //printf("i=%d, j=%d, k=%d, l=%d\n",i,j,k,l);                                                                               
          printf("Quartet %d %d %d %d\n",ppNodeChildren[0][i],ppNodeChildren[0][j],ppNodeChildren[1][k],ppNodeChildren[1][l]);
          // now loop over choices within that quartet
          for (i2=0; i2<seq_counter[ppNodeChildren[0][i]-1]; i2++) {
            for (j2=0; j2<seq_counter[ppNodeChildren[0][j]-1]; j2++) {
              for (k2=0; k2<seq_counter[ppNodeChildren[1][k]-1]; k2++) {
                for (l2=0; l2<seq_counter[ppNodeChildren[1][l]-1]; l2++) {
                   printf("    Sequences are %d %d %d %d\n",ppSp_assign[ppNodeChildren[0][i]-1][i2],ppSp_assign[ppNodeChildren[0][j]-1][j2],ppSp_assign[ppNodeChildren[1][k]-1][k2],ppSp_assign[ppNodeChildren[1][l]-1][l2]);
                }
              }  
            }                                                                                                                         
          }
        }
      }
    }                                                                                                                                 
  }

  }



  /* three left, one right, no parent */
  /* asymmetric, tau3                 */

  if (ppNodeChildren[0][0]>2) {
  for (i=1; i<ppNodeChildren[0][0]+1; i++){
    for (j=i+1; j<ppNodeChildren[0][0]+1; j++) {
      for (k=j+1; k<ppNodeChildren[0][0]+1; k++) {  
        for (l=1; l<ppNodeChildren[1][0]+1; l++) { 
          //printf("i=%d, j=%d, k=%d, l=%d\n",i,j,k,l);
          printf("Quartet %d %d %d %d\n",ppNodeChildren[0][i],ppNodeChildren[0][j],ppNodeChildren[0][k],ppNodeChildren[1][l]);
          // now loop over choices within that quartet
          for (i2=0; i2<seq_counter[ppNodeChildren[0][i]-1]; i2++) {
            for (j2=0; j2<seq_counter[ppNodeChildren[0][j]-1]; j2++) {
              for (k2=0; k2<seq_counter[ppNodeChildren[0][k]-1]; k2++) {
                for (l2=0; l2<seq_counter[ppNodeChildren[1][l]-1]; l2++) {
                   printf("    Sequences are %d %d %d %d\n",ppSp_assign[ppNodeChildren[0][i]-1][i2],ppSp_assign[ppNodeChildren[0][j]-1][j2],ppSp_assign[ppNodeChildren[0][k]-1][k2],ppSp_assign[ppNodeChildren[1][l]-1][l2]);
                }                                                                                                                     
              }                                                                                                                       
            }
          }
        }
      }
    }
  }

  }



  /* one left, three right, no parent */
  /* asymmetric, tau3                 */

  if (ppNodeChildren[1][0]>2) {
  for (i=1; i<ppNodeChildren[0][0]+1; i++){
    for (j=1; j<ppNodeChildren[1][0]+1; j++) {
      for (k=j+1; k<ppNodeChildren[1][0]+1; k++) {
        for (l=k+1; l<ppNodeChildren[1][0]+1; l++) {  
          //printf("i=%d, j=%d, k=%d, l=%d\n",i,j,k,l);
          printf("Quartet %d %d %d %d\n",ppNodeChildren[0][i],ppNodeChildren[1][j],ppNodeChildren[1][k],ppNodeChildren[1][l]);
          // now loop over choices within that quartet
          for (i2=0; i2<seq_counter[ppNodeChildren[0][i]-1]; i2++) {
            for (j2=0; j2<seq_counter[ppNodeChildren[1][j]-1]; j2++) {
              for (k2=0; k2<seq_counter[ppNodeChildren[1][k]-1]; k2++) {
                for (l2=0; l2<seq_counter[ppNodeChildren[1][l]-1]; l2++) {
                   printf("    Sequences are %d %d %d %d\n",ppSp_assign[ppNodeChildren[0][i]-1][i2],ppSp_assign[ppNodeChildren[1][j]-1][j2],ppSp_assign[ppNodeChildren[1][k]-1][k2],ppSp_assign[ppNodeChildren[1][l]-1][l2]);
                }                                                                                                                     
              }
            }
          }
        }
      }
    }
  }

  }


  printf("tau_counter is %d, tau is %f, Estimate of tau is %f\n",tau_counter, tau,tau/tau_counter);
  return tau/tau_counter;

}


double GetDerCompLik(int mynode) {
                
  int i, j, k, l;
  int tax1, tax2, tax3, tax4;
        
  /* one left, one right, two parent */
  /* symmetric, tau1                 */

  /* two left, one right, one parent */
  /* asymmetric, tau2                */

  /* one left, two right, one parent */                                                                                               
  /* asymmetric, tau2                */

  /* two left, two right, no parent */
  /* symmetric, tau3                */

  /* three left, one right, no parent */
  /* asymmetric, tau3                 */

  /* one left, three right, no parent */
  /* asymmetric, tau3                 */
  
  return 0.0;
 
}
