/*** Compute the S and A matrices to compute site pattern probs ***/

void ComputeAandS(double mytheta) {

   int i;
   double m = 4.0/3.0, t = 2.0*theta, mt = 8.0*theta/3.0;

   //compute the C matrix in the symmetric case - indexing starts at 1 to match R code
        // row 1
        for (i=1; i<10; i++) smat[1][i] = 1.0/256;
        // row 2
        smat[2][1]=smat[2][3]=smat[2][5]=smat[2][7] = 3.0/(256*(1+mt));
        smat[2][2]=smat[2][4]=smat[2][6]=smat[2][8]=smat[2][9]=-1.0/(256*(1+mt));
        // row 3
        smat[3][1]=smat[3][2]=smat[3][5]=smat[3][6]=3.0/(256*(1+mt));
        smat[3][3]=smat[3][4]=smat[3][7]=smat[3][8]=smat[3][9]=-1.0/(256*(1+mt));
        // row 4
        smat[4][1]=smat[4][5] = 9.0/(256*pow(1+mt,2));
        smat[4][2]=smat[4][3]=smat[4][6]=smat[4][7] = -3.0/(256*pow(1+mt,2));
        smat[4][4]=smat[4][8]=smat[4][9] = 1.0/(256*pow(1+mt,2));
        // row 5
        smat[5][1]=12.0/(256*(1+mt));
        smat[5][2]=smat[5][3]=smat[5][4]=4.0/(256*(1+mt));
        smat[5][5]=smat[5][6]=smat[5][7]=smat[5][9]=-4.0/(256*(1+mt));
        smat[5][8]=0.0;
        // row 6
        smat[6][1]=24.0/(256*(1+m*t)*(2+mt));
        smat[6][2]=smat[6][4]=smat[6][5]=smat[6][7]=-8.0/(256*(1+mt)*(2+mt));
        smat[6][3]=8.0/(256*(1+mt)*(2+mt));
        smat[6][6]=smat[6][9]=8.0/(256*(1+mt)*(2+mt));
        smat[6][8]=0.0;
        // row 7
        smat[7][1]=24.0/(256*(1+mt)*(2+mt));
        smat[7][2]=8.0/(256*(1+mt)*(2+mt));
        smat[7][3]=smat[7][4]=smat[7][5]=smat[7][6]=-8.0/(256*(1+mt)*(2+mt));
        smat[7][7]=smat[7][9]=8.0/(256*(1+mt)*(2+mt));
        smat[7][8]=0.0;
        // row 8
        smat[8][1]=48.0/(256*(1+m*t)*pow(2+mt,2));
        smat[8][2]=smat[8][3]=smat[8][5]=-16.0/(256*(1+mt)*pow(2+mt,2));
        smat[8][4]=smat[8][6]=smat[8][7]=16.0/(256*(1+mt)*pow(2+mt,2));
        smat[8][8]=0.0;
        smat[8][9]=-16.0/(256*(1+mt)*pow(2+mt,2));
        // row 9
        smat[9][1]=6.0*mt*(4+mt)*(4+3*mt)/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt));
        smat[9][2]=smat[9][3]=-2.0*mt*(4+mt)*(4+3*mt)/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt));
        smat[9][4]=mt*(32+40*mt+10*pow(m,2)*pow(t,2))/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt));
        smat[9][5]=2.0*mt*pow(4+mt,2)/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt));
        smat[9][6]=smat[9][7]=2.0*mt*mt*(4+mt)/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt));
        smat[9][8]=-mt*mt*(4+2*mt)/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt));
        smat[9][9]=2.0*mt*mt*mt/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt));

   //compute the C matrix in the symmetric case - indexing starts at 1 to match R code
        // row 1 
        amat[1][1] = 1.0/256; 
        amat[1][2] = 3.0/((256)*(1+mt)); 
        amat[1][3] = 6.0/(256*(1+mt));
        amat[1][4] = 12.0/(256*(1+mt)*(2+mt)); 
        amat[1][5] = 9.0/(256*(1+mt)); 
        amat[1][6] = 12.0/(256*(1+mt)*(2+mt)); 
        amat[1][7] = 9.0/(256*pow(1+mt,2)); 
        amat[1][8] = 24.0/(256*(1+mt)*(2+mt)); 
        amat[1][9] = 48.0/(256*(1+mt)*pow(2+mt,2)); 
        amat[1][10] = (6*mt*(4+mt)*(4+3*mt))/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt));

        // row 2 
        amat[2][1] = 1.0/256; 
        amat[2][2] = -1.0/(256*(1+mt)); 
        amat[2][3] = 2.0/(256*(1+mt));
        amat[2][4] = -4.0/(256*(1+mt)*(2+mt)); 
        amat[2][5] = 5.0/(256*(1+mt)); 
        amat[2][6] = -4.0/(256*(1+mt)*(2+mt)); 
        amat[2][7] = -3.0/(256*pow(1+mt,2)); 
        amat[2][8] = 8.0/(256*(1+mt)*(2+mt)); 
        amat[2][9] = -16.0/(256*(1+mt)*pow(2+mt,2)); 
        amat[2][10] = -(2*mt*(4+mt)*(4+3*mt))/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt)); 

        // row 3 
        amat[3][1] = 1.0/256; 
        amat[3][2] = 3.0/(256*(1+mt)); 
        amat[3][3] = -2.0/(256*(1+mt)); 
        amat[3][4] = -4.0/(256*(1+mt)*(2+mt)); 
        amat[3][5] = 5.0/(256*(1+mt)); 
        amat[3][6] = 12.0/(256*(1+mt)*(2+mt)); 
        amat[3][7] = -3.0/(256*pow(1+mt,2)); 
        amat[3][8] = -8.0/(256*(1+mt)*(2+mt)); 
        amat[3][9] = -16.0/(256*(1+mt)*pow(2+mt,2)); 
        amat[3][10] = -(2*mt*(4+mt)*(4+3*mt))/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt)); 

        // row 4 
        amat[4][1] = 1.0/256; 
        amat[4][2] = 3.0/(256*(1+mt)); 
        amat[4][3] = 6.0/(256*(1+mt)); 
        amat[4][4] = 12.0/(256*(1+mt)*(2+mt)); 
        amat[4][5] = -3.0/(256*(1+mt)); 
        amat[4][6] = -4.0/(256*(1+mt)*(2+mt)); 
        amat[4][7] = -3.0/(256*pow(1+mt,2)); 
        amat[4][8] = -8.0/(256*(1+mt)*(2+mt)); 
        amat[4][9] = -16.0/(256*(1+mt)*pow(2+mt,2)); 
        amat[4][10] = -(2*mt*(4+mt)*(4+3*mt))/(256*pow(1+mt,2)*pow(2+mt,2)*(3+mt)); 

        // row 5 
        amat[5][1] = 1.0/256; 
        amat[5][2] = 3.0/((256)*(1+mt)); 
        amat[5][3] = -2.0/((256)*(1+mt)); 
        amat[5][4] = -4.0/((256)*(1+mt)*(2+mt)); 
        amat[5][5] = 1.0/((256)*(1+mt)); 
        amat[5][6] = -4.0/((256)*(1+mt)*(2+mt)); 
        amat[5][7] = 9.0/((256)*pow(1+mt,2)); 
        amat[5][8] = -8.0/((256)*(1+mt)*(2+mt)); 
        amat[5][9] = -16.0/((256)*(1+mt)*pow(2+mt,2)); 
        amat[5][10] = (2*mt*pow(4+mt,2))/((256)*pow(1+mt,2)*pow(2+mt,2)*(3+mt)); 

        // row 6 
        amat[6][1] = 1.0/256; 
        amat[6][2] = -1.0/((256)*(1+mt)); 
        amat[6][3] = 2.0/((256)*(1+mt)); 
        amat[6][4] = -4.0/((256)*(1+mt)*(2+mt)); 
        amat[6][5] = 1.0/((256)*(1+mt)); 
        amat[6][6] = -4.0/((256)*(1+mt)*(2+mt)); 
        amat[6][7] = 1.0/((256)*pow(1+mt,2)); 
        amat[6][8] = -8.0/((256)*(1+mt)*(2+mt)); 
        amat[6][9] = 16.0/((256)*(1+mt)*pow(2+mt,2)); 
        amat[6][10] = mt*(2*16+40*mt+(10.0)*mt*mt)/((256)*pow(1+mt,2)*pow(2+mt,2)*(3+mt)); 

        // row 7 
        amat[7][1] = 1.0/256; 
        amat[7][2] = -1.0/((256)*(1+mt)); 
        amat[7][3] = -2.0/((256)*(1+mt));
        amat[7][4] = 4.0/((256)*(1+mt)*(2+mt)); 
        amat[7][5] = 1.0/((256)*(1+mt)); 
        amat[7][6] = 4.0/((256)*(1+mt)*(2+mt)); 
        amat[7][7] = -3.0/((256)*pow(1+mt,2)); 
        amat[7][8] = -8.0/((256)*(1+mt)*(2+mt)); 
        amat[7][9] = 16.0/((256)*(1+mt)*pow(2+mt,2)); 
        amat[7][10] = 2.0*mt*mt*(4+mt)/((256)*pow(1+mt,2)*pow(2+mt,2)*(3+mt)); 

        // row 8 
        amat[8][1] = 1.0/256;
        amat[8][2] = 3.0/((256)*(1+mt)); 
        amat[8][3] = -2.0/((256)*(1+mt)); 
        amat[8][4] = -4.0/((256)*(1+mt)*(2+mt)); 
        amat[8][5] = -3.0/((256)*(1+mt)); 
        amat[8][6] = -4.0/((256)*(1+mt)*(2+mt)); 
        amat[8][7] = -3.0/((256)*pow(1+mt,2)); 
        amat[8][8] = 8.0/((256)*(1+mt)*(2+mt)); 
        amat[8][9] = 16.0/((256)*(1+mt)*pow(2+mt,2)); 
        amat[8][10] = 2.0*mt*mt*(4+mt)/((256)*pow(1+mt,2)*pow(2+mt,2)*(3+mt)); 

        // row 9 
        amat[9][1] = 1.0/256; 
        amat[9][2] = -1.0/((256)*(1+mt)); 
        amat[9][3] = -2.0/((256)*(1+mt)); 
        amat[9][4] = 4.0/((256)*(1+mt)*(2+mt)); 
        amat[9][5] = 1.0/((256)*(1+mt)); 
        amat[9][6] = -4.0/((256)*(1+mt)*(2+mt)); 
        amat[9][7] = 1.0/((256)*pow(1+mt,2)); 
        amat[9][8] = 0.0;
        amat[9][9] = 0.0; 
        amat[9][10] = -mt*mt*(4+2*mt)/((256)*pow(1+mt,2)*pow(2+mt,2)*(3+mt));

        // row 10 
        amat[10][1] = 1.0/256;
        amat[10][2] = -1.0/((256)*(1+mt));
        amat[10][3] = 2.0/((256)*(1+mt)); 
        amat[10][4] = -4.0/((256)*(1+mt)*(2+mt)); 
        amat[10][5] = -3.0/((256)*(1+mt)); 
        amat[10][6] = 4.0/((256)*(1+mt)*(2+mt)); 
        amat[10][7] = 1.0/((256)*pow(1+mt,2)); 
        amat[10][8] = 0.0; 
        amat[10][9] = 0.0; 
        amat[10][10] = -mt*mt*(4+2*mt)/((256)*pow(1+mt,2)*pow(2+mt,2)*(3+mt)); 

        // row 11 
        amat[11][1] = 1.0/256; 
        amat[11][2] = -1.0/((256)*(1+mt)); 
        amat[11][3] = -2.0/((256)*(1+mt)); 
        amat[11][4] = 4.0/((256)*(1+mt)*(2+mt)); 
        amat[11][5] = -3.0/((256)*(1+mt)); 
        amat[11][6] = 4.0/((256)*(1+mt)*(2+mt)); 
        amat[11][7] = 1.0/((256)*pow(1+mt,2)); 
        amat[11][8] = 8.0/((256)*(1+mt)*(2+mt)); 
        amat[11][9] = -16.0/((256)*(1+mt)*pow(2+mt,2)); 
        amat[11][10] = 2.0*mt*mt*mt/((256)*pow(1+mt,2)*pow(2+mt,2)*(3+mt)); 

}

// Function to swap two elements
void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Sorting function for 4 elements
void sort4(int arr[]) {
    // Compare and sort the 4 elements explicitly
    if (arr[0] > arr[1]) swap(&arr[0], &arr[1]);
    if (arr[0] > arr[2]) swap(&arr[0], &arr[2]);
    if (arr[0] > arr[3]) swap(&arr[0], &arr[3]);
    if (arr[1] > arr[2]) swap(&arr[1], &arr[2]);
    if (arr[1] > arr[3]) swap(&arr[1], &arr[3]);
    if (arr[2] > arr[3]) swap(&arr[2], &arr[3]);
}



// Returns weights for a base code (0-3 = unambiguous, 4+ = ambiguous)
// Code below modified from Claude, 5/16/26

BaseWeights get_base_weights(int base_code) {
    BaseWeights bw = {{0.0, 0.0, 0.0, 0.0}};
    switch (base_code) {
        // Unambiguous
        case 0: bw.weight[0] = 1.0; break;  // A
        case 1: bw.weight[1] = 1.0; break;  // C
        case 2: bw.weight[2] = 1.0; break;  // G
        case 3: bw.weight[3] = 1.0; break;  // T
        // 2-fold ambiguities
        case 6:  bw.weight[0]=bw.weight[2]=0.5; break;  // R = A/G
        case 9:  bw.weight[1]=bw.weight[3]=0.5; break;  // Y = C/T
        case 5:  bw.weight[0]=bw.weight[1]=0.5; break;  // M = A/C
        case 10:  bw.weight[2]=bw.weight[3]=0.5; break; // K = G/T
        case 7:  bw.weight[0]=bw.weight[3]=0.5; break;  // W = A/T
        case 8:  bw.weight[1]=bw.weight[2]=0.5; break;  // S = C/G
        // 3-fold ambiguities
        case 14: bw.weight[0]=bw.weight[1]=bw.weight[2]=1.0/3.0; break;  // V = A/C/G
        case 13: bw.weight[0]=bw.weight[1]=bw.weight[3]=1.0/3.0; break;  // H = A/C/T
        case 12: bw.weight[0]=bw.weight[2]=bw.weight[3]=1.0/3.0; break;  // D = A/G/T
        case 11: bw.weight[1]=bw.weight[2]=bw.weight[3]=1.0/3.0; break;  // B = C/G/T
        // 4-fold ambiguity
        case 15: bw.weight[0]=bw.weight[1]=bw.weight[2]=bw.weight[3]=0.25; break;  // N
        default: break;  // gap or unknown: zero weight (site skipped) - this includes case=4
    }
    return bw;
}



/***  Function to search through ppTwoRow to find parent ***/
/***  of a specified node.                               ***/

int find_parent(int target) {

int i, j, parent=0;

 for (j=0; j<2; j++) {
   for (i=0; i<ntaxa; i++) {
     if (ppTwoRow[j][i]==target) {
      parent = i+ntaxa+1;
      break;
     }
   }
 }

 //if (parent==0) { printf("Could not find parent in find_parent for target = %d\n",target); }
 return(parent);

}


/*** Function that returns 1 if a cherry and 0 if not in the complete tree ***/

int IsCherry(int q, int r) {

   if (find_parent(q) == find_parent(r)) return 1;
   else return 0;

}

/*** Function that returns 1 if the first two arguments form a cherry and 0 if not in the quartet tree ***/

int IsCherryQuartet(int q, int r, int w, int x) {

   //printf("in is cherry quartet %d %d %d %d",q,r,w,x);

   if (find_parent(q) == find_parent(r)) return 1; // easy -- q and r are a cherry in the full tree
   else {
      //printf("where I want to be\n");
    if ((ppMatrix[q][r] < ppMatrix[q][w]) && (ppMatrix[q][r] < ppMatrix[r][w]) && (ppMatrix[q][r] < ppMatrix[q][x]) && (ppMatrix[q][r] < ppMatrix[r][x])) return 1;
    else {
       //printf("leaving is cherry quartet\n");
       //printf("%f %f %f %f %f \n",ppMatrix[q][r],ppMatrix[q][w],ppMatrix[r][w],ppMatrix[q][x],ppMatrix[r][x]);
        return 0;
        }
   }

}



/*** Function to count number of cherries in full tree ***/

int CountCherries(int q, int r, int s, int t) {

   int sum = IsCherry(q,r) + IsCherry(q,s) + IsCherry(q,t) + IsCherry(r,s) + IsCherry(r,t) + IsCherry(s,t);
   //printf("number of cherries is %d\n",sum);
   return sum;

}


/***  Function to find the generation of the input node   ***/

int find_gen(int current_node) {

  int gen_counter=0, parent;

    parent = current_node;

    while (parent != ntaxa+1) {

      parent = find_parent(parent);
      gen_counter++;

    }

    return gen_counter;

}


/***  Function to construct the IND matrix                        ***/
/*** Instead of number of nodes, store time of parent in ppMatrix ***/

void make_indmat() {

  int q, r;
  int parent1, parent2, check=0;
  int dist1=0, dist2=1;
  int time_parent;

  for (q=1; q<ntaxa+1; q++) {

    for (r=q+1; r<ntaxa+1; r++) {

      parent1 = find_parent(q);
      parent2 = find_parent(r);

      while (parent1 != ntaxa+1 && check != 1) {

        while (parent2 != ntaxa+1 && check != 1) {

          if (parent1 == parent2) {

            check = 1;
            break;

          }

          parent2 = find_parent(parent2);
          if (check != 1) dist2++;

        }

        if (check == 1) break;

        else { 

	  parent1 = find_parent(parent1);
          dist1++;
          dist2 = 1;
          parent2 = find_parent(r);
          
        }

      }

      //ppMatrix[q][r] = dist1 + dist2;
      //printf("For taxa %d and %d, parent1 is %d, parent2 is %d\n",q,r,parent1,parent2); 
      if (parent1>parent2) time_parent = parent2;
      else time_parent = parent1;
      ppMatrix[q][r] = TimeVec[time_parent];

      if ((parent1 == ntaxa+1 || parent2 == ntaxa+1) && check == 0) {

        //ppMatrix[q][r] = find_gen(q) + find_gen(r) - 1;
	ppMatrix[q][r] = TimeVec[ntaxa+1];
        check = 1;

      }
 
      check = 0;
      dist1 = 0;
      dist2 = 1;

    }

    check = 0;
    dist1 = 0;
    dist2 = 1;


  }

  /* fill the upper diagonal */
  for (q=1; q<ntaxa+1; q++) {
     for (r=q+1; r<ntaxa+1; r++) {
	ppMatrix[r][q] = ppMatrix[q][r];
     }
  }

  /*printf("\n\n ppMatrix is:\n");
  for (q=1; q<ntaxa+1; q++) {
	for (r=q+1; r<ntaxa+1; r++) {
		printf("%f ",ppMatrix[q][r]);
	}
	printf("\n");
  }*/

}


/*** Find the MRCA of two nodes in the complete tree ***/

int FindMRCA(int q, int r) {

  int parent1, parent2, check=0;
  int mrca;

  parent1 = find_parent(q);
  parent2 = find_parent(r);

  while (parent1 != ntaxa+1 && check != 1) {
    while (parent2 != ntaxa+1 && check != 1) {
          if (parent1 == parent2) {
            check = 1;
            break;
          }
          parent2 = find_parent(parent2);
    }
    if (check == 1) break;
    else { 
          parent1 = find_parent(parent1);
          parent2 = find_parent(r);
    }
  }

  if (parent1>parent2) mrca = parent2;
  else mrca = parent1;

  return mrca;

}

/*** Populates ppNodeChildrenLeftQuart row 0 with left children ***/
/*** of node mynode; [0][0] entry is number of children         ***/

void AllMyChildrenLeftQuartLeft(int mynode) {

  int i, j=0, currentpar;
  int leftpar = ppTwoRow[0][mynode-(ntaxa+1)];

  for (i=1; i<ntaxa+1; i++)  {

        currentpar = i;
        while (currentpar != leftpar && currentpar != ntaxa+1) currentpar = find_parent(currentpar);
        if (currentpar!=ntaxa+1) {
                j++;
                ppNodeChildrenLeftQuart[0][j] = i;
        }
  }
  ppNodeChildrenLeftQuart[0][0] = j;

}

/*** Populates ppNodeChildrenLeftQuart row 1 with right children ***/
/*** of node mynode; [1][0] entry is number of children          ***/

void AllMyChildrenRightQuartetLeft(int mynode) {

  int i, j=0, currentpar;
  int rightpar = ppTwoRow[1][mynode-(ntaxa+1)];

  for (i=1; i<ntaxa+1; i++)  {

        currentpar = i;
        while (currentpar != rightpar && currentpar != ntaxa+1) currentpar = find_parent(currentpar);
        if (currentpar != ntaxa+1) {
                j++;
                ppNodeChildrenLeftQuart[1][j] = i;
        }
  }
  ppNodeChildrenLeftQuart[1][0] = j;

}


/*** Populates ppNodeChildrenRightQuart row 0 with left children ***/
/*** of node mynode; [0][0] entry is number of children         ***/
  
void AllMyChildrenLeftQuartRight(int mynode) {
  
  int i, j=0, currentpar;
  int leftpar = ppTwoRow[0][mynode-(ntaxa+1)];
        
  for (i=1; i<ntaxa+1; i++)  {
                
        currentpar = i;
        while (currentpar != leftpar && currentpar != ntaxa+1) currentpar = find_parent(currentpar);
        if (currentpar!=ntaxa+1) {
                j++;
                ppNodeChildrenRightQuart[0][j] = i;
        }
  }  
  ppNodeChildrenRightQuart[0][0] = j;

}


/*** Populates ppNodeChildrenRightQuart row 1 with right children ***/
/*** of node mynode; [1][0] entry is number of children          ***/
  
void AllMyChildrenRightQuartetRight(int mynode) {
  
  int i, j=0, currentpar;
  int rightpar = ppTwoRow[1][mynode-(ntaxa+1)];
        
  for (i=1; i<ntaxa+1; i++)  {
                
        currentpar = i;
        while (currentpar != rightpar && currentpar != ntaxa+1) currentpar = find_parent(currentpar);
        if (currentpar != ntaxa+1) { 
                j++;
                ppNodeChildrenRightQuart[1][j] = i;
        }
  }  
  ppNodeChildrenRightQuart[1][0] = j;

}


/*** This function finds the children of the parent of mynode on the "other" ***/
/*** side; populates row 2 of ppNodeChildrenRightQuart; [2][0] entry is    ***/
/*** number of children                                                      ***/

void AllMyChildrenQuartet(int mynode) {

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
                ppNodeChildrenRightQuart[2][j] = i;
        }

  }
  ppNodeChildrenRightQuart[2][0] = j;
}


/*** Function to find all duplicate quartets and add ***/
/*** their site pattern counts; uses current values  ***/
/*** of ppTwoRowQuart and TimeVecQuart               ***/

void FindDupQuarts(int is_symm, double spcounts[15]) {

  int i, j, k, l, m;
  int mrca, mrca3, parent1, parent2, leftnode, sideind;
  int mytaxarr[4];
  double dupspcounts[15];
 
  /* set ppNodeChildrenLeftQuart and ppNodeChildrenRightQuart to 0 */
  for (i=0; i<ntaxa+1; i++) {
	ppNodeChildrenLeftQuart[0][i] = 0;
	ppNodeChildrenLeftQuart[1][i] = 0;
	ppNodeChildrenRightQuart[0][i] = 0;
	ppNodeChildrenRightQuart[1][i] = 0;
	ppNodeChildrenRightQuart[2][i] = 0;
  }
  /* initialize spcounts */
  for (i=0; i<15; i++) { spcounts[i] = 0.0; dupspcounts[i] = 0.0;}

  
  if (is_symm==2) mrca = FindMRCA(ppTwoRowQuart[0][ntaxa+2],ppTwoRowQuart[0][ntaxa+3]);
  else if (is_symm==1) mrca = FindMRCA(ppTwoRowQuart[1][ntaxa+1],ppTwoRowQuart[0][ntaxa+3]);
       else { 
	printf("At the time of error, the tree and times are:\n");
	for (i=0; i<ntaxa-1; i++) {
                printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
                printf("%f %f\n ",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
                }
        printf("\n");

	printf("Error in FindDupQuarts .... exiting.\n\n"); 
	exit(1);
	}

  // case 1: symmetric
  // fill vectors of descendents
  // LeftQuart = left descendents of left cherry in row 0, right descendents in row 1
  // RightQuart = left descendents of right cherry in row 0, right descendents in row 1

  if (is_symm==2) {

  	parent1 = FindMRCA(ppTwoRowQuart[0][ntaxa+2],ppTwoRowQuart[1][ntaxa+2]);
  	AllMyChildrenLeftQuartLeft(parent1);
  	AllMyChildrenRightQuartetLeft(parent1);

  	parent2 = FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+3]);
  	AllMyChildrenLeftQuartRight(parent2);
  	AllMyChildrenRightQuartetRight(parent2);

  	//for (i=0; i<ntaxa; i++) printf("%d %d %d %d\n",ppNodeChildrenLeftQuart[0][i],ppNodeChildrenLeftQuart[1][i],ppNodeChildrenRightQuart[0][i],ppNodeChildrenRightQuart[1][i]);

   	/* now get counts */
	for (i=0; i<ppNodeChildrenLeftQuart[0][0]; i++) {
	  for (j=0; j<ppNodeChildrenLeftQuart[1][0]; j++) {
	    for (k=0;k<ppNodeChildrenRightQuart[0][0]; k++) {
	      for (l=0; l<ppNodeChildrenRightQuart[1][0]; l++) {
		mytaxarr[0] = ppNodeChildrenLeftQuart[0][i+1];
		mytaxarr[1] = ppNodeChildrenLeftQuart[1][j+1];
		mytaxarr[2] = ppNodeChildrenRightQuart[0][k+1];
		mytaxarr[3] = ppNodeChildrenRightQuart[1][l+1];
		CountQuartetSitePatterns(mytaxarr[0],mytaxarr[1],mytaxarr[2],mytaxarr[3],dupspcounts);
		//CountQuartetSitePatterns(ppNodeChildrenLeftQuart[0][i+1],ppNodeChildrenLeftQuart[1][j+1],ppNodeChildrenRightQuart[0][k+1],ppNodeChildrenRightQuart[1][l+1],dupspcounts);
		for (m=0; m<15; m++) spcounts[m] = spcounts[m]+dupspcounts[m];
		sort4(mytaxarr);
		qvec[(int)subset_to_index(mytaxarr[0]-1,mytaxarr[1]-1,mytaxarr[2]-1,mytaxarr[3]-1)]=1;
		//printf("Marking quartet %d %d %d %d done\n",mytaxarr[0],mytaxarr[1],mytaxarr[2],mytaxarr[3]);
	      }
            }
          }
	}
        //printf("The overall counts are:\n");
	//for (m=0; m<15; m++) printf("%f ",dupspcounts[m]);
	//printf("\n");
  }

  // case 2: asymmetric
  // fill vectors of descendents
  // LeftQuart = descendents on "long" side of root in row 0
  // RightQuart = descendents of left side of cherry in row 0
  // RightQuart = descendents of right side of cherry in row 1
  // RightQuart= descendents of sister of cherry in row 2
  
  else if (is_symm==1) {

    	parent1 = ppTwoRowQuart[1][ntaxa+1];
	parent2 = find_parent(parent1);
	if (parent2 == mrca) {
		if (ppTwoRow[0][mrca-(ntaxa+1)] == parent1) sideind=0;  //left
		else sideind=1; //right
	}
	else {
		while (parent2 != mrca) {

			parent1 = parent2;
			parent2 = find_parent(parent1);
		}
		if (parent1 == ppTwoRow[0][parent2-(ntaxa+1)]) sideind=0; //left
		else sideind=1; //right	
	}

	if (sideind==0) AllMyChildrenLeftQuartLeft(parent2);
	else AllMyChildrenRightQuartetLeft(parent2);   // now the left side of the mrca is filled

	// now fill ppNodeChildrenQuartRight for the other 3 taxa -- cherry first
	parent1 = FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+3]);
	AllMyChildrenLeftQuartRight(parent1);
        AllMyChildrenRightQuartetRight(parent1);

	// figure out which side the cherry was on in the complete tree
	mrca3 = FindMRCA(ppTwoRowQuart[1][ntaxa+2],ppTwoRowQuart[1][ntaxa+3]);
	parent2 = find_parent(parent1);
        if (parent2 == mrca3) AllMyChildrenQuartet(parent1);
        else {
                while (parent2 != mrca3) {

                        parent1 = parent2;
                        parent2 = find_parent(parent1);
                }
		AllMyChildrenQuartet(parent1);
        }

	//for (i=0; i<ntaxa; i++) printf("%d %d %d %d %d\n",ppNodeChildrenLeftQuart[0][i],ppNodeChildrenLeftQuart[1][i],ppNodeChildrenRightQuart[0][i],ppNodeChildrenRightQuart[1][i],ppNodeChildrenRightQuart[2][i]);

	/* now get counts */
        for (i=0; i<ppNodeChildrenLeftQuart[sideind][0]; i++) {
          for (j=0; j<ppNodeChildrenRightQuart[2][0]; j++) {
            for (k=0;k<ppNodeChildrenRightQuart[0][0]; k++) {
              for (l=0; l<ppNodeChildrenRightQuart[1][0]; l++) {
		mytaxarr[0] = ppNodeChildrenLeftQuart[sideind][i+1];
		mytaxarr[1] = ppNodeChildrenRightQuart[2][j+1];
		mytaxarr[2] = ppNodeChildrenRightQuart[0][k+1];
		mytaxarr[3] = ppNodeChildrenRightQuart[1][l+1];
		CountQuartetSitePatterns(mytaxarr[0],mytaxarr[1],mytaxarr[2],mytaxarr[3],dupspcounts);
                //CountQuartetSitePatterns(ppNodeChildrenLeftQuart[sideind][i+1],ppNodeChildrenRightQuart[2][j+1],ppNodeChildrenRightQuart[0][k+1],ppNodeChildrenRightQuart[1][l+1],dupspcounts);
                for (m=0; m<15; m++) spcounts[m] = spcounts[m]+dupspcounts[m];
		sort4(mytaxarr);
                qvec[(int)subset_to_index(mytaxarr[0]-1,mytaxarr[1]-1,mytaxarr[2]-1,mytaxarr[3]-1)]=1;
                //printf("Marking quartet %d %d %d %d done\n",mytaxarr[0],mytaxarr[1],mytaxarr[2],mytaxarr[3]);
              }
            }
          }
        }
        //printf("The overall counts are:\n");
        //for (m=0; m<15; m++) printf("%f ",dupspcounts[m]);
        //printf("\n");

  }

}


/*** Function to pull off quartet trees ***/
/*** and store them in ppTwoRowQuart    ***/
/*** Times are stored in TimeVecQuart   ***/

int GetQuartetTree(int tax1, int tax2, int tax3, int tax4, int OrderedVec[5]) {

   int i, sum=0, leftover1, leftover2;
   int resvec[5];  // element 0 = is_symmetric; elements 1-4 gives the ordering of the taxa
   for (i=0; i<5;i++) resvec[i]=0;

   make_indmat();

   // tax1 and tax2
   if (IsCherryQuartet(tax1,tax2,tax3,tax4) == 1) {
	ppTwoRowQuart[0][ntaxa+1+2] = tax1;
	ppTwoRowQuart[1][ntaxa+1+2] = tax2;
	leftover1 = tax3; leftover2 = tax4;
	sum = 1;
	resvec[1]=3; resvec[2]=4;
   }
   // tax1 and tax3
   if (IsCherryQuartet(tax1,tax3,tax2,tax4) == 1) {
	if (sum == 0) {
		ppTwoRowQuart[0][ntaxa+1+2] = tax1;
        	ppTwoRowQuart[1][ntaxa+1+2] = tax3;
		leftover1 = tax2; leftover2 = tax4;
        	sum = 1; 
		resvec[1]=3; resvec[3]=4;
        }
	else if (sum == 1) {
		ppTwoRowQuart[0][ntaxa+1+1] = tax1;
        	ppTwoRowQuart[1][ntaxa+1+1] = tax3;	
		sum = 2;
		resvec[1]=3; resvec[2]=1; resvec[3]=4; resvec[4]=2;
        } 
    }
   if (sum<2) {
	if (IsCherryQuartet(tax1,tax4,tax2,tax3) == 1) {
		if (sum == 0) {
                	ppTwoRowQuart[0][ntaxa+1+2] = tax1;
                	ppTwoRowQuart[1][ntaxa+1+2] = tax4;
			leftover1 = tax2; leftover2 = tax3;
                	sum = 1;
			resvec[1]=3; resvec[4]=4;
        	}
        	else if (sum == 1) {
                	ppTwoRowQuart[0][ntaxa+1+1] = tax1;
                	ppTwoRowQuart[1][ntaxa+1+1] = tax4; 
                	sum = 2;
			resvec[1]=3; resvec[2]=1; resvec[3]=2; resvec[4]=4;
        	}
    	}
    }
    if (sum<2) {
        if (IsCherryQuartet(tax2,tax3,tax1,tax4) == 1) {
                if (sum == 0) {
                        ppTwoRowQuart[0][ntaxa+1+2] = tax2;
                        ppTwoRowQuart[1][ntaxa+1+2] = tax3;
			leftover1 = tax1; leftover2 = tax4;
			sum = 1;
			resvec[2]=3; resvec[3]=4;
                }
        	else if (sum == 1) {
                        ppTwoRowQuart[0][ntaxa+1+1] = tax2;
                        ppTwoRowQuart[1][ntaxa+1+1] = tax3;
			sum = 2;
			resvec[1]=1; resvec[2]=3; resvec[3]=4; resvec[4]=2;
                }
        }
    }
    if (sum<2) {
        if (IsCherryQuartet(tax2,tax4,tax1,tax3) == 1) {
                if (sum == 0) {
                        ppTwoRowQuart[0][ntaxa+1+2] = tax2;
                        ppTwoRowQuart[1][ntaxa+1+2] = tax4;
			leftover1 = tax1; leftover2 = tax3;
                        sum = 1;
			resvec[2]=3; resvec[4]=4;
                }
                else if (sum == 1) {
                        ppTwoRowQuart[0][ntaxa+1+1] = tax2;
                        ppTwoRowQuart[1][ntaxa+1+1] = tax4;
			sum = 2;
			resvec[1]=1; resvec[2]=3; resvec[3]=2; resvec[4]=4;
                }
        }
    }
    if (sum<2) {
        if (IsCherryQuartet(tax3,tax4,tax1,tax2) == 1) {
                if (sum == 0) {
                        ppTwoRowQuart[0][ntaxa+1+2] = tax3;
                        ppTwoRowQuart[1][ntaxa+1+2] = tax4;
			leftover1 = tax1; leftover2 = tax2;
			sum = 1;
			resvec[3]=3; resvec[4]=4;
                }
                else if (sum == 1) {
                        ppTwoRowQuart[0][ntaxa+1+1] = tax3;
                        ppTwoRowQuart[1][ntaxa+1+1] = tax4;
			sum = 2;
			resvec[1]=1; resvec[2]=2; resvec[3]=3; resvec[4]=4;
                }
        }
    }

    if (sum == 2) {
	ppTwoRowQuart[0][ntaxa+1+0] = ntaxa+2; //6;
	ppTwoRowQuart[1][ntaxa+1+0] = ntaxa+3; //7;
    }
    else if (sum == 1) {
	ppTwoRowQuart[0][ntaxa+1+0] = ntaxa+2; //6;
	ppTwoRowQuart[0][ntaxa+1+1] = ntaxa+3; //7;
	if (ppMatrix[ppTwoRowQuart[0][ntaxa+1+2]][leftover1] < ppMatrix[ppTwoRowQuart[0][ntaxa+1+2]][leftover2]) {
		ppTwoRowQuart[1][ntaxa+1+0] = leftover2;
		ppTwoRowQuart[1][ntaxa+1+1] = leftover1;
		if (resvec[1] == 0) {
			if (tax1==leftover2) {
				resvec[1] = 1;
				if (resvec[2]==0) resvec[2] = 2;
				else if (resvec[3]==0) resvec[3] = 2;
				else if (resvec[4]==0) resvec[4] = 2;
		        }
                        else {
                                resvec[1] = 1;
                                if (resvec[2]==0) resvec[2] = 2;
                                else if (resvec[3]==0) resvec[3] = 2;
                                else if (resvec[4]==0) resvec[4] = 2; 
                        }
                }
         	else if (resvec[2] == 0) {
			if (tax2==leftover2) 	{
				resvec[2] = 1;
                                if (resvec[3]==0) resvec[3] = 2;
                                else if (resvec[4]==0) resvec[4] = 2;
                        }
                        else {
                                resvec[2] = 2;
                                if (resvec[3]==0) resvec[3] = 1;
                                else if (resvec[4]==0) resvec[4] = 1; 
                	}
		}
		else if (resvec[3] == 0) {
			if (tax3==leftover2) { resvec[3] = 1; resvec[4]=2; }
			else { resvec[3] = 2; resvec[4] = 1;}
		}
	}
	else {
		ppTwoRowQuart[1][ntaxa+1+0] = leftover1;
                ppTwoRowQuart[1][ntaxa+1+1] = leftover2;
		if (resvec[1] == 0) {
                        if (tax1==leftover1) {
                                resvec[1] = 1;
                                if (resvec[2]==0) resvec[2] = 2;
                                else if (resvec[3]==0) resvec[3] = 2;
                                else if (resvec[4]==0) resvec[4] = 2;
                        }                       
                        else {
                                resvec[1] = 2;
                                if (resvec[2]==0) resvec[2] = 1;
                                else if (resvec[3]==0) resvec[3] = 1;
                                else if (resvec[4]==0) resvec[4] = 1;
                        }
                }
		else if (resvec[2] == 0) {
                        if (tax2==leftover1)  {
                                resvec[2] = 1;
                                if (resvec[3]==0) resvec[3] = 2;
                                else if (resvec[4]==0) resvec[4] = 2;
                        }
                        else {
                                resvec[2] = 2;
                                if (resvec[3]==0) resvec[3] = 1;
                                else if (resvec[4]==0) resvec[4] = 1;
                        }
                }
		else if (resvec[3] == 0) {
                        if (tax3==leftover1) { resvec[3] = 1; resvec[4] = 2; }
                        else { resvec[3] = 2; resvec[4] = 1; }
                }
	}
    }

   resvec[0] = sum;
   OrderedVec[0] = sum;
   //printf("%d %d %d %d %d\n",resvec[0],resvec[1],resvec[2],resvec[3],resvec[4]); 

   if (resvec[1]==1) OrderedVec[1]=tax1; 
        else if (resvec[1]==2) OrderedVec[2]=tax1;
        else if (resvec[1]==3) OrderedVec[3]=tax1; 
        else if (resvec[1]==4) OrderedVec[4]=tax1; 
        else printf("There was a problem determining order in CountQuartetSitePatterns 1 - %d %d %d %d %d %d %d %d \n\n",tax1,tax2,tax3,tax4,resvec[1],resvec[2],resvec[3],resvec[4]);
    if (resvec[2]==1) OrderedVec[1]=tax2; 
        else if (resvec[2]==2) OrderedVec[2]=tax2;
        else if (resvec[2]==3) OrderedVec[3]=tax2;
        else if (resvec[2]==4) OrderedVec[4]=tax2;
        else printf("There was a problem determining order in CountQuartetSitePatterns 2 -  %d %d %d %d %d %d %d %d \n\n",tax1,tax2,tax3,tax4,resvec[1],resvec[2],resvec[3],resvec[4]);
    if (resvec[3]==1) OrderedVec[1]=tax3;
        else if (resvec[3]==2) OrderedVec[2]=tax3;
        else if (resvec[3]==3) OrderedVec[3]=tax3;
        else if (resvec[3]==4) OrderedVec[4]=tax3;
        else printf("There was a problem determining order in CountQuartetSitePatterns 3 -  %d %d %d %d %d %d %d %d \n\n",tax1,tax2,tax3,tax4,resvec[1],resvec[2],resvec[3],resvec[4]);
    if (resvec[4]==1) OrderedVec[1]=tax4;
        else if (resvec[4]==2) OrderedVec[2]=tax4;
        else if (resvec[4]==3) OrderedVec[3]=tax4;
        else if (resvec[4]==4) OrderedVec[4]=tax4;
        else printf("There was a problem determining order in CountQuartetSitePatterns 4  - %d %d %d %d %d %d %d %d \n\n",tax1,tax2,tax3,tax4,resvec[1],resvec[2],resvec[3],resvec[4]);

  return sum;

}



/*** Function to count site patterns for a specified quartet ***/
/*** adding over all of its lineages.  Counts are returned   ***/
/*** in a 15-dimensional vector called pvec.                 ***/
/*** Ambiguity codes are accommodated by splitting counts    ***/
/*** over site patterns. If include_gaps = 0, sites with     ***/
/*** gaps in at least one taxa are excluded. If = 1, sites   ***/
/*** for which a particular quartet has no gap are included. ***/

void CountQuartetSitePatterns(int tax1, int tax2, int tax3, int tax4, double pvec[15]) {

    int i, j, k, l;
    int i2, j2, k2, l2;
    int ltax1, ltax2, ltax3, ltax4;
    double count_array[16][16], p[15];
    double count_sites=0.0, sumps=0.0, datasum=0.0;
    double site_total=0.0;

    //printf("Ordered taxa are %d %d %d %d\n",tax1,tax2,tax3,tax4);

    for (i=0; i<15; i++) pvec[i]=0.0; 
  
    for (i2=0; i2<seq_counter[tax1-1]; i2++) {
    	for (j2=0; j2<seq_counter[tax2-1]; j2++) {
            for (k2=0; k2<seq_counter[tax3-1]; k2++) {
                for (l2=0; l2<seq_counter[tax4-1]; l2++) {

			count_sites = 0.0;
			for (i=0; i<15; i++) p[i]=0.0; 
			for (i=0; i<16; i++) for (j=0; j<16; j++) count_array[i][j]=0.0;

			//if (verbose==1) printf("Sequences are %d %d %d %d\n, ",ppSp_assign[tax1-1][i2],ppSp_assign[tax2-1][j2],ppSp_assign[tax3-1][k2],ppSp_assign[tax4-1][l2]);

		   for (i=0; i<num_unique; i++) {

			int b1 = ppBase_unique[ppSp_assign[tax1-1][i2]][i];
    			int b2 = ppBase_unique[ppSp_assign[tax2-1][j2]][i];
    			int b3 = ppBase_unique[ppSp_assign[tax3-1][k2]][i];
    			int b4 = ppBase_unique[ppSp_assign[tax4-1][l2]][i];

			// Fast path: no ambiguity codes
    			if (b1 < 4 && b2 < 4 && b3 < 4 && b4 < 4) {
        			count_sites += site_counter[i];
				p[pattern_index[b1*4 + b2][b3*4 + b4]] += site_counter[i];
        			continue;
    			}

    			// Skip sites where any base is a gap/unknown (weight sums to 0)
    			// You can define a max valid code, e.g. 14 for N
    			if (b1 > 14 || b2 > 14 || b3 > 14 || b4 > 14) continue;

			// Slow path: ambiguity codes
			//BaseWeights w1 = get_base_weights(b1);
                        //BaseWeights w2 = get_base_weights(b2);
                        //BaseWeights w3 = get_base_weights(b3);
                        //BaseWeights w4 = get_base_weights(b4);
			double *w1 = base_weight_table[b1];
    			double *w2 = base_weight_table[b2];
    			double *w3 = base_weight_table[b3];
    			double *w4 = base_weight_table[b4];

			double w12[16], w34[16];
    			for (int a = 0; a < 4; a++)
        		  for (int b = 0; b < 4; b++)
            		    w12[a*4+b] = w1[a] * w2[b];
    			for (int c = 0; c < 4; c++)
        		  for (int d = 0; d < 4; d++)
            		    w34[c*4+d] = w3[c] * w4[d];

    			count_sites += site_counter[i];
    			double sc = (double)site_counter[i];
			for (int row = 0; row < 16; row++) {
    			 if (w12[row] == 0.0) continue;
    			  double contrib = sc * w12[row];
    			  for (int col = 0; col < 16; col++) {
			   p[pattern_index[row][col]] += contrib * w34[col];
    			  }
			 }

			}

    			//printf(" \n\n\t\tcount_sites is %f, ",count_sites);
    			// order is: 0 = xxxx; 1 = xxxy; 2 = xxyx; 3 = xyxx; 4 = yxxx
    			// 5 = xyxy; 6 = yxxy; 7 = xxyy; 8 = xyxz; 9 = xyzx; 10 = yxxz
    			// 11 = yxzx; 12 = xxyz; 13 = yzxx; 14 = xyzw

    			sumps = 0.0;
    			for (i=0; i<15; i++) sumps += p[i];
			//printf(" sumps is %f, ",sumps);
    			if (fabs((1.0/count_sites)*sumps - 1.0) > 0.005) {
        			printf("There was a problem counting site patterns ... exiting.");
        			exit(1);
    			}

    			/*if (verbose==1) { 
				printf("\tSite pattern counts are: ");
    				for (i=0; i<15; i++) printf("%f ",p[i]);
    				printf("\n"); 
			}*/
			

			for (i=0; i<15; i++) pvec[i] = pvec[i] + p[i]; 
    
			datasum = datasum + count_sites;
			//printf("\tcount_sites is %f, datasum is %f\n",count_sites,datasum);
 
		    }
		}
	     }
	}

	/*if (verbose==1) {
		printf("\n\n The overall site pattern counts are: ");
		for (i=0; i<15; i++) printf("%f ",pvec[i]);
		printf("\n");
	}*/
	sumps = 0.0;
        for (i=0; i<15; i++) sumps += pvec[i];
	//if (verbose==1) printf("The total number of site patterns is %f, datasum is %f\n",sumps,datasum);
        if (fabs((1.0/datasum)*sumps - 1.0) > 0.005) {
            printf("There was a problem counting overall site patterns ... exiting.");                                                        
            exit(1);
        }
}



double SymmetricQuartetLikelihood(int nn){

  int i, sump=0;
  double m = 4.0/3.0, t=2.0*theta; // t is 2*theta   
  double sbeta[10][1];
  double q[12],check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3, expt1, expt2, expt3;

  // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  //printf("Times are %f %f %f\n",t1,t2,t3);
  expt1 = exp(-2*m*t1);
  expt2 = exp(-2*m*t2);
  expt3 = exp(-2*m*t3);

  sbeta[1][0]=1.0;
  sbeta[2][0]=expt1;
  sbeta[3][0]=expt2;
  sbeta[4][0]=expt1*expt2;
  sbeta[5][0]=expt3;
  sbeta[6][0]=exp(-m*t1)*expt3;
  sbeta[7][0]=exp(-m*t2)*expt3;
  sbeta[8][0]=exp(-m*t1)*exp(-m*t2)*expt3;
  sbeta[9][0]=exp(2*t1/t)*exp(2*t2/t)*exp(-4*t3*(m+1/t));

  // order is: 1 = xxxx; 2 = xxxy = xxyx; 3 = xyxx = yxxx; 4 = xyxy = yxxy
  // 5 = xxyy; 6 = xxyz; 7 = yzxx; 8 = xyxz = yxxz = xyzx = yxzx; 9 = xyzw
  // to get site patterns probs, multiply transpose(smat)*sbeta
  // weighted version is below

  for (i=1; i<12; i++) q[i]=0.0;
    for (i=1; i<10; i++) {
            q[1] += 4*smat[i][1]*sbeta[i][0];
            q[2] += 12*smat[i][2]*sbeta[i][0];//24
            q[3] += 12*smat[i][3]*sbeta[i][0];//24
            q[4] += 12*smat[i][4]*sbeta[i][0];//24
            q[5] += 12*smat[i][5]*sbeta[i][0];
            q[6] += 24*smat[i][6]*sbeta[i][0];
            q[7] += 24*smat[i][7]*sbeta[i][0];
            q[8] += 24*smat[i][8]*sbeta[i][0];//96
            q[9] += 24*smat[i][9]*sbeta[i][0];      
  }

   // compute quartet likelihood
   quart_lik = StoreQuarts[nn]->spprobs[0]*log(q[1]) + StoreQuarts[nn]->spprobs[1]*log(q[2]) + StoreQuarts[nn]->spprobs[2]*log(q[2]) + StoreQuarts[nn]->spprobs[3]*log(q[3]) + StoreQuarts[nn]->spprobs[4]*log(q[3]) 
                    + StoreQuarts[nn]->spprobs[5]*log(q[4]) + StoreQuarts[nn]->spprobs[6]*log(q[4]) + StoreQuarts[nn]->spprobs[7]*log(q[5]) + StoreQuarts[nn]->spprobs[12]*log(q[6]) 
                    + StoreQuarts[nn]->spprobs[13]*log(q[7]) + StoreQuarts[nn]->spprobs[8]*log(q[8]) + StoreQuarts[nn]->spprobs[10]*log(q[8]) + StoreQuarts[nn]->spprobs[9]*log(q[8]) + StoreQuarts[nn]->spprobs[11]*log(q[8]) + StoreQuarts[nn]->spprobs[14]*log(q[9]);
   
  //if (verbose==1) printf("The likelihood of the symmetric quartet is %f\n\n",quart_lik);
  return(quart_lik);

}  


double AsymmetricQuartetLikelihood(int nn){                            

  int i, sump=0;
  double m = 4.0/3.0, t=2.0*theta; // t is 2*theta
  double abeta[11][1];
  double q[12],check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3, expt1, expt2, expt3;


  // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  expt1 = exp(-2.0*m*t1);
  expt2 = exp(-2.0*m*t2);
  expt3 = exp(-2.0*m*t3);

  abeta[1][0]=1.0; 
  abeta[2][0]=expt1;
  abeta[3][0]=expt2; 
  abeta[4][0]=exp(-m*t1)*expt2; 
  abeta[5][0]=expt3; 
  abeta[6][0]=exp(-m*t1)*expt3; 
  abeta[7][0]=expt1*expt3; 
  abeta[8][0]=exp(-m*t2)*expt3; 
  abeta[9][0]=exp(-m*t1)*exp(-m*t2)*expt3;
  abeta[10][0]=exp((2.0/t)*(t1-t2))*expt2*expt3; 

  // order is: xxxx; xxxy = xxyx; xyxx; yxxx; xxyy; xyxy = yxxy; xxyz
  // yzxx; xyxz = xyzx; yxxz = yxzx; xyzw
  // to get site patterns probs, multiply amat*abeta (not transpose)
  // weighted version is below

  for (i=1; i<12; i++) q[i]=0.0;
   for (i=1; i<11; i++) {
            q[1] += 4*amat[1][i]*abeta[i][0]; //4
            q[2] += 12*amat[2][i]*abeta[i][0]; //24
            q[3] += 12*amat[3][i]*abeta[i][0]; //12
            q[4] += 12*amat[4][i]*abeta[i][0]; //12
            q[5] += 12*amat[5][i]*abeta[i][0]; //12
            q[6] += 12*amat[6][i]*abeta[i][0]; //24
            q[7] += 24*amat[7][i]*abeta[i][0]; //24
            q[8] += 24*amat[8][i]*abeta[i][0]; //24
            q[9] += 24*amat[9][i]*abeta[i][0];  //48
            q[10] += 24*amat[10][i]*abeta[i][0]; //48
            q[11] += 24*amat[11][i]*abeta[i][0]; //24
   }

  // compute likelihood
  quart_lik = StoreQuarts[nn]->spprobs[0]*log(q[1]) + StoreQuarts[nn]->spprobs[1]*log(q[2]) + StoreQuarts[nn]->spprobs[2]*log(q[2]) + StoreQuarts[nn]->spprobs[3]*log(q[3]) + StoreQuarts[nn]->spprobs[4]*log(q[4])
                    + StoreQuarts[nn]->spprobs[7]*log(q[5]) + StoreQuarts[nn]->spprobs[5]*log(q[6]) + StoreQuarts[nn]->spprobs[6]*log(q[6]) + StoreQuarts[nn]->spprobs[12]*log(q[7]) + StoreQuarts[nn]->spprobs[13]*log(q[8])
                    + StoreQuarts[nn]->spprobs[8]*log(q[9]) + StoreQuarts[nn]->spprobs[9]*log(q[9]) + StoreQuarts[nn]->spprobs[10]*log(q[10]) + StoreQuarts[nn]->spprobs[11]*log(q[10]) + StoreQuarts[nn]->spprobs[14]*log(q[11]);

  //if (verbose==1) printf("The likelihood of the quartet is %f\n\n",quart_lik);
  return(quart_lik);

}


double GetCompLik() {

  int i, j, k, l, m;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          //printf("Quartet %d %d %d %d\n",i,j,k,l);
	  //printf("\t %llu %d %d\n",subset_to_index(i-1,j-1,k-1,l-1),(int)subset_to_index(i-1,j-1,k-1,l-1),qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]);
	  //printf("\t Have we looked at this quartet yet? %d\n",qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]); 
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

		GetQuartetTree(i,j,k,l,ovec);
  		//printf("ovec is: %d %d %d %d %d\n",ovec[0],ovec[1],ovec[2],ovec[3],ovec[4]);
                FindDupQuarts(ovec[0],duppvec);
            
   		// Fill the information for the current quartet into the struct StoreQuarts
   		for (m=0; m<15; m++) StoreQuarts[num_unique_quarts]->spprobs[m] = duppvec[m];
   		StoreQuarts[num_unique_quarts]->ncherries = ovec[0];
   		StoreQuarts[num_unique_quarts]->t1 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+3])];
   		if (ovec[0]==2) {
        		StoreQuarts[num_unique_quarts]->t2 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+2],ppTwoRowQuart[1][ntaxa+2])];
        		StoreQuarts[num_unique_quarts]->t3 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+2])];
   		}
   		else if (ovec[0]==1) {
        		StoreQuarts[num_unique_quarts]->t2 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+2])];
        		StoreQuarts[num_unique_quarts]->t3 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+1])];
   		}

		num_unique_quarts+=1;

	  }
        }
      }
    }
  }

  for (i=0; i<num_unique_quarts; i++) {

	//printf("%d ",i);
 	//for (j=0; j<15; j++) printf("%d ",StoreQuarts[i]->spprobs[j]);
	//printf("%f %f %f ",*StoreQuarts[i]->t1,*StoreQuarts[i]->t2,*StoreQuarts[i]->t3);
	//printf("%d\n",StoreQuarts[i]->ncherries);

	if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
	//printf("\t %f\n",complik);

  }


  //printf("The composite likelihood is %f; this required %d quartet comps\n",complik,num_unique_quarts);
  return(complik);

}

/* switch specifies which node in the quartet to use */
/* switch = 1, 2, 3 for t1, t2, t3                   */

void SymmetricQuartetLikelihoodDerivatives(int nn, int blswitch, double *symders){

  int i, sump=0;
  double m = 4.0/3.0, t=2.0*theta; // t is 2*theta   
  double sbeta[10][1], sdbeta[10][1], sddbeta[10][1];
  double q[12], dq[12], ddq[12], check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3, expt1, expt2, expt3, emt1, emt2;

  symders[0] = 0.0;  /* storage for first derivative */
  symders[1] = 0.0;  /* storage for second derivative */

  // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  printf("Times are %f %f %f\n",t1,t2,t3);
  expt1 = exp(-2.0*m*t1);
  expt2 = exp(-2.0*m*t2);
  expt3 = exp(-2.0*m*t3);
  emt1 = exp(-m*t1);
  emt2 = exp(-m*t2); 

  sbeta[1][0]=1.0;
  sbeta[2][0]=expt1;
  sbeta[3][0]=expt2;
  sbeta[4][0]=expt1*expt2;
  sbeta[5][0]=expt3;
  sbeta[6][0]=emt1*expt3;
  sbeta[7][0]=emt2*expt3;
  sbeta[8][0]=emt1*emt2*expt3;
  sbeta[9][0]=exp(2*t1/t)*exp(2*t2/t)*exp(-4*t3*(m+1/t));


  if (blswitch == 1) {

  	sdbeta[1][0]=0.0;
  	sdbeta[2][0]=-2.0*m*expt1;
  	sdbeta[3][0]=0.0;
  	sdbeta[4][0]=-2.0*m*expt1*expt2;
  	sdbeta[5][0]=0.0;
  	sdbeta[6][0]=-m*expt1*expt3;
  	sdbeta[7][0]=0.0;
  	sdbeta[8][0]=-m*emt1*emt2*expt3;
  	sdbeta[9][0]=(2.0/t)*exp(2.0*t1/t)*exp(2.0*t2/t)*exp(-4.0*t3*(m+1/t));

  	sddbeta[1][0]=0.0;
  	sddbeta[2][0]=4.0*(m/2.0)*expt1;
  	sddbeta[3][0]=0.0;
  	sddbeta[4][0]=4.0*(m/2)*expt1*expt2;
  	sddbeta[5][0]=0.0;
  	sddbeta[6][0]=m*m*emt1*expt3;
  	sddbeta[7][0]=0.0;
  	sddbeta[8][0]=m*m*emt1*emt2*expt3;
  	sddbeta[9][0]=(4.0/(t*t))*exp(2.0*t1/t)*exp(2.0*t2/t)*exp(-4.0*t3*(m+1/t));

  }

  else if (blswitch == 2) {

	sdbeta[1][0]=0.0;
  	sdbeta[2][0]=0.0;
  	sdbeta[3][0]=-2.0*m*expt2;
  	sdbeta[4][0]=-2.0*m*expt1*expt2;
  	sdbeta[5][0]=0.0;
  	sdbeta[6][0]=0.0;
  	sdbeta[7][0]=-m*emt2*expt3;
  	sdbeta[8][0]=-m*emt1*emt2*expt3;
  	sdbeta[9][0]=(2.0/t)*exp(2.0*t1/t)*exp(2.0*t2/t)*exp(-4.0*t3*(m+1/t));
  
  	sddbeta[1][0]=0.0;
  	sddbeta[2][0]=0.0;
  	sddbeta[3][0]=4.0*m*m*expt2;
  	sddbeta[4][0]=4.0*m*m*expt1*expt2;
  	sddbeta[5][0]=0.0;
  	sddbeta[6][0]=0.0;
  	sddbeta[7][0]=m*m*emt2*expt3;
  	sddbeta[8][0]=m*m*emt1*emt2*expt3;
  	sddbeta[9][0]=(4.0/(t*t))*exp(2.0*t1/t)*exp(2.0*t2/t)*exp(-4.0*t3*(m+1/t));

  }

  else if (blswitch == 3) {

	sdbeta[1][0]=0.0;
  	sdbeta[2][0]=0.0;
  	sdbeta[3][0]=0.0;
  	sdbeta[4][0]=0.0;
  	sdbeta[5][0]=-2.0*m*expt3;
  	sdbeta[6][0]=-2.0*m*emt1*expt3;
  	sdbeta[7][0]=-2.0*m*emt2*expt3;
  	sdbeta[8][0]=-2.0*m*emt1*emt2*expt3;
  	sdbeta[9][0]=(-4.0/t-4.0*m)*exp(2.0*t1/t)*exp(2.0*t2/t)*exp(-4.0*t3*(m+1.0/t)); 
  
  	sddbeta[1][0]=0.0;
  	sddbeta[2][0]=0.0;
  	sddbeta[3][0]=0.0;
  	sddbeta[4][0]=0.0;
  	sddbeta[5][0]=4.0*m*m*expt3;
  	sddbeta[6][0]=4.0*m*m*emt1*expt3;
  	sddbeta[7][0]=4.0*m*m*emt2*expt3;
  	sddbeta[8][0]=4.0*m*m*emt1*emt2*expt3;
  	sddbeta[9][0]=pow(-4.0/t-4.0*m,2)*exp(2.0*t1/t)*exp(2.0*t2/t)*exp(-4.0*t3*(m+1/t));

  }

  else { printf("Invalid argument to derivative function. Exiting.\n\n"); exit(1); }

  // order is: 1 = xxxx; 2 = xxxy = xxyx; 3 = xyxx = yxxx; 4 = xyxy = yxxy
  // 5 = xxyy; 6 = xxyz; 7 = yzxx; 8 = xyxz = yxxz = xyzx = yxzx; 9 = xyzw
  // to get site patterns probs, multiply transpose(smat)*sbeta
  // weighted version is below

  for (i=1; i<12; i++) { q[i]=0.0; dq[i]=0.0; ddq[i]=0.0; }
    for (i=1; i<10; i++) {

	    q[1] += 4*smat[i][1]*sbeta[i][0];
            q[2] += 12*smat[i][2]*sbeta[i][0];//24
            q[3] += 12*smat[i][3]*sbeta[i][0];//24
            q[4] += 12*smat[i][4]*sbeta[i][0];//24
            q[5] += 12*smat[i][5]*sbeta[i][0];
            q[6] += 24*smat[i][6]*sbeta[i][0];
            q[7] += 24*smat[i][7]*sbeta[i][0];
            q[8] += 24*smat[i][8]*sbeta[i][0];//96
            q[9] += 24*smat[i][9]*sbeta[i][0]; 

            dq[1] += 4*smat[i][1]*sdbeta[i][0];
            dq[2] += 12*smat[i][2]*sdbeta[i][0];//24
            dq[3] += 12*smat[i][3]*sdbeta[i][0];//24
            dq[4] += 12*smat[i][4]*sdbeta[i][0];//24
            dq[5] += 12*smat[i][5]*sdbeta[i][0];
            dq[6] += 24*smat[i][6]*sdbeta[i][0];
            dq[7] += 24*smat[i][7]*sdbeta[i][0];
            dq[8] += 24*smat[i][8]*sdbeta[i][0];//96
            dq[9] += 24*smat[i][9]*sdbeta[i][0];     

	    ddq[1] += 4*smat[i][1]*sddbeta[i][0];
            ddq[2] += 12*smat[i][2]*sddbeta[i][0];//24                                                                                 
            ddq[3] += 12*smat[i][3]*sddbeta[i][0];//24
            ddq[4] += 12*smat[i][4]*sddbeta[i][0];//24                                                                                 
            ddq[5] += 12*smat[i][5]*sddbeta[i][0];
            ddq[6] += 24*smat[i][6]*sddbeta[i][0];
            ddq[7] += 24*smat[i][7]*sddbeta[i][0];
            ddq[8] += 24*smat[i][8]*sddbeta[i][0];//96
            ddq[9] += 24*smat[i][9]*sddbeta[i][0]; 
  }

  for (i=1; i<10; i++) printf("%f %f %f\n",sdbeta[i][0],dq[i],ddq[i]);


   // compute first derivative of the quartet likelihood
   symders[0] = StoreQuarts[nn]->spprobs[0]*(1.0/q[1])*dq[1] + StoreQuarts[nn]->spprobs[1]*(1.0/q[2])*dq[2] + StoreQuarts[nn]->spprobs[2]*(1.0/q[2])*dq[2] + StoreQuarts[nn]->spprobs[3]*(1.0/q[3])*dq[3] + StoreQuarts[nn]->spprobs[4]*(1.0/q[3])*dq[3] 
                    + StoreQuarts[nn]->spprobs[5]*(1.0/q[4])*dq[4] + StoreQuarts[nn]->spprobs[6]*(1.0/q[4])*dq[4] + StoreQuarts[nn]->spprobs[7]*(1.0/q[5])*dq[5] + StoreQuarts[nn]->spprobs[12]*(1.0/q[6])*dq[6] 
                    + StoreQuarts[nn]->spprobs[13]*(1.0/q[7])*dq[7] + StoreQuarts[nn]->spprobs[8]*(1.0/q[8])*dq[8] + StoreQuarts[nn]->spprobs[10]*(1.0/q[8])*dq[8] + StoreQuarts[nn]->spprobs[9]*(1.0/q[8])*dq[8] + StoreQuarts[nn]->spprobs[11]*(1.0/q[8])*dq[8] 
		    + StoreQuarts[nn]->spprobs[14]*(1.0/q[9])*dq[9];
   
   // compute second derivative of the quartet likelihood 
   symders[1] = StoreQuarts[nn]->spprobs[0]*((1.0/q[1])*ddq[1]-(1.0/(q[1]*q[1]))*(dq[1]*dq[1])) 
		+ StoreQuarts[nn]->spprobs[1]*((1.0/q[2])*ddq[2]-(1.0/(q[2]*q[2]))*(dq[2]*dq[2])) 
		+ StoreQuarts[nn]->spprobs[2]*((1.0/q[2])*ddq[2]-(1.0/(q[2]*q[2]))*(dq[2]*dq[2])) 
		+ StoreQuarts[nn]->spprobs[3]*((1.0/q[3])*ddq[3]-(1.0/(q[3]*q[3]))*(dq[3]*dq[3])) 
		+ StoreQuarts[nn]->spprobs[4]*((1.0/q[3])*ddq[3]-(1.0/(q[3]*q[3]))*(dq[3]*dq[3]))
                + StoreQuarts[nn]->spprobs[5]*((1.0/q[4])*ddq[4]-(1.0/(q[4]*q[4]))*(dq[4]*dq[4])) 
		+ StoreQuarts[nn]->spprobs[6]*((1.0/q[4])*ddq[4]-(1.0/(q[4]*q[4]))*(dq[4]*dq[4])) 
		+ StoreQuarts[nn]->spprobs[7]*((1.0/q[5])*ddq[5]-(1.0/(q[5]*q[5]))*(dq[5]*dq[5])) 
		+ StoreQuarts[nn]->spprobs[12]*((1.0/q[6])*ddq[6]-(1.0/(q[6]*q[6]))*(dq[6]*dq[6]))
                + StoreQuarts[nn]->spprobs[13]*((1.0/q[7])*ddq[7]-(1.0/(q[7]*q[7]))*(dq[7]*dq[7])) 
		+ StoreQuarts[nn]->spprobs[8]*((1.0/q[8])*ddq[8]-(1.0/(q[8]*q[8]))*(dq[8]*dq[8])) 
		+ StoreQuarts[nn]->spprobs[10]*((1.0/q[8])*ddq[8]-(1.0/(q[8]*q[8]))*(dq[8]*dq[8]))
		+ StoreQuarts[nn]->spprobs[9]*((1.0/q[8])*ddq[8]-(1.0/(q[8]*q[8]))*(dq[8]*dq[8]))
		+ StoreQuarts[nn]->spprobs[11]*((1.0/q[8])*ddq[8]-(1.0/(q[8]*q[8]))*(dq[8]*dq[8]))
                + StoreQuarts[nn]->spprobs[14]*((1.0/q[9])*ddq[9]-(1.0/(q[9]*q[9]))*(dq[9]*dq[9]));

   
  if (verbose==1) printf("The derivatives of the likelihood of the symmetric quartet are %f %f\n\n",symders[0],symders[1]);

}  


void AsymmetricQuartetLikelihoodDerivatives(int nn, int blswitch, double *asymders){                            

  int i, sump=0;
  double m = 4.0/3.0, t=2.0*theta; // t is 2*theta
  double abeta[11][1], adbeta[11][1], addbeta[11][0];
  double q[12], dq[12], ddq[12], check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3, expt1, expt2, expt3, emt1, emt2;

  // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  expt1 = exp(-2.0*m*t1);
  expt2 = exp(-2.0*m*t2);
  expt3 = exp(-2.0*m*t3);
  emt1 = exp(-m*t1);
  emt2 = exp(-m*t2);

  abeta[1][0]=1.0; 
  abeta[2][0]=expt1;
  abeta[3][0]=expt2; 
  abeta[4][0]=emt1*expt2; 
  abeta[5][0]=expt3; 
  abeta[6][0]=emt1*expt3; 
  abeta[7][0]=expt1*expt3; 
  abeta[8][0]=emt2*expt3; 
  abeta[9][0]=emt1*emt2*expt3;
  abeta[10][0]=exp((2.0/t)*(t1-t2))*expt2*expt3; 

  if (blswitch == 1) {

	adbeta[1][0]=0.0;
  	adbeta[2][0]=-2*m*expt1;
 	adbeta[3][0]=0.0;
  	adbeta[4][0]=-m*emt1*expt2;
  	adbeta[5][0]=0.0;
  	adbeta[6][0]=-m*emt1*expt3;
  	adbeta[7][0]=-2*m*expt1*expt3;
  	adbeta[8][0]=0.0;
  	adbeta[9][0]=-m*emt1*emt2*expt3;
  	adbeta[10][0]=(2.0/t)*exp((2.0/t)*(t1-t2))*expt2*expt3;

	addbeta[1][0]=0.0;
  	addbeta[2][0]=4*m*m*expt1;
  	addbeta[3][0]=0.0;
  	addbeta[4][0]=m*m*emt1*expt2;
  	addbeta[5][0]=0.0;
  	addbeta[6][0]=m*m*emt1*expt3;
  	addbeta[7][0]=4*m*m*expt1*expt3;
  	addbeta[8][0]=0.0;
  	addbeta[9][0]=m*m*emt1*emt2*expt3;
  	addbeta[10][0]=(4.0/(t*t))*exp((2.0/t)*(t1-t2))*expt2*expt3;

  }
 
  else if (blswitch == 2) {

        adbeta[1][0]=0.0;
  	adbeta[2][0]=0.0; 
  	adbeta[3][0]=-2*m*expt2;
  	adbeta[4][0]=-2*m*emt1*expt2;
  	adbeta[5][0]=0.0;
  	adbeta[6][0]=0.0;
  	adbeta[7][0]=0.0;
  	adbeta[8][0]=-m*emt2*expt3;
  	adbeta[9][0]=-m*emt1*emt2*expt3;
  	adbeta[10][0]=(-2.0/t-2.0*m)*exp((2.0/t)*(t1-t2))*expt2*expt3;
  
        addbeta[1][0]=0.0;
  	addbeta[2][0]=0.0;
  	addbeta[3][0]=4*m*m*expt2;
  	addbeta[4][0]=4*m*m*emt1*expt2;
  	addbeta[5][0]=0.0;
  	addbeta[6][0]=0.0;
  	addbeta[7][0]=0.0;
  	addbeta[8][0]=m*m*emt2*expt3;
  	addbeta[9][0]=m*m*emt1*emt2*expt3;
  	addbeta[10][0]=((2.0/t+2.0*m)*(2.0/t+2.0*m))*exp((2.0/t)*(t1-t2))*expt2*expt3;

  }

  else if (blswitch == 3) {

        adbeta[1][0]=0.0;
  	adbeta[2][0]=0.0; 
  	adbeta[3][0]=0.0;
  	adbeta[4][0]=0.0;
  	adbeta[5][0]=-2*m*expt3;
  	adbeta[6][0]=-2*m*emt1*expt3;
  	adbeta[7][0]=-2*m*expt1*expt3;
  	adbeta[8][0]=-2*m*emt2*expt3;
  	adbeta[9][0]=-2*m*emt1*emt2*expt3;
  	adbeta[10][0]=(-2*m)*exp((2.0/t)*(t1-t2))*expt2*expt3;
  
        addbeta[1][0]=0.0;
  	addbeta[2][0]=0.0;  
  	addbeta[3][0]=0.0;
  	addbeta[4][0]=0.0;
  	addbeta[5][0]=4*m*m*expt3;
  	addbeta[6][0]=4*m*m*emt1*expt3;
  	addbeta[7][0]=4*m*m*expt1*expt3;
  	addbeta[8][0]=4*m*m*emt2*expt3;
  	addbeta[9][0]=4*m*m*emt1*emt2*expt3;
  	addbeta[10][0]=(4*m*m)*exp((2.0/t)*(t1-t2))*expt2*expt3;

  }

  else { printf("Invalid argument to derivative function. Exiting.\n\n"); exit(1); }

  // order is: xxxx; xxxy = xxyx; xyxx; yxxx; xxyy; xyxy = yxxy; xxyz
  // yzxx; xyxz = xyzx; yxxz = yxzx; xyzw
  // to get site patterns probs, multiply amat*abeta (not transpose)
  // weighted version is below

  for (i=1; i<12; i++) { q[i]=0.0; dq[i]=0.0; ddq[i]=0.0; }
   for (i=1; i<11; i++) {
            q[1] += 4*amat[1][i]*abeta[i][0]; //4
            q[2] += 12*amat[2][i]*abeta[i][0]; //24
            q[3] += 12*amat[3][i]*abeta[i][0]; //12
            q[4] += 12*amat[4][i]*abeta[i][0]; //12
            q[5] += 12*amat[5][i]*abeta[i][0]; //12
            q[6] += 12*amat[6][i]*abeta[i][0]; //24
            q[7] += 24*amat[7][i]*abeta[i][0]; //24
            q[8] += 24*amat[8][i]*abeta[i][0]; //24
            q[9] += 24*amat[9][i]*abeta[i][0];  //48
            q[10] += 24*amat[10][i]*abeta[i][0]; //48
            q[11] += 24*amat[11][i]*abeta[i][0]; //24

	    dq[1] += 4*amat[1][i]*adbeta[i][0]; //4
            dq[2] += 12*amat[2][i]*adbeta[i][0]; //24
            dq[3] += 12*amat[3][i]*adbeta[i][0]; //12
            dq[4] += 12*amat[4][i]*adbeta[i][0]; //12
            dq[5] += 12*amat[5][i]*adbeta[i][0]; //12
            dq[6] += 12*amat[6][i]*adbeta[i][0]; //24
            dq[7] += 24*amat[7][i]*adbeta[i][0]; //24
            dq[8] += 24*amat[8][i]*adbeta[i][0]; //24                                                                                                   
            dq[9] += 24*amat[9][i]*adbeta[i][0];  //48                                                                                                  
            dq[10] += 24*amat[10][i]*adbeta[i][0]; //48                                                                                                 
            dq[11] += 24*amat[11][i]*adbeta[i][0]; //24   

	    ddq[1] += 4*amat[1][i]*addbeta[i][0]; //4
            ddq[2] += 12*amat[2][i]*addbeta[i][0]; //24
            ddq[3] += 12*amat[3][i]*addbeta[i][0]; //12
            ddq[4] += 12*amat[4][i]*addbeta[i][0]; //12
            ddq[5] += 12*amat[5][i]*addbeta[i][0]; //12
            ddq[6] += 12*amat[6][i]*addbeta[i][0]; //24
            ddq[7] += 24*amat[7][i]*addbeta[i][0]; //24
            ddq[8] += 24*amat[8][i]*addbeta[i][0]; //24                                                                                                   
            ddq[9] += 24*amat[9][i]*addbeta[i][0];  //48                                                                                                  
            ddq[10] += 24*amat[10][i]*addbeta[i][0]; //48                                                                                                 
            ddq[11] += 24*amat[11][i]*addbeta[i][0]; //24   

   }

  for (i=1; i<12; i++) printf("%f %f %f %f\n",abeta[i][0],adbeta[i][0],dq[i],ddq[i]);

  // compute first derviative of the quartet likelihood
  asymders[0] = StoreQuarts[nn]->spprobs[0]*(1.0/q[1])*dq[1] 
		+ StoreQuarts[nn]->spprobs[1]*(1.0/q[2])*dq[2] 
		+ StoreQuarts[nn]->spprobs[2]*(1.0/q[2])*dq[2] 
		+ StoreQuarts[nn]->spprobs[3]*(1.0/q[3])*dq[3] 
		+ StoreQuarts[nn]->spprobs[4]*(1.0/q[4])*dq[4]
                + StoreQuarts[nn]->spprobs[7]*(1.0/q[5])*dq[5] 
		+ StoreQuarts[nn]->spprobs[5]*(1.0/q[6])*dq[6] 
		+ StoreQuarts[nn]->spprobs[6]*(1.0/q[6])*dq[6] 
		+ StoreQuarts[nn]->spprobs[12]*(1.0/q[7])*dq[7] 
		+ StoreQuarts[nn]->spprobs[13]*(1.0/q[8])*dq[8]
                + StoreQuarts[nn]->spprobs[8]*(1.0/q[9])*dq[9] 
		+ StoreQuarts[nn]->spprobs[9]*(1.0/q[9])*dq[9] 
		+ StoreQuarts[nn]->spprobs[10]*(1.0/q[10])*dq[10] 
		+ StoreQuarts[nn]->spprobs[11]*(1.0/q[10])*dq[10]
		+ StoreQuarts[nn]->spprobs[14]*(1.0/q[11])*dq[11];

 // compute second derviative of the quartet likelihood
  asymders[1] = StoreQuarts[nn]->spprobs[0]*((1.0/q[1])*ddq[1]-(1.0/(q[1]*q[1]))*(dq[1]*dq[1]))
                + StoreQuarts[nn]->spprobs[1]*((1.0/q[2])*ddq[2]-(1.0/(q[2]*q[2]))*(dq[2]*dq[2]))  
                + StoreQuarts[nn]->spprobs[2]*((1.0/q[2])*ddq[2]-(1.0/(q[2]*q[2]))*(dq[2]*dq[2]))  
                + StoreQuarts[nn]->spprobs[3]*((1.0/q[3])*ddq[3]-(1.0/(q[3]*q[3]))*(dq[3]*dq[3]))  
                + StoreQuarts[nn]->spprobs[4]*((1.0/q[4])*ddq[4]-(1.0/(q[4]*q[4]))*(dq[4]*dq[4]))     
                + StoreQuarts[nn]->spprobs[7]*((1.0/q[5])*ddq[5]-(1.0/(q[5]*q[5]))*(dq[5]*dq[5])) 
                + StoreQuarts[nn]->spprobs[5]*((1.0/q[6])*ddq[6]-(1.0/(q[6]*q[6]))*(dq[6]*dq[6])) 
                + StoreQuarts[nn]->spprobs[6]*((1.0/q[6])*ddq[6]-(1.0/(q[6]*q[6]))*(dq[6]*dq[6])) 
                + StoreQuarts[nn]->spprobs[12]*((1.0/q[7])*ddq[7]-(1.0/(q[7]*q[7]))*(dq[7]*dq[7])) 
                + StoreQuarts[nn]->spprobs[13]*((1.0/q[8])*ddq[8]-(1.0/(q[8]*q[8]))*(dq[8]*dq[8])) 
                + StoreQuarts[nn]->spprobs[8]*((1.0/q[9])*ddq[9]-(1.0/(q[9]*q[9]))*(dq[9]*dq[9])) 
                + StoreQuarts[nn]->spprobs[9]*((1.0/q[9])*ddq[9]-(1.0/(q[9]*q[9]))*(dq[9]*dq[9])) 
                + StoreQuarts[nn]->spprobs[10]*((1.0/q[10])*ddq[10]-(1.0/(q[10]*q[10]))*(dq[10]*dq[10])) 
                + StoreQuarts[nn]->spprobs[11]*((1.0/q[10])*ddq[10]-(1.0/(q[10]*q[10]))*(dq[10]*dq[10]))
                + StoreQuarts[nn]->spprobs[14]*((1.0/q[11])*ddq[11]-(1.0/(q[11]*q[11]))*(dq[11]*dq[11])); 


  if (verbose==1) printf("The derivatives of the likelihood of the asymmetric quartet are %f %f\n\n",asymders[0],asymders[1]);

}


/*** Functions to compute first and second derivatives  ***/
/*** of the composite log likelihood for one node time. ***/
/*** The argument optnode is the integer index of the   ***/
/*** node whose time is being optimized.                ***/

void GetCompLikDerivatives(int optnode) {

  int i, j, k, l, m;
  int ovec[5], bl_ind[nquarts+1];
  int timemap1, timemap2, timemap3;
  double duppvec[15];
  double complik=0.0, bl_ders[2], der1=0.0, der2=0.0;

  bl_ders[0] = bl_ders[1] = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) { qvec[i]=0; bl_ind[i]=0;}

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          printf("Quartet %d %d %d %d\n",i,j,k,l);
          //printf("\t %llu %d %d\n",subset_to_index(i-1,j-1,k-1,l-1),(int)subset_to_index(i-1,j-1,k-1,l-1),qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]);
          //printf("\t Have we looked at this quartet yet? %d\n",qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]); 
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
                //printf("ovec is: %d %d %d %d %d\n",ovec[0],ovec[1],ovec[2],ovec[3],ovec[4]);
                FindDupQuarts(ovec[0],duppvec);
            
                // Fill the information for the current quartet into the struct StoreQuarts
                for (m=0; m<15; m++) StoreQuarts[num_unique_quarts]->spprobs[m] = duppvec[m];
                StoreQuarts[num_unique_quarts]->ncherries = ovec[0];

		timemap1 = FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+3]);
                StoreQuarts[num_unique_quarts]->t1 = &TimeVec[timemap1];
		if (optnode == timemap1) bl_ind[num_unique_quarts] = 1;

                if (ovec[0]==2) {
			timemap2 = FindMRCA(ppTwoRowQuart[0][ntaxa+2],ppTwoRowQuart[1][ntaxa+2]);
                	timemap3 = FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+2]);
                        StoreQuarts[num_unique_quarts]->t2 = &TimeVec[timemap2];
                        StoreQuarts[num_unique_quarts]->t3 = &TimeVec[timemap3];
			if (optnode == timemap2) bl_ind[num_unique_quarts] = 2;
			else if (optnode == timemap3) bl_ind[num_unique_quarts] = 3;
                }
                else if (ovec[0]==1) {
			timemap2 = FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+2]);
			timemap3 = FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+1]);
                        StoreQuarts[num_unique_quarts]->t2 = &TimeVec[timemap2];
                        StoreQuarts[num_unique_quarts]->t3 = &TimeVec[timemap3];
			if (optnode == timemap2) bl_ind[num_unique_quarts] = 2;
                        else if (optnode == timemap3) bl_ind[num_unique_quarts] = 3;
                }

                num_unique_quarts+=1;

          }
        }
      }
    }
  }

  for (i=0; i<num_unique_quarts; i++) {

        printf("%d ",i);
        //for (j=0; j<15; j++) printf("%d ",StoreQuarts[i]->spprobs[j]);
        printf("%f %f %f ",*StoreQuarts[i]->t1,*StoreQuarts[i]->t2,*StoreQuarts[i]->t3);
	printf("Opt node is %d, bl_ind is %d \n",optnode,bl_ind[i]);
        //printf("%d\n",StoreQuarts[i]->ncherries);

        if (StoreQuarts[i]->ncherries==2 && bl_ind[i]!=0) {
		SymmetricQuartetLikelihoodDerivatives(i,bl_ind[i],bl_ders);
		der1+=bl_ders[0];
		der2+=bl_ders[1];
	}
        else if (StoreQuarts[i]->ncherries==1 && bl_ind[i]!=0) {
		AsymmetricQuartetLikelihoodDerivatives(i,bl_ind[i],bl_ders);
		printf("blders: %f %f\n",bl_ders[0],bl_ders[1]);
                der1+=bl_ders[0];
                der2+=bl_ders[1];
	}

  }


  printf("The derviatives of the composite likelihood are %f  %f\n",der1,der2);

}





/*******************************************************/
/*** Functions for multi-allelic SNPs                ***/
/*******************************************************/


/*** Function to count site patterns for a specified quartet ***/
/*** Note that if include-gaps = 0, only - are excluded.     ***/
/*** ? still appear in ppBase_unique. The code below removes ***/
/*** sites with ? or any other ambiguity codes before        ***/
/*** counting site patterns.                                 ***/

double SymmetricQuartetLikelihood_msnp(int nn){

  int i, sump=0;
  int datasum, count_noambigs;
  double m = 4.0/3.0, t=2.0*theta; // t is 2*theta
  double sbeta[10][1];
  double q[12], check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3;

  // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  //printf("Times are %f %f %f\n",t1,t2,t3);

  sbeta[1][0]=1.0;
  sbeta[2][0]=exp(-2*m*t1);
  sbeta[3][0]=exp(-2*m*t2);
  sbeta[4][0]=exp(-2*m*t1)*exp(-2*m*t2);
  sbeta[5][0]=exp(-2*m*t3);
  sbeta[6][0]=exp(-m*t1)*exp(-2*m*t3);
  sbeta[7][0]=exp(-m*t2)*exp(-2*m*t3);
  sbeta[8][0]=exp(-m*t1)*exp(-m*t2)*exp(-2*m*t3);
  sbeta[9][0]=exp(2*t1/t)*exp(2*t2/t)*exp(-4*t3*(m+1/t));  

  // order is: 1 = xxxx; 2 = xxxy = xxyx; 3 = xyxx = yxxx; 4 = xyxy = yxxy
  // 5 = xxyy; 6 = xxyz; 7 = yzxx; 8 = xyxz = yxxz = xyzx = yxzx; 9 = xyzw
  // to get site patterns probs, multiply transpose(smat)*sbeta
  // weighted version is below

  for (i=1; i<12; i++) q[i]=0.0;
    for (i=1; i<10; i++) {
            q[1] += 4*smat[i][1]*sbeta[i][0];
            q[2] += 12*smat[i][2]*sbeta[i][0];
            q[3] += 12*smat[i][3]*sbeta[i][0];
            q[4] += 12*smat[i][4]*sbeta[i][0];
            q[5] += 12*smat[i][5]*sbeta[i][0];
            q[6] += 24*smat[i][6]*sbeta[i][0];
            q[7] += 24*smat[i][7]*sbeta[i][0];
            q[8] += 24*smat[i][8]*sbeta[i][0];
            q[9] += 24*smat[i][9]*sbeta[i][0];      
  }

    
  // compute quartet likelihood
  for (i=0; i<15; i++) sump += StoreQuarts[nn]->spprobs[i];
  prsum = 1.0 - q[1];
  datasum = sump - StoreQuarts[nn]->spprobs[0]; 
  quart_lik = StoreQuarts[nn]->spprobs[1]*log(q[2]) + StoreQuarts[nn]->spprobs[2]*log(q[2]) + StoreQuarts[nn]->spprobs[3]*log(q[3]) + StoreQuarts[nn]->spprobs[4]*log(q[3])
                    + StoreQuarts[nn]->spprobs[5]*log(q[4]) + StoreQuarts[nn]->spprobs[6]*log(q[4]) + StoreQuarts[nn]->spprobs[7]*log(q[5]) + StoreQuarts[nn]->spprobs[12]*log(q[6])
                    + StoreQuarts[nn]->spprobs[13]*log(q[7]) + StoreQuarts[nn]->spprobs[8]*log(q[8]) + StoreQuarts[nn]->spprobs[10]*log(q[8]) + StoreQuarts[nn]->spprobs[9]*log(q[8]) + StoreQuarts[nn]->spprobs[11]*log(q[8]) + StoreQuarts[nn]->spprobs[14]*log(q[9])
                    - datasum*log(prsum);
 
  //if (verbose==1) printf("The likelihood of the symmetric quartet is %f\n\n",quart_lik);
  return(quart_lik);

}
 

double AsymmetricQuartetLikelihood_msnp(int nn){

  int i, sump=0;
  int datasum, count_noambigs;
  double m = 4.0/3.0, t=2.0*theta; // t is 2*theta
  double abeta[12][1];
  double q[12], check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3;

  // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  //printf("Times are %f %f %f\n",t1,t2,t3);

  abeta[1][0]=1.0; 
  abeta[2][0]=exp(-2.0*m*t1);
  abeta[3][0]=exp(-2.0*m*t2); 
  abeta[4][0]=exp(-m*t1)*exp(-2.0*m*t2); 
  abeta[5][0]=exp(-2.0*m*t3); 
  abeta[6][0]=exp(-m*t1)*exp(-2.0*m*t3); 
  abeta[7][0]=exp(-2.0*m*t1)*exp(-2.0*m*t3); 
  abeta[8][0]=exp(-m*t2)*exp(-2.0*m*t3); 
  abeta[9][0]=exp(-m*t1)*exp(-m*t2)*exp(-2.0*m*t3);
  abeta[10][0]=exp((2.0/t)*(t1-t2))*exp(-2.0*m*(t2+t3)); 

  // order is: xxxx; xxxy = xxyx; xyxx; yxxx; xxyy; xyxy = yxxy; xxyz
  // yzxx; xyxz = xyzx; yxxz = yxzx; xyzw
  // to get site patterns probs, multiply amat*abeta (not transpose)
  // weighted version is below

  for (i=1; i<12; i++) q[i]=0.0;
     for (i=1; i<11; i++) {
            q[1] += 4*amat[1][i]*abeta[i][0]; //4
            q[2] += 12*amat[2][i]*abeta[i][0]; //24
            q[3] += 12*amat[3][i]*abeta[i][0]; //12
            q[4] += 12*amat[4][i]*abeta[i][0]; //12
            q[5] += 12*amat[5][i]*abeta[i][0]; //12
            q[6] += 12*amat[6][i]*abeta[i][0]; //24
            q[7] += 24*amat[7][i]*abeta[i][0]; //24
            q[8] += 24*amat[8][i]*abeta[i][0]; //24
            q[9] += 24*amat[9][i]*abeta[i][0];  //48
            q[10] += 24*amat[10][i]*abeta[i][0]; //48
            q[11] += 24*amat[11][i]*abeta[i][0]; //24
   }

   // compute likelihood
   for (i=0; i<15; i++) sump += StoreQuarts[nn]->spprobs[i];
   prsum = 1.0 - q[1];
   datasum = sump - StoreQuarts[nn]->spprobs[0];
   quart_lik = StoreQuarts[nn]->spprobs[1]*log(q[2]) + StoreQuarts[nn]->spprobs[2]*log(q[2]) + StoreQuarts[nn]->spprobs[3]*log(q[3]) + StoreQuarts[nn]->spprobs[4]*log(q[4])
                    + StoreQuarts[nn]->spprobs[7]*log(q[5]) + StoreQuarts[nn]->spprobs[5]*log(q[6]) + StoreQuarts[nn]->spprobs[6]*log(q[6]) + StoreQuarts[nn]->spprobs[12]*log(q[7]) + StoreQuarts[nn]->spprobs[13]*log(q[8])
                    + StoreQuarts[nn]->spprobs[8]*log(q[9]) + StoreQuarts[nn]->spprobs[9]*log(q[9]) + StoreQuarts[nn]->spprobs[10]*log(q[10]) + StoreQuarts[nn]->spprobs[11]*log(q[10]) + StoreQuarts[nn]->spprobs[14]*log(q[11])
                    - datasum*log(prsum);
        
    //if (verbose==1) printf("The likelihood of the quartet is %f\n\n",quart_lik);
    return(quart_lik);

}


double GetCompLik_msnp() {

  int i, j, k, l, m=0;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;

  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          //printf("Quartet %d %d %d %d\n",i,j,k,l);
          //printf("\t %llu %d %d\n",subset_to_index(i-1,j-1,k-1,l-1),(int)subset_to_index(i-1,j-1,k-1,l-1),qvec[(int)subset_to_index$
          //printf("\t Have we looked at this quartet yet? %d\n",qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]);
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {
 
		GetQuartetTree(i,j,k,l,ovec);
                //printf("ovec is: %d %d %d %d %d\n",ovec[0],ovec[1],ovec[2],ovec[3],ovec[4]);
                FindDupQuarts(ovec[0],duppvec);
            
                // Fill the information for the current quartet into the struct StoreQuarts
                for (m=0; m<15; m++) StoreQuarts[num_unique_quarts]->spprobs[m] = duppvec[m];
                StoreQuarts[num_unique_quarts]->ncherries = ovec[0];
                StoreQuarts[num_unique_quarts]->t1 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+3])];
                if (ovec[0]==2) {
                        StoreQuarts[num_unique_quarts]->t2 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+2],ppTwoRowQuart[1][ntaxa+2])];
                        StoreQuarts[num_unique_quarts]->t3 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+2])];
                }
                else if (ovec[0]==1) {
                        StoreQuarts[num_unique_quarts]->t2 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+2])];
                        StoreQuarts[num_unique_quarts]->t3 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+1])];
                }

                num_unique_quarts+=1;

          }
        }
      }
    }
  }

  for (i=0; i<num_unique_quarts; i++) {

        //printf("%d ",i);
        //for (j=0; j<15; j++) printf("%d ",StoreQuarts[i]->spprobs[j]);
        //printf("%f %f %f ",*StoreQuarts[i]->t1,*StoreQuarts[i]->t2,*StoreQuarts[i]->t3);
        //printf("%d\n",StoreQuarts[i]->ncherries);

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_msnp(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood_msnp(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        //printf("\t %f\n",complik);

  }

  //printf("The composite likelihood is %f; this required %d quartet comps\n",complik,num_unique_quarts);
  return(complik);
 
}
  



/*******************************************/
/*** Functions for rate-variation models ***/
/*******************************************/


/*** Compute rate variation model parameters for a given alpha ***/

void GetRateParams() {

   int i, errcode;
   double gg;
   double cutoffs[ncat+1];

   cutoffs[0] = 0.0;
   cutoffs[ncat] = 1000000.0;
   gg = lgamma(ratepar);

   /* first find the cutoff points */
   for (i=1; i<ncat; i++) cutoffs[i] = ppchi2((double)i/ncat,2*ratepar,gg,&errcode)/(2.0*ratepar);
   //for (i=0; i<ncat+1; i++) printf("Category %d has cutoff %f\n",i,cutoffs[i]);

   /* now find mean rates */
   for (i=1; i<ncat+1; i++) rvals[i] = ncat*(gammad(ratepar*cutoffs[i],ratepar+1,&errcode)-gammad(ratepar*cutoffs[i-1],ratepar+1.0,&errcode));

   //printf("means of intervals are:\n");
   //for (i=1; i<ncat+1; i++) printf("%f\n",rvals[i]);

}




/*** Function to compute quartet likelihood for a specified quartet ***/
/*** Note that if include-gaps = 0, only - are excluded.            ***/
/*** ? still appear in ppBase_unique. The code below removes        ***/
/*** sites with ? or any other ambiguity codes before               ***/
/*** counting site patterns.                                        ***/

double SymmetricQuartetLikelihood_ratevar(int nn){

  int i, k, sump=0;
  double m = 4.0/3.0, t=2.0*theta; // t is 2*theta
  double sbeta[10][1];
  double q[12], check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3;

  GetRateParams();

   // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  //printf("Times are %f %f %f\n",t1,t2,t3);

  for (i=1; i<12; i++) q[i]=0.0;
  for (k=0; k<ncat; k++) {

		//printf("for category %d, rate is %f\n",k,rvals[k+1]);

        	sbeta[1][0]=1.0;
        	sbeta[2][0]=exp(-2*m*rvals[k+1]*t1);
        	sbeta[3][0]=exp(-2*m*rvals[k+1]*t2);
        	sbeta[4][0]=exp(-2*m*rvals[k+1]*t1)*exp(-2*m*rvals[k+1]*t2);
        	sbeta[5][0]=exp(-2*m*rvals[k+1]*t3);
        	sbeta[6][0]=exp(-m*rvals[k+1]*t1)*exp(-2*m*rvals[k+1]*t3);
        	sbeta[7][0]=exp(-m*rvals[k+1]*t2)*exp(-2*m*rvals[k+1]*t3);
        	sbeta[8][0]=exp(-m*rvals[k+1]*t1)*exp(-m*rvals[k+1]*t2)*exp(-2*m*rvals[k+1]*t3);
        	sbeta[9][0]=exp(2*rvals[k+1]*t1/t)*exp(2*rvals[k+1]*t2/t)*exp(-4*rvals[k+1]*t3*(m+1/t));

        	// order is: 1 = xxxx; 2 = xxxy = xxyx; 3 = xyxx = yxxx; 4 = xyxy = yxxy
        	// 5 = xxyy; 6 = xxyz; 7 = yzxx; 8 = xyxz = yxxz = xyzx = yxzx; 9 = xyzw
        	// to get site patterns probs, multiply transpose(smat)*sbeta
        	// weighted version is below

        	for (i=1; i<10; i++) {
            		q[1] += (1.0/ncat)*4*smat[i][1]*sbeta[i][0];
            		q[2] += (1.0/ncat)*12*smat[i][2]*sbeta[i][0];
            		q[3] += (1.0/ncat)*12*smat[i][3]*sbeta[i][0];
            		q[4] += (1.0/ncat)*12*smat[i][4]*sbeta[i][0];
            		q[5] += (1.0/ncat)*12*smat[i][5]*sbeta[i][0];
            		q[6] += (1.0/ncat)*24*smat[i][6]*sbeta[i][0];
            		q[7] += (1.0/ncat)*24*smat[i][7]*sbeta[i][0];
            		q[8] += (1.0/ncat)*24*smat[i][8]*sbeta[i][0];
            		q[9] += (1.0/ncat)*24*smat[i][9]*sbeta[i][0];      
        	}

   }

   // compute quartet likelihood
   quart_lik = StoreQuarts[nn]->spprobs[0]*log(q[1]) + StoreQuarts[nn]->spprobs[1]*log(q[2]) + StoreQuarts[nn]->spprobs[2]*log(q[2]) + StoreQuarts[nn]->spprobs[3]*log(q[3]) + StoreQuarts[nn]->spprobs[4]*log(q[3])
                    + StoreQuarts[nn]->spprobs[5]*log(q[4]) + StoreQuarts[nn]->spprobs[6]*log(q[4]) + StoreQuarts[nn]->spprobs[7]*log(q[5]) + StoreQuarts[nn]->spprobs[12]*log(q[6])
                    + StoreQuarts[nn]->spprobs[13]*log(q[7]) + StoreQuarts[nn]->spprobs[8]*log(q[8]) + StoreQuarts[nn]->spprobs[10]*log(q[8]) + StoreQuarts[nn]->spprobs[9]*log(q[8]) + StoreQuarts[nn]->spprobs[11]*log(q[8]) + StoreQuarts[nn]->spprobs[14]*log(q[9]);

  //if (verbose==1) printf("The likelihood of the quartet is %f\n\n",quart_lik);
  return(quart_lik);

}


double AsymmetricQuartetLikelihood_ratevar(int nn){
                
  int i, k, sump=0;
  double m = 4.0/3.0, t=2.0*theta; // t is 2*theta
  double abeta[11][1];
  double q[12], check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3;

  GetRateParams();
                        
   // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  //printf("Times are %f %f %f\n",t1,t2,t3);   

  for (i=1; i<12; i++) q[i]=0.0;
	for (k=0; k<ncat; k++) {

		// printf("for category %d, rate is %f\n",k,rvals[k+1]);

        	abeta[1][0]=1.0; 
        	abeta[2][0]=exp(-2.0*m*rvals[k+1]*t1);
        	abeta[3][0]=exp(-2.0*m*rvals[k+1]*t2); 
        	abeta[4][0]=exp(-m*rvals[k+1]*t1)*exp(-2.0*m*rvals[k+1]*t2); 
        	abeta[5][0]=exp(-2.0*m*rvals[k+1]*t3); 
        	abeta[6][0]=exp(-m*rvals[k+1]*t1)*exp(-2.0*m*rvals[k+1]*t3); 
        	abeta[7][0]=exp(-2.0*m*rvals[k+1]*t1)*exp(-2.0*m*rvals[k+1]*t3); 
        	abeta[8][0]=exp(-m*rvals[k+1]*t2)*exp(-2.0*m*rvals[k+1]*t3); 
        	abeta[9][0]=exp(-m*rvals[k+1]*t1)*exp(-m*rvals[k+1]*t2)*exp(-2.0*m*rvals[k+1]*t3);
        	abeta[10][0]=exp((2.0/t)*(rvals[k+1]*t1-rvals[k+1]*t2))*exp(-2.0*m*(rvals[k+1]*t2+rvals[k+1]*t3)); 

        	// order is: xxxx; xxxy = xxyx; xyxx; yxxx; xxyy; xyxy = yxxy; xxyz
        	// yzxx; xyxz = xyzx; yxxz = yxzx; xyzw
        	// to get site patterns probs, multiply amat*abeta (not transpose)
        	// weighted version is below

        	for (i=1; i<11; i++) {
            		q[1] += (1.0/ncat)*4*amat[1][i]*abeta[i][0]; //4
            		q[2] += (1.0/ncat)*12*amat[2][i]*abeta[i][0]; //24
            		q[3] += (1.0/ncat)*12*amat[3][i]*abeta[i][0]; //12
            		q[4] += (1.0/ncat)*12*amat[4][i]*abeta[i][0]; //12
            		q[5] += (1.0/ncat)*12*amat[5][i]*abeta[i][0]; //12
            		q[6] += (1.0/ncat)*12*amat[6][i]*abeta[i][0]; //24
            		q[7] += (1.0/ncat)*24*amat[7][i]*abeta[i][0]; //24
            		q[8] += (1.0/ncat)*24*amat[8][i]*abeta[i][0]; //24
            		q[9] += (1.0/ncat)*24*amat[9][i]*abeta[i][0];  //48
            		q[10] += (1.0/ncat)*24*amat[10][i]*abeta[i][0]; //48
            		q[11] += (1.0/ncat)*24*amat[11][i]*abeta[i][0]; //24
        	}

  }

  // compute likelihood
  quart_lik = StoreQuarts[nn]->spprobs[0]*log(q[1]) + StoreQuarts[nn]->spprobs[1]*log(q[2]) + StoreQuarts[nn]->spprobs[2]*log(q[2]) + StoreQuarts[nn]->spprobs[3]*log(q[3]) + StoreQuarts[nn]->spprobs[4]*log(q[4])
                    + StoreQuarts[nn]->spprobs[7]*log(q[5]) + StoreQuarts[nn]->spprobs[5]*log(q[6]) + StoreQuarts[nn]->spprobs[6]*log(q[6]) + StoreQuarts[nn]->spprobs[12]*log(q[7]) + StoreQuarts[nn]->spprobs[13]*log(q[8])
                    + StoreQuarts[nn]->spprobs[8]*log(q[9]) + StoreQuarts[nn]->spprobs[9]*log(q[9]) + StoreQuarts[nn]->spprobs[10]*log(q[10]) + StoreQuarts[nn]->spprobs[11]*log(q[10]) + StoreQuarts[nn]->spprobs[14]*log(q[11]);

  //if (verbose==1) printf("The likelihood of the quartet is %f\n\n",quart_lik);
  return(quart_lik);

}



double GetCompLik_ratevar() {
   
  int i, j, k, l, m=0;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  GetRateParams();

  num_unique_quarts = 0;

  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          //printf("Quartet %d %d %d %d\n",i,j,k,l);
          //printf("\t %llu %d %d\n",subset_to_index(i-1,j-1,k-1,l-1),(int)subset_to_index(i-1,j-1,k-1,l-1),qvec[(int)subset_to_index$
          //printf("\t Have we looked at this quartet yet? %d\n",qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]);
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {
        
		GetQuartetTree(i,j,k,l,ovec);
                //printf("ovec is: %d %d %d %d %d\n",ovec[0],ovec[1],ovec[2],ovec[3],ovec[4]);
                FindDupQuarts(ovec[0],duppvec);
            
                // Fill the information for the current quartet into the struct StoreQuarts
                for (m=0; m<15; m++) StoreQuarts[num_unique_quarts]->spprobs[m] = duppvec[m];
                StoreQuarts[num_unique_quarts]->ncherries = ovec[0];
                StoreQuarts[num_unique_quarts]->t1 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+3])];
                if (ovec[0]==2) {
                        StoreQuarts[num_unique_quarts]->t2 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+2],ppTwoRowQuart[1][ntaxa+2])];
                        StoreQuarts[num_unique_quarts]->t3 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+2])];
                }
                else if (ovec[0]==1) {
                        StoreQuarts[num_unique_quarts]->t2 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+2])];
                        StoreQuarts[num_unique_quarts]->t3 = &TimeVec[FindMRCA(ppTwoRowQuart[0][ntaxa+3],ppTwoRowQuart[1][ntaxa+1])];
                }

                num_unique_quarts+=1;

          }
        }
      }
    }
  }

  for (i=0; i<num_unique_quarts; i++) {

        //printf("%d ",i);
        //for (j=0; j<15; j++) printf("%d ",StoreQuarts[i]->spprobs[j]);
        //printf("%f %f %f ",*StoreQuarts[i]->t1,*StoreQuarts[i]->t2,*StoreQuarts[i]->t3);
        //printf("%d\n",StoreQuarts[i]->ncherries);

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_ratevar(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood_ratevar(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        //printf("\t %f\n",complik);

  }

  //printf("The composite likelihood is %f; this required %d quartet comps\n",complik,num_unique_quarts);
  return(complik);

  
}



