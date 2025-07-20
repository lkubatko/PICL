/********************************************************************/
/***  Functions to write trees to files                           ***/
/********************************************************************/

/***  Write saved trees to file  ***/

void write_species_tree(int node, int previous_node) {

  if (ppTwoRow[0][node-(ntaxa+1)]<=ntaxa && ppTwoRow[1][node-(ntaxa+1)]<=ntaxa) {

    //printf("(%d",ppTwoRow[0][node-(ntaxa+1)]);
    printf("(%s",taxname[ppTwoRow[0][node-(ntaxa+1)]-1]);
    printf(":%f",fabs((-1.0)*(TimeVec[ppTwoRow[0][node-(ntaxa+1)]]-TimeVec[node])));
    //printf(",%d",ppTwoRow[1][node-(ntaxa+1)]);
    printf(",%s",taxname[ppTwoRow[1][node-(ntaxa+1)]-1]);
    printf(":%f",fabs((-1.0)*(TimeVec[ppTwoRow[1][node-(ntaxa+1)]]-TimeVec[node])));
    printf(")");
    if (node != previous_node) {
      printf(":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

  if (ppTwoRow[0][node-(ntaxa+1)]<=ntaxa && ppTwoRow[1][node-(ntaxa+1)]>ntaxa) {

    //printf("(%d",ppTwoRow[0][node-(ntaxa+1)]);
    printf("(%s",taxname[ppTwoRow[0][node-(ntaxa+1)]-1]);
    printf(":%f",fabs((-1.0)*(TimeVec[ppTwoRow[0][node-(ntaxa+1)]]-TimeVec[node])));
    printf(",");
    write_species_tree(ppTwoRow[1][node-(ntaxa+1)],node);
    printf(")");
    if (node != previous_node) {
      printf(":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

  if (ppTwoRow[0][node-(ntaxa+1)]>ntaxa && ppTwoRow[1][node-(ntaxa+1)]<=ntaxa) {

    //printf("(%d",ppTwoRow[1][node-(ntaxa+1)]);
    printf("(%s",taxname[ppTwoRow[1][node-(ntaxa+1)]-1]);
    printf(":%f",fabs((-1.0)*(TimeVec[ppTwoRow[1][node-(ntaxa+1)]]-TimeVec[node])));
    printf(",");
    write_species_tree(ppTwoRow[0][node-(ntaxa+1)],node);
    printf(")");
    if (node != previous_node) {
      printf(":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

  if (ppTwoRow[0][node-(ntaxa+1)]>ntaxa && ppTwoRow[1][node-(ntaxa+1)]>ntaxa) {

    printf("(");
    write_species_tree(ppTwoRow[0][node-(ntaxa+1)],node);
    printf(",");
    write_species_tree(ppTwoRow[1][node-(ntaxa+1)],node);
    printf(")");
    if (node != previous_node) {
      printf(":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

}


/***********************************************************/
/*** Write tree to file out (outtree.tre)                ***/
/***********************************************************/

void write_species_tree_out(int node, int previous_node) {

  if (ppTwoRow[0][node-(ntaxa+1)]<=ntaxa && ppTwoRow[1][node-(ntaxa+1)]<=ntaxa) {

    //printf("(%d",ppTwoRow[0][node-(ntaxa+1)]);
    fprintf(out,"(%s",taxname[ppTwoRow[0][node-(ntaxa+1)]-1]);
    fprintf(out,":%f",fabs((-1.0)*(TimeVec[ppTwoRow[0][node-(ntaxa+1)]]-TimeVec[node])));
    //printf(",%d",ppTwoRow[1][node-(ntaxa+1)]);
    fprintf(out,",%s",taxname[ppTwoRow[1][node-(ntaxa+1)]-1]);
    fprintf(out,":%f",fabs((-1.0)*(TimeVec[ppTwoRow[1][node-(ntaxa+1)]]-TimeVec[node])));
    fprintf(out,")");
    if (node != previous_node) {
      fprintf(out,":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

  if (ppTwoRow[0][node-(ntaxa+1)]<=ntaxa && ppTwoRow[1][node-(ntaxa+1)]>ntaxa) {

    //printf("(%d",ppTwoRow[0][node-(ntaxa+1)]);
    fprintf(out,"(%s",taxname[ppTwoRow[0][node-(ntaxa+1)]-1]);
    fprintf(out,":%f",fabs((-1.0)*(TimeVec[ppTwoRow[0][node-(ntaxa+1)]]-TimeVec[node])));
    fprintf(out,",");
    write_species_tree_out(ppTwoRow[1][node-(ntaxa+1)],node);
    fprintf(out,")");
    if (node != previous_node) {
      fprintf(out,":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

  if (ppTwoRow[0][node-(ntaxa+1)]>ntaxa && ppTwoRow[1][node-(ntaxa+1)]<=ntaxa) {

    //printf("(%d",ppTwoRow[1][node-(ntaxa+1)]);
    fprintf(out,"(%s",taxname[ppTwoRow[1][node-(ntaxa+1)]-1]);
    fprintf(out,":%f",fabs((-1.0)*(TimeVec[ppTwoRow[1][node-(ntaxa+1)]]-TimeVec[node])));
    fprintf(out,",");
    write_species_tree_out(ppTwoRow[0][node-(ntaxa+1)],node);
    fprintf(out,")");
    if (node != previous_node) {
      fprintf(out,":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

  if (ppTwoRow[0][node-(ntaxa+1)]>ntaxa && ppTwoRow[1][node-(ntaxa+1)]>ntaxa) {

    fprintf(out,"(");
    write_species_tree_out(ppTwoRow[0][node-(ntaxa+1)],node);
    fprintf(out,",");
    write_species_tree_out(ppTwoRow[1][node-(ntaxa+1)],node);
    fprintf(out,")");
    if (node != previous_node) {
      fprintf(out,":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

}


/***********************************************************/
/*** Write tree to file picltrees.tre                    ***/
/***********************************************************/

void write_species_tree_out_file(int node, int previous_node) {

  if (ppTwoRow[0][node-(ntaxa+1)]<=ntaxa && ppTwoRow[1][node-(ntaxa+1)]<=ntaxa) {

    //printf("(%d",ppTwoRow[0][node-(ntaxa+1)]);
    fprintf(pt,"(%s",taxname[ppTwoRow[0][node-(ntaxa+1)]-1]);
    fprintf(pt,":%f",fabs((-1.0)*(TimeVec[ppTwoRow[0][node-(ntaxa+1)]]-TimeVec[node])));
    //printf(",%d",ppTwoRow[1][node-(ntaxa+1)]);
    fprintf(pt,",%s",taxname[ppTwoRow[1][node-(ntaxa+1)]-1]);
    fprintf(pt,":%f",fabs((-1.0)*(TimeVec[ppTwoRow[1][node-(ntaxa+1)]]-TimeVec[node])));
    fprintf(pt,")");
    if (node != previous_node) {
      fprintf(pt,":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

  if (ppTwoRow[0][node-(ntaxa+1)]<=ntaxa && ppTwoRow[1][node-(ntaxa+1)]>ntaxa) {

    //printf("(%d",ppTwoRow[0][node-(ntaxa+1)]);
    fprintf(pt,"(%s",taxname[ppTwoRow[0][node-(ntaxa+1)]-1]);
    fprintf(pt,":%f",fabs((-1.0)*(TimeVec[ppTwoRow[0][node-(ntaxa+1)]]-TimeVec[node])));
    fprintf(pt,",");
    write_species_tree_out_file(ppTwoRow[1][node-(ntaxa+1)],node);
    fprintf(pt,")");
    if (node != previous_node) {
      fprintf(pt,":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

  if (ppTwoRow[0][node-(ntaxa+1)]>ntaxa && ppTwoRow[1][node-(ntaxa+1)]<=ntaxa) {

    //printf("(%d",ppTwoRow[1][node-(ntaxa+1)]);
    fprintf(pt,"(%s",taxname[ppTwoRow[1][node-(ntaxa+1)]-1]);
    fprintf(pt,":%f",fabs((-1.0)*(TimeVec[ppTwoRow[1][node-(ntaxa+1)]]-TimeVec[node])));
    fprintf(pt,",");
    write_species_tree_out_file(ppTwoRow[0][node-(ntaxa+1)],node);
    fprintf(pt,")");
    if (node != previous_node) {
      fprintf(pt,":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

  if (ppTwoRow[0][node-(ntaxa+1)]>ntaxa && ppTwoRow[1][node-(ntaxa+1)]>ntaxa) {

    fprintf(out,"(");
    write_species_tree_out_file(ppTwoRow[0][node-(ntaxa+1)],node);
    fprintf(pt,",");
    write_species_tree_out_file(ppTwoRow[1][node-(ntaxa+1)],node);
    fprintf(pt,")");
    if (node != previous_node) {
      fprintf(pt,":%f",fabs((-1.0)*(TimeVec[node]-TimeVec[previous_node])));
      }

  }

}



/***********************************************************/
/*** Function to print data in data.phy in PHYLIP format ***/
/***********************************************************/

void print_PHYLIP() {

  int i, j, rr;

  printf("\n\n");

   for (i=0; i<ntaxa; i++) {

     printf("%-10s",psname[i]);
     printf(" ");

     for (j=0; j<num_unique; j++){

       for (rr=0; rr<site_counter[j]; rr++) {

       if (ppBase_unique[i][j]==0) printf("A");
       if (ppBase_unique[i][j]==1) printf("G");
       if (ppBase_unique[i][j]==2) printf("C");
       if (ppBase_unique[i][j]==3) printf("T");
       if (ppBase_unique[i][j]==4) printf("-");
       if (ppBase_unique[i][j]==5) printf("M");
       if (ppBase_unique[i][j]==6) printf("R");
       if (ppBase_unique[i][j]==7) printf("W");
       if (ppBase_unique[i][j]==8) printf("S");
       if (ppBase_unique[i][j]==9) printf("Y");
       if (ppBase_unique[i][j]==10) printf("K");
       if (ppBase_unique[i][j]==11) printf("B");
       if (ppBase_unique[i][j]==12) printf("D");
       if (ppBase_unique[i][j]==13) printf("H");
       if (ppBase_unique[i][j]==14) printf("V");
       if (ppBase_unique[i][j]==15) printf("N");			      

       }

     }

     printf("\n");

   }

   printf("\n\n");

}


void remove_CONSTANT() {

   int i, j, state;
   int total_sites = 0;

   for (i=0; i<num_unique; i++) {

     state = ppBase_unique[0][i];
     j = 0;

     while (ppBase_unique[j][i] == state && j<ntaxa) j+=1; 

     if (j == ntaxa) site_counter[i] = 0;
     
     total_sites = total_sites + site_counter[i];
	
    }

    /*printf("\n");
    for (i=0; i<num_unique; i++) printf("%d ",site_counter[i]);
    printf("\n");*/

    printf("The total number of sites after removal is %d\n\n",total_sites);

}



/**** Code to map quartets to unique ordered integers   ****/
/**** Generated by ChatGPT on 12/30/24 with following:  ****/
/**** Find a one-to-one mapping of sets of 4 integers   ****/
/**** selected from n integers to the numbers 0 through ****/
/**** n choose 4 -1                                     ****/

// Function to compute the binomial coefficient C(n, k)
unsigned long long binomial(int n, int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    unsigned long long res = 1;
    for (int i = 1; i <= k; i++) {
        res = res * (n - i + 1) / i;
    }
    return res;
}

// Function to compute the index of a subset {a, b, c, d}
unsigned long long subset_to_index(int a, int b, int c, int d) {
    return binomial(a, 1) + binomial(b, 2) + binomial(c, 3) + binomial(d, 4);
}

// Function to compute the subset corresponding to a given index
void index_to_subset(unsigned long long index, int *a, int *b, int *c, int *d, int n) {
    *a = 0;
    while (binomial(*a + 1, 1) <= index) (*a)++;
    index -= binomial(*a, 1);

    *b = *a + 1;
    while (binomial(*b + 1, 2) <= index) (*b)++;
    index -= binomial(*b, 2);

    *c = *b + 1;
    while (binomial(*c + 1, 3) <= index) (*c)++;
    index -= binomial(*c, 3);

    *d = *c + 1;
    while (binomial(*d + 1, 4) <= index) (*d)++;
    // No need to subtract for the final element
}

