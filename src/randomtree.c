/***  Function to search through ppTwoRow to find parent ***/
/***  of a specified node.                               ***/
  
int find_parentr(int target) {

int i, j, parent=0;

 for (j=0; j<2; j++) {
   for (i=0; i<ntaxa; i++) {
     if (ppTwoRow[j][i]==target) {
      parent = i+ntaxa+1;
      break;
     }
   }
 }

 if (parent==0) {
 
   printf("Could not find parent in find_parent for target = %d\n",target); 
   exit(0);

}
 
return(parent);

}

/***  Function to find the generation of the input node   ***/
  
int find_genr(int current_node) {

  int gen_counter=0, parent;

    parent = current_node;

    while (parent != ntaxa+1) {

      parent = find_parentr(parent);
      gen_counter++;

    }

    return gen_counter;

}

void build_random_tree(int seed1, int seed2) {

  int i, j;
  int count = 0, found_one = 0;
  int max_geni = 0, place_max = 0;
  int curr_length;
  int *ppNodes[2];
  int *TimesAssigned, *genarray;
  int eLength, eeLength, fLength, ffLength, iLenseq;
  int node1, node2, parent, parent2;
  int new_counter;
  int count_genarray=0, count_gen=0, count_gen2=0;
  int max_gen=0;
  double tmax=1.0;
  double *RandomTimes;
  double scale;

  /* allocate memory */
  
  eLength = (2*ntaxa+1)*sizeof(double);
  eeLength = (2*ntaxa+1)*sizeof(double);
  fLength = (ntaxa-2)*sizeof(double);
  ffLength = (ntaxa-2)*sizeof(int);
  iLenseq = (ntaxa)*sizeof(int);
  
  TimesAssigned = (int*)malloc(ffLength);
  if (TimesAssigned==NULL) printf("Can't memalloc TimesAssigned\n");
  
  RandomTimes = (double*)malloc(fLength);
  if (RandomTimes==NULL) printf("Can't memalloc RandomTimes\n");
  
  genarray = (int*)malloc(eeLength);
  if (genarray==NULL) printf("Can't memalloc genarray\n");
  
  for (i=0; i<2; i++)
    {
      ppNodes[i] = (int*)malloc(iLenseq);
      if (ppNodes[i]==NULL) printf("Can't memalloc ppNodes[%d]\n",i);
    }
  
  /* done memory allocation */

  scale = ((int)ntaxa/2)*theta;

  curr_length=ntaxa;
  new_counter=ntaxa+2;
  seedj = seed1;
  seedk = seed2;
  for (i=0; i<2*ntaxa; i++) TimeVec[i] = 0.0;

  for (i=0; i<curr_length; i++) {
    
    ppNodes[0][i]=i+1;
    ppNodes[1][i]=0;
    
  }
  
  /*  construct the random tree */

  while (count<ntaxa-2) {
    
    /* randomly pick two nodes to join by generating random */
    /* numbers from a uniform distribution                  */
    
    node1 = (int)(ranf()*(curr_length));
    node2 = (int)(ranf()*(curr_length));

    while (node1==node2) node2 = (int)(ranf()*(curr_length));
    
    ppTwoRow[0][new_counter-(ntaxa+1)]=ppNodes[0][node1];
    ppTwoRow[1][new_counter-(ntaxa+1)]=ppNodes[0][node2];
    
    ppNodes[1][node1]=1;
    ppNodes[1][node2]=1;
    
    /* shift */

    found_one = 0;

    for (i=0; i<curr_length; i++ ) {
      
      if (ppNodes[1][i]==1 && found_one==0) {
	
	ppNodes[0][i]=new_counter;
	ppNodes[1][i]=0;
	found_one=1;
	
      }
      
      if (ppNodes[1][i]==1 && found_one==1) found_one=2;
      
      if (found_one==2) {
	
	ppNodes[0][i]=ppNodes[0][i+1];
	
      }
      
    }
    
    curr_length = curr_length-1;
    
    for (i=0; i<curr_length; i++) {
      
      ppNodes[1][i]=0;
      
    }
    
    /* done shift */
    
    count++;
    new_counter++;
    
  }
  
  ppTwoRow[0][0]=ppNodes[0][0];
  ppTwoRow[1][0]=ppNodes[0][1];

  
  /*  equally space times throughout the tree  */

  for (i=0; i<ntaxa;i++) {

    genarray[i] = find_genr(i+1);
    if (genarray[i]>max_gen) {

      parent=i;
      max_gen = genarray[i];

    }

  }

  parent = parent + 1;
  parent2 = parent;

  while (count_genarray<ntaxa) {
    
    while (parent!=ntaxa+1) {
      
      if (TimeVec[parent]!=0.0) break;
      parent = find_parentr(parent);
      count_gen++;

    }
    
    if (parent==ntaxa+1) tmax=1.0;
    else tmax=1.0-TimeVec[parent];
    
    parent = find_parentr(parent2);
    count_gen2 = 1;
    
    while (TimeVec[parent]==0.0 && parent!=ntaxa+1) {
      
      if (TimeVec[parent]==0.0) TimeVec[parent] = scale*(count_gen2*(tmax/((double)count_gen)));
      count_gen2++;
      parent = find_parentr(parent);
      
    }
    
    count_gen = 0;
    count_genarray++;
    genarray[parent2-1]=0;

    /* find_max_gen */

    max_geni = 0;
    place_max = 0;
    for (i=0; i<ntaxa;i++) {

      if (genarray[i]>max_geni) {

	max_geni = genarray[i];
	place_max=i;
	
      }
      
    }

    /* end find_max_gen */
    
    parent = place_max + 1;
    parent2 = parent;
    
  }
  TimeVec[ntaxa+1] = scale;
  
}




void assign_times() {

  int i, j;
  int count = 0;
  int max_geni = 0, place_max = 0;
  int *genarray;
  int eeLength;
  int parent, parent2;
  int count_genarray=0, count_gen=0, count_gen2=0;
  int max_gen=0;
  double tmax=1.0;
  double scale;

  /* allocate memory */
  
  eeLength = (2*ntaxa+1)*sizeof(double);
  
  genarray = (int*)malloc(eeLength);
  if (genarray==NULL) printf("Can't memalloc genarray\n");
  
  /* done memory allocation */

  scale = ((int)ntaxa/2)*theta;

  /*  equally space times throughout the tree  */

  for (i=0; i<ntaxa;i++) {

    genarray[i] = find_genr(i+1);
    if (genarray[i]>max_gen) {

      parent=i;
      max_gen = genarray[i];

    }

  }

  parent = parent + 1;
  parent2 = parent;

  while (count_genarray<ntaxa) {
    
    while (parent!=ntaxa+1) {
      
      if (TimeVec[parent]!=0.0) break;
      parent = find_parentr(parent);
      count_gen++;

    }
    
    if (parent==ntaxa+1) tmax=1.0;
    else tmax=1.0-TimeVec[parent];
    
    parent = find_parentr(parent2);
    count_gen2 = 1;
    
    while (TimeVec[parent]==0.0 && parent!=ntaxa+1) {
      
      if (TimeVec[parent]==0.0) TimeVec[parent] = scale*(count_gen2*(tmax/((double)count_gen)));
      count_gen2++;
      parent = find_parentr(parent);
      
    }
    
    count_gen = 0;
    count_genarray++;
    genarray[parent2-1]=0;

    /* find_max_gen */

    max_geni = 0;
    place_max = 0;
    for (i=0; i<ntaxa;i++) {

      if (genarray[i]>max_geni) {

        max_geni = genarray[i];
        place_max=i;
        
      }
      
    }

    /* end find_max_gen */
    
    parent = place_max + 1;
    parent2 = parent;
    
  }
  TimeVec[ntaxa+1] = scale;

}
