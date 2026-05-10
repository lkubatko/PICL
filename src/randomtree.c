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
  int num_times_assigned;
  double tmax=1.0;
  double *RandomTimes;
  double scale;

  /* allocate memory */
  
  eLength = (2*ntaxa+1)*sizeof(double);
  eeLength = (2*ntaxa+1)*sizeof(double);
  fLength = (ntaxa-2)*sizeof(double);
  ffLength = (2*ntaxa+1)*sizeof(int);
  iLenseq = (ntaxa)*sizeof(int);

  // keep track of which times have been assigned
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

  scale = ((int)ntaxa)*theta;

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

  // assign times to internal nodes
  num_times_assigned = 0;
  for (i=ntaxa+1; i<2*ntaxa; i++) TimesAssigned[i] = 0;
  for (i=0; i<ntaxa+1; i++) TimesAssigned[i] = 1;
	
  while (num_times_assigned < ntaxa-2) {

	for (i=0; i<ntaxa-1; i++) {

		if (TimesAssigned[ppTwoRow[0][i]]==1 && TimesAssigned[ppTwoRow[1][i]]==1 && TimesAssigned[i+ntaxa+1]==0) {

			TimeVec[i+ntaxa+1] = scale*(num_times_assigned+1)/ntaxa;
			TimesAssigned[i+ntaxa+1] = 1;
			num_times_assigned++;

	  	}
 	}

  }
 
  TimeVec[ntaxa+1] = scale;
  
}




/* The following code was developed using ChatGPT on Jan 24, 2026 */

double compute_times(int node){

   int c1, c2;
   double t1, t2, t;
   double scale;

   scale = ((int)ntaxa/2)*theta;

   if (node <= ntaxa) {
	TimeVec[node] = 0.0;
	return 0.0;
   }
 
   c1 = ppTwoRow[0][node-(ntaxa+1)];
   c2 = ppTwoRow[1][node-(ntaxa+1)];
   t1 = compute_times(c1);
   t2 = compute_times(c2);

   t = 0.5 * (1.0 + (t1 > t2 ? t1 : t2));
   TimeVec[node] = scale*t;

  return t;
	
}
