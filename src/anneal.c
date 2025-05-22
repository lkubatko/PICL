/*************************************************************/
/*** Functions to carry out the simulated annealing search ***/
/*************************************************************/

/*** Annealing to search tree space for full data model ***/

void anneal_full() {

   int i, counteri, burnin;
   double curr_lik, prop_lik, max_change=0.0, U;

   counteri = 1;
   burnin = 100;

   /* Step 0: Make sure that temp time vector matches starting time vector */
   for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i];

   /* Step 1: estimate U, the upper bound on the change in likelihood */

   curr_lik = GetCompLik();
   U = curr_lik;
   ci = U/(1+counteri*beta);
   while (counteri<burnin) {
	trbldg();
	prop_lik = GetCompLik();
	if (fabs(prop_lik-curr_lik)>max_change) max_change=fabs(prop_lik-curr_lik);
      	counteri++;
   }
   U = max_change;

  /* Step 2: anneal */
  counteri = 1;
  while (counteri<(int)max_it/2) {
	trbldg();
	ci = U/(1+counteri*beta);
        counteri += 1;
  }
  bl_uphill_full();
  while (counteri<max_it) {
        trbldg();
        ci = U/(1+counteri*beta);
        counteri += 1;
  }

}

/*** Annealing to search tree space for the rate variation model ***/

void anneal_ratevar() {

   int i, counteri, burnin;
   double curr_lik, prop_lik, max_change=0.0, U;

   counteri = 1;
   burnin = 100;

   /* Step 0: Make sure that temp time vector matches starting time vector */
   for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i];  

   /* Step 1: estimate U, the upper bound on the change in likelihood */

   curr_lik = GetCompLik_ratevar();
   U = curr_lik;
   ci = U/(1+counteri*beta);
   while (counteri<burnin) {

        trbldg_ratevar();
        prop_lik = GetCompLik_ratevar();
        if (fabs(prop_lik-curr_lik)>max_change) max_change=fabs(prop_lik-curr_lik);
        counteri++;
   }
   U = max_change;

  /* Step 2: anneal */
  counteri = 1;
  while (counteri<(int)max_it/2) {
        trbldg_ratevar();
        ci = U/(1+counteri*beta);
        counteri += 1;
  }
  bl_uphill_ratevar();
  while (counteri<max_it) {
        trbldg_ratevar();
        ci = U/(1+counteri*beta);
        counteri += 1;
  }

}


/*** Annealing to search tree space for the SNP model ***/

void anneal_msnp() {

   int i, counteri, burnin;
   double curr_lik, prop_lik, max_change=0.0, U;

   counteri = 1;
   burnin = 100;

   /* Step 0: Make sure that temp time vector matches starting time vector */
   for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i];  

   /* Step 1: estimate U, the upper bound on the change in likelihood */

   curr_lik = GetCompLik_msnp();
   U = curr_lik;
   ci = U/(1+counteri*beta);
   while (counteri<burnin) {

        trbldg_msnp();
        prop_lik = GetCompLik_msnp();
        if (fabs(prop_lik-curr_lik)>max_change) max_change=fabs(prop_lik-curr_lik);
        counteri++;
   }
   U = max_change;


  /* Step 2: anneal */
  counteri = 1;
  while (counteri<(int)max_it/2) {
        trbldg_msnp();
        ci = U/(1+counteri*beta);
        counteri += 1;
  }
  bl_uphill_msnp();
  while (counteri<max_it) {
        trbldg_msnp();   
        ci = U/(1+counteri*beta);
        counteri += 1;
  }

}




/*** Annealing for branch lengths in a fixed tree ***/
/*** Full data model                              ***/

void bl_anneal_full() {

  int it, rnode;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;

  int i, j, k, l, m;
  int ovec[5], duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
           if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
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

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        //printf("\t %f\n",complik);

  }
 
  curr_lik = complik;

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        curr_time = TimeVec[rnode+(ntaxa+1)];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(rnode+(ntaxa+1))];
               
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
      
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

		complik = 0.0;
		for (i=0; i<num_unique_quarts; i++) {

        		if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood(i);
        		else complik+=AsymmetricQuartetLikelihood(i);
                	
  		}

                prop_lik = complik;
      
                if (curr_lik > prop_lik) {
                	//printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) TimeVec[rnode+(ntaxa+1)]  = curr_time;
			else curr_lik = prop_lik;
		}
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }
        else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
 
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

		complik = 0.0;
		for (i=0; i<num_unique_quarts; i++) {

        		if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood(i);
        		else complik+=AsymmetricQuartetLikelihood(i);
          	}

                prop_lik = complik;

                if (curr_lik > prop_lik) {
                	//printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                     	if (exp((prop_lik-curr_lik)/ci) < ranf()) TimeVec[rnode+(ntaxa+1)]  = curr_time;
			else curr_lik = prop_lik;
		}
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }

  }

}



/*** Uphill search for branch lengths in a fixed tree ***/
/*** Full data model                                  ***/

void bl_uphill_full() {

  int it, rnode;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;

  int i, j, k, l, m;
  int ovec[5], duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
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

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}

  }

  curr_lik = complik;

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        curr_time = TimeVec[rnode+(ntaxa+1)];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(rnode+(ntaxa+1))];
               
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
      
                prop_time = ranf()*(par_time-child_time);
		if (child_time + prop_time > 0.00001) TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

		complik = 0.0;
		for (i=0; i<num_unique_quarts; i++) {

        		if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood(i);
        		else complik+=AsymmetricQuartetLikelihood(i);
     		}

                prop_lik = complik;
      
                if (curr_lik > prop_lik) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                else curr_lik = prop_lik;
        }
        else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
 
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

		complik = 0.0;
		for (i=0; i<num_unique_quarts; i++) {

        		if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood(i);
        		else complik+=AsymmetricQuartetLikelihood(i);
     		}

                prop_lik = complik;

                if (curr_lik > prop_lik) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                else curr_lik = prop_lik;
        }

  }

}





/*** Annealing for branch lengths in a fixed tree ***/
/*** Multi-allelic SNP data model                 ***/

void bl_anneal_msnp() {

  int it, rnode;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik;

  int i, j, k, l, m;
  int ovec[5], duppvec[15];
  double complik = 0.0, U;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
           if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
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

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_msnp(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood_msnp(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        //printf("\t %f\n",complik);

  }
 
  curr_lik = complik;

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        curr_time = TimeVec[rnode+(ntaxa+1)];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(rnode+(ntaxa+1))];
               
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
      
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

                complik = 0.0;
                for (i=0; i<num_unique_quarts; i++) {

                        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_msnp(i);
                        else complik+=AsymmetricQuartetLikelihood_msnp(i);
                        
                }

                prop_lik = complik;
      
                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                        else curr_lik = prop_lik;
                }
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }

       else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
 
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

                complik = 0.0;
                for (i=0; i<num_unique_quarts; i++) {

                        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_msnp(i);
                        else complik+=AsymmetricQuartetLikelihood_msnp(i);
                }

                prop_lik = complik;

                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                        else curr_lik = prop_lik;
                }
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }

  }

}



/*** Uphill search for branch lengths in a fixed tree ***/
/*** Multi-allelic SNP data model                     ***/

void bl_uphill_msnp() {

int it, rnode;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;

  int i, j, k, l, m;
  int ovec[5], duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
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

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_msnp(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood_msnp(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}

  }

  curr_lik = complik;

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        curr_time = TimeVec[rnode+(ntaxa+1)];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(rnode+(ntaxa+1))];
               
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
      
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

                complik = 0.0;
                for (i=0; i<num_unique_quarts; i++) {

                        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_msnp(i);
                        else complik+=AsymmetricQuartetLikelihood_msnp(i);
                }

                prop_lik = complik;
      
                if (curr_lik > prop_lik) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                else curr_lik = prop_lik;
        }

       else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
 
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

                complik = 0.0;
                for (i=0; i<num_unique_quarts; i++) {

                        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_msnp(i);
                        else complik+=AsymmetricQuartetLikelihood_msnp(i);
                }

                prop_lik = complik;

                if (curr_lik > prop_lik) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                else curr_lik = prop_lik;
        }

  }



}



/*** Annealing for branch lengths in a fixed tree ***/
/*** Rate variation data model                    ***/

void bl_anneal_ratevar() {

  int it, rnode;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;

  int i, j, k, l, m;
  int ovec[5], duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
           if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
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

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_ratevar(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood_ratevar(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        //printf("\t %f\n",complik);

  }
 
  curr_lik = complik;

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        curr_time = TimeVec[rnode+(ntaxa+1)];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(rnode+(ntaxa+1))];
               
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
      
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

                complik = 0.0;
                for (i=0; i<num_unique_quarts; i++) {

                        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_ratevar(i);
                        else complik+=AsymmetricQuartetLikelihood_ratevar(i);
                        
                }

                prop_lik = complik;
      
                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                        else curr_lik = prop_lik;
                }
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }

        else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
 
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

                complik = 0.0;
                for (i=0; i<num_unique_quarts; i++) {

                        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_ratevar(i);
                        else complik+=AsymmetricQuartetLikelihood_ratevar(i);
                }

                prop_lik = complik;

                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                        else curr_lik = prop_lik;
                }
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }

  }

}



/*** Uphill search for branch lengths in a fixed tree ***/
/*** Rate variation data model                        ***/

void bl_uphill_ratevar() {

int it, rnode;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik;

  int i, j, k, l, m;
  int ovec[5], duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
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

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_ratevar(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood_ratevar(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}

  }

  curr_lik = complik;

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        curr_time = TimeVec[rnode+(ntaxa+1)];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(rnode+(ntaxa+1))];
               
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
      
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

                complik = 0.0;
                for (i=0; i<num_unique_quarts; i++) {

                        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_ratevar(i);
                        else complik+=AsymmetricQuartetLikelihood_ratevar(i);
                }

                prop_lik = complik;
      
                if (curr_lik > prop_lik) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                else curr_lik = prop_lik;
        }

       else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];
 
                prop_time = ranf()*(par_time-child_time);
                TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

                complik = 0.0;
                for (i=0; i<num_unique_quarts; i++) {

                        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_ratevar(i);
                        else complik+=AsymmetricQuartetLikelihood_ratevar(i);
                }

                prop_lik = complik;

                if (curr_lik > prop_lik) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                else curr_lik = prop_lik;
        }

  }


}
