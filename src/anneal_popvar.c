/*************************************************************/
/*** Functions to carry out the simulated annealing search ***/
/*************************************************************/

/*** Annealing to search tree space for full data model ***/

void anneal_full_popvar() {

   int i, counteri, burnin, samplesize, startcount;
   int bound1;
   double curr_lik, max_change=0.0, U;
   double xsum=0.0, ysum=0.0, xsumsq=0.0, covarsum=0.0;
   double b0, b1, b1old, R;    

   /* set params */
   counteri = 1;
   burnin = 1000;
   bound1 = floor(log(prob_bound/ntaxa)/log(1.0-1.0/(ntaxa-2)));

   /* Step 0: Make sure that temp time vector matches starting time vector */
   for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i];

   /* Step 1: estimate U, the upper bound on the change in likelihood */

   curr_lik = GetCompLik_popvar();
   U = curr_anneal_lik;
   ci = U/(1+counteri*beta);
   while (counteri<burnin) {

	trbldg_popvar();
	if (fabs(curr_lik-curr_anneal_lik)>max_change) max_change=fabs(curr_lik-curr_anneal_lik);
      	counteri++;

   }
   U = max_change;
   max_cl = curr_anneal_lik;
   b1old = 0.0;
   printf("\nburnin is complete ... \nstarting first annealing ...\n");


  /* Step 2: annealing 1 */
  num_reject = 0;
  curr_anneal_lik = GetCompLik_popvar();
  while (num_reject < bound1 && counteri < burnin+(int)max_it/2) {
	trbldg_popvar();
	ci = U/(1+counteri*beta);
	xsum += log(counteri);
        ysum += log(-1.0*curr_anneal_lik);
        xsumsq += log(counteri)*log(counteri);
        covarsum += log(counteri)*log(-1.0*curr_anneal_lik);
	counteri += 1;

	if (counteri%test_increment == 0) { /* adaptively update beta */
		samplesize = counteri-burnin;
		b1 = (samplesize*covarsum-xsum*ysum)/(samplesize*xsumsq-xsum*xsum);
		R = fabs(b1-b1opt)/fabs(b1old-b1opt);
		if (R>1) {
			if (b1<b1opt) beta = (2*beta)/R; /* b1 below b1opt */
			else if (b1<-1.0*b1opt) beta = R*beta; /* b1 between b1opt and -b1opt */
			else beta = beta/(R);  /* b1 above -b1opt */
		}
		b1old = b1;
		RescaleTree();
	}
  }
  printf("first annealing is complete ...\n");

  /* Step 3: branch length optimization */
  bl_uphill_full_popvar();
  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i]; 

  /* Step 4: annealing 2 */
  printf("starting second annealing ...\n");
  num_reject = 0;
  startcount = counteri;
  while (num_reject < mult_iter*bound1 && counteri < burnin+max_it) {
        trbldg_popvar();
        ci = U/(1+counteri*beta);
	xsum += log(counteri);
        ysum += log(-1.0*curr_anneal_lik);
        xsumsq += log(counteri)*log(counteri);
        covarsum += log(counteri)*log(-1.0*curr_anneal_lik);
        counteri += 1;
	if (counteri%test_increment == 0) { /* adaptively update beta */
		samplesize = counteri-burnin;
        	b1 = (samplesize*covarsum-xsum*ysum)/(samplesize*xsumsq-xsum*xsum);
		R = fabs(b1-b1opt)/fabs(b1old-b1opt);
		if (R>1) {
                        if (b1<b1opt) beta = (2*beta)/R; /* b1 below b1opt */
                        else if (b1<-1.0*b1opt) beta = R*beta; /* b1 between b1opt and -b1opt */
                        else beta = beta/(R);  /* b1 above -b1opt */
                }
		b1old = b1;
		RescaleTree();
        }
  }
  printf("second annealing is complete ...\n");

  /* Step 5: branch length optimization of estimated tree */
  bl_uphill_full_popvar();
  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i]; 

  printf("final optimization is complete ... \nthe total number of iterations used for simulated annealing was %d.\n\n",counteri);
  if (counteri == max_it) printf("WARNING: maximum number of iterations reached without satisfying the stopping criterion.\n\n");

}




/*** Annealing for branch lengths in a fixed tree ***/
/*** Full data model                              ***/

void bl_anneal_full_popvar() {

  int it, rnode;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;

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

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_popvar(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood_popvar(i);
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

        		if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_popvar(i);
        		else complik+=AsymmetricQuartetLikelihood_popvar(i);
                	
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

        		if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_popvar(i);
        		else complik+=AsymmetricQuartetLikelihood_popvar(i);
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

void bl_uphill_full_popvar() {

  int it, rnode;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;

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

        if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_popvar(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood_popvar(i);
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
		if (par_time - child_time > 0.00001 && prop_time > 0.001) TimeVec[rnode+(ntaxa+1)]  =  child_time + prop_time;

		complik = 0.0;
		for (i=0; i<num_unique_quarts; i++) {

        		if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_popvar(i);
        		else complik+=AsymmetricQuartetLikelihood_popvar(i);
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

        		if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_popvar(i);
        		else complik+=AsymmetricQuartetLikelihood_popvar(i);
     		}

                prop_lik = complik;

                if (curr_lik > prop_lik) TimeVec[rnode+(ntaxa+1)]  = curr_time;
                else curr_lik = prop_lik;

        }

  }

  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i];

}
