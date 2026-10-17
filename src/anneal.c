/*************************************************************/
/*** Functions to carry out the simulated annealing search ***/
/*************************************************************/

/*** Annealing to search tree space for full data model ***/

void anneal_full() {

   int i, counteri, burnin, samplesize, startcount;
   int bound1, max_it_bl_switch;
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

   curr_lik = GetCompLik();
   U = curr_anneal_lik;
   ci = U/(1+counteri*beta);
   while (counteri<burnin) {

	trbldg();
	if (fabs(curr_lik-curr_anneal_lik)>max_change) max_change=fabs(curr_lik-curr_anneal_lik);
      	counteri++;

   }
   U = max_change;
   max_cl = curr_anneal_lik;
   b1old = 0.0;
   printf("\n   Burnin is complete.  \n   Starting first annealing ...");

  /* Step 2: annealing 1 */
  num_reject = 0;
  curr_anneal_lik = GetCompLik();
  while (num_reject < bound1 && counteri < burnin+(int)max_it/2) {
	trbldg();
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
  printf(" first annealing is complete.\n");

  /* Step 3: branch length optimization */
  //bl_uphill_full();
  //for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i]; 

  /* Step 4: annealing 2 */
  printf("   Starting second annealing ...");
  num_reject = 0;
  startcount = counteri;
  while (num_reject < mult_iter*bound1 && counteri < burnin+max_it) {

        trbldg();
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
  printf(" second annealing is complete.\n");

  /* Step 5: branch length optimization of estimated tree */
  bl_uphill_full();
  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i]; 

  printf("   Final optimization is complete. \n   The total number of iterations used for simulated annealing was %d.\n\n",counteri);
  if (counteri == max_it) printf("WARNING: maximum number of iterations reached without satisfying the stopping criterion.\n\n");

}

/*** Annealing to search tree space for the rate variation model ***/

void anneal_ratevar() {

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

   OptAlpha();
   curr_lik = GetCompLik_ratevar();
   U = curr_anneal_lik;
   ci = U/(1+counteri*beta);
   while (counteri<burnin) {
        trbldg_ratevar();
        if (fabs(curr_lik-curr_anneal_lik)>max_change) max_change=fabs(curr_lik-curr_anneal_lik);
        counteri++;
   }
   U = max_change;
   max_cl = curr_anneal_lik;
   b1old = 0.0;
   printf("\n   Burnin is complete. \n   Starting first annealing ...");

  /* Step 2: annealing 1 */
  num_reject = 0;
  curr_anneal_lik = GetCompLik_ratevar();
  while (num_reject < bound1 && counteri < burnin+(int)max_it/2) {
        trbldg_ratevar();
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
                RescaleTree_ratevar();
		OptAlpha();
        }
  }
  printf(" first annealing is complete.\n");

  /* Step 3: branch length optimization */
  bl_uphill_ratevar();
  OptAlpha();
  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i];   

  /* Step 4: annealing 2 */
  printf("   Starting second annealing ...");
  num_reject = 0;
  startcount = counteri;
  while (num_reject < mult_iter*bound1 && counteri < burnin+max_it) {
        trbldg_ratevar();
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
                RescaleTree_ratevar();
		OptAlpha();
        }
  }

  printf(" second annealing is complete.\n");
  /* Step 5: branch length optimization of estimated tree */
  bl_uphill_ratevar();
  OptAlpha();
  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i]; 

  printf("   Final optimization is complete. \n   The total number of iterations used for simulated annealing was %d.\n\n",counteri);
  if (counteri == max_it) printf("WARNING: maximum number of iterations reached without satisfying the stopping criterion.\n\n");

}


/*** Annealing to search tree space for the SNP model ***/

void anneal_msnp() {

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

   curr_lik = GetCompLik_msnp();
   U = curr_anneal_lik;
   ci = U/(1+counteri*beta);
   while (counteri<burnin) {

        trbldg_msnp();
        if (fabs(curr_lik-curr_anneal_lik)>max_change) max_change=fabs(curr_lik-curr_anneal_lik);
        counteri++;
   }
   U = max_change;
   max_cl = curr_anneal_lik;
   b1old = 0.0;
   printf("\n   Burnin is complete. \n   Starting first annealing ...");

  /* Step 2: annealing 1 */
  num_reject = 0;
  curr_anneal_lik = GetCompLik_msnp();
  while (num_reject < bound1 && counteri < burnin+(int)max_it/2) {
        trbldg_msnp();
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
                RescaleTree_msnp();
        }
  }
  printf(" first annealing is complete.\n");

 /* Step 3: branch length optimization */
  bl_uphill_msnp();
  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i]; 

  /* Step 4: annealing 2 */
  printf("   Starting second annealing ...");
  num_reject = 0;
  startcount = counteri;
  while (num_reject < mult_iter*bound1 && counteri < burnin+max_it) {
        trbldg_msnp();   
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
                RescaleTree_msnp();
        }

  }
  printf(" second annealing is complete.\n");

  /* Step 5: branch length optimization of estimated tree */
  bl_uphill_ratevar();
  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i]; 

  printf("   Final optimization is complete. \n   The total number of iterations used for simulated annealing was %d.\n\n",counteri);
  if (counteri == max_it) printf("WARNING: maximum number of iterations reached without satisfying the stopping criterion.\n\n");

}


void anneal_genetree() {

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

   curr_lik = GetCompLik_genetree();
   U = curr_anneal_lik;
   ci = U/(1+counteri*beta);
   while (counteri<burnin) {
	
        trbldg_genetree();
        if (fabs(curr_lik-curr_anneal_lik)>max_change) max_change=fabs(curr_lik-curr_anneal_lik);
        counteri++;
   }
   U = max_change;
   max_cl = curr_anneal_lik;
   b1old = 0.0;
   printf("\n   Burnin is complete. \n   Starting first annealing ...");

  /* Step 2: annealing 1 */
  num_reject = 0;
  curr_anneal_lik = GetCompLik_genetree();
  while (num_reject < bound1 && counteri < burnin+(int)max_it/2) {
        trbldg_genetree();
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
                RescaleTree_genetree();
        }
  }
  printf(" first annealing is complete.\n");

  /* Step 3: branch length optimization */
  bl_uphill_genetree();
  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i]; 

  /* Step 4: annealing 2 */
  printf("   Starting second annealing ...");
  num_reject = 0;
  startcount = counteri;
  while (num_reject < mult_iter*bound1 && counteri < burnin+max_it) {
        trbldg_genetree();
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
                RescaleTree_genetree();
        }
  }
  printf(" second annealing is complete.\n");

  /* Step 5: branch length optimization of estimated tree */
  bl_uphill_genetree();
  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i]; 

  printf("   Final optimization is complete. \n   The total number of iterations used for simulated annealing was %d.\n\n",counteri);
  if (counteri == max_it) printf("WARNING: maximum number of iterations reached without satisfying the stopping criterion.\n\n");

}



/*** Annealing for branch lengths in a fixed tree ***/
/*** Full data model                              ***/

void bl_anneal_full() {

  int it, rnode, node, dk, qi;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;
  QuartDeps qd;
  double *quart_lik_cache, *old_vals;

  int i, j, k, l, m;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;
  make_indmat();

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
           if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
                 FindDupQuarts(ovec[0],duppvec,num_unique_quarts+1);

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

  quart_lik_cache = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (i=0; i<num_unique_quarts; i++) {

        if (StoreQuarts[i]->ncherries==2) quart_lik_cache[i] = SymmetricQuartetLikelihood(i);
        else if (StoreQuarts[i]->ncherries==1) quart_lik_cache[i] = AsymmetricQuartetLikelihood(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        complik += quart_lik_cache[i];

  }

  curr_lik = complik;

  qd = BuildQuartDeps(2*ntaxa+1);
  old_vals = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        node = rnode+(ntaxa+1);
        curr_time = TimeVec[node];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(node)];

                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood(qi) : AsymmetricQuartetLikelihood(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                	//printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) {
                                TimeVec[node]  = curr_time;
                                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                        }
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

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood(qi) : AsymmetricQuartetLikelihood(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                	//printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) {
                                TimeVec[node]  = curr_time;
                                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                        }
			else curr_lik = prop_lik;
		}
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }

  }

  free(old_vals);
  FreeQuartDeps(&qd);
  free(quart_lik_cache);

  if (curr_lik > max_cl) {
      max_cl = curr_lik;

      for (i=0; i<ntaxa-1; i++) {
          ppTwoRow_best[0][i] = ppTwoRow[0][i];
          ppTwoRow_best[1][i] = ppTwoRow[1][i];
      }
      for (i=0; i<2*ntaxa+1; i++) {
          TimeVec_best[i] = TimeVec[i];
      }
  }

}



/*** Uphill search for branch lengths in a fixed tree ***/
/*** Full data model                                  ***/

void bl_uphill_full() {

  int it, rnode, node, dk, qi;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;
  QuartDeps qd;
  double *quart_lik_cache, *old_vals;

  int i, j, k, l, m;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;
  make_indmat();

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
                FindDupQuarts(ovec[0],duppvec,num_unique_quarts+1);

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

  /* quart_lik_cache[i] always holds SymmetricQuartetLikelihood(i) or
     AsymmetricQuartetLikelihood(i) evaluated at the *current* TimeVec state -
     i.e. exactly what a fresh call would return right now. Below we only
     ever refresh the entries whose dependencies changed (see QuartDeps
     comment above GetCompLik), and always sum the full cache in ascending
     index order, so the composite log-likelihood is computed with the exact
     same terms in the exact same order as the original full-recompute loop
     - bit-identical, not just mathematically-equivalent. */
  quart_lik_cache = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (i=0; i<num_unique_quarts; i++) {

        if (StoreQuarts[i]->ncherries==2) quart_lik_cache[i] = SymmetricQuartetLikelihood(i);
        else if (StoreQuarts[i]->ncherries==1) quart_lik_cache[i] = AsymmetricQuartetLikelihood(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        complik += quart_lik_cache[i];

  }

  curr_lik = complik;

  qd = BuildQuartDeps(2*ntaxa+1);
  old_vals = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        node = rnode+(ntaxa+1);
        curr_time = TimeVec[node];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(node)];

                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

		if (par_time - child_time > 0.00001 && prop_time > 0.001) TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood(qi) : AsymmetricQuartetLikelihood(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        TimeVec[node]  = curr_time;
                        for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                }
                else curr_lik = prop_lik;
        }
        else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood(qi) : AsymmetricQuartetLikelihood(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        TimeVec[node]  = curr_time;
                        for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                }
                else curr_lik = prop_lik;

        }

  }

  free(old_vals);
  FreeQuartDeps(&qd);
  free(quart_lik_cache);

  if (curr_lik > max_cl) {
      max_cl = curr_lik;

      for (i=0; i<ntaxa-1; i++) {
          ppTwoRow_best[0][i] = ppTwoRow[0][i];
          ppTwoRow_best[1][i] = ppTwoRow[1][i];
      }
      for (i=0; i<2*ntaxa+1; i++) {
          TimeVec_best[i] = TimeVec[i];
      }
  }

  for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i];

}





/*** Annealing for branch lengths in a fixed tree ***/
/*** Multi-allelic SNP data model                 ***/

void bl_anneal_msnp() {

  int it, rnode, node, dk, qi;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik;
  QuartDeps qd;
  double *quart_lik_cache, *old_vals;

  int i, j, k, l, m;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0, U;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;
  make_indmat();

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
           if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
                 FindDupQuarts(ovec[0],duppvec,num_unique_quarts+1);

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

  quart_lik_cache = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

    for (i=0; i<num_unique_quarts; i++) {

        if (StoreQuarts[i]->ncherries==2) quart_lik_cache[i] = SymmetricQuartetLikelihood_msnp(i);
        else if (StoreQuarts[i]->ncherries==1) quart_lik_cache[i] = AsymmetricQuartetLikelihood_msnp(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        complik += quart_lik_cache[i];

  }

  curr_lik = complik;

  qd = BuildQuartDeps(2*ntaxa+1);
  old_vals = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        node = rnode+(ntaxa+1);
        curr_time = TimeVec[node];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(node)];

                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood_msnp(qi) : AsymmetricQuartetLikelihood_msnp(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) {
                                TimeVec[node]  = curr_time;
                                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                        }
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

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood_msnp(qi) : AsymmetricQuartetLikelihood_msnp(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) {
                                TimeVec[node]  = curr_time;
                                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                        }
                        else curr_lik = prop_lik;
                }
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }

  }

  free(old_vals);
  FreeQuartDeps(&qd);
  free(quart_lik_cache);

}



/*** Uphill search for branch lengths in a fixed tree ***/
/*** Multi-allelic SNP data model                     ***/

void bl_uphill_msnp() {

int it, rnode, node, dk, qi;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;
  QuartDeps qd;
  double *quart_lik_cache, *old_vals;

  int i, j, k, l, m;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;
  make_indmat();

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
                FindDupQuarts(ovec[0],duppvec,num_unique_quarts+1);

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

  quart_lik_cache = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (i=0; i<num_unique_quarts; i++) {

        if (StoreQuarts[i]->ncherries==2) quart_lik_cache[i] = SymmetricQuartetLikelihood_msnp(i);
        else if (StoreQuarts[i]->ncherries==1) quart_lik_cache[i] = AsymmetricQuartetLikelihood_msnp(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        complik += quart_lik_cache[i];

  }

  curr_lik = complik;

  qd = BuildQuartDeps(2*ntaxa+1);
  old_vals = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        node = rnode+(ntaxa+1);
        curr_time = TimeVec[node];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(node)];

                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood_msnp(qi) : AsymmetricQuartetLikelihood_msnp(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        TimeVec[node]  = curr_time;
                        for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                }
                else curr_lik = prop_lik;
        }

       else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood_msnp(qi) : AsymmetricQuartetLikelihood_msnp(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        TimeVec[node]  = curr_time;
                        for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                }
                else curr_lik = prop_lik;
        }

  }

  free(old_vals);
  FreeQuartDeps(&qd);
  free(quart_lik_cache);

}



/*** Annealing for branch lengths in a fixed tree ***/
/*** Rate variation data model                    ***/

void bl_anneal_ratevar() {

  int it, rnode, node, dk, qi;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;
  QuartDeps qd;
  double *quart_lik_cache, *old_vals;

  int i, j, k, l, m;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;
  make_indmat();

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
           if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
                 FindDupQuarts(ovec[0],duppvec,num_unique_quarts+1);

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

  quart_lik_cache = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (i=0; i<num_unique_quarts; i++) {

        if (StoreQuarts[i]->ncherries==2) quart_lik_cache[i] = SymmetricQuartetLikelihood_ratevar(i);
        else if (StoreQuarts[i]->ncherries==1) quart_lik_cache[i] = AsymmetricQuartetLikelihood_ratevar(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        complik += quart_lik_cache[i];

  }

  curr_lik = complik;

  qd = BuildQuartDeps(2*ntaxa+1);
  old_vals = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        node = rnode+(ntaxa+1);
        curr_time = TimeVec[node];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(node)];

                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood_ratevar(qi) : AsymmetricQuartetLikelihood_ratevar(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) {
                                TimeVec[node]  = curr_time;
                                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                        }
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

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood_ratevar(qi) : AsymmetricQuartetLikelihood_ratevar(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) {
                                TimeVec[node]  = curr_time;
                                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                        }
                        else curr_lik = prop_lik;
                }
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }

  }

  free(old_vals);
  FreeQuartDeps(&qd);
  free(quart_lik_cache);

}



/*** Uphill search for branch lengths in a fixed tree ***/
/*** Rate variation data model                        ***/

void bl_uphill_ratevar() {

int it, rnode, node, dk, qi;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik;
  QuartDeps qd;
  double *quart_lik_cache, *old_vals;

  int i, j, k, l, m;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;
  make_indmat();

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
                FindDupQuarts(ovec[0],duppvec,num_unique_quarts+1);

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

  quart_lik_cache = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (i=0; i<num_unique_quarts; i++) {

        if (StoreQuarts[i]->ncherries==2) quart_lik_cache[i] = SymmetricQuartetLikelihood_ratevar(i);
        else if (StoreQuarts[i]->ncherries==1) quart_lik_cache[i] = AsymmetricQuartetLikelihood_ratevar(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
        complik += quart_lik_cache[i];

  }

  curr_lik = complik;

  qd = BuildQuartDeps(2*ntaxa+1);
  old_vals = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        node = rnode+(ntaxa+1);
        curr_time = TimeVec[node];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(node)];

                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood_ratevar(qi) : AsymmetricQuartetLikelihood_ratevar(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        TimeVec[node]  = curr_time;
                        for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                }
                else curr_lik = prop_lik;
        }

       else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = (StoreQuarts[qi]->ncherries==2) ? SymmetricQuartetLikelihood_ratevar(qi) : AsymmetricQuartetLikelihood_ratevar(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        TimeVec[node]  = curr_time;
                        for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                }
                else curr_lik = prop_lik;
        }

  }

  free(old_vals);
  FreeQuartDeps(&qd);
  free(quart_lik_cache);

}


void bl_anneal_genetree(){

  int it, rnode, node, dk, qi;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;
  QuartDeps qd;
  double *quart_lik_cache, *old_vals;

  int i, j, k, l, m;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;
  make_indmat();

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
           if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
                FindDupQuarts(ovec[0],duppvec,num_unique_quarts+1);

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

  quart_lik_cache = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (i=0; i<num_unique_quarts; i++) {
        quart_lik_cache[i] = QuartetLikelihood_genetree(i);
        complik += quart_lik_cache[i];
  }

  curr_lik = complik;

  qd = BuildQuartDeps(2*ntaxa+1);
  old_vals = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        node = rnode+(ntaxa+1);
        curr_time = TimeVec[node];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(node)];

                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = QuartetLikelihood_genetree(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) {
                                TimeVec[node]  = curr_time;
                                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                        }
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

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = QuartetLikelihood_genetree(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        //printf("acceptance probability is %f\n",exp((prop_lik-curr_lik)/ci));
                        if (exp((prop_lik-curr_lik)/ci) < ranf()) {
                                TimeVec[node]  = curr_time;
                                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                        }
                        else curr_lik = prop_lik;
                }
                else curr_lik = prop_lik;
                ci = U/(1+it*beta);
        }

  }

  free(old_vals);
  FreeQuartDeps(&qd);
  free(quart_lik_cache);

}


void bl_uphill_genetree() {

  int it, rnode, node, dk, qi;
  double par_time, child_time;
  double curr_time, prop_time;
  double curr_lik, prop_lik, U;
  QuartDeps qd;
  double *quart_lik_cache, *old_vals;

  int i, j, k, l, m;
  int ovec[5];
  double duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;
  make_indmat();

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

                GetQuartetTree(i,j,k,l,ovec);
                FindDupQuarts(ovec[0],duppvec,num_unique_quarts+1);

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

  quart_lik_cache = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (i=0; i<num_unique_quarts; i++) {
        quart_lik_cache[i] = QuartetLikelihood_genetree(i);
        complik += quart_lik_cache[i];
  }

  curr_lik = complik;

  qd = BuildQuartDeps(2*ntaxa+1);
  old_vals = (double*)malloc((num_unique_quarts > 0 ? num_unique_quarts : 1)*sizeof(double));

  for (it=0; it<max_it_bl; it++) {

        /* pick a node at random */
        rnode = floor((ntaxa-1)*ranf());
        node = rnode+(ntaxa+1);
        curr_time = TimeVec[node];
        //printf("random node is %d\n",rnode);
        //printf("time of random node is %f\n",TimeVec[rnode+(ntaxa+1)]);

        if (rnode!=0) {  // non-root node
                par_time = TimeVec[find_parent(node)];

                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                if (child_time + prop_time > 0.00001) TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = QuartetLikelihood_genetree(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        TimeVec[node]  = curr_time;
                        for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                }
                else curr_lik = prop_lik;
        }
        else {  // adjust the time of the root
                par_time = 2*TimeVec[ntaxa+1];
                if (TimeVec[ppTwoRow[0][rnode]]<TimeVec[ppTwoRow[1][rnode]]) child_time = TimeVec[ppTwoRow[1][rnode]];
                else child_time = TimeVec[ppTwoRow[0][rnode]];

                prop_time = ranf()*(par_time-child_time);

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) old_vals[dk-qd.offset[node]] = quart_lik_cache[qd.dep_quart[dk]];

                TimeVec[node]  =  child_time + prop_time;

                for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) {
                        qi = qd.dep_quart[dk];
                        quart_lik_cache[qi] = QuartetLikelihood_genetree(qi);
                }

                prop_lik = 0.0;
                for (i=0; i<num_unique_quarts; i++) prop_lik += quart_lik_cache[i];

                if (curr_lik > prop_lik) {
                        TimeVec[node]  = curr_time;
                        for (dk=qd.offset[node]; dk<qd.offset[node+1]; dk++) quart_lik_cache[qd.dep_quart[dk]] = old_vals[dk-qd.offset[node]];
                }
                else curr_lik = prop_lik;
        }

  }

  free(old_vals);
  FreeQuartDeps(&qd);
  free(quart_lik_cache);

}


