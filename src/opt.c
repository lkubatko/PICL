/**********************************************/
/*** Functions to optimize model parameters ***/
/**********************************************/

/*** Optimize the rate variation parameter (alpha) ***/

void OptAlpha() {

  int i;
  double curr_lik, new_lik;
  float alpha, new_alpha;

  curr_lik = GetCompLik_ratevar();
  //printf("starting lik is %f\n",curr_lik);

  //for (i=0; i<max_it_bl; i++) {
  for (i=0; i<100; i++) {

        alpha = ratepar;
	new_alpha = (double)ratepar*gengam((float)10.0,(float)10.0);
	//new_alpha = (double)genunf((float)ratepar-0.05,(float)ratepar+0.05);
	if (new_alpha>1000.0) new_alpha = 1000.0;  // max = 1000 => no rate var
        ratepar = new_alpha;
        //printf("proposed alpha is %f ",ratepar);
        GetRateParams();
  	new_lik = GetCompLik_ratevar();
  	//printf("New lik is %f\n",new_lik);

        if (new_lik > curr_lik) {
	
	   //printf("\t accepting -- new alpha is %f\n",ratepar);
	   curr_lik = new_lik;

	}	
	else ratepar = alpha;

  }

}


/**********************************************/
/*** Node time optimization                 ***/
/**********************************************/


void OptTimeQuick(int node, int reps) {

  int i;
  double mintime, maxtime, proptime;
  double derivs[2], oldlik, lik;
  
  if (TimeVec[ppTwoRow[0][node-(ntaxa+1)]] > TimeVec[ppTwoRow[1][node-(ntaxa+1)]]) mintime = TimeVec[ppTwoRow[0][node-(ntaxa+1)]];
  else mintime = TimeVec[ppTwoRow[1][node-(ntaxa+1)]];
  if (node!=ntaxa+1) maxtime = TimeVec[parents[node]];
  else maxtime = 2*TimeVec[node];
  printf("node is %d, mintime is %f, maxtime is %f\n",node,mintime,maxtime);

  for (i=0; i<reps; i++) {

        printf("%f ",TimeVec[node]);
	oldlik = GetCompLik();
	GetCompLikDerivatives(node,derivs);
	proptime = TimeVec[node] - derivs[0]/derivs[1];
	if (proptime>mintime && proptime<maxtime) { printf(" yes "); 
	TimeVec[node] = proptime; }
	lik = GetCompLik();

	printf("%f %f %f\n",TimeVec[node],oldlik,lik);

  }

  TimeVec_temp[node] = TimeVec[node];


 /*for (i=0; i<1000; i++) {  
	TimeVec[node] = TimeVec[node]+0.00001;
	printf("%f %f\n",TimeVec[node],GetCompLik());
 }*/


}


void OptAllTimes(int quick, int max, double tol) {

  int i, j;
  double diff, oldlik, newlik;  

  diff = 1000;
  RescaleTree();

  while (diff>tol) {

	oldlik = GetCompLik();

  	for (i=ntaxa+2; i<2*ntaxa; i++) {

		printf("\n\nNode %d:\n",i);	
		OptTimeQuick(i,max);

  	}

  	newlik = GetCompLik();
	diff = newlik-oldlik;
	printf("\t Difference is %f\n",diff);
	RescaleTree();
  }

}
