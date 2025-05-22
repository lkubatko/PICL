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

  for (i=0; i<max_it_bl; i++) {

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
