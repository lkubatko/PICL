/*************************************************/
/*** Function to perform local rearrangements  ***/
/*** and make reject/accept decision           ***/
/*************************************************/

void trbldg_ratevar() {

  int i, j, k, index, target, parent, left_child, right_child, sib, side_target;
  int unsuccess_est, root_ind;
  double min_time, old_lik, curr_lik, prop_lik, ran_num, diff;
  double alpha, new_alpha;
  float ran_est;

  old_lik = curr_anneal_lik;
  root_ind = 0;
 
  /* Randomly choose an internal node as the target node */
  
  target = (floor)(ranf()*(ntaxa-2))+ntaxa+2;
  if (target<ntaxa+2 || target>2*ntaxa) {
    printf("\t\t Warning: There  was a target generating problem in file trbldg.c at line=%d",__LINE__);
    target = ntaxa+2;
  }


  /* Identify parent, right_child, left_child, and sib for     */
  /* the target.  Also let the variable side_target            */
  /* be zero if the target is on the parent's left side        */
  /* and one if the target is on the parent's right side.      */
  /* ppTwoRow_temp and TimeVec_temp contain the current values */
  
  parent = parents[target];
  left_child = ppTwoRow[0][target-(ntaxa+1)];
  right_child = ppTwoRow[1][target-(ntaxa+1)];
  if (ppTwoRow[0][parent-(ntaxa+1)] == target) {
    
    sib = ppTwoRow[1][parent-(ntaxa+1)];
    side_target=0;
    
  }
  
  else {
    
    sib = ppTwoRow[0][parent-(ntaxa+1)];
    side_target = 1;
    
  }

  /* Generate 1, 2, or 3 at random to decide which of the three  */
  /* possible local rearrangements will be made.  1 = no change, */
  /* 2 = move left_child to branch connecting sib and parent,    */
  /* 3  = move right_child to branch connecting sib and parent   */
  
  index=(floor)(ranf()*3)+1;
 
  if (index>3 || index<1) {
    printf("\t\t Warning: There was an index generating problem in file trbldg.c at line=%d",__LINE__);
    index = 1;
  }

  if (index==2) {
   
    ppTwoRow[side_target][parent-(ntaxa+1)]=right_child;
    ppTwoRow[(side_target+1)%2][parent-(ntaxa+1)]=target;
    ppTwoRow[side_target][target-(ntaxa+1)]=left_child;
    ppTwoRow[(side_target+1)%2][target-(ntaxa+1)]=sib;
    parents[right_child]=parent;
    parents[sib]=target;
    
  }
  
  else { 
    if (index==3){
      
      ppTwoRow[side_target][parent-(ntaxa+1)]=left_child;
      ppTwoRow[(side_target+1)%2][parent-(ntaxa+1)]=target;
      ppTwoRow[side_target][target-(ntaxa+1)]=right_child;
      ppTwoRow[(side_target+1)%2][target-(ntaxa+1)]=sib;
      parents[left_child]=parent;
      parents[sib]=target;
      
    }
  }


/* printf("Index was %d, Tree has been re-arranged :) The proposed tree is\n",index);
 for (i=0; i<ntaxa-1; i++) {
        printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
        printf("%f %f\n ",TimeVec[i+ntaxa+1],TimeVec_temp[i+ntaxa+1]);
  }
*/

  
  /*  Generation of new time for the target node */
  if (TimeVec_temp[ppTwoRow[0][target-(ntaxa+1)]] > TimeVec[ppTwoRow[1][target-(ntaxa+1)]]) min_time = TimeVec[ppTwoRow[0][target-(ntaxa+1)]];
  else min_time = TimeVec[ppTwoRow[1][target-(ntaxa+1)]];
  //printf("min_time is %f, max time is %f\n",min_time,TimeVec[parent]);
  diff = TimeVec[parent]-min_time;
  if (diff > 0.0001) {
  	TimeVec[target] = genunf((float)min_time+diff/4,(float)TimeVec[parent]-diff/4);
  	//printf("New time of target is %f\n",TimeVec[target]);

  	/* Propose a new value of the rate parameter */
  	alpha = ratepar;
  	//new_alpha = (double)ratepar*gengam((float)10.0,(float)10.0);
  	//ratepar = new_alpha;
  	//GetRateParams();

  	/* Compute CL of new time and make decision */
  	prop_lik = GetCompLik_ratevar();

  	if (prop_lik > curr_lik) {
	
        	//printf("ACCEPT\n");
		num_reject=0;
		curr_anneal_lik = prop_lik;
		ppTwoRow_temp[0][parent-(ntaxa+1)] = ppTwoRow[0][parent-(ntaxa+1)];
    		ppTwoRow_temp[1][parent-(ntaxa+1)] = ppTwoRow[1][parent-(ntaxa+1)];
    		ppTwoRow_temp[0][target-(ntaxa+1)] = ppTwoRow[0][target-(ntaxa+1)];
     		ppTwoRow_temp[1][target-(ntaxa+1)] = ppTwoRow[1][target-(ntaxa+1)];
    		TimeVec_temp[target] = TimeVec[target];
   		parents_temp[left_child] = parents[left_child];
    		parents_temp[right_child] = parents[right_child];
    		parents_temp[sib] = parents[sib];

  	}
  	else {

		/* make accept/reject decision */
		ran_num = ranf();
		//printf("Random number is %f, accept prob is %f\n",ran_num,exp((prop_lik-curr_anneal_lik)/ci));

		if (ran_num<=exp((prop_lik-curr_anneal_lik)/ci)) { /* accept */
			//printf("ACCEPT lower\n");
			num_reject=0; 
        		curr_anneal_lik = prop_lik;
        		ppTwoRow_temp[0][parent-(ntaxa+1)] = ppTwoRow[0][parent-(ntaxa+1)];
        		ppTwoRow_temp[1][parent-(ntaxa+1)] = ppTwoRow[1][parent-(ntaxa+1)];
        		ppTwoRow_temp[0][target-(ntaxa+1)] = ppTwoRow[0][target-(ntaxa+1)];                                                                                
        		ppTwoRow_temp[1][target-(ntaxa+1)] = ppTwoRow[1][target-(ntaxa+1)];
        		TimeVec_temp[target] = TimeVec[target];
        		parents_temp[left_child] = parents[left_child];
        		parents_temp[right_child] = parents[right_child];
        		parents_temp[sib] = parents[sib];
		}

		else {
			//printf("Reject\n");
			num_reject++;
			curr_anneal_lik = old_lik;
     			ppTwoRow[0][parent-(ntaxa+1)] = ppTwoRow_temp[0][parent-(ntaxa+1)];
    			ppTwoRow[1][parent-(ntaxa+1)] = ppTwoRow_temp[1][parent-(ntaxa+1)];
    			ppTwoRow[0][target-(ntaxa+1)] = ppTwoRow_temp[0][target-(ntaxa+1)];
    			ppTwoRow[1][target-(ntaxa+1)] = ppTwoRow_temp[1][target-(ntaxa+1)];
    			TimeVec[target] = TimeVec_temp[target];
    			parents[left_child] = parents_temp[left_child];
    			parents[right_child] = parents_temp[right_child];
    			parents[sib] = parents_temp[sib];
			ratepar = alpha;
			GetRateParams();
        	}

  	}

	/*printf("Returning with value:\n");
	for (i=0; i<ntaxa-1; i++) {
        	printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
        	printf("%f %f\n ",TimeVec[i+ntaxa+1],TimeVec_temp[i+ntaxa+1]);
  	}*/

  	if (curr_anneal_lik > max_cl) {

		for (i=0; i<ntaxa-1; i++) {
            		ppTwoRow_best[0][i] = ppTwoRow[0][i];
            		ppTwoRow_best[1][i] = ppTwoRow[1][i];
  		}
  		for (i=0; i<2*ntaxa+1; i++) TimeVec_best[i] = TimeVec[i]; 
		max_cl = curr_anneal_lik;
  	}
    }
    else {

        ppTwoRow[0][parent-(ntaxa+1)] = ppTwoRow_temp[0][parent-(ntaxa+1)];
        ppTwoRow[1][parent-(ntaxa+1)] = ppTwoRow_temp[1][parent-(ntaxa+1)];
        ppTwoRow[0][target-(ntaxa+1)] = ppTwoRow_temp[0][target-(ntaxa+1)];
        ppTwoRow[1][target-(ntaxa+1)] = ppTwoRow_temp[1][target-(ntaxa+1)];
        TimeVec[target] = TimeVec_temp[target];
        parents[left_child] = parents_temp[left_child];
        parents[right_child] = parents_temp[right_child];
        parents[sib] = parents_temp[sib];

  }


}


void RescaleTree_ratevar() {

  int i, j, ntries=5;
  float multiplier;
  double prop_lik, alpha, new_alpha;

  for (i=0; i<ntries; i++) {

	/* rescale the tree */
        multiplier = gennor(1.0,0.01);
        for (j=0; j<2*ntaxa+1; j++) TimeVec[j] = multiplier*TimeVec[j];
        prop_lik = GetCompLik_ratevar();

        if (curr_anneal_lik > prop_lik) for (j=0; j<2*ntaxa+1; j++) TimeVec[j] = TimeVec_temp[j];
        else {
                for (j=0; j<2*ntaxa+1; j++) TimeVec_temp[j] = TimeVec[j];
                curr_anneal_lik = prop_lik;
        }

 }

}






















