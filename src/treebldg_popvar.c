/*************************************************/
/*** Function to perform local rearrangements  ***/
/*** and make reject/accept decision           ***/
/*************************************************/

void trbldg_popvar() {

  int i, j, k, index, target, parent, left_child, right_child, sib, side_target, root_ind;
  int unsuccess_est;
  double min_time, old_lik, curr_lik, prop_lik, ran_num;
  double par_time, child_time, prop_time, diff;
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

     //printf("parent is %d, left_child is %d, right child is %d, sib is %d, side_target is %d\n",parent,left_child,right_child,sib,side_target);
  
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


    //printf("Index was %d, Tree has been re-arranged :) The proposed tree is\n",index);
    /*for (i=0; i<ntaxa-1; i++) {
        printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
        printf("%f %f\n ",TimeVec[i+ntaxa+1],TimeVec_temp[i+ntaxa+1]);
    }*/


    /*  Generation of new time for the target node */
    if (TimeVec[ppTwoRow[0][target-(ntaxa+1)]] > TimeVec[ppTwoRow[1][target-(ntaxa+1)]]) min_time = TimeVec[ppTwoRow[0][target-(ntaxa+1)]];
    else min_time = TimeVec[ppTwoRow[1][target-(ntaxa+1)]];
    //printf("min_time is %10.8f, max time is %10.8f\n",min_time,TimeVec[parent]);
    diff = TimeVec[parent]-min_time;
    //printf("diff is %f\n",diff);
    if (diff > 0.0001) {
    	TimeVec[target] = genunf((float)min_time+diff/4,(float)TimeVec[parent]-diff/4);
    	//printf("New time of target %d is %10.8f\n",target,TimeVec[target]);

    	/* Compute CL of new tree and make decision */
    	prop_lik = GetCompLik_popvar();

    	//printf("curr lik is %f, prop lik is %f\n",curr_anneal_lik,prop_lik);

    	if (prop_lik > curr_anneal_lik) {
	
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
		//printf("Random number is %f, ci is %f, accept prob is %f\n",ran_num,ci,exp((prop_lik-curr_anneal_lik)/ci));

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

        	}

    	}

		/*printf("Returning with value:\n");
		for (i=0; i<ntaxa-1; i++) {
        		printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
        		printf("%f %f\n ",TimeVec[i+ntaxa+1],TimeVec_temp[i+ntaxa+1]);
  		}*/

    	//printf("\t curr_anneal_lik is %f, max_cl is %f\n",curr_anneal_lik,max_cl); 
        
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


/* Uphill stochastic search for a new root time */

void RootUpdate() {

   int i, ntries = 5;
   double curr_time, par_time, child_time, prop_time, prop_lik;

  printf("Before update, root time is %f and likelihood is %f, ",TimeVec[ntaxa+1],curr_anneal_lik);

  for (i=0; i<ntries; i++) {

   	curr_time = TimeVec[ntaxa+1];
   	par_time = 2*TimeVec[ntaxa+1];
   	if (TimeVec[ppTwoRow[0][0]]<TimeVec[ppTwoRow[1][0]]) child_time = TimeVec[ppTwoRow[1][0]];
   	else child_time = TimeVec[ppTwoRow[0][0]];
 
   	prop_time = ranf()*(par_time-child_time);
   	TimeVec[ntaxa+1]  =  child_time + prop_time;

   	prop_lik = GetCompLik();

   	if (curr_anneal_lik > prop_lik) TimeVec[ntaxa+1]  = curr_time;
   	else {
		curr_anneal_lik = prop_lik;
		TimeVec_temp[ntaxa+1] = TimeVec[ntaxa+1];
	}

  }

  printf("After update, root time is %f and likelihood is %f\n",TimeVec[ntaxa+1],curr_anneal_lik);

}


void RescaleTree_popvar() {

  int i, j, ntries=5;
  float multiplier;
  double prop_lik;

  for (i=0; i<ntries; i++) {

	multiplier = gennor(1.0,0.01);
	for (j=0; j<2*ntaxa+1; j++) TimeVec[j] = multiplier*TimeVec[j];
	prop_lik = GetCompLik_popvar();

	if (curr_anneal_lik > prop_lik) for (j=0; j<2*ntaxa+1; j++) TimeVec[j] = TimeVec_temp[j];
	else {
		for (j=0; j<2*ntaxa+1; j++) TimeVec_temp[j] = TimeVec[j];
		curr_anneal_lik = prop_lik;
	}

 }

 

}













