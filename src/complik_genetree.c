// Compute transition probabilities for the JC69 model
double tpii(double brlen) { return(0.25+0.75*exp(-4.0*brlen/3.0)); }
double tpij(double brlen) { return(0.25-0.25*exp(-4.0*brlen/3.0)); }


double QuartetLikelihood_genetree(int nn){

  int i, sump=0;
  double q[15],check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3;
  double pii[5], pij[5];  // stores transition probabilities for the 5 branches

  // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  //printf("Times are %f %f %f\n",t1,t2,t3);
  if (StoreQuarts[nn]->ncherries==2) {
  	pii[0] = tpii(t2); pii[1] = tpii(t2); pii[2] = tpii(2*t3-t1-t2); pii[3] = tpii(t1); pii[4] = tpii(t1);
  	pij[0] = tpij(t2); pij[1] = tpij(t2); pij[2] = tpij(2*t3-t1-t2); pij[3] = tpij(t1); pij[4] = tpij(t1);
  }
  else if (StoreQuarts[nn]->ncherries==1) {
	pii[0] = tpii(2*t3-t2); pii[1] = tpii(t2); pii[2] = tpii(t2-t1); pii[3] = tpii(t1); pii[4] = tpii(t1);
        pij[0] = tpij(2*t3-t2); pij[1] = tpij(t2); pij[2] = tpij(t2-t1); pij[3] = tpij(t1); pij[4] = tpij(t1);
  }
  else {
	printf("Problem detecting symmetry in QuartetLikelihood_genetree ... exiting.\n\n");
	exit(1);
  }

  // In the StoreQuarts counts, the order is: 
  // 0 = xxxx; 1 = xxxy; 2 = xxyx; 3 = xyxx; 4 = yxxx
  // 5 = xyxy; 6 = yxxy; 7 = xxyy; 8 = xyxz; 9 = xyzx; 10 = yxxz
  // 11 = yxzx; 12 = xxyz; 13 = yzxx; 14 = xyzw


  // order from C-K 2015 is: 0 = xxxx; 1 = xxxy; 2 = xxyx; 3 = xyxx; 4 = yxxx; 5 = xyxy; 6 = yxxy
  // 7 = xxyy; 8 = xxyz; 9 = yzxx; 10 = xyxz; 11 = yxxz; 12 = xyzx; 13 = yxzx; 14 = xyzw
  // so what is below is re-ordered as above
  // weighted version is below for consistency with the species tree likelihood
  // note that the difference is only in a constant that is not a function of the parameters
  // however, the gene tree likelihood will not match PAUP* for 4 tips unless the weights are removed below

  for (i=1; i<15; i++) q[i]=0.0;
  // xxxx:
  q[0] = 0.25*(pii[0]*pii[1]*pii[2]*pii[3]*pii[4] + 3*pii[0]*pii[1]*pij[2]*pij[3]*pij[4] + 3*pij[0]*pij[1]*pij[2]*pii[3]*pii[4] + 3*pij[0]*pij[1]*pii[2]*pij[3]*pij[4] 
	+ 3*2*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);  
  // xxxy:
  q[1] = 0.25*(pii[0]*pii[1]*pii[2]*pii[3]*pij[4] + pii[0]*pii[1]*pij[2]*pij[3]*pii[4] + pij[0]*pij[1]*pii[2]*pij[3]*pii[4] + 2*pii[0]*pii[1]*pij[2]*pij[3]*pij[4] 
	+ 3*pij[0]*pij[1]*pij[2]*pii[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pij[3]*pii[4] + 2*pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 4*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // xxyx:
  q[2] = 0.25*(pii[0]*pii[1]*pii[2]*pij[3]*pii[4] + pii[0]*pii[1]*pij[2]*pii[3]*pij[4] + pij[0]*pij[1]*pii[2]*pii[3]*pij[4] + 2*pii[0]*pii[1]*pij[2]*pij[3]*pij[4]
	+ 3*pij[0]*pij[1]*pij[2]*pij[3]*pii[4] + 2*pij[0]*pij[1]*pij[2]*pii[3]*pij[4] + 2*pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 4*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // xyxx:
  q[3] = 0.25*(pii[0]*pij[1]*pii[2]*pii[3]*pii[4] + pij[0]*pii[1]*pij[2]*pii[3]*pii[4] + pij[0]*pii[1]*pii[2]*pij[3]*pij[4] + 2*pij[0]*pii[1]*pij[2]*pij[3]*pij[4]
	+ 3*pii[0]*pij[1]*pij[2]*pij[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pii[3]*pii[4] + 2*pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 4*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // yxxx:
  q[4] = 0.25*(pij[0]*pii[1]*pii[2]*pii[3]*pii[4] + pii[0]*pij[1]*pij[2]*pii[3]*pii[4] + pii[0]*pij[1]*pii[2]*pij[3]*pij[4] + 2*pii[0]*pij[1]*pij[2]*pij[3]*pij[4]
	+ 3*pij[0]*pii[1]*pij[2]*pij[3]*pij[4] +2*pij[0]*pij[1]*pij[2]*pii[3]*pii[4] + 2*pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 4*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // xxyy:
  q[7] = 0.25*(pii[0]*pii[1]*pii[2]*pij[3]*pij[4] + pij[0]*pij[1]*pii[2]*pii[3]*pii[4] + pii[0]*pii[1]*pij[2]*pii[3]*pii[4] + 2*pii[0]*pii[1]*pij[2]*pij[3]*pij[4]
	+ 2*pij[0]*pij[1]*pij[2]*pii[3]*pii[4] + 2*pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 7*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // xyxy:
  q[5] = 0.25*(pii[0]*pij[1]*pii[2]*pii[3]*pij[4] + pii[0]*pij[1]*pij[2]*pij[3]*pii[4] + pij[0]*pii[1]*pij[2]*pii[3]*pij[4] + pij[0]*pii[1]*pii[2]*pij[3]*pii[4] 
	+ 2*pii[0]*pij[1]*pij[2]*pij[3]*pij[4] + 2*pij[0]*pii[1]*pij[2]*pij[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pii[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pij[3]*pii[4] 
	+ 2*pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // yxxy:
  q[6] = 0.25*(pii[0]*pij[1]*pii[2]*pij[3]*pii[4] + pii[0]*pij[1]*pij[2]*pii[3]*pij[4] + pij[0]*pii[1]*pij[2]*pij[3]*pii[4] + pij[0]*pii[1]*pii[2]*pii[3]*pij[4] 
	+ 2*pii[0]*pij[1]*pij[2]*pij[3]*pij[4] + 2*pij[0]*pii[1]*pij[2]*pij[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pii[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pij[3]*pii[4] 
	+ 2*pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // xxyz:
  q[12] = 0.25*(pii[0]*pii[1]*pii[2]*pij[3]*pij[4] + pii[0]*pii[1]*pij[2]*pii[3]*pij[4] + pii[0]*pii[1]*pij[2]*pij[3]*pii[4] + pij[0]*pij[1]*pii[2]*pii[3]*pij[4]
	+ pij[0]*pij[1]*pii[2]*pij[3]*pii[4] + pii[0]*pii[1]*pij[2]*pij[3]*pij[4] + pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pij[3]*pii[4]
	+ 2*pij[0]*pij[1]*pij[2]*pii[3]*pij[4] + 5*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // yzxx: 
  q[13] = 0.25*(pij[0]*pij[1]*pii[2]*pii[3]*pii[4] + pii[0]*pij[1]*pij[2]*pii[3]*pii[4] + pii[0]*pij[1]*pii[2]*pij[3]*pij[4] + pij[0]*pii[1]*pij[2]*pii[3]*pii[4]
	+ pij[0]*pii[1]*pii[2]*pij[3]*pij[4] + pij[0]*pij[1]*pij[2]*pii[3]*pii[4] + pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 2*pii[0]*pij[1]*pij[2]*pij[3]*pij[4] + 
	+ 2*pij[0]*pii[1]*pij[2]*pij[3]*pij[4] + 5*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // xyxz:
  q[8] = 0.25*(pii[0]*pij[1]*pii[2]*pii[3]*pij[4] + pii[0]*pij[1]*pij[2]*pij[3]*pii[4] + pij[0]*pii[1]*pij[2]*pii[3]*pij[4] + pij[0]*pii[1]*pii[2]*pij[3]*pij[4]
	+ pij[0]*pii[1]*pij[2]*pij[3]*pii[4] + pij[0]*pij[1]*pii[2]*pij[3]*pii[4] + 2*pii[0]*pij[1]*pij[2]*pij[3]*pij[4] + 2*pij[0]*pij[1]*pij[2]*pii[3]*pij[4] 
	+ pij[0]*pii[1]*pij[2]*pij[3]*pij[4] + pij[0]*pij[1]*pij[2]*pij[3]*pii[4] + pij[0]*pij[1]*pii[2]*pij[3]*pij[4] + 3*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // yxxz:
  q[10] = 0.25*(pii[1]*pij[0]*pii[2]*pii[3]*pij[4] + pii[1]*pij[0]*pij[2]*pij[3]*pii[4] + pij[1]*pii[0]*pij[2]*pii[3]*pij[4] + pij[1]*pii[0]*pii[2]*pij[3]*pij[4] +
	+ pij[1]*pii[0]*pij[2]*pij[3]*pii[4] + pij[1]*pij[0]*pii[2]*pij[3]*pii[4] + 2*pii[1]*pij[0]*pij[2]*pij[3]*pij[4] + 2*pij[1]*pij[0]*pij[2]*pii[3]*pij[4]
	+ pij[1]*pii[0]*pij[2]*pij[3]*pij[4] + pij[1]*pij[0]*pij[2]*pij[3]*pii[4] + pij[1]*pij[0]*pii[2]*pij[3]*pij[4] + 3*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // xyzx:
  q[9] = 0.25*(pii[0]*pij[1]*pii[2]*pii[4]*pij[3] + pii[0]*pij[1]*pij[2]*pij[4]*pii[3] + pij[0]*pii[1]*pij[2]*pii[4]*pij[3] + pij[0]*pii[1]*pii[2]*pij[4]*pij[3]
	+ pij[0]*pii[1]*pij[2]*pij[4]*pii[3] + pij[0]*pij[1]*pii[2]*pij[4]*pii[3] + 2*pii[0]*pij[1]*pij[2]*pij[4]*pij[3] + 2*pij[0]*pij[1]*pij[2]*pii[4]*pij[3]
	+ pij[0]*pii[1]*pij[2]*pij[4]*pij[3] + pij[0]*pij[1]*pij[2]*pij[4]*pii[3] + pij[0]*pij[1]*pii[2]*pij[4]*pij[3] + 3*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // yxzx:
  q[11] = 0.25*(pii[1]*pij[0]*pii[2]*pii[4]*pij[3] + pii[1]*pij[0]*pij[2]*pij[4]*pii[3] + pij[1]*pii[0]*pij[2]*pii[4]*pij[3] + pij[1]*pii[0]*pii[2]*pij[4]*pij[3]
	+ pij[1]*pii[0]*pij[2]*pij[4]*pii[3] + pij[1]*pij[0]*pii[2]*pij[4]*pii[3] + 2*pii[1]*pij[0]*pij[2]*pij[4]*pij[3] + 2*pij[1]*pij[0]*pij[2]*pii[4]*pij[3]
	+ pij[1]*pii[0]*pij[2]*pij[4]*pij[3] + pij[1]*pij[0]*pij[2]*pij[4]*pii[3] + pij[1]*pij[0]*pii[2]*pij[4]*pij[3] + 3*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);
  // xyzw:
  q[14] = 0.25*(pii[0]*pij[1]*pii[2]*pij[3]*pij[4] + pii[0]*pij[1]*pij[2]*pii[3]*pij[4] + pii[0]*pij[1]*pij[2]*pij[3]*pii[4] + pij[0]*pii[1]*pii[2]*pij[3]*pij[4]
	+ pij[0]*pii[1]*pij[2]*pii[3]*pij[4] + pij[0]*pii[1]*pij[2]*pij[3]*pii[4] + pij[0]*pij[1]*pii[2]*pii[3]*pij[4] + pij[0]*pij[1]*pii[2]*pij[3]*pii[4]
	+ pii[0]*pij[1]*pij[2]*pij[3]*pij[4] + pij[0]*pii[1]*pij[2]*pij[3]*pij[4] + pij[0]*pij[1]*pij[2]*pij[3]*pii[4] + pij[0]*pij[1]*pij[2]*pii[3]*pij[4]
	+ 4*pij[0]*pij[1]*pij[2]*pij[3]*pij[4]);

   /*check_qsum = 0.0;
   for (i=0; i<15; i++) {
	printf("%d %f %d\n",i,q[i],StoreQuarts[nn]->spprobs[i]);
	check_qsum += q[i];
   }
   printf("Sum of qs is: %f\n\n",check_qsum);	*/

   // compute quartet likelihood
   quart_lik = StoreQuarts[nn]->spprobs[0]*log(q[0]) + StoreQuarts[nn]->spprobs[1]*log(q[1]) + StoreQuarts[nn]->spprobs[2]*log(q[2]) + StoreQuarts[nn]->spprobs[3]*log(q[3]) + StoreQuarts[nn]->spprobs[4]*log(q[4]) 
                    + StoreQuarts[nn]->spprobs[5]*log(q[5]) + StoreQuarts[nn]->spprobs[6]*log(q[6]) + StoreQuarts[nn]->spprobs[7]*log(q[7]) + StoreQuarts[nn]->spprobs[8]*log(q[8]) 
                    + StoreQuarts[nn]->spprobs[9]*log(q[9]) + StoreQuarts[nn]->spprobs[10]*log(q[10]) + StoreQuarts[nn]->spprobs[11]*log(q[11]) + StoreQuarts[nn]->spprobs[12]*log(q[12]) + StoreQuarts[nn]->spprobs[13]*log(q[13]) + StoreQuarts[nn]->spprobs[14]*log(q[14]);
   
  //if (verbose==1) printf("The likelihood of the quartet with symmetry %d is %f\n\n",StoreQuarts[nn]->ncherries,quart_lik);
  return(quart_lik);

}  


double GetCompLik_genetree() {

  int i, j, k, l, m;
  int ovec[5], duppvec[15];
  double complik = 0.0;

  num_unique_quarts = 0;
  for (i=0; i<nquarts+1; i++) qvec[i]=0;

  for (i=1; i<ntaxa+1; i++){
    for (j=i+1; j<ntaxa+1; j++) {
      for (k=j+1; k<ntaxa+1; k++) {
        for (l=k+1; l<ntaxa+1; l++) {
          //printf("Quartet %d %d %d %d\n",i,j,k,l);
	  //printf("\t %llu %d %d\n",subset_to_index(i-1,j-1,k-1,l-1),(int)subset_to_index(i-1,j-1,k-1,l-1),qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]);
	  //printf("\t Have we looked at this quartet yet? %d\n",qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]); 
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

		GetQuartetTree(i,j,k,l,ovec);
  		//printf("ovec is: %d %d %d %d %d\n",ovec[0],ovec[1],ovec[2],ovec[3],ovec[4]);
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

  complik=0.0;
  for (i=0; i<num_unique_quarts; i++) {

	//printf("%d ",i);
 	//for (j=0; j<15; j++) printf("%d ",StoreQuarts[i]->spprobs[j]);
	//printf("%f %f %f ",*StoreQuarts[i]->t1,*StoreQuarts[i]->t2,*StoreQuarts[i]->t3);
	//printf("%d\n",StoreQuarts[i]->ncherries);

	complik+=QuartetLikelihood_genetree(i);

  }


  //printf("The composite likelihood is %f; this required %d quartet comps\n",complik,num_unique_quarts);
  return(complik);

}


