extern double besselK_c(double nu, double x);

/*** Compute the S and A matrices to compute site pattern probs ***/

void ComputeAandS_popvar(double alpha, double lambda) {

   int i;
   double m = 4.0/3.0;
   double constant = (tgamma(1.0 + alpha) * (12.0 * (lambda * lambda * lambda * lambda * lambda) + lambda * m * (40.0 * (lambda * lambda * lambda) * (1.0 + alpha) + lambda * m * (51.0 * lambda + 31.0 * (3.0 + alpha) * m) * ((2.0 + alpha) * (1.0 + alpha)) + 9.0 * (m * m * m) * ((4.0 + alpha) * (3.0 + alpha) * (2.0 + alpha) * (1.0 + alpha))) + (m * m * m * m * m) * ((5.0 + alpha) * (4.0 + alpha) * (3.0 + alpha) * (2.0 + alpha) * (1.0 + alpha))));
   double sE1 = (lambda * (12.0 * pow(lambda, 4) * tgamma(1.0 + alpha) + lambda * m * (28.0 * pow(lambda, 2) + (2.0 + alpha) * m * (23.0 * lambda + 8.0 * (3.0 + alpha) * m)) * tgamma(2.0 + alpha) + pow(m, 4) * tgamma(5.0 + alpha))) / constant;

   double sE3 = ((lambda * lambda) * (12.0 * (lambda * lambda * lambda) * tgamma(1.0 + alpha) + lambda * m * (16.0 * lambda + 7.0 * (2.0 + alpha) * m) * tgamma(2.0 + alpha) + (m * m * m) * tgamma(4.0 + alpha))) / constant;   

    double sE5 = ((lambda * lambda) * (6.0 * (lambda * lambda * lambda) * tgamma(1.0 + alpha) + lambda * m * (11.0 * lambda + 6.0 * (2.0 + alpha) * m) * tgamma(2.0 + alpha) + (m * m * m) * tgamma(4.0 + alpha)))/ constant;

    double sE7 = ((lambda * lambda * lambda) * (3.0 * (lambda * lambda) + (1.0 + alpha) * m * (4.0 * lambda + (2.0 + alpha) * m)) * tgamma(1.0 + alpha))/ constant;

   //compute the C matrix in the symmetric case - indexing starts at 1 to match R code
        // c0 terms
        for (i=1; i<10; i++) smat[1][i] = 1.0/256;

        // c1 terms
        smat[2][1]=smat[2][3]=smat[2][5]=smat[2][7] = 3.0/256 * sE1;
        smat[2][2]=smat[2][4]=smat[2][6]=smat[2][8]=smat[2][9]= -1.0/256*sE1;

        // c2 terms
        smat[3][1]=smat[3][2]=smat[3][5]=smat[3][6]=3.0/(256)*sE1;
        smat[3][3]=smat[3][4]=smat[3][7]=smat[3][8]=smat[3][9]=-1.0/(256)*sE1;

        // c3 terms
        smat[4][1]=smat[4][5] = 9.0/(256)*sE3;
        smat[4][2]=smat[4][3]=smat[4][6]=smat[4][7] = -3.0/(256)*sE3;
        smat[4][4]=smat[4][8]=smat[4][9] = 1.0/(256)*sE3;

        // c4 terms
        smat[5][1]=12.0/(256)*sE1;
        smat[5][2]=smat[5][3]=smat[5][4]=4.0/(256)*sE1;
        smat[5][5]=smat[5][6]=smat[5][7]=smat[5][9]=-4.0/(256)*sE1;
        smat[5][8]=0.0;

        // c5 terms
        smat[6][1]=24.0/(256)*sE5;
        smat[6][2]=smat[6][4]=smat[6][5]=smat[6][7]=-8.0/(256)*sE5;
        smat[6][3]=8.0/(256)*sE5;
        smat[6][6]=smat[6][9]=8.0/(256)*sE5;
        smat[6][8]=0.0;

        // c6 terms
        smat[7][1]=24.0/(256)*sE5;
        smat[7][2]=8.0/(256)*sE5;
        smat[7][3]=smat[7][4]=smat[7][5]=smat[7][6]=-8.0/(256)*sE5;
        smat[7][7]=smat[7][9]=8.0/(256)*sE5;
        smat[7][8]=0.0;

        // c7 terms 
        smat[8][1]=48.0/(256)*sE7;
        smat[8][2]=smat[8][3]=smat[8][5]=-16.0/(256)*sE7;
        smat[8][4]=smat[8][6]=smat[8][7]=16.0/(256)*sE7;
        smat[8][8]=0.0;
        smat[8][9]=-16.0/(256)*sE7;

        // c8 terms
        for (i=1; i<10; i++) smat[9][i] = 0.0;

   //compute the C matrix in the asymmetric case - indexing starts at 1 to match R code
    double aE1 = (lambda * tgamma(1.0 + alpha) * (12.0 * (lambda * lambda * lambda * lambda) + lambda * m * (28.0 * (lambda * lambda) + (2.0 + alpha) * m * (23.0 * lambda + 8.0 * (3.0 + alpha) * m)) * (1.0 + alpha) + (m * m * m * m) * ((4.0 + alpha) * (3.0 + alpha) * (2.0 + alpha) * (1.0 + alpha)))) / constant;

    double aE3 = ((lambda * lambda) * tgamma(1.0 + alpha) * (6.0 * (lambda * lambda * lambda) + lambda * m * (11.0 * lambda + 6.0 * (2.0 + alpha) * m) * (1.0 + alpha) + (m * m * m) * ((3.0 + alpha) * (2.0 + alpha) * (1.0 + alpha)))) / constant;

    double aE6 = (pow(lambda, 2) * (12.0 * pow(lambda, 3) * tgamma(1.0 + alpha) + lambda * m * (16.0 * lambda + 7.0 * (2.0 + alpha) * m) * tgamma(2.0 + alpha) + pow(m, 3) * tgamma(4.0 + alpha))) / constant;

    double aE8 = ((lambda * lambda * lambda) * (3.0 * (lambda * lambda) + (1.0 + alpha) * m * (4.0 * lambda + (2.0 + alpha) * m)) * tgamma(1.0 + alpha))/ constant;
        for(int j = 1; j <= 11; j++) { amat[1][j] = 1.0 / 256.0; }

    // c1 terms            
    amat[2][1] = amat[2][3] = amat[2][4] = amat[2][5] = amat[2][8] = 3.0 / 256.0 * aE1;
    amat[2][2] = amat[2][6] = amat[2][7] = amat[2][9] = amat[2][10] = amat[2][11] = -1.0 / 256.0 * aE1;

    // c2 terms            
    amat[3][1] = amat[3][4] = 6.0 / 256.0 * aE1;
    amat[3][2] = amat[3][6] = amat[3][10] = 2.0 / 256.0 * aE1;
    amat[3][3] = amat[3][5] = amat[3][7] = amat[3][8] = amat[3][9] = amat[3][11] = -2.0 / 256.0 * aE1;

    // c3 terms 
    amat[4][1] = amat[4][4] = 12.0 / 256.0 * aE3;
    amat[4][2] = amat[4][3] = amat[4][5] = amat[4][6] = amat[4][8] = amat[4][10] = -4.0 / 256.0 * aE3;
    amat[4][7] = amat[4][9] = amat[4][11] = 4.0 / 256.0 * aE3;

    // c4 terms 
    amat[5][1] = 9.0 / 256.0 * aE1;
    amat[5][2] = amat[5][3] = 5.0 / 256.0 * aE1;
    amat[5][4] = amat[5][8] = amat[5][10] = amat[5][11] = -3.0 / 256.0 * aE1;
    amat[5][5] = amat[5][6] = amat[5][7] = amat[5][9] = 1.0 / 256.0 * aE1;

    // c5 terms 
    amat[6][1] = amat[6][3] = 12.0 / 256.0 * aE3;
    amat[6][2] = amat[6][4] = amat[6][5] = amat[6][6] = amat[6][8] = amat[6][9] = -4.0 / 256.0 * aE3;
    amat[6][7] = amat[6][10] = amat[6][11] = 4.0 / 256.0 * aE3;

    // c6 terms 
    amat[7][1] = amat[7][5] = 9.0 / 256.0 * aE6;
    amat[7][2] = amat[7][3] = amat[7][4] = amat[7][7] = amat[7][8] = -3.0 / 256.0 * aE6;
    amat[7][6] = amat[7][9] = amat[7][10] = amat[7][11] = 1.0 / 256.0 * aE6;

    // c7 terms 
    amat[8][1] = 24.0 / 256.0 * aE3;
    amat[8][2] = 8.0 / 256.0 * aE3;
    amat[8][3] = amat[8][4] = amat[8][5] = amat[8][6] = amat[8][7] = -8.0 / 256.0 * aE3;
    amat[8][8] = amat[8][11] = 8.0 / 256.0 * aE3;
    amat[8][9] = amat[8][10] = 0.0;

    // c8 terms 
    amat[9][1] = 48.0 / 256.0 * aE8;
    amat[9][2] = amat[9][3] = amat[9][4] = amat[9][5] = -16.0 / 256.0 * aE8;
    amat[9][6] = amat[9][7] = 16.0 / 256.0 * aE8;
    amat[9][8] = 16.0 / 256.0 * aE8;
    amat[9][9] = amat[9][10] = 0.0;
    amat[9][11] = -16.0 / 256.0 * aE8;
    // c9 terms 
    for (i=1; i<12; i++) amat[10][i] = 0.0;
}

double SymmetricQuartetLikelihood_popvar(int nn){

  int i, sump=0;
  double m = 4.0/3.0;  
  double sbeta[10][1];
  double q[12],check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3;

  // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;
  //printf("Times are %f %f %f\n",t1,t2,t3);

  sbeta[1][0]=1.0;
  sbeta[2][0]=exp(-2*m*t1);
  sbeta[3][0]=exp(-2*m*t2);
  sbeta[4][0]=exp(-2*m*t1)*exp(-2*m*t2);
  sbeta[5][0]=exp(-2*m*t3);
  sbeta[6][0]=exp(-m*t1)*exp(-2*m*t3);
  sbeta[7][0]=exp(-m*t2)*exp(-2*m*t3);
  sbeta[8][0]=exp(-m*t1)*exp(-m*t2)*exp(-2*m*t3);
  sbeta[9][0]=exp(-4*t3*m);

  // Need to get correct entries in smat
  int kappa = 4;

    double branch_lengths = 2 * t3 - t1 - t2;
    double z = 2 * sqrt(2 * branch_lengths * lambda);

    double constant = (tgamma(1.0 + alpha) * (12.0 * (lambda * lambda * lambda * lambda * lambda) + lambda * m * (40.0 * (lambda * lambda * lambda) * (1.0 + alpha) + lambda * m * (51.0 * lambda + 31.0 * (3.0 + alpha) * m) * ((2.0 + alpha) * (1.0 + alpha)) + 9.0 * (m * m * m) * ((4.0 + alpha) * (3.0 + alpha) * (2.0 + alpha) * (1.0 + alpha))) + (m * m * m * m * m) * ((5.0 + alpha) * (4.0 + alpha) * (3.0 + alpha) * (2.0 + alpha) * (1.0 + alpha))));
   
    double sE81 = (pow(2.0, 1.0 + alpha / 2.0) * pow(lambda, 3.0) * pow((-t1 - t2 + 2.0 * t3) * lambda, alpha / 2.0) * m * (2.0 * (t1 + t2 - 2.0 * t3) * lambda * ((-pow(kappa, 2.0)) * lambda - 2.0 * (t1 + t2 - 2.0 * t3) * pow(m, 2.0) + 2.0 * (t1 + t2 - 2.0 * t3) * kappa * pow(m, 2.0)) * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))) - sqrt(2.0) * sqrt(-t1 - t2 + 2.0 * t3) * pow(kappa, 2.0) * pow(lambda, 1.5) * (3.0 + alpha + 2.0 * t1 * m + 2.0 * t2 * m - 4.0 * t3 * m) * besselK_c(-3.0 - alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))))) / constant;          

    double sE84 =(pow(2.0, 2.0 + alpha / 2.0) * pow(lambda, 3.0) * pow((-t1 - t2 + 2.0 * t3) * lambda, alpha / 2.0) * m * (2.0 * (t1 + t2 - 2.0 * t3) * lambda * (2.0 * (t1 + t2 - 2.0 * t3) * pow(m, 2.0) - 2.0 * (t1 + t2 - 2.0 * t3) * kappa * pow(m, 2.0) + pow(kappa, 2.0) * (-lambda + (t1 + t2 - 2.0 * t3) * pow(m, 2.0))) * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))) - sqrt(2.0) * sqrt(-t1 - t2 + 2.0 * t3) * kappa * pow(lambda, 1.5) * (-2.0 * (t1 + t2 - 2.0 * t3) * m + kappa * (3.0 + alpha + 3.0 * t1 * m + 3.0 * t2 * m - 6.0 * t3 * m)) * besselK_c(-3.0 - alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))))) / constant;  
             
    double sE85 = (pow(2.0, 1.0 + alpha / 2.0) * pow(lambda, 3.0) * pow((-t1 - t2 + 2.0 * t3) * lambda, alpha / 2.0) * m * (2.0 * (t1 + t2 - 2.0 * t3) * lambda * ((-pow(kappa, 2.0)) * lambda + 2.0 * (t1 + t2 - 2.0 * t3) * pow(m, 2.0)) * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))) - sqrt(2.0) * sqrt(-t1 - t2 + 2.0 * t3) * kappa * pow(lambda, 1.5) * ((3.0 + alpha) * kappa + 4.0 * (t1 + t2 - 2.0 * t3) * m) * besselK_c(-3.0 - alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))))) / constant;   
              
    double sE86 = (pow(2.0, 2.0 + alpha / 2.0) * -1.0 * pow(-t1 - t2 + 2.0 * t3, (2.0 + alpha) / 2.0) * pow(lambda, 4.0 + alpha / 2.0) * pow(m, 2.0) * (2.0 * (t1 + t2 - 2.0 * t3) * m * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))) - sqrt(2.0) * kappa * sqrt(-((t1 + t2 - 2.0 * t3) * lambda)) * besselK_c(-3.0 - alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))))) /  constant;

    double sE88 = (pow(2.0, 2.0 + alpha / 2.0) * pow(lambda, 6.0) * pow((-t1 - t2 + 2.0 * t3) * lambda, alpha / 2.0) * pow(m, 2.0) * (sqrt(2.0) * kappa * pow((-t1 - t2 + 2.0 * t3) / lambda, 1.5) * besselK_c(3.0 + alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))) + (2.0 * pow(t1 + t2 - 2.0 * t3, 2.0) * (-2.0 + kappa) * m * besselK_c(4.0 + alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda)))) / pow(lambda, 2.0)))/ constant;    
               
    double sE89 = (pow(2.0, 3.0 + alpha / 2.0) * pow(-t1 - t2 + 2.0 * t3, (4.0 + alpha) / 2.0) * pow(lambda, 4.0 + alpha / 2.0) * pow(m, 3.0) * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt(-((t1 + t2 - 2.0 * t3) * lambda))))/ constant;                
    
     // c8 terms
        smat[9][1]=6.0*sE81/256;
        smat[9][2]=smat[9][3]=-2.0*sE81/256;

        smat[9][4]=sE84/256;
        smat[9][5]=2.0*sE85/256;
        smat[9][6]=smat[9][7]=2.0*sE86/256;
        smat[9][8]=-sE88/256;
        smat[9][9]=2.0*sE89/256;
  

  // order is: 1 = xxxx; 2 = xxxy = xxyx; 3 = xyxx = yxxx; 4 = xyxy = yxxy
  // 5 = xxyy; 6 = xxyz; 7 = yzxx; 8 = xyxz = yxxz = xyzx = yxzx; 9 = xyzw
  // to get site patterns probs, multiply transpose(smat)*sbeta
  // weighted version is below

  for (i=1; i<12; i++) q[i]=0.0;
    for (i=1; i<10; i++) {
            q[1] += 4*smat[i][1]*sbeta[i][0];
            q[2] += 12*smat[i][2]*sbeta[i][0];//24
            q[3] += 12*smat[i][3]*sbeta[i][0];//24
            q[4] += 12*smat[i][4]*sbeta[i][0];//24
            q[5] += 12*smat[i][5]*sbeta[i][0];
            q[6] += 24*smat[i][6]*sbeta[i][0];
            q[7] += 24*smat[i][7]*sbeta[i][0];
            q[8] += 24*smat[i][8]*sbeta[i][0];//96
            q[9] += 24*smat[i][9]*sbeta[i][0];      
  }

   // compute quartet likelihood
   quart_lik = StoreQuarts[nn]->spprobs[0]*log(q[1]) + StoreQuarts[nn]->spprobs[1]*log(q[2]) + StoreQuarts[nn]->spprobs[2]*log(q[2]) + StoreQuarts[nn]->spprobs[3]*log(q[3]) + StoreQuarts[nn]->spprobs[4]*log(q[3]) 
                    + StoreQuarts[nn]->spprobs[5]*log(q[4]) + StoreQuarts[nn]->spprobs[6]*log(q[4]) + StoreQuarts[nn]->spprobs[7]*log(q[5]) + StoreQuarts[nn]->spprobs[12]*log(q[6]) 
                    + StoreQuarts[nn]->spprobs[13]*log(q[7]) + StoreQuarts[nn]->spprobs[8]*log(q[8]) + StoreQuarts[nn]->spprobs[10]*log(q[8]) + StoreQuarts[nn]->spprobs[9]*log(q[8]) + StoreQuarts[nn]->spprobs[11]*log(q[8]) + StoreQuarts[nn]->spprobs[14]*log(q[9]);
   
  //if (verbose==1) printf("The likelihood of the symmetric quartet is %f\n\n",quart_lik);
  return(quart_lik);

}  


double AsymmetricQuartetLikelihood_popvar(int nn){                            

  int i, sump=0;
  double m = 4.0/3.0; // t is 2*theta
  double abeta[11][1];
  double q[12],check_qsum, prsum;
  double quart_lik = 0.0;
  double t1, t2, t3;


  // t1 = TimeVecQuart[ntaxa+3]; t2 = TimeVecQuart[ntaxa+2]; t3 = TimeVecQuart[ntaxa+1]
  t1 = *StoreQuarts[nn]->t1;
  t2 = *StoreQuarts[nn]->t2;
  t3 = *StoreQuarts[nn]->t3;

  abeta[1][0]=1.0;
  abeta[2][0]=exp(-2.0*m*t1);
  abeta[3][0]=exp(-2.0*m*t2);
  abeta[4][0]=exp(-m*t1)*exp(-2.0*m*t2);
  abeta[5][0]=exp(-2.0*m*t3);
  abeta[6][0]=exp(-m*t1)*exp(-2.0*m*t3);
  abeta[7][0]=exp(-2.0*m*t1)*exp(-2.0*m*t3);
  abeta[8][0]=exp(-m*t2)*exp(-2.0*m*t3);
  abeta[9][0]=exp(-m*t1)*exp(-m*t2)*exp(-2.0*m*t3);
  abeta[10][0]=exp(-2*m*(t2+t3));

  // Updated amat
  int kappa = 4;

  double constant = (tgamma(1.0 + alpha) * (12.0 * (lambda * lambda * lambda * lambda * lambda) + lambda * m * (40.0 * (lambda * lambda * lambda) * (1.0 + alpha) + lambda * m * (51.0 * lambda + 31.0 * (3.0 + alpha) * m) * ((2.0 + alpha) * (1.0 + alpha)) + 9.0 * (m * m * m) * ((4.0 + alpha) * (3.0 + alpha) * (2.0 + alpha) * (1.0 + alpha))) + (m * m * m * m * m) * ((5.0 + alpha) * (4.0 + alpha) * (3.0 + alpha) * (2.0 + alpha) * (1.0 + alpha))));
  
  double bessel_arg = 2 * sqrt(2) * sqrt((-t1 + t2) * lambda);

    double K2 = besselK_c(2, bessel_arg);
    double K3 = besselK_c(3, bessel_arg);
    double K4 = besselK_c(4, bessel_arg);

    double aE91 = (pow(2.0, 1.0 + alpha / 2.0) * (t1 - t2) * pow(lambda, 4.0) * pow((-t1 + t2) * lambda, alpha / 2.0) * m * (-2.0 * (pow(kappa, 2.0) * lambda + 2.0 * (t1 - t2) * pow(m, 2.0) + 2.0 * (-t1 + t2) * kappa * pow(m, 2.0)) * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda)) + (sqrt(2.0) * pow(kappa, 2.0) * lambda * (3.0 + alpha + 2.0 * t1 * m - 2.0 * t2 * m) * besselK_c(-3.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda))) / sqrt((-t1 + t2) * lambda))) / constant;

    double aE95 = (-1.0 * pow(2.0, 1.0 + alpha / 2.0) * pow(-t1 + t2, (1.0 + alpha) / 2.0) * pow(lambda, (7.0 + alpha) / 2.0) * m * (2.0 * sqrt((-t1 + t2) * lambda) * ((-pow(kappa, 2.0)) * lambda + 2.0 * (t1 - t2) * pow(m, 2.0)) * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda)) + sqrt(2.0) * kappa * lambda * ((3.0 + alpha) * kappa + 4.0 * (t1 - t2) * m) * besselK_c(-3.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda)))) / constant;

    double aE96 = (pow(2.0, 2.0 + alpha / 2.0) * (t1 - t2) * pow(lambda, 4.0) * pow((-t1 + t2) * lambda, alpha / 2.0) * m * (-2.0 * (2.0 * (-t1 + t2) * pow(m, 2.0) + 2.0 * (t1 - t2) * kappa * pow(m, 2.0) + pow(kappa, 2.0) * (lambda + (-t1 + t2) * pow(m, 2.0))) * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda)) + (sqrt(2.0) * kappa * lambda * (2.0 * (-t1 + t2) * m + kappa * (3.0 + alpha + 3.0 * t1 * m - 3.0 * t2 * m)) * besselK_c(-3.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda))) / sqrt((-t1 + t2) * lambda))) / constant;

    double aE97 = (pow(2.0, 2.0 + alpha / 2.0) * (t1 - t2) * pow(lambda, 4.0) * pow((-t1 + t2) * lambda, alpha / 2.0) * pow(m, 2.0) * (2.0 * (t1 - t2) * m * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda)) - sqrt(2.0) * kappa * sqrt((-t1 + t2) * lambda) * besselK_c(-3.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda)))) / constant;

    double aE99 = (pow(2.0, 2.0 + alpha / 2.0) * (t1 - t2) * pow(lambda, 4.0) * pow((-t1 + t2) * lambda, alpha / 2.0) * pow(m, 2.0) * (2.0 * (t1 - t2) * (-2.0 + kappa) * m * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda)) - sqrt(2.0) * kappa * sqrt((-t1 + t2) * lambda) * besselK_c(-3.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda)))) / constant;

    double aE911 = (pow(2.0, 3.0 + alpha / 2.0) * pow((-t1 + t2) / lambda, (4.0 + alpha) / 2.0) * pow(lambda, 6.0 + alpha) * pow(m, 3.0) * besselK_c(-4.0 - alpha, 2.0 * sqrt(2.0) * sqrt((-t1 + t2) * lambda))) / constant;
    
    // c9 terms 
    amat[10][1] = 6.0 / 256.0 * aE91;
    amat[10][2] = amat[10][3] = amat[10][4] = -2.0 / 256.0 * aE91;
    amat[10][5] = 2.0 / 256.0 * aE95;
    amat[10][6] = 1.0 / 256.0 * aE96;
    amat[10][7] = amat[10][8] = 2.0 / 256.0 * aE97;
    amat[10][9] = amat[10][10] = -1.0 / 256.0 * aE99;
    amat[10][11] = 2.0 / 256.0 * aE911;                   


  // order is: xxxx; xxxy = xxyx; xyxx; yxxx; xxyy; xyxy = yxxy; xxyz
  // yzxx; xyxz = xyzx; yxxz = yxzx; xyzw
  // to get site patterns probs, multiply amat*abeta (not transpose) -- changed to transpose
  // weighted version is below

  for (i=1; i<12; i++) q[i]=0.0;
   for (i=1; i<11; i++) {
            q[1] += 4*amat[i][1]*abeta[i][0]; //4
            q[2] += 12*amat[i][2]*abeta[i][0]; //24
            q[3] += 12*amat[i][3]*abeta[i][0]; //12
            q[4] += 12*amat[i][4]*abeta[i][0]; //12
            q[5] += 12*amat[i][5]*abeta[i][0]; //12
            q[6] += 12*amat[i][6]*abeta[i][0]; //24
            q[7] += 24*amat[i][7]*abeta[i][0]; //24
            q[8] += 24*amat[i][8]*abeta[i][0]; //24
            q[9] += 24*amat[i][9]*abeta[i][0];  //48
            q[10] += 24*amat[i][10]*abeta[i][0]; //48
            q[11] += 24*amat[i][11]*abeta[i][0]; //24
   }

  // compute likelihood
  quart_lik = StoreQuarts[nn]->spprobs[0]*log(q[1]) + StoreQuarts[nn]->spprobs[1]*log(q[2]) + StoreQuarts[nn]->spprobs[2]*log(q[2]) + StoreQuarts[nn]->spprobs[3]*log(q[3]) + StoreQuarts[nn]->spprobs[4]*log(q[4])
                    + StoreQuarts[nn]->spprobs[7]*log(q[5]) + StoreQuarts[nn]->spprobs[5]*log(q[6]) + StoreQuarts[nn]->spprobs[6]*log(q[6]) + StoreQuarts[nn]->spprobs[12]*log(q[7]) + StoreQuarts[nn]->spprobs[13]*log(q[8])
                    + StoreQuarts[nn]->spprobs[8]*log(q[9]) + StoreQuarts[nn]->spprobs[9]*log(q[9]) + StoreQuarts[nn]->spprobs[10]*log(q[10]) + StoreQuarts[nn]->spprobs[11]*log(q[10]) + StoreQuarts[nn]->spprobs[14]*log(q[11]);

  //if (verbose==1) printf("The likelihood of the quartet is %f\n\n",quart_lik);
  return(quart_lik);

}


double GetCompLik_popvar() {

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
          //printf("Quartet %d %d %d %d\n",i,j,k,l);
	  //printf("\t %llu %d %d\n",subset_to_index(i-1,j-1,k-1,l-1),(int)subset_to_index(i-1,j-1,k-1,l-1),qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]);
	  //printf("\t Have we looked at this quartet yet? %d\n",qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]); 
          if (qvec[(int)subset_to_index(i-1,j-1,k-1,l-1)]==0) {

		GetQuartetTree(i,j,k,l,ovec);
  		//printf("ovec is: %d %d %d %d %d\n",ovec[0],ovec[1],ovec[2],ovec[3],ovec[4]);
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

  for (i=0; i<num_unique_quarts; i++) {

	//printf("%d ",i);
 	//for (j=0; j<15; j++) printf("%d ",StoreQuarts[i]->spprobs[j]);
	//printf("%f %f %f ",*StoreQuarts[i]->t1,*StoreQuarts[i]->t2,*StoreQuarts[i]->t3);
	//printf("%d\n",StoreQuarts[i]->ncherries);

	if (StoreQuarts[i]->ncherries==2) complik+=SymmetricQuartetLikelihood_popvar(i);
        else if (StoreQuarts[i]->ncherries==1) complik+=AsymmetricQuartetLikelihood_popvar(i);
        else { printf("There was an error forming quartets .... exiting.\n"); exit(1);}
	//printf("\t %f\n",complik);

  }


  //printf("The composite likelihood is %f; this required %d quartet comps\n",complik,num_unique_quarts);
  return(complik);

}


