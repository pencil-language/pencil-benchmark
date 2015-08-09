void extract(int Ne, float image[const restrict static Ne]) {
#pragma scop
#pragma pencil independent
  for (int i=0; i<Ne; i++){
    image[i] = exp(image[i]/255);      /*scale to 0-1 and exponential map*/
  }
#pragma endscop
}

void compress(int Ne, float image[const restrict static Ne]) {
#pragma scop
#pragma pencil independent
  for (int i=0; i<Ne; i++) {
    image[i] = log(image[i])*255; //log-compress image, scale image up from 0-1 to 0-255
  }
#pragma endscop
}

void diffusion(int Nr, int Nc, int Ne, float q0sqr,
    int iN[const restrict static Nr], int iS[const restrict static Nr],
    int jE[const restrict static Nc], int jW[const restrict static Nc],
    float dN[const restrict static Ne], float dS[const restrict static Ne],
    float dE[const restrict static Ne], float dW[const restrict static Ne],
    float image[const restrict static Ne], float c[const restrict static Ne])
{
#pragma scop
#pragma pencil independent
  for (int j=0; j<Nc; j++) {
    for (int i=0; i<Nr; i++) {

      int    k = i + Nr*j;
      float Jc = image[k]; /*current pixel element*/

      /*compute directional derivatives of current pixel*/
      int iNi = iN[i];
      int Nrj = Nr*j;
      int iSi = iS[i];
      int jWj = jW[j];
      int jEj = jE[j];
      int NrjWj = Nr*jWj;
      int NrjEj = Nr*jEj;

      float dNk_l = image[iNi + Nrj] - Jc;
      float dSk_l = image[iSi + Nrj] - Jc;
      float dWk_l = image[i + NrjWj] - Jc;
      float dEk_l = image[i + NrjEj] - Jc;

      dN[k] = dNk_l;
      dS[k] = dSk_l;
      dW[k] = dWk_l;
      dE[k] = dEk_l;

      /*compute normalized discrete gradient magnitude squared from derivates*/
      float G2 = (dNk_l*dNk_l + dSk_l*dSk_l
          + dWk_l*dWk_l + dEk_l*dEk_l) / (Jc*Jc);

      /*compute normalized discrete laplacian*/
      float L = (dNk_l + dSk_l + dWk_l + dEk_l) / Jc;

      // ICOV (equ 31/35)
      float num  = (0.5f*G2) - ((1.0f/16.0f)*(L*L)) ;/*based on gradient and laplacian*/
      float den  = 1 + (.25f*L);									  /*based on laplacian*/
      float qsqr = num/(den*den);

      // diffusion coefficent (equ 33) (every element of IMAGE)
      den = (qsqr-q0sqr) / (q0sqr * (1+q0sqr));/*compute diffusion coefficient*/
      c[k] = 1.0f / (1.0f+den);

      if( c[k] < 0) c[k] = 0;/*clip diffusion to 0-1*/
      if( c[k] > 1) c[k] = 1;
    }
  }
#pragma endscop
}

void divergence(int Nr, int Nc, int Ne, float lambda,
    int iS[const restrict static Nr], int jE[const restrict static Nc],
    float dN[const restrict static Ne], float dS[const restrict static Ne],
    float dE[const restrict static Ne], float dW[const restrict static Ne],
    float image[const restrict static Ne], float c[const restrict static Ne])
{
#pragma scop
#pragma pencil independent
  for (int j=0; j<Nc; j++){
    for (int i=0; i<Nr; i++){
      int k = i + Nr*j;/*current pixel element position*/

      /*compute directional diffusion coefficient*/
      int Nrj = Nr*j;
      int iSi = iS[i];
      int jEj = jE[j];
      int NrjEj = Nr*jEj;

      float cN = c[k];
      float cS = c[iSi + Nrj];
      float cW = c[k];
      float cE = c[i + NrjEj];

      float D = cN*dN[k] + cS*dS[k] + cW*dW[k] + cE*dE[k];/*compute divergence*/
      image[k] = image[k] + 0.25f*lambda*D;/*update image*/
    }
  }
#pragma endscop
}

void srad(int niter, int NeROI, int Nr, int Nc, int Ne, float lambda,
    int iN[const restrict static Nr], int iS[const restrict static Nr],
    int jE[const restrict static Nc], int jW[const restrict static Nc],
    float dN[const restrict static Ne], float dS[const restrict static Ne],
    float dE[const restrict static Ne], float dW[const restrict static Ne],
    float image[const restrict static Ne], float c[const restrict static Ne])
{
  extract(Ne, image);
  for (int iter=0; iter<niter; iter++)
  {
    float sum=0;
    float sum2=0;
    for (int i = 0; i < Nr; i++) {
      for (int j = 0; j < Nc; j++) {
        float tmp = image[i + Nr*j];               /*compute std of image pixels*/
        sum  += tmp;
        sum2 += tmp*tmp;
      }
    }
    float meanROI = sum / NeROI;                      /*mean value of image in ROI*/
    float varROI  = (sum2 / NeROI) - meanROI*meanROI; /*variance*/
    float q0sqr   = varROI / (meanROI*meanROI);       /*standard deviation*/

    diffusion( Nr, Nc, Ne, q0sqr, iN, iS, jE, jW, dN, dS, dE, dW, image, c);
    divergence(Nr, Nc, Ne, lambda,    iS, jE,     dN, dS, dE, dW, image, c);
  }
  compress(Ne, image);
}
