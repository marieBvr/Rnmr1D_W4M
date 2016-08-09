#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

using namespace std;
using namespace Rcpp;

class data_info {
public:
    double  pmin;
    double  pmax;
    int size_c;
    int size_l;
};

double _abs(double x)
{
   if (x<0.0) x=-x;
   return x;
}

int indMaxC(SEXP x, int start, int end)
{
   NumericVector v(x);
   int n1 = (start>0) ? start : 0;
   int n2 = (end<v.size()) ? end : v.size();
   int indx = n1;
   for(int i = n1; i <= n2; i++) if (v[i]>v[indx]) indx=i;
   return indx;
}

int indMinC(SEXP x, int start, int end)
{
   NumericVector v(x);
   int n1 = (start>0) ? start : 0;
   int n2 = (end<v.size()) ? end : v.size();
   int indx = n1;
   for(int i = n1; i <= n2; i++) if (v[i]<v[indx]) indx=i;
   return indx;
}

double maxC(SEXP x, int start, int end)
{
   NumericVector v(x);
   return v[indMaxC(x, start, end)];
}
double minC(SEXP x, int start, int end)
{
   NumericVector v(x);
   return v[indMinC(x, start, end)];
}

// [[Rcpp::export]]
SEXP SDL(SEXP x, double Sigma)
{
   NumericVector X(x);
   int N = X.size();
   NumericVector Out(N);
   double v1, v2;
   #pragma omp parallel for schedule(static) 
   for(int n = 0; n<N; n++)
   {
       v1 = X[n]*X[n]; v2 = Sigma*Sigma; 
       Out[n] = (12.0*v1-v2)/pow(4.0*v1+v2, 3);
   }
   return Out;
}

// ---------------------------------------------------
//  Read / Write the Matrix of spectra wihtin a binary file
// ---------------------------------------------------

// [[Rcpp::export]]
void C_write_pack (SEXP x, double pmin, double pmax, SEXP ff)
{
   // Matrix of spectra : 1 row = 1 spectrum, 1 column = a same value of ppm
   NumericMatrix xx(x);

   std::string fname = as<std::string>(ff); 

   // Header table
   data_info* inforec = new data_info();
   inforec->pmax=pmax;
   inforec->pmin=pmin;
   inforec->size_l=xx.ncol();
   inforec->size_c=xx.nrow()+2;

   // Open the file as binary mode
   std::ofstream outBinFile;
   outBinFile.open(fname.c_str(), std::ios::out | std::ios::binary);

   // write the header table
   outBinFile.write( (char *)inforec, sizeof(data_info) );

   // Buffer
   double* buf = new double[inforec->size_c];

   // write data by column
   buf[0] = buf[xx.nrow()+1] = 0.0;
   for(int i = 0; i<xx.ncol(); i++) {
      for(int k = 1; k<=xx.nrow(); k++) buf[k] = xx(k-1, i);
      outBinFile.write( (char *)buf, (unsigned)(inforec->size_c)*sizeof(double) );
   }

   outBinFile.flush();
   outBinFile.close(); 
}

// [[Rcpp::export]]
SEXP C_read_pack (SEXP ff)
{
   string fname = as<string>(ff); 

   // Header table
   data_info* inforec = new data_info();
   
   // Open the file as binary mode
   ifstream inBinFile;
   inBinFile.open(fname.c_str(), ios::in | ios::binary);
   inBinFile.seekg (0, ios::beg);
   
   // read the header table
   inBinFile.read( (char *)inforec, sizeof(data_info) );
   double pmax = inforec->pmax;
   double pmin = inforec->pmin;
   int ncol = (unsigned)(inforec->size_l);
   int nrow = (unsigned)(inforec->size_c)-2;

   // Matrix of spectra : 1 row = 1 spectrum, 1 column = a same value of ppm
   NumericMatrix M(nrow, ncol);

   // Buffer
   double* buf = new double[inforec->size_c];

   // Read data by column
   buf[0] = buf[nrow+1] = 0.0;
   for(int i = 0; i<ncol; i++) {
       inBinFile.read( (char *)buf, (unsigned)(inforec->size_c)*sizeof(double) );
       for(int k = 1; k<=nrow; k++) M(k-1, i) = buf[k];
   }

   inBinFile.close(); 

   return Rcpp::List::create(_["int"] = M,
                             _["nspec"] = nrow,
                             _["size"] = ncol,
                             _["ppm_min"] = pmin,
                             _["ppm_max"] = pmax );

}

// ---------------------------------------------------
//  Major Peaks finder
// ---------------------------------------------------

// [[Rcpp::export]]
SEXP C_findPeaks(SEXP x, SEXP y, double sig, double DW, double minRatio, double minInt0, double minInt, bool bPlasma)
{
   NumericVector specR(x), V(y);
   List pkList = List::create();
   int TD = specR.size();
   int WN_2 = (int)(0.02/DW+0.5);
   int OFFSET = (int)(0.05*TD+0.5);
   int dNmin = (int)(0.01/DW+0.5);
   int dNmax = (int)(5/DW+0.5);
   
   int iprev = OFFSET;
   bool fRange = false;
   for(int i = (OFFSET+1); i<(TD-OFFSET); i++)
   {
       if ( ( maxC(V, i-WN_2, i+WN_2)-minC(V, i-WN_2, i+WN_2) )< sig ) {
          if (fRange) {
              fRange=false;
              int dN=i-iprev+1;
              if (dN<dNmin || dN>dNmax) continue;
              int N = indMaxC(specR,iprev,i);
              double maxInt=specR[N];
              double Nratio= (double)(N)/(double)(TD);
              if (bPlasma) {
                  if ( maxInt<minInt ) continue;
              } else {
                  if ( Nratio<minRatio && maxInt<minInt0 ) continue;
                  if ( Nratio>minRatio && maxInt<minInt ) continue;
              }
              int dn = bPlasma ? (int)(2*(i-iprev)+0.5) : (int)((i-iprev)/10+0.5);
              List elem = List::create(
                _["N"] = N+1,
                _["maxInt"] = (double)(maxInt),
                _["dppm"] = dN*DW,
                _["n1"] = iprev-dn+1,
                _["n2"] = i+dn+1 
              );
              pkList.push_back(elem);
          }
          iprev=i;
      } else {
          fRange=true;
      }
   }
   return pkList;

}

// ---------------------------------------------------
//  Phasing Routines
// ---------------------------------------------------

// [[Rcpp::export]]
SEXP C_corr_spec_re (SEXP l)
{
   List spec(l);
   NumericVector re = as<NumericVector>(spec["re"]);
   NumericVector im = as<NumericVector>(spec["im"]);
   List proc = as<List>(spec["proc"]);
   double phc0 = as<double>(proc["phc0"]);
   double phc1 = as<double>(proc["phc1"]);
   List pklisk = as<List>(spec["pklist"]);
   int I1 = as<int>(spec["I1"]) - 1;
   List pk1 = as<List>(pklisk[I1]);
   int N1 = as<int>(pk1["N"]) - 1;
   double Isign = (as<bool>(proc["NegPhi"])) ? 1.0 : -1.0 ;
   int TD = re.size();
   double phi;

   NumericVector S_re(TD);
   //NumericVector S_im(TD);
   #pragma omp parallel for schedule(static) 
   for (int i=0; i<TD; i++)
   {
       phi = ( phc0 + (N1-i)*phc1 )*M_PI/180.0;
       S_re[i] = cos(phi)*re[i] + Isign*sin(phi)*im[i];
   }
   return S_re;
}

// [[Rcpp::export]]
SEXP C_corr_spec_im (SEXP l)
{
   List spec(l);
   NumericVector re = as<NumericVector>(spec["re"]);
   NumericVector im = as<NumericVector>(spec["im"]);
   List proc = as<List>(spec["proc"]);
   List pklisk = as<List>(spec["pklist"]);
   double phc0 = as<double>(proc["phc0"]);
   double phc1 = as<double>(proc["phc1"]);
   int I1 = as<int>(spec["I1"]) - 1;
   List pk1 = as<List>(pklisk[I1]);
   int N1 = as<int>(pk1["N"]) - 1;
   double Isign = (as<bool>(proc["NegPhi"])) ? 1.0 : -1.0 ;
   int TD = re.size();
   double phi;

   NumericVector S_im(TD);
   #pragma omp parallel for schedule(static) 
   for (int i=0; i<TD; i++)
   {
       phi = ( phc0 + (N1-i)*phc1 )*M_PI/180.0;
       S_im[i] = cos(phi)*im[i] - Isign*sin(phi)*re[i];
   }
   return S_im;
}

// [[Rcpp::export]]
double C_phc0_fn (SEXP l, double phc0)
{
   List spec(l);
   NumericVector re = as<NumericVector>(spec["re"]);
   NumericVector im = as<NumericVector>(spec["im"]);
   List proc = as<List>(spec["proc"]);
   List pklist = as<List>(spec["pklist"]);
   int I1 = as<int>(spec["I1"]) - 1;
   List pk1 = as<List>(pklist[I1]);
   int n1 = as<int>(pk1["n1"]);
   int n2 = as<int>(pk1["n2"]);
   double phi = phc0*M_PI/180.0;
   double Isign = (as<bool>(proc["NegPhi"])) ? 1.0 : -1.0 ;

   double S=0;
   double V;
   for (int i=(n1-1); i<n2; i++)
   { 
       V = cos(phi)*re[i] + Isign*sin(phi)*im[i];
       if (V<0.0) S += V;
   }
   return _abs(S);
}

// [[Rcpp::export]]
double C_phc1_fn (SEXP l, double phc1)
{
   List spec(l);
   NumericVector re = as<NumericVector>(spec["re"]);
   NumericVector im = as<NumericVector>(spec["im"]);
   List proc = as<List>(spec["proc"]);
   double phc0 = as<double>(proc["phc0"]);
   List param = as<List>(spec["param"]);
   bool plasma = as<bool>(param["PLASMA"]);
   List pklisk = as<List>(spec["pklist"]);
   int I0 = as<int>(spec["I0"]) - 1;
   List pk0 = as<List>(pklisk[I0]);
   int I1 = as<int>(spec["I1"]) - 1;
   List pk1 = as<List>(pklisk[I1]);
   int N1 = as<int>(pk1["N"]) - 1;
   double Isign = (as<bool>(proc["NegPhi"])) ? 1.0 : -1.0 ;

   int n1,n2;
   double S=0.0;
   if (plasma) {
      double V,phi;
      double M=0.0;
      int np=0;
      n1 = as<int>(pk0["n1"]) - 1;
      n2 = (int)(0.5*( as<int>(pk1["n1"]) + as<int>(pk0["n2"]) ) + 0.5) - 1;
      for (int i=n1; i<=n2; i++) { 
         phi = ( phc0 + (N1-i)*phc1 )*M_PI/180.0;
         V = cos(phi)*re[i] + Isign*sin(phi)*im[i];
         if (V<0.0) { M += V; np++; }
      }
      M /= (double)(np-1);
      for (int i=n1; i<=n2; i++) { 
         phi = ( phc0 + (N1-i)*phc1 )*M_PI/180.0; V = cos(phi)*re[i] + Isign*sin(phi)*im[i];
         if (V<0.0) S += (V-M);
      }
   } else {
      double V1,V2,phi;
      n1 = as<int>(pk0["n1"]) - 1;
      n2 = as<int>(pk0["n2"]) - 1;
      phi = ( phc0 + (N1-n1)*phc1 )*M_PI/180.0; V1 = cos(phi)*re[n1] + Isign*sin(phi)*im[n1];
      phi = ( phc0 + (N1-n2)*phc1 )*M_PI/180.0; V2 = cos(phi)*re[n2] + Isign*sin(phi)*im[n2];
      S = V2 - V1;
   }
   return _abs(S);
}

// [[Rcpp::export]]
double C_optim_phc (SEXP s, SEXP p, Function fn)
{
    List spec(s);
    NumericVector phcrange(p);
    NumericVector res(1);
    double phc_1, phc_2, phcm, v1, v2, v;
    double phc_d = 0;
    double prev = 0;
    bool flg = true;
    double tol = 0.0001;
    int cnt = 0;
    int cnt_max = 100 ;
    while (flg)
    {
       phcm = 0.5*(phcrange[0]+phcrange[1]);
       phc_1 = 0.5*(phcrange[0] + phcm);
       phc_2 = 0.5*(phcrange[1] + phcm);
       res = fn(s, phc_1); v1 = (double)res[0];
       res = fn(s, phc_2); v2 = (double)res[0];
       if (_abs(v1)<_abs(v2)) {
           v=v1; phc_d=phc_1; phcrange[1] = phcm ;
       } else {
           v=v2; phc_d=phc_2; phcrange[0] = phcm ;
       }
       if (v==0) { 
           flg = false;
       } else if (_abs((phc_d-prev)/phc_d)<tol) flg = false;
       prev = phc_d;
       cnt++;
       if (cnt>cnt_max) flg = false;
    }
    return phc_d;
}

// ---------------------------------------------------
//  Baseline Correction Routines
// ---------------------------------------------------

// [[Rcpp::export]]
SEXP lowpass1 (SEXP x, double alpha)
{
   NumericVector y1(x);
   NumericVector y2(clone(x));
   int n=y1.size();
   int k;
   for (k=1; k<n; k++) y2[k] = y2[k-1] + alpha * (y1[k] - y2[k-1]);
   return(y2);
}

// [[Rcpp::export]]
double WinMoy (SEXP v, int n1, int n2)
{
    NumericVector specR(v);
    int k;
    double  moy=0.0;
    for (k=n1; k<=n2; k++) moy += specR[k];
    moy /= (double)(n2-n1+1);
    return moy;
}

// [[Rcpp::export]]
SEXP Smooth (SEXP v, int n)
{
    NumericVector specR(v);
    int n1,n2;
    int N = specR.size();
    NumericVector S(N);
    for (int count=0; count<N; count++) {
        n1 = count >= n ? count - n : 0;
        n2 = count <= N - n - 1 ? count + n : N - 1;
        S[count] = WinMoy(specR,n1,n2);
    }
    return S;
}

// [[Rcpp::export]]
void Ajust_LB (SEXP s, SEXP b, int n1, int n2)
{
    NumericVector specR(s), lb(b);
    int k,ni;
    double  a,diff,diff_max,lb_line;

    a=(lb[n2]-lb[n1])/(n2-n1);
    diff_max=0.0; ni=n1;
    for (k=n1; k<n2; k++) {
        lb_line = a*(k-n1)+lb[n1];
        diff = specR[k]< lb_line ? lb_line - specR[k] : 0.0 ;
        if (diff>diff_max) { diff_max=diff; ni=k; }
    }
    if (ni>n1 && ni<n2) {
        a=(specR[ni]-lb[n1])/(ni-n1);
        for (k=n1+1; k<=ni; k++)
            lb[k]=a*(k-n1)+lb[n1];
        Ajust_LB(specR,lb,ni,n2);
    }
    else
        for (k=n1; k<n2; k++)
            lb[k]=a*(k-n1)+lb[n1];
}

// [[Rcpp::export]]
SEXP C_Estime_LB (SEXP l, SEXP s, int istart, int iend)
{
   List spec(l);
   NumericVector specR(s);
   List proc = as<List>(spec["proc"]);
   List param = as<List>(spec["param"]);
   double sig = as<double>(spec["sig"]);
   int S1 = as<int>(param["WINDOWSIZE"]);
   int NEIGH = as<int>(param["NEIGHSIZE"]);
   double NFAC = as<double>(param["NOISELEVEL"]);
   int TD = specR.size();
   int OFFSET = 1;
   int count,n1,n2,k,cnt;
   int ws = TD>32768 ? 2 : 1 ;

   // Create the BL vector initialize with spectrum values
   NumericVector lb(TD), BL(TD), m1(TD), m2(TD);

   // (s1,neigh) = (50,35) => soft, (25,15) => intermediate, (10,5) => hard

   m1=Smooth(specR,S1*ws);
   m2=Smooth(specR,4*ws);

   cnt=0;
   for (count=0; count<TD; count++) {
        if (count<OFFSET || count>(TD-OFFSET) ) {
            lb[count]=m1[count];  n1=count; continue;
        }
        if (count<istart || count>iend ) {
            lb[count]=0;  n1=count; n2=count; cnt=0; continue;
        }
        if ((_abs(m2[count] - m1[count])) <= NFAC*sig) {
            if (cnt==0) n2=count;
            cnt++;
        }
        else {
            if (cnt<NEIGH*ws) { cnt=0; continue; }
            for (k=n2; k<count; k++) lb[k] = m1[k];
            if (n1<n2) {
                Ajust_LB(specR,lb,n1,n2);
            }
            n1=count-1;
            cnt=0;
        }
   }
   if (cnt>0) for (k=n2; k<count; k++) lb[k] = m1[k];
   return(lb);
}

// ---------------------------------------------------
//  Reference Spectrun
// ---------------------------------------------------

// [[Rcpp::export]]
SEXP spec_ref_interval (SEXP x, int istart, int iend)
{
   NumericMatrix VV(x);
   int n_specs = VV.nrow();
   int size_m = iend-istart+1;
   int count,k;

   NumericVector vref(size_m);

   for(count=0; count<size_m; count++) {
        vref[count]=0.0;
        for (k=0; k<n_specs; k++) {
            vref[count] += VV(k, istart+count);
        }
   }
   for (count=0; count<size_m; count++) vref[count] /= (double)(n_specs);

   return(vref);
}

// [[Rcpp::export]]
SEXP C_spec_ref (SEXP x)
{
   NumericMatrix VV(x);
   int count_max = VV.ncol();
   NumericVector vref = spec_ref_interval (x, 0, count_max-1);
   return(vref);
}

// [[Rcpp::export]]
SEXP C_MedianSpec(SEXP x) {
   NumericMatrix VV(x);
   int n_specs = VV.nrow();
   int count_max = VV.ncol();
   int position = n_specs / 2; // Euclidian division
   NumericVector out(count_max);
   for (int j = 0; j < count_max; j++) { 
        NumericVector y = VV(_,j); // Copy column -- original will not be mod
        std::nth_element(y.begin(), y.begin() + position, y.end()); 
        out[j] = y[position];  
   }
   return out;
}

// ---------------------------------------------------
//  Alignment Algorithms
// ---------------------------------------------------

/* Apodisation par "sigmoides symétriques" aux extremites de la zone d'alignement */
void _apodize (SEXP x, int k, int n) {
   NumericMatrix VV(x);
   double lambda=2.0;
   int N=4;
   int i;
   for (i=n-N; i<n; i++) VV(k, i) *= 1.0/(1.0+exp(-lambda*(n-N/2-i)));
   VV(k, n)=0.0;
   for (i=n+1; i<(n+N); i++) VV(k, i) *= 1.0/(1.0+exp(-lambda*(i-n-N/2)));
}

// [[Rcpp::export]]
SEXP C_segment_shifts (SEXP x, int istart, int iend, int decal_max)
{
   NumericMatrix VV(x);
   int n_specs = VV.nrow();
   int size_m = iend-istart+1;
   int i,j,k, ij, decal;
   double som, min_sse, sse;

   // NumericVector vref(size_m);
   NumericVector vk(size_m);
   NumericVector shift_v(n_specs);

   /* Spectre de reference Vref */
   NumericVector vref = spec_ref_interval (x, istart-1, iend-1);

   som=0.0;
   for (i=0; i<size_m; i++) som +=  vref[i];
   for (i=0; i<size_m; i++) vref[i] = vref[i]/som;

   /* Taille de la fenetre de glissement entre les deux massifs Vref et Vk */
   decal= (int)(size_m/3);
   if (decal_max>0 && decal>decal_max) decal = decal_max;

   /* Pour chaque spectre */
   for (k=0; k<n_specs; k++) {

       /* init */
       shift_v[k]=0;

       /* Segment du spectre Vk a aligner */
       som=0.0;
       for (i=0; i<size_m; i++) som +=  VV(k, istart+i-1);
       if (som==0.0) continue;
       for (i=0; i<size_m; i++) vk[i] = VV(k, istart+i-1)/som;

       /* Recherche le shift qui minimise la Somme des Erreurs Quadratiques (SSE) */
       min_sse=DBL_MAX;
       for (j=-decal; j<=decal; j++) {
           sse=0.0;
           for (i=0; i<size_m; i++) {
               ij=i+j;
               sse += 10.0*(vref[i] - vk[ij])*(vref[i] - vk[ij]);
           }
           if (sse<min_sse) { shift_v[k]=j; min_sse=sse; }
       }
   }

   return(shift_v);
}

// [[Rcpp::export]]
int C_align_segment (SEXP x, SEXP s, int istart, int iend)
{
   NumericMatrix VV(x);
   NumericVector shift_v(s);

   int n_specs = VV.nrow();
   int size_m = iend-istart+1;
   int i, k, ij, delta, moy_shift;

   NumericVector vk(size_m);

   /* Calcul le shift moyen */
   moy_shift = 0;
   for (k=0; k<n_specs; k++) moy_shift += shift_v[k];
   moy_shift = (int)(moy_shift/n_specs);

   /* Translate les massifs */
   for (k=0; k<n_specs; k++) {
       delta = shift_v[k]-moy_shift;
       if (delta==0) continue;
       for (i=0; i<size_m; i++) {
           vk[i]=0.0;
           ij=i+delta;
           if ( ij>=0 && ij<size_m)
               vk[i]=VV(k, istart+ij-1);
       }
       for (i=0; i<size_m; i++) VV(k,istart+i-1)=vk[i];
       _apodize(x, k,istart-1);
       _apodize(x, k,iend-1);
   }

   return(moy_shift);
}

// ---------------------------------------------------
//  Binning Algorithms
// ---------------------------------------------------

struct BinData {
   int VREF;
   int n_buckets;
   int inoise_start;
   int inoise_end;
   double delta_ppm;
   double R;
   double ynoise;
   double vnoise;
   double noise_fac;
   double bin_fac;
   double peaknoise_rate;
   double BUCMIN;
};

// [[Rcpp::export]]
double C_noise_estimation(SEXP x, int n1, int n2)
{
   NumericVector V(x);
   double  ym, sum_y, sum_y2, y_noise;
   int i;
   sum_y=sum_y2=0.0;
   for (i=n1; i<n2; i++) {
       sum_y2 += V[i]*V[i];
       sum_y  += V[i];
   }
   ym=_abs(sum_y);
   y_noise = sqrt(( sum_y2 - ym*ym/_abs(n2-n1) )/_abs(n2-n1-1));
   return y_noise;
}

/*-------- Bin Evaluation Criterion (BEC)----------------------------------*/
double bin_value(SEXP x, SEXP v, struct BinData *bdata, int n1, int n2)
{
   NumericMatrix VV(x);
   NumericVector vref(v);
   int n_specs = VV.nrow();
   int i,k;
   double Amax=0.0;
   double vb=0.0;

   if (bdata->VREF==1) {
       for(i=n1; i<=n2; i++) if (vref[i]>Amax) Amax=vref[i];
       vb += pow((Amax-vref[n1])*(Amax-vref[n2]),bdata->R);
   } else {
        for (k=0; k<n_specs; k++) {
            Amax=0.0;
            for(i=n1; i<=n2; i++) if (VV(k,i)>Amax) Amax=VV(k,i);
            vb += pow((Amax-VV(k,n1))*(Amax-VV(k,n2)),bdata->R);
        }
        vb /= n_specs;
   }
   return vb;
}

void save_bucket(SEXP b, SEXP v, struct BinData *bdata, int n1, int n2)
{
   NumericMatrix buckets(b);
   NumericVector vref(v);
   int i, flg;
   while (vref[n1]==0.0) n1++;
   while (vref[n2]==0.0) n2--;
   flg=0;
   for(i=n1; i<=n2; i++)
      if (vref[i]>bdata->peaknoise_rate*bdata->ynoise) { flg=1; break; }
   if (flg==0) return;
   if (C_noise_estimation(vref, n1, n2) < bdata->noise_fac*bdata->ynoise) return;
   if ( ( _abs(n1-n2)*bdata->delta_ppm )<bdata->BUCMIN ) return;
   if ( ( _abs(n1-n2)*bdata->delta_ppm )>1 ) return;
   Rprintf("\tsave bucket [%03d]: range %d - %d\n",bdata->n_buckets+1, n1, n2);
   buckets(bdata->n_buckets,0)=n1;
   buckets(bdata->n_buckets,1)=n2;
   bdata->n_buckets++;
}

int find_buckets(SEXP x, SEXP b, SEXP v, struct BinData *bdata, int n1, int n2)
{
   NumericMatrix buckets(b);
   NumericVector vref(v);
   int count,ncut,nmin;
   double vb1,vb2,vbsum, vbmax;
   vbmax=bin_value(x, vref,bdata,n1,n2);
   nmin=(int)(bdata->BUCMIN/bdata->delta_ppm);
   ncut=0;
   for(count=(n1+nmin); count<(n2-nmin); count++) {
        vb1=bin_value(x, vref,bdata,n1,count);
        vb2=bin_value(x, vref,bdata,count,n2);
        vbsum=vb1+vb2;
        if (vbsum>vbmax && vb1>bdata->vnoise && vb2>bdata->vnoise) {
            vbmax=vbsum;
            ncut=count;
        }
   }
   //Rprintf("find bucket: range %d - %d, vbmax=%f, nmin=%d, ncut=%d\n",n1,n2, vbmax, nmin, ncut);
   if (ncut>0) {
        if (find_buckets(x, buckets,vref,bdata,n1,ncut)==0) save_bucket(buckets,vref,bdata,n1,ncut);
        if (find_buckets(x, buckets,vref,bdata,ncut,n2)==0) save_bucket(buckets,vref,bdata,ncut,n2);
   }
   return ncut;
}

// [[Rcpp::export]]
SEXP C_aibin_buckets(SEXP x, SEXP b, SEXP v, SEXP l, int n1, int n2)
{
   NumericMatrix buckets(b);
   NumericVector vref(v);
   List blist(l);
   struct BinData bdata;
   int i, ret;

   bdata.n_buckets=0;
   bdata.VREF = as<int>(blist["VREF"]);
   bdata.R = as<double>(blist["R"]);
   bdata.noise_fac = as<double>(blist["noise_fac"]);
   bdata.bin_fac = as<double>(blist["bin_fac"]);
   bdata.peaknoise_rate = as<double>(blist["peaknoise_rate"]);
   bdata.BUCMIN = as<double>(blist["BUCMIN"]);
   bdata.delta_ppm = as<double>(blist["dppm"]);
   bdata.ynoise = as<double>(blist["ynoise"]);
   bdata.inoise_start = as<int>(blist["inoise_start"]);
   bdata.inoise_end = as<int>(blist["inoise_end"]);
   bdata.vnoise = bdata.bin_fac*bin_value(x, vref, &bdata, bdata.inoise_start, bdata.inoise_end);
   Rprintf("AIBIN: range %d - %d, vnoise=%f\n",n1,n2, bdata.vnoise);

   ret=find_buckets(x, buckets,vref,&bdata,n1-1,n2-1);
   Rprintf("Returned Value = %d, number of buckets found = %d\n",ret, bdata.n_buckets);
   if (bdata.n_buckets==0) return R_NilValue;

   NumericMatrix M(bdata.n_buckets, 2);
   for (i=0; i<bdata.n_buckets; i++) { M(i,0) = buckets(i,0)+1; M(i,1) = buckets(i,1)+1; }
   return M;
}

// [[Rcpp::export]]
SEXP C_buckets_integrate (SEXP x, SEXP b, int mode)
{
   NumericVector V(x);
   NumericMatrix Buc(b);
   int n_bucs = Buc.nrow();
   int i, m;

   //Matrix of the Buckets' integration : 1 row = 1 spectrum, 1 column = 1 bucket
   NumericVector Vbuc(n_bucs);

   // for each bucket
   for (m=0; m<n_bucs; m++) {
       // for each point
       Vbuc[m] = 0.0;
       for (i=(Buc(m,0)-1); i<(Buc(m,1)-1); i++) Vbuc[m] += 0.5*( V[i] + V[i+1] );
       if (mode == -1) Vbuc[m] /= ( Buc(m,1) - Buc(m,0) + 1 );
       if (mode ==  1) Vbuc[m] *= ( Buc(m,1) - Buc(m,0) + 1 );
   }

   return(Vbuc);
}


// [[Rcpp::export]]
SEXP C_all_buckets_integrate (SEXP x, SEXP b, int mode)
{
   NumericMatrix VV(x);
   NumericMatrix Buc(b);
   int n_specs = VV.nrow();
   int n_bucs = Buc.nrow();
   int k;

   //Matrix of the Buckets' integration : 1 row = 1 spectrum, 1 column = 1 bucket
   NumericMatrix M(n_specs, n_bucs);

   // for each spectrum
   for (k=0; k<n_specs; k++) {
        NumericVector Y = VV(k,_);
        NumericVector Z = C_buckets_integrate( Y, Buc, mode );
        M(k,_) = Z;
   }

   return(M);
}

// [[Rcpp::export]]
SEXP C_buckets_CSN_normalize (SEXP b)
{
   NumericMatrix buckets(b);
   int n_specs = buckets.nrow();
   int n_bucs = buckets.ncol();
   NumericMatrix M(n_specs, n_bucs);
   int k, m;
   double sumB;

   // for each spectrum
   for (k=0; k<n_specs; k++) {
       sumB=0.0;
       // for each bucket
       for (m=0; m<n_bucs; m++) sumB += buckets(k,m);
       for (m=0; m<n_bucs; m++) M(k,m) = 100* buckets(k,m)/sumB;
   }
   return(M);
}



