#include <stdlib.h>
#include "../dr_admm_slp.hpp"
#include "gtest/gtest.h"


#ifdef DOUBLE
#define FTYPE double
#define EPS 1e-8
#endif

#ifdef SINGLE
#define FTYPE float
#define EPS 1e-5
#endif

class SLPTest : public ::testing::Test{
protected:

  SLPTest(){
  }

  virtual ~SLPTest(){
  }

  virtual void SetUp(){
  }

  virtual void TearDown(){
  }
};

TEST_F(SLPTest, Levinson) {

  int n = 5;
  FTYPE mu[]   = {1.0, 2, -0.5, 1.2, -0.2};
  FTYPE b[] = {6.2, 6.7, -1.2, 5.2, 0.2};
  FTYPE x[] = {5.510242403550700, -0.917036531239330, -0.556094230112666, 1.669170365312393, -1.213895527483783};

  FTYPE xh[] = {0.0, 0.0, 0.0, 0.0, 0.0};
  
  levinson(n, mu, b, xh);

  for( int i = 0 ; i < n ; i++)
    ASSERT_NEAR(x[i], xh[i], EPS) << "Vectors x and xh differ at index " << i;

}


TEST_F(SLPTest, FastLevinson) {

  int n = 5;
  FTYPE mu[]   = {1.0, 2, -0.5, 1.2, -0.2};
  FTYPE b[] = {6.2, 6.7, -1.2, 5.2, 0.2};
  FTYPE x[] = {5.510242403550700, -0.917036531239330, -0.556094230112666, 1.669170365312393, -1.213895527483783};
  FTYPE r[] = {-0.669421487603306,  0.559779614325069, -0.293112947658402, -0.330578512396694};

  FTYPE xh[] = {0.0, 0.0, 0.0, 0.0, 0.0};
  
  flevinson * l = new flevinson(n);
  l->update_mu(mu);
  l->solve(b, xh);

  for( int i = 0 ; i < n-1 ; i++ )
    ASSERT_NEAR(r[i], l->r[i], EPS) << "Vectors x and xh differ at index " << i;

  for( int i = 0 ; i < n ; i++)
    ASSERT_NEAR(x[i], xh[i], EPS) << "Vectors x and xh differ at index " << i;

  delete l;
}


TEST_F(SLPTest, S) {

  int n = 5;
  FTYPE x[]   = {1.0, 2, -0.5, 1.2, -0.2};
  FTYPE t = 0.6;
  FTYPE y[] = {0.4, 1.4, 0.0, 0.6, 0.0};
  FTYPE yh[] = {0.0, 0.0, 0.0, 0.0, 0.0};
  
  S(x, n, t, yh);

  for( int i = 0 ; i < n ; i++)
    ASSERT_NEAR(y[i], yh[i], EPS) << "Vectors y and yh differ at index " << i;

}

TEST_F(SLPTest, autocorr) {

  int n = 5;
  int k = 10;
  FTYPE x[]   =  {1.0, 2.0, -0.5, 1.2, -0.2, 0.1, 0.6, 0.8, 0.1, -0.5};
  FTYPE r[] = {8.0, 0.71, 1.74, 1.02, 0.59, 1.12};
  FTYPE rh1[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  FTYPE rh2[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  autocorr(x, k, rh1, n);

  for( int i = 0 ; i < n ; i++)
    ASSERT_NEAR(r[i], rh1[i], EPS) << "Vectors r and rh1 differ at index " << i;

  fastautocorr * f = new fastautocorr(k);
  f->autocorr(x, k, rh2, n);

  for( int i = 0 ; i < n ; i++)
    ASSERT_NEAR(r[i], rh2[i], EPS) << "Vectors r and rh2 differ at index " << i;
  
  delete f;
}

TEST_F(SLPTest, filter) {

  int n = 3;
  int k = 10;
  FTYPE x[]   =  {1.0, 2.0, -0.5, 1.2, -0.2, 0.1, 0.6, 0.8, 0.1, -0.5};
  FTYPE a[] = {2.0, 4.0, -1.0, 1.0};
  FTYPE y[] = {2.0, 8.0, 6.0, -0.6, 6.9, -2.3, 3.0, 3.7, 2.9, -0.8};
  FTYPE yh[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  filter(a, n, x, k, yh);

  for( int i = 0 ; i < k ; i++)
    ASSERT_NEAR(y[i], yh[i], EPS) << "Vectors y and yh differ at index " << i;

}

TEST_F(SLPTest, fftfilter) {

  int n = 3;
  int k = 10;
  FTYPE x[]   =  {1.0, 2.0, -0.5, 1.2, -0.2, 0.1, 0.6, 0.8, 0.1, -0.5};
  FTYPE a[] = {2.0, 4.0, -1.0, 1.0};
  FTYPE y[] = {2.0, 8.0, 6.0, -0.6, 6.9, -2.3, 3.0, 3.7, 2.9, -0.8};
  FTYPE yh[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  fftfilter * h = new fftfilter(n, x, k);

  h->filter(a, yh);
  for( int i = 0 ; i < k ; i++)
    ASSERT_NEAR(y[i], yh[i], EPS) << "Vectors y and yh differ at index " << i;

  h->set_signal(x);

  h->filter(a, yh);
  for( int i = 0 ; i < k ; i++)
    ASSERT_NEAR(y[i], yh[i], EPS) << "Vectors y and yh differ at index " << i;

  delete h;
}

TEST_F(SLPTest, adjointfilter) {

  int l = 4;
  int k = 10;
  FTYPE x[] = {1.0, 2.0, -0.5, 1.2, -0.2, 0.1, 0.6, 0.8, 0.1, -0.5};
  FTYPE y[] = {2.0, 8.0, 6.0, -0.6, 6.9, -2.3, 3.0, 3.7, 2.9, -0.8};
  FTYPE a[] = {18.119999999999997, 33.800000000000004, -0.540000000000000,  17.019999999999996};
  FTYPE ah[] = {0.0, 0.0, 0.0, 0.0};

  adjointfilter(x, y, k, ah, l);

  for( int i = 0 ; i < l ; i++ )
    ASSERT_NEAR(a[i], ah[i], EPS) << "Vectors r and rh differ at index " << i;

}

TEST_F(SLPTest, fftadjointfilter) {

  int l = 4;
  int k = 10;
  FTYPE x[] = {1.0, 2.0, -0.5, 1.2, -0.2, 0.1, 0.6, 0.8, 0.1, -0.5};
  FTYPE y[] = {2.0, 8.0, 6.0, -0.6, 6.9, -2.3, 3.0, 3.7, 2.9, -0.8};
  FTYPE a[] = {18.119999999999997, 33.800000000000004, -0.540000000000000,  17.019999999999996};
  FTYPE ah[] = {0.0, 0.0, 0.0, 0.0};

  fftadjointfilter * h = new fftadjointfilter(l, x, k);

  h->set_signal(x);
  h->filter(y, ah);

  for( int i = 0 ; i < l ; i++ )
    ASSERT_NEAR(a[i], ah[i], EPS) << "Vectors r and rh differ at index " << i;

  h->set_signal(x);
  h->filter(y, ah);
  for( int i = 0 ; i < l ; i++ )
    ASSERT_NEAR(a[i], ah[i], EPS) << "Vectors r and rh differ at index " << i;

  delete h;
}

TEST_F(SLPTest, PQ) {
  
  int m = 2;
  int n = 4;

  FTYPE A[] = {3.578396939725760, 2.769437029884877, -1.349886940156521, 3.034923466331855, 0.725404224946106, -0.063054873189656, 0.714742903826096, -0.204966058299775}; //Column major
  FTYPE b[] = {-0.124144348216312, 1.489697607785465};

  leastnorm *s = new leastnorm(m, n, A, b);
  
  FTYPE x[] = {1.409034489800479,  1.417192413429614,  0.671497133608080, -1.207486922685038};
  FTYPE xh[] = {0.165334060946720, 0.275704771196147, 0.665918539710267, -1.156590843664913};

  s->PQ(x);

  for( int i = 0 ; i < n ; i++)
    ASSERT_NEAR(x[i], xh[i], EPS) << "Vectors x and xh differ at index " << i;

  delete s;
}


TEST_F(SLPTest, DRLEASTNORM) {

  int n = 4;
  int m = 2;

  FTYPE A[] = {3.578396939725760, 2.769437029884877, -1.349886940156521, 3.034923466331855, 0.725404224946106, -0.063054873189656, 0.714742903826096, -0.204966058299775}; //Column major
  FTYPE b[] = {-0.124144348216312, 1.489697607785465};

  leastnorm *s = new leastnorm(m, n, A, b);
  s->set_epsilon(1e-14);

  FTYPE z[] = {0.0, 0.0, 0.0, 0.0};
  FTYPE xs[] = {0.111939235559929,  0.388704676348588,  0.000000000001635,  0.000000000001131};

  dr(s, z);

  for( int i = 0 ; i < n ; i++ )
    ASSERT_NEAR(s->u[i], xs[i], 1e-6) << "Vectors x and xs differ at index " << i;

  delete s;
}

TEST_F(SLPTest, PQSLP) {
  
  int m = 6;
  int n = 3;

  FTYPE x[] = {0, 0.537667139546100,  1.833885014595086, -2.258846861003648, 0, 0};


  slp_dr * s = new slp_dr(m, n, 0, x, 1.0);

  FTYPE z[] = {0.862173320368121,  0.318765239858981, -1.307688296305273, -0.433592022305684, 0.342624466538650, 3.578396939725760, 2.769437029884877, -1.349886940156521, 3.034923466331855};

  FTYPE zh[] = {0.395245850887881, 1.009133017685083, -0.562228058594951, 0,  0.212510706064351,  1.267413106064431,  0.655542517225683, -3.310538560775258, 1.269987085325380};

  s->PQ(z); //in-place projection

  for( int i = 0 ; i < m+n ; i++)
    ASSERT_NEAR(z[i], zh[i], EPS) << "Vectors z and zh differ at index " << i;

  delete s;
}


TEST_F(SLPTest, PQADMMSLP) {
  
  int m = 6;
  int n = 3;

  FTYPE x[] = {0, 0.537667139546100,  1.833885014595086, -2.258846861003648, 0, 0};


  slp_admm *s = new slp_admm(m, n, 0, x, 2.0);

  FTYPE z[] = {0.862173320368121,  0.318765239858981, -1.307688296305273, -0.433592022305684, 0.342624466538650, 3.578396939725760, 2.769437029884877, -1.349886940156521, 3.034923466331855};

  FTYPE zh[] = { -0.707964004654025, -1.763775978426263, 0.315651764053395, 0.537667139546100, 2.024209505237052, -1.135522378939491, 0.732832242625563, -2.281484436154914, 0.356504498201137};

  s->PQ(z); //in-place projection

  for( int i = 0 ; i < m+n ; i++)
    ASSERT_NEAR(z[i], zh[i], EPS) << "Vectors z and zh differ at index " << i;

  delete s;
}


TEST_F(SLPTest, PQSLPFAST) {
  
  int m = 6;
  int n = 3;

  FTYPE x[] = {0, 0.537667139546100,  1.833885014595086, -2.258846861003648, 0, 0};


  slp_dr * s = new slp_dr(m, n, 1, x, 1.0);

  FTYPE z[] = {0.862173320368121,  0.318765239858981, -1.307688296305273, -0.433592022305684, 0.342624466538650, 3.578396939725760, 2.769437029884877, -1.349886940156521, 3.034923466331855};

  FTYPE zh[] = {0.395245850887881, 1.009133017685083, -0.562228058594951, 0,  0.212510706064351,  1.267413106064431,  0.655542517225683, -3.310538560775258, 1.269987085325380};

  s->PQ(z); //in-place projection

  for( int i = 0 ; i < m+n ; i++)
    ASSERT_NEAR(z[i], zh[i], EPS) << "Vectors z and zh differ at index " << i;

  delete s;
}

TEST_F(SLPTest, PQADMMSLPFAST) {
  
  int m = 6;
  int n = 3;

  FTYPE x[] = {0, 0.537667139546100,  1.833885014595086, -2.258846861003648, 0, 0};


  slp_admm *s = new slp_admm(m, n, 1, x, 2.0);

  FTYPE z[] = {0.862173320368121,  0.318765239858981, -1.307688296305273, -0.433592022305684, 0.342624466538650, 3.578396939725760, 2.769437029884877, -1.349886940156521, 3.034923466331855};

  FTYPE zh[] = { -0.707964004654025, -1.763775978426263, 0.315651764053395, 0.537667139546100, 2.024209505237052, -1.135522378939491, 0.732832242625563, -2.281484436154914, 0.356504498201137};

  s->PQ(z); //in-place projection

  for( int i = 0 ; i < m+n ; i++)
    ASSERT_NEAR(z[i], zh[i], EPS) << "Vectors z and zh differ at index " << i;

  delete s;
}

/* Matlab code to construct a result
m = 20;
n = 5;

s = filter([1, 0.8, -0.7, 0.5], 1, randn(m-n,1));
X1 = [0; s; zeros(n-1, 1)];
x = [s; zeros(n, 1)];
X = toeplitz(X1, zeros(1, n));
gamma = 1.0;

cvx_begin
cvx_quiet(true)
variables alpha(n)
minimize norm(x-X*alpha, 1)+gamma*norm(alpha, 1)
cvx_end
*/

TEST_F(SLPTest, DRSLP) {

  int m = 20;
  int n = 5;
  FTYPE gamma = 1.0;
  FTYPE X1[] = {0,  -0.532011376808821,   1.256494493216122,   
                0.842361493336700,  -2.627876731706770,   0.355005750528082,
                -1.843010108681678,  -0.875113723440593,   0.038084513542038,
                0.858277827027098,   1.070762671208747,  -2.471519081963466,
                1.691440144780230,   3.137882003152485,  -0.896409948579252,
                -1.752407618513786, 0, 0, 0, 0};

  slp_dr * s = new slp_dr(m, n, 0, X1, gamma);
  s->set_epsilon(1e-14);

  FTYPE z[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  FTYPE alpha[] = {0.000000000000900,  -0.529968471732926,  -0.000000000024857,  -0.000000000055337,  -0.000000000029813};

  dr(s, z);

  for( int i = 0 ; i < n ; i++ )
    ASSERT_NEAR(s->u[i], alpha[i], 1e-3) << "Vectors x and xs differ at index " << i;

  delete s;
}


TEST_F(SLPTest, ADMMSLP) {

  int m = 20;
  int n = 5;
  FTYPE gamma = 1.0;
  FTYPE X1[] = {0,  -0.532011376808821,   1.256494493216122,   
                0.842361493336700,  -2.627876731706770,   0.355005750528082,  
                -1.843010108681678,  -0.875113723440593,   0.038084513542038,
                0.858277827027098,   1.070762671208747,  -2.471519081963466,
                1.691440144780230,   3.137882003152485,  -0.896409948579252,
                -1.752407618513786, 0, 0, 0, 0};

  slp_admm * s = new slp_admm(m, n, 0, X1, gamma);

  FTYPE y[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  FTYPE u[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  FTYPE alphah[] = {0.000000000000900,  -0.529968471732926,  -0.000000000024857,  -0.000000000055337,  -0.000000000029813};

  s->set_epsilon(1e-14);
  s->set_kmax(4000);

  admm(s, y, u);
  
  for( int i = 0 ; i < n ; i++ )
    ASSERT_NEAR(alphah[i], (1/gamma)*y[i], 1e-4) << "Vectors x and xs differ at index " << i;

  delete s;
}


TEST_F(SLPTest, DRSLPFAST) {

  int m = 20;
  int n = 5;
  FTYPE gamma = 1.0;
  FTYPE X1[] = {0,  -0.532011376808821,   1.256494493216122,   0.842361493336700,  -2.627876731706770,   0.355005750528082,  -1.843010108681678,  -0.875113723440593,   0.038084513542038,   0.858277827027098,   1.070762671208747,  -2.471519081963466,   1.691440144780230,   3.137882003152485,  -0.896409948579252,  -1.752407618513786, 0, 0, 0, 0};

  slp_dr *s = new slp_dr(m, n, 1, X1, gamma);
  s->set_epsilon(1e-14);

  slp_dr *ss = new slp_dr(m, n, 0, X1, gamma);
  ss->set_epsilon(1e-14);

  FTYPE z[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  FTYPE zz[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  FTYPE alpha[] = {0.000000000000900,  -0.529968471732926,  -0.000000000024857,  -0.000000000055337,  -0.000000000029813};

  dr(s, z);
  dr(ss, zz);

  for( int i = 0 ; i < n ; i++ )
    ASSERT_NEAR(s->u[i], alpha[i], 1e-3) << "Vectors x and xs differ at index " << i;

  for( int i = 0 ; i < n ; i++ )
    ASSERT_NEAR(ss->u[i], s->u[i], EPS) << "Vectors u and us differ at index " << i;

  delete s;
  delete ss;
}

TEST_F(SLPTest, ADMMSLPFAST) {

  int m = 20;
  int n = 5;
  FTYPE gamma = 1.0;
  FTYPE X1[] = {0,  -0.532011376808821,   1.256494493216122,   0.842361493336700,  -2.627876731706770,   0.355005750528082,  -1.843010108681678,  -0.875113723440593,   0.038084513542038,   0.858277827027098,   1.070762671208747,  -2.471519081963466,   1.691440144780230,   3.137882003152485,  -0.896409948579252,  -1.752407618513786, 0, 0, 0, 0};

  slp_admm * s = new slp_admm(m, n, 1, X1, gamma);
  slp_admm * ss = new slp_admm(m, n, 0, X1, gamma);

  FTYPE y[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  FTYPE u[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  FTYPE yy[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  FTYPE uu[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  FTYPE alphah[] = {0.000000000000900,  -0.529968471732926,  -0.000000000024857,  -0.000000000055337,  -0.000000000029813};

  s->set_epsilon(1e-14);
  s->set_kmax(4000);

  ss->set_epsilon(1e-14);
  ss->set_kmax(4000);

  admm(s, y, u);
  admm(ss, yy, uu);

  for( int i = 0 ; i < n ; i++ )
    ASSERT_NEAR(alphah[i], (1/gamma)*y[i], 1e-4) << "Vectors x and xs differ at index " << i;

  for( int i = 0 ; i < n ; i++ )
    ASSERT_NEAR((1/gamma)*y[i], (1/gamma)*yy[i], EPS) << "Vectors x and xs differ at index " << i;

  delete s;
  delete ss;
}

int main(int argc, char **argv){
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
