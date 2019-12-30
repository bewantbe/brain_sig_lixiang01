// 
// mex -lfftw3 stft_wnd_c.cpp
// CXXFLAGS='-O2 -DNDEBUG' mkoctfile stft_wnd_c.cpp
// For test
// CXXFLAGS='-O2 -march=native -fopenmp -std=c++14' LDFLAGS='-march=native -fopenmp' mkoctfile --mex stft_wnd_c.cpp
// CXXFLAGS='-O2 -fno-omit-frame-pointer -DNDEBUG' mkoctfile --mex stft_wnd_c.cpp

#include <numeric>  // For std::itoa
#include <vector>
#include <fftw3.h>
#include <eigen3/Eigen/Dense>

// fftw-wisdom -v -o wisdom rof32 rof64 rof128 rof256 rof512 rof1024 rof2048 rof4096 rof8192 rof16384 rof65536

template <typename Derived>
void stft_complex(
  Eigen::MatrixBase<Derived> &xtf_r,
  Eigen::MatrixBase<Derived> &xtf_i, 
  const Eigen::MatrixXd &x,
  const Eigen::VectorXd &wnd,
  size_t sz_hop, size_t sz_fft)
{
  size_t n_var = x.cols();
  size_t n_hop = xtf_r.cols();
  size_t sz_dat = wnd.rows();
  Eigen::MatrixXd u(sz_dat, n_var);

  fftw_complex *in, *out;
  fftw_plan p;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sz_fft);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sz_fft);  
  p = fftw_plan_dft_1d(sz_fft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  for (int l = 0; l < n_hop; l++) {
    //u = x.block(l*sz_hop, 0, sz_dat, n_var).array().colwise() * wnd.array();
    u = x.block(l*sz_hop, 0, sz_dat, n_var).array() * wnd.array();
    for (int j = 0; j < sz_dat; j++) {
      in[j][0] = u(j,0);
      in[j][1] = 0;
    }
    fftw_execute(p);
    double *pr = xtf_r.col(l).data();
    double *pi = xtf_i.col(l).data();
    for (int j = 0; j < sz_fft; j++) {
      pr[j] = out[j][0];
      pi[j] = out[j][1];
    }
  }
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}

template <typename Derived>
void stft_real(
  Eigen::MatrixBase<Derived> &xtf_r,
  Eigen::MatrixBase<Derived> &xtf_i, 
  const Eigen::MatrixXd &x,
  const Eigen::VectorXd &wnd,
  size_t sz_hop, size_t sz_fft)
{
  size_t n_var = x.cols();
  size_t n_hop = xtf_r.cols();
  size_t sz_dat = wnd.rows();
  Eigen::MatrixXd u(sz_dat, n_var);

  double *in;
  fftw_complex *out;
  fftw_plan p;
  in  = (double*) fftw_malloc(sizeof(double) * sz_fft);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sz_fft);  
  p = fftw_plan_dft_r2c_1d(sz_fft, in, out, FFTW_ESTIMATE);

  for (int l = 0; l < n_hop; l++) {
    //u = x.block(l*sz_hop, 0, sz_dat, n_var).array().colwise() * wnd.array();
    u = x.block(l*sz_hop, 0, sz_dat, n_var).array() * wnd.array();
    memcpy(in, u.data(), sizeof(double) * sz_dat);
    fftw_execute(p);
    double *pr = xtf_r.col(l).data();
    double *pi = xtf_i.col(l).data();
    pr[0] = out[0][0];
    pi[0] = out[0][1];
    for (int j = 1; j <= sz_fft/2; j++) {
      pr[j] = out[j][0];
      pi[j] = out[j][1];
      pr[sz_fft-j] = out[j][0];
      pi[sz_fft-j] =-out[j][1];
    }
  }

  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}

#ifdef MATLAB_MEX_FILE
#include <mex.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs < 2 || nrhs > 4) {
    mexErrMsgTxt("stft_wnd(x, wnd [, sz_hop [, sz_fft]])");
  }

  size_t sz_dat = mxGetNumberOfElements(prhs[1]);
  size_t sz_fft = sz_dat;
  size_t sz_hop = sz_dat/2;

  // read parameters
  if (nrhs >= 3 && !mxIsEmpty(prhs[2])) {
    sz_hop = (size_t)mxGetScalar(prhs[2]);
  }
  if (nrhs >= 4 && !mxIsEmpty(prhs[3])) {
    sz_fft = (size_t)mxGetScalar(prhs[3]);
  }

  // wnd should be 1-dim
  const Eigen::Map<const Eigen::VectorXd, Eigen::Aligned>
    wnd(mxGetPr(prhs[1]), sz_dat);
  
  // no zero padding
  // Eigen defaults to column-major
  const Eigen::Map<const Eigen::MatrixXd, Eigen::Aligned>
    x(mxGetPr(prhs[0]), mxGetM(prhs[0]), mxGetN(prhs[0]));
  
  size_t n_hop = (x.rows() - sz_dat) / sz_hop + 1;
  
  // Parameter sane test
  if (x.cols() > 1 || wnd.cols() > 1) {
    mexErrMsgTxt("x and wnd must be column vector.\n");
  }
  
  // output parameter
  plhs[0] = mxCreateDoubleMatrix(sz_fft, n_hop, mxCOMPLEX);
  
  Eigen::Map<Eigen::MatrixXd, Eigen::Aligned>
    xtf_r(mxGetPr(plhs[0]), mxGetM(plhs[0]), mxGetN(plhs[0]));
  Eigen::Map<Eigen::MatrixXd, Eigen::Aligned>
    xtf_i(mxGetPi(plhs[0]), mxGetM(plhs[0]), mxGetN(plhs[0]));

  int ret = fftw_import_wisdom_from_filename("wisdom");
  
  stft_real(xtf_r, xtf_i, x, wnd, sz_hop, sz_fft);
//  stft_complex(xtf_r, xtf_i, x, wnd, sz_hop, sz_fft);

  if (nlhs >= 2) {
    plhs[1] = mxCreateDoubleMatrix(1, n_hop, mxREAL);
    double *tc = mxGetPr(plhs[1]);
    for (int k = 0; k < n_hop; k++)
      tc[k] = sz_hop * k + sz_dat/2.0;
  }

  if (nlhs >= 3) {
    plhs[2] = mxCreateDoubleMatrix(1, sz_fft, mxREAL);
    double *fc = mxGetPr(plhs[2]);
    for (int k = 0; k < (sz_fft+1)/2; k++)
      fc[k] = (double)k / sz_fft;
    for (int k = (sz_fft+1)/2; k < sz_fft; k++)
      fc[k] = (double)(k - sz_fft) / sz_fft;
  }
}

#else

#include <octave/oct.h>

void stft_real_o(
    ComplexMatrix &xtf,
    const Matrix &x,
    const Matrix &wnd,
    size_t sz_hop, size_t sz_fft)
{
//  size_t n_var = x.cols();
  size_t n_hop = xtf.cols();
  size_t sz_dat = wnd.rows();
  const double *wp = wnd.fortran_vec();

  double *in = (double*) fftw_malloc(sizeof(double) * sz_fft);
  fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sz_fft);
  fftw_plan p = fftw_plan_dft_r2c_1d(sz_fft, in, out, FFTW_ESTIMATE);

  for (size_t l = 0; l < n_hop; l++) {
    const double *xp = x.fortran_vec() + l*sz_hop;
    for (size_t j = 0; j < sz_dat; j++) {
      in[j] = xp[j] * wp[j];
    }
    fftw_execute(p);
    std::complex<double> *p = xtf.fortran_vec() + l*sz_fft;
    p[0].real(out[0][0]);
    p[0].imag(out[0][1]);
    for (size_t j = 1; j <= sz_fft/2; j++) {
      p[j].real(out[j][0]);
      p[j].imag(out[j][1]);
      p[sz_fft-j].real( out[j][0]);
      p[sz_fft-j].imag(-out[j][1]);
    }
  }
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}

#define HELP_STR "stft_wnd_c(x, wnd [, sz_hop [, sz_fft]])"

DEFUN_DLD (stft_wnd_c, args, nargout, HELP_STR)
{
  if (args.length() < 2 || args.length() > 4) {
    print_usage ();
    return octave_value_list ();
  }
  // skip error check
  
  Matrix x   = args(0).matrix_value();
  Matrix wnd = args(1).matrix_value();

  // Parameter sane test
  if (x.cols() > 1 || wnd.cols() > 1) {
    printf("x and wnd must be column vector.\n");
    return octave_value_list ();
  }
  
  size_t sz_dat = wnd.rows();
  size_t sz_hop = sz_dat / 2;
  size_t sz_fft = sz_dat;

  // read parameters
  if (args.length() >= 3 && !args(2).is_empty()) {
    sz_hop = args(2).uint64_value();
  }
  if (args.length() >= 4 && !args(3).is_empty()) {
    sz_fft = args(3).uint64_value();
  }
  
  size_t n_hop = (x.rows() - sz_dat) / sz_hop + 1;

  // output argument
  ComplexMatrix xtf(sz_fft, n_hop);

  fftw_import_wisdom_from_filename("wisdom");

  // STFT  
  stft_real_o(xtf, x, wnd, sz_hop, sz_fft);

  // central time of each spectrum frames
  Matrix tc(1, n_hop);
  std::iota(tc.fortran_vec(), tc.fortran_vec() + n_hop, 0);
  tc = tc * sz_hop + sz_dat/2.0;

  // frequency points
  Matrix fqs(1, sz_fft);
  std::iota(fqs.fortran_vec(), fqs.fortran_vec() + (sz_fft/2+1), 0);
  std::iota(fqs.fortran_vec() + (sz_fft/2+1), fqs.fortran_vec() + sz_fft, -floor((sz_fft-1)/2.0));
  fqs /= sz_fft;

  return ovl(xtf, tc, fqs);
}

#endif

// vim: set expandtab shiftwidth=2 softtabstop=2:
