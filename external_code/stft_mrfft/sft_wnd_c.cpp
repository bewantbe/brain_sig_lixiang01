// CXXFLAGS='-O2 -DNDEBUG' mkoctfile sft_wnd_c.cpp

#include <fftw3.h>
#include <octave/oct.h>
#include <eigen3/Eigen/Dense>

typedef std::complex<double> Complex;

ComplexNDArray sft_real(
    const Matrix &x,
    const Matrix &wnd,
    size_t sz_hop, size_t sz_fft,
    int normalize_mode)
{
  size_t x_len = x.rows();
  size_t sz_dat = wnd.rows();
  size_t n_hop = (x_len - sz_dat) / sz_hop + 1;
  size_t n_var = x.cols();
  
  // output parameters
  ComplexNDArray aveS(dim_vector(sz_fft, n_var, n_var), 0.0);
  ComplexMatrix Jk(sz_fft, n_var);

  // prepare FFT
  fftw_import_wisdom_from_filename("wisdom");
  double *in = (double*) fftw_malloc(sizeof(double) * sz_fft);
  fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sz_fft);
  fftw_plan plan = fftw_plan_dft_r2c_1d(sz_fft, in, out, FFTW_ESTIMATE);

  const double *wp = wnd.fortran_vec();

  for (int l = 0; l < n_hop; l++) {
    // windowed FFT of each data slice
    // Jk = fft(bsxfun(@times, wnd, X(bg+1:bg+sz_fft,:,i_trial)));
    for (int id_var = 0; id_var < n_var; id_var++) {
      const double *xp = x.fortran_vec() + l*sz_hop + id_var*x_len;
      for (int j = 0; j < sz_dat; j++) {
        in[j] = xp[j] * wp[j];
      }
      fftw_execute(plan);
      Complex *pjk = Jk.fortran_vec() + sz_fft*id_var;
      pjk[0].real(out[0][0]);
      pjk[0].imag(out[0][1]);
      for (int j = 1; j <= sz_fft/2; j++) {
        pjk[j].real(out[j][0]);
        pjk[j].imag(out[j][1]);
        pjk[sz_fft-j].real(out[j][0]);
        pjk[sz_fft-j].imag(-out[j][1]);
      }
    }
    // FFT to power spectrum
    // S(:, chan1, chan2) = Jk(:,chan1).*conj(Jk(:,chan2));
    for (size_t chan1 = 0; chan1 < n_var; chan1++) {
      Eigen::Map<Eigen::ArrayXcd> S(
        aveS.fortran_vec() + sz_fft*chan1 + sz_fft*n_var*chan1, sz_fft);
      const Eigen::Map<const Eigen::ArrayXcd> J1(
        Jk.fortran_vec() + sz_fft*chan1, sz_fft);
      S += J1 * J1.conjugate();
      for (size_t chan2 = chan1 + 1; chan2 < n_var; chan2++) {
        Eigen::Map<Eigen::ArrayXcd> S12(
          aveS.fortran_vec() + sz_fft*chan1 + sz_fft*n_var*chan2, sz_fft);
        const Eigen::Map<const Eigen::ArrayXcd> J2(
          Jk.fortran_vec() + sz_fft*chan2, sz_fft);
        S12 += J1 * J2.conjugate();
      }
    }
  }
  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);
  
  double qw = 0, sw = 0;
  for (int j = 0; j < sz_dat; j++) {
    qw += wp[j] * wp[j];
    sw += wp[j];
  }

  double fact = qw;  
  if (normalize_mode == 1) {
    fact = sw*sw / 4;
  }

  aveS /= n_hop * fact;

  for (size_t chan1 = 0; chan1 < n_var; chan1++) {
    for (size_t chan2 = chan1 + 1; chan2 < n_var; chan2++) {
      Complex *pS12 = aveS.fortran_vec() + sz_fft*chan1 + sz_fft*n_var*chan2;
      Complex *pS21 = aveS.fortran_vec() + sz_fft*chan2 + sz_fft*n_var*chan1;
      for (size_t j = 0; j < sz_fft; j++) {
        pS21[j] += std::conj(pS12[j]);
      }
    }
  }

  return aveS;
}

#define HELP_STR "aveS = stft_wnd_c(x, wnd [, sz_hop [, sz_fft]])"

DEFUN_DLD (sft_wnd_c, args, nargout, HELP_STR)
{
  if (args.length() < 2 || args.length() > 5) {
    print_usage ();
  }
  // skip error check
  
  Matrix x   = args(0).matrix_value();
  Matrix wnd = args(1).matrix_value();

  size_t sz_fft = wnd.rows();
  size_t sz_hop = wnd.rows()/2;
  int normalize_mode = 0;    // 0: statistics, 1: audio
  
  if (args.length() >= 3) {
    if (!args(2).is_empty()) {
      sz_hop = args(2).int64_value();
    }
    if (args.length() >= 4) {
      if (!args(3).is_empty()) {
        sz_fft = args(3).int64_value();
      }
      if (args.length() >= 5) {
        std::string spectrum_normalize = args(4).string_value();
        if (spectrum_normalize == "statistics") {
          normalize_mode = 0;
        } else if (spectrum_normalize == "audio") {
          normalize_mode = 1;
        }
      }
    }
  }

  // fqs = ifftshift((0:sz_fft-1)-floor(sz_fft/2))'/sz_fft;
  Matrix fqs(sz_fft, 1);
  std::iota(fqs.fortran_vec(), fqs.fortran_vec() + (sz_fft+1)/2, 0);
  std::iota(fqs.fortran_vec() + (sz_fft+1)/2, fqs.fortran_vec() + sz_fft, -floor(sz_fft/2.0));
  fqs /= sz_fft;
  
  return ovl(sft_real(x, wnd, sz_hop, sz_fft, normalize_mode), fqs);
}

// vim: set expandtab shiftwidth=2 softtabstop=2:
