#include "DFA.h"

#include "Constants.h"
#include <unsupported/Eigen/FFT>
#include <iostream>

void DFA::periodogram(const Vector& data, Vector& per, VectorC& dft, const bool useFft)
{
  const unsigned int len = data.size();
  const unsigned int halfLen = len / 2 + 1;
  per.resize(halfLen);
  dft.resize(halfLen);

  if (useFft)
    fft(data, dft);
  else
    ft(data, dft);

  dft = dft.segment(0, halfLen).eval();
  dft *= std::sqrt(1.0 / (2.0 * Pi * len));
  per = dft.cwiseAbs();
  per = per.cwiseProduct(per);
}

Vector DFA::periodogram(const Vector& data, const bool useFft)
{
  Vector ret;
  VectorC dft;
  periodogram(data, ret, dft, useFft);
  return ret;
}

DFA::FIR DFA::designBP(const Vector& weightFunction, const unsigned int filterLength, const ScalarType cutoff1, const ScalarType cutoff2, const int lag, const ScalarType lambda, const ScalarType eta, const bool i1, const bool i2)
{
  if (!(0 <= cutoff1 && cutoff1 <= 1 && cutoff1 < cutoff2 && cutoff2 <= 1))
  {
    std::cerr << "Invalid cutoff frequencies: please assure that 0 <= cutoff1 < cutoff2 <= 1" << std::endl;
    return FIR{ Vector::Zero(0), VectorC::Zero(0) };
  }

  const unsigned int len = weightFunction.size();
  const unsigned int K = len - 1;

  // Design target filter transfer function
  const unsigned int cutoff1Idx = static_cast<unsigned int>(std::ceil(cutoff1 * K));
  const unsigned int cutoff2Idx = static_cast<unsigned int>(std::floor(cutoff2 * K));

  Vector gamma = Vector::Zero(len);
  gamma.segment(cutoff1Idx, cutoff2Idx - cutoff1Idx + 1) = Vector::Ones(cutoff2Idx - cutoff1Idx + 1);

  Vector weightH = Vector::Zero(len);
  if (cutoff2Idx < len)
  {
    unsigned int k = 0;
    for (; k < cutoff2Idx; ++k)
      weightH[k] = weightFunction[k];
    --k;
    for (unsigned int etaBase = 1; etaBase <= len - cutoff2Idx; ++etaBase)
      weightH[k + etaBase] = weightFunction[k + etaBase] * std::pow(etaBase, eta);
  }
  else
  {
    for (unsigned int i = 0; i < len; ++i)
      weightH[i] = weightFunction[i];
  }

  // First order filter restriction: assign a large weight to frequency zero
  weightH[0] = i1 ? std::max(weightH[0], 1.0e+10) : weightH[0];

  MatrixC X, Xy;
  X.resize(len, filterLength);
  Xy.resize(len, filterLength);
  ScalarTypeC cplxFac = ScalarTypeC(0.0, -1.0) * static_cast<ScalarType>(lag * Pi);

  for (unsigned int row = 0; row <= K; ++row)
  {
    const ScalarTypeC cplxW = std::exp(cplxFac * static_cast<ScalarType>(static_cast<unsigned int>(row) / K));
    Xy(row, 0) = cplxW;
    X(row, 0) = cplxW * std::sqrt(weightH[row]);
  }

  if (i2)
  {
    // Second order restriction: time shift in frequency zero vanishes
    const unsigned int L = filterLength - 1;
    X = X.block(0, 0, X.rows(), L).eval();
    Xy = Xy.block(0, 0, Xy.rows(), L).eval();

    for (unsigned int col = 1; col < L; ++col)
    {
      const ScalarType phiFac1 = (col - lag) * Pi;
      const ScalarType phiFac2 = (L - lag) * Pi;
      for (unsigned int row = 0; row <= K; ++row)
      {
        const ScalarType phi1 = phiFac1 * row / K;
        const ScalarType phi2 = phiFac2 * row / K;
        const ScalarType wReal = std::cos(phi1) - (static_cast<ScalarType>(col) / L) * std::cos(phi2);
        const ScalarType wImag = std::sin(phi1) - (static_cast<ScalarType>(col) / L) * std::sin(phi2);
        Xy(row, col) = ScalarTypeC(wReal, wImag);
        X(row, col) = ScalarTypeC(wReal, std::sqrt(1.0 + gamma[row] * lambda) * wImag) * std::sqrt(weightH[row]);
      }
    }
  }
  else
  {
    for (unsigned int col = 1; col < filterLength; ++col)
    {
      const ScalarType phiFac = (col - lag) * Pi;
      for (unsigned int row = 0; row <= K; ++row)
      {
        const ScalarType phi = phiFac * row / K;
        const ScalarType wReal = std::cos(phi);
        const ScalarType wImag = std::sin(phi);
        Xy(row, col) = ScalarTypeC(wReal, wImag);
        X(row, col) = ScalarTypeC(wReal, std::sqrt(1.0 + gamma[row] * lambda) * wImag) * std::sqrt(weightH[row]);
      }
    }
  }

  // Solve system for filter coefficients
  Matrix XReal = X.real();
  Matrix XImag = X.imag();
  Matrix XTX = XReal.transpose() * XReal + XImag.transpose() * XImag;
  Vector coeffs = XTX.colPivHouseholderQr().solve(Xy.real().transpose() * gamma.cwiseProduct(weightH));

  // Build filter coefficients
  if (i2)
  {
    // Last weight is a function of the previous ones through the second order restriction
    Vector i2Coeffs;
    i2Coeffs.resize(filterLength);
    const unsigned int L = filterLength - 1;
    ScalarType sum = 0.0;
    for (unsigned int n = 0; n < L; ++n)
    {
      i2Coeffs[n] = coeffs[n];
      sum += n * coeffs[n];
    }
    i2Coeffs[L] = -sum / L;
    coeffs = i2Coeffs;
  }

  // Build transfer function
  const VectorC cplxCoeffs = coeffs;
  VectorC trf;
  trf.resize(len);
  trf[0] = ScalarTypeC(coeffs.sum(), 0);
  for (unsigned int k = 1; k <= K; ++k)
  {
    VectorC cplxW;
    cplxW.resize(filterLength);
    cplxFac = ScalarTypeC(0.0, 1.0) * static_cast<ScalarType>(Pi * k / K);
    for (unsigned int n = 0; n < filterLength; ++n)
      cplxW[n] = std::exp(static_cast<ScalarType>(n) * cplxFac);

    trf[k] = cplxW.dot(cplxCoeffs);
  }

  return FIR { coeffs, trf };
}

DFA::FIR DFA::designLP(const Vector& weightFunction, const unsigned filterLength, const ScalarType cutoff, const int lag, const ScalarType lambda, const ScalarType eta, const bool i1, const bool i2)
{
  return designBP(weightFunction, filterLength, 0, cutoff, lag, lambda, eta, i1, i2);
}

DFA::FIR DFA::designHP(const Vector& weightFunction, const unsigned filterLength, const ScalarType cutoff, const int lag, const ScalarType lambda, const ScalarType eta, const bool i1, const bool i2)
{
  return designBP(weightFunction, filterLength, cutoff, 1, lag, lambda, eta, i1, i2);
}

void DFA::ft(const Vector& in, VectorC& out)
{
  const unsigned int len = in.size();
  out.resize(len);

  for (unsigned int k = 0; k < len; k++)
  {
    ScalarTypeC sum(0.0, 0.0);
    for (unsigned int j = 0; j < len; j++)
      sum += in[j] * std::exp(ScalarTypeC(0.0, -2.0 * static_cast<int>(j * k) * Pi / len));
    out[k] = sum;
  }
}

VectorC DFA::ft(const Vector& in)
{
  VectorC ret;
  ft(in, ret);
  return ret;
}

void DFA::ift(const VectorC& in, Vector& out)
{
  const unsigned int len = in.size();
  out.resize(len);

  const ScalarType scaleFac = 1.0 / len;
  for (unsigned int k = 0; k < len; k++)
  {
    ScalarTypeC sum(0.0, 0.0);
    for (unsigned int j = 0; j < len; j++)
      sum += in[j] * std::exp(ScalarTypeC(0.0, 2.0 * static_cast<int>(j * k) * Pi / len));
    out[k] = sum.real() * scaleFac;
  }
}

Vector DFA::ift(const VectorC& in)
{
  Vector ret;
  ift(in, ret);
  return ret;
}

void DFA::fft(const Vector& in, VectorC& out)
{
  Eigen::FFT<ScalarType> fft;
  fft.fwd(out, in);
}

VectorC DFA::fft(const Vector& in)
{
  VectorC ret;
  fft(in, ret);
  return ret;
}

void DFA::ifft(const VectorC& in, Vector& out)
{
  Eigen::FFT<ScalarType> fft;
  fft.inv(out, in);
}

Vector DFA::ifft(const VectorC& in)
{
  Vector ret;
  ifft(in, ret);
  return ret;
}

void DFA::apply(const Vector& data, const Vector& filter, Vector& signal)
{
  const unsigned int len = data.size();
  const unsigned int firLen = filter.size();

  if (len < firLen)
  {
    signal.resize(0);
    return;
  }

  signal = Vector::Zero(len);
  for (unsigned int i = firLen - 1; i < len; i++)
  {
    ScalarType sum = 0.0;
    for (unsigned int l = 0; l < firLen; l++)
      sum += filter[l] * data[i - l];
    signal[i] = sum;
  }
}

Vector DFA::apply(const Vector& data, const Vector& filter)
{
  Vector ret;
  apply(data, filter, ret);
  return ret;
}
