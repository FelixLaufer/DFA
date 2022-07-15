/** Direct Filter Approach
* Adapted from Marc Wildi, 2013: https://blog.zhaw.ch/sef/files/eco2_script.pdf
* Implementation is similar to the original R code provided here: 
* https://blog.zhaw.ch/sef/2010/10/05/analytical-dfa-a-review-of-code-versions-examples-and-parameters/
**/

#ifndef _DFA_H_
#define _DFA_H_

#include "EigenTypes.h"

class DFA
{
public:
  struct FIR
  {
    Vector coefficients;
    VectorC transferFunction;
  };

  static void periodogram(const Vector& data, Vector& per, VectorC& dft, const bool useFft = false);
  static Vector periodogram(const Vector& data, const bool useFft = false);

  /** DFA computation
  * filterLength: desired length of the computed filter
  * cutoff1, cutoff2: cutoff frequencies [0, PI] of the desired target bandpass (cutoff1 = 0 implies a lowpass design)
  * lag: lag = 0 implies real-time filtering, lag = filterLength/2 implies symmetric filtering
  * lambda: emphasizes phase artifacts in the customized criterion / timeliness
  * eta: emphasizes noise-suppression / smoothness
  * i1, i2: allow for filter constraints in frequency zero:
  *         i1 = true implies A(0) = Gamma(0) (important for trend extraction)
  *         i2 = true implies phi(0) = 0 (might be important for critical point extraction)
  **/
  static FIR designBP(const Vector& weightFunction, const unsigned int filterLength, const ScalarType cutoff1, const ScalarType cutoff2, const int lag = 0, const ScalarType lambda = 0, const ScalarType eta = 0, const bool i1 = false, const bool i2 = false);
  static FIR designLP(const Vector& weightFunction, const unsigned int filterLength, const ScalarType cutoff, const int lag = 0, const ScalarType lambda = 0, const ScalarType eta = 0, const bool i1 = false, const bool i2 = false);
  static FIR designHP(const Vector& weightFunction, const unsigned int filterLength, const ScalarType cutoff, const int lag = 0, const ScalarType lambda = 0, const ScalarType eta = 0, const bool i1 = false, const bool i2 = false);

  static void ft(const Vector& in, VectorC& out);
  static VectorC ft(const Vector& in);
  static void ift(const VectorC& in, Vector& out);
  static Vector ift(const VectorC& in);
  static void fft(const Vector& in, VectorC& out);
  static VectorC fft(const Vector& in);
  static void ifft(const VectorC& in, Vector& out);
  static Vector ifft(const VectorC& in);
  static void apply(const Vector& data, const Vector& filter, Vector& signal);
  static Vector apply(const Vector& data, const Vector& filter);
};

#endif
