# DFA - Direct Filtering Approach

Design customized FIR (finite impulse response) high, low and band pass filters that fit your use-case and data.
In (real-time) signal extraction there is an unavoidable trilemma between accurateness, timeliness and smootheness of a filter's response.
The direct filtering approach (DFA) by Marc Wildi allows to control the mutual interaction of those three criterions using two weighting parameters.
This is an efficient C++ implementation of the original R code of DFA provided by Marc Wildi; see also: https://www.researchgate.net/publication/319521573_Direct_Filter_Approach_DFA

## Usage
```cpp
  // low-pass filter design
  const unsigned int len = 20; // number of filter coeffs/tabs
  const ScalarType cutoff = 0.1; // normalized cutoff frequency
  const ScalarType lag = 0; // filter lag 0 implies "zero" lag, i.e. real time filtering
  const ScalarType lambda = 0; // lambda emphasizes timeliness vs. accurateness, 0 implies a the MSE solution, try to increase it!
  const ScalarType eta = 0; // eta emphasizes smoothing vs. accurateness, 0 implies a the MSE solution try to increase it!
  const Vector per = DFA::periodogram(data); // compute the periodogram, i.e. an estimate of the signal's power spectry density
  const DFA::FIR fir = DFA::designLP(per, len, cutoff, lag, lambda, eta); // design the low-pass filter

  // offline filtering
  const Vector offlineSignal = DFA::apply(data, fir.coefficients);

  // online filtering
  Vector onlineSignal = Vector::Zero(data.size());
  DFA::OnlineFIR filt = DFA::OnlineFIR(fir.coefficients);
  for (unsigned int i = 0; i < data.size(); ++i)
    onlineSignal[i] = filt.process(data[i]);
```

## Requires
- Eigen3
