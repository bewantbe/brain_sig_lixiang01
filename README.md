# Usage

Modify and run `run_freq_band.m` in Matlab(or GNU Octave) to get alpha/beta/gamma band of spikes from Intan measured "amplifier" signal.

## Possible improvement

* Replace the naive spike-sorting algorithm to a sophisticated algorithm (`spike_sorting_naive.m`).
  - Current naive algorithm: high-pass (>200Hz) then thresholding.
* Automate processing of full data set (`run_freq_band.m`, see also `file_looper.m`).
* Save resulting plot(pictures) to file.
* Relate the spike data to the behavioral data.
