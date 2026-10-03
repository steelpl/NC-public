# Executed synthetic experiment results

This report was generated from completed CSV outputs. It contains no global or ISMN results. The experiment contains **36 scenarios**, **500 Monte Carlo attempts per scenario/method**, and **72000 covariance-estimator attempts** in total. The sensitivity summary contains 30 rows with 100 realizations per setting.

All main merging results apply the same WA SPD safeguard to TC, the repository TCH variant, NC raw, and NC PSD. NC raw denotes the current optimizer before the native PSD feasibility map; it is not an unprotected matrix inversion. TCH denotes the repository penalty optimizer's first reconstructed covariance output, not its fourth fitted covariance output.

Each realization uses 1095 training samples and 1095 independently generated test samples. The recorded seed is 20260917. The configured SNR values are 0.1, 3, 10 dB, corresponding to signal/mean-error-power ratios 1.02329, 1.99526, 10. These SNR settings are sensitivity scenarios, not claims about typical operational retrieval SNR.

## Denominators and interpretation

Estimation and merge failure rates use every attempted realization. Indefiniteness rates use only finite covariance outputs (`NFiniteCovariances`), and are NA if that denominator is zero. Accuracy means are conditional on available estimates and are accompanied by valid counts. Equal-weight and known-covariance controls are not covariance estimators; covariance-estimation metrics for them are NA. The paired comparison table uses only realizations where both compared methods are valid. Confidence intervals quantify Monte Carlo uncertainty across independent realizations, not independent daily samples.

## Paired temporal comparison

Positive differences below mean larger seasonal than stationary error. Read the CI and common-valid count together; a small change does not establish universal immunity to seasonality.

| SNR dB | Base error corr. | Signal-error corr. | Method | Seasonal minus stationary test MSE [95% CI] | Common merge n |
|---|---|---|---|---|---|
| 0.1 | 0 | 0 | TC | -7.70769e-05 [-0.000161806, 7.6526e-06] | 500 |
| 0.1 | 0 | 0 | TCH | 8.42139e-06 [-2.39339e-05, 4.07767e-05] | 500 |
| 0.1 | 0 | 0 | NC_PSD | -0.000161388 [-0.000497492, 0.000174716] | 500 |
| 0.1 | 0 | 0.3 | TC | -0.00100096 [-0.00255852, 0.000556603] | 500 |
| 0.1 | 0 | 0.3 | TCH | -0.000854787 [-0.00243566, 0.000726086] | 500 |
| 0.1 | 0 | 0.3 | NC_PSD | -0.00129655 [-0.00287027, 0.00027717] | 500 |
| 0.1 | 0.4 | 0 | TC | -0.000355549 [-0.000639202, -7.18969e-05] | 500 |
| 0.1 | 0.4 | 0 | TCH | -1.97053e-05 [-0.000113064, 7.36532e-05] | 500 |
| 0.1 | 0.4 | 0 | NC_PSD | -4.12668e-05 [-0.00154877, 0.00146624] | 500 |
| 0.1 | 0.4 | 0.3 | TC | -0.00118826 [-0.0032611, 0.000884582] | 500 |
| 0.1 | 0.4 | 0.3 | TCH | -0.000690756 [-0.00285659, 0.00147508] | 500 |
| 0.1 | 0.4 | 0.3 | NC_PSD | -0.000914718 [-0.00326806, 0.00143863] | 500 |
| 0.1 | 0.99 | 0 | TC | NA [NA, NA] | 0 |
| 0.1 | 0.99 | 0 | TCH | -0.0163207 [-0.0271852, -0.00545625] | 500 |
| 0.1 | 0.99 | 0 | NC_PSD | -0.00248392 [-0.00911743, 0.00414958] | 500 |
| 0.1 | 0.99 | 0.3 | TC | NA [NA, NA] | 0 |
| 0.1 | 0.99 | 0.3 | TCH | -0.000519633 [-0.0162237, 0.0151845] | 500 |
| 0.1 | 0.99 | 0.3 | NC_PSD | -1.3651e-05 [-0.00658421, 0.0065569] | 500 |
| 3 | 0 | 0 | TC | -2.88773e-05 [-5.98819e-05, 2.12733e-06] | 500 |
| 3 | 0 | 0 | TCH | 4.52518e-06 [-8.00014e-06, 1.70505e-05] | 500 |
| 3 | 0 | 0 | NC_PSD | -8.57932e-05 [-0.000349541, 0.000177955] | 500 |
| 3 | 0 | 0.3 | TC | -0.000510655 [-0.00130959, 0.00028828] | 500 |
| 3 | 0 | 0.3 | TCH | -0.00043918 [-0.00124973, 0.000371371] | 500 |
| 3 | 0 | 0.3 | NC_PSD | -0.000643356 [-0.00148728, 0.000200572] | 500 |
| 3 | 0.4 | 0 | TC | -0.000171155 [-0.00030201, -4.02994e-05] | 500 |
| 3 | 0.4 | 0 | TCH | -9.21119e-06 [-5.05641e-05, 3.21417e-05] | 500 |
| 3 | 0.4 | 0 | NC_PSD | -0.00035447 [-0.000963043, 0.000254103] | 500 |
| 3 | 0.4 | 0.3 | TC | -0.000613741 [-0.00168029, 0.000452804] | 500 |
| 3 | 0.4 | 0.3 | TCH | -0.000367948 [-0.00147665, 0.000740753] | 500 |
| 3 | 0.4 | 0.3 | NC_PSD | -0.000487247 [-0.00164926, 0.000674765] | 500 |
| 3 | 0.99 | 0 | TC | NA [NA, NA] | 0 |
| 3 | 0.99 | 0 | TCH | -0.00111021 [-0.00234532, 0.000124897] | 500 |
| 3 | 0.99 | 0 | NC_PSD | 0.000400497 [-0.00194311, 0.0027441] | 500 |
| 3 | 0.99 | 0.3 | TC | NA [NA, NA] | 0 |
| 3 | 0.99 | 0.3 | TCH | -0.00370398 [-0.00955872, 0.00215075] | 500 |
| 3 | 0.99 | 0.3 | NC_PSD | -0.00263259 [-0.0137042, 0.00843906] | 500 |
| 10 | 0 | 0 | TC | -2.89284e-06 [-6.98167e-06, 1.19599e-06] | 500 |
| 10 | 0 | 0 | TCH | 1.09699e-06 [-6.69318e-07, 2.8633e-06] | 500 |
| 10 | 0 | 0 | NC_PSD | -1.30303e-06 [-2.76216e-05, 2.50155e-05] | 500 |
| 10 | 0 | 0.3 | TC | -0.000100743 [-0.000260215, 5.87297e-05] | 500 |
| 10 | 0 | 0.3 | TCH | -8.77338e-05 [-0.000249402, 7.39345e-05] | 500 |
| 10 | 0 | 0.3 | NC_PSD | -9.06687e-05 [-0.000252284, 7.09466e-05] | 500 |
| 10 | 0.4 | 0 | TC | -2.35247e-05 [-4.04753e-05, -6.57398e-06] | 500 |
| 10 | 0.4 | 0 | TCH | -4.71409e-07 [-6.16788e-06, 5.22506e-06] | 500 |
| 10 | 0.4 | 0 | NC_PSD | -3.28809e-05 [-0.000144833, 7.90708e-05] | 500 |
| 10 | 0.4 | 0.3 | TC | -0.000119636 [-0.000333574, 9.43021e-05] | 500 |
| 10 | 0.4 | 0.3 | TCH | -7.66861e-05 [-0.000297284, 0.000143912] | 500 |
| 10 | 0.4 | 0.3 | NC_PSD | -0.000107574 [-0.000339227, 0.000124079] | 500 |
| 10 | 0.99 | 0 | TC | NA [NA, NA] | 0 |
| 10 | 0.99 | 0 | TCH | -4.50092e-05 [-0.00018778, 9.77619e-05] | 500 |
| 10 | 0.99 | 0 | NC_PSD | 0.000891099 [0.000147441, 0.00163476] | 500 |
| 10 | 0.99 | 0.3 | TC | NA [NA, NA] | 0 |
| 10 | 0.99 | 0.3 | TCH | -3.58118e-05 [-0.000964611, 0.000892987] | 500 |
| 10 | 0.99 | 0.3 | NC_PSD | -0.000184129 [-0.000862956, 0.000494698] | 500 |

## Primary NC versus comparator on common valid realizations

These are NC PSD minus comparator differences. Negative values favor NC PSD for that metric and subset; positive values favor the comparator. They do not compare successful NC estimates with failed comparator attempts.

At the middle tested SNR (3 dB), using Seasonal truth and signal-error correlation 0:

- At base error correlation 0, `NC_PSD_minus_TC` gives a test-MSE difference of 0.00142286 [0.00117547, 0.00167025] on 500 common valid realizations.
- At base error correlation 0, `NC_PSD_minus_TCH` gives a test-MSE difference of 0.00142439 [0.00117646, 0.00167231] on 500 common valid realizations.
- At base error correlation 0.4, `NC_PSD_minus_TC` gives a test-MSE difference of 0.00595591 [0.00546064, 0.00645117] on 500 common valid realizations.
- At base error correlation 0.4, `NC_PSD_minus_TCH` gives a test-MSE difference of 0.00462542 [0.00411938, 0.00513145] on 500 common valid realizations.
- At base error correlation 0.99, `NC_PSD_minus_TC` gives a test-MSE difference of NA [NA, NA] on 0 common valid realizations.
- At base error correlation 0.99, `NC_PSD_minus_TCH` gives a test-MSE difference of -0.047986 [-0.0523312, -0.0436408] on 500 common valid realizations.

| Truth | SNR dB | Base corr. | Signal-error corr. | Contrast | Test MSE difference [95% CI] | Common n |
|---|---|---|---|---|---|---|
| seasonal_multiscale | 0.1 | 0 | 0 | NC_PSD_minus_TC | 3.58044e-05 [-0.000186002, 0.00025761] | 500 |
| seasonal_multiscale | 0.1 | 0 | 0 | NC_PSD_minus_TCH | 8.25666e-05 [-0.000130534, 0.000295667] | 500 |
| seasonal_multiscale | 0.1 | 0 | 0 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 0.1 | 0 | 0.3 | NC_PSD_minus_TC | 0.000294055 [0.000105162, 0.000482947] | 500 |
| seasonal_multiscale | 0.1 | 0 | 0.3 | NC_PSD_minus_TCH | -0.000460298 [-0.000706557, -0.000214039] | 500 |
| seasonal_multiscale | 0.1 | 0 | 0.3 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 0.1 | 0.4 | 0 | NC_PSD_minus_TC | 0.0101899 [0.00895289, 0.0114269] | 500 |
| seasonal_multiscale | 0.1 | 0.4 | 0 | NC_PSD_minus_TCH | 0.00594336 [0.0046575, 0.00722921] | 500 |
| seasonal_multiscale | 0.1 | 0.4 | 0 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 0.1 | 0.4 | 0.3 | NC_PSD_minus_TC | 0.0117678 [0.0108474, 0.0126883] | 500 |
| seasonal_multiscale | 0.1 | 0.4 | 0.3 | NC_PSD_minus_TCH | 0.000945353 [-0.000180986, 0.00207169] | 500 |
| seasonal_multiscale | 0.1 | 0.4 | 0.3 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 0.1 | 0.99 | 0 | NC_PSD_minus_TC | NA [NA, NA] | 0 |
| seasonal_multiscale | 0.1 | 0.99 | 0 | NC_PSD_minus_TCH | -0.350843 [-0.363091, -0.338595] | 500 |
| seasonal_multiscale | 0.1 | 0.99 | 0 | NC_PSD_minus_NC_raw | 0.0150687 [0.013838, 0.0162993] | 500 |
| seasonal_multiscale | 0.1 | 0.99 | 0.3 | NC_PSD_minus_TC | NA [NA, NA] | 0 |
| seasonal_multiscale | 0.1 | 0.99 | 0.3 | NC_PSD_minus_TCH | -0.551878 [-0.565661, -0.538095] | 500 |
| seasonal_multiscale | 0.1 | 0.99 | 0.3 | NC_PSD_minus_NC_raw | 0.00740331 [0.00667788, 0.00812874] | 500 |
| seasonal_multiscale | 3 | 0 | 0 | NC_PSD_minus_TC | 0.00142286 [0.00117547, 0.00167025] | 500 |
| seasonal_multiscale | 3 | 0 | 0 | NC_PSD_minus_TCH | 0.00142439 [0.00117646, 0.00167231] | 500 |
| seasonal_multiscale | 3 | 0 | 0 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 3 | 0 | 0.3 | NC_PSD_minus_TC | 0.00236335 [0.00213011, 0.0025966] | 500 |
| seasonal_multiscale | 3 | 0 | 0.3 | NC_PSD_minus_TCH | 0.00200825 [0.00177247, 0.00224404] | 500 |
| seasonal_multiscale | 3 | 0 | 0.3 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 3 | 0.4 | 0 | NC_PSD_minus_TC | 0.00595591 [0.00546064, 0.00645117] | 500 |
| seasonal_multiscale | 3 | 0.4 | 0 | NC_PSD_minus_TCH | 0.00462542 [0.00411938, 0.00513145] | 500 |
| seasonal_multiscale | 3 | 0.4 | 0 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 3 | 0.4 | 0.3 | NC_PSD_minus_TC | 0.00249704 [0.0021152, 0.00287889] | 500 |
| seasonal_multiscale | 3 | 0.4 | 0.3 | NC_PSD_minus_TCH | -0.00218669 [-0.00265377, -0.00171961] | 500 |
| seasonal_multiscale | 3 | 0.4 | 0.3 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0 | NC_PSD_minus_TC | NA [NA, NA] | 0 |
| seasonal_multiscale | 3 | 0.99 | 0 | NC_PSD_minus_TCH | -0.047986 [-0.0523312, -0.0436408] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0 | NC_PSD_minus_NC_raw | 0.0154825 [0.0147947, 0.0161702] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0.3 | NC_PSD_minus_TC | NA [NA, NA] | 0 |
| seasonal_multiscale | 3 | 0.99 | 0.3 | NC_PSD_minus_TCH | -0.142049 [-0.15252, -0.131578] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0.3 | NC_PSD_minus_NC_raw | 0.00989484 [0.0091054, 0.0106843] | 500 |
| seasonal_multiscale | 10 | 0 | 0 | NC_PSD_minus_TC | 2.67315e-05 [5.48412e-06, 4.79789e-05] | 500 |
| seasonal_multiscale | 10 | 0 | 0 | NC_PSD_minus_TCH | 2.47758e-05 [3.77069e-06, 4.57809e-05] | 500 |
| seasonal_multiscale | 10 | 0 | 0 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 10 | 0 | 0.3 | NC_PSD_minus_TC | 5.32681e-05 [3.29339e-05, 7.36023e-05] | 500 |
| seasonal_multiscale | 10 | 0 | 0.3 | NC_PSD_minus_TCH | -5.58854e-06 [-2.96513e-05, 1.84742e-05] | 500 |
| seasonal_multiscale | 10 | 0 | 0.3 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 10 | 0.4 | 0 | NC_PSD_minus_TC | 0.000193065 [0.000119091, 0.000267038] | 500 |
| seasonal_multiscale | 10 | 0.4 | 0 | NC_PSD_minus_TCH | 0.000135457 [6.05416e-05, 0.000210373] | 500 |
| seasonal_multiscale | 10 | 0.4 | 0 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 10 | 0.4 | 0.3 | NC_PSD_minus_TC | 0.000154651 [9.79081e-05, 0.000211395] | 500 |
| seasonal_multiscale | 10 | 0.4 | 0.3 | NC_PSD_minus_TCH | -0.00046008 [-0.000536483, -0.000383676] | 500 |
| seasonal_multiscale | 10 | 0.4 | 0.3 | NC_PSD_minus_NC_raw | 0 [0, 0] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0 | NC_PSD_minus_TC | NA [NA, NA] | 0 |
| seasonal_multiscale | 10 | 0.99 | 0 | NC_PSD_minus_TCH | 0.00822819 [0.00745825, 0.00899812] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0 | NC_PSD_minus_NC_raw | 0.00425361 [0.00399522, 0.004512] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0.3 | NC_PSD_minus_TC | NA [NA, NA] | 0 |
| seasonal_multiscale | 10 | 0.99 | 0.3 | NC_PSD_minus_TCH | -0.0038427 [-0.00459671, -0.00308869] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0.3 | NC_PSD_minus_NC_raw | 0.00165893 [0.00151329, 0.00180457] | 500 |

## Severe error-correlation stress and PSD cost

Error correlations at or above 0.9 are deliberately severe stress settings. The native PSD map and the WA SPD floor are distinct steps. In the raw-versus-PSD contrast both estimates pass through the same WA safeguard; its MSE difference therefore measures the additional native PSD map, not a comparison against an unprotected inverse. Numerical admissibility does not imply recovery of the physical generating covariance or improved merging accuracy.

| Truth | SNR | Base corr. | Signal-error corr. | Method | Attempts | Finite covariances | Estimation failure % | Indefinite % | PSD repair % | Mean theta scale | Test MSE [95% CI] | Valid merge n |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| seasonal_multiscale | 0.1 | 0.99 | 0 | TC | 500 | 0 | 100 | NA | NA | NA | NA [NA, NA] | 0 |
| seasonal_multiscale | 0.1 | 0.99 | 0 | TCH | 500 | 500 | 0 | 100 | NA | NA | 1.43924 [1.42613, 1.45235] | 500 |
| seasonal_multiscale | 0.1 | 0.99 | 0 | NC_raw | 500 | 500 | 0 | 100 | 0 | 1 | 1.07333 [1.06632, 1.08034] | 500 |
| seasonal_multiscale | 0.1 | 0.99 | 0 | NC_PSD | 500 | 500 | 0 | 0 | 100 | 0.999051 | 1.0884 [1.08058, 1.09621] | 500 |
| seasonal_multiscale | 0.1 | 0.99 | 0.3 | TC | 500 | 0 | 100 | NA | NA | NA | NA [NA, NA] | 0 |
| seasonal_multiscale | 0.1 | 0.99 | 0.3 | TCH | 500 | 500 | 0 | 100 | NA | NA | 1.60984 [1.59505, 1.62462] | 500 |
| seasonal_multiscale | 0.1 | 0.99 | 0.3 | NC_raw | 500 | 500 | 0 | 100 | 0 | 1 | 1.05056 [1.04402, 1.0571] | 500 |
| seasonal_multiscale | 0.1 | 0.99 | 0.3 | NC_PSD | 500 | 500 | 0 | 0 | 100 | 0.999714 | 1.05796 [1.05099, 1.06493] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0 | TC | 500 | 0 | 100 | NA | NA | NA | NA [NA, NA] | 0 |
| seasonal_multiscale | 3 | 0.99 | 0 | TCH | 500 | 500 | 0 | 100 | NA | NA | 0.630073 [0.625024, 0.635122] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0 | NC_raw | 500 | 500 | 0 | 100 | 0 | 1 | 0.566605 [0.563134, 0.570075] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0 | NC_PSD | 500 | 500 | 0 | 0 | 100 | 0.998407 | 0.582087 [0.578187, 0.585987] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0.3 | TC | 500 | 0 | 100 | NA | NA | NA | NA [NA, NA] | 0 |
| seasonal_multiscale | 3 | 0.99 | 0.3 | TCH | 500 | 500 | 0 | 100 | NA | NA | 0.664081 [0.658446, 0.669717] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0.3 | NC_raw | 500 | 500 | 0 | 82.2 | 0 | 1 | 0.512137 [0.503731, 0.520544] | 500 |
| seasonal_multiscale | 3 | 0.99 | 0.3 | NC_PSD | 500 | 500 | 0 | 0 | 82.2 | 0.999348 | 0.522032 [0.513056, 0.531008] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0 | TC | 500 | 0 | 100 | NA | NA | NA | NA [NA, NA] | 0 |
| seasonal_multiscale | 10 | 0.99 | 0 | TCH | 500 | 500 | 0 | 100 | NA | NA | 0.108925 [0.1082, 0.10965] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0 | NC_raw | 500 | 500 | 0 | 100 | 0 | 1 | 0.1129 [0.112178, 0.113621] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0 | NC_PSD | 500 | 500 | 0 | 0 | 100 | 0.999236 | 0.117153 [0.116251, 0.118055] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0.3 | TC | 500 | 0 | 100 | NA | NA | NA | NA [NA, NA] | 0 |
| seasonal_multiscale | 10 | 0.99 | 0.3 | TCH | 500 | 500 | 0 | 100 | NA | NA | 0.114267 [0.1134, 0.115133] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0.3 | NC_raw | 500 | 500 | 0 | 99.4 | 0 | 1 | 0.108765 [0.108124, 0.109406] | 500 |
| seasonal_multiscale | 10 | 0.99 | 0.3 | NC_PSD | 500 | 500 | 0 | 0 | 99.4 | 0.999791 | 0.110424 [0.109702, 0.111146] | 500 |

## Initialization and eigenvalue-floor sensitivity

The figure uses initialization `default_spectral` for its floor sweep and floor `1e-08` for its initialization comparison. The table retains every tested combination. Stability of this numerical experiment does not resolve physical nonidentifiability. The two-stage PSD map can increase the off-diagonal loss, so more starts do not guarantee a better post-map objective or better prediction.

| Setting | Relative floor | n | Mean loss before/after PSD | Test MSE [95% CI] | PSD repair % | Floor % | Mean weight L1 |
|---|---|---|---|---|---|---|---|
| one_spectral | 1e-10 | 100 | 0.184992 / 0.225866 | 0.537787 [0.531498, 0.544075] | 84 | 84 | 2.49399 |
| one_spectral | 1e-08 | 100 | 0.184992 / 0.225866 | 0.537787 [0.531498, 0.544075] | 84 | 84 | 2.49399 |
| one_spectral | 1e-06 | 100 | 0.184992 / 0.225866 | 0.537786 [0.531498, 0.544075] | 84 | 84 | 2.49396 |
| one_spectral | 0.0001 | 100 | 0.184992 / 0.225866 | 0.537737 [0.531449, 0.544024] | 84 | 84 | 2.49041 |
| one_spectral | 0.01 | 100 | 0.184992 / 0.225866 | 0.53345 [0.527223, 0.539676] | 84 | 90 | 2.17794 |
| default_spectral | 1e-10 | 100 | 0.184992 / 0.225866 | 0.537787 [0.531498, 0.544075] | 84 | 84 | 2.49399 |
| default_spectral | 1e-08 | 100 | 0.184992 / 0.225866 | 0.537787 [0.531498, 0.544075] | 84 | 84 | 2.49399 |
| default_spectral | 1e-06 | 100 | 0.184992 / 0.225866 | 0.537786 [0.531498, 0.544075] | 84 | 84 | 2.49396 |
| default_spectral | 0.0001 | 100 | 0.184992 / 0.225866 | 0.537737 [0.531449, 0.544024] | 84 | 84 | 2.49041 |
| default_spectral | 0.01 | 100 | 0.184992 / 0.225866 | 0.53345 [0.527223, 0.539676] | 84 | 90 | 2.17794 |
| six_spectral | 1e-10 | 100 | 0.184992 / 0.225866 | 0.537787 [0.531498, 0.544075] | 84 | 84 | 2.49399 |
| six_spectral | 1e-08 | 100 | 0.184992 / 0.225866 | 0.537787 [0.531498, 0.544075] | 84 | 84 | 2.49399 |
| six_spectral | 1e-06 | 100 | 0.184992 / 0.225866 | 0.537786 [0.531498, 0.544075] | 84 | 84 | 2.49396 |
| six_spectral | 0.0001 | 100 | 0.184992 / 0.225866 | 0.537737 [0.531449, 0.544024] | 84 | 84 | 2.49041 |
| six_spectral | 0.01 | 100 | 0.184992 / 0.225866 | 0.53345 [0.527223, 0.539676] | 84 | 90 | 2.17794 |
| random8 | 1e-10 | 100 | 0.182588 / 0.257577 | 0.546743 [0.540402, 0.553084] | 100 | 100 | 2.82272 |
| random8 | 1e-08 | 100 | 0.182588 / 0.257577 | 0.546743 [0.540402, 0.553084] | 100 | 100 | 2.82272 |
| random8 | 1e-06 | 100 | 0.182588 / 0.257577 | 0.546743 [0.540402, 0.553084] | 100 | 100 | 2.82267 |
| random8 | 0.0001 | 100 | 0.182588 / 0.257577 | 0.546673 [0.540333, 0.553013] | 100 | 100 | 2.81783 |
| random8 | 0.01 | 100 | 0.182588 / 0.257577 | 0.540803 [0.534509, 0.547097] | 100 | 100 | 2.40236 |
| shift_low | 1e-10 | 100 | 0.184631 / 0.231356 | 0.537762 [0.53142, 0.544105] | 93 | 93 | 2.63489 |
| shift_low | 1e-08 | 100 | 0.184631 / 0.231356 | 0.537762 [0.53142, 0.544105] | 93 | 93 | 2.63488 |
| shift_low | 1e-06 | 100 | 0.184631 / 0.231356 | 0.537762 [0.531419, 0.544104] | 93 | 93 | 2.63484 |
| shift_low | 0.0001 | 100 | 0.184631 / 0.231356 | 0.537705 [0.531364, 0.544046] | 93 | 93 | 2.63076 |
| shift_low | 0.01 | 100 | 0.184631 / 0.231356 | 0.532837 [0.526607, 0.539067] | 93 | 100 | 2.27154 |
| shift_high | 1e-10 | 100 | 0.184436 / 0.264525 | 0.539336 [0.532908, 0.545764] | 97 | 97 | 2.93134 |
| shift_high | 1e-08 | 100 | 0.184436 / 0.264525 | 0.539336 [0.532908, 0.545764] | 97 | 97 | 2.93134 |
| shift_high | 1e-06 | 100 | 0.184436 / 0.264525 | 0.539335 [0.532907, 0.545763] | 97 | 97 | 2.93129 |
| shift_high | 0.0001 | 100 | 0.184436 / 0.264525 | 0.539262 [0.532836, 0.545688] | 97 | 97 | 2.92621 |
| shift_high | 0.01 | 100 | 0.184436 / 0.264525 | 0.533123 [0.526824, 0.539421] | 97 | 100 | 2.48656 |

## Identifiability and numerical-optimization qualifications

The exported p=2 and p=3 tables contain 3 and 2 candidate decompositions, respectively. The p=2 example permits distinct PSD residuals with the same loss; its candidates imply different loadings, so known equal gains would supply additional identifying information. The p=3 common-mode example demonstrates covariance ambiguity while its unit-sum optimal weights remain unchanged. Thus nonidentifiability does not force weight changes in every special case.

For a compatible p=3 covariance with nonnegative TC residual variances, the TC signal products can supply a zero-off-diagonal-loss feasible candidate. A positive final NC loss may therefore include numerical optimization approximation. The projected-subgradient stopping rule is not a certificate of global optimality, and a measured NC/TC accuracy difference must not be presented as a theorem of intrinsic statistical inferiority.

The stationary/seasonal comparison uses common error realizations when signal-error correlation is zero. With unit gains, truth cancels in pairwise differences; nearly identical TCH fitted variance results can therefore be structural. Its first-output full-covariance reconstruction can still vary with Q, and severe near-collinearity can amplify optimizer roundoff sensitivity.

All results describe the supplied synthetic generator. Do not claim universally superior accuracy, universal failure-free behavior, identified physical error covariances, or completed global validation from these outputs.
