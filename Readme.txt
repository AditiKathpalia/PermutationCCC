%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		Permutation Compression Complexity Causality (PCCC) Toolbox
				Version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This is a Permutation Compression Complexity Causality or PCCC (MATLAB) toolbox which uses 'Compression Complexity Causality (CCC)' and 'Effort to Compress (ETC)', published in the following works:

1. Kathpalia, A., & Nagaraj, N. (2019). Data-based intervention approach for Complexity-Causality measure. PeerJ Computer Science, 5, e196.
2. Nagaraj, N., Balasubramanian, K., & Dey, S. (2013). A new complexity measure for time series analysis and classification. The European Physical Journal Special Topics, 222(3), 847-860.


This toolbox can be used standalone after including it on the MATLAB path. It is meant for research purposes only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Run 'Main_PCCC_surr_analysis.m' for a demo using coupled time-series from Rossler systems in 'Rossler_data_demo.mat'

List of files and their functionality:
1.  'Main_PCCC_surr_analysis.m' loads Rossler data and estimates PCCC between the coupled processes in both directions and determines whether PCCC in either direction is significant.
2. 'perm_binning_std.m' is used to bin/ symbolize the time series using ordinal patterns embedding.
3. 'Partition.m' is used for equidistant binning/ symbolization of the time series.
4. 'CCC_binned_seqs.m' estimates CCC between pre-binned sequences. This is called PCCC when the cause is binned using permutation/ ordinal patterns. This is what is being done in the demo.
5. 'AAFTsur.m' generates surrogates of given data using the Amplitude Adjusted Fourier Transform method.
6.  'stationary_bootstrap.m' generates surrogates of given data using the stationary bootstrap method.
7. 'ETC_1D.m' and 'ETC_2D.m' are employed by the 'CCC_binned_seqs.m' to compute CCC. These functions compute individual and joint ETC values respectively.
