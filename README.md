# heartrate
This folder contains two methods of detecting heartbeats. xcorr_detector.m uses cross correlation and cwt_detector.m uses a continuous wavelet transform. 

1. xcorr_detector.m takes in a .wav file (containing the heartbeat vibrations) and template .wav file (containing a template heartbeat) and starts by applying a lowpass equiripple filter with a passband frequency of 200Hz and stopband frequency of 300Hz. Both files are then resampled at 1/10 the original sampling rate to increase the speed of the algorithm. The algorithm finds the max cross correlation between the template and the heartbeat file and finds the peaks in cross correlation and where they occur. 

...The find_rate function takes in the inputs pks (a vector containing the height of the peaks), locs (a vector containing the location of those peaks), rs_beat (the resampled and filtered heartbeat file), and rs (the sampling rate that the file is resampled to). Its output is a matrix such that the first column contains the heartrates, the second column contains the times at which they occur, and the third column contains the ratio of hte height of the peak and the root mean square of the surrounding signal, which can be interpreted as a measure of confidence of the heartbeat. 

...The find_similarity function takes in the inputs real (a matrix that contains the times at which heartbeats occur in the second column), rate (the output of find_rate), rs_beat, and rs (which are the same as those found in find_rate). It uses a normalized cross correlation to output the similarity (out of 1) between "real" and "rate".

...The find_variability function takes in rate (the output of find_rate), rs_beat, and rs (which are the same as those found in find_rate). It outputs the standard deviation per second of the heartrates. 

2. cwt_detector.m also filters and resamples the heartbeat .wav file and uses a template to find the continuous wavelet transform between the template and heartbeat file. 
