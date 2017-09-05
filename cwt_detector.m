%created by Jeff Chow to detect heartrates and heartrate variability from
%vibration sensors
close all
clear
[beat,Fs] = audioread('clear_heart_rate.wav');
[tmp,Fs] = audioread('single_beat.wav');

Hd = fdesign.lowpass('Fp,Fst,Ap,Ast',200,...
    400,1,60,Fs); %lowpass filter w/ passband frequency 200Hz & stopband frequency 300Hz
d = design(Hd,'equiripple'); % equiripple filter
ft_beat = filter(d,beat);
ft_tmp = filter(d,tmp);

rs_beat = resample(ft_beat, 1, 10); %resample file at 1/10 sampling rate
rs_tmp = resample(ft_tmp, 1, 10);

[wt,f] = cwt(rs_beat, rs_tmp, 4410); 
[~,idx125] = min(abs(f-125));
csf = wt(idx125,:);
plot(abs(csf))