close all
clear
load adapted_wavelet.mat 
[beat,Fs] = audioread('clear_heart_rate.wav');
[tmp,Fs] = audioread('single_beat.wav');
%{
d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.01,0.02,1,60);
Hd = design(d,'equiripple');
ft_beat = filter(Hd,beat);

ft_tmp = filter(Hd,tmp);
%}
[wt,f] = cwt(beat, Y, 44100); 
[~,idx125] = min(abs(f-125));
csf = wt(idx125,:);
plot(abs(csf))