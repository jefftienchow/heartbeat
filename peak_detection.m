%created by Jeff Chow to detect heartrates and heartrate variability from
%vibration sensors
close all
clear
[beat,Fs] = audioread('vibration_sensor.wav');
[tmp,Fs] = audioread('single_beat.wav');
text = importdata('clear_label.txt');

Hd = fdesign.lowpass('Fp,Fst,Ap,Ast',200,...
    400,1,60,Fs); %lowpass filter w/ passband frequency 200Hz & stopband frequency 300Hz
d = design(Hd,'equiripple'); % equiripple filter
ft_beat = filter(d,beat);
ft_tmp = filter(d,tmp);

rs_beat = resample(ft_beat, 1, 10); %resample file at 1/10 sampling rate
rs_tmp = resample(ft_tmp, 1, 10);
rs=3333.3; %sampling frequency after file is resampled

stepsize=2.5;
strt=1;
stp=strt+length(rs_tmp);
acc=[];
S=0;

for i=1:length(rs_beat)/(length(rs_tmp)/stepsize)-1 
    S=S+1;
    z = xcorr(rs_beat(strt:stp),rs_tmp);
    [cur_max,index] = max(z);
    acc(i,1) = cur_max;
    acc(i,2) = index;
    stp = round(1+(i+1)*length(rs_tmp)/stepsize);
    strt = round(1+i*length(rs_tmp)/stepsize);
end
%takes the max of the cross correlation taken from segments of the file of
%of length equal to the size of the template; this is incremented by the
%size of the template divided by 4

[peak1,tmp_locs1] = findpeaks(acc(:,1)); 
%peak1 is the height of the peak in the max cross correlations; tmp_locs
%are the indices at which they occur

locs1 = []*length(tmp_locs1);
for i = 1:length(peak1)
    locs1(i) = tmp_locs1(i)*length(rs_tmp)/stepsize + (acc(tmp_locs1(i),2))-length(rs_tmp/stepsize);
end
%locs1 contains the locations of the peaks in cross correlations in 
%relation to the whole file

[peak2,tmp_locs2] = findpeaks(peak1);
locs2 = []*length(tmp_locs2);
for i = 1:length(peak2)
    locs2(i) = round(locs1(tmp_locs2(i)));
end

%locs2 contains the locations of the peaks of peak1

real = []; %matrix containing the heartbeats detected by hand 
for i=1:(length(text)-1)
    real(i,1) = text(i+1)-text(i);
    real(i,2) = text(i);
end

rate = find_rate(peak2,locs2,rs_beat,rs);
similarity = find_similarity(real,rate,rs_beat,rs);
variability = find_variability(rate,rs_beat,rs);

%***the following is for the purpose of graphing only***
x_val = []*length(rate(:,2));
for i = 1:length(rate(:,2))
    x_val(i,1) = rate(i,2);
    x_val(i,2) = .2;
end

time =[]*length(rs_beat);
for i=1:length(rs_beat)
    time(i) = i/rs;
end
location = []*length(locs1);
for i =1:length(locs1)
    location(i) = locs1(i)/rs;
end
variability_time = []*length(variability);
for i = 1:length(variability)
    variability_time(i) = i;
end

figure
subplot(2,1,1)
scatter(x_val(:,1),x_val(:,2),'r')
hold on;
plot(time,rs_beat)
hold on;
plot(location,peak1)
hold off
legend('heartbeat','raw data','cross correlation')
subplot(2,1,2)
plot(variability_time, variability)

function y = find_rate(pks,locs,rs_beat,rs)
    %returns the matrix "rate" where 
    %rate(:,1) contains the heartrates
    %rate(:,2) contains the times at which the heartbeats occur
    %rate(:,3) contains the confidence level measured in the ratio of
    %the height of the peak and the rootmean square of the surrounding 
    %signal
    rl_pk_loc = [];
    rl_pk = [];
    S = 0;
    cur = 0;
    R = 1;
    U = 1;
    max_pks = [];
    for i = 1:(length(rs_beat))/(.14*rs)-1
        list_pks = [];
        T = 0;
        
        while locs(U)<(i*.14*rs)
            U = U+1;
            T = T+1;
            list_pks(T) = pks(U);
        end
        %if the peak is within the next .14sec, append it to the list of
        %peaks
        if isempty(list_pks)==1
            max_pks(i) = 0;
        else
            max_pks(i) = max(list_pks);
        end
        %take the max of the list of peaks
        if i>10
            cur = mean(max_pks((i-10):i))/2;
            
        end
        
        while locs(R)<(i*.14*rs)
            R = R+1;
            if pks(R) > cur
                S = S+1;
                rl_pk_loc(S) = locs(R);
                rl_pk(S) = pks(R);
            end
        end
        %if the peak is above the current threshold, append it to the list
        %of real peaks
    end
    S=1;
    heartbeat_loc = [];
    
    for i = 1:(length(rl_pk_loc)-1)
        cur_max = rl_pk(i);
        index = i;
        T=1;
        while (rl_pk_loc(T+i)-(rl_pk_loc(i)))/rs < .1 & T+i<length(rl_pk_loc)
            
            if rl_pk(T+i)>cur_max
                cur_max = rl_pk(T+i);
                index = T+i;
            end
            %find the highest peak within the next .1sec of every real peak
            
            T = T+1;
        end  
        
        if S == 1
            heartbeat_loc(1,1) = rl_pk_loc(index);
            heartbeat_loc(1,2) = rl_pk(index);
            S=S+1;
        else
            
            if heartbeat_loc(S-1,1)~=rl_pk_loc(index)
                
                heartbeat_loc(S,1) = rl_pk_loc(index);
                heartbeat_loc(S,2) = rl_pk(index);
                S=S+1;
            end    
        end    
        %add these highest peaks to the list of heartbeats, ignoring
        %repeats
        
    end
    rate = [];
    S = 0;
    temp = 0;
   
    for i = 1:(length(heartbeat_loc(:,1))-1)

        if (heartbeat_loc(i+1,1)-heartbeat_loc(i,1))/rs + temp< .1
            temp = temp+(heartbeat_loc(i+1,1)-heartbeat_loc(i,1))/rs;
        %if the heartbeats are within .1sec of each other, ignore the 
        %peak until the time between hearbeats time is greater than or
        %equal to .1sec
        else
           S = S+1;
           rate(S,1) = (heartbeat_loc(i+1,1) - heartbeat_loc(i,1))/rs + temp;
           rate(S,2) = heartbeat_loc(i+1,1)/rs;
           temp = 0;
           
           if heartbeat_loc(i+1,1)-350<0||heartbeat_loc(i+1,1)+350>length(rs_beat)
               rate(S,3) = 0;
               disp(heartbeat_loc(i+1,1))
           else
               rate(S,3) = heartbeat_loc(i+1,2)/rms(rs_beat(heartbeat_loc(i+1,1)-350:heartbeat_loc(i+1,1)+350));
           end
        end
        
    end
    y = rate;

end

function similarity = find_similarity(real,rate,rs_beat,rs)
    %returns the similarity (out of 1) between two different measures of
    %the heart rate
    S = 1;
    R = 1;
    list_xcorr = []*length(rs_beat)/rs;
    for i = 1:(length(rs_beat))/rs
        real_range = 0;
        rate_range = 0;
        
        while real(S,2)<i & S<length(real(:,2))
            real_range = real_range+1;
            S = S+1;
        end
        %keep track of the number of heartbeats detected in one second
        
        while rate(R,2)<i & R<length(rate(:,2))
            rate_range = rate_range+1;
            R = R+1;
        end    
        if S==length(real(:,2)) || R==length(rate(:,2))
            break
        end

        if real_range<rate_range
            z = normxcorr2(real((S-real_range):S),rate((R-rate_range):R));
            
        else
            z = normxcorr2(rate((R-rate_range):R),real((S-real_range):S));
        end
        cur_max = max(z);
        %take the max cross correlation between the heartbeats detected
        %within the last second
        list_xcorr(i) = cur_max;
    end
    similarity = mean(list_xcorr);
end

function deviation = find_variability(rate,rs_beat,rs)
    %returns the standard deviation per second of the heart rates
    seconds = length(rs_beat)/rs;
    S = 0;
    T = 1;
    deviation = []*seconds;
    for i = 1:(seconds)
        S = S+1;
        R = 1;
        list_rates = [];
        while rate(T,2)<i & T<length(rate(:,2))
            list_rates(R) = rate(T,1);
            R = R+1;
            T = T+1;
        end    
        deviation(S) = std(list_rates);
    end    
end
