close all
clear
[beat,Fs] = audioread('clear_heart_rate.wav');
[tmp,Fs] = audioread('single_beat.wav');
text = importdata('clear_label.txt');

d = designfilt('lowpassfir', ...
    'PassbandFrequency',0.01,'StopbandFrequency',0.02, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod','equiripple');
ft_beat = filter(d,beat);
ft_tmp = filter(d,tmp);

rs_beat = resample(ft_beat, 1, 10);
rs_tmp = resample(ft_tmp, 1, 10);


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

[pks,tmp_locs] = findpeaks(acc(:,1));
locs = []*length(tmp_locs);
for i = 1:length(pks)
    locs(i) = tmp_locs(i)*length(rs_tmp)/stepsize + (acc(tmp_locs(i),2))-length(rs_tmp/stepsize);
end
real = [];
for i=1:(length(text)-1)
    
    if text(i+1)-text(i)<.135
        real(i,1) = .135;
        real(i,2) = text(i);
    elseif text(i+1)-text(i)>.146
        real(i,1) = .146;
        real(i,2) = text(i);
    else
        real(i,1) = text(i+1)-text(i);
        real(i,2) = text(i);
    end
end

rate = find_rate(pks,locs,rs_beat,4100);


x_val = []*length(rate(:,2));
for i = 1:length(rate(:,2))
    x_val(i,1) = rate(i,2)*4100;
    x_val(i,2) = .2;
end
figure
subplot(2,1,1)
scatter(x_val(:,1),x_val(:,2),'r')
xlim([1 20000])

hold on;
plot(rs_beat)
xlim([1 20000])

subplot(2,1,2)
plot(locs,pks)
xlim([1 20000])
    


disp(find_similarity(real,rate,rs_beat,4100))

function y = find_rate(pks,locs,beat,rs)
    rl_pk_loc = [];
    rl_pk = [];
    S = 0;
    cur = 0;
    R = 1;
    U = 1;
    debugger=0;
    for i = 1:(length(beat))/(.14*rs)
        list_pks = [];
        T = 0;
        
        while locs(U)<(i*.14*rs)
            U = U+1;
            T = T+1;
            list_pks(T) = pks(U);
            
        end
        
        max_pks(i) = max(list_pks);
        
        
        if i>10
            cur = mean(max_pks((i-10):i))/1.75;
            if debugger<20
                debugger=debugger+1;
                disp(cur)
            end
        end
        
        while locs(R)<(i*.14*rs)
            R = R+1;
            if pks(R) > cur
                S = S+1;
                rl_pk_loc(S) = locs(R);
                rl_pk(S) = pks(R);
            end
        end
    end
    rate = [];
    S = 0;
    temp = 0;
    for i = 1:(length(rl_pk_loc)-1)

        if (rl_pk_loc(i+1)-rl_pk_loc(i))/rs + temp< .13
            temp = temp+(rl_pk_loc(i+1) - rl_pk_loc(i))/rs;
        
        %{
        elseif (rl_pk_loc(i+1) - rl_pk_loc(i))/4410 + temp> .16
            K=1;
            while ((rl_pk_loc(i+1) - rl_pk_loc(i))/4410 + temp)/K>.16
                K=K+1;
            end
            for j = 1:K
                S = S+1;
                rate(S) = ((rl_pk_loc(i+1) - rl_pk_loc(i))/4410 + temp)/K;
            end
            temp = 0;
%}
        
            else

                S = S+1;

                rate(S,1) = (rl_pk_loc(i+1) - rl_pk_loc(i))/rs + temp;
                rate(S,2) = rl_pk_loc(i)/rs;

                temp = 0;
        end
    end
    y = rate;

end
   
function similarity = find_similarity(real,rate,beat,rs)
    S = 1;
    R = 1;
    list_xcorr = []*length(beat)/rs;
    for i = 1:(length(beat))/rs
        real_range = 0;
        rate_range = 0;
        while real(S,2)<i & S<length(real(:,2))
            real_range = real_range+1;
            S = S+1;
            
        end
        while rate(R,2)<i & R<length(real(:,2))
            rate_range = rate_range+1;
            R = R+1;

        end    
        if S==length(real(:,2)) | R==length(rate(:,2))
            break
        end
      
        z = xcorr(real((S-real_range):S),rate((R-rate_range):R));
        mex = max(z);
        list_xcorr(i) = mex;
    end
    similarity = mean(list_xcorr);
end