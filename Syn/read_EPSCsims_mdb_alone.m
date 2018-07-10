
clear
clc

fbase = 'fig11_PFCapic';
%WTD1_Nov8IR3a
% gAMPA = 0.00021;
% tau1 = 2.7006;
% tau2 = 2.1829;

%WTD2_Nov8IR2b
% gAMPA = 0.000228;
% tau1 = 6.2513;
% tau2 = 1.9087;

%HETD1_Nov7IR3a
% gAMPA = 0.0002030;
% tau1 = 4.532;
% tau2 = 1.7101;

%HETD2_Apr20IR3a
gAMPA = 0.000225;
tau1 = 4.5559;
tau2 = 1.8262;

writeSmry = 0;

inbase = sprintf('%s_tR%.4f_tF%.4f_gAMP%.7f',fbase,tau1,tau2,gAMPA);

figttl = sprintf('PFC test: t_1 %.4f t_2 %.4f gAMP %.7f',tau1,tau2,gAMPA);
makeDistPlot = 1;
TipsOnly = 1;

[t,v]=readNRNbin_Vclamp(inbase,0);

ylim([-0.05 0])
txt_fname = sprintf('%s_dist.txt',inbase);
[dat] = dlmread(txt_fname);
nSyn = size(dat,1);

tmp_ras=importdata('test_raster.txt');

sTimes = tmp_ras(:,1);
spk_ind=tmp_ras(:,2);
l_spks=length(sTimes);

idx=zeros(l_spks-1,2);
PSCbase=zeros(l_spks-1,1);

v_EPSP = zeros(60000,l_spks-1);
peak = [];
amp = [];
rise = [];
decay = [];
hfw = [];
% figure(1)
emp_spks=[1 l_spks];
emp_spks2=[];
for k=2:l_spks-1
idx(k,:) = [min(find(t >= sTimes(k)))  max(find(t < sTimes(k+1)))];
% idx(k,2)-idx(k,1);
% figure(1)
% plot(v(idx(k,1):idx(k,2)))
% hold on;
tstep = [0 : idx(k,2)-idx(k,1)-1];
t_EPSP = t(1+tstep);


    idx_EPSP = idx(k,1);
    tmp = v(idx(k,1)+1:idx(k,2));
    [mn,mnI]=min(tmp);
    
    v_EPSP(1+tstep,k-1) = tmp-tmp(1);
%     plot(v_EPSP(:,k-1))
%     hold on;
    peak(end+1,:) = [t_EPSP(mnI)  abs(mn-tmp(1))];
    
    [amp_tmp,rise_tmp,decay_tmp,hfw_tmp]=analyze_EPSC(t_EPSP,tmp); 
    if (decay_tmp~=0)
        amp(end+1)=amp_tmp;
        rise(end+1)=rise_tmp;
        decay(end+1)=decay_tmp;
        hfw(end+1)=hfw_tmp;
    else
        emp_spks=[emp_spks k];
        emp_spks2=[emp_spks2 k];
    end
end

v_EPSP(30000:end,:)=[];

% nSyn = size(v_EPSP,1);
fprintf('Found %d EPSCs\n',length(amp)-length(find(amp==0)));
fprintf('Mean amp\t%.2f\n',mean(nonzeros(amp)*1e3));
fprintf('Mean rise\t%.2f\n',mean(nonzeros(rise)));
fprintf('Mean decay\t%.2f\n',mean(nonzeros(decay)));
fprintf('Mean 1/2 width\t%.2f\n',mean(nonzeros(hfw)));
mean_EPSP=zeros(length(v_EPSP),1);

targ_spks=find(dat(:,2)>=0.67*max(dat(:,2)));

for i=1:length(v_EPSP)
mean_EPSP(i) = mean(v_EPSP(i,targ_spks));
end

figure(1)
l_EPSP=length(v_EPSP);
t_plot=0:0.025:0.025*(l_EPSP-1);
for i=1:length(find(dat(:,2)>=0.67*max(dat(:,2))))
    if (ismember(i,emp_spks2)==0)
        if(max(v_EPSP(:,i)*1000)<0.5)
        plot(t_plot,v_EPSP(:,i)*1000,'k',t_plot,mean_EPSP*1000,'r');
        hold on;
        end
    end
end
plot(t_plot,mean_EPSP*1000,'r');
ttl = sprintf('%s:  superimposing all EPSCs',figttl);
title('superimposed EPSCs');
xlabel('time (ms)');
xlim([0 50])
ylim([-15 1])
ylabel('EPSC (pA)');
fname = sprintf('%s_meanEPSCond.fig',inbase);
% figure(2)
% 
% distan=zeros(nSyn,1);
% for i=1:length(dat)
%     distan(i)=sqrt(dat(i,2)^2+dat(i,3)^2);
% end
% ave_persyn=zeros(nSyn,1);
% spk_indtmp=spk_ind;
% spk_indtmp(emp_spks)=[];
% spks_syn=[spk_indtmp amp'];
% for i=1:length(dat)-2
%     tmp_spks=find(spks_syn(:,1)==i);
%     ave_persyn(i)=spks_syn(tmp_spks,2);
% end
% for i=1:length(dat)
% if ((ave_persyn(i)~=0)&&(dat(i,2)>=0.67*max(dat(:,2))))
% plot(dat(i,2),ave_persyn(i)*1000,'ro')
% hold on;
% end
% end
% title('Synapse location and EPSC amplitudes');
% xlabel('Distance from soma (microns)');
% ylabel('Max EPSC amplitude (pA)');


% subplot(1,3,3);
% xbar=[.5:.1:60];
% n_elements=histc(peak(:,2)*1000,xbar);
% npts=max(size(peak));
% c_elements = cumsum(n_elements);
% 
% plot(xbar,c_elements/npts,'o-');

% xlim([0 60]);

% title('Cumulative frequency histograms');
% xlabel('EPSC amplitude (pA)');
% ylabel('cumulative (%)');
% avamp=mean(amp)*1e3;
% tr=mean(rise);
% td=mean(decay);
% avhfw=mean(hfw);
% fprintf('Found %d EPSCs\n',length(amp));
% fprintf('Mean amp\t%.2f\n',mean(amp)*1e3);
% fprintf('Mean rise\t%.2f\n',mean(rise));
% fprintf('Mean decay\t%.2f\n',mean(decay));
% fprintf('Mean 1/2 width\t%.2f\n',mean(hfw));

figure(2)
plot(t_plot,mean_EPSP*1000,'r');
ttl = sprintf('%s:  superimposing all EPSCs',figttl);
title('superimposed EPSCs');
xlabel('time (ms)');
xlim([0 40])
ylim([-18 1])
ylabel('EPSC (pA)');
