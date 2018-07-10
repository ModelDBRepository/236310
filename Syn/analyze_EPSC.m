function [amp, rise, decay, hfwidth, rise10_90, tc_dcy] = analyze_EPSC(t,cur)
% function [amp, rise, hfwidth, rise10_90] = analyze_EPSC(t,cur)
%
%   assume we are only sent time vs current data of a single EPSC.

%
%  Find the 'easy' stats: amplitude and half-width.
%
[mn, mnI] = min(cur);
[mx, mxI] = max(cur);
Ishift = max(cur(1),cur(end))-cur';

amp = abs(cur(mxI)-cur(mnI));
hfmx = (mn+mx)/2;
hfidx = find(cur<hfmx);
hfbds = [hfidx(1) hfidx(end)];
hfwidth = abs(t(hfbds(2))-t(hfbds(1)));
%plot(t,Ishift)
%
%  now fit exponential.  This is hard because the data are negative and
%  near zero, so linearization doesn't work well. 
%

% Only fit to the part of the curve where the EPSC is 'active':
% assume EPSC starts when when abs(dI/dt) > 1e-5 and ends when abs(dI/dt) <
% 1e-5 for the last time. 
%

[dI]=get_dVdt(t,cur);
% plot(dI)
[EPidx]=find(abs(dI)>1e-5);
% fprintf('This is the length of EPidx %i\n',length(EPidx));
if (length(EPidx)<1)
    rise=0;
    decay=0;
    rise10_90=0;
    amp=0;
    return;
else
    EPst  = EPidx(1);
EPend = EPidx(end);

rise = t(mnI) - t(EPst);
[mn10,idx10] = min(find(abs(Ishift)>=0.1*amp));
[mn90,idx90] = min(find(abs(Ishift)>=0.9*amp));
rise10_90 = t(mn90)-t(mn10);

decay = t(EPend)-t(mnI);
EPbds = [EPst mnI EPend];


% now:  decay time constant as fit to exponential.
tidx = mnI:EPend;


l=length(t(tidx));
t_tmp=t(tidx);
Ishift_tmp=Ishift(tidx);
% plot(t_tmp,Ishift_tmp);
if (size(Ishift_tmp)>0)
    e_decay=Ishift_tmp(1)*exp(-1);
    [I_crit idx_crit]=min(find(Ishift_tmp<=e_decay));
    rise = rise10_90;
    decay = I_crit*0.025;
else
    rise = 0;
    decay = 0;
    amp=0;
    hfwidth=0;
end

return;
end
function [dVdt] = get_dVdt(tvec, vvec)
% function [dVdt,dVall] = get_dVdt(tvec, vvec)
%
%   get_dVdt    use the central difference formula to estimate the
%               derivative of the second vector with respect to the first. 
%
%   INPUT       TVEC, vector of times
%               VVEC, vector of membrane potential, for example.
%   OUTPUT      dVdt, the derivative of VVEC with respect to TVEC.
%
%   Christina Weaver, christina.weaver@mssm.edu, July 2005
%

for i = 2:length(tvec)-1 
    dVdt(i) = (vvec(i+1)-vvec(i-1))/(tvec(i+1)-tvec(i-1));
end;


