clear all; close all; 

% Daily power spectral density binned by solar elevation 
cd ..
load 'data/BlueCIbasis.mat'
    
% Handle for month
dv = datevec(B.time); mo = dv(:,2); 

% Quartiles for daily CI and daily categorical CI
[nfreq,ntime,ns] = size(B.sm);
% First examine counts, which should be independent of frequency
sn = {'n','dd','d','all'};
pcount = squeeze(B.ct(1,:,:)); pcount = pcount';
pcount(4,:) = sum(pcount(1:3,:));
        
figure(1); clf; 

mincount = [380*12 100*12 500*12 1200*12]; % night, dd, day, total; x12 to account for switch from 1-min to 5-sec LTSA
        
for S = 1:4;
    eval(['subplot(41' int2str(S) ')']);
    plot(B.time,pcount(S,:),'-o');
    ylabel('Observation count');
    axis tight; datetick('x'); hold on;
    xl = get(gca,'Xlim'); plot(xl,mincount(S)+[0 0],'r--');
    idx = find(pcount(S,:) < mincount(S));
    spct = round(10000*numel(idx)/ntime)/100;
    title([sn{S} ': ' num2str(spct) '% of days < minimum threshold']);
end

% Get daily call indices by SE category, screened
clear DCI; DCI.time = B.time; DCI.sn = sn; 
for S = 1:3;
    sm = B.sm(:,:,S); ct = B.ct(:,:,S);
    % remove data for days with insufficient sampling
    xcl = find(ct(1,:) < mincount(S)); ct(:,xcl) = NaN;
    % Place into structure
    clear L; L.time = B.time; L.freq = B.freq; L.ltsa = sm./ct;
    % Get the daily call index for this SE category
    c = call_index(L); DCI.blue(S,:) = c.blue;
end
% Get daily call index for all data, regardless of SE category
sm = sum(B.sm,3); ct = sum(B.ct,3);
S = 4;
xcl = find(ct(1,:) < mincount(S)); ct(:,xcl) = NaN;
clear L; L.time = B.time; L.freq = B.freq; L.ltsa = sm./ct;
c = call_index(L); DCI.blue(S,:) = c.blue;
q = DCI.blue([1 3],:); % isolate night and day
q = q - min(q(:));  % scale range above minimum
S = 5; DCI.blue(S,:) = q(1,:)./q(2,:);

% Plot the daily values by category
sn = {'n','dd','d','all','ratio'};
DCI.sn = sn;

figure(2); clf; 
for S = 1:5;
    eval(['subplot(23' int2str(S) ')']);
    plot(DCI.time,DCI.blue(S,:),'o');
    ng = numel(find(~isnan(DCI.blue(S,:))));
    title([sn{S} ': ' int2str(ng) ' good indices']);
    datetick('x'); axis tight;
end

% Save dates and CI values as .csv
a=cellstr(datestr(DCI.time));
t=cell2table(a);
writetable(t,'ci_daily_seq_dates3.csv')
csvwrite('ci_daily_seq3.csv',DCI.blue)