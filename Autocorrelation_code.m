%% Autocorrelation code

% this code estimates the periodicity of a given trace by finding the
% cross-correlation values of the trace against itself and taking the mean
% of the distance betwen the peaks in the reulting output which  have
% passed a bootstrapping threshold.

% The data we inputted is an array of two generation mVenus-Hes1 traces which
% had a linear trend removed.

%Hes1_detrended_time_traces = "Insert data here";


Autocorrelation_periodicity = nan(1,length(Hes1_detrended_time_traces));

% loop through all traces
for cell_index = 1:length(Hes1_detrended_time_traces)
    
    % Assign y to the trace vector
    y = Hes1_detrended_time_traces;
    
    % Assign t to the time vecotr
    t = 1:length(Hes1_detrended_time_traces);
    
    % run autocorrelation function on y
    % "cr" is a vector of autocorrolation values and "lag" is a time vecotr
    % associated with the lags
    
    [cr,lag]=xcorr(y);
    
    %% boostrap to get confidence bounds from noise
    
    % This performs autocorrelation (AC) on 100 different random permutaions of
    % the vector y and provides a threshold boundary of autocorrelation
    % values. We also decide here how stringent our threshold will be, that
    % is, how many standard deviations away from the randomised AC values
    % must our threshold be for a peak in the original AC function to be
    % considered nonrandom. E.g. in this case beloiw we choose 1.
    
    BootMat=[]; 
    for k=1:100
        % bootstrap
        kidx=randperm(numel(y));
        randvec=y(kidx);
        [cboot,lag]=xcorr(randvec);
        BootMat=[BootMat cboot(:)];
    end
    lag=lag*(t(2)-t(1));
    sd=std(BootMat');
    how_many_standard_deviations = 1;
    m=-how_many_standard_deviations*sd; m(lag==0)=NaN;
    M=how_many_standard_deviations*sd; M(lag==0)=NaN;
    
   % We now ask which peaks in "cr" surpass our threshold and thus can be
   % accepted as a nonrandom peqak in the autocorrelation function.
    
    [pks,locs]=findpeaks(cr,'MinPeakProminence',0.10*max(cr));
    
    peaks_and_locations=[pks,locs];
    
    threshold_passing_peaks_and_locs = [];
    
    for i = 1:length(locs)       
        if peaks_and_locations(i,2) >= length(Hes1_detrended_time_traces{cell_index}) &&...
                peaks_and_locations(i,1) > M(locs(i))
            
            threshold_passing_peaks_and_locs(end+1,:) = peaks_and_locations(i,:);           
        end       
    end
    
    % We also check if we didn't have any "passing" peaks.
    
    if ~isempty(threshold_passing_peaks_and_locs)
        threshold_passing_peaks_and_locs(:,2)=(threshold_passing_peaks_and_locs(:,2)-length(Hes1_detrended_time_traces{cell_index}));
        
        Peaks_Above_threshold = [0;threshold_passing_peaks_and_locs(:,2)]./4;
        
        
        % We take the mean of all the distances between the peaks in the AC
        % function which passed our theshold. We consider this value our
        % periodicity for this two-generation trace
        
        Autocorrelation_periodicity(cell_index) = nanmean(diff(Peaks_Above_threshold));
    end
end


