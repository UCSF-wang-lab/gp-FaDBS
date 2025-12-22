function STFT = calcRCS_STFT(aligned_data,gapFillType,windowLength,percentOverlap,nfft)
%% 
% A = calcRCS_STFT(aligned_data,[],0.1,0.5,[]);
% A = calcRCS_STFT(aligned_data,[],[],0.996,[]); % Legacy


%%
if ~exist('gapFillType','var') || isempty(gapFillType)
    gapFillType = 'blank';
end

if ~exist('windowLength','var') || isempty(windowLength)
    windowLength = inf;
end

if ~exist('percentOverlap','var') || isempty(percentOverlap)
    percentOverlap = 0.5;
end

if ~exist('nfft','var') || isempty(nfft)
    nfft = [];
end

%% Left RCS
if isfield(aligned_data,'left_LFP_table')
    % Spectrogram hyperparameters
    left_sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);
    
    if isinf(windowLength)
        WINDOW = left_sr;
    elseif windowLength > 0 && windowLength <= 1
        WINDOW = left_sr*windowLength;
        WINDOW = WINDOW + mod(WINDOW,2);
    else
        WINDOW = left_sr*windowLength;
        WINDOW = WINDOW + mod(WINDOW,2);
    end
    
    NOVERLAP = round(WINDOW*percentOverlap);
    
    if isempty(nfft)
        nfft = 2^nextpow2(left_sr);
    end
    
    chan_tag_inds = cellfun(@(x) contains(x,'chan'),aligned_data.DeviceSettings.Left.timeDomainSettings.Properties.VariableNames);
    chan_col_names = aligned_data.DeviceSettings.Left.timeDomainSettings.Properties.VariableNames(chan_tag_inds);
    left_chan_names = cellfun(@(x) aligned_data.DeviceSettings.Left.timeDomainSettings.(x){end},chan_col_names,'UniformOutput',false);
    
    left_spect = {};
    left_spect_freq = {};
    left_spect_time = {};
    left_PSD = {};
    remove_ind = [];
    for i = 1:length(left_chan_names)
        same_chan = cellfun(@(x) strcmp(left_chan_names{i},x),left_chan_names(1:i-1));
        if sum(same_chan) == 0  && ~any(isnan(aligned_data.left_LFP_table.(['key',num2str(i-1)]))) % Not a duplicate channel recording
            [data,time] = addEmptyData(aligned_data.left_taxis,aligned_data.left_LFP_table.(['key',num2str(i-1)]),left_sr,gapFillType);
            [left_spect{end+1},left_spect_freq{end+1},left_spect_time{end+1},left_PSD{end+1}]=spectrogram(data,WINDOW,NOVERLAP,nfft,left_sr);
        else
            remove_ind = [remove_ind,i];
        end
    end
    left_chan_names(remove_ind) = [];
    
    STFT.Left.Values = left_spect;
    STFT.Left.Time = left_spect_time;
    STFT.Left.Freq_Values = left_spect_freq;
    STFT.Left.PSD = left_PSD;
    STFT.Left.Chan_Names = left_chan_names;
end

%% Right RCS
if isfield(aligned_data,'right_LFP_table')
    % Spectrogram hyperparameters
    right_sr = aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate(end);
    
    if isinf(windowLength)
        WINDOW = right_sr;
    elseif windowLength > 0 && windowLength <= 1
        WINDOW = right_sr*windowLength;
        WINDOW = WINDOW + mod(WINDOW,2);
    else
        WINDOW = right_sr*windowLength;
        WINDOW = WINDOW + mod(WINDOW,2);
    end
    
    NOVERLAP = round(WINDOW*percentOverlap);
    
    if isempty(nfft)
        nfft = 2^nextpow2(right_sr);
    end
    
    chan_tag_inds = cellfun(@(x) contains(x,'chan'),aligned_data.DeviceSettings.Right.timeDomainSettings.Properties.VariableNames);
    chan_col_names = aligned_data.DeviceSettings.Right.timeDomainSettings.Properties.VariableNames(chan_tag_inds);
    right_chan_names = cellfun(@(x) aligned_data.DeviceSettings.Right.timeDomainSettings.(x){end},chan_col_names,'UniformOutput',false);
    
    right_spect = {};
    right_spect_freq = {};
    right_spect_time = {};
    right_PSD = {};
    remove_ind = [];
    for i = 1:length(right_chan_names)
        same_chan = cellfun(@(x) strcmp(right_chan_names{i},x),right_chan_names(1:i-1));
        if sum(same_chan) == 0  && ~any(isnan(aligned_data.left_LFP_table.(['key',num2str(i-1)]))) % Not a duplicate channel recording
            [data,time] = addEmptyData(aligned_data.right_taxis,aligned_data.right_LFP_table.(['key',num2str(i-1)]),right_sr,gapFillType);
            [right_spect{end+1},right_spect_freq{end+1},right_spect_time{end+1},right_PSD{end+1}]=spectrogram(data,WINDOW,NOVERLAP,nfft,right_sr);
        else
            remove_ind = [remove_ind,i];
        end
    end
    right_chan_names(remove_ind) = [];
    
    STFT.Right.Values = right_spect;
    STFT.Right.Time = right_spect_time;
    STFT.Right.Freq_Values = right_spect_freq;
    STFT.Right.PSD = right_PSD;
    STFT.Right.Chan_Names = right_chan_names;
end