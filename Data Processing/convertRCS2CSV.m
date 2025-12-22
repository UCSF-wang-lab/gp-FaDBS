function convertRCS2CSV(file_name,save_path,save_name)
% Processes the aligned data files from RC+S and converts them to a csv to
% use with Tanner Dixon's RCS simulation code
% struct.
% Author:   Kenneth Louie
% Date:     10-24-2022

% A filepath to the aligned data file needs to be passed in
% If it exist load in the file.
if ~exist('file_name','var') || isempty(file_name)
    error('File to convert must be provided.')
else
    if isempty(fileparts(file_name))
        error('Filename needs to be the entire path including directories.');
    else
        load(file_name);
    end
end

[left_data,right_data] = simplifyRCSDataTable(aligned_data);
[left_settings,right_settings] = extractRCSDeviceSettings(aligned_data);
% gait_events = aligned_data.gait_events;

% Save CSV files to save path
if ~exist('save_name','var') || isempty(save_name)
    % Default saving file name
    save_filename = generateFileName(file_name,aligned_data);
end

if ~exist('save_path','var') || isempty(save_path)
    save_path = fileparts(file_name);
end

save_name = fullfile(save_path,save_filename);

if ~isempty(left_data)
    writetable(left_data,[save_name,'LEFT.csv']);
end

if ~isempty(right_data)
    writetable(right_data,[save_name,'RIGHT.csv']);
end

if ~isempty(left_settings)
    writetable(left_settings,[save_name,'LEFT_SETTINGS.csv']);
end

if ~isempty(right_settings)
    writetable(right_settings,[save_name,'RIGHT_SETTINGS.csv']);
end

% if ~isempty(gait_events)
%     writetable(gait_events,[save_name,'Gait_Events.csv']);
% end

end

function varargout = simplifyRCSDataTable(data_table)
% Processes the aligned data files from RC+S and converts them to a csv to
% use with Tanner Dixon's RCS simulation code.


% Create seperate tables for left and right sides if detected
if sum(contains(fields(data_table),'left')) > 0
    out_table_left = table(data_table.left_taxis,...
        data_table.left_LFP_table.key0,...
        data_table.left_LFP_table.key1,...
        data_table.left_LFP_table.key2,...
        data_table.left_LFP_table.key3,...
        'VariableNames',{'time','key0','key1','key2','key3'});
else
    out_table_left = [];
end

if sum(contains(fields(data_table),'right')) > 0
    out_table_right = table(data_table.right_taxis,...
        data_table.right_LFP_table.key0,...
        data_table.right_LFP_table.key1,...
        data_table.right_LFP_table.key2,...
        data_table.right_LFP_table.key3,...
        'VariableNames',{'time','key0','key1','key2','key3'});
else
    out_table_right = [];
end

varargout{1} = out_table_left;
varargout{2} = out_table_right;

end


function varargout = extractRCSDeviceSettings(data)
% Extract device settings from the aligned data

if isfield(data.DeviceSettings,'Left')
    setting_table_left = getSettingTable(data.DeviceSettings.Left);
else
    setting_table_left = [];
end

if isfield(data.DeviceSettings,'Right')
    setting_table_right = getSettingTable(data.DeviceSettings.Right);
else
    setting_table_right = [];
end

varargout{1} = setting_table_left;
varargout{2} = setting_table_right;
end

function out_table = getSettingTable(device_settings,adaptive_settings)
out_table = [];

% amplitude gains for the recording pairs
amp_fields = fields(device_settings.metaData.ampGains);
temp = [];
for i = 1:length(amp_fields)
    temp = [temp,num2str(device_settings.metaData.ampGains.(amp_fields{i})),'|'];
end
temp = temp(1:end-1);
amp_gains{1} = temp;

% sampling frequency
fs_td = unique(device_settings.timeDomainSettings.samplingRate);

% fft size, update interval, bit shift
fft_size = device_settings.fftSettings.fftConfig.size;
fft_interval = device_settings.fftSettings.fftConfig.interval;

% UGLYYYYYYYYYY
if contains(device_settings.fftSettings.fftConfig.bandFormationConfig,'0')
    fft_bitshift = 0;
elseif contains(device_settings.fftSettings.fftConfig.bandFormationConfig,'1')
    fft_bitshift = 1;
elseif contains(device_settings.fftSettings.fftConfig.bandFormationConfig,'2')
    fft_bitshift = 2;
elseif contains(device_settings.fftSettings.fftConfig.bandFormationConfig,'3')
    fft_bitshift = 3;
elseif contains(device_settings.fftSettings.fftConfig.bandFormationConfig,'4')
    fft_bitshift = 4;
elseif contains(device_settings.fftSettings.fftConfig.bandFormationConfig,'5')
    fft_bitshift = 5;
elseif contains(device_settings.fftSettings.fftConfig.bandFormationConfig,'6')
    fft_bitshift = 6;
else
    fft_bitshift = 7;
end

% power band info
if ~isempty(device_settings.powerSettings)
    temp = [device_settings.powerSettings.powerBands(1).lowerBound';device_settings.powerSettings.powerBands(1).upperBound'];
    temp = temp(:)';      % slick way to interleave values
    temp = sprintf('%.4f|',temp);
    temp = temp(1:end-1);
    power_bands{1} = temp;
else
    power_bands{1} = [];
end

% adaptive settings
if exist('adaptive_settings','var')
    fprintf('cool');
else
    update_rate = nan;
    subtract_vector = nan;
    multiply_vector = nan;
    dual_threshold = nan;
    threshold = nan;
    onset = nan;
    termination = nan;
    blank_duration = nan;
    blank_both = nan;
    target_amp = nan;
    rise_time = 999;
    fall_time = 999;
end

out_table = table(amp_gains,fs_td,fft_size,fft_interval,fft_bitshift,power_bands,...
    update_rate,subtract_vector,multiply_vector,dual_threshold, threshold,onset,...
    termination,blank_duration,blank_both,target_amp,rise_time,fall_time,...
    'VariableNames',{'amp_gains','fs_td','FFT_size','FFT_interval','FFT_bitshift','power_bands',...
    'update_rate','subtract_vector','multiply_vector','dual_threshold','threshold','onset',...
    'termination','blank_duration','blank_both','target_amp','rise_time','fall_time'});
end

function file_name = generateFileName(loaded_file_name,data)
parts = strsplit(loaded_file_name,'/');
file_name = parts{end};
file_name = strrep(file_name,'w_Gait_Events_Ken.mat','');
end
