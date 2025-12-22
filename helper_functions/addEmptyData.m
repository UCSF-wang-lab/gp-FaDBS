function varargout = addEmptyData(time_vec,LFP_signal,sampling_freq,type)
out_data = LFP_signal;
out_time = time_vec;

% Find time gaps
time_diff = diff(time_vec);
ind = find(time_diff > 1/sampling_freq + 1e-6);

cum_new_val_added = 0;
new_data_vals = cell(1,length(ind));
new_time_vals = cell(1,length(ind));
for i = 1:length(ind)
    new_time_vals{i} = (time_vec(ind(i))+1/sampling_freq:1/sampling_freq:time_vec(ind(i)+1)-1/sampling_freq)';
    if strcmp(type,'blank')
        new_data_vals{i} = zeros(length(new_time_vals{i}),1);
    end
    
    out_time = [out_time(1:ind(i)+cum_new_val_added);new_time_vals{i};out_time(ind(i)+cum_new_val_added+1:end)];
    out_data = [out_data(1:ind(i)+cum_new_val_added);new_data_vals{i};out_data(ind(i)+cum_new_val_added+1:end)];
    
    cum_new_val_added = cum_new_val_added + length(new_time_vals{i});
end

varargout{1} = out_data;
varargout{2} = out_time;
end