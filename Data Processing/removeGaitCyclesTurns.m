function gait_events_turns_removed = removeGaitCyclesTurns(xsens_data,gait_events,filter_threshold)
% Check to see if last row of xsens data is nan
if isnan(xsens_data{end,2})
    xsens_data(end,:) = [];
end

% Filter and convert to degrees
[b,a]=butter(4,(1.5/30));
pelvis_data = abs(filtfilt(b,a,xsens_data.Pelvis_angVelZ)*57.3);

% Find all threshold crossings in both directions
pos_ind = find((pelvis_data(2:end) > filter_threshold) & (pelvis_data(1:end-1) < filter_threshold));
neg_ind = find((pelvis_data(2:end) < filter_threshold) & (pelvis_data(1:end-1) > filter_threshold));

% Filter through the positive and negative threshold crossings to remove
% erroneous crossings
n_ind_check = 60;
min_cross_dist = 150;
remove_ind = [];
for i = 1:length(pos_ind)
    if pos_ind(i)+n_ind_check <= length(pelvis_data)
        if min(pelvis_data(pos_ind(i)+1:pos_ind(i)+n_ind_check))<filter_threshold
            remove_ind(end+1) = i;
        end
    else
        if mean(pelvis_data(pos_ind(i)+1:end))<filter_threshold
            remove_ind(end+1) = i;
        end
    end
end
pos_ind(remove_ind) = [];
remove_ind = find(diff(pos_ind)<min_cross_dist);
pos_ind(remove_ind+1) = [];

remove_ind = [];
for i = 1:length(neg_ind)
    if neg_ind(i) + n_ind_check <= length(pelvis_data)
        if max(pelvis_data(neg_ind(i):neg_ind(i)+n_ind_check))>filter_threshold
            if (sum(pelvis_data(neg_ind(i):neg_ind(i)+n_ind_check)>filter_threshold)/length(pelvis_data(neg_ind(i):neg_ind(i)+n_ind_check)))>0.3
                remove_ind(end+1) = i;
            end
        end
    else
        if max(pelvis_data(neg_ind(i):end))>filter_threshold
            if (sum(pelvis_data(neg_ind(i):end)>filter_threshold)/length(pelvis_data(neg_ind(i):end)))>0.3
                remove_ind(end+1) = i;
            end
        end
    end
end
neg_ind(remove_ind) = [];
remove_ind = find(diff(neg_ind)<min_cross_dist);
neg_ind(remove_ind+1) = [];

% Combine start and end times
turn_times = [];
for i = 1:min([length(pos_ind),length(neg_ind)])
    ind_diff = neg_ind-pos_ind(i);
    matching_neg_ind = find(ind_diff>0,1,'first');
    if ~isempty(matching_neg_ind)
        turn_times(i,1) = xsens_data.Time(pos_ind(i));
        turn_times(i,2) = xsens_data.Time(neg_ind(matching_neg_ind));
    end
end

remove_ind = [];
for i = 1:height(gait_events)
    for j = 1:size(turn_times,1)
        if any(gait_events{i,:}>=turn_times(j,1) & gait_events{i,:} <= turn_times(j,2))
            remove_ind(end+1) = i;
        end
    end
end

gait_events_turns_removed = gait_events;
gait_events_turns_removed(remove_ind,:) = [];

end