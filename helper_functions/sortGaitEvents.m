function sorted_gait_events = sortGaitEvents(gait_events,cycle_start_event)
if strcmp(gait_events.Properties.VariableNames{1},cycle_start_event)
    sorted_gait_events = gait_events;
else
    gait_event_order = [];
    switch cycle_start_event
        case 'LHS'
            gait_event_order = {'LHS','RTO','RHS','LTO'};
        case 'LTO'
            gait_event_order = {'LTO','LHS','RTO','RHS'};
        case 'RHS'
            gait_event_order = {'RHS','LTO','LHS','RTO'};
        case 'RTO'
            gait_event_order = {'RTO','RHS','LTO','LHS'};
    end
    
    shift_ind = find(cellfun(@(x) strcmp(x,cycle_start_event),gait_events.Properties.VariableNames))-1;
    
    sorted_gait_events = nan(1+height(gait_events),4);
    for i = 1:length(gait_event_order)-shift_ind
        sorted_gait_events(2:end,i) = gait_events.(gait_event_order{i});
    end
    
    for j = length(gait_event_order)-(shift_ind-1):length(gait_event_order)
        sorted_gait_events(1:end-1,j) = gait_events.(gait_event_order{j});
    end

    end_of_table = false;
    count = 1;
    while ~end_of_table
        for k = 2:size(sorted_gait_events,2)
            event_diff = sorted_gait_events(count,k)-sorted_gait_events(count,k-1);
            if event_diff > 1.5
                X = sorted_gait_events(1:count,:);
                Y = sorted_gait_events(count+1:end,:);
                Z = X(end,:);

                X(end,k:end) = nan;
                Z(1:k-1) = nan;

                sorted_gait_events = [X;Z;Y];
            end
        end
        if count + 1 > size(sorted_gait_events,1)
            end_of_table = true;
        end
        count = count + 1;
    end
    
    sorted_gait_events = array2table(sorted_gait_events,'VariableNames',gait_event_order);
end

% Remove any lines with all NAN
remove_ind = [];
for i = 1:height(sorted_gait_events)
    if sum(isnan(sorted_gait_events{i,:})) == 4
        remove_ind = [remove_ind,i];
    end
end
sorted_gait_events(remove_ind,:) = [];
end