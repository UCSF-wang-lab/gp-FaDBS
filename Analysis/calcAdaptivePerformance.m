function performanceOutput = calcAdaptivePerformance(data,phase,hemisphere_true_state,save_output)
% Input variables
%   phase [=] The gait phase of interest for the detector to work on.
%             Current valid options: cSwing (contralateral swing)
%
%   hemisphere_true_state [=] A 2 x 1 cell array that indicates the true
%                             state for the desired gait phase of interest.
%                             First input is the true state for the left
%                             hemisphere and the second is for the right
%                             hemisphere. Put NaN to ignore a hemisphere.

if ~exist("save_output",'var') || isempty(save_output)
    save_output = false;
end

% Variables to hold confusion matrix values
left_cMatrix = [0,0];
right_cMatrix = [0,0];

if length(data)==1 % Passed in a single structure of aligned data
    left_cMatrix_byFile = zeros([1,2]);
    right_cMatrix_byFile = zeros([1,2]);
else % Passed in a list of files to parse into data
    left_cMatrix_byFile = zeros([length(data),2]);
    right_cMatrix_byFile = zeros([length(data),2]);
end

for i = 1:length(data)
    % Load data
    if length(data)==1
        aligned_data = data;
    else
        load(fullfile(data(i).folder,data(i).name));
    end

    % Adjust gait events so you don't have to worry about event wrapping
    if isfield(aligned_data,'left_adaptive_taxis')
        gel = adjustGaitEvents(aligned_data.gait_events,phase,'left');
        
        % Average gait phase duration. Use this value to fill in missing
        % gait events if number of marked gait events is >70%
        if (sum(isnan(gel{:,1})) + sum(isnan(gel{:,2})))/(2*height(gel)) < 0.3
            mean_gp_dur = mean(gel{:,2}-gel{:,1},'omitnan');

            for j= 1:height(gel)
                if isnan(gel{j,1}) && ~isnan(gel{j,2})
                    gel{j,1} = gel{j,2}-mean_gp_dur;
                elseif ~isnan(gel{j,1}) && isnan(gel{j,2})
                    gel{j,2} = gel{j,1} + mean_gp_dur;
                end
            end
        end
    else
        gel = nan;
    end

    if isfield(aligned_data,'right_adaptive_taxis')
        ger = adjustGaitEvents(aligned_data.gait_events,phase,'right');

        % Average gait phase duration.
        if (sum(isnan(ger{:,1})) + sum(isnan(ger{:,2})))/(2*height(ger)) < 0.3
            mean_gp_dur = mean(ger{:,2}-ger{:,1},'omitnan');

            for j= 1:height(ger)
                if isnan(ger{j,1}) && ~isnan(ger{j,2})
                    ger{j,1} = ger{j,2}-mean_gp_dur;
                elseif ~isnan(ger{j,1}) && isnan(ger{j,2})
                    ger{j,2} = ger{j,1} + mean_gp_dur;
                end
            end
        end
    else
        ger = nan;
    end

    % Calculate confusion matrix for this file
    cMatrix_file = calcAccuracyTimings(aligned_data,gel,ger,hemisphere_true_state);

    % Add cMatrix to left and right side appropriately
    if isfield(cMatrix_file,'Left')
        left_cMatrix = left_cMatrix + cMatrix_file.Left;
        left_cMatrix_byFile(i,:) = cMatrix_file.Left;
        left_baseline_byFile(i) = cMatrix_file.Left_baseline_accuracy;
    end

    if isfield(cMatrix_file,'Right')
        right_cMatrix = right_cMatrix + cMatrix_file.Right;
        right_cMatrix_byFile(i,:) = cMatrix_file.Right;
        right_baseline_byFile(i) = cMatrix_file.Left_baseline_accuracy;
    end
end

if exist('left_cMatrix','var')
    performanceOutput.Left.Durations = left_cMatrix;
    performanceOutput.Left_by_file.Durations = left_cMatrix_byFile;
    performanceOutput.Left_baseline_by_file = left_baseline_byFile;

    performanceOutput.Left.Percentages = left_cMatrix./sum(left_cMatrix,2);
    performanceOutput.Left_by_file.Percentages = left_cMatrix_byFile./sum(left_cMatrix_byFile,2);
end

if exist('right_cMatrix','var')
    performanceOutput.Right.Durations = right_cMatrix;
    performanceOutput.Right_by_file.Durations = right_cMatrix_byFile;
    performanceOutput.Right_baseline_by_file = right_baseline_byFile;

    performanceOutput.Right.Percentages = right_cMatrix./sum(right_cMatrix,2);
    performanceOutput.Right_by_file.Percentages = right_cMatrix_byFile./sum(right_cMatrix_byFile,2);
end

% Create table
variableNames = {'Side','TP (dur)','TN (dur)','Accuracy'};
performanceTable = table({'Left';'Right'},...
    [left_cMatrix(1);right_cMatrix(1)],[left_cMatrix(2);right_cMatrix(2)],...
    [(left_cMatrix(1)+left_cMatrix(2))/sum(left_cMatrix);(right_cMatrix(1)+right_cMatrix(2))/sum(right_cMatrix)],...
    'VariableNames',variableNames);

performanceOutput.Table = performanceTable;

if save_output
    % create output variable name
    if length(data)==1
        % NOTHING RIGHT NOW
    else
        [folder,filename] = fileparts(fullfile(data(1).folder,data(1).name));
    end

    filename_parts = strsplit(filename,'_');
    aDBS_visit_number = filename_parts{find(~cellfun(@isempty,regexp(filename_parts,'aDBS[0-9][0-9][0-9]')))};
    setting_num_ind = find(contains(filename_parts,'Setting'));
    setting_str = [filename_parts{setting_num_ind},'_',filename_parts{setting_num_ind+1}];
    root_folder = fileparts(folder);
    root_folder = [root_folder,'/Analysis Data/aDBS'];
    
    save_filename = fullfile(root_folder,[aDBS_visit_number,'_',setting_str,'_Accuracies.csv']);
    writetable(performanceTable,save_filename);
end

end


function gait_events = adjustGaitEvents(gait_events_table,phase,hemisphere)

% depending on the phase type, pick a different starting gait phase event.
% will add more phases in the future
switch phase
    case 'cSwing'
        if strcmpi(hemisphere,'left')
            gait_events = sortGaitEvents(gait_events_table,'RTO');
        elseif strcmpi(hemisphere,'right')
            gait_events = sortGaitEvents(gait_events_table,'LTO');
        end
end
end



function cMatrix = calcAccuracyTimings(aligned_data,gel,ger,hemisphere_true_state)
% Setup output variable
% Order: [true positive,true negative]
% durations
left_dur = [0,0];
right_dur = [0,0];

% Calc left hemisphere state accuracy
if istable(gel)
    % Grab state detector current state and calculate if the first data
    % point is in contralateral swing phase
    A = aligned_data.left_adaptive_taxis(1) >= gel{:,1};
    B = aligned_data.left_adaptive_taxis(1) <= gel{:,2};
    C = find(A&B,1);

    if ~isempty(C)
        prev_gEvent = 1;
    else
        prev_gEvent = 0;
    end
    prev_state = aligned_data.left_adaptive_table.CurrentAdaptiveState(1);

    % Loop through the rest of the data points
    for i = 2:height(aligned_data.left_adaptive_taxis)
        A = aligned_data.left_adaptive_taxis(i) >= gel{:,1};
        B = aligned_data.left_adaptive_taxis(i) <= gel{:,2};
        C = find(A&B,1);

        % curr data point is within a gait phase of interest
        if ~isempty(C)  
            if prev_gEvent == 1 % Previous time point was in gait event of interest
                if strcmp(prev_state,hemisphere_true_state{1})
                    left_dur(1) = left_dur(1) + (aligned_data.left_adaptive_taxis(i)-aligned_data.left_adaptive_taxis(i-1));
                end
            elseif prev_gEvent == 0 % Previous time point was not in gait event of interest
                gait_event_start = gel{C,1};
                time_epoch_start_2_event = gait_event_start-aligned_data.left_adaptive_taxis(i-1);
                time_event_2_epoch_end = aligned_data.left_adaptive_taxis(i)-gait_event_start;
                if strcmp(prev_state,hemisphere_true_state{1})
                    left_dur(1) = left_dur(1) + time_event_2_epoch_end;
                elseif ~strcmp(prev_state,hemisphere_true_state{1})
                    left_dur(2) = left_dur(2) + time_epoch_start_2_event;
                end
            end

            % Update prev state variables
            prev_gEvent = 1;
            prev_state = aligned_data.left_adaptive_table.CurrentAdaptiveState(i);
        else    % data point is not within a gait phase of interest
            if prev_gEvent == 1
                A = aligned_data.left_adaptive_taxis(i-1) >= gel{:,1};
                B = aligned_data.left_adaptive_taxis(i-1) <= gel{:,2};
                C = find(A&B,1);
                gait_event_start = gel{C,2};
                time_epoch_start_2_event = gait_event_start-aligned_data.left_adaptive_taxis(i-1);
                time_event_2_epoch_end = aligned_data.left_adaptive_taxis(i)-gait_event_start;
                if strcmp(prev_state,hemisphere_true_state{1})
                    left_dur(1) = left_dur(1) + time_epoch_start_2_event;
                elseif ~strcmp(prev_state,hemisphere_true_state{1})
                    left_dur(2) = left_dur(2) + time_event_2_epoch_end;
                end
            elseif prev_gEvent == 0
                if ~strcmp(prev_state,hemisphere_true_state{1})
                    left_dur(2) = left_dur(2) + (aligned_data.left_adaptive_taxis(i)-aligned_data.left_adaptive_taxis(i-1));
                end
            end

            % Update prev state variables
            prev_gEvent = 0;
            prev_state = aligned_data.left_adaptive_table.CurrentAdaptiveState(i);
        end
    end
end

% Calc right hemisphere state accuracy
if istable(ger)
    % Grab state detector current state and calculate if the first data
    % point is in contralateral swing phase
    A = aligned_data.right_adaptive_taxis(1) >= ger{:,1};
    B = aligned_data.right_adaptive_taxis(1) <= ger{:,2};
    C = find(A&B,1);

    if ~isempty(C)
        prev_gEvent = 1;
    else
        prev_gEvent = 0;
    end
    prev_state = aligned_data.right_adaptive_table.CurrentAdaptiveState(1);
    
    % Loop through the rest of the data points
    for i = 2:height(aligned_data.right_adaptive_taxis)
        A = aligned_data.right_adaptive_taxis(i) >= ger{:,1};
        B = aligned_data.right_adaptive_taxis(i) <= ger{:,2};
        C = find(A&B,1);

        % curr data point is within a gait phase of interest
        if ~isempty(C)  
            if prev_gEvent == 1 % Previous time point was in gait event of interest
                if strcmp(prev_state,hemisphere_true_state{2})
                    right_dur(1) = right_dur(1) + (aligned_data.right_adaptive_taxis(i)-aligned_data.right_adaptive_taxis(i-1));
                end
            elseif prev_gEvent == 0 % Previous time point was not in gait event of interest
                gait_event_start = ger{C,1};
                time_epoch_start_2_event = gait_event_start-aligned_data.right_adaptive_taxis(i-1);
                time_event_2_epoch_end = aligned_data.right_adaptive_taxis(i)-gait_event_start;
                if strcmp(prev_state,hemisphere_true_state{2})
                    right_dur(1) = right_dur(1) + time_event_2_epoch_end;
                elseif ~strcmp(prev_state,hemisphere_true_state{2})
                    right_dur(2) = right_dur(2) + time_epoch_start_2_event;
                end
            end

            % Update prev state variables
            prev_gEvent = 1;
            prev_state = aligned_data.right_adaptive_table.CurrentAdaptiveState(i);
        else    % data point is not within a gait phase of interest
            if prev_gEvent == 1
                A = aligned_data.right_adaptive_taxis(i-1) >= ger{:,1};
                B = aligned_data.right_adaptive_taxis(i-1) <= ger{:,2};
                C = find(A&B,1);
                gait_event_start = ger{C,2};
                time_epoch_start_2_event = gait_event_start-aligned_data.right_adaptive_taxis(i-1);
                time_event_2_epoch_end = aligned_data.right_adaptive_taxis(i)-gait_event_start;
                if strcmp(prev_state,hemisphere_true_state{2})
                    right_dur(1) = right_dur(1) + time_epoch_start_2_event;
                elseif ~strcmp(prev_state,hemisphere_true_state{2})
                    right_dur(2) = right_dur(2) + time_event_2_epoch_end;
                end
            elseif prev_gEvent == 0
                if ~strcmp(prev_state,hemisphere_true_state{2})
                    right_dur(2) = right_dur(2) + (aligned_data.right_adaptive_taxis(i)-aligned_data.right_adaptive_taxis(i-1));
                end
            end

            % Update prev state variables
            prev_gEvent = 0;
            prev_state = aligned_data.right_adaptive_table.CurrentAdaptiveState(i);
        end
    end
end

% calculate baseline accuracy
if istable(gel)
    left_baseline = 0;

    A = find(gel{:,1}>=aligned_data.left_adaptive_taxis(1),1,'first');
    B = find(gel{:,2}>=aligned_data.left_adaptive_taxis(1),1,'first');
    C = find(gel{:,3}>=aligned_data.left_adaptive_taxis(1),1,'first');
    D = find(gel{:,4}>=aligned_data.left_adaptive_taxis(1),1,'first');
    table_start = min([A,B,C,D]);

    A = find(gel{:,1}<=aligned_data.left_adaptive_taxis(end),1,'last');
    B = find(gel{:,2}<=aligned_data.left_adaptive_taxis(end),1,'last');
    C = find(gel{:,3}<=aligned_data.left_adaptive_taxis(end),1,'last');
    D = find(gel{:,4}<=aligned_data.left_adaptive_taxis(end),1,'last');
    table_end = max([A,B,C,D]);

    for i = table_start:table_end
        if gel{i,1} >= 0 && ~isnan(gel{i,2})
            left_baseline = left_baseline + (gel{i,2}-gel{i,1});
        elseif gel{i,1} < 0 && gel{i,2} > 0
            left_baseline = left_baseline + gel{i,2};
        end
    end
end

if istable(ger)
    right_baseline = 0;

    A = find(ger{:,1}>=aligned_data.right_adaptive_taxis(1),1,'first');
    B = find(ger{:,2}>=aligned_data.right_adaptive_taxis(1),1,'first');
    C = find(ger{:,3}>=aligned_data.right_adaptive_taxis(1),1,'first');
    D = find(ger{:,4}>=aligned_data.right_adaptive_taxis(1),1,'first');
    table_start = min([A,B,C,D]);

    A = find(ger{:,1}<=aligned_data.right_adaptive_taxis(end),1,'last');
    B = find(ger{:,2}<=aligned_data.right_adaptive_taxis(end),1,'last');
    C = find(ger{:,3}<=aligned_data.right_adaptive_taxis(end),1,'last');
    D = find(ger{:,4}<=aligned_data.right_adaptive_taxis(end),1,'last');
    table_end = max([A,B,C,D]);

    for i = table_start:table_end
        if ger{i,1} >= 0 && ~isnan(ger{i,2})
            right_baseline = right_baseline + (ger{i,2}-ger{i,1});
        elseif ger{i,1} < 0 && ger{i,2} > 0
            right_baseline = right_baseline + ger{i,2};
        end
    end
end

% create output
if sum(left_dur == 0) ~= 2
    cMatrix.Left = left_dur;
    cMatrix.Left_baseline_accuracy = left_baseline/aligned_data.left_adaptive_taxis(end);
end

if sum(right_dur == 0) ~= 2
    cMatrix.Right = right_dur;
    cMatrix.Right_baseline_accuracy = right_baseline/aligned_data.right_adaptive_taxis(end);
end
end