function [gait_metrics_table,cadence,gait_speed] = calcGaitMetrics(filename,varargin)
% Calculates spatial and temporal gait metrics from Delsys or Xsens data
%
% INPUTS:   
%           filesname           [=] Matrix of gait events. Must be 4
%
% OPTIONAL INPUTS:
%
%           level_type          [=] "none"/"single". Changes the way gait
%                                   events are automatically detected and
%                                   how gait metrics are calculated.
%                                   Default is "single".
%           visit_name          [=] String to denote visit type. Examples:
%                                   preop, preprogramming, 
%                                   preprogrammingUnilat,
%                                   preprogrammingBilat, dbsOpt,
%                                   dbsOptUnilat,dbsOptBilat,aDBS.
%
%           trial_num           [=] Trial number if this task is repeated.
%           
%           med_state           [=] String input to denote medication state
%                                   of data. Can be on/off/low/na. Default
%                                   is na.
%
%           stim_state          [=] String input to denote stimulation
%                                   state. Can be clinical/aDBS/off/na.
%
% OUTPUTS:  gait_metrics_table  [=] Table of gait events by gait cycle.
%
% Author:   Kenneth Louie
% Date:     06/27/2023
%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nargin/2
    switch varargin{i*2-1}
        case 'xsens_filename'
            xsens_filename = varargin{i*2};
        case 'delsys_filename'
            delsys_filename = varargin{i*2};
        case 'gait_event_filename'
            gait_event_filename = varargin{i*2};
        case 'level_type'
            level_type = varargin{i*2};
        case 'visit_name'
            visit_name = varargin{i*2};
        case 'med_state'
            med_state = varargin{i*2};
        case 'stim_state'
            stim_state = varargin{i*2};
        case 'trial_num'
            trial_num = varargin{i*2};
        case 'subjectID'
            subjectID = varargin{i*2};
        case 'save_data'
            save_data = varargin{i*2};
    end
    
end

%% Load data
if ~isempty(filename)
    load(filename);
    xsens_data = aligned_data.Xsens;
    delsys_data = aligned_data.Delsys;
    gait_events = aligned_data.gait_events;
    
    if isfield(aligned_data,'med_condition')
        med_state = aligned_data.med_condition;
    else
        med_state = 'na';
    end
    
    if isfield(aligned_data,'stim_condition')
        stim_state = aligned_data.stim_condition;
    else
        stim_state = 'na';
    end

    if isfield(aligned_data,'trial_num')
        trial_num = aligned_data.trial_num;
    end
end

if (~exist('xsens_data','var') || isempty(xsens_data)) && exist('xsens_filename','var')
    try
        xsens_data = readtable(xsens_filename);
    catch
        error('Xsens file does not exist, data is empty, or wrong filename.');
    end
end

if (~exist('delsys_data','var') || isempty(delsys_data)) && exist('delsys_filename','var')
    try
        temp_data = load(delsys_filename);
        delsys_data = temp_data.out_struct;
    catch
        error('Delsys file does not exist, data is empty, or wrong filename.');
    end
end

if (~exist('gait_events','var') || isempty(gait_events)) && exist('gait_event_filename','var')
    try
        gait_events = readtable(gait_event_filename);
    catch
        error('gait event file does not exist, data is empty, or wrong filename.');
    end
end

%% Set warnings if xsens dataset is missing
if ~exist('xsens_data','var') && exist('delsys_data','var')
    warning(sprintf('Gait Metric Calculation. \n No Xsens data detected. \n Delsys data detected \n Only temporal gait metrics will be calculated.'));
end

%% Default values
if ~exist('level_type','var') || isempty(level_type)
    level_type = 'single';
end

if ~exist('visit_name','var') || isempty(visit_name)
    visit_name = 'na';
end

if ~exist('med_state','var') || isempty(med_state)
    med_state = 'na';
end

if ~exist('stim_state','var') || isempty(stim_state)
    stim_state = 'na';
end

if ~exist('trial_num','var') || isempty(trial_num)
    trial_num = 1;
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'na';
end

if ~exist('adaptive_setting_num','var')
    adaptive_setting_num = nan;
end

if ~exist('save_data','var') || isempty(save_data)
    save_data = false;
end

%% Calculate gait matrics and build output table
% out_table = array2table(nan(n_gait_events*2,11));
% out_table.Properties.VariableNames = {'Gait_Cycle','Side',...
%     'Step_Time','Step_Length'};
% out_table.Gait_Cycle = repelem(1:n_gait_events,2)';
% out_table.Side = repmat(['R';'L'],n_gait_events,1);

if exist('xsens_data','var') && ~isempty(xsens_data)
    [gait_metrics_table] = getMetrics(gait_events,xsens_data,level_type,single_direction,direction_turn_time_range);
elseif ~exist('xsens_data','var') && exist('delsys_data','var') && ~isempty(delsys_data)
    [gait_metrics_table] = getMetrics(gait_events,delsys_data,level_type,single_direction,direction_turn_time_range);
end

SubjectID = repmat({subjectID},height(gait_metrics_table),1);
VisitName = repmat({visit_name},height(gait_metrics_table),1);
MedState = repmat({med_state},height(gait_metrics_table),1);
StimState = repmat({stim_state},height(gait_metrics_table),1);
TrialNum = repmat(trial_num,height(gait_metrics_table),1);

if ~isempty(adaptive_setting_num)
    SettingNum = repmat(adaptive_setting_num,height(gait_metrics_table),1);
    gait_metrics_table = addvars(gait_metrics_table,SubjectID,VisitName,SettingNum,MedState,StimState,TrialNum,'Before','GaitCycle');
else
    gait_metrics_table = addvars(gait_metrics_table,SubjectID,VisitName,MedState,StimState,TrialNum,'Before','GaitCycle');
end



if save_data
    % gait metrics table
    [folder,file] = fileparts(filename);
    folder = strrep(folder,"Aligned Data","Analysis Data/Gait_Metrics");
    if ~exist(folder,'dir')
        mkdir(folder);
    end

    idx = strfind(file,"w");
    file = [file(1:idx-1),'gait_metrics.csv'];

    writetable(gait_metrics_table,fullfile(folder,file));
end

end

function [gait_metrics_table] = getMetrics(gait_events,foot_data,level_type)
% [LHS,RTO,RHS,LTO]
% step length, width, time
% stride length, time

if istable(foot_data)
    % Order gait events
    gait_events_turns_removed = removeGaitCyclesTurns(foot_data,gait_events,30);
    gait_events_ordered = sortGaitEvents(gait_events_turns_removed,'LHS');
    
    % Determine which gait cycles to consider
    gc_start_search = find(gait_events_ordered.LHS>=foot_data.Time(1),1,'first');
    gc_end_search = find(gait_events_ordered.LHS<=foot_data.Time(end),1,'last');
    gait_events_ordered_trim = gait_events_ordered(gc_start_search:gc_end_search,:);
end

if ~isempty(gait_events_ordered_trim)
    n_gait_events = size(gait_events_ordered_trim,1);
    step_length = nan(n_gait_events*2,1);
    step_time = nan(n_gait_events*2,1);

    if strcmp(level_type,'none')
    end

    if strcmp(level_type,'single')
        % [LHS,RTO,RHS,LTO]
        % Step metrics using left leg as the line of progression and vector projections
        for i = 1:height(gait_events_ordered_trim)-1
            if (~isnan(gait_events_ordered_trim.RHS(i)) && ~isnan(gait_events_ordered_trim.LHS(i+1))) && ((gait_events_ordered_trim.LHS(i+1)-gait_events_ordered_trim.RHS(i))<1)
                [~,left_foot_pos1_ind] = min(abs(foot_data.Time-gait_events_ordered_trim.LHS(i)));
                [~,right_foot_pos1_ind] = min(abs(foot_data.Time-gait_events_ordered_trim.RHS(i)));
                [~,left_foot_pos2_ind] = min(abs(foot_data.Time-gait_events_ordered_trim.LHS(i+1)));

                left_foot_pos1 = [foot_data.LeftFoot_PosX(left_foot_pos1_ind),foot_data.LeftFoot_PosY(left_foot_pos1_ind)];
                right_foot_pos1 = [foot_data.RightFoot_PosX(right_foot_pos1_ind),foot_data.RightFoot_PosY(right_foot_pos1_ind)];
                left_foot_pos2 = [foot_data.LeftFoot_PosX(left_foot_pos2_ind),foot_data.LeftFoot_PosY(left_foot_pos2_ind)];

                % Right step
                left_foot_seg = left_foot_pos2-left_foot_pos1;
                left_right_foot_seg = right_foot_pos1-left_foot_pos1;
                projected_point = left_foot_pos1 + dot(left_right_foot_seg,left_foot_seg)/dot(left_foot_seg,left_foot_seg)*left_foot_seg;
                step_length(i*2) = sqrt(sum((projected_point-left_foot_pos1).^2));

                % Left Step length
                step_length((i*2)-1) = stride_length((i*2)-1)-step_length(i*2);
            end
        end

        % Step metrics using pelvis as the line of progression

        % Temporal metrics
        % Left
        step_time(1:2:(n_gait_events-1)*2) = gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.RHS(1:end-1);     % Left step time = time from right heel strike to left heel-strike
        
        % Right
        step_time(2:2:(n_gait_events-1)*2) = gait_events_ordered_trim.RHS(1:end-1)-gait_events_ordered_trim.LHS(1:end-1);     % Left step time = time from right heel strike to left heel-strike
    end

    % Filter values that are outside what is normal
    step_length(step_length>0.70*mean(stride_length,'omitnan') | step_length < 0) = nan;
    step_time(step_time>1 | step_time < 0) = nan;

    side = repmat({'L';'R'},n_gait_events,1);
    gait_cycle = repelem(1:n_gait_events,2)';

    % Create table
    gait_metrics_table = table(gait_cycle,side,step_length,step_time,stride_length,stride_time,swing_time,stance_time,dst,'VariableNames',{'GaitCycle','Side','StepLength','StepTime'});
else
    error("No Gait Events within this file.");
end
end