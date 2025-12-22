 %% Load data
P2_data = '/Users/USER/Documents/P2_dbs_opt_data.csv';
P3_data = '/Users/USER/Documents/P3_dbs_opt_data.csv';

folder_list = {P2_data,P3_data};
gait_event_order = {'LHS','RTO','RHS','LTO'};
colors = CBMap('SwingPhases');
                 
%% Process and Plot Data
fig_hand = figure('Name','Aggregate PSD');
leg_hand = [];
ax_hand = [];
inset_ax = [];
legend_items = [];

custom_xticks = [10,30,50];

for i = 1:length(folder_list)
    patient_files = dir(folder_list{i});
    PSD_gait_events.Left = {};
    PSD_gait_events.Right = {};
    
    for j = 1:length(patient_files)
        load(fullfile(patient_files(j).folder,patient_files(j).name));
        signal_analysis_data = calcRCS_STFT(aligned_data,[],1,0.9,[]);
        gait_events_sorted = sortGaitEvents(aligned_data.gait_events,'LTO');
        
        nleft = height(aligned_data.left_LFP_table);
        
        if isfield(aligned_data,'right_LFP_table')
            nright = height(aligned_data.right_LFP_table);
        else
            nright = [];
        end
        
        if ~isempty(nright)
            gait_event_start_ind = find(gait_events_sorted.LTO>min([aligned_data.left_taxis(1),aligned_data.right_taxis(1)]),1,'first');
            gait_event_end_ind = find(gait_events_sorted.RHS<max([aligned_data.left_taxis(nleft),aligned_data.right_taxis(nright)]),1,'last');
            gait_event_range = [gait_event_start_ind,gait_event_end_ind];
        else
            gait_event_start_ind = find(gait_events_sorted.LTO>aligned_data.left_taxis(1),1,'first');
            gait_event_end_ind = find(gait_events_sorted.RHS<aligned_data.left_taxis(nleft),1,'last');
            gait_event_range = [gait_event_start_ind,gait_event_end_ind];
        end
        
        if isfield(signal_analysis_data,'Left')
            for m = 1:length(signal_analysis_data.Left.Chan_Names)
                if ~isfield(PSD_gait_events.Left,'Left_Swing')
                    PSD_gait_events.Left.Left_Swing = cell(1,length(signal_analysis_data.Left.Chan_Names));
                    PSD_gait_events.Left.Right_Swing = cell(1,length(signal_analysis_data.Left.Chan_Names));
                end
                
                for n = gait_event_range(1):gait_event_range(2)
                    left_swing_event_start_time = gait_events_sorted.LTO(n);
                    left_swing_event_end_time = gait_events_sorted.LHS(n);
                    left_swing_event_end_time = left_swing_event_start_time + (left_swing_event_end_time-left_swing_event_start_time)/2;
                    if ~isnan(left_swing_event_start_time) && ~isnan(left_swing_event_end_time)
                        [~,min_ind] = min(abs(signal_analysis_data.Left.Time{m}-left_swing_event_start_time));
                        [~,max_ind] = min(abs(signal_analysis_data.Left.Time{m}-left_swing_event_end_time));
                        power_values = signal_analysis_data.Left.PSD{m}(:,min_ind:max_ind);
                        
                        if ~any(isinf(power_values)|isnan(power_values)|power_values==0, 'all')
                            PSD_gait_events.Left.Left_Swing{m}(:,end+1) = mean(power_values,2);
                        end
                    end
                    
                    right_swing_event_start_time = gait_events_sorted.RTO(n);
                    right_swing_event_end_time = gait_events_sorted.RHS(n);
                    right_swing_event_end_time = right_swing_event_start_time + (right_swing_event_end_time-right_swing_event_start_time)/2;
                    if ~isnan(right_swing_event_start_time) && ~isnan(right_swing_event_end_time)
                        [~,min_ind] = min(abs(signal_analysis_data.Left.Time{m}-right_swing_event_start_time));
                        [~,max_ind] = min(abs(signal_analysis_data.Left.Time{m}-right_swing_event_end_time));
                        power_values = signal_analysis_data.Left.PSD{m}(:,min_ind:max_ind);
                        
                        if ~any(isinf(power_values)|isnan(power_values)|power_values==0, 'all')
                            PSD_gait_events.Left.Right_Swing{m}(:,end+1) = mean(power_values,2);
                        end
                    end
                end
            end
        end
        
        if isfield(signal_analysis_data,'Right')
            for m = 1:length(signal_analysis_data.Right.Chan_Names)
                if ~isfield(PSD_gait_events.Right,'Left_Swing')
                    PSD_gait_events.Right.Left_Swing = cell(1,length(signal_analysis_data.Right.Chan_Names));
                    PSD_gait_events.Right.Right_Swing = cell(1,length(signal_analysis_data.Right.Chan_Names));
                end
                
                for n = gait_event_range(1):gait_event_range(2)
                    left_swing_event_start_time = gait_events_sorted.LTO(n);
                    left_swing_event_end_time = gait_events_sorted.LHS(n);
                    left_swing_event_end_time = left_swing_event_start_time + (left_swing_event_end_time-left_swing_event_start_time)/2;
                    if ~isnan(left_swing_event_start_time) && ~isnan(left_swing_event_end_time)
                        [~,min_ind] = min(abs(signal_analysis_data.Right.Time{m}-left_swing_event_start_time));
                        [~,max_ind] = min(abs(signal_analysis_data.Right.Time{m}-left_swing_event_end_time));
                        power_values = signal_analysis_data.Right.PSD{m}(:,min_ind:max_ind);
                        
                        if ~any(isinf(power_values)|isnan(power_values)|power_values==0, 'all')
                            PSD_gait_events.Right.Left_Swing{m}(:,end+1) = mean(power_values,2);
                        end
                    end
                    
                    right_swing_event_start_time = gait_events_sorted.RTO(n);
                    right_swing_event_end_time = gait_events_sorted.RHS(n);
                    right_swing_event_end_time = right_swing_event_start_time + (right_swing_event_end_time-right_swing_event_start_time)/2;
                    if ~isnan(right_swing_event_start_time) && ~isnan(right_swing_event_end_time)
                        [~,min_ind] = min(abs(signal_analysis_data.Right.Time{m}-right_swing_event_start_time));
                        [~,max_ind] = min(abs(signal_analysis_data.Right.Time{m}-right_swing_event_end_time));
                        power_values = signal_analysis_data.Right.PSD{m}(:,min_ind:max_ind);
                        
                        if ~any(isinf(power_values)|isnan(power_values)|power_values==0, 'all')
                            PSD_gait_events.Right.Right_Swing{m}(:,end+1) = mean(power_values,2);
                        end
                    end
                end
            end
        end
    end
    
    % Plot 
    figure(fig_hand);
    switch(i)
        case 1
            ax_hand(end+1) = subplot(1,3,1);
            plot(signal_analysis_data.Left.Freq_Values{1},mean(10*log10(PSD_gait_events.Left.Left_Swing{2}),2),'Color',colors.Left_Swing,'LineWidth',1.15,'DisplayName','Swing Phase Left');
            hold on;
            plot(signal_analysis_data.Left.Freq_Values{1},mean(10*log10(PSD_gait_events.Left.Right_Swing{2}),2),'Color',colors.Right_Swing,'LineWidth',1.15,'DisplayName','Swing Phase Right');
            title('GPi');
            ylabel('db/Hz');
            xlabel('Frequency (Hz)');
        case 2
            ax_hand(end+1) = subplot(1,3,2);
            plot(signal_analysis_data.Left.Freq_Values{1},mean(10*log10(PSD_gait_events.Right.Left_Swing{3}),2),'Color',colors.Left_Swing,'LineWidth',1.15,'DisplayName','Swing Phase Left');
            hold on;
            plot(signal_analysis_data.Left.Freq_Values{1},mean(10*log10(PSD_gait_events.Right.Right_Swing{3}),2),'Color',colors.Right_Swing,'LineWidth',1.15,'DisplayName','Swing Phase Right');
            title('M1');
            ylabel('db/Hz');
            xlabel('Frequency (Hz)');
            
            ax_hand(end+1) = subplot(1,3,3);
            leg_hand(1) = plot(signal_analysis_data.Left.Freq_Values{1},mean(10*log10(PSD_gait_events.Right.Left_Swing{4}),2),'Color',colors.Left_Swing,'LineWidth',1.15,'DisplayName','Swing Phase Left');
            hold on;
            leg_hand(2) = plot(signal_analysis_data.Left.Freq_Values{1},mean(10*log10(PSD_gait_events.Right.Right_Swing{4}),2),'Color',colors.Right_Swing,'LineWidth',1.15,'DisplayName','Swing Phase Right');
            title('PM');
            ylabel('db/Hz');
            xlabel('Frequency (Hz)');
    end
end