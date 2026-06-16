%% EXAMPLE 1
load('/Users/USER/Documents/aDBS_exampleFigure3A.mat');

legend_elements = [];

figure;
ax(1) = subplot(2,1,2);
hold on;
gait_events_sorted = sortGaitEvents(aligned_data.gait_events,'LTO');
for i = 1:height(gait_events_sorted)
    if (gait_events_sorted.LTO(i)>0 && gait_events_sorted.LHS(i)>0) && (gait_events_sorted.LTO(i)<=aligned_data.right_adaptive_taxis(end) && gait_events_sorted.LHS(i)<=aligned_data.right_adaptive_taxis(end))
        if isempty(legend_elements)
            legend_elements(1) = fill([gait_events_sorted.LTO(i),gait_events_sorted.LHS(i),gait_events_sorted.LHS(i),gait_events_sorted.LTO(i)],[0,0,12000,12000],[98,179,196]./255,'EdgeColor','none');
        else
            fill([gait_events_sorted.LTO(i),gait_events_sorted.LHS(i),gait_events_sorted.LHS(i),gait_events_sorted.LTO(i)],[0,0,12000,12000],[98,179,196]./255,'EdgeColor','none');
        end
    end
end
hold on;
plot(aligned_data.right_adaptive_taxis,aligned_data.right_adaptive_table.Ld0_output,'-','color','#ED2812','linewidth',1.75);
yline(mode(aligned_data.right_adaptive_table.Ld0_lowThreshold),'--k','linewidth',0.75);

A = cellfun(@(x)contains(x,'1'),aligned_data.right_adaptive_table.CurrentAdaptiveState);
B = A(1:end-1) == A(2:end);
adaptive_time = [];
adaptive_state = [];
insert_inds = find(~B);

raw_data_ind = 1;
for i = 1:length(insert_inds)
    adaptive_time = [adaptive_time;aligned_data.right_adaptive_taxis(raw_data_ind:insert_inds(i));aligned_data.right_adaptive_taxis(insert_inds(i)+1)];
    adaptive_state = [adaptive_state;A(raw_data_ind:insert_inds(i));A(insert_inds(i))];
    raw_data_ind = insert_inds(i)+1;
end
plot(adaptive_time,adaptive_state.*12000,'-k','linewidth',1.25);
ylim([-10,13000]);
title('Right Hemisphere');

%%%%%

ax(2) = subplot(2,1,1);
hold on;
gait_events_sorted = sortGaitEvents(aligned_data.gait_events,'RTO');
for i = 1:height(gait_events_sorted)
    if (gait_events_sorted.RTO(i)>0 && gait_events_sorted.RHS(i)>0) && (gait_events_sorted.RTO(i)<=aligned_data.left_adaptive_taxis(end) && gait_events_sorted.RHS(i)<=aligned_data.left_adaptive_taxis(end))
        if length(legend_elements) == 1
            legend_elements(2) = fill([gait_events_sorted.RTO(i),gait_events_sorted.RHS(i),gait_events_sorted.RHS(i),gait_events_sorted.RTO(i)],[0,0,70000,70000],[114,97,170]./255,'EdgeColor','none','DisplayName','Right Leg Swing');
        else
            fill([gait_events_sorted.RTO(i),gait_events_sorted.RHS(i),gait_events_sorted.RHS(i),gait_events_sorted.RTO(i)],[0,0,70000,70000],[114,97,170]./255,'EdgeColor','none');
        end
        
    end
end
hold on;
legend_elements(3) = plot(aligned_data.left_adaptive_taxis,aligned_data.left_adaptive_table.Ld0_output,'-','color','#ED2812','linewidth',1.75);
legend_elements(4) = yline(mode(aligned_data.left_adaptive_table.Ld0_lowThreshold),'--k','linewidth',0.75);

A = cellfun(@(x)contains(x,'0'),aligned_data.left_adaptive_table.CurrentAdaptiveState);
B = A(1:end-1) == A(2:end);
adaptive_time = [];
adaptive_state = [];
insert_inds = find(~B);

raw_data_ind = 1;
for i = 1:length(insert_inds)
    adaptive_time = [adaptive_time;aligned_data.left_adaptive_taxis(raw_data_ind:insert_inds(i));aligned_data.left_adaptive_taxis(insert_inds(i)+1)];
    adaptive_state = [adaptive_state;A(raw_data_ind:insert_inds(i));A(insert_inds(i))];
    raw_data_ind = insert_inds(i)+1;
end
legend_elements(5) = plot(adaptive_time,adaptive_state.*70000,'-k','linewidth',1.25);
ylim([-1,75000])
title('Left Hemisphere');

leg_hand = legend(legend_elements,'Left Leg Swing','Right Leg Swing','Biomarker Power','Threshold','Adaptive State (0/1)','NumColumns',2,'Location','southoutside');
leg_hand.FontSize = 5;
leg_hand.Position(2) = 0.04;
leg_hand.ItemTokenSize = [15,2];

linkaxes(ax,'x');

xlabel(ax(1),'Time (s)');
ylabel(ax,'A.U.');

% Set subplot sizes
set(ax(1),'Position',[0.12,0.22,0.83,0.325]);
set(ax(2),'Position',[0.12,0.62,0.83,0.325]);