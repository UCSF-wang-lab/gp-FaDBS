% Script to extract each patient's gait phase biomarker from aggregated
% simulated RCS power table. Relies on tables that were created using
% aggregateRCSSimSpecData.

files = {'/Users/USER/Documents/P1/Data_Aggregates/P1_AggregateSpecData_fs500_nfft_256.csv';...
    '/Users/USER/Documents/P2/Data_Aggregates/P2_AggregateSpecData_fs500_nfft_256.csv';...
    '/Users/USER/Documents/P3/Data_Aggregates/P3_AggregateSpecData_fs500_nfft_256.csv';...
    '/Users/USER/Documents/P4/Data_Aggregates/P4_AggregateSpecData_fs500_nfft_256_2024-06-12.csv';...
    '/Users/USER/Documents/P5/Data_Aggregates/P5_AggregateSpecData_fs500_nfft_1024.csv'};

left_inds = [59,60;...
    40,41;...
    15,18;...
    19,21;...
    253,258];

left_fft_size = [256,256,256,256,1024];

left_fft_int = [256,256,256,205,410];

left_bitshift = [5,7,3,3,1];

right_inds = [102,111;...
    105,109;...
    14,18;...
    nan,nan;...
    140,142];

right_fft_size = [256,256,256,nan,1024];

right_fft_int = [256,205,256,nan,410];

right_bitshift = [1,0,3,nan,2];

extracted_table = [];

for i = 1:length(files)
    data = readtable(files{i});

    if ~isnan(left_bitshift(i))
        if i < 4
            filter_inds_left = contains(data.VisitName,'Bilateral') & strcmp(data.Side,'left') & data.Bitshift == left_bitshift(i) & data.NFFT == left_fft_size(i) & data.FFTInterval == left_fft_int(i);
            part_table = data(filter_inds_left,[1:10]);
            part_table.PowerBand = sum(data{filter_inds_left,left_inds(i,:)},2);
            extracted_table = [extracted_table;part_table];
        elseif i == 4
            filter_inds_left = contains(data.VisitName,'Unilateral') & strcmp(data.Side,'left') & data.Bitshift == left_bitshift(i) & data.NFFT == left_fft_size(i) & data.FFTInterval == left_fft_int(i);
            part_table = data(filter_inds_left,[1:10]);
            part_table.PowerBand = sum(data{filter_inds_left,left_inds(i,:)},2);
            extracted_table = [extracted_table;part_table];
        elseif i == 5
            filter_inds_left = strcmp(data.Side,'left') & data.Bitshift == left_bitshift(i) & data.NFFT == left_fft_size(i) & data.FFTInterval == left_fft_int(i);
            part_table = data(filter_inds_left,[1:10]);
            part_table.PowerBand = sum(data{filter_inds_left,left_inds(i,:)},2);
            extracted_table = [extracted_table;part_table];
        end
    end

    if ~isnan(right_bitshift(i))
        if i < 4
            filter_inds_right = contains(data.VisitName,'Bilateral') & strcmp(data.Side,'right') & data.Bitshift == right_bitshift(i) & data.NFFT == right_fft_size(i) & data.FFTInterval == right_fft_int(i);
            part_table = data(filter_inds_right,[1:10]);
            part_table.PowerBand = sum(data{filter_inds_right,right_inds(i,:)},2);
            extracted_table = [extracted_table;part_table];
        elseif i == 4
            filter_inds_right = contains(data.VisitName,'Unilateral') & strcmp(data.Side,'right') & data.Bitshift == right_bitshift(i) & data.NFFT == right_fft_size(i) & data.FFTInterval == right_fft_int(i);
            part_table = data(filter_inds_right,[1:10]);
            part_table.PowerBand = sum(data{filter_inds_right,right_inds(i,:)},2);
            extracted_table = [extracted_table;part_table];
        elseif i == 5
            filter_inds_right = strcmp(data.Side,'right') & data.Bitshift == right_bitshift(i) & data.NFFT == right_fft_size(i) & data.FFTInterval == right_fft_int(i);
            part_table = data(filter_inds_right,[1:10]);
            part_table.PowerBand = sum(data{filter_inds_right,right_inds(i,:)},2);
            extracted_table = [extracted_table;part_table];
        end
    end
end