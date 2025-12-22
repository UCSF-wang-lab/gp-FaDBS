function aggregateTables = aggregateRCSSimSpecData(patientFolder)
% This script goes through each visit folder for the selected patient,
% look for converted data, and assign a class to each converted data
% point. The data in the converted data folder contains simulated RCS power
% values, which is the output of Tanner Dixon's RCS raw LFP conversion
% function. 

patientFolder = '/Users/USER/Documents/P1';
folderList = dir(patientFolder);
subjectID = patientFolder(find(patientFolder=='/',1,'last')+1:end);
aggregateTables = [];

% Go through all of the patient's folders
for i = 1:length(folderList)
    % Check to make sure the Rover folder is not included
    if ~strcmp(folderList(i).name,'Rover')
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('Parsing folder: %s\n',folderList(i).name);
        
        % Grab all the spectrogram and mat files
        spec_files = dir(fullfile(folderList(i).folder,folderList(i).name,'*bitshift*'));
        mat_files = dir(fullfile(folderList(i).folder,folderList(i).name,'Aligned Data'));
        
        if ~isempty(spec_files)            
            % Loop through mat files and associated spec files and fill in a
            % table
            for j = 1:length(mat_files)
                % Load in mat file
                fprintf('\nReading in mat file\n %s\n\n',fullfile(folderList(1).folder,folderList(i).name,mat_files{j}));
                load(fullfile(folderList(1).folder,folderList(i).name,mat_files{j}));
                
                % Check which spec files are associated with this file
                file_name_elements = strsplit(mat_files{j},'_');
                specPatientState = sprintf('%s_stim_%s_med',aligned_data.stim_condition,aligned_data.med_condition);
                specFileMatch = cellfun(@(x)contains(x,specPatientState,'IgnoreCase',true),{spec_files.name});

                if ~isempty(find(cellfun(@(x)contains(x,'Trial'),file_name_elements)))   
                    trial_ind = find(cellfun(@(x)contains(x,'Trial'),file_name_elements))+1;
                    trial_match = cellfun(@(x)contains(x,sprintf('Trial_%d',str2num(file_name_elements{trial_ind})),'IgnoreCase',true),{spec_files.name});
                    specFileMatch = and(specFileMatch,trial_match);
                end
                
                if ~isempty(find(cellfun(@(x)contains(x,'Turn'),file_name_elements)))
                    turn_ind = find(cellfun(@(x)contains(x,'Turn'),file_name_elements))-1;
                    turn_match = cellfun(@(x)contains(x,sprintf('%s_Turn',file_name_elements{turn_ind}),'IgnoreCase',true),{spec_files.name});
                    specFileMatch = and(specFileMatch,turn_match);
                end
                
                if ~isempty(find(cellfun(@(x)contains(x,'Set'),file_name_elements)))
                    set_ind = find(cellfun(@(x)contains(x,'Set'),file_name_elements))+1;
                    set_match = cellfun(@(x)contains(x,sprintf('Set_%d',str2num(file_name_elements{set_ind})),'IgnoreCase',true),{spec_files.name});
                    specFileMatch = and(specFileMatch,set_match);
                end

                if ~isempty(find(cellfun(@(x)contains(x,'Part'),file_name_elements)))
                    part_ind = find(cellfun(@(x)contains(x,'Part'),file_name_elements))+1;
                    part_match = cellfun(@(x)contains(x,sprintf('Part_%d',str2num(file_name_elements{part_ind})),'IgnoreCase',true),{spec_files.name});
                    specFileMatch = and(specFileMatch,part_match);
                end

                
                specFileNames = {spec_files(specFileMatch).name};
                
                % Used to keep track of spec params. Reduces additional
                % computations
                prev_fs = 0;
                prev_nfft = 0;
                prev_fft_int = 0;
                prev_hemisphere = 0;
                
                % Create table of all spec files that are related to each
                % other
                for k = 1:length(specFileNames)
                    fprintf('Associated rcs sim file: %s\n',specFileNames{k});
                    
                    % Load in spectral data
                    specData = readtable(fullfile(folderList(1).folder,folderList(i).name,specFileNames{k}));
                    
                    % Grab current spec parameters and hemisphere
                    specParam = strsplit(specFileNames{k},'_');
                    fs = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'fs'),specParam))+1)));
                    nfft = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'nfft'),specParam))+1)));
                    fft_int = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'fftint'),specParam))+1)));
                    overlap = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'overlap'),specParam))+1)));
                    bitshift = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'fftbitshift'),specParam))+1)));
                    
                    if contains(specFileNames{k},'LEFT')
                        hemisphere = "left";
                    else
                        hemisphere = "right";
                    end
                    
                    % Check to see if this calculation was done previously
                    % but for a different bitshift or hemisphere. If it
                    % wasn't, calculate it for the new fft settings.
                    sameSpecParams = 0;
                    if ~sameSpecParams
                        state = zeros(height(specData),1);
                        
                        % Update values
                        prev_fs = fs;
                        prev_nfft = nfft;
                        prev_fft_int = fft_int;
                        prev_hemisphere = hemisphere;
                        
                        % New parameters, calculate the correct states base
                        % on gait events.
                        filteredGaitEvents = removeGaitCyclesTurns(aligned_data.Xsens,aligned_data.gait_events,30);
                        sortedGaitEvents = sortGaitEvents(filteredGaitEvents,'RTO');
                        if ~isempty(sortedGaitEvents)    
                            for m = 1:height(specData)
                                % Specific for swing phase. Future will change
                                % to make this generic.
                                if strcmp(hemisphere,'left')
                                    cToeOff = (specData.time(m)>=sortedGaitEvents.RTO) & (specData.time(m)<=sortedGaitEvents.RHS);
                                    iToeOff = (specData.time(m)>=sortedGaitEvents.LTO) & (specData.time(m)<=sortedGaitEvents.LHS);
                                    if sum(cToeOff) == 1
                                        state(m) = 1;
                                    elseif sum(iToeOff) == 1
                                        state(m) = 2;
                                    end
                                else
                                    cToeOff = (specData.time(m)>=sortedGaitEvents.LTO) & (specData.time(m)<=sortedGaitEvents.LHS);
                                    iToeOff = (specData.time(m)>=sortedGaitEvents.RTO) & (specData.time(m)<=sortedGaitEvents.RHS);
                                    if sum(cToeOff) == 1
                                        state(m) = 1;
                                    elseif sum(iToeOff) == 1
                                        state(m) = 2;
                                    end
                                end
                            end
                        end
                    end
                    
                    % Create table with states and values. Always
                    % recalculate as verification on the RCS recording
                    % contacts are the same.                
                    includeColumns = ~cellfun(@isempty,cellfun(@(x)regexp(x,'_[4-9]_|_[0-5][0-9]_'),specData.Properties.VariableNames,'UniformOutput',false));
                    currTable = specData(:,includeColumns);
                    
                    % Other columns
                    subjectIDVec = repmat({subjectID},height(currTable),1);                   
                    visitNameVec = repmat({'dbsOptBilateral'},height(currTable),1);
                    
                    medStateVec = repmat({aligned_data.med_condition},height(currTable),1);
                    stimStateVec = repmat({aligned_data.stim_condition},height(currTable),1);
                    hemisphereVec = repmat({hemisphere},height(currTable),1);
                    fsVec = repmat(fs,height(currTable),1);
                    nfftVec = repmat(nfft,height(currTable),1);
                    overlapVec = repmat(overlap,height(currTable),1);
                    fftIntVec = repmat(fft_int,height(currTable),1);
                    bitshiftVec = repmat(bitshift,height(currTable),1);
                    extraTable = table(subjectIDVec,visitNameVec,medStateVec,stimStateVec,hemisphereVec,fsVec,nfftVec,overlapVec,fftIntVec,bitshiftVec,state,...
                        'VariableNames',{'SubjectID','VisitName','MedState','StimState','Side','FS','NFFT','OverlapPercent','FFTInterval','Bitshift','Class'});
                    outTable = [extraTable,currTable];
                    
                    aggregateTableName = sprintf('%s_AggregateSpecData_fs%d_nfft%d',subjectID,fs,nfft);
                    
                    if isempty(aggregateTables) || sum(cellfun(@(x) strcmp(x,aggregateTableName),fields(aggregateTables))) == 0
                        aggregateTables.(aggregateTableName) = outTable;
                    elseif sum(cellfun(@(x) strcmp(x,aggregateTableName),fields(aggregateTables))) == 1
                        aggregateTables.(aggregateTableName) = [aggregateTables.(aggregateTableName);outTable];
                    end
                end
            end
        end
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    end
end
end