function MP_saveLatent(dataIndex, model_path)

% load the fitting results, separate them into individual sessions for
% later linear regression

animalFolder = unique(dataIndex.LogFilePath);
for ii = 1:length(animalFolder)
    Ind = strfind(animalFolder{ii},filesep);
    startInd = Ind(end);
    animalList{ii} = animalFolder{ii}(startInd:end);
end

for ii = 1:length(animalList)
    filename1 = fullfile(model_path,[animalList{ii},'_FQRPE.mat']);
    if exist(filename1)
        load(filename1);
        display('Now loading model FQ_RPE');
        stats_new = struct();
        % iterate within animals? reassign latent variables into each session
        iterate_Index = find(strcmp(dataIndex.LogFilePath, animalFolder{ii}));
        
        for jj = 1:length(iterate_Index)
            savefilepath = dataIndex.BehPath{iterate_Index(jj)};
            behFile = fullfile(savefilepath, [dataIndex.LogFileName{iterate_Index(jj)}(1:end-4),'_FQRPElatentV.mat']);
            
            field = fieldnames(stats_sim);
            for kk = 1:length(field)-2
                start = 1+ sum(stats_sim.sLength(1:jj-1));
                en = start + stats_sim.sLength(jj) - 1;
                if length(stats_sim.(field{kk})) > 1 && ~strcmp(field{kk},'sLength')
                    stats_new.(field{kk}) = stats_sim.(field{kk})(start:en);
                end
                
                %stats_new.beta = stats_sim.beta;
            end
            stats_new.alpha = stats_sim.alpha;
            stats_new.beta = stats_sim.beta;
            save(behFile,'stats_new');
        end
    end
    filename2 = fullfile(model_path,[animalList{ii},'_FQRPECK.mat']);
    if exist(filename2)
        load(filename2);
        display('Now loading model FQ_RPE_CK');
        stats_new = struct();
        % iterate within animals? reassign latent variables into each session
        iterate_Index = find(strcmp(dataIndex.LogFilePath, animalFolder{ii}));
        
        for jj = 1:length(iterate_Index)
            savefilepath = dataIndex.BehPath{iterate_Index(jj)};
            behFile = fullfile(savefilepath, [dataIndex.LogFileName{iterate_Index(jj)}(1:end-4),'_FQRPECKlatentV.mat']);
            
            field = fieldnames(stats_sim);
            for kk = 1:length(field)-4
                start = 1+ sum(stats_sim.sLength(1:jj-1));
                en = start + stats_sim.sLength(jj) - 1;
                if length(stats_sim.(field{kk})) > 1 && ~strcmp(field{kk},'sLength')
                    stats_new.(field{kk}) = stats_sim.(field{kk})(start:en);
                end
                
                %stats_new.beta = stats_sim.beta;
            end
            stats_new.alpha = stats_sim.alpha;
            stats_new.beta = stats_sim.beta;
            stats_new.alphac = stats_sim.alphac;
            stats_new.betac = stats_sim.betac;
            
            save(behFile,'stats_new');
        end
    end
        
end
end
    
    