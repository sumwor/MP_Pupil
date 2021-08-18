function [ dirs, expData  ] = expData_bandit(data_dir)
%%% expData_bandit %%%
%PURPOSE: Create data structure for imaging tiff files and behavioral log files
%AUTHOR:  AC Kwan, 180429
%
%INPUT ARGUMENTS
%   data_dir:    The base directory to which the raw data are stored.  
%
%OUTPUT VARIABLES
%   dirs:        The subfolder structure within data_dir to work with
%   expData:     Info regarding each experiment

dirs.data = fullfile(data_dir,'data');
dirs.analysis = fullfile(data_dir,'analysis');
dirs.summary = fullfile(data_dir,'summary');

%%
l = 1;
subdir{l}='1802'; l=l+1;

expData=[];
i=1; 
for k=1:numel(subdir)
    cd(fullfile(dirs.data,subdir{k}));  % go to the data directory
    flist=rdir(fullfile('**','*.log'));      % find list of files with .log extension, including subfolders
    for j=1:numel(flist)
        if flist(j).isdir==0    %if it is not a folder
            expData(i).sub_dir = subdir{k};
            expData(i).animalID = k;
            expData(i).logfile = flist(j).name;
            expData(i).logfilenum = j;
            expData(i).onefolder = true;    %many .log files are located in one folder
            i=i+1;
        end
    end  
end