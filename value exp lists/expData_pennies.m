function [ dirs, expData  ] = expData_pennies(data_dir)
%%% expData_pennies %%%
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
subdir{l}='760'; l=l+1;
subdir{l}='761'; l=l+1;
subdir{l}='770'; l=l+1;
subdir{l}='772'; l=l+1;
subdir{l}='862'; l=l+1;
subdir{l}='863'; l=l+1;
subdir{l}='866'; l=l+1;
subdir{l}='867'; l=l+1;
subdir{l}='868'; l=l+1;
subdir{l}='880'; l=l+1;
subdir{l}='882'; l=l+1;
subdir{l}='883'; l=l+1;
% subdir{1} = 'saline';
% subdir{l}='0.1'; l=l+1;
% subdir{l}='1'; l=l+1;
% subdir{l}='862A2'; l=l+1;
% subdir{l}='863A1'; l=l+1;
% subdir{l}='863A2'; l=l+1;
% subdir{l}='866A1'; l=l+1;
% subdir{l}='866A2'; l=l+1;
% subdir{l}='867A1'; l=l+1;
% subdir{l}='867A2'; l=l+1;
% subdir{l}='868A1'; l=l+1;
% subdir{l}='868A2'; l=l+1;
% subdir{l}='880A1'; l=l+1;
% subdir{l}='880A2'; l=l+1;
% subdir{l}='881A1'; l=l+1;
% subdir{l}='881A2'; l=l+1;
% subdir{l}='882A1'; l=l+1;
% subdir{l}='882A2'; l=l+1;
% subdir{l}='883A1'; l=l+1;
% subdir{l}='883A2'; l=l+1;
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