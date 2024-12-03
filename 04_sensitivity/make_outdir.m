%% CREATE OUTPUT DIRECTORY WITH TIMESTAMP
% This function makes a directory 'runs/yyyy-mm-dd/' in the current working directory (if needed),
% adds a subdirectory 'runs/yyyy-mm-dd/runX-THHMM-addname/', and creates nested subdirectories 
% yyyy-mm-dd = today's date
% X = run number today
% THHMM = time in hours and minutes
% Example: runs/2021-01-02/run3-T0456/ = 3rd run of the day at 4:56am on January 2, 2021

function outdir = make_outdir(addname)

%% CREATE PARENT DIRECTORIES, IF NEEDED
% create 'runs' directory and 'runs/[date]' subdirectory for today

% check that runs folder exists, and make it if it doesn't
if isfolder('runs') == 0
    mkdir('runs')
end

% define current date and time
today = datestr(now, 'yyyy-mm-dd');
time = datestr(now, 'HHMM');

% set parent directory within runs to today's date
parent_dir = fullfile('runs', today);

% make runs subfolder for current date (parent directory), if needed
if isfolder(parent_dir) == 0
    mkdir(parent_dir);
end

%% DETERMINE WHICH RUN NUMBER IS NEXT

% initialize a container for existing run numbers, default to 0
run_nums = [0];

% list directories in the parent directory
dirlist = dir(parent_dir);

% extract run number from each directory and add to vector
for i = 1:length(dirlist)
    filename = dirlist(i).name;
    
    if contains(filename, 'run')
        nums = regexp(filename, '\d+', 'match');
        run_nums = [run_nums, str2num(nums{1})];
    end
    
end

% find the maximum run number and increment 1 -> new run number
suffix = max(run_nums) + 1;

%% CREATE NEXT OUTPUT DIRECTORY

% name subfolder in today's results for this run number with timestamp
% outdir format: runs/yyyy-mm-dd/run1-THHMM-addname
base_outdir = fullfile('runs', today, 'run');
outdir = strcat(base_outdir, int2str(suffix), '-T', time, '-', addname);    

% make output directory
mkdir(outdir);      % main output directory
mkdir(fullfile(outdir, 'inputs'))   % subdir for inputs
mkdir(fullfile(outdir, 'results'))  % subdir for results
mkdir(fullfile(outdir, 'figures'))  % subdir for figures

end