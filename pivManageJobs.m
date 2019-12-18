function [im1out,im2out,pivParOut] = pivManageJobs(im1,im2,pivParIn)
% pivDistributeTreatment - divide the treatment of image sequence into several jobs and outputs settings for the
% first untreated job.
%
% Usage:
% [im1out,im2out,pivParOut] = pivManageJobs(im1,im2,pivParIn)
%
% Inputs:
%    im1,im2 ... list of images to be treated (cell containing paths to files with first and second image in image
%                pair)
%    pivParIn ... (struct) parameters defining the evaluation. Following fields are considered:
%       jmParallelJobs ... defines, to how many jobs is the treatment of present task distributed. Used by
%          pivManageJobs.m.
%       jmLockExpirationTime ... maximum age of a lock file (in seconds). If a lock file is older than this limit, 
%          it is regarded as non-existent. 
%       anTargetPath ... specified folder, to which files containing results of PIV analysis of each image
%          pair are stored. 
%       anPairsOnly ... if true and if .anOnDrive is true, the subroutine pivAnalyzeImageSequence.m processes only 
%          image pairs, but it does not produce final file with sequence data. If this option is set, much less
%          memory is required for processing.
% Outputs:
%    im1out,im2out ... list of images to be treated in the job which should start now (cell containing paths to 
%          files with first and second image in image pair)
%    pivParOut ... (struct) parameters defining the evaluation by the next jobs. Most fields are copied from 
%          pivParIn with these exceptions:
%       jmLockFile ... (this field is added). Contains the name of a file, which pivAnalyzeImageSequence should 
%          update every time when it treats an image par. 
%       anPairsOnly ... This parameer is set to true, if job to be started is not the last job.
%
%
% How this subroutine should be used:
%
%    This subroutine allows to distribute treatment of an image sequence to several independently running Matlab
%    instances (they can run on the same computer, or on several different computers is they are sharing the same
%    data folder via network). By this, faster treatment of an image sequence can be achieved.
%      The usage is following (see also example_06b_Sequence_multiprocessor.m):
%        1. Prepare and debug the treatment as for the use on a single machine
%        2. Modify the program in the following way:
%             - set the parameter pivPar.jmParallelJobs to the number of matlab instances N, which will be treating
%               the sequence
%             - let this subroutine to modify image lists im1 and im2, and parameters pivPar, by calling it just 
%               before command pivAnalyzeImageSequence
%        3. Start N Matlab instances and run the treatment in each of them.
%        4. Not all Matlab instances will "remember" the results of the treatment (most of them will run with
%           parameter anPairsOnly set). To get result of all the sequence, run the treatment program once more.
%        5. If something goes wrong, go to the output folder (specified by anTargetPath) and erase files
%           "Joblist.mat" and all files with filename in form "lock*.lck").
%
%
% How this subroutine works:
%     1. It assures that only one instance of this subroutine is running at the same time. Therefore, it creates a
%        lock file (file "lock_Master_Wait.lck" in the output folder). If this file exists in the output folder, if
%        waits until this file is deleted, or until it gets expired (older then jmLockExpirationTime). When this
%        subroutine is finishing, it erases this lock file. (Parts 2 - 3 of the code below)
%     2. Subroutine checks, if the treatment of image sequence is already distributed in several jobs, or if this
%        distribution should be done (typically when this subroutine is run for the first time). The distribution
%        should be done if
%          - file JobList.mat is not present in the output folder
%          - if some of the lock files in the output folder is expired
%        If the treatment is already distributed into several jobs, read information about jobs (from file
%        JobList.mat), save it variable JobList, and skip the next step.
%        (Parts 5 - 6 of the code below)
%     3. Distribute the treatment of image sequence inseveral jobs. For this purpose, 
%          - compare the output files present in the output folder with expected list of files, which should be
%            present after all the treatment
%          - check, which jobs are running (non-expired lock file "lock_Job_XXX.lck" is present in output folder,
%            where XXX is the job number) and identify, how many jobs still should be started
%          - distribute evenly the image pairs to jobs; save information about images treated by each job to
%            variable JobList; save there also information that this particular job should start
%          - Note: to Job no. 1, all image pairs are always attributed. This job will hence treat all image pairs
%            (reading results of previous treatment by another jobs from output files, if they are available)
%        (Parts 7a, 7b and 7c of the code below)
%     4. To outputs of subroutine, assign (if there is some job to be started)
%          - image lists of image pairs, which should be treated by the next job (job with highest number, marked 
%            as "should start" in variable JobList)
%          - pivParIn; if job number is bigger than 1, set anPairsOnly to true;
%          - before leaving, create lock file "lock_Job_XXX.lck", where XXX is job number to be started as next
%        If there is no job to be started and file "lock_Job_001.lck" does not exist, assign to output images for 
%        job no. 1 (containing all image pairs). The subsequent treatment hence would read all available results.
%        If the "lock_Job_001.lck" exists and is not expired, then output empty lists (treatment will finish).
%          
%
%        
%
%%
% This subroutine is a part of
%
% =========================================
%               PIVsuite
% =========================================
%
% PIVsuite is a set of subroutines intended for processing of data acquired with PIV (particle image
% velocimetry) within Matlab environment.
%
% Written by Jiri Vejrazka, Institute of Chemical Process Fundamentals, Prague, Czech Republic
%
% For the use, see files example_XX_xxxxxx.m, which acompany this file. PIVsuite was tested with
% Matlab 8.2 (R2013b).
%
% In the case of a bug, please, contact me: vejrazka (at) icpf (dot) cas (dot) cz
%
%
% Requirements:
%     Image Processing Toolbox
%         (required only if pivPar.smMethod is set to 'gaussian')
%
%     inpaint_nans.m
%         subroutine by John D'Errico, available at http://www.mathworks.com/matlabcentral/fileexchange/4551
%
%     smoothn.m
%         subroutine by Damien Garcia, available at
%         http://www.mathworks.com/matlabcentral/fileexchange/274-smooth
%
% Credits:
%    PIVsuite is a redesigned version of a part of PIVlab software [3], developped by W. Thielicke and
%    E. J. Stamhuis. Some parts of this code are copied or adapted from it (especially from its
%    piv_FFTmulti.m subroutine).
%
%    PIVsuite uses 3rd party software:
%        inpaint_nans.m, by J. D'Errico, [2]
%        smoothn.m, by Damien Garcia, [5]
%
% References:
%   [1] Adrian & Whesterweel, Particle Image Velocimetry, Cambridge University Press, 2011
%   [2] John D'Errico, inpaint_nans subroutine, http://www.mathworks.com/matlabcentral/fileexchange/4551
%   [3] W. Thielicke and E. J. Stamhuid, PIVlab 1.31, http://pivlab.blogspot.com
%   [4] Raffel, Willert, Wereley & Kompenhans, Particle Image Velocimetry: A Practical Guide. 2nd edition,
%       Springer, 2007
%   [5] Damien Garcia, smoothn subroutine, http://www.mathworks.com/matlabcentral/fileexchange/274-smooth
%


pivParIn = pivParams([],pivParIn,'DefaultsJobManagement');

%% 1 Check the existance of output folder. If it does not exist, try to create it. If unsuccessful, give a
% message.

if ~pivParIn.anOnDrive
    error('pivManageJobs.m: Job management can be used only with option anOnDrive = true.');
end

if ~exist(pivParIn.anTargetPath,'dir')
    try
        mkdir(pivParIn.anTargetPath);
    catch  %#ok<*CTCH>
        error('pivDistributeTreatment: Target folder does not exist. Failed to create it.');
    end
end


%% 2 Check the presence of lock_Master_Wait file. If exists or expired, wait.
pause(0.5+0.5*rand);
auxt0 = tic;
auxLastEcho = -inf;
auxLastPing = -inf;
auxEcho = false;
while fileAge([pivParIn.anTargetPath, '/lock_Master_Wait.lck'])<pivParIn.jmLockExpirationTime
    auxEcho = true;
    if toc(auxt0)-auxLastEcho > 60
        fprintf('\nJob distribution file is locked. Waiting for unlocking...');
        auxLastEcho = toc(auxt0);
    end
    if toc(auxt0)-auxLastPing > 10
        fprintf(' (%4.1f s)...',toc(auxt0));
        auxLastPing = toc(auxt0);
    end
    pause(0.5+0.5*rand);
end
if auxEcho, fprintf('\n'); end


%% 3 Create the lock file with _Wait status.
if exist([pivParIn.anTargetPath, '/lock_Master_Wait.lck'],'file')
    delete([pivParIn.anTargetPath, '/lock_Master_Wait.lck']);
end
fMasterLock = fopen([pivParIn.anTargetPath, '/lock_Master_Wait.lck'],'w');
fprintf(fMasterLock,'Distributing/identifying tasks.');
fclose(fMasterLock);


%% 5 Read JobList, if the file exists, or initialize this variable:
RedistributeJobs = false;
if exist([pivParIn.anTargetPath, '/JobList.mat'],'file')
    load([pivParIn.anTargetPath, '/JobList.mat']);
    delete([pivParIn.anTargetPath, '/JobList.mat']);
    try
        if numel(JobList.ShouldStart)~=pivParIn.jmParallelJobs     %#ok<NODEF>
            clear('JobList');
        end
    catch
    end
end
if ~exist('JobList','var')
    RedistributeJobs = true;
    JobList.lockFiles = cell(pivParIn.jmParallelJobs,1);
    JobList.ShouldStart = true(pivParIn.jmParallelJobs,1);
    JobList.im1 = cell(pivParIn.jmParallelJobs,1);
    JobList.im2 = cell(pivParIn.jmParallelJobs,1);
end


%% 6 Test, if jobs should be redistributed:
auxI = find(~JobList.ShouldStart);
if ~RedistributeJobs
    for kk = 1:numel(auxI)
        if (~exist(JobList.lockFiles{auxI(kk)},'file') || ...
                fileAge(JobList.lockFiles{auxI(kk)})>pivParIn.jmLockExpirationTime) && ...
                ~exist([JobList.lockFiles{auxI(kk)}(1:end-4) '_Finished.lck'],'file')
            RedistributeJobs = true;
        end
    end
end

%% 7 Distribute jobs (if neccessary): a) find, which image pairs should be treated; b) distribute jobs to jobs with
%    non-existent or expired lock file

if RedistributeJobs
    %% 7a Get list of image pairs, for which result files are missing
    % get list of existing output files
    fprintf('Distributing image pairs to %d jobs... ',pivParIn.jmParallelJobs);
    tic;
    aux = dir([pivParIn.anTargetPath, '/piv*.mat']);
    % renew lock file (dir command might take some time)
    fMasterLock = fopen([pivParIn.anTargetPath, '/lock_Master_Wait.lck'],'w');
    fprintf(fMasterLock,'Distributing/identifying tasks.');
    fclose(fMasterLock);
    outputList = cell(numel(aux),1);
    for ki = 1:numel(aux)
        outputList{ki,1} = aux(ki).name;
    end
    % add image pairs queued for processing to the list of existing jobs
    for kk = 2:pivParIn.jmParallelJobs;     % start at 2 (Job001 treat always all image pairs)
        % add images only if corresponding lock file is not expired
        if fileAge([pivParIn.anTargetPath, '/lock_Job_' num2str(kk,'%03d') '.lck'])< pivParIn.jmLockExpirationTime
            % lock file for kkth job is not expired
            JobList.ShouldStart(kk) = false;
            auxim1 = JobList.im1{kk};
            auxim2 = JobList.im2{kk};
            auxList = cell(numel(auxim1),1);
            for ki = 1:numel(auxim1)
                [~,filename1] = treatImgPath(auxim1{ki});
                [~,filename2] = treatImgPath(auxim2{ki});
                auxList{ki,1} = ['piv_',filename1, '_', filename2, '.mat'];
            end
            outputList = appendCells(outputList,auxList);
        else
            % lock file for kkth job is expired
            JobList.ShouldStart(kk) = true;
            if exist([pivParIn.anTargetPath, '/lock_Job_' num2str(kk,'%03d') '.lck'],'file')
                delete([pivParIn.anTargetPath, '/lock_Job_' num2str(kk,'%03d') '.lck']);
            end
        end
    end
    if numel(outputList)>0, outputList = sort(outputList); end;
    % get list of required output files
    requiredList = cell(numel(im1),1);
    for ki = 1:numel(im1)
        [~,filename1] = treatImgPath(im1{ki});
        [~,filename2] = treatImgPath(im2{ki});
        requiredList{ki,1} = ['piv_',filename1, '_', filename2, '.mat'];
    end
    [requiredList,auxSortOrder] = sort(requiredList);
    % check presence of required files in the list of existing/queued files (complicated for better speed...)
    MissingI = ones(size(requiredList));
    if numel(outputList)>0
        auxPrevious = 0;
        for ki = 1:numel(requiredList)
            auxFound = false; % first check, if it is not close to previously found file
            for kj = 1:10
                try
                    if strcmpi(requiredList{ki,1},outputList{auxPrevious+kj})
                        auxPrevious = auxPrevious + kj;
                        auxFound = true;
                        MissingI(ki) = 0;
                        break
                    end
                catch
                end
            end
            if auxFound, continue; end
            for kj = 1:numel(outputList)
                if strcmpi(requiredList{ki,1},outputList{kj})
                    auxPrevious = kj;
                    MissingI(ki) = 0;
                    break
                end
            end
            % update lock file sometimes
            if ki/200 == round(ki/200)
                fMasterLock = fopen([pivParIn.anTargetPath, '/lock_Master_Wait.lck'],'w');
                fprintf(fMasterLock,'Distributing/identifying tasks.');
                fclose(fMasterLock);
            end
        end
    end
    % keep indices of image pairs with missing results
    MissingI = auxSortOrder(logical(MissingI));
    MissingPairs = numel(MissingI);
    % keep only missing files in lists of images
    MissingIm1 = cell(numel(MissingI),1);
    MissingIm2 = cell(numel(MissingI),1);
    for kk = 1:numel(MissingI)
        MissingIm1{kk} = im1{MissingI(kk)};
        MissingIm2{kk} = im2{MissingI(kk)};
    end
    % If only few image pairs should be treated, limit number of starting jobs:
    if MissingPairs < pivParIn.jmParallelJobs
        JobList.ShouldStart(MissingPairs+1:end) = false;
    end
    % Now we have a list of image pairs, for which a result file is missing and their treatment is not queued among
    % active jobs.
    % Check whether Job001 (treating all image pairs) is active
    if ~exist([pivParIn.anTargetPath, '/lock_Job_001.lck'], 'file')
        JobList.ShouldStart(1) = 1;
    elseif fileAge([pivParIn.anTargetPath, '/lock_Job_001.lck']) > pivParIn.jmLockExpirationTime
        delete([pivParIn.anTargetPath, '/lock_Job_001.lck']);
        JobList.ShouldStart(1) = 1;
    end
    % Get number of jobs, among which the treatment of remaining image pairs will be distributed
    %% 7b Distribute image pairs to jobs
    I = find(JobList.ShouldStart);
    if pivParIn.anPairsOnly
        auxFirstNormalJob = 1;
    else 
        auxFirstNormalJob = 2;
    end
    for kk = numel(I):-1:auxFirstNormalJob
        auxStartI = ceil(1+(kk-1)*(MissingPairs-1)/numel(I));
        auxStopI = floor(1+kk*(MissingPairs-1)/numel(I));
        auxIm1 = cell(auxStopI-auxStartI,1);
        auxIm2 = cell(auxStopI-auxStartI,1);
        for ki = auxStartI:auxStopI
            auxIm1{ki-auxStartI+1} = MissingIm1{ki};
            auxIm2{ki-auxStartI+1} = MissingIm2{ki};
        end
        JobList.lockFiles{I(kk)} = [pivParIn.anTargetPath, '/lock_Job_' num2str(I(kk),'%03d') '.lck'];
        JobList.im1{I(kk)} = auxIm1;
        JobList.im2{I(kk)} = auxIm2;
        if exist(JobList.lockFiles{I(kk)},'file')
            delete(JobList.lockFiles{I(kk)});
        end
        if exist([JobList.lockFiles{I(kk)}(1:end-4) '_Finished.lck'],'file')
            delete([JobList.lockFiles{I(kk)}(1:end-4) '_Finished.lck']);
        end
        
    end
        
        
    %% 7c If Job001 is to be started and ~anPairsOnly, attribute to it all image pairs
    if JobList.ShouldStart(1) && ~pivParIn.anPairsOnly
        auxIm1 = cell(numel(im1),1);
        auxIm2 = cell(numel(im1),1);
        for ki = 1:numel(im1)
            auxIm1{ki} = im1{ki};
            auxIm2{ki} = im2{ki};
        end
        JobList.lockFiles{1} = [pivParIn.anTargetPath, '/lock_Job_001.lck'];
        JobList.im1{1} = auxIm1;
        JobList.im2{1} = auxIm2;
        if exist(JobList.lockFiles{1},'file')
            delete(JobList.lockFiles{1});
        end
        if exist([JobList.lockFiles{1}(1:end-4) '_Finished.lck'],'file')
            delete([JobList.lockFiles{1}(1:end-4) '_Finished.lck']);
        end
    end
    fprintf('Finished in %5.2f sec.\n',toc);
else  % end of "if RedistributeJobs"
    fprintf('Job distribution read from file. \n');
end   % end of "if RedistributeJobs"


%% 8. Start last unstarted task from the list
I = find(JobList.ShouldStart);
if numel(I)>0
    I = I(end);
    % output required data
    pivParOut = pivParIn;
    pivParOut.seqJobNumber = I;
    pivParOut.jmLockFile = JobList.lockFiles{I};
    if I>1, pivParOut.anPairsOnly = true; end
    im1out = JobList.im1{I};
    im2out = JobList.im2{I};
    % create lock file for the given task
    fLock = fopen(JobList.lockFiles{I},'w');
    fprintf(fLock,'Starting task.');
    fclose(fLock);
    % update master lock file
    JobList.ShouldStart(I) = false;
    save([pivParIn.anTargetPath, '/JobList.mat'],'JobList');
elseif fileAge([pivParIn.anTargetPath, '/lock_Job_001.lck'])> pivParIn.jmLockExpirationTime
    % nothing should be started: run the treatment, only if lockFile{1} does not exist:
    pivParOut = pivParIn;
    pivParOut.seqJobNumber = 1;
    pivParOut.jmLockFile = [pivParIn.anTargetPath, '/lock_Job_001.lck'];
    im1out = im1;
    im2out = im2;
    % create lock file for the given task
    fLock = fopen(pivParOut.jmLockFile,'w');
    fprintf(fLock,'Starting task.');
    fclose(fLock);
else
    pivParOut = [];
    im1out = [];
    im2out = [];
    fprintf('Treatment of all image pairs is attributed to jobs. No treatment will occur. \n');
    fprintf('    (User might want to erase lock files in the output folder.)\n');
end

% remove the lock file with _Wait status.
pause(0.5);
delete([pivParIn.anTargetPath, '/lock_Master_Wait.lck']);

end


%% local functions
function [imgNo, filename, folder] = treatImgPath(path)
% separate the path to get the folder, filename, and number if contained in the name
filename = '';
imgNo = [];
folder = '';
if numel(path)>0
    path = path(end:-1:1);
    I = find((path=='\')|(path=='/'));
    I = I(1);
    Idot = find(path=='.');
    Idot = Idot(1);
    try
        folder = path(I+1:end);
        folder = folder(end:-1:1);
    catch  %#ok<CTCH>
        folder = '';
        I = length(path)+1;
    end
    try
        filename = path(Idot+1:I-1);
        filename = filename(end:-1:1);
    catch  %#ok<CTCH>
        filename = '';
    end
    try
        aux = regexp(filename,'[0-9]');
        aux = filename(aux);
        imgNo = str2double(aux);
    catch  %#ok<CTCH>
        imgNo = [];
    end
end
end


%% local function
function [out] = appendCells(in1,in2)
out = cell(numel(in1)+numel(in2),1);
for ki = 1:numel(in1)
    out{ki,1} = in1{ki};
end
for ki = 1:numel(in2)
    out{ki + numel(in1),1} = in2{ki};
end
end


%% local function
function [age] = fileAge(filename)
try
    aux = dir(filename);
    if numel(aux)==0
        age = Inf;
    else
        age = (now-aux.datenum)*24*3600;
    end
catch
    age = Inf;
end
end