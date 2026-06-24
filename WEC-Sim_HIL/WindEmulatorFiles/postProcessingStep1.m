% post processing script to retrieve data from the target

clearvars; close all; clc;

%% === parameters =========================================================
tgName = 'BlueSpeedgoat';
%tgName = 'EGIBaseline2';

mdlName = 'windEmulatorStep4';
dataDir = fullfile('c:','simulink_raw');
processedDir = fullfile('c:','simulink_processed');
matFileName = 'windEmulatorStep4';
% =========================================================================

%% === test Speedgoat connection ==========================================
tg = slrealtime(tgName);

try 
   tg.connect
catch ME
   fprintf('\n*** Target %s not connected. Stopping program. Check connection.\n',tgName)
   fprintf('\n*** Matlab error \n %s \n\n',ME.getReport)   
   return  
end

if tg.isConnected 
   fprintf('\n*** Target %s is connected at IP address %s. \n\n',tg.TargetSettings.name,tg.TargetSettings.address)
end

ipAddress = tg.TargetSettings.address;
% =========================================================================  

%% === copy data from target ==============================================
newFolder = fullfile(dataDir,mdlName);
if exist(newFolder, 'dir')
    rmdir(newFolder,'s') % clear old data
end
mkdir(newFolder);

system(['pscp -r slrt@', ipAddress, ':/home/slrt/applications/', mdlName, '/* ' ,newFolder])
% =========================================================================  

%% === Import Logged Data into MATLAB and view in Simulation Data Inspector
importLogData(newFolder)
Simulink.sdi.view;
% =========================================================================  

%% === Export to .MAT file from Simulation Data Inspector =================
runIDs = Simulink.sdi.getAllRunIDs; % get run ID
runID = runIDs(end);

rtRun = Simulink.sdi.getRun(runID); % get data for last run
SignalData = rtRun.export;

timestamp = rtRun.DateCreated;
timestamp.Format = 'yyyy_MM_dd_HH_mm_ss';
timeStr = char(timestamp); 
dataFileName = [matFileName,'_',timeStr,'.mat'];

Simulink.sdi.exportRun(runID,'to','file','filename',dataFileName); % export to .mat
% =========================================================================  

%% === copy data to unique path for later processing ======================
if ~exist(processedDir,'dir')
    mkdir(processedDir)
end
movefile(dataFileName,processedDir)

fprintf('Wrote file %s to %s',dataFileName,processedDir);
% =========================================================================  







