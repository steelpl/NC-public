function summary = Run_All_Synthetic(outputRoot, options)
%RUN_ALL_SYNTHETIC Sequential synthetic experiments with isolated run outputs.
%   Run_All_Synthetic; % dependency preflight only; no calculations
%   summary = Run_All_Synthetic([],struct('mode','run'));
%   summary = Run_All_Synthetic([],struct('nMC',3,'nSensitivity',2, ...
%       'nLegacy',2,'runLegacy',true,'mode','run')); % SOFTWARE SMOKE CHECK ONLY
% Defaults: 500 main, 100 sensitivity, 500 legacy replicates.
% Existing nonempty output directories are NEVER overwritten. No auto-resume.
% No global datasets are read by this entry point.
if nargin<1, outputRoot=[]; end
if nargin<2 || isempty(options), options=struct(); end
opts=resolveOptions(options);
codeRoot=fileparts(mfilename('fullpath'));
if ~strcmp(opts.mode,'run')
    summary=existingWorkflow(codeRoot,outputRoot,opts); return
end
previousDir=pwd; previousPath=path; previousRNG=rng;
previousVisible=get(groot,'defaultFigureVisible');
environmentGuard=onCleanup(@() restoreEnvironment(previousDir,previousPath,previousRNG,previousVisible)); %#ok<NASGU>
if isempty(outputRoot)
    outputRoot=fullfile(codeRoot,'results',['synthetic_' datestr(now,'yyyymmdd_HHMMSS_FFF')]);
end
outputRoot=absolutePath(outputRoot);
prepareEmptyOutput(outputRoot);
cd(codeRoot);
addpath(codeRoot,'-begin');
if opts.invisibleFigures, set(groot,'defaultFigureVisible','off'); end
oldDiary=get(0,'Diary'); oldDiaryFile=get(0,'DiaryFile');
diaryGuard=onCleanup(@() restoreDiary(oldDiary,oldDiaryFile)); %#ok<NASGU>
diary(fullfile(outputRoot,'execution.log')); diary on;
fprintf('Synthetic integrated run: %s\n',outputRoot);
fprintf('Main=%d, sensitivity=%d, legacy=%d, runLegacy=%d\n', ...
    opts.nMC,opts.nSensitivity,opts.nLegacy,opts.runLegacy);

names={'01_core_checks','02_legacy_TC_NC','03_legacy_TCH_NC', ...
    '04_legacy_sensitivity','05_legacy_WA','06_revision_synthetic', ...
    '07_identifiability','08_output_validation','09_publication_artifacts'};
stages=repmat(struct('Stage','','Status','pending','StartedUTC','', ...
    'FinishedUTC','','Seconds',NaN,'OutputDirectory','','Message',''),numel(names),1);
for k=1:numel(names), stages(k).Stage=names{k}; end
manifest=struct('runType','synthetic','status','running','completed',false, ...
    'scientificResultsReady',false,'createdUTC',utcNow(),'finishedUTC','', ...
    'codeRoot',codeRoot,'outputRoot',outputRoot,'options',opts, ...
    'matlabVersion',version,'computer',computer,'toolboxes',ver, ...
    'productionReplicationCounts',opts.nMC>=500 && opts.nSensitivity>=100 && (~opts.runLegacy || opts.nLegacy>=500), ...
    'actualGlobalDataExecuted',false,'errorIdentifier','','errorMessage','');
current=0; stageTimer=tic;
try
    required={'NC','WA','TC','TCH'};
    for k=1:numel(required)
        expected=fullfile(codeRoot,[required{k} '.m']);
        assert(strcmpi(which(required{k}),expected),'RunAllSynthetic:ShadowedCore', ...
            'Unexpected %s selected: %s',required{k},which(required{k}));
    end
    assert(exist('Test_Revision_Core','file')==2 && exist('run_revision_synthetic','file')==2 ...
        && exist('run_revision_identifiability','file')==2,'RunAllSynthetic:MissingDependency', ...
        'Install the revision directory alongside this entry point.');
    if opts.makePublicationOutputs
        assert(exist('Build_Synthetic_Revision_Artifacts','file')==2, ...
            'RunAllSynthetic:MissingPostprocessor','Missing integrated/synthetic postprocessor.');
    end
    if opts.runLegacy
        assert(strcmpi(which('EeeTGEN'),fullfile(codeRoot,'EeeTGEN.m')), ...
            'RunAllSynthetic:ShadowedCore','Unexpected EeeTGEN selected: %s',which('EeeTGEN'));
        assert(strcmpi(which('WA_Submitted'),fullfile(codeRoot,'WA_Submitted.m')), ...
            'RunAllSynthetic:MissingBaseline','Missing or shadowed baseline WA estimator.');
        assert(exist('boxplot','file')==2,'RunAllSynthetic:LegacyToolbox', ...
            'Legacy WA plotting needs boxplot (Statistics Toolbox). Set runLegacy=false to run the base-MATLAB revision suite.');
    end
    manifest.sourceFiles=snapshotSources(codeRoot,outputRoot,opts);
    writeState(outputRoot,manifest,stages);
    mainDir=fullfile(outputRoot,names{6});
    sensitivityDir=fullfile(outputRoot,names{7});
    legacyScripts={'Sim_TC_vs_NC.m','Sim_TCH_vs_NC.m','Test_NC_Sensitivity.m','Sim_WA_Comparison.m'};
    for current=1:numel(names)
        folder=fullfile(outputRoot,names{current});
        stages(current).OutputDirectory=folder;
        if (ismember(current,2:5) && ~opts.runLegacy) || (current==9 && ~opts.makePublicationOutputs)
            stages(current).Status='skipped'; stages(current).Message='Disabled by explicit run options.';
            writeState(outputRoot,manifest,stages); continue
        end
        stages(current).Status='running'; stages(current).StartedUTC=utcNow();
        stageTimer=tic; writeState(outputRoot,manifest,stages);
        fprintf('\n[%d/%d] %s\n',current,numel(names),names{current});
        switch current
            case 1
                Test_Revision_Core(folder);
            case {2,3,4,5}
                runLegacy(codeRoot,legacyScripts{current-1},folder,opts.nLegacy);
            case 6
                run_revision_synthetic(opts.nMC,folder);
            case 7
                run_revision_identifiability(opts.nSensitivity,folder);
            case 8
                mkdir(folder);
                audit=validateOutputs(mainDir,sensitivityDir,opts);
                save(fullfile(folder,'validation.mat'),'audit');
                writeJSON(fullfile(folder,'validation.json'),audit);
                manifest.validation=audit;
            case 9
                artifacts=Build_Synthetic_Revision_Artifacts(mainDir,sensitivityDir,folder);
                manifest.publicationArtifacts=artifacts;
        end
        stages(current).Seconds=toc(stageTimer);
        stages(current).Status='complete'; stages(current).FinishedUTC=utcNow();
        writeState(outputRoot,manifest,stages);
    end
    for k=1:numel(manifest.sourceFiles)
        assert(strcmp(fileHash(manifest.sourceFiles(k).sourcePath),manifest.sourceFiles(k).sha256), ...
            'RunAllSynthetic:SourceChanged','Source code changed during execution: %s',manifest.sourceFiles(k).sourcePath);
    end
    manifest.completed=true; manifest.status='complete'; manifest.finishedUTC=utcNow();
    manifest.scientificResultsReady=manifest.productionReplicationCounts && opts.makePublicationOutputs;
    manifest.readinessNote=['Counts and software checks alone do not establish superiority, physical identifiability, ', ...
        'or final manuscript readiness. Review the saved figures, tables and methodological limitations.'];
    writeState(outputRoot,manifest,stages);
    writeRunIndex(outputRoot,manifest,stages);
    fprintf('\nCOMPLETE: %s\nScientific-output readiness flag: %d\n',outputRoot,manifest.scientificResultsReady);
    if ~manifest.productionReplicationCounts
        fprintf('Reduced-replicate SOFTWARE CHECK ONLY: do not use these estimates as the submission experiment.\n');
    end
catch ME
    manifest.completed=false; manifest.scientificResultsReady=false;
    manifest.status='failed'; manifest.errorIdentifier=ME.identifier;
    manifest.errorMessage=ME.message; manifest.finishedUTC=utcNow();
    if current>0
        stages(current).Status='failed'; stages(current).Message=ME.message;
        stages(current).Seconds=toc(stageTimer); stages(current).FinishedUTC=utcNow();
    end
    writeState(outputRoot,manifest,stages);
    fid=fopen(fullfile(outputRoot,'ERROR.txt'),'w','n','UTF-8');
    if fid>=0, fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off')); fclose(fid); end
    fprintf(2,'FAILED: %s. Partial results retained in %s\n',ME.message,outputRoot);
    rethrow(ME)
end
summary=manifest;
end

function opts=resolveOptions(input)
assert(isstruct(input) && isscalar(input),'options must be a scalar struct.');
opts=struct('nMC',500,'nSensitivity',100,'nLegacy',500,'runLegacy',true, ...
    'makePublicationOutputs',true,'invisibleFigures',true,'mode','preflight');
for field=fieldnames(input)'
    assert(isfield(opts,field{1}),'Unknown option: %s',field{1}); opts.(field{1})=input.(field{1});
end
assert(ismember(opts.mode,{'run','preflight','audit','postprocess'}),'Unknown workflow mode.');
for field={'nMC','nSensitivity','nLegacy'}
    x=opts.(field{1}); assert(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=2&&x==fix(x), ...
        '%s must be an integer >=2.',field{1});
end
for field={'runLegacy','makePublicationOutputs','invisibleFigures'}
    assert(islogical(opts.(field{1}))&&isscalar(opts.(field{1})),'%s must be logical.',field{1});
end
end

function summary=existingWorkflow(root,runDir,options)
oldDir=pwd; oldPath=path; oldRng=rng;
guard=onCleanup(@()restoreEnvironment(oldDir,oldPath,oldRng,get(groot,'defaultFigureVisible'))); %#ok<NASGU>
cd(root); addpath(root,'-begin');
if strcmp(options.mode,'preflight')
    for name={'NC','WA','WA_Submitted','TC','TCH','run_revision_synthetic','run_revision_identifiability','Build_Synthetic_Revision_Artifacts','Test_Revision_Core'}
        assert(strcmpi(which(name{1}),fullfile(root,[name{1} '.m'])),'Missing or shadowed %s.',name{1});
    end
    summary=struct('status','PREFLIGHT_ONLY','executionCompleted',false,'codeRoot',root);
    fprintf('PASS: synthetic dependencies resolve to %s. No experiments executed.\\n',root); return
end
assert(~isempty(runDir),'Provide the existing completed run directory.');
runDir=absolutePath(runDir);
manifest=jsondecode(fileread(fullfile(runDir,'run_manifest.json')));
mainDir=fullfile(runDir,'06_revision_synthetic'); sensitivityDir=fullfile(runDir,'07_identifiability');
summary=validateOutputs(mainDir,sensitivityDir,manifest.options);
summary.originalRunCompleted=manifest.completed;
summary.mode=options.mode; summary.sourceRun=runDir;
if strcmp(options.mode,'postprocess')
    output=fullfile(runDir,['postprocess_' datestr(now,'yyyymmdd_HHMMSS_FFF')]);
    prepareEmptyOutput(output);
    summary.publication=Build_Synthetic_Revision_Artifacts(mainDir,sensitivityDir,output);
    writeJSON(fullfile(output,'postprocess_manifest.json'),summary);
end
end

function prepareEmptyOutput(folder)
if isfolder(folder)
    entries=dir(folder); entries=entries(~ismember({entries.name},{'.','..'}));
    assert(isempty(entries),'RunAllSynthetic:OutputExists','Output folder must be new or empty: %s',folder);
else
    assert(~isfile(folder),'Output path is a file: %s',folder);
    [ok,msg]=mkdir(folder); assert(ok,'Cannot create output folder: %s',msg);
end
end

function rows=snapshotSources(root,out,opts)
files={'Run_All_Synthetic.m','NC.m','WA.m','TC.m','TCH.m', ...
    'Test_Revision_Core.m','run_revision_synthetic.m', ...
    'run_revision_identifiability.m'};
if opts.makePublicationOutputs
    files=[files {'Build_Synthetic_Revision_Artifacts.m','JoHFigureStyle.m', ...
        'Draw_JoH_Graphical_Abstract.m'}];
end
if opts.runLegacy
    files=[files {'EeeTGEN.m','WA_Submitted.m','Sim_TC_vs_NC.m','Sim_TCH_vs_NC.m','Test_NC_Sensitivity.m','Sim_WA_Comparison.m'}];
end
rows=repmat(struct('relativePath','','sourcePath','','sha256',''),numel(files),1);
for k=1:numel(files)
    source=fullfile(root,files{k}); target=fullfile(out,'source_snapshot',files{k});
    assert(isfile(source),'Missing source file: %s',source);
    if ~isfolder(fileparts(target)), mkdir(fileparts(target)); end
    copyfile(source,target);
    rows(k)=struct('relativePath',files{k},'sourcePath',source,'sha256',fileHash(source));
end
end

function runLegacy(root,name,folder,nMC)
mkdir(folder);
source=fileread(fullfile(root,name));
pattern='(?m)^nMC\s*=\s*\d+\s*;';
assert(numel(regexp(source,pattern,'match'))==1,'Expected one top-level nMC assignment: %s',name);
effective=regexprep(source,pattern,sprintf('nMC = %d;',nMC));
script=fullfile(folder,name); fid=fopen(script,'w','n','UTF-8'); assert(fid>=0);
fprintf(fid,'%s',effective); fclose(fid);
fid=fopen(fullfile(folder,'INTERPRETATION.md'),'w','n','UTF-8'); assert(fid>=0);
fprintf(fid,['# Original experiment design with current NC\n\nSource: %s; nMC=%d. The copied script changes only the nMC assignment. ', ...
    'Current NC.m (default PSD feasibility) is used; WA_Submitted retains the pre-revision WA rule. This rerun does not reproduce the pre-revision NC estimator. Original gain assumptions and population/in-sample ', ...
    'evaluation remain unchanged. This is not the independent-test, unit-gain revision experiment. ', ...
    'Legacy figures retain their original failure definition; use the revision outputs for PSD/WA diagnostics.\n'],name,nMC);
fclose(fid);
previous=pwd; guard=onCleanup(@() cd(previous)); %#ok<NASGU>
oldFigures=findall(groot,'Type','figure');
figureGuard=onCleanup(@() closeCreatedFigures(oldFigures)); %#ok<NASGU>
cd(folder); executeLegacyScript(script);
assert(isfile(fullfile(folder,'workspace.mat')),'Legacy workspace was not saved.');
saved=load(fullfile(folder,'workspace.mat'),'nMC');
assert(isfield(saved,'nMC') && saved.nMC==nMC,'Legacy replicate override was not applied.');
expectedNames={'Sim_TC_vs_NC.m','Sim_TCH_vs_NC.m','Test_NC_Sensitivity.m','Sim_WA_Comparison.m'};
expectedPNG={'Sim_TC_vs_NC_results.png','Sim_TCH_vs_NC_results.png', ...
    'Test_NC_Sensitivity_results.png','Sim_WA_Comparison.png'};
whichScript=find(strcmp(name,expectedNames),1);
imageFile=fullfile(folder,'figures',expectedPNG{whichScript});
assert(isfile(imageFile),'Missing native figure: %s',imageFile);
imageInfo=imfinfo(imageFile); assert(imageInfo.Width>0 && imageInfo.Height>0,'Invalid native PNG.');
allFigures=findall(groot,'Type','figure'); newFigures=allFigures(~ismember(allFigures,oldFigures));
assert(~isempty(newFigures),'Legacy script produced no figure: %s',name);
for k=1:numel(newFigures)
    base=fullfile(folder,sprintf('legacy_figure_%02d',k));
    exportgraphics(newFigures(k),[base '.pdf'],'ContentType','vector');
    savefig(newFigures(k),[base '.fig']);
end
close(newFigures);
end

function executeLegacyScript(scriptPath)
% Script CLEAR acts only in this isolated workspace, not in the orchestrator.
try
    run(scriptPath);
    save('workspace.mat','-v7.3');
catch legacyError
    try
        save('workspace_failed.mat','-v7.3');
    catch
        % Preserve the original scientific/plotting exception if saving fails.
    end
    rethrow(legacyError)
end
end

function closeCreatedFigures(oldFigures)
figures=findall(groot,'Type','figure');
figures=figures(~ismember(figures,oldFigures));
if ~isempty(figures), close(figures); end
end

function audit=validateOutputs(mainDir,sensitivityDir,opts)
R=readtable(fullfile(mainDir,'replicate_metrics.csv'),'TextType','string');
S=readtable(fullfile(mainDir,'scenario_summary.csv'),'TextType','string');
P=readtable(fullfile(mainDir,'paired_seasonal_minus_stationary.csv'),'TextType','string');
I=readtable(fullfile(sensitivityDir,'multistart_floor_replicates.csv'),'TextType','string');
assert(height(R)==36*6*opts.nMC && height(S)==36*6,'Unexpected main output row counts.');
keys=string(R.Scenario)+"_"+string(R.Method); groups=unique(keys);
assert(numel(groups)==36*6,'Missing scenario/method combinations.');
for k=1:numel(groups), assert(sum(keys==groups(k))==opts.nMC,'Unequal replicate support.'); end
assert(~hasError(R.EstimatorError) && ~hasError(R.MergeError),'Unexpected estimator/merge exceptions in output.');
finiteMerge=isfinite(R.MergeMSE); weights=R.Weight1+R.Weight2+R.Weight3;
assert(all(abs(weights(finiteMerge)-1)<1e-8),'Weights do not sum to one.');
assert(all(R.MergeMSE(finiteMerge)>=0),'Negative merge MSE.');
psd=string(R.Method)=="NC_PSD";
assert(all(isfinite(R.MinEigenvalue(psd))) && all(R.MinEigenvalue(psd)>=-1e-8),'Invalid NC PSD output.');
assert(all(R.NonfiniteEstimate(psd)==0 & R.NegativeVariance(psd)==0),'Invalid NC variance output.');
assert(all(abs(R.SampleTruthVariance-1)<1e-10),'Truth variance normalization differs.');
paired=ismember(string(P.Method),["equal_weight","oracle"]) & P.SignalErrorCorrelation==0;
assert(all(abs(P.MergeMSEDifference(paired))<1e-12),'Paired synthetic control mismatch.');
assert(height(I)==opts.nSensitivity*6*5,'Unexpected sensitivity row count.');
assert(all(isfinite(I.HeldoutMergeMSE)) && all(I.HeldoutMergeMSE>=0),'Invalid sensitivity merging output.');
for k=1:height(S)
    use=R.Scenario==S.Scenario(k) & string(R.Method)==string(S.Method(k));
    assert(S.NMerge(k)==sum(isfinite(R.MergeMSE(use))) && ...
        S.NVarianceEstimates(k)==sum(isfinite(R.VarianceMSE(use))),'Summary denominators differ.');
end
C=readtable(fullfile(sensitivityDir,'p2_exact_nonuniqueness.csv'));
assert(all(C.MinEigenvalue>0) && all(abs(C.OffDiagonalLoss)<1e-12),'Invalid p2 counterexample.');
audit=struct('status','PASS','mainRows',height(R),'scenarios',36, ...
    'replicatesPerScenario',opts.nMC,'finiteMerges',sum(finiteMerge), ...
    'NCPSDEstimates',sum(psd),'sensitivityRows',height(I),'pairedControlChecks',sum(paired));
fprintf('PASS: output consistency audit (%d rows, %d scenarios).\n',height(R),36);
end

function yes=hasError(value)
if isnumeric(value), yes=any(isfinite(value));
else, text=string(value); yes=any(~ismissing(text) & strlength(strtrim(text))>0); end
end

function writeState(folder,manifest,rows)
stages=struct2table(rows);
save(fullfile(folder,'run_manifest.mat'),'manifest','stages');
writeJSON(fullfile(folder,'run_manifest.json'),manifest);
writetable(stages,fullfile(folder,'stage_status.csv'));
end

function writeRunIndex(folder,manifest,stages)
fid=fopen(fullfile(folder,'RESULTS.md'),'w','n','UTF-8'); assert(fid>=0);
guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Synthetic run outputs\n\nCompleted: %s\n\n',manifest.finishedUTC);
fprintf(fid,'Main=%d; sensitivity=%d; legacy=%d. Full-count flag=%d.\n\n', ...
    manifest.options.nMC,manifest.options.nSensitivity,manifest.options.nLegacy,manifest.productionReplicationCounts);
fprintf(fid,'| Stage | Status | Seconds | Folder |\n|---|---|---:|---|\n');
for k=1:numel(stages)
    fprintf(fid,'| %s | %s | %.2f | %s |\n',stages(k).Stage,stages(k).Status,stages(k).Seconds,stages(k).Stage);
end
fprintf(fid,['\nRead stage 09 for CSV tables, figures and the synthetic results report. ', ...
    'Stages 02-05 retain legacy designs; stages 06-07 provide the added revision experiments. ', ...
    'Reduced-replicate runs are software checks only. No global data were executed.\n']);
end

function writeJSON(file,value)
fid=fopen(file,'w','n','UTF-8'); assert(fid>=0,'Cannot write %s',file);
guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end

function hash=fileHash(file)
fid=fopen(file,'rb'); assert(fid>=0); bytes=fread(fid,Inf,'*uint8'); fclose(fid);
digest=java.security.MessageDigest.getInstance('SHA-256'); digest.update(typecast(bytes,'int8'));
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end

function folder=absolutePath(value)
assert((ischar(value)&&isrow(value)) || (isstring(value)&&isscalar(value)), 'Output path must be text.');
folder=char(java.io.File(char(value)).getCanonicalPath());
end

function text=utcNow
text=char(datetime('now','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss'Z'"));
end

function restoreDiary(state,file)
diary off;
set(0,'DiaryFile',file);
if strcmpi(state,'on'), diary on; end
end

function restoreEnvironment(folder,searchPath,rngState,visible)
path(searchPath); cd(folder); rng(rngState); set(groot,'defaultFigureVisible',visible);
end
