function report = Export_Public_Figures(runDir, outputDir)
%EXPORT_PUBLIC_FIGURES Export a compact gallery from a completed production run.
% Run Figure_Supplementary_NC_PSD_Diagnostic from runDir first.
% No covariance fits are performed here. Existing nonempty galleries are rejected.
% PNG/PDF and small numeric summaries are public; workspaces remain in results/.
root=fileparts(mfilename('fullpath'));
if nargin<2 || isempty(outputDir), outputDir=fullfile(root,'figures'); end
manifest=jsondecode(fileread(fullfile(runDir,'run_manifest.json')));
assert(manifest.completed && manifest.productionReplicationCounts, ...
    'NCGallery:Incomplete','A completed production run is required.');
assert(manifest.options.nMC==500 && manifest.options.nSensitivity==100 && ...
    manifest.options.nLegacy==500 && manifest.options.runLegacy, ...
    'NCGallery:Counts','This gallery requires the full default design.');
if isfolder(outputDir)
    items=dir(outputDir); items=items(~ismember({items.name},{'.','..'}));
    assert(isempty(items),'NCGallery:Exists','Use a new or empty gallery folder.');
else
    mkdir(outputDir);
end
main=fullfile(runDir,'09_publication_artifacts');
names={'Sim_TC_vs_NC_results','Sim_TCH_vs_NC_results', ...
    'Test_NC_Sensitivity_results','Sim_WA_Comparison', ...
    'Figure_R1_Variance_RMSE_alpha_0','Figure_R1_Variance_RMSE_alpha_0p3', ...
    'Figure_R2_Heldout_Merging','Figure_R3_Stress_Feasibility_and_PSD_Cost', ...
    'Figure_R4_Initialization_and_Eigenfloor','Figure_R4b_Initialization', ...
    'Figure_R5_Merging_Controls','Figure_R6_Paired_Seasonal_Change', ...
    'Graphical_Abstract_NC','FigS_NC_PSD_Diagnostic'};
legacyStages={'02_legacy_TC_NC','03_legacy_TCH_NC','04_legacy_sensitivity','05_legacy_WA'};
files={}; layouts=struct('figure',{},'legends',{},'axes',{});
for k=1:numel(names)
    if k<=4
        source=fullfile(runDir,legacyStages{k},'legacy_figure_01.fig');
    elseif k==14
        source=fullfile(runDir,'figures','Supplementary_NC_PSD_Diagnostic','figures',[names{k} '.fig']);
    else
        source=fullfile(main,[names{k} '.fig']);
    end
    assert(isfile(source),'NCGallery:Missing','Missing figure: %s',source);
    if k==13
        % Preserve the verified original exports. Re-exporting a reopened FIG
        % can clip long text in MATLAB R2025b even when its PNG/PDF are correct.
        for ext={'.png','.pdf'}
            copyfile(fullfile(main,[names{k} ext{1}]),fullfile(outputDir,[names{k} ext{1}]));
        end
        layouts(end+1)=struct('figure',names{k},'legends',0,'axes',1); %#ok<AGROW>
        files=[files,{[names{k} '.png'],[names{k} '.pdf']}]; %#ok<AGROW>
        fprintf('Gallery %d/%d: %s (original MATLAB export)\n',k,numel(names),names{k});
        continue
    end
    f=openfig(source,'invisible'); guard=onCleanup(@()close(f));
    axesList=findall(f,'Type','axes');
    for ax=axesList'
        if isprop(ax,'Toolbar') && ~isempty(ax.Toolbar), ax.Toolbar.Visible='off'; end
        disableDefaultInteractivity(ax);
        if k~=13
            ax.FontName='Arial'; ax.FontSize=12;
            ax.Title.FontName='Arial'; ax.Title.FontSize=13;
            if k<=4, ax.Title.String=replace(string(ax.Title.String),' | ',newline); end
            ax.XLabel.FontName='Arial'; ax.XLabel.FontSize=12;
            ax.YLabel.FontName='Arial'; ax.YLabel.FontSize=12;
        end
    end
    legends=findall(f,'Type','legend');
    [~,order]=sort(arrayfun(@(lg)numel(lg.String),legends),'descend');
    legends=legends(order);
    % Keep the richest key when other panels repeat its method labels.
    keys={};
    for lg=legends'
        labels=regexprep(string(lg.String),'\s+(IQR|median)$','','ignorecase');
        key=strjoin(sort(unique(labels)),'|');
        if any(strcmp(keys,key)), delete(lg); else
            keys{end+1}=char(key); %#ok<AGROW>
            if k~=13, lg.FontName='Arial'; lg.FontSize=11; end
        end
    end
    drawnow;
    exportgraphics(f,fullfile(outputDir,[names{k} '.png']), ...
        'Resolution',600,'BackgroundColor','white','Padding',12);
    exportgraphics(f,fullfile(outputDir,[names{k} '.pdf']), ...
        'ContentType','vector','BackgroundColor','white','Padding',12);
    layouts(end+1)=struct('figure',names{k},'legends',numel(findall(f,'Type','legend')), ...
        'axes',numel(axesList)); %#ok<AGROW>
    files=[files,{[names{k} '.png'],[names{k} '.pdf']}]; %#ok<AGROW>
    fprintf('Gallery %d/%d: %s\n',k,numel(names),names{k});
    clear guard
end
summaryDir=fullfile(outputDir,'tables'); mkdir(summaryDir);
for p=dir(fullfile(main,'Table_*.csv'))'
    copyfile(fullfile(p.folder,p.name),fullfile(summaryDir,p.name));
end
diagDir=fullfile(runDir,'figures','Supplementary_NC_PSD_Diagnostic','tables');
for p=dir(fullfile(diagDir,'*.csv'))'
    copyfile(fullfile(p.folder,p.name),fullfile(summaryDir,p.name));
end
copyfile(fullfile(runDir,'01_core_checks','core_regression_results.csv'),summaryDir);
copyfile(fullfile(main,'synthetic_results.md'),fullfile(outputDir,'synthetic_results.md'));
configuration=jsondecode(fileread(fullfile(runDir,'06_revision_synthetic','run_configuration.json')));
writeJSON(fullfile(outputDir,'experiment_configuration.json'),configuration.config);
status=readtable(fullfile(runDir,'stage_status.csv'),'Delimiter',',', ...
    'ReadVariableNames',true,'VariableNamingRule','preserve','TextType','string');
writetable(status(:,{'Stage','Status','Seconds'}),fullfile(summaryDir,'stage_timings.csv'));
sources=rmfield(manifest.sourceFiles,'sourcePath');
for k=1:numel(sources), sources(k).relativePath=strrep(sources(k).relativePath,'\','/'); end
report=struct('runName',getName(runDir),'generatedUTC',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
    'productionReplicationCounts',true,'mainReplicates',500,'originalDesignReplicates',500, ...
    'sensitivityReplicates',100,'diagnosticPixels',1000,'matlabVersion',manifest.matlabVersion, ...
    'figures',{names},'formats',{{'600 dpi PNG (graphical abstract: original 300 dpi)','vector PDF'}},'sourceFiles',sources, ...
    'layout',layouts,'stageStatus',{cellstr(status.Status)},'generatedFiles',{files});
report.exportSources=struct('file',{'Export_Public_Figures.m','Figure_Supplementary_NC_PSD_Diagnostic.m'}, ...
    'sha256',{fileHash(fullfile(root,'Export_Public_Figures.m')),fileHash(fullfile(root,'Figure_Supplementary_NC_PSD_Diagnostic.m'))});
writeJSON(fullfile(outputDir,'run_metadata.json'),report);
fid=fopen(fullfile(outputDir,'README.md'),'w','n','UTF-8'); assert(fid>=0); g=onCleanup(@()fclose(fid));
fprintf(fid,'# MATLAB-generated NC figures\n\nFull current-code execution: 500 main and original-design replicates; 100 sensitivity replicates; 1,000 supplementary diagnostic pixels. MATLAB %s.\n\n',manifest.matlabVersion);
fprintf(fid,'All figures were calculated/exported in MATLAB. PNGs use 600 dpi, except the graphical abstract retains its original 300-dpi MATLAB export; PDFs are vector exports. Numeric data and plot series are unchanged during gallery export; duplicate identical legends and axes toolbars are removed. The graphical abstract retains its original layout.\n\n');
fprintf(fid,'These are outputs of the current NC implementation, not a replacement of the earlier submitted manuscript results. Re-running the original designs with current NC can change their values. The PSD diagnostic uses the estimation samples; paired merging uses independent test samples.\n\n');
fprintf(fid,'Reproduce with `Run_All_Synthetic([],struct(''mode'',''run''))`, then change to the completed run directory and run `Figure_Supplementary_NC_PSD_Diagnostic`. Call `Export_Public_Figures(runDir,newEmptyGalleryDirectory)` to rebuild this gallery. Full workspaces and editable FIG files remain locally under `results/`. Compact tables are in `tables/`; seeds/settings and source hashes are in the JSON files.\n\n');
for k=1:numel(names)
    fprintf(fid,'## %s\n\n[PDF](%s.pdf)\n\n![%s](%s.png)\n',strrep(names{k},'_',' '),names{k},strrep(names{k},'_',' '),names{k});
    if k<numel(names), fprintf(fid,'\n'); end
end
end

function name=getName(folder)
[~,name]=fileparts(folder);
end
function writeJSON(file,data)
fid=fopen(file,'w','n','UTF-8'); assert(fid>=0); g=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(data,'PrettyPrint',true));
end

function hash=fileHash(file)
fid=fopen(file,'rb'); assert(fid>=0); bytes=fread(fid,Inf,'*uint8'); fclose(fid);
digest=java.security.MessageDigest.getInstance('SHA-256'); digest.update(typecast(bytes,'int8'));
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
