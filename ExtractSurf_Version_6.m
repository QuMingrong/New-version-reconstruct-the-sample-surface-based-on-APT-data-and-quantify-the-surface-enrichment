%% ExtractSurf_Version_6.m
% HEA surface peeling + pre-analysis + Ir-centered statistics + layer-selectable 1NN
% 功能：
%  - 选 HEA 点云，拟合上表面并按法向距离剥离多层
%  - 导出各层点云和 HEA 统计
%  - 1 nm 富集 / KDE / Posterior / E_s(r)
%  - kNN 距离频率（计数） + PPP baseline
%  - 1NN^s（中心元素 -> 其它 HEA 元素）
%  - 1NN multi vs union（多元素 vs 合并 baseline）
%  - New: GUI 中可选 “中心原子层数 + 近邻原子层数”

clear; clc;

%% 0) Select files
[aptFile, aptPath] = uigetfile({'*.csv'}, 'Select APT point cloud file (apt.csv)');
if isequal(aptFile,0), error('Cancelled.'); end
[mzwFile, mzwPath] = uigetfile({'*.xlsx;*.csv'}, 'Select m/z window table (APT_mz_windows.xlsx)');
if isequal(mzwFile,0), error('Cancelled.'); end
outDefault = fullfile(aptPath, 'HEA_surface_output.xlsx');
[outFile, outPath] = uiputfile({'*.xlsx'}, 'Save surface Excel as...', outDefault);
if isequal(outFile,0), error('Cancelled.'); end
outXLSX = fullfile(outPath, outFile);

%% 1) Parameters (will be set after Pre-analysis for h)
h_user     = NaN;   % XY grid step (nm) suggested by Pre-analysis
q_user     = 0.98;  % quantile for upper envelope
Nmin_user  = 20;    % min HEA points per XY cell
layer_order = "topmost";  % or 'closest'
HEAset      = ["Ir","Ru","Rh","Pt","Pd"];   % global HEA list

%% 2) Load data + map m/z -> species
T = readtable(fullfile(aptPath, aptFile), ...
    'VariableNamingRule','preserve', 'TextType','string');
cn = matlab.lang.makeUniqueStrings(strtrim(string(T.Properties.VariableNames)));
T.Properties.VariableNames = cn;
xCol = find(ismember(upper(cn), ["X","CX","COORDX"]), 1);
yCol = find(ismember(upper(cn), ["Y","CY","COORDY"]), 1);
zCol = find(ismember(upper(cn), ["Z","CZ","COORDZ"]), 1);
mzCol = find(contains(lower(cn),["mz","m/z","mass","tocharge"]), 1);
assert(~(isempty(xCol)||isempty(yCol)||isempty(zCol)||isempty(mzCol)), ...
    'X/Y/Z/mz not detected.');
T.X  = double(T{:,xCol});
T.Y  = double(T{:,yCol});
T.Z  = double(T{:,zCol});
T.mz = double(T{:,mzCol});

[~,~,extW] = fileparts(mzwFile);
if strcmpi(extW,'.csv')
    W = readtable(fullfile(mzwPath, mzwFile), ...
        'VariableNamingRule','preserve', 'TextType','string');
else
    W = readtable(fullfile(mzwPath, mzwFile),'FileType','spreadsheet', ...
        'VariableNamingRule','preserve','TextType','string');
end
W.Properties.VariableNames = matlab.lang.makeUniqueStrings( ...
    strtrim(string(W.Properties.VariableNames)));
assert(all(ismember(["Species","mz_min","mz_max"], ...
    string(W.Properties.VariableNames))), ...
    'Window table must have: Species, mz_min, mz_max');

label = assignSpeciesFromWindows(T.mz, W);
T.Species = categorical(label);

%% 3) Pre-analysis (KNN & NND) to suggest XY step h, with plots
[h_user, method_tag, preStats] = preAnalysisDialog(T, HEAset, outXLSX); %#ok<NASGU>
if isnan(h_user)
    error('Cancelled at Pre-analysis.');
end

%% 4) Parameters dialog
prompt = {'XY grid step h (nm):', ...
          'Quantile q (0–1, 1=strict max):', ...
          'Min HEA points / XY cell (Nmin):', ...
          'Layer-1 definition (topmost/closest):'};
defAns = {sprintf('%.3f', h_user), sprintf('%.2f', q_user), ...
          sprintf('%d', Nmin_user), char(layer_order)};
answ = inputdlg(prompt,'HEA-only Surface Parameters',1,defAns);
if isempty(answ), error('Cancelled.'); end
h    = str2double(answ{1});
q    = str2double(answ{2});
Nmin = str2double(answ{3});
layer_order = lower(strtrim(answ{4}));
if ~ismember(layer_order,["topmost","closest"])
    layer_order = "topmost";
end

%% 5) Build HEA-only upper envelope (surface) + normals
Th = T(ismember(string(T.Species), HEAset), :);
assert(~isempty(Th), 'No HEA (Ir/Ru/Rh/Pt/Pd) points after filtering.');
XYZ = [Th.X, Th.Y, Th.Z];

xb = floor(min(XYZ(:,1))/h)*h : h : ceil(max(XYZ(:,1))/h)*h;
yb = floor(min(XYZ(:,2))/h)*h : h : ceil(max(XYZ(:,2))/h)*h;
nx = numel(xb)-1; ny = numel(yb)-1;

zq  = nan(nx,ny);  nxG = nan(nx,ny);  nyG = nan(nx,ny);  nzG = nan(nx,ny);
for ix = 1:nx
  for iy = 1:ny
    in = XYZ(:,1)>=xb(ix) & XYZ(:,1)<xb(ix+1) & ...
         XYZ(:,2)>=yb(iy) & XYZ(:,2)<yb(iy+1);
    if nnz(in) < Nmin, continue; end
    z_local   = XYZ(in,3);
    zq(ix,iy) = quantile(z_local, q);       % upper envelope per cell
    Xi = [XYZ(in,1), XYZ(in,2), ones(nnz(in),1)];
    ai = Xi \ z_local;                      % plane fit z ~ ax+by+c
    n  = [-ai(1), -ai(2), 1]; n = n/norm(n);
    nxG(ix,iy)=n(1); nyG(ix,iy)=n(2); nzG(ix,iy)=n(3);
  end
end
[Xg, Yg] = ndgrid(xb(1:end-1)+h/2, yb(1:end-1)+h/2);
Zg = zq;

%% 6) Signed-normal distance for ALL atoms; per-cell HEA ranking
FnZ = griddedInterpolant({xb(1:end-1)+h/2, yb(1:end-1)+h/2}, Zg,  'nearest','nearest');
Fnx= griddedInterpolant({xb(1:end-1)+h/2, yb(1:end-1)+h/2}, nxG, 'nearest','nearest');
Fny= griddedInterpolant({xb(1:end-1)+h/2, yb(1:end-1)+h/2}, nyG, 'nearest','nearest');
Fnz= griddedInterpolant({xb(1:end-1)+h/2, yb(1:end-1)+h/2}, nzG, 'nearest','nearest');

x=T.X; y=T.Y; z=T.Z;
zsurf=FnZ(x,y); nxv=Fnx(x,y); nyv=Fny(x,y); nzv=Fnz(x,y);
valid = isfinite(zsurf)&isfinite(nxv)&isfinite(nyv)&isfinite(nzv);

nvec = [nxv(valid), nyv(valid), nzv(valid)];
nvec = nvec ./ vecnorm(nvec,2,2);
sref = [x(valid), y(valid), zsurf(valid)];
d    = sum( ([x(valid),y(valid),z(valid)] - sref) .* nvec, 2 );

nxg = size(Zg,1); nyg = size(Zg,2);
ix_all = max(1, min(nxg, floor((x - xb(1))/h)+1));
iy_all = max(1, min(nyg, floor((y - yb(1))/h)+1));
Ivalid = find(valid);

[cellOrderByD, maxRankAvailable] = buildCellOrder(T, Ivalid, ix_all, iy_all, d, HEAset, layer_order);

% 构建 L{1..K_pre} 层 mask
K_pre = min(5, maxRankAvailable);
L = repmat({false(height(T),1)}, 1, K_pre);
for k = 1:K_pre
    L{k} = mask_from_rank(cellOrderByD, k, height(T));
end
L_top5 = false(height(T),1);
for k=1:K_pre, L_top5 = L_top5 | L{k}; end

% 把关键变量放到 base 里供 GUI 使用
assignin('base','T',T);
assignin('base','outXLSX',outXLSX);
assignin('base','HEAset',HEAset);
assignin('base','layer_order',string(layer_order));
assignin('base','cellOrderByD',cellOrderByD);
assignin('base','maxRankAvailable',maxRankAvailable);
assignin('base','L',L);
assignin('base','L_top5',L_top5);
assignin('base','Xg',Xg); assignin('base','Yg',Yg); assignin('base','Zg',Zg);

%% 7) Export default layers (1..5 or up to available) + Z/D stats
colHex = struct('Ir','#FF0000','Pt','#0000FF','Ru','#FF00FF','Pd','#00CCFF','Rh','#000000');

maskGrid = isfinite(Zg);
writetable(table(Xg(maskGrid),Yg(maskGrid),Zg(maskGrid), ...
    'VariableNames',{'X','Y','Z'}), ...
    outXLSX,'Sheet','SurfaceGrid');

Kexp = K_pre;
Layer = (1:Kexp)'; Count=zeros(Kexp,1); MeanZ=nan(Kexp,1); MedianZ=nan(Kexp,1);
StdZ=nan(Kexp,1); MinZ=nan(Kexp,1); MaxZ=nan(Kexp,1);
d_full = nan(height(T),1); d_full(Ivalid)=d;
MeanD=nan(Kexp,1); MedianD=nan(Kexp,1); StdD=nan(Kexp,1);

for k=1:Kexp
    Mk = mask_from_rank(cellOrderByD, k, height(T));
    Tk = T(Mk & ismember(string(T.Species),HEAset), {'X','Y','Z','Species'});
    Count(k)=height(Tk);
    if Count(k)>0
        MeanZ(k)=mean(Tk.Z,'omitnan'); MedianZ(k)=median(Tk.Z,'omitnan');
        StdZ(k)=std(Tk.Z,0,'omitnan');  MinZ(k)=min(Tk.Z);  MaxZ(k)=max(Tk.Z);
    end
    dk = d_full(Mk & ismember(string(T.Species),HEAset));
    MeanD(k)=mean(dk,'omitnan'); MedianD(k)=median(dk,'omitnan'); StdD(k)=std(dk,0,'omitnan');

    % 元素 + ColorHex
    colorHexCol = strings(height(Tk),1);
    for i=1:height(Tk)
        s = char(string(Tk.Species(i)));
        if isfield(colHex,s), colorHexCol(i)=colHex.(s); else, colorHexCol(i)="#808080"; end
    end
    Tk.ColorHex = colorHexCol;
    writetable(Tk, outXLSX, 'Sheet', sprintf('SurfaceLayer%d',k));

    if height(Tk)>0
        S = groupsummary(Tk,'Species');
        S.Properties.VariableNames = {'Species','GroupCount'};
        S.Fraction = S.GroupCount / sum(S.GroupCount);
    else
        S = table(categorical([]),[],[], ...
            'VariableNames',{'Species','GroupCount','Fraction'});
    end
    writetable(S, outXLSX, 'Sheet', sprintf('Layer%d_Summary',k));
end
writetable(table(Layer,Count,MeanZ,MedianZ,StdZ,MinZ,MaxZ), ...
    outXLSX, 'Sheet','LayerZ_Stats');
writetable(table(Layer,MeanD,MedianD,StdD), ...
    outXLSX, 'Sheet','LayerD_Stats');

%% 8) Plots — surface and layers (简单 3D 可视化)
figure('Color','w');
surf(Xg,Yg,Zg,'EdgeColor','none','FaceAlpha',0.95);
axis equal tight; view(3); camlight headlight; lighting gouraud;
xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (nm)');
title('HEA-only Surface (mesh only)'); colorbar;

for k=1:Kexp
    Mk = mask_from_rank(cellOrderByD, k, height(T));
    Tk = T(Mk & ismember(string(T.Species),HEAset), :);
    figure('Color','w'); hold on;
    legendNames = {};
    for s = HEAset
        in = Tk.Species == categorical(s);
        if any(in)
            c = hex2rgb(colHex.(char(s)));
            scatter3(Tk.X(in),Tk.Y(in),Tk.Z(in),10,'filled', ...
                     'MarkerFaceColor',c,'MarkerEdgeColor','none');
            legendNames{end+1} = char(s); %#ok<AGROW>
        end
    end
    axis equal tight; grid on; view(3);
    xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (nm)');
    title(sprintf('Surface Layer %d (HEA only)', k));
    if ~isempty(legendNames), legend(legendNames,'Location','northeastoutside'); end
end

Ntop = min(2, maxRankAvailable);
Mtop = false(height(T),1);
for k=1:Ntop, Mtop = Mtop | mask_from_rank(cellOrderByD,k,height(T)); end
Ts = T(Mtop & ismember(string(T.Species),HEAset),:);
figure('Color','w'); hold on;
legendNames = {};
for s = HEAset
    in = Ts.Species == categorical(s);
    if any(in)
        c = hex2rgb(colHex.(char(s)));
        scatter3(Ts.X(in),Ts.Y(in),Ts.Z(in),8,'filled', ...
                 'MarkerFaceColor',c,'MarkerEdgeColor','none');
        legendNames{end+1}=char(s);
    end
end
axis equal tight; grid on; view(3);
xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (nm)');
title(sprintf('Surface Layers 1-%d (HEA only, element-colored)', Ntop));
if ~isempty(legendNames), legend(legendNames,'Location','northeastoutside'); end

%% 9) Open interactive app (layer-selectable enrichment + kNN + 1NN)
createSurfaceApp();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ===================== Helper: Pre-analysis dialog =====================
function [h_suggest, method_tag, S] = preAnalysisDialog(T, HEAset, outXLSX)
    c3 = 4*pi/3;
    h_suggest = NaN;
    method_tag = "";
    S = struct();

    ui0.f = uifigure('Name','Pre-analysis: KNN & NND (Suggest XY step h)',...
                     'Position',[80 80 1060 660]);

    % ---- Top controls ----
    uilabel(ui0.f,'Text','Universe','Position',[20 620 80 20]);
    ui0.universe = uidropdown(ui0.f,'Items',{'HEA only','All species'}, ...
        'Value','HEA only','Position',[90 615 160 28]);

    uilabel(ui0.f,'Text','Max samples (subsample)','Position',[270 620 160 20]);
    ui0.maxN = uieditfield(ui0.f,'numeric','Limits',[1000 Inf],...
        'Value',250000,'Position',[430 615 110 28]);

    uilabel(ui0.f,'Text','Method to set h','Position',[560 620 130 20]);
    ui0.method = uidropdown(ui0.f,'Items',...
        {'KNN: min(ρ^{-1/3}, median R5)', ...
         'NND: MLE ( (c3·mean(R1^3))^{1/3} )', ...
         'NND: Median ( median(R1)/0.549 )', ...
         'Custom'}, ...
        'Value','KNN: min(ρ^{-1/3}, median R5)',...
        'Position',[680 615 300 28]);

    uilabel(ui0.f,'Text','Custom h (nm)','Position',[20 585 100 20]);
    ui0.hCustom = uieditfield(ui0.f,'numeric','Limits',[0.02 5],...
        'Value',0.35,'Enable','off','Position',[120 580 90 28]);

    ui0.btnRun = uibutton(ui0.f,'Text','Compute / Refresh',...
        'Position',[230 580 150 30],'ButtonPushedFcn',@onCompute);

    ui0.btnAccept = uibutton(ui0.f,'Text','Accept & Continue',...
        'Position',[390 580 150 30],'ButtonPushedFcn',@onAccept);

    ui0.btnExport = uibutton(ui0.f,'Text','Export Pre-analysis to Excel',...
        'Position',[550 580 220 30],'ButtonPushedFcn',@onExport);

    ui0.btnCancel = uibutton(ui0.f,'Text','Cancel',...
        'Position',[780 580 90 30],'ButtonPushedFcn',@(src,evt) delete(ui0.f));

    % ---- Tables / text ----
    uilabel(ui0.f,'Text','KNN summary','FontWeight','bold',...
        'Position',[20 540 120 20]);
    ui0.tabK = uitable(ui0.f,'Position',[20 410 510 120], ...
        'ColumnName',{'Metric','Value'}, 'ColumnEditable',[false false]);

    uilabel(ui0.f,'Text','NND summary','FontWeight','bold',...
        'Position',[20 380 120 20]);
    ui0.tabN = uitable(ui0.f,'Position',[20 250 510 120], ...
        'ColumnName',{'Metric','Value'}, 'ColumnEditable',[false false]);

    uilabel(ui0.f,'Text','Suggested h (nm)','FontWeight','bold',...
        'Position',[20 220 120 20]);
    ui0.hShow = uieditfield(ui0.f,'text','Editable','off',...
        'Position',[140 215 130 28]);

    % ---- Axes ----
    uilabel(ui0.f,'Text','R5 distribution (KNN)','FontWeight','bold',...
        'Position',[560 540 200 20]);
    ui0.axK = uiaxes(ui0.f,'Position',[560 340 470 200]);
    title(ui0.axK, 'R5 (nm) — histogram + KDE');
    xlabel(ui0.axK,'R5 (nm)'); ylabel(ui0.axK,'Density');

    uilabel(ui0.f,'Text','R1 (NND) distribution','FontWeight','bold',...
        'Position',[560 300 200 20]);
    ui0.axN = uiaxes(ui0.f,'Position',[560 100 470 200]);
    title(ui0.axN, 'R1 (nm) — histogram + KDE');
    xlabel(ui0.axN,'R1 (nm)'); ylabel(ui0.axN,'Density');

    ui0.method.ValueChangedFcn = @(src,evt) ...
        set(ui0.hCustom,'Enable', ternary(strcmp(src.Value,'Custom'),'on','off'));

    drawnow;
    onCompute();   % 先算一次

    % 关键：用 waitfor 等 figure 被 delete
    waitfor(ui0.f);

    tmp = [];
    if isappdata(0,'PREANALYSIS_RESULT')
        tmp = getappdata(0,'PREANALYSIS_RESULT');
        rmappdata(0,'PREANALYSIS_RESULT');
    end
    if ~isempty(tmp)
        h_suggest = tmp.h_suggest;
        method_tag = tmp.method_tag;
        S = tmp.S;
    end

    % ---- 内部函数 ----
    function onCompute(~,~)
        if strcmp(ui0.universe.Value, 'HEA only')
            mask = ismember(string(T.Species), HEAset);
        else
            mask = true(height(T),1);
        end
        X = [double(T.X(mask)), double(T.Y(mask)), double(T.Z(mask))];
        X = X(all(isfinite(X),2),:);
        X = unique(X,'rows','stable');
        N = size(X,1);

        if N < 2000
            maxN = N;
        else
            maxN = min(N, max(1000, round(ui0.maxN.Value)));
        end
        if N > maxN
            idx = randperm(N, maxN);
        else
            idx = 1:N;
        end
        Xs = X(idx,:);
        Ns = size(Xs,1);

        Mdl = createns(Xs,'NSMethod','kdtree');
        [~, D6] = knnsearch(Mdl, Xs, 'K', 6);
        R5 = D6(:,6);
        medR5 = median(R5,'omitnan');

        [~, D2] = knnsearch(Mdl, Xs, 'K', 2);
        R1 = D2(:,2);
        meanR1 = mean(R1,'omitnan');
        medR1  = median(R1,'omitnan');

        try
            [~, Vh] = convhulln(X);
        catch
            mins = min(X,[],1);
            maxs = max(X,[],1);
            Vh   = prod(maxs - mins);
        end
        rho = size(X,1) / max(Vh, eps);
        s_voxel = rho^(-1/3);

        h_knn     = min(s_voxel, medR5);
        s_nnd_mle = (c3 * mean(R1.^3,'omitnan'))^(1/3);
        s_nnd_med = medR1 / 0.549;

        S = struct('Universe', ui0.universe.Value, ...
                   'N_total', size(X,1), ...
                   'N_sample', Ns, ...
                   'rho', rho, ...
                   's_voxel', s_voxel, ...
                   'medianR5', medR5, ...
                   'h_knn', h_knn, ...
                   'meanR1', meanR1, ...
                   'medianR1', medR1, ...
                   's_nnd_mle', s_nnd_mle, ...
                   's_nnd_med', s_nnd_med);

        Ktbl = {
            'Total points (selected)',   num2str(S.N_total);
            'Sampled points (KD-tree)',  num2str(S.N_sample);
            'Density ρ (atoms/nm^3)',    sprintf('%.4f', S.rho);
            'Voxel scale ρ^{-1/3} (nm)', sprintf('%.3f', S.s_voxel);
            'median R5 (nm)',            sprintf('%.3f', S.medianR5);
            'Suggested h (KNN)',         sprintf('%.3f', S.h_knn);
            };
        ui0.tabK.Data = Ktbl;

        Ntbl = {
            'mean R1 (nm)',     sprintf('%.3f', S.meanR1);
            'median R1 (nm)',   sprintf('%.3f', S.medianR1);
            's_{NND,MLE} (nm)', sprintf('%.3f', S.s_nnd_mle);
            's_{NND,Median} (nm)', sprintf('%.3f', S.s_nnd_med);
            };
        ui0.tabN.Data = Ntbl;

        cla(ui0.axK); hold(ui0.axK,'on');
        if ~isempty(R5)
            histogram(ui0.axK, R5, 'Normalization','pdf',...
                'EdgeColor','none','FaceAlpha',0.35);
            try
                [f,xx] = ksdensity(R5);
                plot(ui0.axK, xx, f, 'LineWidth', 1.8);
            catch
            end
            xline(ui0.axK, medR5,'--','LineWidth',1.2);
            xline(ui0.axK, s_voxel,'-.','LineWidth',1.2);
            legend(ui0.axK, {'hist (pdf)','KDE','median R5','ρ^{-1/3}'},...
                   'Location','northeast');
        end
        hold(ui0.axK,'off');

        cla(ui0.axN); hold(ui0.axN,'on');
        if ~isempty(R1)
            histogram(ui0.axN, R1, 'Normalization','pdf',...
                'EdgeColor','none','FaceAlpha',0.35);
            try
                [f2,xx2] = ksdensity(R1);
                plot(ui0.axN, xx2, f2, 'LineWidth', 1.8);
            catch
            end
            xline(ui0.axN, medR1,'--','LineWidth',1.2);
            xline(ui0.axN, s_nnd_mle,'-.','LineWidth',1.2);
            legend(ui0.axN, {'hist (pdf)','KDE','median R1','s_{NND,MLE}'},...
                   'Location','northeast');
        end
        hold(ui0.axN,'off');

        curMethod = ui0.method.Value;
        switch curMethod
            case 'KNN: min(ρ^{-1/3}, median R5)'
                h_now = h_knn; method_tag = 'KNN';
            case 'NND: MLE ( (c3·mean(R1^3))^{1/3} )'
                h_now = s_nnd_mle; method_tag = 'NND_MLE';
            case 'NND: Median ( median(R1)/0.549 )'
                h_now = s_nnd_med; method_tag = 'NND_Median';
            otherwise
                h_now = ui0.hCustom.Value; method_tag = 'Custom';
        end
        ui0.hShow.Value = sprintf('%.3f', h_now);

        setappdata(ui0.f,'PRE_S',S);
        setappdata(ui0.f,'PRE_hNow',h_now);
        setappdata(ui0.f,'PRE_method',curMethod);
        drawnow;
    end

    function onAccept(~,~)
        S  = getappdata(ui0.f,'PRE_S');
        hN = getappdata(ui0.f,'PRE_hNow');
        mth= getappdata(ui0.f,'PRE_method');
        setappdata(0,'PREANALYSIS_RESULT', ...
            struct('h_suggest',hN,'method_tag',mth,'S',S));
        if isvalid(ui0.f), delete(ui0.f); end
    end

    function onExport(~,~)
        try
            Slocal  = getappdata(ui0.f,'PRE_S');
            if isempty(Slocal), return; end
            Stab = struct2table(Slocal,'AsArray',true);
            writetable(Stab, outXLSX, 'Sheet','PreAnalysis_Summary');
            uialert(ui0.f, 'Pre-analysis summary exported to Excel.', 'Export');
        catch ME
            uialert(ui0.f, sprintf('Export failed:\n%s', ME.message), 'Export Error');
        end
    end
end

%% ================== Helper functions (mapping/colors/etc.) =============
function label = assignSpeciesFromWindows(mz, W)
    sp = string(W.Species);
    mm = [double(W.mz_min), double(W.mz_max)];
    mz = mz(:);
    label = strings(size(mz));
    speciesList = unique(sp,'stable');
    for k = 1:numel(speciesList)
        s = speciesList(k); rows = find(sp==s);
        mask = false(size(mz));
        for r = 1:numel(rows)
            lo = mm(rows(r),1); hi = mm(rows(r),2);
            mask = mask | (mz>=lo & mz<=hi);
        end
        assignable = (label=="") & mask;
        label(assignable) = s;
    end
    label(label=="") = "Unknown";
end

function rgb = hex2rgb(hex)
    hex = char(hex); 
    if hex(1) == '#', hex = hex(2:end); end
    rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))]/255;
end

function [cellOrderByD, maxRankAvailable] = buildCellOrder(T, Ivalid, ix_all, iy_all, d, HEAset, layer_order)
    d_min = -Inf;
    isHEA_all   = ismember(string(T.Species), HEAset);
    isHEA_valid = isHEA_all(Ivalid);
    keep = isHEA_valid & (d >= d_min);

    nx = max(ix_all); ny = max(iy_all);
    cellID_valid = sub2ind([nx,ny], ix_all(Ivalid), iy_all(Ivalid));

    uCells = unique(cellID_valid(keep));
    G_idx = []; G_d = []; G_cell = []; G_rank = [];

    for c = reshape(uCells,1,[])
        J = find(cellID_valid == c & keep);
        if isempty(J), continue; end
        switch lower(layer_order)
            case 'closest'
                [~,ord] = sort(d(J), 'ascend');
            otherwise
                [~,ord] = sort(d(J), 'descend');
        end
        J = J(ord);
        ranks = (1:numel(J))';
        G_idx  = [G_idx;  Ivalid(J)];            %#ok<AGROW>
        G_d    = [G_d;    d(J)];                 %#ok<AGROW>
        G_cell = [G_cell; repmat(c,numel(J),1)]; %#ok<AGROW>
        G_rank = [G_rank; ranks];                %#ok<AGROW>
    end

    cellOrderByD.idxGlobal = G_idx;
    cellOrderByD.d         = G_d;
    cellOrderByD.cellID    = G_cell;
    cellOrderByD.rank      = G_rank;

    maxRankAvailable = 0;
    if ~isempty(G_rank), maxRankAvailable = max(G_rank); end
end

function M = mask_from_rank(cellOrderByD, k, N)
    M = false(N,1);
    if isempty(fieldnames(cellOrderByD)), return; end
    rows = (cellOrderByD.rank == k);
    if any(rows), M(cellOrderByD.idxGlobal(rows)) = true; end
end

function [r, b] = fcc_shell_params(d1, Smax, epsR)
    Smax = min(Smax,5);
    a = sqrt(2)*d1;
    rFCC = [a/sqrt(2), a, a*sqrt(3/2), a*sqrt(2), a*sqrt(5/2)];
    r = rFCC(1:Smax).';
    nextR = a*sqrt(3);
    b = zeros(Smax,1);
    if Smax>=2
        b(1) = (r(1)+r(2))/2;
    else
        b(1) = r(1)+0.5*(nextR-r(1));
    end
    for s=2:Smax-1
        b(s) = (r(s)+r(s+1))/2;
    end
    b(end) = (r(end)+nextR)/2;
    b = b + epsR;
end

function [tbl, pairCounts] = kde_pdf_pairs(T, centerMask, neighborsMask, neighList, r_min, r_max, bw, nGrid)
    neighList = string(neighList(:)');
    centers = find(centerMask);
    tbl = table(string.empty, [], [], [], ...
        'VariableNames', {'NeighborSpecies','r','pdf_kde','Npairs'});
    pairCounts = zeros(1, numel(neighList));
    if isempty(centers), return; end

    XYZ=[T.X,T.Y,T.Z];
    nbrIdx=find(neighborsMask);
    if isempty(nbrIdx), return; end
    Mdl=createns(XYZ(nbrIdx,:), 'NSMethod','kdtree');

    D = cell(numel(neighList),1);
    idxCell = rangesearch(Mdl, XYZ(centers,:), r_max);
    for ii=1:numel(centers)
        c=centers(ii); glb=nbrIdx(idxCell{ii}); glb(glb==c)=[]; 
        if isempty(glb), continue; end
        dist=vecnorm(XYZ(glb,:)-XYZ(c,:),2,2);
        sp=string(T.Species(glb));
        sel=(dist>=r_min & dist<=r_max); dist=dist(sel); sp=sp(sel);
        for k=1:numel(neighList)
            if ~isempty(dist), D{k} = [D{k}; dist(sp==neighList(k))]; end %#ok<AGROW>
        end
    end

    r_grid = linspace(r_min, r_max, nGrid).';
    for k=1:numel(neighList)
        d = D{k}; pairCounts(k)=numel(d);
        if isempty(d)
            f=zeros(nGrid,1); Npairs=0;
        else
            try
                [f,~]=ksdensity(d, r_grid, 'Bandwidth', bw, ...
                                'Support', [r_min r_max], ...
                                'BoundaryCorrection','reflection');
            catch
                d_ext=[2*r_min-d; d; 2*r_max-d];
                [f,~]=ksdensity(d_ext, r_grid, 'Bandwidth', bw, ...
                                'Support', [r_min r_max]);
            end
            A=trapz(r_grid,f); if A>0, f=f/A; end
            Npairs=numel(d);
        end
        tmp=table(repmat(neighList(k),numel(r_grid),1), r_grid, f(:), ...
                  repmat(Npairs,numel(r_grid),1), ...
                  'VariableNames',{'NeighborSpecies','r','pdf_kde','Npairs'});
        tbl=[tbl; tmp]; %#ok<AGROW>
    end
end

%% ============== kNN frequency helpers (PPP baselines; counts) ==========
function [centers, empCounts, baseCounts, edges, meta] = knn_frequency_with_ppp(XYZ_pool, XYZ_centers, k, binSpec)
    P = XYZ_pool(all(isfinite(XYZ_pool),2),:);
    P = unique(P,'rows','stable');
    Q = XYZ_centers(all(isfinite(XYZ_centers),2),:);

    if size(P,1) < (k+1)
        error('Neighbor pool has too few points for k = %d.', k);
    end
    if isempty(Q), error('No centers to query.'); end

    try
        [~, V] = convhulln(P);
    catch
        mins = min(P,[],1); maxs = max(P,[],1);
        V = prod(maxs - mins);
    end
    rho = size(P,1) / max(V, eps);

    Mdl = createns(P, 'NSMethod','kdtree');
    [~, Dk] = knnsearch(Mdl, Q, 'K', k+1);
    Rk = Dk(:, k+1);

    mins = min(P,[],1); maxs = max(P,[],1);
    edge_clear = min( [Q - mins, maxs - Q], [], 2 );
    keep = Rk <= edge_clear;
    Rk_kept = Rk(keep);
    if isempty(Rk_kept)
        error('All R_k removed by edge correction. Try smaller k or wider ROI.');
    end

    switch lower(binSpec.mode)
        case 'fd'
            iq = iqr(Rk_kept);
            bw = 2*iq / max(numel(Rk_kept)^(1/3), 1);
            if ~isfinite(bw) || bw<=0
                nb = max(10, min(80, round(sqrt(numel(Rk_kept)))));
                edges = linspace(min(Rk_kept), max(Rk_kept), nb+1);
            else
                nb = max(10, ceil((max(Rk_kept)-min(Rk_kept))/bw));
                edges = linspace(min(Rk_kept), max(Rk_kept), nb+1);
            end
        otherwise
            nb = max(10, binSpec.nbins);
            edges = linspace(min(Rk_kept), max(Rk_kept), nb+1);
    end
    histCounts = histcounts(Rk_kept, edges);
    empCounts  = histCounts;
    centers = (edges(1:end-1) + edges(2:end)) / 2;

    rho_loc = rho;
    mu = @(r) rho_loc * (4*pi/3) .* (r.^3);
    cdf_edges = zeros(size(edges));
    for i=1:numel(edges)
        m = mu(edges(i));
        s = 0; eg = exp(-m); pow = 1;
        for n=0:k-1
            if n>0, pow = pow * m; end
            s = s + eg * (pow / factorial(n));
        end
        cdf_edges(i) = 1 - s;
    end
    massProb   = diff(cdf_edges);
    baseCounts = massProb * numel(Rk_kept);

    meta = struct('rho',rho,'V',V, 'kept_fraction',mean(keep), ...
                  'Nq',numel(Rk),'Nq_kept',numel(Rk_kept));
end

function CS = crossSpeciesKNNfreq_withPPP(T, centersMask, neighList, roiMask, k, binSpec)
    neighList = string(neighList(:)');
    idxC = find(centersMask & roiMask);
    if isempty(idxC), error('No centers.'); end

    Q = [double(T.X(idxC)), double(T.Y(idxC)), double(T.Z(idxC))];

    P_all = [double(T.X(roiMask)), double(T.Y(roiMask)), double(T.Z(roiMask))];
    if size(P_all,1) < k+1, error('ROI has too few points for k=%d.', k); end
    try
        [~, V] = convhulln(P_all);
    catch
        mins = min(P_all,[],1); maxs = max(P_all,[],1);
        V = prod(maxs - mins);
    end

    mins = min(P_all,[],1); maxs = max(P_all,[],1);
    edge_clear = min([Q - mins, maxs - Q], [], 2);

    D_all = [];
    rho = zeros(1,numel(neighList));
    Ns  = zeros(1,numel(neighList));
    Neff = zeros(1,numel(neighList));
    D_map = cell(1,numel(neighList));

    for j = 1:numel(neighList)
        s = neighList(j);
        idxS = find(roiMask & string(T.Species) == s);
        Ps = [double(T.X(idxS)), double(T.Y(idxS)), double(T.Z(idxS))];
        Ns(j) = size(Ps,1);
        if Ns(j) < k
            D_map{j} = []; rho(j) = 0; Neff(j) = 0;
            continue;
        end
        rho(j) = Ns(j) / max(V, eps);

        Mdl = createns(Ps, 'NSMethod','kdtree');
        [~, Dk] = knnsearch(Mdl, Q, 'K', k);
        Rk = Dk(:, k);
        keep = isfinite(Rk) & (Rk <= edge_clear);
        Rk = Rk(keep);

        D_map{j} = Rk;
        D_all = [D_all; Rk]; %#ok<AGROW>
        Neff(j) = numel(Rk);
    end
    if isempty(D_all), error('No valid center-neighbor pairs after edge correction.'); end

    switch lower(binSpec.mode)
        case 'fd'
            iq = iqr(D_all);
            bw = 2*iq / max(numel(D_all)^(1/3), 1);
            if ~isfinite(bw) || bw<=0
                nb = max(10, min(80, round(sqrt(numel(D_all)))));
            else
                nb = max(10, ceil((max(D_all)-min(D_all))/bw));
            end
        otherwise
            nb = max(10, binSpec.nbins);
    end
    edges = linspace(min(D_all), max(D_all), nb+1);
    centers = (edges(1:end-1)+edges(2:end))/2;

    S = numel(neighList); B = numel(centers);
    empCounts  = zeros(B,S);
    baseCounts = zeros(B,S);

    for j = 1:S
        d = D_map{j};
        if isempty(d), continue; end
        hc = histcounts(d, edges);
        empCounts(:,j) = hc(:);

        rhos = rho(j);
        mu = @(r) rhos*(4*pi/3).*(r.^3);
        cdf_edges = zeros(size(edges));
        for ii=1:numel(edges)
            m = mu(edges(ii));
            ssum = 0; eg = exp(-m); pow = 1;
            for n=0:k-1
                if n>0, pow = pow * m; end
                ssum = ssum + eg * (pow / factorial(n));
            end
            cdf_edges(ii) = 1 - ssum;
        end
        massProb = diff(cdf_edges);
        baseCounts(:,j) = massProb(:) * Neff(j);
    end

    CS.centers     = centers;
    CS.empCounts   = empCounts;
    CS.baseCounts  = baseCounts;
    CS.edges       = edges;
    CS.species     = neighList;
    CS.meta = struct('V',V,'rho',rho,'Ns',Ns,'Neff',Neff,'k',k,'Ncenters',numel(idxC));
end

%% =================== Interactive App (with layer-selectable 1NN) ========
function createSurfaceApp()
    persistent ui;
    if ~isempty(ui) && isvalid(ui.f)
        figure(ui.f);
        return;
    end

    T            = evalin('base','T');
    outXLSX      = evalin('base','outXLSX');
    HEAset       = evalin('base','HEAset');
    layer_order  = evalin('base','layer_order'); %#ok<NASGU>
    cellOrderByD = evalin('base','cellOrderByD');
    maxRankAvail = evalin('base','maxRankAvailable');
    L            = evalin('base','L');

    Ntot = height(T);

    ui.f = uifigure('Name','Surface Layers, Enrichment & kNN Frequency',...
                    'Position',[100 100 1200 760]);

    % -------------------- 顶部：表面 / 层数 & 富集 --------------------

    uilabel(ui.f,'Text','Layer selection mode','Position',[20 720 150 22]);
    ui.layerMode = uidropdown(ui.f,'Items',{'Union 1..N','Only layer N'}, ...
        'Value','Union 1..N','Position',[170 716 120 28]);

    uilabel(ui.f,'Text','Surface / layer N','Position',[300 720 120 22]);
    ui.Nsurf = uieditfield(ui.f,'numeric','Limits',[1 max(1,maxRankAvail)], ...
        'RoundFractionalValues','on','Value',min(3,maxRankAvail), ...
        'Position',[420 716 60 28]);

    uilabel(ui.f,'Text','Radius r (nm)','Position',[500 720 90 22]);
    ui.R = uieditfield(ui.f,'numeric','Limits',[0.1 5],'Value',1.0,...
        'Position',[590 716 60 28]);

    ui.btn = uibutton(ui.f,'Text','Compute 1 nm enrichment',...
        'Position',[660 716 180 32], 'ButtonPushedFcn',@onComputeEnrichment);

    % 新增：中心层数 / 近邻层数
    uilabel(ui.f,'Text','Center layers (e.g. 2-3)',...
        'Position',[20 690 170 22]);
    ui.centerLayers = uieditfield(ui.f,'text','Value','1',...
        'Position',[190 686 80 28]);

    uilabel(ui.f,'Text','Neighbor layers (e.g. 2-5)',...
        'Position',[290 690 180 22]);
    ui.neighLayers = uieditfield(ui.f,'text','Value','1-3',...
        'Position',[470 686 80 28]);

    uilabel(ui.f,'Text','Center element','Position',[570 690 90 22]);
    avail = intersect(HEAset, string(categories(removecats(T.Species))),'stable');
    if isempty(avail)
        avail = string(categories(removecats(T.Species)));
    end
    ui.center = uidropdown(ui.f,'Items',cellstr(avail),...
        'Value',char(avail(1)),'Position',[660 686 80 28]);

    uilabel(ui.f,'Text','Neighbor mode','Position',[750 690 100 22]);
    ui.mode = uidropdown(ui.f,'Items',{'HEA (exclude center)',...
                                       'All species (exclude center)'}, ...
        'Value','HEA (exclude center)','Position',[850 686 200 28]);

    uilabel(ui.f,'Text','Baseline','Position',[20 660 80 22]);
    ui.base = uidropdown(ui.f,'Items',{'Global baseline','Surface baseline'},...
        'Value','Global baseline','Position',[90 656 150 28]);

    ui.btnExportCenters = uibutton(ui.f,'Text','Export centers of selection',...
        'Position',[260 656 200 28],'ButtonPushedFcn',@onExportCenters);

    uilabel(ui.f,'Text','Enrichment within r (nm)','FontWeight','bold',...
        'Position',[20 630 220 20]);
    ui.tabRes = uitable(ui.f,'Position',[20 418 1140 210], ...
        'ColumnName',{'Neighbor','PairsWithinR','P_withinR',...
                      'pi_baseline','Enrichment E','Tendency'}, ...
        'ColumnEditable',false(1,6));

    uilabel(ui.f,'Text','Per-layer Z stats (HEA only)','FontWeight','bold',...
        'Position',[20 388 220 20]);
    ui.tabZ = uitable(ui.f,'Position',[20 258 560 120], ...
        'ColumnName',{'Layer','Count','MeanZ','MedianZ','StdZ','MinZ','MaxZ'}, ...
        'ColumnEditable',false(1,7));
    Kshow = min(10, maxRankAvail);
    [Layer, Count, MeanZ, MedianZ, StdZ, MinZ, MaxZ] = ...
        layerStats(T, cellOrderByD, HEAset, Kshow);
    ui.tabZ.Data = table(Layer, Count, MeanZ, MedianZ, StdZ, MinZ, MaxZ);

    ui.info = uitextarea(ui.f,'Position',[600 258 560 120],'Editable','off');

    uilabel(ui.f,'Text','Distance-based analysis (KDE / Posterior / Enrichment)',...
        'FontWeight','bold','Position',[20 228 400 20]);
    ui.btnPlotKDE = uibutton(ui.f,'Text','Plot KDE/Posterior/Enrichment',...
        'Position',[420 224 230 28], 'ButtonPushedFcn',@onPlotKDEPosterior);

    % ----------------------- 底部：kNN & 1NN 工具 ---------------------
    p = uipanel(ui.f,'Title','kNN Frequency (empirical vs PPP baseline) & 1NN tools',...
        'Position',[20 10 1140 210]);

    uilabel(p,'Text','Mode','Position',[10 176 40 20]);
    ui.knnMode = uidropdown(p,'Items',{'Geometric (pooled)',...
                                       'Center-conditioned by neighbor species'}, ...
        'Value','Geometric (pooled)','Position',[60 172 230 28]);

    uilabel(p,'Text','Use current surface selection as ROI',...
        'Position',[310 176 220 20]);
    ui.useSurfROI = uicheckbox(p,'Value',true,'Position',[525 176 20 20]);

    uilabel(p,'Text','Neighbor universe','Position',[10 148 120 20]);
    ui.knnUni = uidropdown(p,'Items',{'HEA only','All species'},...
        'Value','HEA only','Position',[130 144 140 28]);

    uilabel(p,'Text','Centers','Position',[290 148 60 20]);
    ui.centerMode = uidropdown(p,'Items',{'All in ROI','Specific species'},...
        'Value','All in ROI','Position',[350 144 130 28]);
    ui.centerSp   = uidropdown(p,'Items',...
        cellstr(string(categories(removecats(T.Species)))), ...
        'Value',char(avail(1)), 'Position',[490 144 120 28],'Enable','off');

    ui.knnMode.ValueChangedFcn    = @onKnnModeChanged;
    ui.centerMode.ValueChangedFcn = @onCenterModeChanged;

    uilabel(p,'Text','k (nearest neighbor order)','Position',[620 148 170 20]);
    ui.kOrder = uieditfield(p,'numeric','Limits',[1 50],'Value',1,...
        'RoundFractionalValues','on','Position',[790 144 60 28]);

    uilabel(p,'Text','Max centers','Position',[870 148 90 20]);
    ui.knnMaxN = uieditfield(p,'numeric','Limits',[500 Inf],'Value',2e5,...
        'Position',[960 144 100 28]);

    uilabel(p,'Text','Bins','Position',[10 116 40 20]);
    ui.binMode = uidropdown(p,'Items',{'Freedman–Diaconis','Fixed number'},...
        'Value','Freedman–Diaconis','Position',[60 112 170 28]);
    uilabel(p,'Text','nbins','Position',[240 116 50 20]);
    ui.nbins = uieditfield(p,'numeric','Limits',[10 200],'Value',60,...
        'Position',[290 112 70 28]);

    ui.btnKNNfreq = uibutton(p,'Text','Compute kNN Frequency',...
        'Position',[380 112 150 28],'ButtonPushedFcn',@onComputeKNNfreq);
    ui.btnExportKNN = uibutton(p,'Text','Export to Excel',...
        'Position',[540 112 120 28],'ButtonPushedFcn',@onExportKNN);

    ui.btnPlot4 = uibutton(p,'Text','Plot 1NN^s (4 figs)',...
        'Position',[670 112 150 28],'ButtonPushedFcn',@onPlotNN1PerSpecies);

    ui.btnNN1MultiUnion = uibutton(p,'Text','1NN multi vs union',...
        'Position',[830 112 180 28],'ButtonPushedFcn',@onPlotNN1MultiVsUnion);

    ui.axKNN = uiaxes(p,'Position',[10 10 1120 90]);
    title(ui.axKNN,''); xlabel(ui.axKNN,''); ylabel(ui.axKNN,'');

    setappdata(ui.f,'UISTRUCT',ui);

    % ------------ 内部小工具：Build surface mask ------------
    function [surfaceMask, baseName] = buildSurfaceMask(ui_local)
        Ns = max(1, min(maxRankAvail, round(ui_local.Nsurf.Value)));
        switch ui_local.layerMode.Value
            case 'Only layer N'
                surfaceMask = mask_from_rank(cellOrderByD, Ns, Ntot);
                baseName = sprintf('SurfOnlyL%d', Ns);
            otherwise
                surfaceMask = false(Ntot,1);
                for kk=1:Ns
                    surfaceMask = surfaceMask | mask_from_rank(cellOrderByD, kk, Ntot);
                end
                baseName = sprintf('Surf1to%d', Ns);
        end
    end

    % 把 "2-3,5" 解析成层号
    function layers = parseLayerString(str, Kmax)
        s = strtrim(str);
        if isempty(s)
            layers = 1:Kmax; return;
        end
        s = strrep(s,',',' ');
        s = strrep(s,';',' ');
        parts = strsplit(s);
        layers = [];
        for ii = 1:numel(parts)
            tok = strtrim(parts{ii});
            if isempty(tok), continue; end
            if contains(tok,'-')
                ab = strsplit(tok,'-');
                if numel(ab) ~= 2, continue; end
                a = str2double(ab{1}); b = str2double(ab{2});
                if isnan(a) || isnan(b), continue; end
                lo = min(a,b); hi = max(a,b);
                layers = [layers, lo:hi]; %#ok<AGROW>
            else
                v = str2double(tok);
                if ~isnan(v), layers = [layers, v]; end %#ok<AGROW>
            end
        end
        layers = unique(layers);
        layers = layers(layers>=1 & layers<=Kmax);
        if isempty(layers), layers = 1:Kmax; end
    end

    function [centerLayerMask, neighborLayerMask, cLayers, nLayers] = getLayerMasks()
        cLayers = parseLayerString(ui.centerLayers.Value, maxRankAvail);
        nLayers = parseLayerString(ui.neighLayers.Value, maxRankAvail);

        centerLayerMask   = false(Ntot,1);
        neighborLayerMask = false(Ntot,1);
        for kk = cLayers
            if kk<=numel(L), centerLayerMask   = centerLayerMask   | L{kk}; end
        end
        for kk = nLayers
            if kk<=numel(L), neighborLayerMask = neighborLayerMask | L{kk}; end
        end
    end

    % ------------------ 富集 & KDE 回调 -------------------
    function onComputeEnrichment(src,~)
        f  = ancestor(src,'figure'); ui = getappdata(f,'UISTRUCT');

        [surfaceMask, baseName] = buildSurfaceMask(ui);
        [centerLayerMask, neighborLayerMask, cL, nL] = getLayerMasks();

        centerA = string(ui.center.Value);
        centerMask = surfaceMask & centerLayerMask & (string(T.Species) == centerA);

        modeStr = ui.mode.Value;
        baseStr = ui.base.Value;

        if contains(modeStr,'HEA')
            neighList = setdiff(HEAset, centerA, 'stable');
            neighborsUniverseMask = ismember(string(T.Species), neighList);
        else
            allSp = string(categories(removecats(T.Species)));
            neighList = setdiff(allSp, centerA, 'stable');
            neighborsUniverseMask = ismember(string(T.Species), neighList);
        end

        neighborsMask = surfaceMask & neighborLayerMask & neighborsUniverseMask;

        R = ui.R.Value;
        [ResTbl, metaTxt] = enrichment_within_radius(T, centerMask, neighborsMask,...
                            neighList, R, surfaceMask, baseStr);

        metaTxt = sprintf('%s\nCenter layers: %s | Neighbor layers: %s',...
            metaTxt, mat2str(cL), mat2str(nL));

        ui.tabRes.Data = ResTbl;
        ui.info.Value = sprintf("%s\nSelection mode: %s", ...
            metaTxt, ui.layerMode.Value);

        try
            writetable(ResTbl, outXLSX, 'Sheet', ...
                sprintf('%s_Enrich_%s', baseName, char(centerA)));
        catch
        end

        onPlotKDEPosterior(ui.btnPlotKDE,[]);
        setappdata(f,'UISTRUCT',ui);
    end

    function onPlotKDEPosterior(src,~)
        f  = ancestor(src,'figure'); ui = getappdata(f,'UISTRUCT');

        [surfaceMask, ~] = buildSurfaceMask(ui);
        [centerLayerMask, neighborLayerMask, ~, ~] = getLayerMasks();

        centerA = string(ui.center.Value);
        centerMask = surfaceMask & centerLayerMask & (string(T.Species) == centerA);

        modeStr = ui.mode.Value;
        baseStr = ui.base.Value;

        if contains(modeStr,'HEA')
            neighList = setdiff(HEAset, centerA, 'stable');
            neighborsUniverseMask = ismember(string(T.Species), neighList);
        else
            allSp = string(categories(removecats(T.Species)));
            neighList = setdiff(allSp, centerA, 'stable');
            neighborsUniverseMask = ismember(string(T.Species), neighList);
        end
        neighborsMask = surfaceMask & neighborLayerMask & neighborsUniverseMask;

        if ~any(centerMask) || ~any(neighborsMask), return; end

        r_min = 0.22; r_max = 1.00; nGrid = 400; bw = 0.03;
        [tblKDE, ~] = kde_pdf_pairs(T, centerMask, neighborsMask, ...
                                    neighList, r_min, r_max, bw, nGrid);

        r  = unique(tblKDE.r); nR = numel(r); M  = numel(neighList);
        F  = zeros(nR, M);
        for kk = 1:M
            rows = (tblKDE.NeighborSpecies == neighList(kk));
            F(:,kk) = tblKDE.pdf_kde(rows);
        end

        if strcmp(baseStr,'Surface baseline')
            baseMask = neighborsMask & surfaceMask;
        else
            baseMask = neighborsMask;
        end
        spBase = string(T.Species(baseMask));
        spBase = spBase(ismember(spBase, neighList));
        counts = arrayfun(@(s) sum(spBase == s), neighList);
        if sum(counts)==0
            pi_s = ones(M,1)/M;
        else
            pi_s = counts(:)/sum(counts);
        end

        g = F * pi_s; g(g<=0) = eps;
        P = (F .* pi_s.') ./ g;
        E = P ./ pi_s.';

        r_plot = [0; r];
        F_plot = [zeros(1,M); F];
        P_plot = [zeros(1,M); P];
        E_plot = [zeros(1,M); E];

        [rFCC, ~] = fcc_shell_params(0.25, 5, 0.04);

        % PDF
        figure('Color','w'); hold on;
        for ii = 1:numel(rFCC)
            xline(rFCC(ii),'--','Color',[.75 .75 .75]);
        end
        xline(r_min,':','Color',[.6 .6 .6],'HandleVisibility','off');
        for kk = 1:M
            plot(r_plot, F_plot(:,kk),'LineWidth',1.8,...
                'DisplayName',char(neighList(kk)));
        end
        grid on; box on;
        xlabel('r (nm)'); ylabel('Probability density (1/nm)');
        title(sprintf('Pair-distance PDF (KDE) — Center: %s', centerA));
        legend('Location','northeastoutside'); xlim([0 1.0]);

        % Posterior
        figure('Color','w'); hold on;
        for ii = 1:numel(rFCC)
            xline(rFCC(ii),'--','Color',[.75 .75 .75]);
        end
        xline(r_min,':','Color',[.6 .6 .6],'HandleVisibility','off');
        for kk = 1:M
            plot(r_plot, P_plot(:,kk),'LineWidth',1.8,...
                'DisplayName',sprintf('%s (P)',char(neighList(kk))));
        end
        grid on; box on;
        xlabel('r (nm)'); ylabel('P(neighbor = s | r)'); ylim([0 1]);
        title(sprintf('Posterior around %s (%s baseline)', centerA, lower(baseStr)));
        legend('Location','northeastoutside'); xlim([0 1.0]);

        % Enrichment
        figure('Color','w'); hold on; yline(1,'k:');
        for ii = 1:numel(rFCC)
            xline(rFCC(ii),'--','Color',[.75 .75 .75]);
        end
        xline(r_min,':','Color',[.6 .6 .6],'HandleVisibility','off');
        for kk = 1:M
            plot(r_plot, E_plot(:,kk),'LineWidth',1.8,...
                'DisplayName',sprintf('%s (E)',char(neighList(kk))));
        end
        grid on; box on;
        xlabel('r (nm)'); ylabel(sprintf('Enrichment E_s(r) (%s baseline)',lower(baseStr)));
        title(sprintf('Enrichment vs %s baseline — Center: %s',lower(baseStr),centerA));
        legend('Location','northeastoutside'); xlim([0 1.0]);

        setappdata(f,'UISTRUCT',ui);
    end

    % ------------------- kNN 部分 -----------------------
    function onKnnModeChanged(~,~)
        f  = ancestor(ui.f,'figure'); ui2 = getappdata(f,'UISTRUCT'); %#ok<NASGU>
        ui.centerMode.Value = 'Specific species';
        ui.centerSp.Enable  = 'on';
        setappdata(f,'UISTRUCT',ui);
    end

    function onCenterModeChanged(~,~)
        f  = ancestor(ui.f,'figure'); ui2 = getappdata(f,'UISTRUCT'); %#ok<NASGU>
        ui.centerSp.Enable = ternary(strcmp(ui.centerMode.Value,'Specific species'),'on','off');
        setappdata(f,'UISTRUCT',ui);
    end

    function onComputeKNNfreq(src,~)
        f  = ancestor(src,'figure'); ui = getappdata(f,'UISTRUCT');

        [surfaceMask, baseName] = buildSurfaceMask(ui);
        roi0 = ternary(ui.useSurfROI.Value, surfaceMask, true(Ntot,1));

        [centerLayerMask, neighborLayerMask, cL, nL] = getLayerMasks();
        roiMask = roi0 & neighborLayerMask;

        if strcmp(ui.knnUni.Value,'HEA only')
            poolMaskBase = roiMask & ismember(string(T.Species), HEAset);
            neighListAll = setdiff(intersect(HEAset, string(categories(removecats(T.Species)))),"",'stable');
        else
            poolMaskBase = roiMask;
            neighListAll = setdiff(string(categories(removecats(T.Species))),"",'stable');
        end

        binSpec = struct();
        if strcmp(ui.binMode.Value,'Freedman–Diaconis')
            binSpec.mode = 'fd';
        else
            binSpec.mode = 'nbins';
            binSpec.nbins = max(10, round(ui.nbins.Value));
        end
        k = max(1, round(ui.kOrder.Value));
        maxCenters = max(500, round(ui.knnMaxN.Value));

        colHex = struct('Ir','#FF0000','Pt','#0000FF','Ru','#FF00FF',...
                        'Pd','#00CCFF','Rh','#000000');

        switch ui.knnMode.Value
            case 'Geometric (pooled)'
                XYZ_pool = [double(T.X(poolMaskBase)),...
                            double(T.Y(poolMaskBase)),...
                            double(T.Z(poolMaskBase))];

                if strcmp(ui.centerMode.Value,'Specific species')
                    sp = string(ui.centerSp.Value);
                    centersMask = roiMask & centerLayerMask & (string(T.Species) == sp);
                    centerTag = char(sp);
                else
                    centersMask = roiMask & centerLayerMask;
                    centerTag = 'All';
                end
                XYZ_centers = [double(T.X(centersMask)),...
                               double(T.Y(centersMask)),...
                               double(T.Z(centersMask))];
                if size(XYZ_centers,1) > maxCenters
                    XYZ_centers = XYZ_centers(randperm(size(XYZ_centers,1), maxCenters),:);
                end

                try
                    [cens, empCounts, baseCounts, edges, meta] = ...
                        knn_frequency_with_ppp(XYZ_pool, XYZ_centers, k, binSpec);
                catch ME
                    uialert(ui.f, sprintf('kNN frequency failed:\n%s', ME.message), 'kNN Frequency'); 
                    return;
                end

                cla(ui.axKNN); hold(ui.axKNN,'on');
                stairs(ui.axKNN, cens, empCounts, 'LineWidth',1.8, 'DisplayName','Empirical');
                stairs(ui.axKNN, cens, baseCounts, 'LineWidth',1.8,'LineStyle','--',...
                       'DisplayName','PPP baseline');
                grid(ui.axKNN,'on'); box(ui.axKNN,'on');
                xlabel(ui.axKNN, sprintf('R_{%d} (nm)', k)); ylabel(ui.axKNN, 'Counts');
                title(ui.axKNN, sprintf('kNN distance frequency — ROI=%s | Centers=%s | Universe=%s',...
                    ternary(ui.useSurfROI.Value,'Surface sel.','All points'), centerTag, ui.knnUni.Value));
                legend(ui.axKNN,'Location','northeast'); hold(ui.axKNN,'off');

                setappdata(f,'KNN_LAST', struct('mode','pooled','centers',cens,...
                    'empCounts',empCounts,'baseCounts',baseCounts,'edges',edges,...
                    'k',k,'centerTag',centerTag,'universe',ui.knnUni.Value,...
                    'useSurf',ui.useSurfROI.Value,'meta',meta,'baseName',baseName,...
                    'centerLayers',cL,'neighLayers',nL));

            otherwise   % Center-conditioned by neighbor species
                spC = string(ui.centerSp.Value);
                centersMaskAll = (string(T.Species) == spC) & centerLayerMask;
                centersMask = centersMaskAll & roiMask;   % ROI + 层 + 元素

                if ~any(centersMask)
                    uialert(ui.f, sprintf('No %s centers in current ROI + layers.', char(spC)),...
                        'kNN Frequency'); 
                    return;
                end
                neighList = setdiff(neighListAll, spC, 'stable');

                idxC = find(centersMask);
                if numel(idxC) > maxCenters
                    idxC = idxC(randperm(numel(idxC), maxCenters));
                end
                centersMaskLimited = false(Ntot,1);
                centersMaskLimited(idxC) = true;

                try
                    CS = crossSpeciesKNNfreq_withPPP(T, centersMaskLimited, neighList,...
                            roiMask, k, binSpec);
                catch ME
                    uialert(ui.f, sprintf('kNN (by species) failed:\n%s', ME.message),...
                        'kNN Frequency'); 
                    return;
                end

                cla(ui.axKNN); hold(ui.axKNN,'on');
                for jj = 1:numel(CS.species)
                    s = CS.species(jj);
                    c = [0 0 0];
                    if isfield(colHex, char(s)), c = hex2rgb(colHex.(char(s))); end
                    stairs(ui.axKNN, CS.centers, CS.empCounts(:,jj),...
                        'LineWidth',1.8,'Color',c,'DisplayName',char(s));
                    stairs(ui.axKNN, CS.centers, CS.baseCounts(:,jj),...
                        'LineWidth',1.8,'Color',c,'LineStyle','--','HandleVisibility','off');
                end
                grid(ui.axKNN,'on'); box(ui.axKNN,'on');
                xlabel(ui.axKNN, sprintf('R_{%d}^{(s)} (nm)', k));
                ylabel(ui.axKNN, 'Counts');
                title(ui.axKNN, sprintf('Nearest-%d distance by neighbor species — Center=%s | ROI=%s',...
                    k, char(spC), ternary(ui.useSurfROI.Value,'Surface sel.','All')));
                legend(ui.axKNN,'Location','northeastoutside'); hold(ui.axKNN,'off');

                setappdata(f,'KNN_LAST', struct('mode','bySpecies','centers',CS.centers,...
                    'empCounts',CS.empCounts,'baseCounts',CS.baseCounts,...
                    'edges',CS.edges,'k',k,'centerTag',char(spC),'species',CS.species,...
                    'universe',ui.knnUni.Value,'useSurf',ui.useSurfROI.Value,...
                    'meta',CS.meta,'baseName',baseName,...
                    'centerLayers',cL,'neighLayers',nL));
        end
        setappdata(f,'UISTRUCT',ui);
    end

    function onExportKNN(src,~)
        f = ancestor(src,'figure');
        S = getappdata(f,'KNN_LAST');
        if isempty(S)
            uialert(f,'No kNN frequency results to export.','Export kNN'); 
            return;
        end
        switch S.mode
            case 'pooled'
                Tbl = table(S.centers(:), S.empCounts(:), S.baseCounts(:), ...
                    'VariableNames', {'Rk_center_nm','Empirical_count','PPP_baseline_count'});
                sh = sprintf('KNNfreq_%s_%s_%s_k%d', S.baseName, S.centerTag,...
                    ternary(S.useSurf,'Surf','All'), S.k);
            otherwise
                Tbl = table(S.centers(:), 'VariableNames', {'Rk_center_nm'});
                for jj = 1:numel(S.species)
                    sj = char(S.species(jj));
                    Tbl.(sprintf('Emp_%s_cnt', sj)) = S.empCounts(:,jj);
                    Tbl.(sprintf('PPP_%s_cnt', sj)) = S.baseCounts(:,jj);
                end
                sh = sprintf('KNNfreq_bySp_%s_%s_%s_k%d', S.baseName, S.centerTag,...
                    ternary(S.useSurf,'Surf','All'), S.k);
        end
        sh = regexprep(sh,'[^A-Za-z0-9_]','_');
        try
            writetable(Tbl, outXLSX, 'Sheet', sh);
            metaT = struct2table(S.meta,'AsArray',true);
            writetable(metaT, outXLSX, 'Sheet', sprintf('%s_meta', sh));
            uialert(f, sprintf('Exported to sheet:\n%s', sh), 'Export kNN');
        catch ME
            uialert(f, sprintf('Export failed:\n%s', ME.message), 'Export kNN');
        end
    end

    function onExportCenters(src,~)
        f  = ancestor(src,'figure'); ui = getappdata(f,'UISTRUCT');
        [surfaceMask, baseName] = buildSurfaceMask(ui);
        [centerLayerMask, ~, cL, ~] = getLayerMasks();
        centerA = string(ui.center.Value);
        centerMask = surfaceMask & centerLayerMask & (string(T.Species) == centerA);
        if ~any(centerMask)
            uialert(ui.f,'No centers under current selection.','Export centers'); 
            return;
        end
        C = T(centerMask, {'X','Y','Z','Species'});
        try
            writetable(C, outXLSX, 'Sheet', ...
                sprintf('%s_Centers_%s_L%s', baseName, char(centerA), mat2str(cL)));
            uialert(ui.f, sprintf('Exported %d centers to sheet.', height(C)), 'Export centers');
        catch ME
            uialert(ui.f, sprintf('Failed to write centers:\n%s', ME.message), 'Export centers');
        end
        setappdata(f,'UISTRUCT',ui);
    end

    % -------------------- 1NN 相关图 ---------------------
    function onPlotNN1PerSpecies(src,~)
        f  = ancestor(src,'figure');
        ui = getappdata(f,'UISTRUCT');

        [surfaceMask, baseName] = buildSurfaceMask(ui);
        roi0 = ternary(ui.useSurfROI.Value, surfaceMask, true(Ntot,1));
        [centerLayerMask, neighborLayerMask, cL, ~] = getLayerMasks();
        roiMask = roi0 & neighborLayerMask;

        if strcmp(ui.centerMode.Value,'Specific species')
            spC = string(ui.centerSp.Value);
        else
            spC = string(ui.center.Value);
        end

        centersMask = roiMask & centerLayerMask & (string(T.Species) == spC);
        idxC = find(centersMask);
        if isempty(idxC)
            uialert(ui.f, sprintf('No %s centers in current ROI + layers.', char(spC)),...
                '1NN^s plots');
            return;
        end

        maxCenters = max(500, round(ui.knnMaxN.Value));
        if numel(idxC) > maxCenters
            idxC = idxC(randperm(numel(idxC), maxCenters));
        end
        centersMaskLimited = false(Ntot,1);
        centersMaskLimited(idxC) = true;

        neighList = setdiff(HEAset, spC, 'stable');
        if isempty(neighList)
            uialert(ui.f,'No HEA neighbors available.','1NN^s plots');
            return;
        end

        binSpec = struct();
        if strcmp(ui.binMode.Value,'Freedman–Diaconis')
            binSpec.mode = 'fd';
        else
            binSpec.mode = 'nbins';
            binSpec.nbins = max(10, round(ui.nbins.Value));
        end

        try
            CS = crossSpeciesKNNfreq_withPPP(T, centersMaskLimited, neighList,...
                    roiMask, 1, binSpec);   % k = 1
        catch ME
            uialert(ui.f, sprintf('Failed to compute 1NN^s:\n%s', ME.message),...
                '1NN^s plots');
            return;
        end

        colHex = struct('Ir','#FF0000','Pt','#0000FF','Ru','#FF00FF',...
                        'Pd','#00CCFF','Rh','#000000');

        for j = 1:numel(CS.species)
            s = CS.species(j);
            if any(strcmpi(s, HEAset)) && ~strcmpi(s, spC)
                c = [0 0 0];
                if isfield(colHex, char(s)), c = hex2rgb(colHex.(char(s))); end
                figure('Color','w'); hold on;
                stairs(CS.centers, CS.empCounts(:,j),'LineWidth',1.8,...
                    'Color',c,'DisplayName','Empirical');
                stairs(CS.centers, CS.baseCounts(:,j),'LineWidth',1.8,...
                    'Color',c,'LineStyle','--','DisplayName','PPP baseline');
                grid on; box on;
                xlabel('R_1^{(s)} (nm)'); ylabel('Counts');
                title(sprintf('1NN to %s around %s — ROI=%s | Universe=HEA | C layers=%s', ...
                    char(s), char(spC), ternary(ui.useSurfROI.Value,'Surface sel.','All'), ...
                    mat2str(cL)));
                legend('Location','northeast'); hold off;
            end
        end

        try
            Tbl = table(CS.centers(:), 'VariableNames', {'R1_center_nm'});
            for j = 1:numel(CS.species)
                sj = char(CS.species(j));
                Tbl.(sprintf('Emp_%s_cnt', sj)) = CS.empCounts(:,j);
                Tbl.(sprintf('PPP_%s_cnt', sj)) = CS.baseCounts(:,j);
            end
            sh = sprintf('R1_bySp_%s_%s_%s', baseName, char(spC), ...
                ternary(ui.useSurfROI.Value,'Surf','All'));
            sh = regexprep(sh,'[^A-Za-z0-9_]','_');
            writetable(Tbl, outXLSX, 'Sheet', sh);
        catch
        end

        setappdata(f,'UISTRUCT',ui);
    end

    function onPlotNN1MultiVsUnion(src,~)
        f  = ancestor(src,'figure');
        ui = getappdata(f,'UISTRUCT');

        [surfaceMask, baseName] = buildSurfaceMask(ui);
        roi0 = ternary(ui.useSurfROI.Value, surfaceMask, true(Ntot,1));
        [centerLayerMask, neighborLayerMask, cL, nL] = getLayerMasks();
        roiMask = roi0 & neighborLayerMask;

        if strcmp(ui.centerMode.Value,'Specific species')
            spC = string(ui.centerSp.Value);
        else
            spC = string(ui.center.Value);
        end

        centersMask = roiMask & centerLayerMask & (string(T.Species) == spC);
        idxC = find(centersMask);
        if isempty(idxC)
            uialert(ui.f, sprintf('No %s centers in current ROI + layers.', char(spC)),...
                '1NN multi vs union');
            return;
        end
        maxCenters = max(500, round(ui.knnMaxN.Value));
        if numel(idxC) > maxCenters
            idxC = idxC(randperm(numel(idxC), maxCenters));
        end
        centersMaskLimited = false(Ntot,1);
        centersMaskLimited(idxC) = true;

        neighList = setdiff(HEAset, spC, 'stable');
        neighList = neighList(ismember(neighList,...
            string(categories(removecats(T.Species)))));
        if isempty(neighList)
            uialert(ui.f,'No HEA neighbors available in current ROI.','1NN multi vs union');
            return;
        end

        binMode = ternary(strcmp(ui.binMode.Value,'Freedman–Diaconis'),'fd','fixed');
        nbins   = max(10, round(ui.nbins.Value));

        Q = [double(T.X(centersMaskLimited)),...
             double(T.Y(centersMaskLimited)),...
             double(T.Z(centersMaskLimited))];
        Ncent = size(Q,1);

        nbrUnionMask = roiMask & ismember(string(T.Species), neighList);
        idxU = find(nbrUnionMask);
        if numel(idxU) < 2
            uialert(ui.f,'Too few neighbours in ROI.','1NN multi vs union');
            return;
        end
        P_union = [double(T.X(idxU)), double(T.Y(idxU)), double(T.Z(idxU))];
        M_union = createns(P_union,'NSMethod','kdtree');
        [~, D_all] = knnsearch(M_union, Q, 'K', 1);
        R1_all = D_all(:,1);

        S = numel(neighList);
        R1_s = cell(S,1);
        allR = R1_all;
        for j = 1:S
            s = neighList(j);
            maskS = roiMask & (string(T.Species) == s);
            idxS  = find(maskS);
            if numel(idxS) < 2
                R1_s{j} = NaN(Ncent,1);
                continue;
            end
            P_s = [double(T.X(idxS)), double(T.Y(idxS)), double(T.Z(idxS))];
            M_s = createns(P_s,'NSMethod','kdtree');
            [~, D_s] = knnsearch(M_s, Q, 'K', 1);
            R1_s{j} = D_s(:,1);
            allR = [allR; D_s(:,1)]; %#ok<AGROW>
        end

        allR = allR(isfinite(allR));
        if numel(allR) < 10
            uialert(ui.f,'Too few valid R1 samples.','1NN multi vs union');
            return;
        end

        switch lower(binMode)
            case 'fd'
                iq = iqr(allR);
                h  = 2*iq / max(numel(allR)^(1/3),1);
                if ~isfinite(h) || h <= 0
                    nb = max(10, round(sqrt(numel(allR))));
                else
                    nb = max(10, ceil((max(allR)-min(allR))/h));
                end
            otherwise
                nb = max(10, nbins);
        end
        edges   = linspace(min(allR), max(allR), nb+1);
        centersR1 = (edges(1:end-1) + edges(2:end))/2;

        H_union = histcounts(R1_all(isfinite(R1_all)), edges);
        H_s     = zeros(nb, S);
        for j = 1:S
            rj = R1_s{j};
            rj = rj(isfinite(rj));
            if isempty(rj), continue; end
            H_s(:,j) = histcounts(rj, edges);
        end

        colHex = struct('Ir','#FF0000','Pt','#0000FF','Ru','#FF00FF',...
                        'Pd','#00CCFF','Rh','#000000');

        figure('Color','w'); hold on;
        for j = 1:S
            s = neighList(j);
            c = [0 0 0];
            if isfield(colHex, char(s)), c = hex2rgb(colHex.(char(s))); end
            stairs(centersR1, H_s(:,j),'LineWidth',1.8,...
                   'Color',c,'DisplayName',char(s));
        end
        stairs(centersR1, H_union,'k--','LineWidth',2.0,...
               'DisplayName','All neighbours');

        grid on; box on;
        xlabel(sprintf('R_1 (nm) from %s', spC));
        ylabel('Counts');
        title(sprintf('1NN around %s in ROI (%d centers) | C layers=%s | N layers=%s',...
            spC, Ncent, mat2str(cL), mat2str(nL)));
        legend('Location','northeastoutside');
        hold off;

        try
            Tbl = table(centersR1(:), 'VariableNames', {'R1_center_nm'});
            for j = 1:S
                sj = char(neighList(j));
                Tbl.(sprintf('Emp_%s_cnt', sj)) = H_s(:,j);
            end
            Tbl.All_neigh_cnt = H_union(:);
            sh = sprintf('R1_multi_vs_union_%s_%s_%s', baseName, char(spC),...
                ternary(ui.useSurfROI.Value,'Surf','All'));
            sh = regexprep(sh,'[^A-Za-z0-9_]','_');
            writetable(Tbl, outXLSX, 'Sheet', sh);
        catch
        end

        setappdata(f,'UISTRUCT',ui);
    end
end

function [Layer, Count, MeanZ, MedianZ, StdZ, MinZ, MaxZ] = layerStats(T, cellOrderByD, HEAset, Kshow)
    N = height(T);
    Layer = (1:Kshow)'; Count=zeros(Kshow,1); MeanZ=nan(Kshow,1); MedianZ=nan(Kshow,1);
    StdZ=nan(Kshow,1); MinZ=nan(Kshow,1); MaxZ=nan(Kshow,1);
    for k=1:Kshow
        Mk = mask_from_rank(cellOrderByD,k,N);
        Tk = T(Mk & ismember(string(T.Species),HEAset),:);
        Count(k)=height(Tk);
        if Count(k)>0
            MeanZ(k)=mean(Tk.Z,'omitnan'); MedianZ(k)=median(Tk.Z,'omitnan');
            StdZ(k)=std(Tk.Z,0,'omitnan');  MinZ(k)=min(Tk.Z);  MaxZ(k)=max(Tk.Z);
        end
    end
end

function [ResTbl, metaTxt] = enrichment_within_radius(T, centerMask, neighborsMask, neighList, R, surfaceMask, baseStr)
    neighList = string(neighList(:)');
    centers = find(centerMask);
    ResTbl = table();
    if isempty(centers)
        metaTxt = "No centers."; return;
    end

    XYZ = [T.X, T.Y, T.Z];
    nbrIdx = find(neighborsMask);
    if isempty(nbrIdx)
        metaTxt = "No neighbors under current mode."; return;
    end

    Mdl = createns(XYZ(nbrIdx,:), 'NSMethod','kdtree');
    K = numel(neighList);
    counts = zeros(1,K);
    idxCell = rangesearch(Mdl, XYZ(centers,:), R);
    for ii=1:numel(centers)
        c = centers(ii);
        glb = nbrIdx(idxCell{ii});
        glb(glb==c) = [];
        if isempty(glb), continue; end
        sp = string(T.Species(glb));
        for k=1:K
            counts(k) = counts(k) + sum(sp == neighList(k));
        end
    end
    totInR = sum(counts);
    if totInR==0
        metaTxt = sprintf('No neighbors found within %.3f nm.', R);
        ResTbl = table(neighList', zeros(K,1), zeros(K,1), zeros(K,1), zeros(K,1), strings(K,1), ...
            'VariableNames', {'Neighbor','PairsWithinR','P_withinR','pi_baseline','Enrichment','Tendency'});
        return;
    end
    P_within = counts / totInR;

    if strcmp(baseStr,'Surface baseline')
        baseMask = neighborsMask & surfaceMask;
    else
        baseMask = neighborsMask;
    end
    spBase = string(T.Species(baseMask));
    spBase = spBase(ismember(spBase, neighList));
    baseCounts = arrayfun(@(s) sum(spBase == s), neighList);
    if sum(baseCounts)==0
        pi = ones(1,K)/K;
    else
        pi = baseCounts / sum(baseCounts);
    end

    E = P_within ./ (pi + eps);

    tol = 0.05; tend = strings(K,1);
    for k=1:K
        if E(k) > 1+tol
            tend(k)="Aggregate";
        elseif E(k) < 1-tol
            tend(k)="Repel";
        else
            tend(k)="Neutral";
        end
    end

    ResTbl = table(neighList', counts', P_within', pi', E', tend, ...
        'VariableNames', {'Neighbor','PairsWithinR','P_withinR','pi_baseline','Enrichment','Tendency'});

    metaTxt = sprintf(['Centers: %d (in selected surface)\n' ...
                       'Radius: %.3f nm\n' ...
                       'Neighbors universe: %s\n' ...
                       'Baseline: %s\n' ...
                       'Total neighbor pairs within R: %d'], ...
                       numel(centers), R, strjoin(neighList,','), baseStr, totInR);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
