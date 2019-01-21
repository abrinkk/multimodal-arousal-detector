%% Check sampling frequencies
clear all; close all;
dirIndex = paths;
data_paths;
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

%% WSC
wsc_sigs = {};
wsc_fs = [];
for i = 1:length(f_wsc_edf)
    fprintf('Processsing WSC %.0f/%.0f\n',i,length(f_wsc_edf));
    try
        % Load and preprocess data
        hdr = loadHDR(filepath(p_wsc_edf,f_wsc_edf{i}));
        wsc_sigs = ([wsc_sigs; hdr.label]);
        wsc_fs = ([wsc_fs; hdr.fs]);
    catch me
        disp(me.message);
    end
end
for i = 1:length(f_wsc2_edf)
    fprintf('Processsing WSC2 %.0f/%.0f\n',i,length(f_wsc2_edf));
    try
        % Load and preprocess data
        hdr = loadHDR(filepath(p_wsc2_edf,f_wsc2_edf{i}));
        wsc_sigs = ([wsc_sigs; hdr.label]);
        wsc_fs = ([wsc_fs; hdr.fs]);
    catch me
        disp(me.message);
    end
end
uniqueSigs = unique(wsc_sigs);
dUS = 1:length(uniqueSigs);
dS = zeros(size(wsc_sigs));
dFs = cell(size(dUS));

for i = 1:length(dUS)
    dS(ismember(wsc_sigs,uniqueSigs{i})) = dUS(i);
    dFs{i} = wsc_fs(dS == i);
end

h = figure;
h.Position(3:4) = [600 600];
centerfig(h);
subplot(2,1,1);
boxplot2(dFs)
ylabel('Fs [Hz]')
set(gca,'XTick',1:length(dUS));
set(gca,'XTickLabel',uniqueSigs);
set(gca,'XTickLabelRotation',90);
xlabel('Channel');
box on;
grid minor;
subplot(2,1,2);
counts = hist(dS,dUS);
bar(dUS,counts./length(f_wsc_edf));
set(gca,'XTick',1:length(dUS));
set(gca,'XTickLabel',uniqueSigs);
set(gca,'XTickLabelRotation',90);
xlabel('Channel');
ylabel('Ratio of Availability');
xlim([0 length(dUS)+1]);
ylim([0 1]);
box on;
grid minor;
