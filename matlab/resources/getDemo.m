function [age, bmi, sex] = getDemo(files,ftype)
% Read table data
T_mros = readtable('E:\mros\datasets\mros-visit1-dataset-0.3.0.csv');
T_cfs = readtable('E:\cfs\datasets\cfs-visit5-dataset-0.4.0.csv');
T_wsc2 = load('E:\meta');
T_wsc2 = T_wsc2.tab;
T_ssc = load('E:\APOE_demographics');
T_ssc = T_ssc.t;
T_wsc = readtable('E:\WSC_PLM_ data_all.xlsx');
T_wsc.Record_ID = cellstr([char(T_wsc.SUBJ_ID) repmat('_',size(T_wsc,1),1) num2str(T_wsc.VISIT_NUMBER)]);
T_ssc.Record_ID = cellstr([repmat('SSC_',size(T_ssc,1),1) num2str(T_ssc.patid) repmat('_',size(T_ssc,1),1) num2str(T_ssc.visitsequence)]);

age = zeros(size(files));
bmi = zeros(size(files));
sex = zeros(size(files));
for i = 1:length(files)
    try
        switch ftype(i)
            case 1
                idx = (contains(lower(T_mros.nsrrid),files{i}(13:end-4)));
                age(i) = T_mros.vsage1(idx);
                bmi_t = str2num(cell2mat(T_mros.hwbmi(i)));
                if ~isempty(bmi_t)
                    bmi(i) = bmi_t;
                end
                sex(i) = 2;
            case 2
                idx = (((T_cfs.nsrrid) == str2num(files{i}(12:end-4))));
                age(i) = T_cfs.age(idx);
                bmi(i) = T_cfs.bmi(idx);
                sex(i) = T_cfs.SEX(idx) + 1;
            case 3
                idx = (contains((T_wsc2.Properties.RowNames),files{i}(1:end-4)));
                age(i) = T_wsc2.Age(idx);
                bmi(i) = T_wsc2.BMI(idx);
                sex(i) = (T_wsc2.Sex{idx} == 'M')*2 + (T_wsc2.Sex{idx} == 'F');
            case 4
                idx = (contains((T_wsc.Record_ID),files{i}(1:end-11)));
                age(i) = T_wsc.age(idx);
                bmi(i) = T_wsc.bmi(idx);
                sex(i) = (T_wsc.SEX{idx} == 'M')*2 + (T_wsc.SEX{idx} == 'F');
        end
    end
end
end

