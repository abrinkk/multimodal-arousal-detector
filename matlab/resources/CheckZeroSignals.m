%% Check error in files

clear all; close all;
des_fs = 128;
dirIndex = paths;
f_train = dir(filepath(dirIndex.Data,'train','*.txt'));
f_test = dir(filepath(dirIndex.Data,'test','*.txt'));
f_val = dir(filepath(dirIndex.Data,'val','*.txt'));
f_other = dir(filepath(dirIndex.Data,'other','*.txt'));
f_train = {f_train.name};
f_test = {f_test.name};
f_val = {f_val.name};
f_other = {f_other.name};
f_all = [f_train f_test f_val f_other];

%% Check #Zeros in files

num_zero_all = zeros(size(f_all));
folder = [repelem({'train'},length(f_train)) ...
    repelem({'test'},length(f_test)) ...
    repelem({'val'},length(f_val)) ...
    repelem({'other'},length(f_other))];
time = 5;

for i = 1:length(f_all)
    tic;
    FID = fopen(filepath(dirIndex.Data,folder{i},f_all{i}),'r');
    C = textscan(FID,'%.2f','Delimiter',',');
    fclose(FID);
    num_zero_all(i) = sum(C{1} == 0);
    time_temp = toc;
    time = time*0.9 + time_temp*0.1;
    fprintf('File number: %.0f. Time remaining: %.1f min\n',i,(length(f_all) - i)*time/60);
end

%% Analyze file zero distribution
num_zero_all_sort = sort(num_zero_all);

figure
plot(num_zero_all_sort)

thresh = 1000000;

zero_files = char(f_all(num_zero_all> thresh) );