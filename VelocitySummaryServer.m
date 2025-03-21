folder_data = '/home/guevara/RESULTS/ResIlariaPatSubjFilter_SummaryVelocity/';
folder_save = '/home/guevara/RESULTS/Temporary_Results/';
dir_folder_data = dir(folder_data);
C = cell(length(dir_folder_data)-2,1);
for i = 3:length(dir_folder_data)
    folder_freq = [folder_data,dir_folder_data(i).name];
    dir_folder_freq = dir([folder_data,dir_folder_data(i).name]);
    % Find folders of patients only:
    M = zeros(length(dir_folder_freq)-2,4); % M: col1 = patient or subject
    % col2 = central freq, col3 = mean v, col4 = std vel
    if strfind(dir_folder_data(i).name,'patient')
        M(:,1) = 1; % code: patient = 1; subject = 2
    else
        M(:,1) = 2; % code: patient = 1; subject = 2
    end
    for j = 3:length(dir_folder_freq)
        load([folder_freq,'/',dir_folder_freq(j).name, '/Velocity']);
        file_name = dir_folder_freq(j).name;
        I = strfind(file_name,'_');
        f1 = str2double(file_name(I(1)+1:I(2)-1));
        f2 = str2double(file_name(I(2)+1:end));
        fc = (f1+f2)/2; % central freq of the filter
        M(j-2,2) = fc;
        M(j-2,3) = v_corr_mean;
        M(j-2,4) = v_corr_std;
    end
    C{i-2} = M ;
end
save([folder_save,'velocity_summary_old_data_again'],"C");
