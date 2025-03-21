folder_data = '/home/guevara/RESULTS/Results_Ilaria_Patients_Subjects_filter/';
folder_save = '/home/guevara/RESULTS/ResIlariaPatSubjFilter_SummaryVelocity/';
dir_folder_data = dir(folder_data);
for i = 3:length(dir_folder_data)
    folder_patient_subj = [folder_data,dir_folder_data(i).name,'/' ];
    folder_patient_subj_save = [folder_save,dir_folder_data(i).name,'/' ];
    mkdir(folder_patient_subj_save);
    dir_folder_patient_subj = dir(folder_patient_subj);
    for j = 3:length(dir_folder_patient_subj)
        folder_freq = [folder_patient_subj,dir_folder_patient_subj(j).name,'/' ];
        folder_freq_save = [folder_patient_subj_save,dir_folder_patient_subj(j).name,'/' ];
        mkdir(folder_freq_save);
        dir_folder_freq = dir(folder_freq);
        V = [];
        D = [];
        for k = 3:length(dir_folder_freq)
            load([folder_freq dir_folder_freq(k).name]);
            V = [V;velocity];
            num_steps = zeros(length(d_pieces),1);
            for m = 1: length(d_pieces)
                num_steps(m) = sum(d_pieces{m}>0);
            end
            D = [D;num_steps];
        end
        V1 = V(V>0);% velocity without zero and without NaN
        D1 = D(D>0);
        L = 2; % typical distance between electrodes (in cm)
        N = 6; % number of interpolation points betwen electrodes
        l = L/N; % distance between interpolation points, in cm
        fs = 250;
        Ts = 1/fs;
        v_unit = (l/Ts)/100;% unit velocity (jumping one pixel in one sampling time), in meters
        v_corrected = V1*v_unit; % velocity in meters per second
        v_corr_mean = mean(v_corrected);
        v_corr_std = std(v_corrected);
        save([folder_freq_save,'Velocity'],"v_unit","v_corrected","L","N",...
            "l","fs","D1","v_corr_std","v_corr_mean");
    end
end



