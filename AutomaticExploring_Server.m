%% Open dataset
% addpath(genpath('/home/guevara/toolboxes/eeglab2023.0'));
% I eliminated the prevoius line, with addpath, becuase it doesn't work,
% it gives errors: I think is because eeglab is already on the server
% addpath(genpath('/home/guevara/RESULTS'));
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
folder_data = '/home/guevara/dati_Ramon';
folder_save = '/home/guevara/RESULTS/Results_Ilaria_Patients_Subjects/no_filter/';
data_name = 'patient_006_EC.set'; % patient or subject name
save_name = 'patient_006_EC';
EEG = pop_loadset('filename',data_name,'filepath',folder_data);
mkdir([folder_save save_name] );
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% here please load existing dataset
% (e.g., C:\Users\gueva\Dropbox\...
% RamonSamirCorbetta\Data_from_Ilaria_January_2024\P04_BP4Hz.set)
% use P04_BP4Hz.set because it contains a clear pattern
% keep the fdt files

%% Visualize data for exploration;
% take all the data

data = EEG.data;
data = data';
fs = EEG.srate;
dt = 1/fs;
n = EEG.pnts;
t_total = n*dt;
t = linspace(0,t_total,n);
figure; plot(t,data);
n_cut = length(data);% take all the data
data1 = data(1:n_cut,:);
t1 = t(1:n_cut);
% figure;
% plot(t1,data1);xlabel('time(s)');
% figure; plot(t1,data1(:,1));xlabel('time(s)');% plot one channel


%% Video of my own topoplot showing maximum values;
duration_break = 30; % in seconds
num_runs = round(t_total/duration_break); % number of breaks (it gives problems if too large)
delta_t_borders = 1; % eliminate data on the borders to avoid border problems 
t_list = linspace(0+delta_t_borders,t_total-delta_t_borders,num_runs);
for k = 1:num_runs-1
    % initial and final times:
    t_initial = t_list(k);
    t_final = t_list(k+1);
    mu = 5; % threshold for acceptance of "jump"
    gamma = 3; % threshold for acceptance of "wave" (minimal number of jumps)
    L = 2; % typical distance between electrodes (in cm)
    N = 6; % number of interpolation points between electrodes
    l = L/N; % distance between interpolation points
    numtops = round((t_final-t_initial)*fs);% number of topoplots;
    time_observation = linspace(t_initial,t_final,numtops);
    sample_observation = round(time_observation*fs);
    row = zeros(numtops,1);
    column = zeros(numtops,1);
    % number of interpolation points
    %x = zeros % coordinates of interpolation points
    for i=1:numtops
        %figure;
        [h, grid_or_val, plotrad_or_grid, xmesh, ymesh] = topoplot(data1(sample_observation(i),:),EEG.chanlocs,...
            'emarker','.','electrodes','on','noplot','on');
        ymax = length(ymesh);
        intensity = grid_or_val;
        [maxint,ind] = max(intensity(:));
        [row(i),col(i)] = ind2sub(size(intensity),ind);
        %imagesc(flip(intensity,1)); hold on;
        % the flip is to do image in the correct orientation: nose up
        %plot(col(i),ymax-row(i),'*r');
        % same for ymax -row
        %pause(0.2);
    end
    t = linspace(t_initial,t_final,numtops);
    dt = t(2)-t(1);
    t1 = t(1:end-1);
    d = sqrt((diff(col')).^2 + (diff(ymax - row)).^2);
    %...........................................
    % Coded distances (no jump, small jump, large jump):
    % code for no movement: 0
    index_0 = find(d==0);
    % code for short movement: 1
    index_1 = find(d<mu & d>0);
    % code for large movement: 2
    index_2 = find(d>=mu);
    coded_d = zeros(size(d));
    coded_d(index_0) = 0;
    coded_d(index_1) = 1;
    coded_d(index_2) = 2;
    coded_d_padded = [0;coded_d;0];% it has zeros at the beginning and at the end
    d_padded = [0;d;0];
    t_padded = [t(1)-dt;t';t(end)+dt];
    index_2_in_padded = coded_d_padded == 2;
    positions_2_in_padded = find(coded_d_padded==2);
    num_2_in_coded_d_padded = sum(index_2_in_padded);
    coded_d_pieces = cell(num_2_in_coded_d_padded,1);
    d_pieces = cell(num_2_in_coded_d_padded,1);
    t_pieces = cell(num_2_in_coded_d_padded,1);
    % Example coded_d_padded: 012201202010200002022211110
    for i = 1:num_2_in_coded_d_padded-1
        coded_d_pieces {i} = coded_d_padded(positions_2_in_padded(i)+1:positions_2_in_padded(i+1)-1);
    end
    coded_d_pieces{end} = coded_d_padded(positions_2_in_padded(end)+1:end);
    for i = 1:num_2_in_coded_d_padded-1
        d_pieces {i} = d_padded(positions_2_in_padded(i)+1:positions_2_in_padded(i+1)-1);
    end
    d_pieces{end} = d_padded(positions_2_in_padded(end)+1:end);
    for i = 1:num_2_in_coded_d_padded-1
        t_pieces {i} = t_padded(positions_2_in_padded(i)+1:positions_2_in_padded(i+1)-1);
    end
    t_pieces{end} = t_padded(positions_2_in_padded(end)+1:end);
    velocity = zeros(num_2_in_coded_d_padded,1);
    for i = 1:num_2_in_coded_d_padded
        velocity(i) = sum(d_pieces{i})/length(d_pieces{i});
    end
    % Figures:
%     figure;
%     plot(t1,d);hold on;
%     plot(t1(index_1),d(index_1),'.r');
%     figure;
%     imagesc(flip(intensity,1));hold on;
%     plot(col,ymax-row,'k','LineWidth',3);
%     plot(col(index_1),ymax-row(index_1),'.r','LineWidth',3);
%     plot(col(1),ymax-row(1),'*r','MarkerSize',12);
%     figure;
%     for i=1:num_2_in_coded_d_padded
%         plot(t_pieces{i},d_pieces{i},'*');hold on;
%     end
   xmax = length(xmesh);
    save([folder_save save_name '/' num2str(t_initial) '_' num2str(t_final) '.mat'],...
        'velocity','coded_d_pieces','d','numtops','index_1', ...
        't','t1','t_initial','t_final','t_pieces','d_pieces','row','col',...
        'xmax','ymax','-v7.3');
    % save('prova.mat, 'var1','var2','-v7.3') (from Ilaria)
end

% copiar esto en un editor de texto y salvar en el server:

% #!/bin/bash
% #SBATCH --mail-user guevara.erra@gmail.com
% #SBATCH --mail-type ALL
% #SBATCH --time 5:00:00
% #SBATCH --ntasks 1
% #SBATCH --mem 128G
% #SBATCH -o AutomaticExploring_Server.txt
% #SBATCH -p brains
% matlab < AutomaticExploring_Server.m

% COMO USAR: 
% 1) abrir un terminal command window (command prompt, abrir del buscador window), poner
% ssh brain01 (o ssh brain02)
% y te va a pedir el password (el mismo de siempre no gmail)
% Si quiero trasferir files de mi compu al brain01 o 02 uso el programa pnclogin 
% que esta en el desktop compu chiquita fisica, hay que apretar boto login 
% (el login ya esta automatico) 
% y aparece una terminal con dos pantallas donde puedo tranferir de una a la otra.
% Para correr los programas de matlab desde brain01 o 02 poner en el terminal: 
% matlab -r nombre_del_script (por ejemplo matlab -r prueba)



% 2) Para correr el programa, desde brain01 o brain02 (en pantalla terminal command window)
% sbatch programa (para correrlo, hay que antes salir de matlab)
% squeue (ver la lista de programas que corren)
% scancel (si quiero cancelar lo que estoy haciendo)
% Example:
% guevara@brain01:~$ scancel 80441
% guevara@brain01:~$ sbatch AutomaticExploringServer.sh

% Para crear un file sh, crear ants un file txt en nodepad, y desde
% pnclogin cambiar su nombre poniendo .sh y eliminando .txt
% O mas simple usar AutomaticExplorinServer.sh y cambiarlo y salvarlo

% 3) Para usar viendo como en matlab, poner en el terminal:
% screen 
% matlab
% write matlab code here (escribir el codigo de matlab)
% Por ejemplo, escribir: AutomaticExploring_Server (sin el .m)

% 5) Si hay un problema tipo infinit loop, dice Paolo
% dos2unix file_matlab.m
% poi lo rilanci.

