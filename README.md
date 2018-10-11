
# Matlab-MEG-analysis
A Matlab exam taken with grade 20/22


%% Load the files
% set directory 
```
clear all
close all
clc
od = cd('D:\Matlab\Exam Data')

load data.mat
load rFG.mat
load time.mat
```

%% Question 1 

% number of trials can be obtained by looking at the rows of bubbles, which
     %is also the number of stimuli
% number of responses can be obrained by looking at rows in data.resp
% number of meg data points can be obrained by looking at data.meg rows
```
nsj = length(data); %number of subjects

% Check if there are missin (NaN) values in the data of each participant
for i = 1:nsj 
    if find(isnan(data(i).bubs))
        disp('NaNs in bubs')
    elseif find(isnan(data(i).meg))
        disp('NaNs in meg')
    elseif find(isnan(data(i).resp))
        disp('NaNs in resp')
    else 
        disp('No NaNs')
    end
end
```
%%   
% Create a matrix where each column is a measure(meg, stim, response size)
      % and each row is a participant
```
      
nums = zeros(5, 3); % empty matrix to be used by the for loop

for i = 1:nsj

    % number of stimuli and trials
    nums(i, 2) = size(data(i).bubs, 1);
    % number of behavioural responses
    nums(i, 3) = size(data(i).resp, 2);
    % number of brain responses
    nums(i, 1) = size(data(i).meg, 1);
end
```
% create names 
```
sizes_names = {'Stimuli/Trials' 'Behavioural responses' 'Brain responses(MEG)'};

% Plot the results
bAxes = figure
bar(nums)
colormap hot
xlabel('Participant number')
ylabel('Count')
colormap(bAxes, hot)
legend(sizes_names)
title('Count of available trials and responses for each participant')
```
% Did all subjects complete the same number of trials and is data recorded for each trial? 
% ANSWER: Each subject has the same number of trials, stimuli, brain data.
% But the number of data is different between subjects. 

%% Question 2 

% To find the participant percentage for each category of the 3 responses,
% we can first write a logical statement, finding the length of each
% response for each subject in a loop.
% Then obtain the sum for each participant, and divide each category by
% this sum.
```
for ii = 1:nsj
    
    %0 = 'nuns', 1 = 'Voltaire' and '2' = 'don't know'
    nr_resp_nuns = length(find(data(ii).resp == 0));
    nr_resp_V = length(find(data(ii).resp == 1));
    nr_resp_dn = length(find(data(ii).resp == 2)); 

    % I use a structure to store my results, so later I can have the names
             % of each response (nuns, Voltaire = V, don't know = dn) that I summed 
    
    % get the sum for each participant separately 
    resp_struct(ii).sum_each = sum([nr_resp_nuns, nr_resp_V, nr_resp_dn]);
    
    % write a field for each response, storing the percentages
    resp_struct(ii).perc_nuns = nr_resp_nuns/[resp_struct(ii).sum_each] * 100;
    resp_struct(ii).perc_V = nr_resp_V/[resp_struct(ii).sum_each] * 100;
    resp_struct(ii).perc_dn = nr_resp_dn/[resp_struct(ii).sum_each] * 100;
    
    % If the percentages are correct, sum should give 100
    resp_struct(ii).should_be_100 = sum([resp_struct(ii).perc_nuns, resp_struct(ii).perc_V, resp_struct(ii).perc_dn]);
    
end
```
% matrix to plot
```
resp_perc = [resp_struct.perc_nuns; resp_struct.perc_V; resp_struct.perc_dn]'; 

resp_names = {'Nuns' 'Voltaire' 'Don`t know'};

figure('Name', 'Q2', 'position', [0 0 700 500])
bar(resp_perc, 'stacked')
xlabel('Participant number')
ylabel('Percentage')
ylim([0 100])
colormap hot
title('Percentage responses for the three categories for each subject')
leg = legend(resp_names, 'Location', 'westoutside')
title(leg, 'Response type')
```
% By the plot it can be see overall how much of each option, each
% participant has chosen.

%% Question 3
% Firs I extract the data where the response is Voltaire or Nun
% (separately) for each participant. Then, separately for Nun or Voltaire,
% I run the pixels for each participant for all trials thought each spatial
% frequency. Finally, I sum these resulting matrices for each spatial
% frequency and plot. I do the same for Voltaire and plot.


% index the trials, where the response was nuns
```
for i = 1:nsj 
    all_trials = data(i).bubs; % extract each participant bubs
    nuns = data(i).resp == 0;  % nuns were coded with 0
    nun_trials(i).pixels = all_trials((find(nuns)), :);
end

% index the trials, where the response was Voltaire
for i = 1:nsj 
    all_trials = data(i).bubs;
    Vol = data(i).resp == 1;  % Voltaire was coded with 1
    V_trials(i).pixels = all_trials((find(Vol)), :);
end
```
%% Q#3 Sum Pixels (for Nuns)
% Sum pixels - For every subject - every spatial frequency - sf1, sf2, sf3, sf4, sf5
% set working directory 
```
cd([od '\utils'])

% Starting with a matrix of defined size makes Matlab faster (#Christoph_tips)
final_nun_sf1 = zeros(64, 64, 5);
   
for i = 1:5
    
  trials_n_per_sj = length(nun_trials(i).pixels(:,1));
  nun1 = zeros(64, 64, trials_n_per_sj);

    % loop taking each trial of each subject and running sf1 on each trial
      for ii = 1:trials_n_per_sj;
          % take each trial of each subject
          sj_trial = nun_trials(i).pixels(ii, :);
          % nun1 matrix for sf1, 64x64xTrials, 3rd dimension is the number
                     % of trials for each subject 
          nun1(:, :, ii) = sf1(sj_trial);
          nun2(:, :, ii) = sf2(sj_trial);
          nun3(:, :, ii) = sf3(sj_trial);
          nun4(:, :, ii) = sf4(sj_trial);
          nun5(:, :, ii) = sf5(sj_trial);
          
      end 
      
      % 3D matrix 64x64x5, all trials summed, for each subject
      final_nun_sf1(:, :, i) = sum(nun1, 3); 
      final_nun_sf2(:, :, i) = sum(nun2, 3); 
      final_nun_sf3(:, :, i) = sum(nun3, 3); 
      final_nun_sf4(:, :, i) = sum(nun4, 3); 
      final_nun_sf5(:, :, i) = sum(nun5, 3); 
end
%% Q#3 Same for Voltaire

final_volt_sf1 = zeros(64, 64, 5);
for i = 1:5
    
  trials_n_per_sj = length(V_trials(i).pixels(:,1));
  volt1 = zeros(64, 64, trials_n_per_sj);

    % loop taking each trial of each subject and running sf1 on each trial
      for ii = 1:trials_n_per_sj;
          % take each trial of each subject
          sj_trial = V_trials(i).pixels(ii, :);
          % nun1 matrix for sf1, 64x64xTrials, 3rd dimension is the number
                     % of trials for each subject 
          volt1(:, :, ii) = sf1(sj_trial);
          volt2(:, :, ii) = sf2(sj_trial);
          volt3(:, :, ii) = sf3(sj_trial);
          volt4(:, :, ii) = sf4(sj_trial);
          volt5(:, :, ii) = sf5(sj_trial);
          
      end 
      
      % 3D matrix 64x64x5, all trials summed, for each subject
      final_volt_sf1(:, :, i) = sum(volt1, 3); 
      final_volt_sf2(:, :, i) = sum(volt2, 3); 
      final_volt_sf3(:, :, i) = sum(volt3, 3); 
      final_volt_sf4(:, :, i) = sum(volt4, 3); 
      final_volt_sf5(:, :, i) = sum(volt5, 3); 
end
```

%% Q#3 Plot NUNS
% I have two plots for Nuns and Voltaire trials, across participants and
% spatial frequencies. I could plot either the pixels superimposed on
% Dali's image, or plain pixel maps. I decide to show both. For nuns I
% superimposed the pixels on the Dali's image with the provided function
% plot_pixel_effect. For Voltaire trials, I use imagesc to show the pixels
% (and to show that I can name the axes and put separate titles, as taught in lectures).

```
figure('Name', 'Q3Nuns', 'position', [0 0 1500 900])
for s = 1:5
    pos = [1:5; 6:10; 11:15; 16:20; 21:25];
    n = pos(s, :);
        subplot(5, 5, n(1))   % There is a subplot for each sf and each participant
        plot_pixel_effect(final_nun_sf1(:, :, s))
        subplot(5, 5, n(2))
        plot_pixel_effect(final_nun_sf2(:, :, s))
        subplot(5, 5, n(3))
        plot_pixel_effect(final_nun_sf3(:, :, s))
        subplot(5, 5, n(4))
        plot_pixel_effect(final_nun_sf4(:, :, s))
        subplot(5, 5, n(5))
        plot_pixel_effect(final_nun_sf5(:, :, s))
    
end
gcf, 
suptitle('Nuns Trials')
```
% It can be seen that there is fine detail with sf1 and less detail with
% sf5.


%% Q#3 Plot Voltaire with pixels
```
figure('Name', 'Q3Voltaire', 'position', [0 0 1500 900])
for s = 1:5
    pos = [1:5; 6:10; 11:15; 16:20; 21:25];
    n = pos(s, :);
        subplot(5, 5, n(1)) 
        imagesc(final_volt_sf1(:, :, s)), axis image  
        title(['Sj 1 with sf' num2str(s)])
        xlabel('n of pixels')
        ylabel('n of pixels')
        subplot(5, 5, n(2))
        imagesc(final_volt_sf2(:, :, s)), axis image
        title(['Sj 2 with sf' num2str(s)])
        xlabel('n of pixels')
        ylabel('n of pixels')
        subplot(5, 5, n(3))
        imagesc(final_volt_sf3(:, :, s)), axis image
        title(['Sj 3 with sf' num2str(s)])
        ylabel('n of pixels')
        xlabel('n of pixels')
        subplot(5, 5, n(4))
        imagesc(final_volt_sf4(:, :, s)), axis image
        title(['Sj 4 with sf' num2str(s)])
        xlabel('n of pixels')
        ylabel('n of pixels')
        subplot(5, 5, n(5))
        imagesc(final_volt_sf5(:, :, s)), axis image
        title(['Sj 5 with sf' num2str(s)])
        xlabel('n of pixels')
        ylabel('n of pixels')
    
end
```
%% Question 4 Find the sum across participants

% For this task we need to sum the matrices of all subjects across spatial
% frequencies. E.g subj1_sf1 + subj2_sp1 and so on...This is what I do in
% this section.

% 3rd dimension of final_nun_sf1 has 5 tabs, and represents each participant 
% Sum across the 3rd dimension to get a matrix representing all
% participants.
% For Nuns
```
all_data(1).nun = sum(final_nun_sf1, 3);
all_data(2).nun = sum(final_nun_sf2, 3);
all_data(3).nun = sum(final_nun_sf3, 3);
all_data(4).nun = sum(final_nun_sf4, 3);
all_data(5).nun = sum(final_nun_sf5, 3);
   
% For Voltaire
all_data(1).volt = sum(final_volt_sf1, 3);
all_data(2).volt = sum(final_volt_sf2, 3);
all_data(3).volt = sum(final_volt_sf3, 3);
all_data(4).volt = sum(final_volt_sf4, 3);
all_data(5).volt = sum(final_volt_sf5, 3);   
  ```
%% Q#4 Plot the nuns and voltaire one next to the other
```
pos = (1:2:10); % prepare a variable to set position
labels_all = [64 32 16 8 4]; % add labels
sp_hand1 = figure('Name', 'Q4', 'position', [0 0 500 750])
for s = 1:5
    
    for n = pos(1, s);
        subplot(5, 2, n)
        imagesc(all_data(s).volt), axis image
        title(['Voltaire trials with sf' num2str(s)])
        xlabel('n of pixels')
        ylabel('n of pixels')
        subplot(5, 2, n+1)
        imagesc(all_data(s).nun), axis image
        title(['Nuns trials with sf' num2str(s)])
        colormap jet
        xlabel('n of pixels')
        ylabel('n of pixels')

    end
    
end
```

%% Question 5 
% Here I find the maximum of each of the summed matrices from Q#4. 
% Devide each summed matrix by its maximum value to get values between 0
% and 1. Then for Nuns and Voltaire exctract the pixels bigger than 0.9,
% and make everything else 0.
% Make an empty matrix the same size as the original one to get colour red
% in the 3d space [ 1 0 0 ]. For Voltaire would flip the 3d dimension to
% get blue [ 0 0 1 ].
% 
```
% Find the maximum of the summed matrices in each spatial frequnecy
      % separately for voltaire trials and nuns trials
n = length(all_data);

for i = 1:n
    max1 = [all_data(i).nun];
    all_data(i).nun_max = max(max1(:));
    max2 = [all_data(i).volt];
    all_data(i).volt_max = max(max2(:));
end

% devide each summed matrix by its maximum value
for i = 1:n
    
   all_data(i).one_matr_n = all_data(i).nun ./ all_data(i).nun_max
   all_data(i).one_matr_v = all_data(i).volt ./ all_data(i).volt_max
   
end
```
%% Q#5  Extract the HOT points for nuns and Voltaire where the value is bigger than 0.9
```
for i = 1:n
   % Extract the matrices for nuns
   hot_sf_nuns = [all_data(i).one_matr_n];
   hot_nun = hot_sf_nuns > 0.9;  % here are some hot nuns
   % index the hot nuns from the matrix (leaving all_data else 0)
   hot_sf_nuns(~hot_nun) = 0; 
   % write a new field in the structure
   all_data(i).hot_nuns = hot_sf_nuns;
   
   % Do the same for hot Voltaire
   hot_sf_v = [all_data(i).one_matr_v];
   hot_v = hot_sf_v > 0.9; 
   hot_sf_v(~hot_v) = 0; 
   all_data(i).colr_v = hot_sf_v;
   
end
%% Plot Q#5 

   % Make a 3D matrix for Voltaire
 for i = 1:5
     %get the matrix with maximums 
   hot_sf_v = [all_data(i).one_matr_n];
   
   dim2 = zeros(size(hot_sf_v, 1), size(hot_sf_v, 2), 2);
   all_data(i).image_3d_v = cat(3, dim2, hot_sf_v);
 end
 
  % Make a 3D matrix for nuns
for i = 1:5
   hot_sf_n = [all_data(i).one_matr_v];
   dim2 = zeros(size(hot_sf_n, 1), size(hot_sf_n, 2), 2);
   all_data(i).image_3d_n = cat(3, hot_sf_n, dim2);
end

  % Make a 3D matrix for both nuns and Voltaire
for i = 1:5
    hot_sf_n = [all_data(i).one_matr_n];
   dim2 = zeros(size(hot_sf_n, 1), size(hot_sf_n, 2), 2);
   all_data(i).image_3d_n = cat(3, hot_sf_n, dim2);
   
   hot_sf_v = [all_data(i).one_matr_v];
   all_data(i).image_3d_v = cat(3, dim2, hot_sf_v);

   dim1 = zeros(size(hot_sf_n, 1), size(hot_sf_n, 2), 1);
   hot_sf_n3d = cat(3, hot_sf_n, dim1);
   all_data(i).image_3d_both = cat(3, hot_sf_n3d, hot_sf_v);
   
end
```
%% Q#5 Plot them together
```
pos1 = [1 4 7 10 13; 2 5 8 11 14; 3 6 9 12 15];

%Define colorbars

col_v = all_data.image_3d_v;
col_v = unique(col_v);
col_3d = zeros(size(col_v, 1), 2);
col_bar_n = cat(2, col_v, col_3d);

figure('position', [0 0 700 900])
for s = 1:5
    
    for n = pos1(1, s);
        subplot(5, 3, n)
        imagesc(all_data(s).image_3d_n), axis image
        title(['Nuns hot points sf' num2str(s)])
        colorbar
        colormap(col_bar_n)
        subplot(5, 3, n+1)
        imagesc(all_data(s).image_3d_v), axis image
        title(['Voltaire hot points sf' num2str(s)])
        subplot(5, 3, n+2)
        imagesc(all_data(s).image_3d_both), axis image
        title(['Both hot points sf' num2str(s)])
        
    end 
end
```
%% Q#6
% Here I extract the nun or Voltaire trials, and index the 3D MEG matrix
% with either of these trials. Then I find the median across all 66 time
% points for each subject and write it in a new 3D matrix.
% Then plot the results. 
```
nsj = length(data);
% Extract nun
% index the trials, where the response was nuns

% Matrix to be written on with medians.
meg_dat_nun_median = zeros(5, 66, 7);
for i = 1:nsj 
    for ii = 1:7
        
    % Find the nun trials
    nun_trials = find(data(i).resp == 0);
    % define the participant number
    meg_dat = data(i).meg; 
    % Extract each voxel to work with
    meg_dat_vox = meg_dat(:, :, ii);
    % index the nun trials, and leave all columns
    nun_dat = meg_dat_vox(nun_trials, :);
    % median across the 66 time points and write (participant, 66, voxel)
    meg_dat_nun_median(i, :, ii) = median(nun_dat);
    end
end
    
% Matrix to be written on with medians.
meg_dat_volt_median = zeros(5, 66, 7);
for i = 1:nsj 
    for ii = 1:7
        
    % Find the VOLTAIRE trials
    v_trials = find(data(i).resp == 1);
    % define the participant number
    meg_dat = data(i).meg; 
    % Extract each voxel to work with
    meg_dat_vox = meg_dat(:, :, ii);
    % index the nun trials, and leave all columns
    volt_dat = meg_dat_vox(v_trials, :);
    % median across the 66 time points and write (participant, 66, voxel)
    meg_dat_volt_median(i, :, ii) = median(volt_dat);
    end
end

```
% Plot the resuls
```
figure('position', [0 0 2000 400])
for vox = 1:7
    for i = 1:5
        pos = [1 2 3 4 5 6 7];
        subplot(1, 7, pos(vox))
        plot1 = plot(meg_dat_nun_median(i,:,vox), 'r');
        xlabel('Time Points')
        ylabel('MEG signal')
        xlim([0 66])
        hold on
        plot2 = plot(meg_dat_volt_median(i,:,vox), 'b');
        
        
    end
end
suptitle('MEG signla across time, for each voxel and all participants, blue is Voltaire, red is Nuns')
```


%% Q#7 
% Make a new folder and print each Figure. Keep the figure from previous
% sections open, so that there is something available to print. 

mkdir 'New_Folder'
cd('O:\GU_2017\MATLAB\Course\Exam 2017-18\Exam 2017-18\Exam Data\utils\New_Folder')

% Save Q2 Plot
print('Q2', '-dtiff')

% Save Q3 Plot
print('Q3Nuns', '-depsc')
print('Q3Voltaire', '-depsc')

% Save Q4 Plot
print('Q4', '-djpeg')

