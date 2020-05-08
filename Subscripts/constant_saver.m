%% File containing constants the pipeline requires

%% Clean up
clearvars
close all
Paths
%% Define the constants

% allocate the structure
constants = struct([]);

% threshold percentile of the traces to use, utilized by Stage2 (which then
% generates file per fish for many other scripts)
constants(1).perc = 40;

% minimum number of stimuli to threshold active traces, also utilized in
% Stage2 as above
constants(1).stimuli = 7;

% correlation threshold to merge clusters (above this, they get merged)
constants(1).correlation = 0.8;

% minimum number of traces to keep a cluster
constants(1).num_thres = 100;

% p value for significance of the cluster correlations
constants(1).pval = 0.05;
%% Define the af labels
af_labels = struct([]);
n=1;
af_labels(n).name = 'SO 1-2';
af_labels(n).number = '11';
n = n+1;
af_labels(n).name = 'SFGS 1-2';
af_labels(n).number = '12';
n = n+1;
af_labels(n).name = 'SFGS 3-4';
af_labels(n).number = '13';
n = n+1;
af_labels(n).name = 'SFGS 5-6';
af_labels(n).number = '14';
n = n+1;
af_labels(n).name = 'SGC';
af_labels(n).number = '15';
n = n+1;
af_labels(n).name = 'SAC/SPV';
af_labels(n).number = '16';
n = n+1;
af_labels(n).name = 'AF4';
af_labels(n).number = '4';
n = n+1;
af_labels(n).name = 'AF5';
af_labels(n).number = '6';
n = n+1;
af_labels(n).name = 'AF6';
af_labels(n).number = '5';
n = n+1;
af_labels(n).name = 'AF7';
af_labels(n).number = '7';
n = n+1;
af_labels(n).name = 'AF8';
af_labels(n).number = '8';
n = n+1;
af_labels(n).name = 'AF9d';
af_labels(n).number = '9';
n = n+1;
af_labels(n).name = 'AF9v';
af_labels(n).number = '3';
n = n+1;


% list of AFs contained in every data set, utilized by several stages
% including Stage3 and others
%determined manually by Clemens
af_list = cell(39,1);

af_list{1} = cell(2,1);
af_list{1}{1} = 'A3_pre';
af_list{1}{2} = cell(4,1);
af_list{1}{2}{1} = [12 13 14 15 16];
af_list{1}{2}{2} = [12 13 14 15 16 9 8 3 5 4];
af_list{1}{2}{3} = [4 6 3 5];
af_list{1}{2}{4} = [6 3 4 5];

af_list{2} = cell(2,1);
af_list{2}{1} = 'A4_pre';
af_list{2}{2} = cell(4,1);
af_list{2}{2}{1} = [12 13 14 15 16];
af_list{2}{2}{2} = [9 8 3 5];
af_list{2}{2}{3} = [9 6 3 4 5];
af_list{2}{2}{4} = [6 3 4 5];

af_list{3} = cell(2,1);
af_list{3}{1} = 'A4_post';
af_list{3}{2} = cell(3,1);
af_list{3}{2}{1} = [12 13 14 15 16];
af_list{3}{2}{2} = [9 8 3 6 5];
af_list{3}{2}{3} = [6 3 4 5];

af_list{4} = cell(2,1);
af_list{4}{1} = 'A5_pre';
af_list{4}{2} = cell(3,1);
af_list{4}{2}{1} = [12 13 14 15 16];
af_list{4}{2}{2} = [9 8 3 5];
af_list{4}{2}{3} = [6 5 4];
% 
af_list{5} = cell(2,1);
af_list{5}{1} = 'A5_post';
af_list{5}{2} = cell(4,1);
af_list{5}{2}{1} = [12 13 14 15 16];
af_list{5}{2}{2} = [9 8 3 4 5];
af_list{5}{2}{3} = [12 13 14 15 16];
af_list{5}{2}{4} = [3 4 5 6];

af_list{6} = cell(2,1);
af_list{6}{1} = 'A6_pre';
af_list{6}{2} = cell(5,1);
af_list{6}{2}{1} = [12 14 15];
af_list{6}{2}{2} = [12 13 14 15 16];
af_list{6}{2}{3} = [12 13 14 15 16];
af_list{6}{2}{4} = [9 8 3 5];
af_list{6}{2}{5} = [6 5 4];
% 
af_list{7} = cell(2,1);
af_list{7}{1} = 'D2_pre';
af_list{7}{2} = cell(5,1);
af_list{7}{2}{1} = [12 13 14 15];
af_list{7}{2}{2} = [12 13 14 15 16];
af_list{7}{2}{3} = [12 13 14 15 16];
af_list{7}{2}{4} = [9 8 5 3 6 4];
af_list{7}{2}{5} = [6 5 4 3];

af_list{8} = cell(2,1);
af_list{8}{1} = 'Imax10_pre';
af_list{8}{2} = cell(3,1);
af_list{8}{2}{1} = [11 12 13 14 15];
af_list{8}{2}{2} = [12 13 14 15 16 9];
af_list{8}{2}{3} = [16 9 8 5 6 4 3];

af_list{9} = cell(2,1);
af_list{9}{1} = 'Imax10_control';
af_list{9}{2} = cell(4,1);
af_list{9}{2}{1} = [11 12 13 14 15];
af_list{9}{2}{2} = [12 13 14 15 16];
af_list{9}{2}{3} = [16 9 8];
af_list{9}{2}{4} = [9 5 3 6 4];

af_list{10} = cell(2,1);
af_list{10}{1} = 'Imax11_pre';
af_list{10}{2} = cell(5,1);
af_list{10}{2}{1} = [11 12 13 14 15];
af_list{10}{2}{2} = [12 13 14 15 16 9];
af_list{10}{2}{3} = [9 8 3 12 13 16];
af_list{10}{2}{4} = [9 8 3 6 4 12 13 14 16 5];
af_list{10}{2}{5} = [4 6];

af_list{11} = cell(2,1);
af_list{11}{1} = 'Imax12_pre';
af_list{11}{2} = cell(5,1);
af_list{11}{2}{1} = [11 12 13 14 15 16];
af_list{11}{2}{2} = [12 13 14 15 16];
af_list{11}{2}{3} = [9 8];
af_list{11}{2}{4} = [9 8 3 5 6];
af_list{11}{2}{5} = [6 5 4];

af_list{12} = cell(2,1);
af_list{12}{1} = 'Imax6_pre';
af_list{12}{2} = cell(3,1);
af_list{12}{2}{1} = [12 13 14 15 16];
af_list{12}{2}{2} = [9 8 3 5 4];
af_list{12}{2}{3} = [6 3 4 5];

af_list{13} = cell(2,1);
af_list{13}{1} = 'Imax7_pre';
af_list{13}{2} = cell(3,1);
af_list{13}{2}{1} = [12 13 14 15 16];
af_list{13}{2}{2} = [9 8 3 5 6 4];
af_list{13}{2}{3} = [6 5 4 3];

af_list{14} = cell(2,1);
af_list{14}{1} = 'Imax8_pre';
af_list{14}{2} = cell(4,1);
af_list{14}{2}{1} = [11 12 13 14 15 16];
af_list{14}{2}{2} = [9 8 3 5];
af_list{14}{2}{3} = [6 3 4 5];
af_list{14}{2}{4} = [6 4];

af_list{15} = cell(2,1);
af_list{15}{1} = 'Imax9_pre';
af_list{15}{2} = cell(4,1);
af_list{15}{2}{1} = [11 12 13 14 15 16];
af_list{15}{2}{2} = [9 8 3];
af_list{15}{2}{3} = [9 8 3 6 5];
af_list{15}{2}{4} = [6 5 4];

af_list{16} = cell(2,1);
af_list{16}{1} = 'Imax1_pre';
af_list{16}{2} = cell(3,1);
af_list{16}{2}{1} = [12 13 14 15 16];
af_list{16}{2}{2} = [9 8 6 3 4 5];
af_list{16}{2}{3} = [6 3 5];

af_list{17} = cell(2,1);
af_list{17}{1} = 'Imax2_pre';
af_list{17}{2} = cell(3,1);
af_list{17}{2}{1} = [12 13 14 15 16];
af_list{17}{2}{2} = [9 8 6 3 4 5];
af_list{17}{2}{3} = [6 3 4 5];

af_list{18} = cell(2,1);
af_list{18}{1} = 'Imax3_pre';
af_list{18}{2} = cell(4,1);
af_list{18}{2}{1} = [12 13 14 15 16];
af_list{18}{2}{2} = [9 8 6 3 4 5];
af_list{18}{2}{3} = [6 3 4 5];
af_list{18}{2}{4} = [6 3 4 5];

af_list{19} = cell(2,1);
af_list{19}{1} = 'Imax4_pre';
af_list{19}{2} = cell(4,1);
af_list{19}{2}{1} = [12 13 14 15 16];
af_list{19}{2}{2} = [9];
af_list{19}{2}{3} = [9 8 6 3 4 5];
af_list{19}{2}{4} = [6 5 4];

af_list{20} = cell(2,1);
af_list{20}{1} = 'Imax5_pre';
af_list{20}{2} = cell(6,1);
af_list{20}{2}{1} = [12 13 14 15];
af_list{20}{2}{2} = [12 13 14 15 16];
af_list{20}{2}{3} = [9 8 3 5];
af_list{20}{2}{4} = [6 3 4 5];
af_list{20}{2}{5} = [6 3 4 5];
af_list{20}{2}{6} = [6 3 4 5];

% % af_list{21} = cell(2,1);
% % af_list{21}{1} = 'Fish1_pre';
% % af_list{21}{2} = cell(3,1);
% % af_list{21}{2}{1} = 10;
% % af_list{21}{2}{2} = [10 9 8 6 5 4 3];
% % af_list{21}{2}{3} = [6 5 4];
% % 
% % af_list{22} = cell(2,1);
% % af_list{22}{1} = 'Fish3_pre';
% % af_list{22}{2} = cell(5,1);
% % af_list{22}{2}{1} = 10;
% % af_list{22}{2}{2} = [10 9 8];
% % af_list{22}{2}{3} = [10 9 8 5 4 3];
% % af_list{22}{2}{4} = [6 5 4 3];
% % af_list{22}{2}{5} = [6 5 4 3];

af_list{23} = cell(2,1);
af_list{23}{1} = 'D2_post';
af_list{23}{2} = cell(4,1);
af_list{23}{2}{1} = [12 13 14 15 16];
af_list{23}{2}{2} = [9 8 5 4 3];
af_list{23}{2}{3} = [6 5 4 3];
af_list{23}{2}{4} = [6 4];

af_list{24} = cell(2,1);
af_list{24}{1} = 'Imax3_post';
af_list{24}{2} = cell(2,1);
af_list{24}{2}{1} = [12 13 14 15 16];
af_list{24}{2}{2} = [9 8 6 3 4 5];

af_list{25} = cell(2,1);
af_list{25}{1} = 'Imax6_post';
af_list{25}{2} = cell(2,1);
af_list{25}{2}{1} = [12 13 14 15 16];
af_list{25}{2}{2} = [9 8 3 4 5];

af_list{26} = cell(2,1);
af_list{26}{1} = 'Imax8_post';
af_list{26}{2} = cell(5,1);
af_list{26}{2}{1} = [12 13 14 15 16];
af_list{26}{2}{2} = [12 13 14 15 16];
af_list{26}{2}{3} = [9 8 3 4 5];
af_list{26}{2}{4} = [6 3 4 5];
af_list{26}{2}{5} = [6 4];

af_list{27} = cell(2,1);
af_list{27}{1} = 'Imax11_control';
af_list{27}{2} = cell(3,1);
af_list{27}{2}{1} = [11 12 13 14 15 16];
af_list{27}{2}{2} = [13 14 15 16 9 8];
af_list{27}{2}{3} = [9 8 6 3 4 5];

af_list{28} = cell(2,1);
af_list{28}{1} = 'Imax12_post';
af_list{28}{2} = cell(4,1);
af_list{28}{2}{1} = [11 12 13 14 15];
af_list{28}{2}{2} = [12 13 14 15 16];
af_list{28}{2}{3} = [9 8 6 3 4 5];
af_list{28}{2}{4} = [6 3 4 5];

af_list{29} = cell(2,1);
af_list{29}{1} = 'A6_post';
af_list{29}{2} = cell(2,1);
af_list{29}{2}{1} = [12 13 14 15 16];
af_list{29}{2}{2} = [9 8 6 3 4 5];

af_list{30} = cell(2,1);
af_list{30}{1} = 'X2_pre';
af_list{30}{2} = cell(3,1);
af_list{30}{2}{1} = [11 12 13 14 15 16];
af_list{30}{2}{2} = [9 8 7 6 5 4 3];
af_list{30}{2}{3} = [6 5 4];

af_list{31} = cell(2,1);
af_list{31}{1} = 'X2_control';
af_list{31}{2} = cell(3,1);
af_list{31}{2}{1} = [11 12 13 14 15 16];
af_list{31}{2}{2} = [9 8 6 5 4 3];
af_list{31}{2}{3} = [6 5 4];

af_list{32} = cell(2,1);
af_list{32}{1} = 'X3_pre';
af_list{32}{2} = cell(3,1);
af_list{32}{2}{1} = [11 12 13 14 15 16];
af_list{32}{2}{2} = [9 8 7 6 5 4 3];
af_list{32}{2}{3} = [6 5 4];

af_list{33} = cell(2,1);
af_list{33}{1} = 'X3_control';
af_list{33}{2} = cell(3,1);
af_list{33}{2}{1} = [11 12 13 14 15 16];
af_list{33}{2}{2} = [9 8 7 6 5 4 3];
af_list{33}{2}{3} = [6 5 4];

af_list{34} = cell(2,1);
af_list{34}{1} = 'X4_pre';
af_list{34}{2} = cell(2,1);
af_list{34}{2}{1} = [12 13 14 15 16];
af_list{34}{2}{2} = [9 8 7 6 3 4 5];

af_list{35} = cell(2,1);
af_list{35}{1} = 'X4_control';
af_list{35}{2} = cell(3,1);
af_list{35}{2}{1} = [12 13 14 15 16];
af_list{35}{2}{2} = [9 8 7 6 3 5 4];
af_list{35}{2}{3} = [6];

af_list{36} = cell(2,1);
af_list{36}{1} = 'X5_pre';
af_list{36}{2} = cell(2,1);
af_list{36}{2}{1} = [12 13 14 15 16];
af_list{36}{2}{2} = [9 8 7 3 5 4 6];

af_list{37} = cell(2,1);
af_list{37}{1} = 'X5_control';
af_list{37}{2} = cell(2,1);
af_list{37}{2}{1} = [12 13 14 15 16];
af_list{37}{2}{2} = [9 4 8 7 6 5 3];

af_list{38} = cell(2,1);
af_list{38}{1} = 'X6_pre';
af_list{38}{2} = cell(3,1);
af_list{38}{2}{1} = [12 13 14 15 16];
af_list{38}{2}{2} = [12 13 14 15 16];
af_list{38}{2}{3} = [9 8 7 6 4 3 5];

af_list{39} = cell(2,1);
af_list{39}{1} = 'X6_control';
af_list{39}{2} = cell(2,1);
af_list{39}{2}{1} = [12 13 14 15 16];
af_list{39}{2}{2} = [9 8 7 6 4 3 5];


af_list = af_list(~cellfun('isempty', af_list'));
%% Save the aforelisted constants

%save path
% save_path = 'E:\Behavioral data\Matlab\AF_proc\Clemens_suite\20170827_Software_pipeline\subfunctions and scripts\pipeline_constants.mat';
save_path = constants_path;
save(save_path,'constants','af_list','af_labels')