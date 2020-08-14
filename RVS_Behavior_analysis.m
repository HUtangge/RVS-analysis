%%
clear all

% add analysis pathway
project_dir = '';
log_dir = fullfile(project_dir, ['raw_data/logfiles']);

subj_ID = {'01', '02', '03', '04', '09', '10', '11', '12', '13', '14',...
           '15', '16', '17', '18', '19', '20', '21', '22', '23', '24',...
           '25', '26', '27', '28'};
       
valid_runs = repmat(1:8, [length(subj_ID) 1]);
valid_runs(2,:) =  [1 2 3 4 5 6 8 9];
valid_runs(6,:) =  [2 3 4 5 6 7 8 9];
valid_runs(8,:) =  [1 2 4 5 6 7 8 9];
valid_runs(10,:) = [1 2 3 4 5 6 8 9];

%% creat matrix for all the subjects
Trial_data = [];
for i = 1:length(subj_ID)
    infiles = fullfile(log_dir, ['*' 'sj_' subj_ID{i} '_*.mat']);
    files{i} = dir(infiles);
    Subj_data = [];
    for j = 1:length(files{i})
        filename = files{i}(j,1).name;
        run = char(regexp(filename, '_', 'split'));
        run = str2num(run(5));
        if files{i}(j,1).bytes > 4000 & ismember(run, valid_runs(i,:))
            load(filename);            
            Subj_data = [Subj_data [log_RVS.Design; log_RVS.responses]];
        end
    end   
    Trial_data(:,:,i) = Subj_data;
end
%% Some description statistics
perf_rehearsal      = sum(Trial_data(16,:)==1) /length(Trial_data(16,:));% performance 
rt_rehearsal        = sum(Trial_data(17,Trial_data(17,:)~=8888)) / length(Trial_data(17,Trial_data(17,:)~=8888));% mean responce time
perf_6_rehearsal    = sum(Trial_data(16,Trial_data(10,:)==7500)==1) / sum(Trial_data(10,:)==7500);% performance for 6 rehearsal
perf_8_rehearsal    = sum(Trial_data(16,Trial_data(10,:)==10000)==1) / sum(Trial_data(10,:)==10000);% performance for 8 rehearsal
perf_10_rehearsal   = sum(Trial_data(16,Trial_data(10,:)==12500)==1) / sum(Trial_data(10,:)==12500);% performance for 10 rehearsal
perf_4_stimulation  = sum(Trial_data(16,Trial_data(6,:)==5000)==1) / sum(Trial_data(6,:)==5000);% performance for 4 stimulation
perf_6_stimulation  = sum(Trial_data(16,Trial_data(6,:)==7500)==1) / sum(Trial_data(6,:)==7500);% performance for 6 stimulation
perf_8_stimulation  = sum(Trial_data(16,Trial_data(6,:)==10000)==1) / sum(Trial_data(6,:)==10000);% performance for 8 stimulation
%%
data =[];
for i = 1:24%length(subj_ID)
    perf_4_6    = sum((Trial_data(16,(Trial_data(6,:,i)==5000 & Trial_data(10,:,i)==7500),i)==1)) / sum((Trial_data(6,:,i)==5000) & (Trial_data(10,:,i)==7500));
    perf_4_8    = sum((Trial_data(16,(Trial_data(6,:,i)==5000 & Trial_data(10,:,i)==10000),i)==1)) / sum((Trial_data(6,:,i)==5000) & (Trial_data(10,:,i)==10000));
    perf_4_10    = sum((Trial_data(16,(Trial_data(6,:,i)==5000 & Trial_data(10,:,i)==12500),i)==1)) / sum((Trial_data(6,:,i)==5000) & (Trial_data(10,:,i)==12500));
    perf_6_6    = sum((Trial_data(16,(Trial_data(6,:,i)==7500 & Trial_data(10,:,i)==7500),i)==1)) / sum((Trial_data(6,:,i)==7500) & (Trial_data(10,:,i)==7500));
    perf_6_8    = sum((Trial_data(16,(Trial_data(6,:,i)==7500 & Trial_data(10,:,i)==10000),i)==1)) / sum((Trial_data(6,:,i)==7500) & (Trial_data(10,:,i)==10000));
    perf_6_10    = sum((Trial_data(16,(Trial_data(6,:,i)==7500 & Trial_data(10,:,i)==12500),i)==1)) / sum((Trial_data(6,:,i)==7500) & (Trial_data(10,:,i)==12500));
    perf_8_6    = sum((Trial_data(16,(Trial_data(6,:,i)==10000 & Trial_data(10,:,i)==7500),i)==1)) / sum((Trial_data(6,:,i)==10000) & (Trial_data(10,:,i)==7500));
    perf_8_8    = sum((Trial_data(16,(Trial_data(6,:,i)==10000 & Trial_data(10,:,i)==10000),i)==1)) / sum((Trial_data(6,:,i)==10000) & (Trial_data(10,:,i)==10000));
    perf_8_10    = sum((Trial_data(16,(Trial_data(6,:,i)==10000 & Trial_data(10,:,i)==12500),i)==1)) / sum((Trial_data(6,:,i)==10000) & (Trial_data(10,:,i)==12500));
    data = [data; [perf_4_6 perf_4_8 perf_4_10 perf_6_6 perf_6_8 perf_6_10 perf_8_6 perf_8_8 perf_8_10]];
end

%% Overview of the behavioral results
total_perf = data(:);
mean(total_perf)
std(total_perf)/sqrt(23)

%% reaction time 
rt_all =[];
for i = 1:24%length(subj_ID)
    perf_4_6    = sum(Trial_data(17,(Trial_data(6,:,i)==5000 & Trial_data(10,:,i)==7500 & Trial_data(17,:,i) ~= 8888),i)) / length(Trial_data(17,(Trial_data(6,:,i)==5000 & Trial_data(10,:,i)==7500 & Trial_data(17,:,i) ~= 8888),i));
    perf_4_8    = sum(Trial_data(17,(Trial_data(6,:,i)==5000 & Trial_data(10,:,i)==10000 & Trial_data(17,:,i) ~= 8888),i)) / length(Trial_data(17,(Trial_data(6,:,i)==5000 & Trial_data(10,:,i)==10000 & Trial_data(17,:,i) ~= 8888),i));
    perf_4_10    = sum(Trial_data(17,(Trial_data(6,:,i)==5000 & Trial_data(10,:,i)==12500 & Trial_data(17,:,i) ~= 8888),i)) / length(Trial_data(17,(Trial_data(6,:,i)==5000 & Trial_data(10,:,i)==12500 & Trial_data(17,:,i) ~= 8888),i));
    perf_6_6    = sum(Trial_data(17,(Trial_data(6,:,i)==7500 & Trial_data(10,:,i)==7500 & Trial_data(17,:,i) ~= 8888),i)) / length(Trial_data(17,(Trial_data(6,:,i)==7500 & Trial_data(10,:,i)==7500 & Trial_data(17,:,i) ~= 8888),i));
    perf_6_8    = sum(Trial_data(17,(Trial_data(6,:,i)==7500 & Trial_data(10,:,i)==10000 & Trial_data(17,:,i) ~= 8888),i)) / length(Trial_data(17,(Trial_data(6,:,i)==7500 & Trial_data(10,:,i)==10000 & Trial_data(17,:,i) ~= 8888),i));
    perf_6_10    = sum(Trial_data(17,(Trial_data(6,:,i)==7500 & Trial_data(10,:,i)==12500 & Trial_data(17,:,i) ~= 8888),i)) / length(Trial_data(17,(Trial_data(6,:,i)==7500 & Trial_data(10,:,i)==12500 & Trial_data(17,:,i) ~= 8888),i));
    perf_8_6    = sum(Trial_data(17,(Trial_data(6,:,i)==10000 & Trial_data(10,:,i)==7500 & Trial_data(17,:,i) ~= 8888),i)) / length(Trial_data(17,(Trial_data(6,:,i)==10000 & Trial_data(10,:,i)==7500 & Trial_data(17,:,i) ~= 8888),i));
    perf_8_8    = sum(Trial_data(17,(Trial_data(6,:,i)==10000 & Trial_data(10,:,i)==10000 & Trial_data(17,:,i) ~= 8888),i)) / length(Trial_data(17,(Trial_data(6,:,i)==10000 & Trial_data(10,:,i)==10000 & Trial_data(17,:,i) ~= 8888),i));
    perf_8_10    = sum(Trial_data(17,(Trial_data(6,:,i)==10000 & Trial_data(10,:,i)==12500 & Trial_data(17,:,i) ~= 8888),i)) / length(Trial_data(17,(Trial_data(6,:,i)==10000 & Trial_data(10,:,i)==12500 & Trial_data(17,:,i) ~= 8888),i));
    rt_all = [rt_all; [perf_4_6 perf_4_8 perf_4_10 perf_6_6 perf_6_8 perf_6_10 perf_8_6 perf_8_8 perf_8_10]];
end

%%
xvalues = {'7.5','10','12.5'};
yvalues = {'5','7.5','10'};
h = heatmap(xvalues,yvalues,data);

h.Title = 'Performance in Different Trials';
h.XLabel = 'Duration of Rehearsal';
h.YLabel = 'Duration of Stimulation';

%% Performance based on stimulus
perf_1 = sum(Trial_data(16,(Trial_data(11,:)==1))==1) / length(Trial_data(16,(Trial_data(11,:)==1))==1);
perf_2 = sum(Trial_data(16,(Trial_data(11,:)==2))==1) / length(Trial_data(16,(Trial_data(11,:)==2))==1);
perf_3 = sum(Trial_data(16,(Trial_data(11,:)==3))==1) / length(Trial_data(16,(Trial_data(11,:)==3))==1);
perf_4 = sum(Trial_data(16,(Trial_data(11,:)==4))==1) / length(Trial_data(16,(Trial_data(11,:)==4))==1);

perfc_1 = sum(Trial_data(16,(Trial_data(12,:)==1))==1) / length(Trial_data(16,(Trial_data(12,:)==1))==1);
perfc_2 = sum(Trial_data(16,(Trial_data(12,:)==2))==1) / length(Trial_data(16,(Trial_data(12,:)==2))==1);
perfc_3 = sum(Trial_data(16,(Trial_data(12,:)==3))==1) / length(Trial_data(16,(Trial_data(12,:)==3))==1);
perfc_4 = sum(Trial_data(16,(Trial_data(12,:)==4))==1) / length(Trial_data(16,(Trial_data(12,:)==4))==1);

%% Plotting
a=bar([rt', rt_c', rt_ic']);
ylim([250 1200]);
xlim([1 8]);
xlabel('Block number');
ylabel('Reaction time (ms)');
title('Mean Reaction Time of Different Blocks');
legend('mean rt', 'rt-correct trials', 'rt-incorrect trials');

hold on 
rt_mean = sum(Trial_data(17,Trial_data(17,:)~=8888)) / length(Trial_data(17,Trial_data(17,:)~=8888));
line([0 9], [rt_mean rt_mean]);
legend('mean rt', 'correct trials rt', 'incorrect trials rt', 'grand mean');
%h.DisplayName = 'mean';

hold on
errorbar_groups([rt', rt_c', rt_ic']',[sd', sd_c', sd_ic']');
xlim([0 9]);


%% Statistics 
[~,~,stats] = anova1(rt_data');
[c,m,h,nms] = multcompare(stats);

%% Plotting the performance
plot(perf, 'b-o');
ylim([60 100]);
xlim([0 9]);
xlabel('Run number', 'FontSize', 16);
ylabel('Performance (%)', 'FontSize', 16);
title('Performance of Different Runs', 'FontSize', 18);
hold on 
perf_mean = 100* (sum(Trial_data(16,:,j)==1) /length(Trial_data(16,:,j)));
h = line([0 9], [perf_mean perf_mean], 'DisplayName', 'grand mean');
legend(h);

boxplot(Subj_perf);
hold on
avg_rt = mean(rt_all);
SEM = std(rt_all)/sqrt(23);
CI = SEM*1.96;

for i = 1 :8
    x = ones(24,1).*i;
    scatter(x, rt_all(:,i));
    hold on
    h = errorbar(i, avg_rt(i), CI(i), 'k', 'linestyle', 'none', 'linewidth', 28);
    errorbarT(h,0.05,1.5);
    hold on    
end
plot(avg_rt);
