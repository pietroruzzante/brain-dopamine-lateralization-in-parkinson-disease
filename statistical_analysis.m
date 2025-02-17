clc; close all; clear all;

%% Initial Dataset

% Load data
data = readtable("Patient_Master.csv");

% Data inspection
disp('INITIAL DATASET')
head(data)

%% Data selection and preparation

% Variable selection
data = data(:,{'ENROLL_AGE','COHORT','GENETICS','BIRTHDT','ETHNICITY','SEX',...
               'HANDED','SXDT','EDUCYRS','NP1PTOT','NP1RTOT','NP2PTOT','NP3TOT',...
               'NP4TOT','DATSCAN_CAUDATE_R','DATSCAN_CAUDATE_L','DATSCAN_PUTAMEN_R',...
               'DATSCAN_PUTAMEN_L','DATSCAN', 'DATSCAN_QUALITY_RATING','LEDD',...
               'MCATOT','NHY','HTCM','WGTKG'});

% Variable categorization
data.COHORT = double(categorical(data.COHORT, {'PD','HC','SWEDD','Prodromal'}));
data.GENETICS = double(categorical(data.GENETICS, {'PINK1','PARKIN','SRDC','HPSM',...
                                                   'RBD','LRRK2','SNCA','GBA'}));
data.ETHNICITY = double(categorical(data.ETHNICITY, {'ASIAN','BLACK','HAWOPI',...
                                                     'INDALS','UNKNOWN','WHITE', ...
                                                     'NOTSPECIFIED'}));
data.SEX = double(categorical(data.SEX, {'Female','Male'}));
data.HANDED = double(categorical(data.HANDED, {'Right','Left','Mixed'}));
data.NP4TOT = cellfun(@(x) ifelse_categorization(strcmp(x, 'NA'), NaN, str2double(x)), data.NP4TOT);

% Age at symptoms calculation
[BirthMonth, BirthYear] = strtok(data.BIRTHDT, '/');
BirthYear = cellfun(@(x) str2double(x(2:end)), BirthYear);
BirthMonth = cellfun(@str2double, BirthMonth);
[SXMonth, SXYear] = strtok(data.SXDT, '/');
SXYear = cellfun(@(x) str2double(x(2:end)), SXYear);
SXMonth = cellfun(@str2double, SXMonth);
data.SXAGE = (SXYear - BirthYear) .* 12 + (SXMonth - BirthMonth);
data.SXAGE = data.SXAGE ./ 12;

% Clear
clear BirthMonth BirthYear SXMonth SXYear 
data = removevars(data,{'BIRTHDT', 'SXDT'});

disp('CATEGORIZED DATA')
head(data)

%% Patient selection

% Remove cohort SWEDD and Prodromal
swedd_mask = (data.COHORT == 3);
data(swedd_mask,:) = [];
prodr_mask = (data.COHORT == 4);
data(prodr_mask,:) = [];
clear swedd_mask prodr_mask

% Remove patients with uncomplete DATSCAN
not_completed_mask = (data.DATSCAN == 0);
data(not_completed_mask,:) = [];
data = removevars(data, {'DATSCAN'});
clear not_completed_mask

% Remove patients with low quality DATSCAN
low_quality_mask = (data.DATSCAN_QUALITY_RATING == 3);
data(low_quality_mask,:) = [];
data = removevars(data, {'DATSCAN_QUALITY_RATING'});
clear low_quality_mask

% Remove left handed patients
left_handed_mask = (data.HANDED == 2);
data(left_handed_mask,:) = [];
mixed_handed_mask = (data.HANDED == 3);
data(mixed_handed_mask,:) = [];
data = removevars(data, {'HANDED'});
clear left_handed_mask mixed_handed_mask

disp('DATASET after patient selection')
head(data)

%% PD - HC splitting

P_idx = (data.COHORT==1);
data_p = data(P_idx, :);
data_p = removevars(data_p,{'COHORT'});

H_idx = (data.COHORT==2);
data_h = data(H_idx, :);
data_h = removevars(data_h,{'COHORT','GENETICS', 'SXAGE'});

% Clear
clear H_idx P_idx

%% Variable NaN percentage

% PD patients
nanCount_p = zeros([size(data_p,2),1]);
for i=1:size(data_p,2)
    tmp = data_p{:,i};
    nanCount_p(i) = sum(isnan(tmp),1)/size(data_p,1);
end
nanCount_p = nanCount_p.*100;

% % nanCount on Variables for PD (%)
% figure
% stem(1:length(nanCount_p), nanCount_p)
% xticks(1:size(data_p,2));
% xticklabels(data_p.Properties.VariableNames);
% title('nanCount on Variables for PD (%)')

% HC patients
nanCount_h = zeros([size(data_h,2),1]);
for i=1:size(data_h,2)
    tmp = data_h{:,i};
    nanCount_h(i) = sum(isnan(tmp),1)/size(data_h,1);
end
nanCount_h = nanCount_h.*100;

% % nanCount on Variables for HC (%)
% figure
% stem(1:length(nanCount_h), nanCount_h)
% xticks(1:size(data_p,2));
% xticklabels(data_h.Properties.VariableNames);
% title('nanCount on Variables for HC (%)')

% Clear
data_p = removevars(data_p, {'NP4TOT','LEDD'});
data_h = removevars(data_h, {'NP4TOT','LEDD'});

clear i nanCount_p nanCount_h tmp

%% Patients NaN percentage

% For PD patients
nanCount_p = sum(ismissing(data_p), 2);  % Count the number of NaNs in each row

% nanCount on Patients for PD (#)
figure;
stem(1:length(nanCount_p), nanCount_p);
title('nanCount on Patients for PD (#)');

% Eliminate rows with more than 2 NaN values
data_p(nanCount_p > 0, :) = [];

% For HC patients
nanCount_h = sum(ismissing(data_h), 2);  % Count the number of NaNs in each row

% nanCount on Patients for HC (#)
figure;
stem(1:length(nanCount_h), nanCount_h);
title('nanCount on Patients for HC (#)');

% Eliminate rows with more than 2 NaN values
data_h(nanCount_h > 0, :) = [];

% Clear
clear nanCount_p nanCount_h 

%% Outliers removal

% Outliers calculation based on single areas
datscan_h = data_h(:,{'DATSCAN_CAUDATE_L','DATSCAN_CAUDATE_R','DATSCAN_PUTAMEN_L','DATSCAN_PUTAMEN_R'});
datscan_p = data_p(:,{'DATSCAN_CAUDATE_L','DATSCAN_CAUDATE_R','DATSCAN_PUTAMEN_L','DATSCAN_PUTAMEN_R'});

idx_h_outliers = outlier_rows(datscan_h);
idx_p_outliers = outlier_rows(datscan_p);

% Outliers removal
data_h = data_h(~idx_h_outliers,:);
data_p = data_p(~idx_p_outliers,:);

% % Boxplots of DATSCAN distributions
% figure
% subplot(211)
% boxplot(table2array(datscan_h))
% title('HC DATSCAN distribution')
% % xticklabels(datscan_h.Properties.VariableNames);
% subplot(212)
% boxplot(table2array(datscan_p))
% title('PD DATSCAN distribution')
% % xticklabels(datscan_p.Properties.VariableNames);

% Clear
clear idx_p_outliers idx_h_outliers

%% Lateralization Index

data_p.LAT_IDX_CAUDATE = (data_p.DATSCAN_CAUDATE_R - data_p.DATSCAN_CAUDATE_L)...
                        ./(data_p.DATSCAN_CAUDATE_R + data_p.DATSCAN_CAUDATE_L);
data_h.LAT_IDX_CAUDATE = (data_h.DATSCAN_CAUDATE_R - data_h.DATSCAN_CAUDATE_L)...
                        ./(data_h.DATSCAN_CAUDATE_R + data_h.DATSCAN_CAUDATE_L);
data_p.LAT_IDX_PUTAMEN = (data_p.DATSCAN_PUTAMEN_R - data_p.DATSCAN_PUTAMEN_L)...
                        ./(data_p.DATSCAN_PUTAMEN_R + data_p.DATSCAN_PUTAMEN_L);
data_h.LAT_IDX_PUTAMEN = (data_h.DATSCAN_PUTAMEN_R - data_h.DATSCAN_PUTAMEN_L)...
                        ./(data_h.DATSCAN_PUTAMEN_R + data_h.DATSCAN_PUTAMEN_L);

datscan_h = data_h(:,{'DATSCAN_CAUDATE_L','DATSCAN_CAUDATE_R','DATSCAN_PUTAMEN_L', ...
                      'DATSCAN_PUTAMEN_R','LAT_IDX_CAUDATE','LAT_IDX_PUTAMEN'});
datscan_p = data_p(:,{'DATSCAN_CAUDATE_L','DATSCAN_CAUDATE_R','DATSCAN_PUTAMEN_L', ...
                      'DATSCAN_PUTAMEN_R','LAT_IDX_CAUDATE','LAT_IDX_PUTAMEN'});

writetable(datscan_h(:,{'LAT_IDX_CAUDATE', 'LAT_IDX_PUTAMEN'}), 'datscan_h.csv');
writetable(datscan_p(:,{'LAT_IDX_CAUDATE', 'LAT_IDX_PUTAMEN'}), 'datscan_p.csv');

%% Data exploration

covariates_p = {'ENROLL_AGE','GENETICS','ETHNICITY','SEX','EDUCYRS','HTCM','WGTKG','SXAGE'};
covariates_h = {'ETHNICITY','SEX','EDUCYRS','HTCM','WGTKG'};

% Descriptive statistics for PD patients
means = varfun(@(x) mean(x, 'omitnan'), data_p, 'InputVariables', covariates_p);
medians = varfun(@(x) median(x, 'omitnan'), data_p, 'InputVariables', covariates_p);
std_devs = varfun(@(x) std(x, 'omitnan'), data_p, 'InputVariables', covariates_p);
minimums = varfun(@(x) min(x, [], 'omitnan'), data_p, 'InputVariables', covariates_p);
maximums = varfun(@(x) max(x, [], 'omitnan'), data_p, 'InputVariables', covariates_p);
percentiles_25 = varfun(@(x) prctile(x, 25), data_p, 'InputVariables', covariates_p);
percentiles_75 = varfun(@(x) prctile(x, 75), data_p, 'InputVariables', covariates_p);

% disp('Mean:')
% disp(means)
% disp('Median:')
% disp(medians)
% disp('Standard Deviation:')
% disp(std_devs)
% disp('Minimums:')
% disp(minimums)
% disp('Maximums:')
% disp(maximums)
% disp('Percentile 25:')
% disp(percentiles_25)
% disp('Percentile 75:')
% disp(percentiles_75)

descriptiveStatistics_p = table(means,medians,std_devs,minimums,maximums,percentiles_25, ...
                                percentiles_75,'VariableNames',{'means','medians', ...
                                'std_devs','minimums','maximums','percentiles_25', ...
                                'percentiles_75'});

% Descriptive statistics for HC patients
means = varfun(@(x) mean(x, 'omitnan'), data_h, 'InputVariables', covariates_h);
medians = varfun(@(x) median(x, 'omitnan'), data_h, 'InputVariables', covariates_h);
std_devs = varfun(@(x) std(x, 'omitnan'), data_h, 'InputVariables', covariates_h);
minimums = varfun(@(x) min(x, [], 'omitnan'), data_h, 'InputVariables', covariates_h);
maximums = varfun(@(x) max(x, [], 'omitnan'), data_h, 'InputVariables', covariates_h);
percentiles_25 = varfun(@(x) prctile(x, 25), data_h, 'InputVariables', covariates_h);
percentiles_75 = varfun(@(x) prctile(x, 75), data_h, 'InputVariables', covariates_h);

% disp('Mean:')
% disp(means)
% disp('Median:')
% disp(medians)
% disp('Standard Deviation:')
% disp(std_devs)
% disp('Minimums:')
% disp(minimums)
% disp('Maximums:')
% disp(maximums)
% disp('Percentile 25:')
% disp(percentiles_25)
% disp('Percentile 75:')
% disp(percentiles_75)

descriptiveStatistics_h = table(means,medians,std_devs,minimums,maximums,percentiles_25, ...
                                percentiles_75,'VariableNames',{'means','medians', ...
                                'std_devs','minimums','maximums','percentiles_25', ...
                                'percentiles_75'});

% Clear
clear means medians std_devs minimums maximums percentiles_25 percentiles_75

%% Graphs
% 
% % Continuous Covariates
% continuousCovariates = {'ENROLL_AGE', 'EDUCYRS', 'HTCM', 'WGTKG','SXAGE'};
% 
% % % Continuous covariates for PD
% % for i = 1:length(continuousCovariates)
% %     figure;
% %
% %     % Histogram
% %     subplot(1, 2, 1);
% %     histogram(data_p.(continuousCovariates{i}), 'Normalization', 'probability', 'DisplayStyle', 'stairs');
% %     title(['PD: Histogram of ', continuousCovariates{i}]);
% %
% %     % Boxplot
% %     subplot(1, 2, 2);
% %     boxplot(data_p.(continuousCovariates{i}));
% %     title(['PD: Boxplot of ', continuousCovariates{i}]);
% % end
% 
% % % Continuous covariates for HC
% % for i = 1:length(continuousCovariates)
% %     figure;
% %
% %     % Histogram
% %     subplot(1, 2, 1);
% %     histogram(data_h.(continuousCovariates{i}), 'Normalization', 'probability', 'DisplayStyle', 'stairs');
% %     title(['HC: Histogram of ', continuousCovariates{i}]);
% %
% %     % Boxplot
% %     subplot(1, 2, 2);
% %     boxplot(data_h.(continuousCovariates{i}));
% %     title(['HC: Boxplot of ', continuousCovariates{i}]);
% % end
% 
% % Categorical covariates
% categoricalCovariates = {'GENETICS', 'ETHNICITY', 'SEX'};
% 
% % % Categorical covariates for PD
% % for i = 1:length(categoricalCovariates)
% %     figure;
% %
% %     % Conteggio di frequenza
% %     [counts, categories] = histcounts(categorical(data_p.(categoricalCovariates{i})));
% %
% %     bar(counts);
% %     set(gca, 'XTickLabel', categories);
% %     title(['PD: Occurrences of ', categoricalCovariates{i}]);
% %
% %     % disp(['PD: Occurrences of ', categoricalVariables{i} ':']);
% %     % disp(table(categories', counts', 'VariableNames', {'Category', 'Count'}));
% % end
% 
% % % Categorical covariates for HC
% % for i = 2:length(categoricalCovariates)
% %     figure;
% %
% %     % Conteggio di frequenza
% %     [counts, categories] = histcounts(categorical(data_h.(categoricalCovariates{i})));
% %
% %     bar(counts);
% %     set(gca, 'XTickLabel', categories);
% %     title(['HC: Occurrences of ', categoricalCovariates{i}]);
% %
% %     % disp(['HC: Occurrences of ', categoricalVariables{i} ':']);
% %     % disp(table(categories', counts', 'VariableNames', {'Category', 'Count'}));
% % end
% 
% % Clear
% clear categories counts continuousCovariates categoricalCovariates

%% Gaussianity test

% Caudate
[~, p_CL_HC] = lillietest(datscan_h.DATSCAN_CAUDATE_L);
[~, p_CR_HC] = lillietest(datscan_h.DATSCAN_CAUDATE_R);
[~, p_CLatIdx_HC] = lillietest(datscan_h.LAT_IDX_CAUDATE);
[~, p_CL_PD] = lillietest(datscan_p.DATSCAN_CAUDATE_L);
[~, p_CR_PD] = lillietest(datscan_p.DATSCAN_CAUDATE_R);
[~, p_CLatIdx_PD] = lillietest(datscan_p.LAT_IDX_CAUDATE);

fprintf('\n')
fprintf('Gaussianity Test Results for Healthy Controls (Left Caudate): p-value = %.4f\n', p_CL_HC);
fprintf('Gaussianity Test Results for Healthy Controls (Right Caudate): p-value = %.4f\n', p_CR_HC);
fprintf('Gaussianity Test Results for Healthy Controls (Lat. Idx Caudate): p-value = %.4f\n', p_CLatIdx_HC);
fprintf('Gaussianity Test Results for PD Patients (Left Caudate): p-value = %.4f\n', p_CL_PD);
fprintf('Gaussianity Test Results for PD Patients (Right Caudate): p-value = %.4f\n', p_CR_PD);
fprintf('Gaussianity Test Results for PD Patients (Lat. Idx Caudate): p-value = %.4f\n', p_CLatIdx_PD);

% Putamen
[~, p_PL_HC] = lillietest(datscan_h.DATSCAN_PUTAMEN_L);
[~, p_PR_HC] = lillietest(datscan_h.DATSCAN_PUTAMEN_R);
[~, p_PLatIdx_HC] = lillietest(datscan_h.LAT_IDX_PUTAMEN);
[~, p_PL_PD] = lillietest(datscan_p.DATSCAN_PUTAMEN_L);
[~, p_PR_PD] = lillietest(datscan_p.DATSCAN_PUTAMEN_R);
[~, p_PLatIdx_PD] = lillietest(datscan_p.LAT_IDX_PUTAMEN);

fprintf('\n')
fprintf('Gaussianity Test Results for Healthy Controls (Left Putamen): p-value = %.4f\n', p_PL_HC);
fprintf('Gaussianity Test Results for Healthy Controls (Left Putamen): p-value = %.4f\n', p_PR_HC);
fprintf('Gaussianity Test Results for Healthy Controls (Lat. Idx Putamen): p-value = %.4f\n', p_PLatIdx_HC);
fprintf('Gaussianity Test Results for PD Patients (Left Putamen): p-value = %.4f\n', p_PL_PD);
fprintf('Gaussianity Test Results for PD Patients (Right Putamen): p-value = %.4f\n', p_PR_PD);
fprintf('Gaussianity Test Results for PD Patients (Lat. Idx Putamen): p-value = %.4f\n', p_PLatIdx_PD);

% Clear
clear p_CL_HC p_CL_PD p_CLatIdx_HC p_CLatIdx_PD p_CR_HC p_CR_PD p_PL_HC p_PL_PD ...
      p_PLatIdx_HC p_PLatIdx_PD p_PR_HC p_PR_PD

%% Descriptive Statistics for DATSCAN variables

fprintf('\n')
fprintf('Descriptive Statistics for Healthy Controls\n');
fprintf('Mean (SD) Caudate L: %.2f (%.2f)\n', mean(datscan_h.DATSCAN_CAUDATE_L), std(datscan_h.DATSCAN_CAUDATE_L));
fprintf('Mean (SD) Caudate R: %.2f (%.2f)\n', mean(datscan_h.DATSCAN_CAUDATE_R), std(datscan_h.DATSCAN_CAUDATE_R));
fprintf('Mean (SD) Putamen L: %.2f (%.2f)\n', mean(datscan_h.DATSCAN_PUTAMEN_L), std(datscan_h.DATSCAN_PUTAMEN_L));
fprintf('Mean (SD) Putamen R: %.2f (%.2f)\n', mean(datscan_h.DATSCAN_PUTAMEN_R), std(datscan_h.DATSCAN_PUTAMEN_R));
fprintf('Median (SD) Lat. Idx Caudate: %.2f (%.2f)\n', median(datscan_h.LAT_IDX_CAUDATE), std(datscan_h.LAT_IDX_CAUDATE));
fprintf('Median (SD) Lat. Idx Putamen: %.2f (%.2f)\n', median(datscan_h.LAT_IDX_PUTAMEN), std(datscan_h.LAT_IDX_PUTAMEN));

fprintf('\n')
fprintf('Descriptive Statistics for PD Patients\n');
fprintf('Mean (SD) Caudate L: %.2f (%.2f)\n', mean(datscan_p.DATSCAN_CAUDATE_L), std(datscan_p.DATSCAN_CAUDATE_L));
fprintf('Mean (SD) Caudate R: %.2f (%.2f)\n', mean(datscan_p.DATSCAN_CAUDATE_R), std(datscan_p.DATSCAN_CAUDATE_R));
fprintf('Mean (SD) Putamen L: %.2f (%.2f)\n', mean(datscan_p.DATSCAN_PUTAMEN_L), std(datscan_p.DATSCAN_PUTAMEN_L));
fprintf('Mean (SD) Putamen R: %.2f (%.2f)\n', mean(datscan_p.DATSCAN_PUTAMEN_R), std(datscan_p.DATSCAN_PUTAMEN_R));
fprintf('Median (SD) Lat. Idx Caudate: %.2f (%.2f)\n', median(datscan_p.LAT_IDX_CAUDATE), std(datscan_p.LAT_IDX_CAUDATE));
fprintf('Median (SD) Lat. Idx Putamen: %.2f (%.2f)\n', median(datscan_p.LAT_IDX_PUTAMEN), std(datscan_p.LAT_IDX_PUTAMEN));

%% SignTest
zero_median_pvalue = [];
s = [];

s(1) = skewness((datscan_h.LAT_IDX_CAUDATE));
s(2) = skewness((datscan_h.LAT_IDX_PUTAMEN));
s(3) = skewness((datscan_p.LAT_IDX_CAUDATE));
s(4) = skewness((datscan_p.LAT_IDX_PUTAMEN));

zero_median_pvalue(1) = signrank((datscan_h.LAT_IDX_CAUDATE));
zero_median_pvalue(2) = signrank((datscan_h.LAT_IDX_PUTAMEN));
zero_median_pvalue(3) = signrank((datscan_p.LAT_IDX_CAUDATE));
zero_median_pvalue(4) = signrank((datscan_p.LAT_IDX_PUTAMEN));


%% Correlation Matrix

% HC
corr_var_h = {'ENROLL_AGE','ETHNICITY','SEX','EDUCYRS','HTCM','WGTKG','LAT_IDX_CAUDATE',...
              'LAT_IDX_PUTAMEN'};
correlationMatrix_h = corr(data_h{:,corr_var_h}, 'rows', 'pairwise');

% % Heatmap
% figure
% heatmap(correlationMatrix_h, 'XData', corr_var_h, 'YData',  corr_var_h, 'Colormap', jet, 'ColorLimits', [-1, 1]);
% title('Correlation Map for Healthy controls')

% PD
corr_var_p = {'ENROLL_AGE','GENETICS','ETHNICITY','SEX','EDUCYRS','HTCM','WGTKG',...
              'SXAGE','LAT_IDX_CAUDATE','LAT_IDX_PUTAMEN'};
correlationMatrix_p = corr(data_p{:,corr_var_p}, 'rows', 'pairwise');

% % Heatmap
figure
heatmap(correlationMatrix_p, 'XData',  corr_var_p, 'YData',  corr_var_p, 'Colormap', jet, 'ColorLimits', [-1, 1]);
title('Correlation Map for Parkinson disease')

%% VIF

% PD
inv_corr_matrix_p = inv(correlationMatrix_p);
for i = 1:size(corr_var_p,2)
    values(i) = inv_corr_matrix_p(i,i);
end
VIF_p = table(corr_var_p',values');
VIF_p = renamevars(VIF_p,{'Var1','Var2'},{'Vars','VIF'});
low_VIF_mask = (VIF_p.VIF<10);
VIF_p(low_VIF_mask,:) = [];

disp('High VIF variables for PD: ')
if size(VIF_p,1)>0
    disp(VIF_p)
else
    disp('None')
end

% Clear
clear inv_corr_matrix_p vars values low_VIF_mask

% HC
inv_corr_matrix_h = inv(correlationMatrix_h);
vars = data_h.Properties.VariableNames;
for i = 1:size(corr_var_h,2)
    values(i) = inv_corr_matrix_h(i,i);
end
VIF_h = table(corr_var_h',values');
VIF_h = renamevars(VIF_h,{'Var1','Var2'},{'Vars','VIF'});
low_VIF_mask = (VIF_h.VIF<10);
VIF_h(low_VIF_mask,:) = [];

disp('High VIF variables for HC: ')
if size(VIF_h,1)>0
    disp(VIF_h)
else
    disp('None')
end

% Clear
clear inv_corr_matrix_h vars values low_VIF_mask corr_var_h corr_var_p VIF_h VIF_p 

% Clear
covariates_p(end) = [];

%% Regression

% Criteria
threshold_FPval = 0.05;

%% Covariates regression for HC
var_datscan = datscan_h(:,{'LAT_IDX_CAUDATE','LAT_IDX_PUTAMEN'}).Properties.VariableNames;
tbl_cov = data_h(:,covariates_h);

% HC linear fitting of covariates
results = struct('DependentVar',{},'Covariates',{},'R2',{},'pvalues',{},'Fstat',{}, ...
                 'FPval',{},'AIC',{},'BIC',{});

for i = 1:length(var_datscan)
    dep_var = var_datscan{i};

    for k = 1:min(length(covariates_h),3)
        combos = nchoosek(1:length(covariates_h), k);

        for j = 1:size(combos, 1)
            cov_idx = combos(j, :);
            cov_vars = covariates_h(cov_idx);

            formula = strcat(dep_var, ' ~ ', strjoin(cov_vars, ' + '));
            data_fit = [datscan_h(:, dep_var), tbl_cov(:, cov_vars)];

            mdl = fitlm(data_fit, formula);

            % % Results visualization
            % disp(mdl);
            % 
            % % R²
            % R2 = mdl.Rsquared.Ordinary;
            % disp(['R²: ', num2str(R2)]);
            % 
            % % p-values of coefficients
            % p_values = mdl.Coefficients.pValue;
            % disp('Coefficients p_values:');
            % disp(p_values);
            % 
            % % % Standardized residuals graph
            % % figure;
            % % plotResiduals(mdl, 'probability');
            % % title(['ProbPlot of Res for ' formula]);
            % 
            % % F-statistic and p-value
            % F_stat = mdl.ModelFitVsNullModel.Fstat;
            % p_value_F = mdl.ModelFitVsNullModel.Pvalue;
            % disp(['F-statistic: ', num2str(F_stat), ', p-value: ', num2str(p_value_F)]);

            results(end+1).DependentVar = dep_var;
            results(end).Covariates = cov_vars;
            results(end).R2 = mdl.Rsquared.Ordinary;
            results(end).pvalues = mdl.Coefficients.pValue;
            results(end).Fstat = mdl.ModelFitVsNullModel.Fstat;
            results(end).FPval = mdl.ModelFitVsNullModel.Pvalue;
            results(end).AIC = mdl.ModelCriterion.AIC;
            results(end).BIC = mdl.ModelCriterion.BIC;
        end
    end
end
covariatesLinFit_h = struct2table(results);

selected_h_lin = covariatesLinFit_h(covariatesLinFit_h.FPval < threshold_FPval, :);
fit_coeff_mask = ones(size(selected_h_lin,1),1);
for i = 1:size(selected_h_lin,1)
    fit_coeff = cell2mat(selected_h_lin{i,'pvalues'});
    if any(fit_coeff>=0.06)
        fit_coeff_mask(i) = 0;
    end
end
idx = find(fit_coeff_mask == 0);
selected_h_lin(idx,:) = [];
clear fit_coeff_mask

% HC non-linear fitting of covariates
results = struct('DependentVar',{},'Covariates',{},'R2',{},'pvalues',{},'Fstat',{}, ...
                 'FPval',{},'AIC',{},'BIC',{});

for i = 1:length(var_datscan)
    dep_var = var_datscan{i};

    for k = 1:min(length(covariates_h), 3)
        combos = nchoosek(1:length(covariates_h), k);

        for j = 1:size(combos, 1)
            cov_idx = combos(j, :);
            cov_vars = covariates_h(cov_idx);

            % Estrazione delle variabili dipendenti e indipendenti
            y = datscan_h.(dep_var);
            X = tbl_cov{:, cov_vars};

            % Definizione dinamica della funzione di modello non lineare
            if k == 1
                modelFun = @(b, x) b(1) + b(2)*x(:, 1) + b(3)*x(:, 1).^2;
                beta0 = [1; 1; 1];
            elseif k == 2
                modelFun = @(b, x) b(1) + b(2)*x(:, 1) + b(3)*x(:, 1).^2 + b(4)*x(:, 2) + b(5)*x(:, 2).^2;
                beta0 = [1; 1; 1; 1; 1];
            elseif k == 3
                modelFun = @(b, x) b(1) + b(2)*x(:, 1) + b(3)*x(:, 1).^2 + b(4)*x(:, 2) + b(5)*x(:, 2).^2 + b(6)*x(:, 3) + b(7)*x(:, 3).^2;
                beta0 = [1; 1; 1; 1; 1; 1; 1];
            end

            % Adattamento del modello non lineare
            [beta, R, J, CovB, MSE] = nlinfit(X, y, modelFun, beta0);

            % Calcolo delle previsioni e della statistica R²
            yfit = modelFun(beta, X);
            SSresid = sum((y - yfit).^2);
            SStotal = (length(y) - 1) * var(y);
            R2 = 1 - SSresid / SStotal;

            % p-values (approx)
            se = sqrt(diag(CovB));
            tstat = beta ./ se;
            p_values = 2 * (1 - tcdf(abs(tstat), length(y) - length(beta)));

            % Calcolo della F-statistica e del p-value associato (approx)
            F_stat = ((SStotal - SSresid) / (length(beta) - 1)) / (SSresid / (length(y) - length(beta)));
            p_value_F = 1 - fcdf(F_stat, length(beta) - 1, length(y) - length(beta));

            % Calcolo di AIC e BIC
            AIC = length(y) * log(SSresid / length(y)) + 2 * length(beta);
            BIC = length(y) * log(SSresid / length(y)) + log(length(y)) * length(beta);

            results(end+1).DependentVar = dep_var;
            results(end).Covariates = cov_vars;
            results(end).R2 = R2;
            results(end).pvalues = p_values;
            results(end).Fstat = F_stat;
            results(end).FPval = p_value_F;
            results(end).AIC = AIC;
            results(end).BIC = BIC;
        end
    end
end
covariatesNLinFit_h = struct2table(results);

selected_h_nlin = covariatesNLinFit_h(covariatesNLinFit_h.FPval < threshold_FPval, :);
fit_coeff_mask = ones(size(selected_h_nlin,1),1);
for i = 1:size(selected_h_nlin,1)
    fit_coeff = cell2mat(selected_h_nlin{i,'pvalues'});
    if any(fit_coeff>=0.06)
        fit_coeff_mask(i) = 0;
    end
end
idx = find(fit_coeff_mask == 0);
selected_h_nlin(idx,:) = [];
clear fit_coeff_mask

%% Covariates regression for PD
var_datscan = datscan_p(:,{'LAT_IDX_CAUDATE','LAT_IDX_PUTAMEN'}).Properties.VariableNames;
tbl_cov = data_p(:,covariates_p);

% PD linear fitting of covariates
results = struct('DependentVar',{},'Covariates',{},'R2',{},'pvalues',{},'Fstat',{}, ...
                 'FPval',{},'AIC',{},'BIC',{});

for i = 1:length(var_datscan)
    dep_var = var_datscan{i};

    for k = 1:min(length(covariates_p),3)
        combos = nchoosek(1:length(covariates_p), k);

        for j = 1:size(combos, 1)
            cov_idx = combos(j, :);
            cov_vars = covariates_p(cov_idx);

            formula = strcat(dep_var, ' ~ ', strjoin(cov_vars, ' + '));
            data_fit = [datscan_p(:, dep_var), tbl_cov(:, cov_vars)];

            mdl = fitlm(data_fit, formula);

            % % Results visualization
            % disp(mdl);
            % 
            % % R²
            % R2 = mdl.Rsquared.Ordinary;
            % disp(['R²: ', num2str(R2)]);
            % 
            % % p-values of coefficients
            % p_values = mdl.Coefficients.pValue;
            % disp('Coefficients p_values:');
            % disp(p_values);
            % 
            % % % Standardized residuals graph
            % % figure;
            % % plotResiduals(mdl, 'probability');
            % % title(['ProbPlot of Res for ' formula]);
            % 
            % % F-statistic and p-value
            % F_stat = mdl.ModelFitVsNullModel.Fstat;
            % p_value_F = mdl.ModelFitVsNullModel.Pvalue;
            % disp(['F-statistic: ', num2str(F_stat), ', p-value: ', num2str(p_value_F)]);

            results(end+1).DependentVar = dep_var;
            results(end).Covariates = cov_vars;
            results(end).R2 = mdl.Rsquared.Ordinary;
            results(end).pvalues = mdl.Coefficients.pValue;
            results(end).Fstat = mdl.ModelFitVsNullModel.Fstat;
            results(end).FPval = mdl.ModelFitVsNullModel.Pvalue;
            results(end).AIC = mdl.ModelCriterion.AIC;
            results(end).BIC = mdl.ModelCriterion.BIC;
        end
    end
end
covariatesLinFit_p = struct2table(results);

selected_p_lin = covariatesLinFit_p(covariatesLinFit_p.FPval < threshold_FPval, :);
fit_coeff_mask = ones(size(selected_p_lin,1),1);
for i = 1:size(selected_p_lin,1)
    fit_coeff = cell2mat(selected_p_lin{i,'pvalues'});
    if any(fit_coeff>=0.06)
        fit_coeff_mask(i) = 0;
    end
end
idx = find(fit_coeff_mask == 0);
selected_p_lin(idx,:) = [];
clear fit_coeff_mask

% PD non-linear fitting of covariates
results = struct('DependentVar',{},'Covariates',{},'R2',{},'pvalues',{},'Fstat',{}, ...
                 'FPval',{},'AIC',{},'BIC',{});

for i = 1:length(var_datscan)
    dep_var = var_datscan{i};

    for k = 1:min(length(covariates_p), 3)
        combos = nchoosek(1:length(covariates_p), k);

        for j = 1:size(combos, 1)
            cov_idx = combos(j, :);
            cov_vars = covariates_p(cov_idx);

            % Estrazione delle variabili dipendenti e indipendenti
            y = datscan_p.(dep_var);
            X = tbl_cov{:, cov_vars};

            % Definizione dinamica della funzione di modello non lineare
            if k == 1
                modelFun = @(b, x) b(1) + b(2)*x(:, 1) + b(3)*x(:, 1).^2;
                beta0 = [1; 1; 1];
            elseif k == 2
                modelFun = @(b, x) b(1) + b(2)*x(:, 1) + b(3)*x(:, 1).^2 + b(4)*x(:, 2) + b(5)*x(:, 2).^2;
                beta0 = [1; 1; 1; 1; 1];
            elseif k == 3
                modelFun = @(b, x) b(1) + b(2)*x(:, 1) + b(3)*x(:, 1).^2 + b(4)*x(:, 2) + b(5)*x(:, 2).^2 + b(6)*x(:, 3) + b(7)*x(:, 3).^2;
                beta0 = [1; 1; 1; 1; 1; 1; 1];
            end

            % Adattamento del modello non lineare
            [beta, R, J, CovB, MSE] = nlinfit(X, y, modelFun, beta0);

            % Calcolo delle previsioni e della statistica R²
            yfit = modelFun(beta, X);
            SSresid = sum((y - yfit).^2);
            SStotal = (length(y) - 1) * var(y);
            R2 = 1 - SSresid / SStotal;

            % p-values (approx)
            se = sqrt(diag(CovB));
            tstat = beta ./ se;
            p_values = 2 * (1 - tcdf(abs(tstat), length(y) - length(beta)));

            % Calcolo della F-statistica e del p-value associato (approx)
            F_stat = ((SStotal - SSresid) / (length(beta) - 1)) / (SSresid / (length(y) - length(beta)));
            p_value_F = 1 - fcdf(F_stat, length(beta) - 1, length(y) - length(beta));

            % Calcolo di AIC e BIC
            AIC = length(y) * log(SSresid / length(y)) + 2 * length(beta);
            BIC = length(y) * log(SSresid / length(y)) + log(length(y)) * length(beta);

            results(end+1).DependentVar = dep_var;
            results(end).Covariates = cov_vars;
            results(end).R2 = R2;
            results(end).pvalues = p_values;
            results(end).Fstat = F_stat;
            results(end).FPval = p_value_F;
            results(end).AIC = AIC;
            results(end).BIC = BIC;
        end
    end
end
covariatesNLinFit_p = struct2table(results);

selected_p_nlin = covariatesNLinFit_p(covariatesNLinFit_p.FPval < threshold_FPval, :);
fit_coeff_mask = ones(size(selected_p_nlin,1),1);
for i = 1:size(selected_p_nlin,1)
    fit_coeff = cell2mat(selected_p_nlin{i,'pvalues'});
    if any(fit_coeff>=0.06)
        fit_coeff_mask(i) = 0;
    end
end
idx = find(fit_coeff_mask == 0);
selected_p_nlin(idx,:) = [];
clear fit_coeff_mask

%% Symptoms regression for PD 
symptoms = {'NP1PTOT','NP1RTOT','NP2PTOT','NP3TOT','MCATOT','NHY'};
var_datscan = datscan_p(:,{'LAT_IDX_CAUDATE','LAT_IDX_PUTAMEN'}).Properties.VariableNames;
tbl_cov = data_p(:,symptoms);

% PD linear fitting of covariates
results = struct('DependentVar',{},'Covariates',{},'R2',{},'pvalues',{},'Fstat',{}, ...
                 'FPval',{},'AIC',{},'BIC',{});

for i = 1:length(var_datscan)
    dep_var = var_datscan{i};

    for k = 1:min(length(symptoms),3)
        combos = nchoosek(1:length(symptoms), k);

        for j = 1:size(combos, 1)
            cov_idx = combos(j, :);
            cov_vars = symptoms(cov_idx);

            formula = strcat(dep_var, ' ~ ', strjoin(cov_vars, ' + '));
            data_fit = [datscan_p(:, dep_var), tbl_cov(:, cov_vars)];

            mdl = fitlm(data_fit, formula);

            % % Results visualization
            % disp(mdl);
            % 
            % % R²
            % R2 = mdl.Rsquared.Ordinary;
            % disp(['R²: ', num2str(R2)]);
            % 
            % % p-values of coefficients
            % p_values = mdl.Coefficients.pValue;
            % disp('Coefficients p_values:');
            % disp(p_values);
            % 
            % % % Standardized residuals graph
            % % figure;
            % % plotResiduals(mdl, 'probability');
            % % title(['ProbPlot of Res for ' formula]);
            % 
            % % F-statistic and p-value
            % F_stat = mdl.ModelFitVsNullModel.Fstat;
            % p_value_F = mdl.ModelFitVsNullModel.Pvalue;
            % disp(['F-statistic: ', num2str(F_stat), ', p-value: ', num2str(p_value_F)]);

            results(end+1).DependentVar = dep_var;
            results(end).Covariates = cov_vars;
            results(end).R2 = mdl.Rsquared.Ordinary;
            results(end).pvalues = mdl.Coefficients.pValue;
            results(end).Fstat = mdl.ModelFitVsNullModel.Fstat;
            results(end).FPval = mdl.ModelFitVsNullModel.Pvalue;
            results(end).AIC = mdl.ModelCriterion.AIC;
            results(end).BIC = mdl.ModelCriterion.BIC;
        end
    end
end
symptomsLinFit_p = struct2table(results);

selected_symp_lin = symptomsLinFit_p(symptomsLinFit_p.FPval < threshold_FPval, :);
fit_coeff_mask = ones(size(selected_symp_lin,1),1);
for i = 1:size(selected_symp_lin,1)
    fit_coeff = cell2mat(selected_symp_lin{i,'pvalues'});
    if any(fit_coeff>=0.06)
        fit_coeff_mask(i) = 0;
    end
end
idx = find(fit_coeff_mask == 0);
selected_symp_lin(idx,:) = [];
clear fit_coeff_mask

% PD non-linear fitting of covariates
results = struct('DependentVar',{},'Covariates',{},'R2',{},'pvalues',{},'Fstat',{}, ...
                 'FPval',{},'AIC',{},'BIC',{});

for i = 1:length(var_datscan)
    dep_var = var_datscan{i};

    for k = 1:min(length(symptoms), 3)
        combos = nchoosek(1:length(symptoms), k);

        for j = 1:size(combos, 1)
            cov_idx = combos(j, :);
            cov_vars = symptoms(cov_idx);

            % Estrazione delle variabili dipendenti e indipendenti
            y = datscan_p.(dep_var);
            X = tbl_cov{:, cov_vars};

            % Definizione dinamica della funzione di modello non lineare
            if k == 1
                modelFun = @(b, x) b(1) + b(2)*x(:, 1) + b(3)*x(:, 1).^2;
                beta0 = [1; 1; 1];
            elseif k == 2
                modelFun = @(b, x) b(1) + b(2)*x(:, 1) + b(3)*x(:, 1).^2 + b(4)*x(:, 2) + b(5)*x(:, 2).^2;
                beta0 = [1; 1; 1; 1; 1];
            elseif k == 3
                modelFun = @(b, x) b(1) + b(2)*x(:, 1) + b(3)*x(:, 1).^2 + b(4)*x(:, 2) + b(5)*x(:, 2).^2 + b(6)*x(:, 3) + b(7)*x(:, 3).^2;
                beta0 = [1; 1; 1; 1; 1; 1; 1];
            end

            % Adattamento del modello non lineare
            [beta, R, J, CovB, MSE] = nlinfit(X, y, modelFun, beta0);

            % Calcolo delle previsioni e della statistica R²
            yfit = modelFun(beta, X);
            SSresid = sum((y - yfit).^2);
            SStotal = (length(y) - 1) * var(y);
            R2 = 1 - SSresid / SStotal;

            % p-values (approx)
            se = sqrt(diag(CovB));
            tstat = beta ./ se;
            p_values = 2 * (1 - tcdf(abs(tstat), length(y) - length(beta)));

            % Calcolo della F-statistica e del p-value associato (approx)
            F_stat = ((SStotal - SSresid) / (length(beta) - 1)) / (SSresid / (length(y) - length(beta)));
            p_value_F = 1 - fcdf(F_stat, length(beta) - 1, length(y) - length(beta));

            % Calcolo di AIC e BIC
            AIC = length(y) * log(SSresid / length(y)) + 2 * length(beta);
            BIC = length(y) * log(SSresid / length(y)) + log(length(y)) * length(beta);

            results(end+1).DependentVar = dep_var;
            results(end).Covariates = cov_vars;
            results(end).R2 = R2;
            results(end).pvalues = p_values;
            results(end).Fstat = F_stat;
            results(end).FPval = p_value_F;
            results(end).AIC = AIC;
            results(end).BIC = BIC;
        end
    end
end
symptomsNLinFit_p = struct2table(results);

selected_symp_nlin = symptomsNLinFit_p(symptomsNLinFit_p.FPval < threshold_FPval, :);
fit_coeff_mask = ones(size(selected_symp_lin,1),1);
for i = 1:size(selected_symp_nlin,1)
    fit_coeff = cell2mat(selected_symp_nlin{i,'pvalues'});
    if any(fit_coeff>=0.06)
        fit_coeff_mask(i) = 0;
    end
end
idx = find(fit_coeff_mask == 0);
selected_symp_nlin(idx,:) = [];
clear fit_coeff_mask idx AIC beta beta0 BIC combos correlationMatrix_h correlationMatrix_p cov_idx cov_vars covariates_h covariates_p covariatesLinFit_h ...
covariatesLinFit_p covariatesNLinFit_p covariatesNLinFit_h CovB data data_fit 

%% ML model

testRatio = 0.2;

% Aggiungi la variabile COHORT
data_h_ml = addvars(data_h, zeros(size(data_h,1),1), 'NewVariableNames', 'COHORT');
data_p_ml = addvars(data_p, ones(size(data_p,1),1), 'NewVariableNames', 'COHORT');

% Seleziona le variabili desiderate
selectedVariables = {'LAT_IDX_PUTAMEN', 'COHORT'};
inputVariables = {'LAT_IDX_PUTAMEN'};

% Sottoseleziona le variabili dai dati
data_h_ml = data_h_ml(:,selectedVariables);
data_p_ml = data_p_ml(1:170,selectedVariables);

% Unisci i dati
data_ml = vertcat(data_h_ml, data_p_ml);

% Mescola i dati
rng('default')
idx = randperm(height(data_ml));
shuffled_data_ml = data_ml(idx, :);

% Calcola le dimensioni dei set di train e test
numObservations = height(shuffled_data_ml);
numTest = round(numObservations * testRatio);
numTrain = numObservations - numTest;

% Dividi i dati in train e test set
train = shuffled_data_ml(1:numTrain, :);
test = shuffled_data_ml(numTrain+1:end, :);

% Estrai le variabili di input e target per il train e test set
X_train = table2array(train(:,inputVariables));
y_train = table2array(train(:, 'COHORT'));
X_test = table2array(test(:,inputVariables));
y_test = table2array(test(:, 'COHORT'));

%% training

% clear idx test train shuffled_data_ml
% 
cvp = cvpartition(y_train, 'KFold', 5);
% 
Mdl = fitcauto(X_train, y_train, ...
    'OptimizeHyperparameters', 'all', ...
    'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations', 50, ...
                                                 'CVPartition', cvp));

% Predizioni del modello sul set di test
[y_pred_test, scores] = predict(Mdl, X_test);

% Calcolo delle metriche di prestazione sul set di test
confMat_test = confusionmat(y_test, y_pred_test);
TP = confMat_test(2,2);
TN = confMat_test(1,1);
FP = confMat_test(1,2);
FN = confMat_test(2,1);

precision_test = TP/(TP+FP);
sensitivity_test = TP/(TP+FN);
specificity_test = TN/(TN+FP);

accuracy = (TP + TN) / (TP + TN + FP + FN);

disp(['Test Accuracy: ', num2str(accuracy)]);
disp(['Test Precision: ', num2str(precision_test)]);
disp(['Test Sensitivity: ', num2str(sensitivity_test)]);
disp(['Test Specificity: ', num2str(specificity_test)]);

%% 


% Calcolare la curva ROC
[X_roc, Y_roc, T_roc, AUC] = perfcurve(y_test, scores(:,2), 1);

J = Y_roc + (1 - X_roc) - 1;
% Trovare l'indice della soglia ottimale
[~, idx] = max(J);

% Determinare la threshold ottimale
threshold = T_roc(idx);

% Visualizzare la curva ROC e la threshold ottimale
figure;
plot(X_roc, Y_roc);
hold on;
plot(X_roc(idx), Y_roc(idx), 'ro');
text(X_roc(idx), Y_roc(idx), ['Threshold = ' num2str(threshold)], 'VerticalAlignment', 'bottom');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curve');
legend(['AUC = ' num2str(AUC)], 'Optimal Threshold', 'Location', 'southeast');
hold off;
