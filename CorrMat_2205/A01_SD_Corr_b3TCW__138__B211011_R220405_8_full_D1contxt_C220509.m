%% Contents of this cord
% 1. Create folder for storing analysis data
% 2. Data import
    % (Preprocess)
    % rules of variance (name)
% 3-1. TomPosi (TP, Tom+) 
% 3-2. TomNega (All, TN, Tom-)
    % 3-2-2. TomNega (first 30 cells = the same number of TomPosi)
% 3-3. All (=Tom Posi + Tom Nega)
        % (4-1. Create Shuffle data of Tom Posi (Shuffle type 1, raw
        % dF/across time and cells)
% 4-2. Create Shuffled data of Tom Posi (Shuffle type 2, raw dF across time, shuffle in each cells, randperm)
% 4-3. Create Shuffled data (Shuffle type 3, raster(binary) across time and cells)
        % (4-4 = 4-2. Create Shuffled data (Shuffle type 4, raster(binary)
        % across time, shuffle in each cells))
        % (4-5 (4-1 circshift).　Create Shuffled data (Shuffle type 5,
        % raster(binary) across time, shuffle in each cells))
% 4-6 (4-2 circshift).　Create Shuffled data (Shuffle type 5, raster(binary) across time, shuffle in each cells)
        % (4-7 (4-3 circshift).　Create Shuffled data (Shuffle type 5,
        % raster(binary) across time, shuffle in each cells))
        % (4-8 (4-4 circshift).　Create Shuffled data (Shuffle type 5,
        % raster(binary) across time, shuffle in each cells))
% 5-1. Random sampling of TomNega (random sample cells, the same number of TomPosi)

%% 1. Create folder for storing analysis data
% 1-1.フォルダ名の作成
%https://jp.mathworks.com/matlabcentral/answers/479642-
cd 'C:\Users\Inokuchi_Lab\Documents\MATLAB\Asai_';%folder作成場所の指定, https://jp.mathworks.com/help/matlab/ref/cd.html

mouseID = char('b3TCW__138__');%フォルダ名, データ保存名, https://jp.mathworks.com/help/matlab/ref/char.html
dinfo = char('R220405_8_full_D1context_');%フォルダ名
anltype = char('cellular_SD_');
ftime = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)

mkdir([mouseID,dinfo,anltype,ftime]) % フォルダの作成
cd ([mouseID,dinfo,anltype,ftime]);%今作成したフォルダに移動
pwd;%https://jp.mathworks.com/help/matlab/ref/pwd.html?searchHighlight=pwd&s_tid=srchtitle_pwd_1

%% 2. Data import

%csv(Excel) fileをMatlab右の空白(Workspace)にドラッグ > Numeric matrixで、time(1列目)とcell ID (1or2行目)を含んだ状態でimport
Data_imp = vD1Contextsort;%rename
% Data_raw = Data_imp(503:end,2:end);%remove time and Tom+/- labels & crop 12000 frames
Data_raw = Data_imp(1:end,2:end);%remove time and Tom+/- labels & crop 12000 frames

%% Preprocess
% Binning%CorrMatの中にbining/skipを作成

%% rules of variance (name)
%Dataはmatrix(time, cells)に使う
%Tom +(Posi) = TP, Tom -(Nega) = TN
%なるべくfunctionの名前(mean,sum,..)で、操作順に

%% 3-1. TomPosi (TP, Tom+)

% set parameter
SD_thres = 3;%binaryの時のSDのthreshold
% data setup
Num_TomPosi = 30;%manualで入力, データごとに値を変える必要あり
dFData_TP = Data_raw(:,1:Num_TomPosi);%extract Tom+ cells/remove Tom- cells
% basic infomation (mean/SD) of this data > binarization based on SD
% threshpold
df_TP_mean = mean(dFData_TP);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
df_TP_sd= std(dFData_TP);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1

sdData_TP = (dFData_TP - df_TP_mean)./df_TP_sd;%各データの細胞ごとのSD
SDbinaryData_TP = sdData_TP>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1

% %     % save the binarized data for correlation matrix using SD-binarized data
% %     SDbinaryData_TP_Time = [Data_imp(503:end,1),SDbinaryData_TP];%add time
% %     %dlmwrite([mouseID,dinfo,'TP_3SDbinary_Time.txt'], SDbinaryData_TP_Time, 'delimiter','\t');%
% %     writematrix([mouseID,dinfo,'TP_3SDbinary_Time.txt'], SDbinaryData_TP_Time, 'delimiter','\t');%

% results from binarized df matrix
TP_binarySD_SyncCells= sum(SDbinaryData_TP, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP =size(dFData_TP, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
    Timeframe_dF_TP = size(dFData_TP,1);%Time(frame)数
TP_SyncCells_Ratio = TP_binarySD_SyncCells./NumCell_dF_TP;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

TP_binarySD_Max = max(SDbinaryData_TP);%活動していない細胞の有無確認 (複数sessionまたいだ動画でTakekawa_systemを使ったため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
TP_Cell_Totalactive = sum(SDbinaryData_TP, 1);%各細胞の活動量 (3SD以上のframe数)
TP_Cell_Totalactive_Mean = mean(TP_Cell_Totalactive);%(TP)全細胞の平均活動量

%Cal of CorrMat > sum, sum_sum > figure(imagesc)
%Correlation matrix from Alan's(original name = out_matrix.m)
%How to use: CorrMat_Alan(binned_datamatrix(time, neurons),shift_window(4))
%     shift_cormat = 4;%parameter setupに移動
% % shift_cormat = 19;%4;
% % TP_CorrMat_Alan(dFData_TP,shift_cormat);toc;%biningなし、biningしたデータを入れるなら使える
tic;TP_CorrMat_Alan_skip = CorrMat_Alan_skip(dFData_TP, 19,4);toc
tic;TP_CorrMat_Alan_bin = CorrMat_Alan_bin(dFData_TP, 5,4);toc
% % writematrix(TP_CorrMat_Alan, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
%Cal of CorrMat_Sum_sum
% % TP_CorrMat_Alan_sum = sum(TP_CorrMat_Alan);
% % TP_CorrMat_sum_sum = sum(TP_CorrMat_Alan_sum);

TP_CorrMat_Alan_skip_sum = sum(TP_CorrMat_Alan_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる,https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum%20nan&s_tid=srchtitle_sum%20nan_1
TP_CorrMat_Alan_skip_sum_sum = sum(TP_CorrMat_Alan_skip_sum);

TP_CorrMat_Alan_bin_sum = sum(TP_CorrMat_Alan_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TP_CorrMat_Alan_bin_sum_sum = sum(TP_CorrMat_Alan_bin_sum);


% Ave_List_TP_CorrMat_sum_sum = mean(TP_CorrMat_sum_sum);%Shuffleの名残
% SD_List_TP_CorrMat_sum_sum = std(TP_CorrMat_sum_sum);
% SEM_List_TP_CorrMat_sum_sum = SD_List_TP_CorrMat_sum_sum./sqrt(iter_TPshuf2);

% % writematrix(TP_CorrMat_sum_sum, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_SumSum.txt'], 'delimiter','\t');%

figure;
imagesc(TP_CorrMat_Alan_skip);
figure;
imagesc(TP_CorrMat_Alan_bin);
dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%一応もう一度時間をとる
saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');

%CorrMat(Asai ver.)
tic;
TP_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_TP,19,4);
toc
writematrix(TP_CorrMat_skip, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_skip.txt'], 'delimiter','\t');%
tic;
TP_CorrMat_bin = CorrMat_Alan_As220506_bin(dFData_TP,5,4);%
toc
writematrix(TP_CorrMat_bin, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_bin.txt'], 'delimiter','\t');%

TP_CorrMat_skip_sum = sum(TP_CorrMat_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TP_CorrMat_skip_sum_sum = sum(TP_CorrMat_skip_sum);

TP_CorrMat_bin_sum = sum(TP_CorrMat_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TP_CorrMat_bin_sum_sum = sum(TP_CorrMat_bin_sum);

figure;
imagesc(TP_CorrMat_skip);
colorbar;caxis([0 0.01]);
title("Tom+")
figure;
imagesc(TP_CorrMat_bin);
colorbar;caxis([0 0.01]);
title("Tom+")

% saveas(gcf, [mouseID,'_CorrMat_As',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_CorrMat_As',dtime], 'tif');



%% 3-2. TomNega (All, TN, Tom-)

% % set parameter
% SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
% data setup
dFData_TN = Data_raw(:, 31:end);%extract Tom- cells/remove Tom+ cells

% basic infomation (mean/SD) of this data > binarization based on SD
% threshpold
df_TN_mean = mean(dFData_TN);%各細胞のmean, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
df_TN_sd = std(dFData_TN);%各細胞のSD, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1

sdData_TN = (dFData_TN - df_TN_mean)./df_TN_sd;%各データの細胞ごとのSD
SDbinaryData_TN = sdData_TN>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1

% %     % save the binarized data for correlation matrix using SD-binarized data
% %     SDbinaryData_TN_Time = [Data_imp(503:end,1),SDbinaryData_TN];%add time
% %     %dlmwrite([mouseID,dinfo,'TN_SDbinary_Time.txt'], SDbinaryData_TN_Time, 'delimiter','\t');%
% %     writematrix([mouseID,dinfo,'TN_SDbinary_Time.txt'], SDbinaryData_TN_Time, 'delimiter','\t');%

%Export dFData_TN for correlation matrix()
dFData_TN_time = [Data_imp(503:end,1),dFData_TN];
writematrix(dFData_TN_time,[mouseID,dinfo,'TN_dF_time_12544fr.txt'],'Delimiter','\t');


% results from binarized df matrix
TN_binarySD_SyncCells = sum(SDbinaryData_TN, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TN = size(dFData_TN,2);%the number of cells
    Timeframe_dF_TN = size(dFData_TN,1);%Time(frame)数
TN_SyncCells_Ratio = TN_binarySD_SyncCells./NumCell_dF_TN;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
    
% TN_binarySD_Max = max(SDbinaryData_TN);%活動していない細胞の有無確認 (複数sessionまたいだ動画でTakekawa_systemを使ったため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
TN_Cell_Totalactive = sum(SDbinaryData_TN, 1);%各細胞の活動量 (3SD以上のframe数)
TN_Cell_Totalactive_Mean = mean(TN_Cell_Totalactive);%(TN)全細胞の平均活動量

%Cal of CorrMat > sum, sum_sum > figure(imagesc)
% % TP_CorrMat_Alan(dFData_TP,shift_cormat);toc;%biningなし、biningしたデータを入れるなら使える
tic;TN_CorrMat_Alan_skip = CorrMat_Alan_skip(dFData_TN, 19,4);toc
tic;TN_CorrMat_Alan_bin = CorrMat_Alan_bin(dFData_TN, 5,4);toc
% % writematrix(TP_CorrMat_Alan, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
%Cal of CorrMat_Sum_sum
% % TP_CorrMat_Alan_sum = sum(TP_CorrMat_Alan);
% % TP_CorrMat_sum_sum = sum(TP_CorrMat_Alan_sum);

TN_CorrMat_Alan_skip_sum = sum(TN_CorrMat_Alan_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TN_CorrMat_skip_sum_sum = sum(TN_CorrMat_Alan_skip_sum);

TN_CorrMat_Alan_bin_sum = sum(TN_CorrMat_Alan_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TN_CorrMat_bin_sum_sum = sum(TN_CorrMat_Alan_bin_sum);


% Ave_List_TP_CorrMat_sum_sum = mean(TP_CorrMat_sum_sum);%Shuffleの名残
% SD_List_TP_CorrMat_sum_sum = std(TP_CorrMat_sum_sum);
% SEM_List_TP_CorrMat_sum_sum = SD_List_TP_CorrMat_sum_sum./sqrt(iter_TPshuf2);

% % writematrix(TP_CorrMat_sum_sum, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_SumSum.txt'], 'delimiter','\t');%

figure;
imagesc(TN_CorrMat_Alan_skip);

figure;
imagesc(TN_CorrMat_Alan_bin);
dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%一応もう一度時間をとる
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');

%CorrMat(Asai ver.)
tic;
TN_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_TN,19,4);
toc
% writematrix(TN_CorrMat_skip, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_skip.txt'], 'delimiter','\t');%
tic;
TN_CorrMat_bin = CorrMat_Alan_As220506_bin(dFData_TN,5,4);%
toc
% writematrix(TN_CorrMat_bin, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_bin.txt'], 'delimiter','\t');%

TN_CorrMat_skip_sum = sum(TN_CorrMat_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TN_CorrMat_skip_sum_sum = sum(TN_CorrMat_skip_sum);

TN_CorrMat_bin_sum = sum(TN_CorrMat_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TN_CorrMat_bin_sum_sum = sum(TN_CorrMat_bin_sum);

figure;
imagesc(TN_CorrMat_skip);
colorbar;caxis([0 0.01]);
figure;
imagesc(TN_CorrMat_bin);
colorbar;caxis([0 0.01]);



%細胞の活動度がある程度似た状態で同期性を評価しないとアンフェア(Data_SumCellをある程度揃えたい)

% %% 3-2-2. TomNega (first 30 cells = the same number of TomPosi)
% %細胞の活動度がある程度似た状態で同期性を評価しないとアンフェア(Data_SumCellをある程度揃えたい)
% %randomにTom Negaからsamplingする
% 
% % % set parameter
% % SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
% % data setup
% dfData_TN30 = Data_raw(:, 31:60);%extract first 30 Tom- cells/remove Tom+ cells
% 
% % basic infomation (mean/SD) of this data > binarization based on SD
% % threshpold
% df_TN30_mean =  mean(dfData_TN30);%各細胞のmean, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
% df_TN30_sd = std(dfData_TN30);%各細胞のSD, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% 
% sdData_TN30 = (dfData_TN30 - df_TN30_mean)./df_TN30_sd;%各データの細胞ごとのSD
% SDbinaryData_TN30 = sdData_TN30>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% 
% % %     % save the binarized data for correlation matrix using SD-binarized data
% % %     SDbinaryData_TN30_Time = [Data_imp(503:end,1),SDbinaryData_TN30];%add time
% % %     %dlmwrite([mouseID,dinfo,'TN30_SDbinary_Time.txt'], SDbinaryData_TN30_Time, 'delimiter','\t');%
% % %     writematrix([mouseID,dinfo,'TN30_SDbinary_Time.txt'], SDbinaryData_TN30_Time, 'delimiter','\t');%
% 
% % results from binarized df matrix
% TN30_binarySD_SyncCells = sum(SDbinaryData_TN30, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%     NumCell_dF_TN30 = size(dfData_TN30,2);
%     Timeframe_dF_TN30 = size(dfData_TN30,1);%Time(frame)数
% TN30_SyndCells_Ratio = TN30_binarySD_SyncCells./NumCell_dF_TN30;
% 
% % TN30_binarySD_Max = max(SDbinaryData_TN30);%活動していない細胞の有無確認 (複数sessionまたいだ動画でTakekawa_systemを使ったため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
% TN30_Cell_Totalactive = sum(SDbinaryData_TN30, 1);%各細胞の活動量 (3SD以上のframe数)
% TN30_Cell_Totalactive_mean = mean(TN30_Cell_Totalactive);%(TN)全細胞の平均活動量


%% 3-3. All (=Tom Posi + Tom Nega)

% % set parameter
% SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
% data setup
dFData_All = Data_raw;

% basic infomation (mean/SD) of this data > binarization based on SD
% threshpold
df_All_mean = mean(dFData_All);%各細胞のmean, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
df_All_sd = std(dFData_All);%各細胞のSD, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1

sdData_All = (dFData_All - df_All_mean)./df_All_sd;%各データの細胞ごとのSD
SDbinaryData_All = sdData_All>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1

% %     % save the binarized data for correlation matrix using SD-binarized data
% %     SDbinaryData_All_Time = [Data_imp(503:end,1),SDbinaryData_All];%add time
% %     %dlmwrite([mouseID,dinfo,'TN_SDbinary_Time.txt'], SDbinaryData_TN_Time, 'delimiter','\t');%
% %     writematrix([mouseID,dinfo,'All_SDbinary_Time.txt'], SDbinaryData_All_Time, 'delimiter','\t');%

%Export dFData_All for correlation matrix()
% dFData_All_time = [Data_imp(503:end,1),dFData_All];
% writematrix(dFData_All_time,[mouseID,dinfo,'All_dF_time_12544fr.txt'],'Delimiter','\t');

% results from binarized df matrix
All_binarySD_SyncCells = sum(SDbinaryData_All, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_All = size(dFData_All,2);
    Timeframe_dF_All = size(dFData_All,1);%Time(frame)数
All_SyncCella_Ratio = All_binarySD_SyncCells./NumCell_dF_All;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

All_binarySD_Max = max(SDbinaryData_All);%活動していない細胞の有無確認 (複数sessionまたいだ動画でTakekawa_systemを使ったため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
All_Cell_Totalactive = sum(SDbinaryData_All, 1);%各細胞の活動量 (3SD以上のframe数)
All_Cell_Totalactive_Mean = mean(All_Cell_Totalactive);%全細胞の平均活動量


%Cal of CorrMat > sum, sum_sum > figure(imagesc)
% % TP_CorrMat_Alan(dFData_TP,shift_cormat);toc;%biningなし、biningしたデータを入れるなら使える
    %tic;All_CorrMat_Alan_skip = CorrMat_Alan_skip(dFData_All, 19,4);toc
    %tic;All_CorrMat_Alan_bin = CorrMat_Alan_bin(dFData_All, 5,4);toc
% % writematrix(TP_CorrMat_Alan, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
%Cal of CorrMat_Sum_sum
% % TP_CorrMat_Alan_sum = sum(TP_CorrMat_Alan);
% % TP_CorrMat_sum_sum = sum(TP_CorrMat_Alan_sum);

    %All_CorrMat_Alan_skip_sum = sum(All_CorrMat_Alan_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
    %All_CorrMat_skip_sum_sum = sum(All_CorrMat_Alan_skip_sum);

    %All_CorrMat_Alan_bin_sum = sum(All_CorrMat_Alan_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
    %All_CorrMat_bin_sum_sum = sum(All_CorrMat_Alan_bin_sum);


% Ave_List_TP_CorrMat_sum_sum = mean(TP_CorrMat_sum_sum);%Shuffleの名残
% SD_List_TP_CorrMat_sum_sum = std(TP_CorrMat_sum_sum);
% SEM_List_TP_CorrMat_sum_sum = SD_List_TP_CorrMat_sum_sum./sqrt(iter_TPshuf2);

% % writematrix(TP_CorrMat_sum_sum, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_SumSum.txt'], 'delimiter','\t');%

    %figure;
    %imagesc(All_CorrMat_Alan_skip);
    %figure;
    %imagesc(All_CorrMat_Alan_bin);
    %dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%一応もう一度時間をとる
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');

%CorrMat(Asai ver.)
tic;
All_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_All,19,4);%25hr計算にかかって、すべてNaN
toc
% writematrix(TN_CorrMat_skip, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_skip.txt'], 'delimiter','\t');%
tic;
All_CorrMat_bin = CorrMat_Alan_As220506_bin(dFData_All,5,4);%
toc
% writematrix(TN_CorrMat_bin, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_bin.txt'], 'delimiter','\t');%

All_CorrMat_skip_sum = sum(All_CorrMat_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
All_CorrMat_skip_sum_sum = sum(All_CorrMat_skip_sum);

All_CorrMat_bin_sum = sum(All_CorrMat_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
All_CorrMat_bin_sum_sum = sum(All_CorrMat_bin_sum);

figure;
imagesc(All_CorrMat_skip);
colorbar;caxis([0 0.01]);
title('All cells')
figure;
imagesc(All_CorrMat_bin);
colorbar;caxis([0 0.01]);
title('All cells')



% % % % %% 4-1. Create Shuffle data of Tom Posi (Shuffle type 1, raw dF/across time and cells)
% % % % % % %一応動くが、データを集めるには向いていない。shuffle data matrixを作成することに集中したやり方 > 後で、大幅変更必要
% % 
% % % Prepare concatenated data for shuffling
% % dfData_TP_preShuf = reshape(dFData_TP,[],1);%Shuffle用にデータを1次元にする
% % L_dfData_TP_preShuf = size(dfData_TP_preShuf, 1);%データの全長、Shuffleデータを生成するために
% % 
% % % parameter
% % iter_TPShuf1 = 10000;%NumCell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)
% % 
% % % Create random ID sets for Shuffle
% % TPShuffle1_IDs = zeros(L_dfData_TP_preShuf,iter_TPShuf1);
% % for i = 1:iter_TPShuf1
% %     TPShuffle1_IDs(:,i)= randperm(L_dfData_TP_preShuf,L_dfData_TP_preShuf);
% % end
% % 
% % % Create random Shuffled data sets (replace raw data according to random IDs)
% % dfData_TP_postShuf1 = zeros(L_dfData_TP_preShuf,iter_TPShuf1);
% % for i = 1:iter_TPShuf1
% %     dfData_TP_postShuf1(:,i) = dfData_TP_preShuf(TPShuffle1_IDs(:,i));
% % 
% % end
% % 
% % % save the shuffled data
% % ftime = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
% % mkdir([mouseID,dinfo,'_TPshuf1_',ftime]) % フォルダの作成
% % cd ([mouseID,dinfo,'_TPshuf1_',ftime]);%今作成したフォルダに移動
% % for i = 1:iter_TPShuf1
% %     dfData_TPShuf1 = reshape(dfData_TP_postShuf1(:,i),size(dFData_TP));
% %     dlmwrite([mouseID,dinfo,'_TPshuf1_',num2str(i),'.txt'], dfData_TPShuf1, 'delimiter','\t');%
% % end
% % cd ..\
% % 
% % % save the shuffled data for correlation matrix
% % mkdir([mouseID,dinfo,'_TPshuf1_Time',ftime]) % フォルダの作成
% % cd ([mouseID,dinfo,'_TPshuf1_Time',ftime]);%今作成したフォルダに移動
% % for i = 1:iter_TPShuf1
% %     Data_TPShuf1_Time = [Data_imp(503:end,1),dfData_TPShuf1];%add time
% %     dlmwrite([mouseID,dinfo,'_TPshuf1_',num2str(i),'_Time.txt'], Data_TPShuf1_Time, 'delimiter','\t');%
% % end
% % cd ..\
% % % %     
% % 
% % % このやり方だとデータを生成するだけならいいが、計算のためにデータを読み込まないといけないので、4-2/3のようにした方がいいだろう

%% 4-2. Create Shuffled data of Tom Posi (Shuffle type 2, raw dF across time, shuffle in each cells, randperm)
%細胞ごとにshuffleしているので、mean/SDは変わらない
% Preparation/setup >> Create random IDs for Shuffle > random data & save > cal > summarize >>
% save

% % % Folders for saving the shuffled data
% % ftime = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
% % mkdir([mouseID,dinfo,'_TPshuf2_',ftime]) % フォルダの作成
% % % % cd ([mouseID,dinfo,'_shuf2_',ftime]);%今作成したフォルダに移動
% % mkdir([mouseID,dinfo,'_TPshuf2_time',ftime]) % フォルダの作成

%set parameters
%for shuffle
iter_TPshuf2 = 5000;%1000:p=0.01を検出するためには最低でもこれくらい必要なのでは？ Num_Time_frames_TP/10;Num_Cell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)
% NumCell_dF_TP =size(dFData_TP, 2);
% for cal (for...end内)
% SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
%for correlation matrix
% shift_cormat = 4;


% % % data setup(shuffleしないversionの名残、消してOK、同じ内容を含んでいる)
% % Timeframe_dF_TPShuf2 = Timeframe_dF_TP;%=size(dFData_TP,1);%Time(frame)数
% % NumCell_dF_TPShuf2= NumCell_dF_TP;%=size(dFData_TP,2)

% prepare matrix for shuffled data and IDs
TPShuf2_IDs = zeros(Timeframe_dF_TP,NumCell_dF_TP);
dFData_TPShuf2 = zeros(Timeframe_dF_TP,NumCell_dF_TP);

    NumCell_TPShuf2 = size(dFData_TPShuf2,2);%=size(dFData_TP,2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
    Timeframe_dF_TPShuf2 = size(dFData_TPShuf2,1);%=size(dFData_TP,1);%Time(frame)数%Time(frame)数

% 4-2-. Matrix for summarizing repeated data
List_df_TPShuf2_mean = zeros(iter_TPshuf2,NumCell_TPShuf2);
List_TPShuf2_binarySD_SyncCells = zeros(iter_TPshuf2,Timeframe_dF_TPShuf2);%for..endの開始前に配置
List_TPShuf2_SyncCells_Ratio = zeros(iter_TPshuf2,Timeframe_dF_TPShuf2);%for..endの開始前に配置
List_TPShuf2_binarySD_Max = zeros(iter_TPshuf2,NumCell_TPShuf2);
List_TPShuf2_Cell_Totalactive = zeros(iter_TPshuf2,NumCell_TPShuf2);

% List_TPShuf2_CorrMatAlan_sum = zeros(iter_TPshuf2,Timeframe_dF_TP-shift_cormat);


for n = 1:iter_TPshuf2
    % Create random IDs for Shuffle in each cells (repeat for total cells)
    for i = 1:NumCell_TPShuf2
        TPShuf2_IDs(:,i)= randperm(Timeframe_dF_TPShuf2,Timeframe_dF_TPShuf2);
    end
    % Create Shuffled data (replace raw data according to random IDs)
    for j = 1:NumCell_TPShuf2
        dFData_TP_cell_preShuf = dFData_TP(:,j);
        dFData_TPShuf2(:,j) = dFData_TP_cell_preShuf(TPShuf2_IDs(:,j));
    
%         dFData_TPShuf2_Time = [Data_imp(503:end,1),dFData_TPShuf2];%add time
    end

% %     % save the shuffled data (skip, 1000で2GBくらい)
% %     cd ([mouseID,dinfo,'_TPshuf2_',ftime]);%今作成したフォルダに移動
% %     dlmwrite([mouseID,dinfo,'_TPshuf2_',num2str(n),'.txt'], dfData_TPShuf2, 'delimiter','\t');%
% %     cd ..\
% %     cd ([mouseID,dinfo,'_TPshuf2_time',ftime]);%今作成したフォルダに移動
% %     dlmwrite([mouseID,dinfo,'_TPshuf2_',num2str(n),'_Time.txt'], dFData_TPShuf2_Time, 'delimiter','\t');%
% %     cd ..\   

    %%%%% cal. (the same flow as raw data, 3-1/2/3) %%%%%%

    % % set parameter (for...endの前に移動);SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
    % data setup (skip);dfData_TPShuf2 = dfData_TPShuf2;

    % basic infomation (mean/SD) of this data > binarization based on SD
    % threshpold
    df_TPShuf2_mean = mean(dFData_TPShuf2);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_TPShuf2_sd= std(dFData_TPShuf2);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1

    sdData_TPShuf2 = (dFData_TPShuf2 - df_TPShuf2_mean)./df_TPShuf2_sd;%各データの細胞ごとのSD
    SDbinaryData_TPShuf2 = sdData_TPShuf2>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1

% %         % save the binarized data for correlation matrix using SD-binarized data
% %         SDbinaryData_TPShuf2_Time = [Data_imp(503:end,1),SDbinaryData_TPShuf2];%add time
% %         % dlmwrite([mouseID,dinfo,'TPShuf2_3SDbinary_Time.txt'], SDbinaryData_TPShuf2_Time, 'delimiter','\t');%
% %         writematrix([mouseID,dinfo,'TPShuf2_SDbinary_Time.txt'], SDbinaryData_TPShuf2_Time, 'delimiter','\t');%

    % results from binarized df matrix
    TPShuf2_binarySD_SyncCells= sum(SDbinaryData_TPShuf2, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    %細胞数、frame数は既に取得済み
    TPShuf2_SyncCells_Ratio = TPShuf2_binarySD_SyncCells./NumCell_TPShuf2;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    TPShuf2_binarySD_Max = max(SDbinaryData_TPShuf2);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
    TPShuf2_Cell_Totalactive = sum(SDbinaryData_TPShuf2, 1);%各細胞の活動量 (3SD以上のframe数)
    TPShuf2_Cell_Totalactive_Mean = mean(TPShuf2_Cell_Totalactive);%全細胞の平均活動量

    % Summarize 
    List_df_TPShuf2_mean(n,:) = df_TPShuf2_mean;%image%TPshuf2ではTPと変わらないようにshuffleしている
    List_TPShuf2_binarySD_SyncCells(n,:) = TPShuf2_binarySD_SyncCells;%
    List_TPShuf2_SyncCells_Ratio(n,:) = TPShuf2_SyncCells_Ratio;%
    List_TPShuf2_binarySD_Max(n,:)= TPShuf2_binarySD_Max;
    List_TPShuf2_Cell_Totalactive(n,:) = TPShuf2_Cell_Totalactive;%%TPshuf2ではTPと変わらないようにshuffleしている
    
%     %Correlation matrix from Alan's(original name = out_matrix.m)
%     %How to use: CorrMat_Alan(binned_datamatrix(time, neurons),shift_window(4))
% %     shift_cormat = 4;%parameter setupに移動
%     CorrMat_Alan = CorrMat_Alan(dFData_TPShuf2,shift_cormat);
%     writematrix(CorrMat_Alan, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
% %     figure;%全shuffleをfigureにしているとデータ量が大変なので、最後のデータだけ一例として保存してみる(for...endの後にある)。
% %     imagesc(CorrMat);
% % 
% %     dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
% %     %一応もう一度時間をとる
% %     saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');
%     CorrMat_sum = sum(CorrMat_Alan);
%     List_TPShuf2_CorrMatAlan_sum(n,:) = CorrMat_sum;
end

%Cal of CorrMat > sum, sum_sum > figure(imagesc)
% % TP_CorrMat_Alan(dFData_TP,shift_cormat);toc;%biningなし、biningしたデータを入れるなら使える
tic;TPShuf2_CorrMat_Alan_skip = CorrMat_Alan_skip(dFData_TPShuf2, 19,4);toc
tic;TPShuf2_CorrMat_Alan_bin = CorrMat_Alan_bin(dFData_TPShuf2, 5,4);toc
% % writematrix(TP_CorrMat_Alan, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
%Cal of CorrMat_Sum_sum
% % TP_CorrMat_Alan_sum = sum(TP_CorrMat_Alan);
% % TP_CorrMat_sum_sum = sum(TP_CorrMat_Alan_sum);

TPShuf2_CorrMat_Alan_skip_sum = sum(TPShuf2_CorrMat_Alan_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TPShuf2_CorrMat_Alan_skip_sum_sum = sum(TPShuf2_CorrMat_Alan_skip_sum);

TPShuf2_CorrMat_Alan_bin_sum = sum(TPShuf2_CorrMat_Alan_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TPShuf2_CorrMat_Alan_bin_sum_sum = sum(TPShuf2_CorrMat_Alan_bin_sum);


% Ave_List_TP_CorrMat_sum_sum = mean(TP_CorrMat_sum_sum);%Shuffleの名残
% SD_List_TP_CorrMat_sum_sum = std(TP_CorrMat_sum_sum);
% SEM_List_TP_CorrMat_sum_sum = SD_List_TP_CorrMat_sum_sum./sqrt(iter_TPshuf2);

% % writematrix(TP_CorrMat_sum_sum, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_SumSum.txt'], 'delimiter','\t');%

figure;
imagesc(TPShuf2_CorrMat_Alan_skip);
colorbar;caxis([0 0.01]);
title("Tom+ Shuffle")
figure;
imagesc(TPShuf2_CorrMat_Alan_bin);
colorbar;caxis([0 0.01]);
title("Tom+ Shuffle")
dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%一応もう一度時間をとる
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');

%CorrMat(Asai ver.)
tic;
TPShuf2_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_TPShuf2,19,4);
toc
% writematrix(TN_CorrMat_skip, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_skip.txt'], 'delimiter','\t');%
tic;
TPShuf2_CorrMat_bin = CorrMat_Alan_As220506_bin(dFData_TPShuf2,5,4);%
toc
% writematrix(TN_CorrMat_bin, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_bin.txt'], 'delimiter','\t');%

TPShuf2_CorrMat_skip_sum = sum(TPShuf2_CorrMat_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TPShuf2_CorrMat_skip_sum_sum = sum(TPShuf2_CorrMat_skip_sum);

TPShuf2_CorrMat_bin_sum = sum(TPShuf2_CorrMat_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TPShuf2_CorrMat_bin_sum_sum = sum(TPShuf2_CorrMat_bin_sum);

figure;
imagesc(TPShuf2_CorrMat_skip);
colorbar;caxis([0 0.01]);%https://jp.mathworks.com/help/matlab/ref/caxis.html
title("Tom+ Shuffle")
    
figure;
imagesc(TPShuf2_CorrMat_bin);
colorbar;caxis([0 0.01]);%https://jp.mathworks.com/help/matlab/ref/caxis.html
title("Tom+ Shuffle")
    



% cal. Average(mean)/SD/SEM for graph
Ave_List_df_TPShuf2_mean = mean(List_df_TPShuf2_mean);
SD_Lisr_df_TPShuf2_mean = std(List_df_TPShuf2_mean);
SEM_List_df_TPShuf2_mean = SD_Lisr_df_TPShuf2_mean./sqrt(iter_TPshuf2);

Ave_List_TPShuf2_binarySD_SyncCells = mean(List_TPShuf2_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
SD_List_TPShuf2_binarySD_SyncCells = std(List_TPShuf2_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
SEM_List_TPShuf2_binarySD_SyncCells = SD_List_TPShuf2_binarySD_SyncCells./sqrt(iter_TPshuf2);%全体平均は、イベント回数をそろえているので意味がない(差がない)

Ave_List_TPShuf2_SyncCells_Ratio = mean(List_TPShuf2_SyncCells_Ratio);
SD_List_TPShuf2_SyncCells_Ratio = std(List_TPShuf2_SyncCells_Ratio);
SEM_List_TPShuf2_SyncCells_Ratio = SD_List_TPShuf2_SyncCells_Ratio./sqrt(iter_TPshuf2);

Ave_List_TPShuf2_binarySD_Max = mean(List_TPShuf2_binarySD_Max);
SD_List_TPShuf2_binarySD_Max = std(List_TPShuf2_binarySD_Max);
SEM_List_TPShuf2_binarySD_Max = SD_List_TPShuf2_binarySD_Max./sqrt(iter_TPshuf2);

Ave_List_TPShuf2_Cell_Totalactive = mean(List_TPShuf2_Cell_Totalactive);
SD_List_TPShuf2_Cell_Totalactive = std(List_TPShuf2_Cell_Totalactive);
SEM_List_TPShuf2_Cell_Totalactive = SD_List_TPShuf2_Cell_Totalactive./sqrt(iter_TPshuf2);

% List_TPShuf2_CorrMat_sum_sum = sum(List_TPShuf2_CorrMatAlan_sum);
% Ave_List_TPShuf2_CorrMat_sum_sum = mean(List_TPShuf2_CorrMat_sum_sum);
% SD_List_TPShuf2_CorrMat_sum_sum = std(List_TPShuf2_CorrMat_sum_sum);
% SEM_List_TPShuf2_CorrMat_sum_sum = SD_List_TPShuf2_CorrMat_sum_sum./sqrt(iter_TPshuf2);

% save the summary and average/SD/SEM
writematrix(List_df_TPShuf2_mean,[mouseID,dinfo,'_TPshuf2_',num2str(n),'_dF_mean.txt'], 'delimiter','\t');%
writematrix(List_TPShuf2_binarySD_SyncCells,[mouseID,dinfo,'_TPshuf2_',num2str(n),'_binarySD_SyncCells.txt'], 'delimiter','\t');%
writematrix(List_TPShuf2_SyncCells_Ratio, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_SyncCell_Ratio.txt'], 'delimiter','\t');%
writematrix(List_TPShuf2_binarySD_Max, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_binarySD_Max.txt'], 'delimiter','\t');%
writematrix(List_TPShuf2_Cell_Totalactive, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_Cell_Totalactive.txt'], 'delimiter','\t');%

% writematrix(List_TPShuf2_CorrMat_sum_sum, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_CorrMat_SumSum.txt'], 'delimiter','\t');%

% % % for comparison between original TP vs shuffled TP
% % Summary_df_TPShuf2_TP_mean = [df_TP_mean;List_df_TPShuf2_mean];%

% Create figures
% parameters for figure
xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);

%figure (Average)
figure;
% xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
plot(xc,df_TP_mean)
hold on
plot(xc,Ave_List_df_TPShuf2_mean);
legend
title('dF average');

figure;
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
plot(xt,TP_binarySD_SyncCells)
hold on
plot(xt,Ave_List_TPShuf2_binarySD_SyncCells);
legend
title('Synchronous active cells');

figure;
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
plot(xt,TP_SyncCells_Ratio)
hold on
plot(xt,Ave_List_TPShuf2_SyncCells_Ratio);
legend
title('Synchronous active cells(%)');

figure;
% xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
plot(xc,TP_binarySD_Max)
hold on
plot(xc,Ave_List_TPShuf2_binarySD_Max);
legend
title('Active(1)/inactive(0)');

figure;
% xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
plot(xc,TP_Cell_Totalactive);
hold on
plot(xc,Ave_List_TPShuf2_Cell_Totalactive);
legend
title('The number of Activity');

%corrmat
figure;
imagesc(CorrMat_Alan);
dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%一応もう一度時間をとる
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');
    

% % %figure (individual)
% % figure;
% % % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% % plot(xc,df_TP_mean)
% % hold on
% % plot(xc,List_df_TPShuf2_mean);
% % legend
% % title('Average dF in each cells');
% % 
% % figure;
% % % xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
% % plot(xt,TP_binarySD_SyncCells)
% % hold on
% % plot(xt,List_TPShuf2_binarySD_SyncCells);
% % legend
% % title('Synchronous active cells');
% % 
% % figure;
% % % xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
% % plot(xt,TP_SyncCells_Ratio)
% % hold on
% % plot(xt,List_TPShuf2_SyncCells_Ratio);
% % legend
% % title('Synchronous active cells(%)');
% % 
% % figure;
% % % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% % plot(xc,TP_binarySD_Max)
% % hold on
% % plot(xc,List_TPShuf2_binarySD_Max);
% % legend
% % title('Active(1)/inactive(0)');
% % 
% % figure;
% % % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% % plot(xc,TP_Cell_Totalactive);
% % hold on
% % plot(xc,List_TPShuf2_Cell_Totalactive);
% % legend
% % title('The number of Activity');

%statistics
%スミルノフ検定 (Statistics and Machine Learning Toolbox が必要)
%使い方を考える (Khaledは、全マウスのshuffle dataとreal dataを比較)
% [h_1,p_1,ks2stat_1] = kstest2(df_TP_mean,Ave_List_df_TPShuf2_mean,'Alpha',0.01)%[h,p,ks2stat] = kstest2(___);
% [h_2,p_2,ks2stat_2] = kstest2(TP_binarySD_SyncCells,Ave_List_TPShuf2_binarySD_SyncCells,'Alpha',0.01)%[h,p,ks2stat] = kstest2(___);



% %% 4-3. Create Shuffled data (Shuffle type 3, raster(binary) across TPShift6time and cells)
% %SDでbinarizeしてからshuffleしているので、totla eventの回数は変わらない
% % Preparation/setup > Create random IDs for Shuffle > random data > reshape & save > cal > summarize >
% % save
% 
% % % % Folders for saving the shuffled data
% % % ftime = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
% % % mkdir([mouseID,dinfo,'_TPshuf3_',ftime]) % フォルダの作成
% % % % % cd ([mouseID,dinfo,'_TPshuf3_',ftime]);%今作成したフォルダに移動
% % % mkdir([mouseID,dinfo,'_TPshuf3_time',ftime]) % フォルダの作成
% 
% %set parameters
% %Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)
% iter_TPshuf3 = 5;%1000:p=0.01を検出するためには最低でもこれくらい必要なのでは？ Num_Time_frames_TP/10;Num_Cell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)
% % NumCell_dF_TP =size(dFData_TP, 2);
% % for cal (for...end内)
% % SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
% 
% % % % data setup(shuffleしないversionの名残、消してOK、同じ内容を含んでいる)
% % % Timeframe_dF_TPShuf2 = Timeframe_dF_TP;%=size(dFData_TP,1);%Time(frame)数
% % % NumCell_dF_TPShuf2= NumCell_dF_TP;%=size(dFData_TP,2)
% 
% % prepare matrix for shuffled data and IDs
% %%%%%TPShuf3_IDs = zeros(Timeframe_dF_TP,NumCell_dF_TP);
% dfData_TPShuf3 = zeros(Timeframe_dF_TP,NumCell_dF_TP);
% 
%     NumCell_TPShuf3 = size(dfData_TPShuf3,2);%=size(dFData_TP,2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
%     Timeframe_dF_TPShuf3 = size(dfData_TPShuf3,1);%=size(dFData_TP,1);%Time(frame)数%Time(frame)数
% 
% % 4-3-. Matrix for summarizing repeated data
% List_df_TPShuf3_mean = zeros(iter_TPshuf3,NumCell_TPShuf3);
% List_TPShuf3_binarySD_SyncCells = zeros(iter_TPshuf3,Timeframe_dF_TPShuf3);%for..endの開始前に配置
% List_TPShuf3_SyncCells_Ratio = zeros(iter_TPshuf3,Timeframe_dF_TPShuf3);%for..endの開始前に配置
% List_TPShuf3_binarySD_Max = zeros(iter_TPshuf3,NumCell_TPShuf3);
% List_TPShuf3_Cell_Totalactive = zeros(iter_TPshuf3,NumCell_TPShuf3);
%     
% 
% % Reshape binarData for shuffle    
% binarData_TP_preShuf_conc = reshape(SDbinaryData_TP,[],1);%Shuffle用にデータを1次元にする
% L_binarData_TP_preShuf_conc = size(binarData_TP_preShuf_conc, 1);%データの全長、Shuffleデータを生成するために
% % Create random ID sets for Shuffle
% TPShuf3_IDs = zeros(L_binarData_TP_preShuf_conc,iter_TPshuf3);
% for i = 1:iter_TPshuf3
%     TPShuf3_IDs(:,i)= randperm(L_binarData_TP_preShuf_conc,L_binarData_TP_preShuf_conc);
% end
% % Create random Shuffled data sets (replace raw data according to random IDs)
% sdData_TP_postShuf3 = zeros(L_binarData_TP_preShuf_conc,iter_TPshuf3);
% for i = 1:iter_TPshuf3
%     sdData_TP_postShuf3(:,i) = binarData_TP_preShuf_conc(TPShuf3_IDs(:,i));
% end
% 
% % reshape (> save the shuffled data) > cal. > summarize
% % Create folders for saving shuffled data
% ftime = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
% mkdir([mouseID,dinfo,'_4_3_TPshuf3_',ftime]) % フォルダの作成
% cd ([mouseID,dinfo,'_4_3_TPshuf3_',ftime]);%今作成したフォルダに移動
% % Reshape & save
% for i = 1:iter_TPshuf3
%     % reshape shuffled binarized Data to original shape (Time,cell)
%     binarData_TPShuf3 = reshape(sdData_TP_postShuf3(:,i),size(dFData_TP));
% % %     writematrix(binarData_TPShuf3, [mouseID,dinfo,'_TPshuf3_binarData_',num2str(i),'.txt'], 'delimiter','\t');%
% % % 
% % % % binarData_Time for CorrMat
% % % mkdir([mouseID,dinfo,'_TPshuf3_Time',ftime]) % フォルダの作成
% % % cd ([mouseID,dinfo,'_TPshuf3_Time',ftime]);%今作成したフォルダに移動
% % % for i = 1:iter_Shuf3
% % %     binarData_TPShuf3_Time = [Data_imp(503:end,1),binarData_TPShuf3];%add time
% % %     writematrix(binarData_TPShuf3_Time, [mouseID,dinfo,'_TPshuf3_binarData_',num2str(i),'_Time.txt'], 'delimiter','\t');%
% % % end
% % % cd ..\
% 
%     %%%%% cal. (the same flow as raw data, 3-1/2/3) %%%%%%
%     % % set parameter (for...endの前に移動);SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
%     % data setup (skip);dfData_TPShuf2 = dfData_TPShuf2;
% 
%     % basic infomation (mean/SD) of this data > binarization based on SD
%     % threshpold
%     binar_TPShuf3_mean = mean(binarData_TPShuf3);%各細胞のmean binary, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
% %     df_TPShuf3_sd= std(binarData_TPShuf3);%各細胞のSD binary, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% 
% %     sdData_TPShuf3 = (binarData_TPShuf3 - df_TPShuf3_mean)./df_TPShuf3_sd;%各データの細胞ごとのSD
% %     SDbinaryData_TPShuf3 = sdData_TPShuf3>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% 
% 
%     % results from binarized df matrix
%     TPShuf3_binarySD_SyncCells= sum(binarData_TPShuf3, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%     %細胞数、frame数は既に取得済み
%     TPShuf3_SyncCells_Ratio = TPShuf3_binarySD_SyncCells./NumCell_TPShuf3;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
% 
%     TPShuf3_binarySD_Max = max(binarData_TPShuf3);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
%     TPShuf3_Cell_Totalactive = sum(binarData_TPShuf3, 1);%各細胞の活動量 (3SD以上のframe数)
%     TPShuf3_Cell_Totalactive_Mean = mean(TPShuf3_Cell_Totalactive);%全細胞の平均活動量
% 
%     % Summarize 
%     List_binar_TPShuf3_mean(i,:) = binar_TPShuf3_mean;%%TPshuf2ではTPと変わらないようにshuffleしている
%     List_TPShuf3_binarySD_SyncCells(i,:) = TPShuf3_binarySD_SyncCells;%
%     List_TPShuf3_SyncCells_Ratio(i,:) = TPShuf3_SyncCells_Ratio;%
%     List_TPShuf3_binarySD_Max(i,:)= TPShuf3_binarySD_Max;
%     List_TPShuf3_Cell_Totalactive(i,:) = TPShuf3_Cell_Totalactive;%%TPshuf2ではTPと変わらないようにshuffleしている
%       
% end
% cd ..\
% 
% % cal. Average(mean)/SD/SEM for graph
% Ave_List_binar_TPShuf3_mean = mean(List_binar_TPShuf3_mean);
% SD_List_binar_TPShuf3_mean = std(List_binar_TPShuf3_mean);
% SEM_List_binar_TPShuf3_mean = SD_List_binar_TPShuf3_mean./sqrt(iter_TPshuf3);
% 
% Ave_List_TPShuf3_binarySD_SyncCells = mean(List_TPShuf3_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
% SD_List_TPShuf3_binarySD_SyncCells = std(List_TPShuf3_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
% SEM_List_TPShuf3_binarySD_SyncCells = SD_List_TPShuf3_binarySD_SyncCells./sqrt(iter_TPshuf3);%全体平均は、イベント回数をそろえているので意味がない(差がない)
% 
% Ave_List_TPShuf3_SyncCells_Ratio = mean(List_TPShuf3_SyncCells_Ratio);
% SD_List_TPShuf3_SyncCells_Ratio = std(List_TPShuf3_SyncCells_Ratio);
% SEM_List_TPShuf3_SyncCells_Ratio = SD_List_TPShuf3_SyncCells_Ratio./sqrt(iter_TPshuf3);
% 
% Ave_List_TPShuf3_binarySD_Max = mean(List_TPShuf3_binarySD_Max);
% SD_List_TPShuf3_binarySD_Max = std(List_TPShuf3_binarySD_Max);
% SEM_List_TPShuf3_binarySD_Max = SD_List_TPShuf3_binarySD_Max./sqrt(iter_TPshuf3);
% 
% Ave_List_TPShuf3_Cell_Totalactive = mean(List_TPShuf3_Cell_Totalactive);
% SD_List_TPShuf3_Cell_Totalactive = std(List_TPShuf3_Cell_Totalactive);
% SEM_List_TPShuf3_Cell_Totalactive = SD_List_TPShuf3_Cell_Totalactive./sqrt(iter_TPshuf3);
% 
% % save the summary and average/SD/SEM
% writematrix(List_binar_TPShuf3_mean,[mouseID,dinfo,'_TPshuf3_',num2str(i),'_binar_mean.txt'], 'delimiter','\t');%
% writematrix(List_TPShuf3_binarySD_SyncCells,[mouseID,dinfo,'_TPshuf3_',num2str(i),'_binarySD_SyncCells.txt'], 'delimiter','\t');%
% writematrix(List_TPShuf3_SyncCells_Ratio, [mouseID,dinfo,'_TPshuf3_',num2str(i),'_SyncCell_Ratio.txt'], 'delimiter','\t');%
% writematrix(List_TPShuf3_binarySD_Max, [mouseID,dinfo,'_TPshuf3_',num2str(i),'_binarySD_Max.txt'], 'delimiter','\t');%
% writematrix(List_TPShuf3_Cell_Totalactive, [mouseID,dinfo,'_TPshuf3_',num2str(i),'_Cell_Totalactive.txt'], 'delimiter','\t');%
% % % % for comparison between original TP vs shuffled TP
% % % Summary_df_TPShuf2_TP_mean = [df_TP_mean;List_df_TPShuf2_mean];%
% 
% % Create figures
% % parameters for figure
% xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
% 
% %figure (Average)
% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,df_TP_mean)%%%%
% hold on
% plot(xc,Ave_List_binar_TPShuf3_mean);
% legend
% title('dF average');
% 
% figure;
% % xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
% plot(xt,TP_binarySD_SyncCells)
% hold on
% plot(xt,Ave_List_TPShuf3_binarySD_SyncCells);
% legend
% title('Synchronous active cells');
% 
% figure;
% % xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
% plot(xt,TP_SyncCells_Ratio)
% hold on
% plot(xt,Ave_List_TPShuf3_SyncCells_Ratio);
% legend
% title('Synchronous active cells(%)');
% 
% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,TP_binarySD_Max)
% hold on
% plot(xc,Ave_List_TPShuf3_binarySD_Max);
% legend
% title('Active(1)/inactive(0)');
% 
% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,TP_Cell_Totalactive);
% hold on
% plot(xc,Ave_List_TPShuf3_Cell_Totalactive);
% legend
% title('The number of Activity');
% 
% %figure (individual)
% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,df_TP_mean)
% hold on
% plot(xc,List_binar_TPShuf3_mean);
% legend
% title('Average dF in each cells');
% 
% figure;
% % xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
% plot(xt,TP_binarySD_SyncCells)
% hold on
% plot(xt,List_TPShuf3_binarySD_SyncCells);
% legend
% title('Synchronous active cells');
% 
% figure;
% % xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
% plot(xt,TP_SyncCells_Ratio)
% hold on
% plot(xt,List_TPShuf3_SyncCells_Ratio);
% legend
% title('Synchronous active cells(%)');
% 
% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,TP_binarySD_Max)
% hold on
% plot(xc,List_TPShuf3_binarySD_Max);
% legend
% title('Active(1)/inactive(0)');
% 
% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,TP_Cell_Totalactive);
% hold on
% plot(xc,List_TPShuf3_Cell_Totalactive);
% legend
% title('The number of Activity');
% 
% %statistics
% %スミルノフ検定
% %trial
% [h_1,p_1,ks2stat_1] = kstest2(df_TP_mean,Ave_List_binar_TPShuf3_mean,'Alpha',0.01)%[h,p,ks2stat] = kstest2(___);
% [h_2,p_2,ks2stat_2] = kstest2(TP_binarySD_SyncCells,Ave_List_TPShuf3_binarySD_SyncCells,'Alpha',0.01)%[h,p,ks2stat] = kstest2(___);

%% 4-4. Create Shuffled data (Shuffle type 4, raster(binary) across time, shuffle in each cells)
% shuffle (type2で十分、活動回数はkeepしたままshuffleできている)

%% 4-5 (4-1 circshift).　Create Shuffled data (Shuffle type 5, raster(binary) across time, shuffle in each cells)

%% 4-6 (4-2 circshift).　Create Shuffled data (Shuffle type 5, raster(binary) across time, shuffle in each cells)
%細胞ごとにshiftしているので、mean/SDは変わらない
% Preparation/setup >> Create random IDs for Shuffle > random data & save > cal > summarize >>
% save
%この計算で比較する前に、全細胞の活動頻度がそれほど高くないことを確認した方がいい (あまり活動頻度が高いと意味がない)

% Folders for saving the shifted data
ftime = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
mkdir([mouseID,dinfo,'_TPShift6_',ftime]) % フォルダの作成
% % cd ([mouseID,dinfo,'_shuf2_',ftime]);%今作成したフォルダに移動
mkdir([mouseID,dinfo,'_TPShift6_time',ftime]) % フォルダの作成

%set parameters
iter_TPShift6 = 5000;%5000-10000;%1000:p=0.01を検出するためには最低でもこれくらい必要なのでは？ Num_Time_frames_TP/10;Num_Cell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)
% NumCell_dF_TP =size(dFData_TP, 2);
% for cal (for...end内)
% SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
%for correlation matrix
% shift_cormat = 19;

% % % data setup

% prepare matrix for shifted data and IDs
TPShift6_IDs = zeros(iter_TPShift6,NumCell_dF_TP);%shiftの履歴として一応残しておく
dFData_TPShift6 = dFData_TP;%zeros(Timeframe_dF_TP,NumCell_dF_TP);

    NumCell_TPShift6 = size(dFData_TPShift6,2);%=size(dFData_TP,2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
    Timeframe_dF_TPShift6 = size(dFData_TPShift6,1);%=size(dFData_TP,1);%Time(frame)数%Time(frame)数

% 4-2-. Matrix for summarizing repeated data
List_df_TPShift6_mean = zeros(iter_TPShift6,NumCell_TPShift6);
List_TPShift6_binarySD_SyncCells = zeros(iter_TPShift6,Timeframe_dF_TPShift6);%for..endの開始前に配置
List_TPShift6_SyncCells_Ratio = zeros(iter_TPShift6,Timeframe_dF_TPShift6);%for..endの開始前に配置
List_TPShift6_binarySD_Max = zeros(iter_TPShift6,NumCell_TPShift6);
List_TPShift6_Cell_Totalactive = zeros(iter_TPShift6,NumCell_TPShift6);
   


% create shifted data > cal. > list
for n = 1:iter_TPShift6%shift回数
dFData_TPShift6 = dFData_TP;%re-set;

    % Create random IDs for Shuffle in each cells (repeat for total cells)
    for i = 1:NumCell_TPShift6%shift the time of all cells one by one
        random_TPshift6 = randperm(Timeframe_dF_TPShift6,1);%randperm(%TPの時間,1);
        dFData_TPShift6(:,i)= circshift(dFData_TPShift6(:,i),random_TPshift6);%create shifted data
        TPShift6_IDs(n,i) = random_TPshift6;%shiftの履歴として一応残しておく
    end
    
%     % for CorrMat; ただ5000回も自分でCorrMatをかけるのは現実的ではないので、Matlabでやった方がいいだろう(AlanのCorrMatは200msずつbinningしたのち、1sec分ずつcorrelationを算出していたはず)
%     dfData_TPShift6_Time = [Data_imp(503:end,1),dfData_TPShift6];%add time
%     % save the shuffled data (skip, 1000で2GBくらい)
%     cd ([mouseID,dinfo,'_TPShift6_',ftime]);%今作成したフォルダに移動
%     writematrix(dfData_TPShift6, [mouseID,dinfo,'_TPShift6_',num2str(n),'.txt'], 'delimiter','\t');%
%     cd ..\
%     cd ([mouseID,dinfo,'_TPShift6_time',ftime]);%今作成したフォルダに移動
%     writematrix(dfData_TPShift6_Time, [mouseID,dinfo,'_TPShift6_',num2str(n),'_Time.txt'], 'delimiter','\t');%
%     cd ..\   

    %%%%% cal. (the same flow as raw data, 3-1/2/3) %%%%%%
% % set parameter (for...endの前に移動);SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
    % data setup (skip);dfData_TPShift6

    % basic infomation (mean/SD) of this data > binarization based on SD
    % threshpold
    df_TPShift6_mean = mean(dFData_TPShift6);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_TPShift6_sd= std(dFData_TPShift6);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1

    sdData_TPShift6 = (dFData_TPShift6 - df_TPShift6_mean)./df_TPShift6_sd;%各データの細胞ごとのSD
    SDbinaryData_TPShift6 = sdData_TPShift6>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1

% %         % save the binarized data for correlation matrix using SD-binarized data
% %         SDbinaryData_TPShift6_Time = [Data_imp(503:end,1),SDbinaryData_TPShift6];%add time
% %         % dlmwrite([mouseID,dinfo,'TPShift6_3SDbinary_Time.txt'], SDbinaryData_TPShift6_Time, 'delimiter','\t');%
% %         writematrix([mouseID,dinfo,'TPShift6_SDbinary_Time.txt'], SDbinaryData_TPShift6_Time, 'delimiter','\t');%

    % results from binarized df matrix
    TPShift6_binarySD_SyncCells= sum(SDbinaryData_TPShift6, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    %細胞数、frame数は既に取得済み
    TPShift6_SyncCells_Ratio = TPShift6_binarySD_SyncCells./NumCell_TPShift6;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    TPShift6_binarySD_Max = max(SDbinaryData_TPShift6);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
    TPShift6_Cell_Totalactive = sum(SDbinaryData_TPShift6, 1);%各細胞の活動量 (3SD以上のframe数)
    TPShift6_Cell_Totalactive_Mean = mean(TPShift6_Cell_Totalactive);%全細胞の平均活動量

    % Summarize 
    List_df_TPShift6_mean(n,:) = df_TPShift6_mean;%%TPShift6ではTPと変わらないようにshuffleしている
    List_TPShift6_binarySD_SyncCells(n,:) = TPShift6_binarySD_SyncCells;%
    List_TPShift6_SyncCells_Ratio(n,:) = TPShift6_SyncCells_Ratio;%
    List_TPShift6_binarySD_Max(n,:)= TPShift6_binarySD_Max;
    List_TPShift6_Cell_Totalactive(n,:) = TPShift6_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
      
     %Correlation matrix from Alan's(original name = out_matrix.m)
    %How to use: CorrMat_Alan(binned_datamatrix(time, neurons),shift_window(4))
%     shift_cormat = 4;%parameter setupに移動

%     CorrMat_Alan = CorrMat_Alan(dfData_TPShift6,shift_cormat);
%     writematrix(CorrMat_Alan, [mouseID,dinfo,'_TPshift6_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
% %     tic;
% %     CorrMat_skip = CorrMat_Alan_As220506_skip(dfData_TPShift6,19,4);
% %     toc
% %     writematrix(CorrMat_skip, [mouseID,dinfo,'_TPshift6_',num2str(n),'_CorrMat_skip.txt'], 'delimiter','\t');%
% %     tic;
% %     CorrMat_bin = CorrMat_Alan_As220506_bin(dfData_TPShift6,5,4);%TNrand
% %     toc
% %     writematrix(CorrMat_bin, [mouseID,dinfo,'_TPshift6_',num2str(n),'_CorrMat_bin.txt'], 'delimiter','\t');%
%     figure;%全shuffleをfigureにしているとデータ量が大変なので、最後のデータだけ一例として保存してみる(for...endの後にある)。
%     imagesc(CorrMat);
% 
%     dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%     %一応もう一度時間をとる
%     saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');
% %     CorrMat_sum = sum(CorrMat_Alan);
% %     List_TPShift6_CorrMatAlan_sum(n,:) = CorrMat_sum;

end

%Cal of CorrMat > sum, sum_sum > figure(imagesc)
% % TP_CorrMat_Alan(dFData_TP,shift_cormat);toc;%biningなし、biningしたデータを入れるなら使える
tic;TPShift6_CorrMat_Alan_skip = CorrMat_Alan_skip(dFData_TPShift6, 19,4);toc
tic;TPShift6_CorrMat_Alan_bin = CorrMat_Alan_bin(dFData_TPShift6, 5,4);toc
% % writematrix(TP_CorrMat_Alan, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
%Cal of CorrMat_Sum_sum
% % TP_CorrMat_Alan_sum = sum(TP_CorrMat_Alan);
% % TP_CorrMat_sum_sum = sum(TP_CorrMat_Alan_sum);

TPShift6_CorrMat_Alan_skip_sum = sum(TPShift6_CorrMat_Alan_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TPShift6_CorrMat_Alan_skip_sum_sum = sum(TPShift6_CorrMat_Alan_skip_sum);

TPShift6_CorrMat_Alan_bin_sum = sum(TPShift6_CorrMat_Alan_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TPShift6_CorrMat_Alan_bin_sum_sum = sum(TPShift6_CorrMat_Alan_bin_sum);


% Ave_List_TP_CorrMat_sum_sum = mean(TP_CorrMat_sum_sum);%Shuffleの名残
% SD_List_TP_CorrMat_sum_sum = std(TP_CorrMat_sum_sum);
% SEM_List_TP_CorrMat_sum_sum = SD_List_TP_CorrMat_sum_sum./sqrt(iter_TPshuf2);

% % writematrix(TP_CorrMat_sum_sum, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_SumSum.txt'], 'delimiter','\t');%

figure;
imagesc(TPShift6_CorrMat_Alan_skip);
colorbar;caxis([0 0.01]);
title("Tom+ Shift")
figure;
imagesc(TPShift6_CorrMat_Alan_bin);
colorbar;caxis([0 0.01]);
title("Tom+ Shift")
dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%一応もう一度時間をとる
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');

%CorrMat(Asai ver.)
tic;
TPShift6_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_TPShift6,19,4);
toc
% writematrix(TN_CorrMat_skip, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_skip.txt'], 'delimiter','\t');%
tic;
TPShift6_CorrMat_bin = CorrMat_Alan_As220506_bin(dFData_TPShift6,5,4);%
toc
% writematrix(TN_CorrMat_bin, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_bin.txt'], 'delimiter','\t');%

TPShift6_CorrMat_skip_sum = sum(TPShift6_CorrMat_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TPShift6_CorrMat_skip_sum_sum = sum(TPShift6_CorrMat_skip_sum);

TPShift6_CorrMat_bin_sum = sum(TPShift6_CorrMat_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TPShift6_CorrMat_bin_sum_sum = sum(TPShift6_CorrMat_bin_sum);

figure;
imagesc(TPShift6_CorrMat_skip);
colorbar;caxis([0 0.01]);%https://jp.mathworks.com/help/matlab/ref/caxis.html
title("Tom+ Shift")

figure;
imagesc(TPShift6_CorrMat_bin);
colorbar;caxis([0 0.01]);%https://jp.mathworks.com/help/matlab/ref/caxis.html
title("Tom+ Shift")  



% cal. Average(mean)/SD/SEM for graph
Ave_List_df_TPShift6_mean = mean(List_df_TPShift6_mean);
SD_Lisr_df_TPShift6_mean = std(List_df_TPShift6_mean);
SEM_List_df_TPShift6_mean = SD_Lisr_df_TPShift6_mean./sqrt(iter_TPShift6);

Ave_List_TPShift6_binarySD_SyncCells = mean(List_TPShift6_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
SD_List_TPShift6_binarySD_SyncCells = std(List_TPShift6_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
SEM_List_TPShift6_binarySD_SyncCells = SD_List_TPShift6_binarySD_SyncCells./sqrt(iter_TPShift6);%全体平均は、イベント回数をそろえているので意味がない(差がない)

Ave_List_TPShift6_SyncCells_Ratio = mean(List_TPShift6_SyncCells_Ratio);
SD_List_TPShift6_SyncCells_Ratio = std(List_TPShift6_SyncCells_Ratio);
SEM_List_TPShift6_SyncCells_Ratio = SD_List_TPShift6_SyncCells_Ratio./sqrt(iter_TPShift6);

Ave_List_TPShift6_binarySD_Max = mean(List_TPShift6_binarySD_Max);
SD_List_TPShift6_binarySD_Max = std(List_TPShift6_binarySD_Max);
SEM_List_TPShift6_binarySD_Max = SD_List_TPShift6_binarySD_Max./sqrt(iter_TPShift6);

Ave_List_TPShift6_Cell_Totalactive = mean(List_TPShift6_Cell_Totalactive);
SD_List_TPShift6_Cell_Totalactive = std(List_TPShift6_Cell_Totalactive);
SEM_List_TPShift6_Cell_Totalactive = SD_List_TPShift6_Cell_Totalactive./sqrt(iter_TPShift6);

% save the summary and average/SD/SEM
% check first session
writematrix(List_df_TPShift6_mean,[mouseID,dinfo,'_TPShift6_',num2str(n),'_dF_mean.txt'], 'delimiter','\t');%
writematrix(List_TPShift6_binarySD_SyncCells,[mouseID,dinfo,'_TPShift6_',num2str(n),'_binarySD_SyncCells.txt'], 'delimiter','\t');%
writematrix(List_TPShift6_SyncCells_Ratio, [mouseID,dinfo,'_TPShift6_',num2str(n),'_SyncCell_Ratio.txt'], 'delimiter','\t');%
writematrix(List_TPShift6_binarySD_Max, [mouseID,dinfo,'_TPShift6_',num2str(n),'_binarySD_Max.txt'], 'delimiter','\t');%
writematrix(List_TPShift6_Cell_Totalactive, [mouseID,dinfo,'_TPShift6_',num2str(n),'_Cell_Totalactive.txt'], 'delimiter','\t');%
% % % for comparison between original TP vs shuffled TP
% % Summary_df_TPShift6_TP_mean = [df_TP_mean;List_df_TPShift6_mean];%




% Create figures
% parameters for figure
xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);

%figure (Average)
% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% % plot(xc,df_TP_mean)
% bar(xc,df_TP_mean);
% hold on
% % plot(xc,Ave_List_df_TPShift6_mean);
% bar(xc,Ave_List_df_TPShift6_mean);
% legend
% title('dF average');
figure;
bar_val = [df_TP_mean;Ave_List_df_TPShift6_mean];
bar(xc,bar_val);
legend('raw','shuffle');
title('dF average');

figure;
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
plot(xt,TP_binarySD_SyncCells)
hold on
plot(xt,Ave_List_TPShift6_binarySD_SyncCells);
legend
title('Synchronous active cells');

    for frame = [1, 1000, 2000, 3000, 4000, 5000, 6000] 
        figure;
        histogram(List_TPShift6_binarySD_SyncCells(:,frame));
        hold on
        %histogram(TP_binarySD_SyncCells(1,:));
        bar(TP_binarySD_SyncCells(frame,:),2500,0.01);%
        legend('shuffle','raw')
        title('Synchronous active cells');
    end

figure;
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
plot(xt,TP_SyncCells_Ratio,'LineWidth',1);
hold on
plot(xt,Ave_List_TPShift6_SyncCells_Ratio,'LineWidth',1);
legend('raw','shuffle');
title('Synchronous active cells(%)');

% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,TP_binarySD_Max)
% hold on
% plot(xc,Ave_List_TPShift6_binarySD_Max);
% legend
% title('Active(1)/inactive(0)');
figure;
bar_val = [TP_binarySD_Max;Ave_List_TPShift6_binarySD_Max];
bar(xc,bar_val);
legend('raw','shuffle');
title('Active(1)/inactive(0)');


% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,TP_Cell_Totalactive_Mean);
% hold on
% plot(xc,Ave_List_TPShift6_Cell_Totalactive);
% legend
% title('The number of Activity');
figure;
bar_val = [TP_Cell_Totalactive;Ave_List_TPShift6_Cell_Totalactive];
bar(xc,bar_val);
legend('raw','shuffle');
title('The number of Activity');





% % %figure (individual)
% % figure;
% % % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% % plot(xc,df_TP_mean)
% % hold on
% % plot(xc,List_df_TPShift6_mean);
% % legend
% % title('Average dF in each cells');
% % 
figure;
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
plot(xt,TP_binarySD_SyncCells)
hold on
plot(xt,List_TPShift6_binarySD_SyncCells(1:3,:),'LineWidth',0.75);
legend('raw','shuffle1')
title('Synchronous active cells');

% % figure;
% % % xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
% % plot(xt,TP_SyncCells_Ratio)
% % hold on
% % plot(xt,List_TPShift6_SyncCells_Ratio);
% % legend
% % title('Synchronous active cells(%)');
% % 
% % figure;
% % % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% % plot(xc,TP_binarySD_Max)
% % hold on
% % plot(xc,List_TPShift6_binarySD_Max);
% % legend
% % title('Active(1)/inactive(0)');
% % 
% % figure;
% % % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% % plot(xc,TP_Cell_Totalactive);
% % hold on
% % plot(xc,List_TPShift6_Cell_Totalactive);
% % legend
% % title('The number of Activity');

%% 4-7 (4-3 circshift).　Create Shuffled data (Shuffle type 5, raster(binary) across time, shuffle in each cells)

%% 4-8 (4-4 circshift).　Create Shuffled data (Shuffle type 5, raster(binary) across time, shuffle in each cells)
%%%%%%%%%%%
%% 5-1. Random sampling of TomNega (random sample cells, the same number of TomPosi)

% 5-1-0. save the shuffled data
ftime = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
mkdir([mouseID,dinfo,'_TNrand_',ftime]) % フォルダの作成
% % cd ([mouseID,dinfo,'_shuf2_',ftime]);%今作成したフォルダに移動
mkdir([mouseID,dinfo,'_TNrand_time',ftime]) % フォルダの作成

tic;
% 5-1-1. Setting parameters
% for creating randomized TN data
dFData_TN = Data_raw(:, 31:end);%extract Tom- cells/remove Tom+ cells
NumCell_TN = size(dFData_TN,2);%the number of cells
Timeframe_dF_TN = size(dFData_TN,1);%Time(frame)数

dFData_TNrand = zeros(Timeframe_dF_TN,NumCell_dF_TP);
    NumCell_TNrand = size(dFData_TNrand,2);%the number of cells,(TNrandはTPに数をそろえているので、TPの数で代用してもいいが、念のため)
    Timeframe_dF_TNrand = size(dFData_TNrand,1);%Time(frame)数

% the number of shuffle/random sampling
iter_TNrand = 150;%5000;%1000:p=0.01を検出するためには最低でもこれくらい必要なのでは？ Num_Time_frames_TP/10;Num_Cell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)

% paremeters for cal.
SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
% data setup (skip);Data_TNrand


% 5-1-2. Matrix for summarizing repeated data
List_df_TNrand_mean = zeros(iter_TNrand,NumCell_TNrand);
List_TNrand_binarySD_SyncCells = zeros(iter_TNrand,Timeframe_dF_TNrand);%for..endの開始前に配置
List_TNrand_SyncCells_Ratio = zeros(iter_TNrand,Timeframe_dF_TNrand);%for..endの開始前に配置
List_TNrand_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
List_TNrand_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);

List_TNrand_CorrMat_skip_sum = zeros(iter_TNrand,3131);%とりあえず手入力
List_TNrand_CorrMat_skip_sum_sum = zeros(iter_TNrand,1);

% 5-1-3. Create data sets of rondomized Tom Nega > cal. mean etc. >
% Summarize the results
for i = 1:iter_TNrand
    % Create random IDs for Shuffle
    TNrand_IDs = randperm(NumCell_TN,NumCell_TN);
    
    % Create Shuffled TN data
    TNrandIDs_DataTN = [TNrand_IDs; dFData_TN];%concatenate shuffled ID with raw data of Tom Nega
    TNrandIDs_DataTN_sort = sortrows(TNrandIDs_DataTN')';%sort, sortrawsで列ごとにsortするために、転置(')してshuffle IDを(1列目から)1行目に持ってきて、sort(rows)したのちに、再転置('), https://jp.mathworks.com/help/matlab/ref/double.sortrows.html
    Data_TNShuf = TNrandIDs_DataTN_sort(2:end,:);%remove shuffled IDs
    dFData_TNrand = Data_TNShuf(:,1:NumCell_dF_TP);
    
%     Data_TNrand_Time
 
%     % save the shuffled data
%     cd ([mouseID,dinfo,'_TNrand_',ftime]);%今作成したフォルダに移動
%     writematrix(Data_TNrand,[mouseID,dinfo,'_TNrand',num2str(NumCell_TP),'_',num2str(i),'.txt'],'delimiter','\t');%
%     cd ..\
%     cd ([mouseID,dinfo,'_TNrand_time',ftime]);%今作成したフォルダに移動
%     writematrix(Data_TNrand_Time,[mouseID,dinfo,'_TNrand',num2str(NumCell_TP),'_',num2str(i),'.txt'],'delimiter','\t');%
%     cd ..\

 %%%%% cal. (the same flow as raw data, 3-1/2/3) %%%%%%

    % basic infomation (mean/SD) of this data > binarization based on SD
    % threshold
    df_TNrand_mean = mean(dFData_TNrand);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_TNrand_sd= std(dFData_TNrand);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1

    sdData_TNrand = (dFData_TNrand - df_TNrand_mean)./df_TNrand_sd;%各データの細胞ごとのSD
    SDbinaryData_TNrand = sdData_TNrand>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1

    % results from binarized df matrix
    TNrand_binarySD_SyncCells= sum(SDbinaryData_TNrand, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    %細胞数、frame数は既に取得済み
    TNrand_SyncCells_Ratio = TNrand_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    TNrand_binarySD_Max = max(SDbinaryData_TNrand);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
    TNrand_Cell_Totalactive = sum(SDbinaryData_TNrand, 1);%各細胞の活動量 (3SD以上のframe数)
    TNrand_Cell_Totalactive_Mean = mean(TNrand_Cell_Totalactive);%全細胞の平均活動量

    % Summarize results from randomized data
    List_df_TNrand_mean(i,:) = df_TNrand_mean;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
    List_TNrand_binarySD_SyncCells(i,:) = TNrand_binarySD_SyncCells;%
    List_TNrand_SyncCells_Ratio(i,:) = TNrand_SyncCells_Ratio;%
    List_TNrand_binarySD_Max(i,:)= TNrand_binarySD_Max;
    List_TNrand_Cell_Totalactive(i,:) = TNrand_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
      
     %Correlation matrix from Alan's(original name = out_matrix.m)
%CorrMat(Asai ver.)
i
tic;
TNrand_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_TNrand,19,4);
toc
writematrix(TNrand_CorrMat_skip, [mouseID,dinfo,'_TP_',num2str(i),'_CorrMat_skip.txt'], 'delimiter','\t');%
%時間短縮のためにbinは省略
%     tic;
%     TNrand_CorrMat_bin = CorrMat_Alan_As220506_bin(dFData_TNrand,5,4);%
%     toc
%     % writematrix(TNrand_CorrMat_bin, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_bin.txt'], 'delimiter','\t');%

TNrand_CorrMat_skip_sum = sum(TNrand_CorrMat_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
TNrand_CorrMat_skip_sum_sum = sum(TNrand_CorrMat_skip_sum);
%時間短縮のためにbinは省略
% TNrand_CorrMat_bin_sum = sum(TNrand_CorrMat_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
% TNrand_CorrMat_bin_sum_sum = sum(TNrand_CorrMat_bin_sum);

List_TNrand_CorrMat_skip_sum(i,:) = TNrand_CorrMat_skip_sum;
List_TNrand_CorrMat_skip_sum_sum(i,:) = TNrand_CorrMat_skip_sum_sum;

end

%Cal of CorrMat > sum, sum_sum > figure(imagesc)
% % TP_CorrMat_Alan(dFData_TP,shift_cormat);toc;%biningなし、biningしたデータを入れるなら使える
%     tic;TNrand_CorrMat_Alan_skip = CorrMat_Alan_skip(dFData_TNrand, 19,4);toc
%     tic;TNrand_CorrMat_Alan_bin = CorrMat_Alan_bin(dFData_TNrand, 5,4);toc
%   % writematrix(TP_CorrMat_Alan, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
%   Cal of CorrMat_Sum_sum
% % TP_CorrMat_Alan_sum = sum(TP_CorrMat_Alan);
% % TP_CorrMat_sum_sum = sum(TP_CorrMat_Alan_sum);

    %TNrand_CorrMat_Alan_skip_sum = sum(TNrand_CorrMat_Alan_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
    %TNrand_CorrMat_Alan_skip_sum_sum = sum(TNrand_CorrMat_Alan_skip_sum);
    
    %TNrand_CorrMat_Alan_bin_sum = sum(TNrand_CorrMat_Alan_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
    %TNrand_CorrMat_Alan_bin_sum_sum = sum(TNrand_CorrMat_Alan_bin_sum);



% % writematrix(TP_CorrMat_sum_sum, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_SumSum.txt'], 'delimiter','\t');%

%     figure;
%     imagesc(TNrand_CorrMat_Alan_skip);
%     figure;
%     imagesc(TNrand_CorrMat_Alan_bin);
% dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%一応もう一度時間をとる
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');

%CorrMat(Asai ver.) for last sifted data
% tic;
% TNrand_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_TNrand,19,4);
% toc
% % writematrix(TN_CorrMat_skip, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_skip.txt'], 'delimiter','\t');%
% tic;
% TNrand_CorrMat_bin = CorrMat_Alan_As220506_bin(dFData_TNrand,5,4);%
% toc
% % writematrix(TN_CorrMat_bin, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_bin.txt'], 'delimiter','\t');%
% 
% TNrand_CorrMat_skip_sum = sum(TNrand_CorrMat_skip,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
% TNrand_CorrMat_skip_sum_sum = sum(TNrand_CorrMat_skip_sum);
% 
% TNrand_CorrMat_bin_sum = sum(TNrand_CorrMat_bin,'omitnan');%nanをomitしてsumを算出。omitしないと結果がNaNになる
% TNrand_CorrMat_bin_sum_sum = sum(TNrand_CorrMat_bin_sum);

figure;
imagesc(TNrand_CorrMat_skip);
colorbar;caxis([0 0.01]);%https://jp.mathworks.com/help/matlab/ref/caxis.html
title("Tom- random")
% 時間短縮のためなし
% figure;
% imagesc(TNrand_CorrMat_bin);
% colorbar;caxis([0 0.01]);%https://jp.mathworks.com/help/matlab/ref/caxis.html
% title("Tom- random")    


 % cal. Average(mean)/SD/SEM for graph
Ave_List_df_TNrand_mean = mean(List_df_TNrand_mean);
SD_Lisr_df_TNrand_mean = std(List_df_TNrand_mean);
SEM_List_df_TNrand_mean = SD_Lisr_df_TNrand_mean./sqrt(iter_TNrand);

Ave_List_TNrand_binarySD_SyncCells = mean(List_TNrand_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
SD_List_TNrand_binarySD_SyncCells = std(List_TNrand_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
SEM_List_TNrand_binarySD_SyncCells = SD_List_TNrand_binarySD_SyncCells./sqrt(iter_TNrand);%全体平均は、イベント回数をそろえているので意味がない(差がない)

Ave_List_TNrand_SyncCells_Ratio = mean(List_TNrand_SyncCells_Ratio);
SD_List_TNrand_SyncCells_Ratio = std(List_TNrand_SyncCells_Ratio);
SEM_List_TNrand_SyncCells_Ratio = SD_List_TNrand_SyncCells_Ratio./sqrt(iter_TNrand);

Ave_List_TNrand_binarySD_Max = mean(List_TNrand_binarySD_Max);
SD_List_TNrand_binarySD_Max = std(List_TNrand_binarySD_Max);
SEM_List_TNrand_binarySD_Max = SD_List_TNrand_binarySD_Max./sqrt(iter_TNrand);

Ave_List_TNrand_Cell_Totalactive = mean(List_TNrand_Cell_Totalactive);
SD_List_TNrand_Cell_Totalactive = std(List_TNrand_Cell_Totalactive);
SEM_List_TNrand_Cell_Totalactive = SD_List_TNrand_Cell_Totalactive./sqrt(iter_TNrand);

%CorrMat
Ave_List_TNrand_CorrMat_skip_sum = mean(List_TNrand_CorrMat_skip_sum);
SD_List_TNrand_CorrMat_skip_sum = std(List_TNrand_CorrMat_skip_sum);
SEM_List_TNrand_CorrMat_skip_sum = SD_List_TNrand_CorrMat_skip_sum./sqrt(iter_TNrand);

Ave_List_TNrand_CorrMat_skip_sum_sum = mean(List_TNrand_CorrMat_skip_sum_sum);
SD_List_TNrand_CorrMat_skip_sum_sum = std(List_TNrand_CorrMat_skip_sum_sum);
SEM_List_TNrand_CorrMat_skip_sum_sum = SD_List_TNrand_CorrMat_skip_sum_sum./sqrt(iter_TNrand);


% save the summary and average/SD/SEM
stime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%一応もう一度時間をとる

% check first session
writematrix(List_df_TNrand_mean,[mouseID,dinfo,'_TNrand_',num2str(i),'_dF_mean_List_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_TNrand_binarySD_SyncCells,[mouseID,dinfo,'_TNrand_',num2str(i),'_binarySD_SyncCells_List_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_TNrand_SyncCells_Ratio, [mouseID,dinfo,'_TNrand_',num2str(i),'_SyncCell_Ratio_List_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_TNrand_binarySD_Max, [mouseID,dinfo,'_TNrand_',num2str(i),'_binarySD_Max_List_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_TNrand_Cell_Totalactive, [mouseID,dinfo,'_TNrand_',num2str(i),'_Cell_Totalactive_List_',stime,'.txt'], 'delimiter','\t');%

writematrix(List_TNrand_CorrMat_skip_sum, [mouseID,dinfo,'_TNrand_',num2str(i),'_CorrMat_skip_sum_List_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_TNrand_CorrMat_skip_sum_sum, [mouseID,dinfo,'_TNrand_',num2str(i),'_CorrMat_skip_sum_sum_List_',stime,'.txt'], 'delimiter','\t');%
% % % for comparison between original TP vs randomized TN

%% 6. Comparison of synchrony (correlation matrix and Sync cells)

diff_CorrMat = TP_CorrMat_bin_sum - TNrand_CorrMat_bin_sum;
figure;plot(diff_CorrMat)


    %CorrMat_Alan_sum

%CorrMat_sum_sum


%     %CorrMat_Alan_sum_sun ()
%     x_corr_sum =  categorical({'Tom+','Tom+Shuffle', 'Tom+Shift', 'Tom-random'});
%     x_corr_sum = reordercats(x_corr_sum,{'Tom+','Tom+Shuffle','Tom+Shift','Tom-random'});%reorderしないとアルファベット順になるよう
% 
%     y_corrAlanskip_sum = [TP_CorrMat_Alan_skip_sum_sum,TPShuf2_CorrMat_Alan_skip_sum_sum, TPShift6_CorrMat_Alan_skip_sum_sum,TNrand_CorrMat_Alan_skip_sum_sum];
% %     figure;bar(y_corrAlanskip_sum)
%     figure;bar(x_corr_sum,y_corrAlanskip_sum)
%     title("Summation of Corrlation Matrix")
%     
%     y_corrAlanbin_sum = [TP_CorrMat_Alan_bin_sum_sum,TPShuf2_CorrMat_Alan_bin_sum_sum, TPShift6_CorrMat_Alan_bin_sum_sum,TNrand_CorrMat_Alan_bin_sum_sum];
% %     figure;bar(y_corrAlanbin_sum)
%     figure;bar(x_corr_sum,y_corrAlanskip_sum)
%     title("Summation of Corrlation Matrix")

%CorrMat_sum_sum
%https://jp.mathworks.com/help/matlab/ref/bar.html
x_corr_sum =  categorical({'Tom+','Tom+Shuffle', 'Tom+Shift', 'Tom-random'});
x_corr_sum = reordercats(x_corr_sum,{'Tom+','Tom+Shuffle','Tom+Shift','Tom-random'});%reorderしないとアルファベット順になるよう
y_corrskip_sum2 = [TP_CorrMat_skip_sum_sum,TPShuf2_CorrMat_skip_sum_sum, TPShift6_CorrMat_skip_sum_sum,TNrand_CorrMat_skip_sum_sum];
% figure;bar(y_corrskip_sum)
figure;bar(x_corr_sum,y_corrskip_sum2)
title("Summation of Corrlation Matrix")

y_corrbin_sum = [TP_CorrMat_bin_sum_sum,TPShuf2_CorrMat_bin_sum_sum, TPShift6_CorrMat_bin_sum_sum,TNrand_CorrMat_bin_sum_sum];
% figure;bar(y_corrbin_sum)
figure;bar(x_corr_sum,y_corrbin_sum)
title("Summation of Corrlation Matrix")


%% 7. Comparison of basic property (firing rate etc.)

% The number of active frames (= )


%% Create figures
% parameters for figure
xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);

%figure (Average)
figure;
% xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,df_TP_mean)
bar(xc,df_TP_mean);
hold on
% plot(xc,Ave_List_df_TNrand_mean);
bar(xc,Ave_List_df_TNrand_mean);
legend
title('dF average');
figure;
bar_val = [df_TP_mean;Ave_List_df_TNrand_mean];
bar(xc,bar_val);
legend('raw','shuffle');
title('dF average');

figure;
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
plot(xt,TP_binarySD_SyncCells)
hold on
plot(xt,Ave_List_TNrand_binarySD_SyncCells);
legend
title('Synchronous active cells');

    for frame = [1, 1000, 2000, 3000, 4000, 5000, 6000] 
        figure;
        histogram(List_TNrand_binarySD_SyncCells(:,frame));
        hold on
        %histogram(TP_binarySD_SyncCells(1,:));
        bar(TP_binarySD_SyncCells(frame,:),2500,0.01);%
        legend('shuffle','raw')
        title('Synchronous active cells');
    end

figure;
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
plot(xt,TP_SyncCells_Ratio,'LineWidth',1);
hold on
plot(xt,Ave_List_TNrand_SyncCells_Ratio,'LineWidth',1);
legend('raw','shuffle');
title('Synchronous active cells(%)');

% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,TP_binarySD_Max)
% hold on
% plot(xc,Ave_List_TNrand_binarySD_Max);
% legend
% title('Active(1)/inactive(0)');
figure;
bar_val = [TP_binarySD_Max;Ave_List_TNrand_binarySD_Max];
bar(xc,bar_val);
legend('raw','shuffle');
title('Active(1)/inactive(0)');


% figure;
% % xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
% plot(xc,TP_Cell_Totalactive_Mean);
% hold on
% plot(xc,Ave_List_TNrand_Cell_Totalactive);
% legend
% title('The number of Activity');
figure;
bar_val = [TP_Cell_Totalactive;Ave_List_TNrand_Cell_Totalactive];
bar(xc,bar_val);
legend('raw','shuffle');
title('The number of Activity');
 
toc

%%%%%%%%%%
%% Memo
% % フォルダ名の作成
% %https://jp.mathworks.com/matlabcentral/answers/479642-
% fname = char(datetime('now','Format','yyyyMMdd'));%https://jp.mathworks.com/help/matlab/ref/char.html
% % fname =
% %     '20190910'
% mkdir(fname) % フォルダの作成
% 
% %
% t = datetime('now', 'Format', 'yyyyMMddHHmmSS');%M:month, m:minute, https://jp.mathworks.com/help/matlab/ref/datetime.html?searchHighlight=datetime&s_tid=srchtitle_datetime_1#d123e285675
% DateString = datestr(t);

% フォルダ移動
%https://jp.mathworks.com/help/matlab/ref/cd.html
% cd 'C:\Users\Inokuchi_Lab\Documents\MATLAB\Asai_';%folder作成場所の指定, https://jp.mathworks.com/help/matlab/ref/cd.html
% 
% mouseID = char('b3TCW__138__');%フォルダ名, データ保存名, https://jp.mathworks.com/help/matlab/ref/char.html
% dinfo = char('R220405_8_full_D1context_');%フォルダ名
% anltype = char('cellular_SD_');
% ftime = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
% 
% mkdir([mouseID,dinfo,anltype,ftime]) % フォルダの作成
% cd ([mouseID,dinfo,anltype,ftime]);%今作成したフォルダに移動
% pwd;%https://jp.mathworks.com/help/matlab/ref/pwd.html?searchHighlight=pwd&s_tid=srchtitle_pwd_1

%shuffuled data/random data
%https://jp.mathworks.com/help/matlab/ref/randperm.html
%r1 = randperm(8,8);これを使ってrandomにする
%https://www.dogrow.net/octave/blog30/
%全データを１列にしてから、randpermを使って、shuffleして、もとの細胞数に戻す
% a = rand(1,10);
% b = randperm(10);
% c = a(b);

% % %220427-28 Circshift/randperm/fftshift/
% % %https://jp.mathworks.com/help/matlab/ref/circshift.html
% % %circshiftをデータ全体に使うと、行列がすべて同じようにshiftするので、細胞のIDは変わっても全体で平均などを取ると変わらないだろう
% % %circshiftをデータ全体に使うと、行列がすべて同じようにshiftするので、同期性は変わらないだろう
% % % > 細胞ごとに別々の値でshufleする必要がある
% % %circshift(data,x,dim) > dataをxだけ下(右、dim = 2)にシフト
% % 
% % %失敗例: 同じ値でshiftすると、位置がスライドされるだけでshuffleされていない
% % A_tri = (1:10)'
% % A_cir1 = circshift(A_tri,3)
% % A_cir2 = circshift(A_cir1,3)
% % %失敗例: 同じ値でshiftすると、位置がスライドされるだけでshuffleされていない
% % B_tri = (1:10)
% % B_circ = circshift(B_tri,4)
% % B_circ2 = circshift(B_circ,4)
% % %失敗例:同じ方向のshiftは繰り返してもshiftするだけで、shuffleはされない
% % B_circ2 = circshift(B_circ,2)
% % 
% % A_trial = [1:10; 11:20; 21:30; 31:40; 41:50;51:60;61:70;71:80;81:90;91:100]
% % A_circ1 = circshift(A_trial,3)
% % A_circ1 = circshift(A_trial(:,2),3)
% % %成功例:上書きしていく > 7列目だけshiftさせられる > これを全細胞に異なるshiftで行う
% % A_trial(:,7) = circshift(A_trial(:,7),2)
% % 
% % %失敗例:並び変えた3列目だけ、新しい３列のzero行列に反映される
% % A_circ2(:,3) = circshift(A_trial(:,3),4);
% % 
% % %こんな感じのfor...endでやれば、いいのかも
% % A_trial = [1:10; 11:20; 21:30; 31:40; 41:50;51:60;61:70;71:80;81:90;91:100];%仮想raw data
% % for i = 1:7%cell
% %     %i = 2;2列目だけshiftさせる場合
% %     for j = 1:10%shift回数
% %         k = randperm(10,1);
% %     
% %         A_trial(:,i) = circshift(A_trial(:,i),k)
% %     end
% % end
% % % shiftのpatternを用意するには...こっちの方がいいだろう
% % A_trial = [1:10; 11:20; 21:30; 31:40; 41:50;51:60;61:70;71:80;81:90;91:100];%仮想raw data
% % for j = 1:10%shift回数
% %     %i = 2;2列目だけshiftさせる場合
% %     for i = 1:7%cell
% %         k = randperm(10,1);
% %     
% %         A_trial(:,i) = circshift(A_trial(:,i),k)
% %         %このshift dataでのデータ計算・抽出
% %     end
% % end

% % % Correlation matrix (bin > correlation coefficiency with time window)
% %corrcoef/corr/corr2
% %https://jp.mathworks.com/help/matlab/ref/corrcoef.html
% % [R,P] = corrcoef(A);
% % [R,P] = corrcoef(A(t:t+5,:),A(t:t+5,:));%matrixの5frameずつcorrelationを計算
% 
% A_trial = randn(50,3);%matrix(50,3)作成
% [R_trial,P_trial] = corrcoef(A_trial);%各列でcorrelation算出
% [R_trial2,P_trial2] = corrcoef(A_trial(1:5,:));%各列間でcorrelation算出し、その総当たりを出す > time windowごとにtrimしたmatrixを作り、その総当たりcorrelationを算出(corrMatの１列分の算出)次のtime windowでも繰り返す、、、というのがいいかも
% 
% A_trial_celltime = A_trial';
% A_trial_timewindow = reshape(A_trial_celltime(),time_window);
% %corr2, https://jp.mathworks.com/help/images/ref/corr2.html?searchHighlight=corr2&s_tid=srchtitle_corr2_1
% corr2(A_trial(1:5,:),A_trial(i:i+4,:))
% corr2(A_trial(1:5,:),A_trial(1:5,:))
% %corr2は、２次元配置されたデータに対するcorr matrixなだけで、計算自体はただのcorrelationなので、
% %Alanのcorrelation matrixの計算の方が同期性を強調していていい気がする。

% %factorial = !(階乗), 
% combination = factorial(m_X)/factorial(2)/factorial(m_X-2);%Asai_220502, nC2

% %function
% function ave = average(x)
%     ave = sum(x(:))/numel(x);
% end
% %これをave.mで保存すると、ave(x)でaverageを算出できる (保存名と関数名が一致していなければならない)
% %https://jp.mathworks.com/help/matlab/ref/function.html

% %度数分布表の作成
% %tabulate() > 度数(value)、数(count)、割合(percent);https://jp.mathworks.com/help/stats/tabulate.html
% table_TNrand_SyncCells1 = tabulate(List_TNrand_binarySD_SyncCells(:,1));
% table_TNrand_SyncCells2 = tabulate(List_TNrand_binarySD_SyncCells(:,2));
% table_TNrand_SyncCells3 = tabulate(List_TNrand_binarySD_SyncCells(:,3));
% table_TNrand_SyncCells4 = tabulate(List_TNrand_binarySD_SyncCells(:,4));
% table_TNrand_SyncCells5 = tabulate(List_TNrand_binarySD_SyncCells(:,5));
% table_TNrand_SyncCells6 = tabulate(List_TNrand_binarySD_SyncCells(:,6));
% table_TNrand_SyncCells7 = tabulate(List_TNrand_binarySD_SyncCells(:,7));
% table_TNrand_SyncCells8 = tabulate(List_TNrand_binarySD_SyncCells(:,8));
% table_TNrand_SyncCells9 = tabulate(List_TNrand_binarySD_SyncCells(:,9));
% table_TNrand_SyncCells10 = tabulate(List_TNrand_binarySD_SyncCells(:,10));
% sum(table_TNrand_SyncCells(:,2));%全データ数
% figure;plot(table_TNrand_SyncCells(:,3));%度数分布の割合のplot
% 
% figure;bar(table_TNrand_SyncCells10(:,3))
% hold on
% bar(table_TP_SyncCells(:,3))
% xticklabels(0:1:9)

%table作成
%https://jp.mathworks.com/help/matlab/matlab_prog/create-a-table.html

%NaNを含むSum
%https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum%20nan&s_tid=srchtitle_sum%20nan_1
sum(TP_CorrMat_skip,'omitnan');
TP_CorrMat_skip_sum = sum(TP_CorrMat_skip,'omitnan');%のように使う

%計算時間の計測
%https://jp.mathworks.com/help/matlab/ref/toc.html?searchHighlight=toc&s_tid=srchtitle_toc_1
% tic;%ストップウォッチスタート
% toc%ストップウォッチストップ

%切り捨て/切り上げ/四捨五入
% floor();%正の数なら切り捨てに使える(負の数の場合は、切り上げになるだろう)
% ceil();%負の数の場合は、切り捨てになるだろう(正の数なら切り上げに使えるだろう)
% fix();%正でも負でも切り捨て
% %https://jp.mathworks.com/help/matlab/ref/floor.html

%figure関連
%colorbarの範囲設定
%caxis;https://jp.mathworks.com/help/matlab/ref/caxis.html
%colorbar;caxis([0 0.1]);%のように使うcmin = 0;cmax = 0.1;caxis([cmin
%cmax])のようにすると変更が楽かも
%透明度
%https://jp.mathworks.com/help/matlab/creating_plots/add-transparency-to-graphics-objects.html

%bargraph
%https://jp.mathworks.com/help/matlab/ref/bar.html
% x_corr_sum = [Tom+ Tom+Shuffle Tom+Shift Tom-random];
% y_corrskip_sum = [TP_CorrMat_skip_sum_sum,TPShuf2_CorrMat_skip_sum_sum, TPShift6_CorrMat_skip_sum_sum,TNrand_CorrMat_skip_sum_sum];
% figure;bar(y_corrskip_sum)
% figure;bar(x_corr_sum,y_corrskip_sum)
% title("Summation of Corrlation Matrix")

%save関連
% %行列データの保存
% writematrix();%Matlab2019aから導入, https://jp.mathworks.com/help/matlab/ref/writematrix.html?searchHighlight=writematrix&s_tid=srchtitle_writematrix_1
% %dlmwrite([mouseID,dinfo,'TP_3SDbinary_Time.txt'], SDbinaryData_TP_Time, 'delimiter','\t');%従来の保存
% %writematrix(SDbinaryData_TP_Time,[mouseID,dinfo,'TP_3SDbinary_Time.txt'], 'delimiter','\t');%
% writetable();
% https://jp.mathworks.com/help/matlab/ref/writetable.html#d123e1491145

%Matlab画面の色設定
%https://jp.mathworks.com/help/matlab/matlab_env/color-settings.html

%% Trial

% % %Summary_Trial (繰り返し計算結果の格納のtrial: 
% % %パラメータ設定 > 格納するmatrix =zerosの作成 > for..格納..end > mean/max... > )
% % iter_Trial = 10;
% % Num_Cell_Trial = 449;
% % Summary_Trial = zeros(iter_Trial,Num_Cell_Trial);
% % for i = 1:iter_Trial
% %     % Create random IDs for Shuffle
% %     Summary_Trial(i,:) = randperm(Num_Cell_Trial,Num_Cell_Trial);
% % end
% % Mean_Summary_Trial = mean(Summary_Trial);
% % SD_Summary_Trial = std(Summary_Trial);
% % SEM_Summary_Trial = SD_Summary_Trial./sqrt(iter_Trial);
% % Max_Summary_Trial = max(Summary_Trial);

% % %Number of iterration (繰り返し計算回数の検討: データポイント (細胞数) が449だと10/100はばらつき大きすぎる、最低でも1000くらいやるとデータが安定している
% % iter_Trial1 = 10;
% % iter_Trial2 = 100;
% % iter_Trial3 = 1000;%最大で+-5%(10%)理想的な平均からずれる
% % iter_Trial4 = 10000;%最大で+-1%(2%)理想的な平均からずれる。10000がいいだろう。
% % 
% % Num_Cell_Trial = 450;
% % Summary_Trial1 = zeros(iter_Trial1,Num_Cell_Trial);
% % Summary_Trial2 = zeros(iter_Trial2,Num_Cell_Trial);
% % Summary_Trial3 = zeros(iter_Trial3,Num_Cell_Trial);
% % Summary_Trial4 = zeros(iter_Trial4,Num_Cell_Trial);
% % for i = 1:iter_Trial1
% %     % Create random IDs for Shuffle
% %     Summary_Trial1(i,:) = randperm(Num_Cell_Trial,Num_Cell_Trial);
% % end
% % for i = 1:iter_Trial2
% %     % Create random IDs for Shuffle
% %     Summary_Trial2(i,:) = randperm(Num_Cell_Trial,Num_Cell_Trial);
% % end
% % for i = 1:iter_Trial3
% %     % Create random IDs for Shuffle
% %     Summary_Trial3(i,:) = randperm(Num_Cell_Trial,Num_Cell_Trial);
% % end
% % for i = 1:iter_Trial4
% %     % Create random IDs for Shuffle
% %     Summary_Trial4(i,:) = randperm(Num_Cell_Trial,Num_Cell_Trial);
% % end
% % 
% % Mean_Summary_Trial1 = mean(Summary_Trial1);
% % Mean_Summary_Trial2 = mean(Summary_Trial2);
% % Mean_Summary_Trial3 = mean(Summary_Trial3);
% % Mean_Summary_Trial4 = mean(Summary_Trial4);
% % 
% % SD_Summary_Trial1 = std(Summary_Trial1);
% % SD_Summary_Trial2 = std(Summary_Trial2);
% % SD_Summary_Trial3 = std(Summary_Trial3);
% % SD_Summary_Trial4 = std(Summary_Trial4);
% % 
% % % SEM_Summary_Trial1 = SD_Summary_Trial1./iter_Trial1;
% % % SEM_Summary_Trial2 = SD_Summary_Trial2./iter_Trial2;
% % % SEM_Summary_Trial3 = SD_Summary_Trial3./iter_Trial3;
% % % SEM_Summary_Trial4 = SD_Summary_Trial4./iter_Trial4;
% % SEM_Summary_Trial1 = SD_Summary_Trial1./sqrt(iter_Trial1);
% % SEM_Summary_Trial2 = SD_Summary_Trial2./sqrt(iter_Trial2);
% % SEM_Summary_Trial3 = SD_Summary_Trial3./sqrt(iter_Trial3);
% % SEM_Summary_Trial4 = SD_Summary_Trial4./sqrt(iter_Trial4);
% % 
% % figure;
% % x_tr = linspace(1,Num_Cell_Trial,Num_Cell_Trial);%https://jp.mathworks.com/help/matlab/ref/linspace.html?searchHighlight=linspace&s_tid=srchtitle_linspace_1#d123e770731
% % plot(x_tr,Mean_Summary_Trial1,'Color',[1 0 0])
% % hold on
% % plot(x_tr,Mean_Summary_Trial2,'Color',[0 1 0])
% % % plot(x_tr,Mean_Summary_Trial3,'Color',[0 0 1])
% % % plot(x_tr,Mean_Summary_Trial4,'Color',[0.3 0.3 0.3])
% % err = SEM_Summary_Trial2;
% % errorbar (x_tr,Mean_Summary_Trial3,err,'CapSize',0.5,'LineWidth',0.25,'Color',[0.3 0.3 0.3]);
% % hold off
% % legend

% % %220428 Number of iterration (randpermの回数比較) (繰り返し計算回数の検討: データポイント (細胞数) が449だと10/100はばらつき大きすぎる、最低でも1000くらいやるとデータが安定している
% % %1000ではやや回数がすくないかも、shuffleの結果のばらつきが大きい =
% shuffleの結果により有意でないということもあるだろう。5000回か10000回
% % it_Trial1 = 10;
% % it_Trial2 = 100;
% % it_Trial3 = 1000;%
% % it_Trial4 = 10000;%
% % it_Trial5 = 5000;%
% % 
% % Data_Trial1 = [1:10000];
% % % Data_Trial1 = [1:1.9:10000;2:2.9:10000;3:3.9:10000;4:4.9:10000;4:4.9:10000;5:5.9:10000;6,6.9:10000];
% % % Summary_Trial2 = zeros(iter_Trial2,Num_Cell_Trial);
% % % Summary_Trial3 = zeros(iter_Trial3,Num_Cell_Trial);
% % % Summary_Trial4 = zeros(iter_Trial4,Num_Cell_Trial);
% % shift_tria_ids = zeros(it_Trial4,1);
% % %shift_tria_ids = zeros(it_Trial4,size(Data_Trial1));
% % 
% % for i = 1:it_Trial4
% %     k = randperm(10000,1)
% % %     for j = 1:Num_Cell_Trial-10
% %     % Create random IDs for Shuffle
% % shift_tria_ids(i) = k;    
% %     Data_Trial1 = circshift(Data_Trial1,k);
% % %      A_trial(:,i) = circshift(A_trial(:,i),k)
% % %     end
% % end
% % figure;
% % plot(Data_Trial1)
% % ideal_mean = (1+10000)/2
% % mean_tria_shift_10 = mean(shift_tria_ids(1:10,1))
% % SD_tria_shift_10 = std(shift_tria_ids(1:10,1))
% % sem_tria_shift_10 = SD_tria_shift./sqrt(10)
% % err_tria_shift_10 = (mean_tria_shift_10 - ideal_mean)/ideal_mean*100
% % 
% % mean_tria_shift_100 = mean(shift_tria_ids(1:100,1))
% % SD_tria_shift_100 = std(shift_tria_ids(1:100,1))
% % sem_tria_shift_100 = SD_tria_shift./sqrt(100)
% % err_tria_shift_100 = (mean_tria_shift_100 - ideal_mean)/ideal_mean*100
% % 
% % mean_tria_shift_1000 = mean(shift_tria_ids(1:1000,1))
% % SD_tria_shift_1000 = std(shift_tria_ids(1:1000,1))
% % sem_tria_shift_1000 = SD_tria_shift./sqrt(1000)
% % err_tria_shift_1000 = (mean_tria_shift_1000 - ideal_mean)/ideal_mean*100
% % 
% % mean_tria_shift_10000 = mean(shift_tria_ids(1:10000,1))
% % SD_tria_shift_10000 = std(shift_tria_ids(1:10000,1))
% % sem_tria_shift_10000 = SD_tria_shift./sqrt(10000)
% % err_tria_shift_10000 = (mean_tria_shift_10000 - ideal_mean)/ideal_mean*100
% % 
% % mean_tria_shift_5000 = mean(shift_tria_ids(1:5000,1))
% % SD_tria_shift_5000 = std(shift_tria_ids(1:5000,1))
% % sem_tria_shift_5000 = SD_tria_shift./sqrt(5000)
% % err_tria_shift_5000 = (mean_tria_shift_5000 - ideal_mean)/ideal_mean*100
% % 
% % x_tria_bar = 1:5
% % tria_data_bar = [mean_tria_shift_10 mean_tria_shift_100 mean_tria_shift_1000 mean_tria_shift_5000 mean_tria_shift_10000]
% % tria_err_bar = [sem_tria_shift_10 sem_tria_shift_100 sem_tria_shift_1000 sem_tria_shift_5000 sem_tria_shift_10000]
% % 
% % bar(x_tria_bar,tria_data_bar)                
% % hold on
% % er = errorbar(x_tria_bar,tria_data_bar,tria_err_bar,tria_err_bar);    
% % er.Color = [0 0 0];                            
% % er.LineStyle = 'none';  
% % hold off
% % 
% % error_p = [err_tria_shift_10 err_tria_shift_100 err_tria_shift_1000 err_tria_shift_5000 err_tria_shift_10000]

% % %220428 randperm vs circshift (fftshift/flip) (10000,1)
% % %circshift/fftshiftでは、時間のスライドのみしかできない。flipでは時間が反転するため、生物らしさが少しおかしくなるかも
% % Data_Trial1 = [1:10000];
% % size_Data_Trial1 = size(Data_Trial1,2);
% % IDs_randp = randperm(size_Data_Trial1,size_Data_Trial1);
% % rand_shift = randperm(size_Data_Trial1,1);
% % IDs_circsh = circshift(Data_Trial1,rand_shift);
% % fft_shift = fftshift(Data_Trial1);
% % flip_data = flip(Data_Trial1);
% % circ_fft_shift = fftshift(IDs_circsh);
% % circ_flip = flip(IDs_circsh);
% % 
% % figure;
% % %plot(Data_Trial1,Data_Trial1);
% % scatter(Data_Trial1,Data_Trial1,1);
% % figure;
% % %scatter(Data_Trial1,IDs_randp);
% % scatter(Data_Trial1,IDs_randp,1);%size=1(default=36)にしないと大きすぎてよくわからない
% % % figure;
% % % plot(Data_Trial1,IDs_randp);%線が見えない
% % figure;
% % scatter(Data_Trial1,IDs_circsh,1);%size=1(default=36)にしないと大きすぎてよくわからない
% % 
% % figure;
% % scatter(Data_Trial1,fft_shift,1);%size=1(default=36)にしないと大きすぎてよくわからない
% % figure;
% % scatter(Data_Trial1,flip_data,1);%size=1(default=36)にしないと大きすぎてよくわからない
% % figure;
% % scatter(Data_Trial1,circ_fft_shift,1);%size=1(default=36)にしないと大きすぎてよくわからない
% % figure;
% % scatter(Data_Trial1,circ_fft_shift,1);%size=1(default=36)にしないと大きすぎてよくわからない
% % figure;
% % scatter(Data_Trial1,circ_flip,1);%size=1(default=36)にしないと大きすぎてよくわからない

% % %220428 randperm vs circshift (fftshift/flip) (5000,3)
%3細胞を想定したデータを作成したが、思ったようにはshuffleできなかった > circshiftで細胞ごとにランダムに時間をずらし、同期性を評価する
Data_Trial1 = [1:5000;1:5000;1:5000]';
size_Data_Trial1 = size(Data_Trial1,2);
sz_Data_Trial1 = size(Data_Trial1,1);

IDs_randp = randperm(size_Data_Trial1,size_Data_Trial1);
rand_shift = randperm(size_Data_Trial1,1);
IDs_circsh = circshift(Data_Trial1,rand_shift);
fft_shift = fftshift(Data_Trial1);
flip_data = flip(Data_Trial1);
circ_fft_shift = fftshift(IDs_circsh);
circ_flip = flip(IDs_circsh);

for i = 1:size_Data_Trial1%shift the time of all cells one by one
        random_Trail1 = randperm(sz_Data_Trial1,1);%randperm(%TPの時間,1);
        Data_Trial1(:,i)= circshift(Data_Trial1(:,i),random_Trail1);%create shifted data
        Data_Trial1_IDs(1,i) = random_Trail1;%shiftの履歴として一応残しておく
end
 
x_Trial1 = linspace(1,sz_Data_Trial1,sz_Data_Trial1);
% figure;
% plot(x_Trial1,Data_Trial1);
% legend
% figure;
% plot(Data_Trial1,x_Trial1);
% legend


% figure
% scatter(x_Trial1,Data_Trial1(:,1),3,"yellow","*","MarkerFaceAlpha",0.8);
% hold on
% scatter(x_Trial1,Data_Trial1(:,2),1,"magenta","+","MarkerFaceAlpha",0.8);
% scatter(x_Trial1,Data_Trial1(:,3),1,"cyan","o","MarkerFaceAlpha",0.8);
% hold off
% legend
figure
scatter(Data_Trial1(:,1),x_Trial1,3,"yellow","*","MarkerFaceAlpha",0.8);
hold on
scatter(Data_Trial1(:,2),x_Trial1,1,"magenta","+","MarkerFaceAlpha",0.8);
scatter(Data_Trial1(:,3),x_Trial1,1,"cyan","o","MarkerFaceAlpha",0.8);
hold off
legend


figure;
%scatter(Data_Trial1,IDs_randp);
scatter(x_Trial1,IDs_randp,1);%size=1(default=36)にしないと大きすぎてよくわからない
% figure;
% plot(Data_Trial1,IDs_randp);%線が見えない
figure;
scatter(x_Trial1,IDs_circsh,1);%size=1(default=36)にしないと大きすぎてよくわからない

figure;
scatter(x_Trial1,fft_shift,1);%size=1(default=36)にしないと大きすぎてよくわからない
figure;
scatter(x_Trial1,flip_data,1);%size=1(default=36)にしないと大きすぎてよくわからない
figure;
scatter(x_Trial1,circ_fft_shift,1);%size=1(default=36)にしないと大きすぎてよくわからない
figure;
scatter(x_Trial1,circ_fft_shift,1);%size=1(default=36)にしないと大きすぎてよくわからない
figure;
scatter(x_Trial1,circ_flip,1);%size=1(default=36)にしないと大きすぎてよくわからない


%figure type 1 (plot)
%figure;
% x_tri = 1:Num_Cell_Trial;,
% plot(x_tri,Mean_Summary_Trial,'-',x_tri,Summary_Trial,'--')
% plot(x_tri,Mean_Summary_Trial,'-',x_tri,Summary_Trial,'--','Color',[157 168 223]/255);
%errorbarの追加, https://jp.mathworks.com/help/matlab/ref/errorbar.html

%figure type 2 (linspace + color change)
x_tr = linspace(1,Num_Cell_Trial,Num_Cell_Trial);%https://jp.mathworks.com/help/matlab/ref/linspace.html?searchHighlight=linspace&s_tid=srchtitle_linspace_1#d123e770731
ax = axes;
hold on
plot(x_tr,Mean_Summary_Trial,x_tr,Summary_Trial)
hold off

ax.ColorOrder = [0.8 0.8 0.9; 0.2 0.2 0.8];

%figure (errorbar)
figure;
x_tr = linspace(1,Num_Cell_Trial,Num_Cell_Trial);%https://jp.mathworks.com/help/matlab/ref/linspace.html?searchHighlight=linspace&s_tid=srchtitle_linspace_1#d123e770731
% plot(x_tr,Mean_Summary_Trial,x_tr,Summary_Trial,'LineWidth',1,'Color',[0.1 0.1 0.1])
plot(x_tr,Mean_Summary_Trial,'LineWidth',0.75,'Color',[0.1 0.1 0.1])
hold on

err = SEM_Summary_Trial;
errorbar (x_tri,Mean_Summary_Trial,err,'CapSize',1,'LineWidth',0.5,'Color',[0.8 0.8 0.9]);
hold off
% hold on; https://jp.mathworks.com/help/matlab/ref/hold.html#d123e582639

% % %CorrMat並列処理を試したが、計算時間は短縮されなかった。またはエラーのため、使わない
%error
% % tic;
% % TP_CorrMat_skip_par = CorrMat_Alan_As220506_skip_par1(dFData_TP,19,4);
% % toc
    % エラー: ファイル: CorrMat_Alan_As220506_skip_par1.m、行: 29、列: 9
    % スライス化された変数 'CorrMat_Alan_As220506_skip_par1' をインデックス付けする場合、for ループ変数 'j' の範囲は、正の定数または変数の行ベクトル
    % でなければなりません。詳細については、MATLAB の並列 for ループ、"スライス化された変数をもつ入れ子の for ループ" を参照してください。
% % tic;
% % TP_CorrMat_skip_par = CorrMat_Alan_As220506_skip_par2(dFData_TP,19,4);
% % toc
    % 使い方によるエラー CorrMat_Alan_As220506_skip_par2
    % Index in position 1 exceeds array bounds. Index must not exceed 12544.
% % tic;
% % TP_CorrMat_bin_par1 = CorrMat_Alan_As220506_bin_par1(dFData_TP,5,4);%
% % toc
% % tic;
% % TP_CorrMat_bin_par2 = CorrMat_Alan_As220506_bin_par2(dFData_TP,5,4);%
% % toc
% % figure;
% % imagesc(TP_CorrMat_skip_par);
% % colorbar;caxis([0 0.01]);
% % figure;
% % imagesc(TP_CorrMat_bin_par);
% % colorbar;caxis([0 0.01]);

