%% PCA (principal component analysis)

%PCA (principal components analysis)
%https://jp.mathworks.com/help/stats/pca.html?searchHighlight=pca&s_tid=srchtitle_pca_1https://jp.mathworks.com/help/stats/pca.html?searchHighlight=pca&s_tid=srchtitle_pca_1
%[coeff_TP,score_TP,latent_TP,tsquared_TP,explained_TP,mu_TP] = pca(dFData_TP);%メインの計算はこれだけ。しかも計算もすぐ終わる
    %計算にかけるdataは(縦にtime, 横にneuronのmatrix)
    % coeff: 各PCにおけるneuronの寄与度 (縦にneuron, 横にPC1,2,,,,)
    % score: 各PCにおける各timeでの値(各timeにおける各PCの寄与度、軸を取り直したときの値) (縦にneuron, 横にPC1,2,,,,)
    % latent: 使っていない
    % tsquared: ホテリングのT2乗値、外れ値の検定に使えるらしいが使っていない
    % explained: 全体への各PCの寄与度
    % mu: 使っていない
    % この順番でしか出力できなそうなので、explainが欲しいときはcoeff,score,latent,tsquared,explainedが出力される


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


%PCA (principal components analysis)
%https://jp.mathworks.com/help/stats/pca.html?searchHighlight=pca&s_tid=srchtitle_pca_1https://jp.mathworks.com/help/stats/pca.html?searchHighlight=pca&s_tid=srchtitle_pca_1


[coeff_TPSh6,score_TPSh6,latent_TPSh6,tsquared_TPSh6,explained_TPSh6,mu_TPSh6] = pca(dFData_TPShift6);
[coeff_TNrand,score_TNrand,latent_TNrand,tsquared_TNrand,explained_TNrand,mu_TNrand] = pca(dFData_TNrand);
[coeff1_TNrand,score1_TNrand,latent1_TNrand,tsquared1_TN,explained1_TNrand] = pca(dFData_TNrand);
[coeff_TN,score_TN,latent_TN,tsquared_TN,explained_TNrand,mu_TN] = pca(dFData_TN);
[coeff_TNrand,score_TNrand,latent_TNrand,tsquared_TNrand,explained_TNrand,mu_TNrand] = pca(dFData_TNrand);
[coeff_TN,score_TN,latent_TN,tsquared_TN,explained_TN,mu_TN] = pca(dFData_TN);


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
Data_raw = Data_imp(503:end,2:end);%remove time and Tom+/- labels & crop 12000 frames

%% rules of variance (name)
%Dataはmatrix(time, cells)に使う
%Tom +(Posi) = TP, Tom -(Nega) = TN
%なるべくfunctionの名前(mean,sum,..)で、操作順に

%% 3-1. TomPosi (TP, Tom+)

    % set parameter
%     SD_thres = 3;%binaryの時のSDのthreshold
% data setup
Num_TomPosi = 30;%manualで入力, データごとに値を変える必要あり
dFData_TP = Data_raw(:,1:Num_TomPosi);%extract Tom+ cells/remove Tom- cells
%     % basic infomation (mean/SD) of this data > binarization based on SD
%     % threshpold
    df_TP_mean = mean(dFData_TP);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_TP_sd= std(dFData_TP);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    
    sdData_TP = (dFData_TP - df_TP_mean)./df_TP_sd;%各データの細胞ごとのSD(z-score)
    SDbinaryData_TP = sdData_TP>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1

    %scaling data (0-1)
    df_TP_min = min(dFData_TP);
    df_TP_max = max(dFData_TP);
    scData_TP = (dFData_TP - df_TP_min)./(df_TP_max - df_TP_min);

% %     % save the binarized data for correlation matrix using SD-binarized data
% %     SDbinaryData_TP_Time = [Data_imp(503:end,1),SDbinaryData_TP];%add time
% %     writematrix(SDbinaryData_TP_Time,[mouseID,dinfo,'TP_3SDbinary_Time.txt'],'delimiter','\t');%

%     % results from binarized df matrix
%     TP_binarySD_SyncCells= sum(SDbinaryData_TP, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%         NumCell_dF_TP =size(dFData_TP, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
%         Timeframe_dF_TP = size(dFData_TP,1);%Time(frame)数
%     TP_SyncCells_Ratio = TP_binarySD_SyncCells./NumCell_dF_TP;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
%     
%     TP_binarySD_Max = max(SDbinaryData_TP);%活動していない細胞の有無確認 (複数sessionまたいだ動画でTakekawa_systemを使ったため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
%     TP_Cell_Totalactive = sum(SDbinaryData_TP, 1);%各細胞の活動量 (3SD以上のframe数)
%     TP_Cell_Totalactive_Mean = mean(TP_Cell_Totalactive);%(TP)全細胞の平均活動量

% PCA(svd)
% tic;
[coeff_TP,score_TP,latent_TP,tsquared_TP,explained_TP,mu_TP] = pca(dFData_TP);
[coeff_sdTP,score_sdTP,latent_sdTP,tsquared_sdTP,explained_sdTP,mu_sdTP] = pca(sdData_TP);
[coeff_scTP,score_scTP,latent_scTP,tsquared_scTP,explained_scTP,mu_scTP] = pca(scData_TP);

% Save the results
    %stime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
writematrix(coeff_TP,[mouseID,dinfo,'TP_PCA_coeff.txt'],'delimiter','\t');%
writematrix(score_TP,[mouseID,dinfo,'TP_PCA_score.txt'],'delimiter','\t');%
writematrix(explained_TP,[mouseID,dinfo,'TP_PCA_explained.txt'],'delimiter','\t');%
writematrix(latent_TP,[mouseID,dinfo,'TP_PCA_latent.txt'],'delimiter','\t');%

writematrix(coeff_sdTP,[mouseID,dinfo,'TPsd_PCA_coeff.txt'],'delimiter','\t');%
writematrix(score_sdTP,[mouseID,dinfo,'TPsd_PCA_score.txt'],'delimiter','\t');%
writematrix(explained_sdTP,[mouseID,dinfo,'TPsd_PCA_explained.txt'],'delimiter','\t');%
writematrix(latent_sdTP,[mouseID,dinfo,'TPsd_PCA_latent.txt'],'delimiter','\t');%

writematrix(coeff_scTP,[mouseID,dinfo,'TPsc_PCA_coeff.txt'],'delimiter','\t');%
writematrix(score_scTP,[mouseID,dinfo,'TPsc_PCA_score.txt'],'delimiter','\t');%
writematrix(explained_scTP,[mouseID,dinfo,'TPsc_PCA_explained.txt'],'delimiter','\t');%
% toc

figure;
plot(coeff_TP())

%% 5-1. Random sampling of TomNega (random sample cells, the same number of TomPosi)

% 5-1-0. save the shuffled data
stime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
% ftime = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
mkdir([mouseID,dinfo,'_TNrand_PCA',stime]) % フォルダの作成
% % cd ([mouseID,dinfo,'_shuf2_',ftime]);%今作成したフォルダに移動
mkdir([mouseID,dinfo,'_TNrandsd_PCA',stime]) % フォルダの作成
mkdir([mouseID,dinfo,'_TNrandsc_PCA',stime]) % フォルダの作成

tic;
% 5-1-1. Setting parameters
% for creating randomized TN data
dFData_TN = Data_raw(:, 31:end);%extract Tom- cells/remove Tom+ cells
NumCell_TN = size(dFData_TN,2);%the number of cells
Timeframe_dF_TN = size(dFData_TN,1);%Time(frame)数

dFData_TNrand = zeros(Timeframe_dF_TN,NumCell_dF_TP);
    NumCell_TNrand = size(dFData_TNrand,2);%the number of cells,(TNrandはTPに数をそろえているので、TPの数で代用してもいいが、念のため)
    Timeframe_dF_TNrand = size(dFData_TNrand,1);%Time(frame)数

    tic
% the number of shuffle/random sampling
iter_TNrand = 10000;%5000;%1000:p=0.01を検出するためには最低でもこれくらい必要なのでは？ Num_Time_frames_TP/10;Num_Cell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)

% paremeters for cal.
SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
% data setup (skip);Data_TNrand


% 5-1-2. Matrix for summarizing repeated data
List_df_TNrand_mean = zeros(iter_TNrand,NumCell_TNrand);
List_TNrand_binarySD_SyncCells = zeros(iter_TNrand,Timeframe_dF_TNrand);%for..endの開始前に配置
List_TNrand_SyncCells_Ratio = zeros(iter_TNrand,Timeframe_dF_TNrand);%for..endの開始前に配置
List_TNrand_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
List_TNrand_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);

List_df_TNrand_pcaEx = zeros(NumCell_TNrand,iter_TNrand);
List_sd_TNrand_pcaEx = zeros(NumCell_TNrand,iter_TNrand);
List_sc_TNrand_pcaEx = zeros(NumCell_TNrand,iter_TNrand);

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

    %scaling data (0-1)
%     df_TP_min = min(dFData_TP);
%     df_TP_max = max(dFData_TP);
%     scData_TP = (dFData_TP - df_TP_min)./(df_TP_max - df_TP_min);

    dF_TNrand_min = min(dFData_TNrand);
    dF_TNrand_max = max(dFData_TNrand);
    scData_TNrand = (dFData_TNrand - dF_TNrand_min)./(dF_TNrand_max - dF_TNrand_min);

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


    [coeff_TNrand,score_TNrand,latent_TNrand,tsquared_TNrand,explained_TNrand,mu_TNrand] = pca(dFData_TNrand);
    [coeff_sdTNrand,score_sdTNrand,latent_sdTNrand,tsquared_sdTNrand,explained_sdTNrand,mu_sdTNrand] = pca(sdData_TNrand);
    [coeff_scTNrand,score_scTNrand,latent_scTNrand,tsquared_scTNrand,explained_scTNrand,mu_scTNrand] = pca(scData_TNrand);
    
    cd ([mouseID,dinfo,'_TNrand_PCA',stime]) % フォルダの作成
% % cd ([mouseID,dinfo,'_shuf2_',ftime]);%今作成したフォルダに移動
    writematrix(coeff_TNrand, [mouseID,dinfo,'_TNrand_',num2str(i),'_pcaCoeff.txt'], 'delimiter','\t');%
    writematrix(score_TNrand, [mouseID,dinfo,'_TNrand_',num2str(i),'_pcaScore.txt'], 'delimiter','\t');%
    List_df_TNrand_pcaEx(:,i) = explained_TNrand;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
    cd ..\

    cd ([mouseID,dinfo,'_TNrandsd_PCA',stime]) % フォルダの作成

    writematrix(coeff_sdTNrand, [mouseID,dinfo,'_TNrandsd_',num2str(i),'_pcaCoeff.txt'], 'delimiter','\t');%
    writematrix(score_sdTNrand, [mouseID,dinfo,'_TNrandsd_',num2str(i),'_pcaScore.txt'], 'delimiter','\t');%
    List_sd_TNrand_pcaEx(:,i) = explained_sdTNrand;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
    cd ..\
    
    cd ([mouseID,dinfo,'_TNrandsc_PCA',stime]) % フォルダの作成
    writematrix(coeff_scTNrand, [mouseID,dinfo,'_TNrandsc',num2str(i),'_pcaCoeff.txt'], 'delimiter','\t');%
    writematrix(score_scTNrand, [mouseID,dinfo,'_TNrandsc_',num2str(i),'_pcaScore.txt'], 'delimiter','\t');%
    List_sc_TNrand_pcaEx(:,i) = explained_scTNrand;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
    cd ..\

end

toc

% dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-

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


% save the summary and average/SD/SEM
stime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%一応もう一度時間をとる

% check first session
writematrix(List_df_TNrand_mean,[mouseID,dinfo,'_TNrand_',num2str(i),'_dF_mean_List_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_TNrand_binarySD_SyncCells,[mouseID,dinfo,'_TNrand_',num2str(i),'_binarySD_SyncCells_List_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_TNrand_SyncCells_Ratio, [mouseID,dinfo,'_TNrand_',num2str(i),'_SyncCell_Ratio_List_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_TNrand_binarySD_Max, [mouseID,dinfo,'_TNrand_',num2str(i),'_binarySD_Max_List_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_TNrand_Cell_Totalactive, [mouseID,dinfo,'_TNrand_',num2str(i),'_Cell_Totalactive_List_',stime,'.txt'], 'delimiter','\t');%

writematrix(List_df_TNrand_pcaEx, [mouseID,dinfo,'_TNrand_',num2str(i),'_PCA_explain_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_sd_TNrand_pcaEx, [mouseID,dinfo,'_TNrandsd_',num2str(i),'_PCA_explain_',stime,'.txt'], 'delimiter','\t');%
writematrix(List_sc_TNrand_pcaEx, [mouseID,dinfo,'_TNrandsc_',num2str(i),'_PCA_explain_',stime,'.txt'], 'delimiter','\t');%


%%

x = 1:30;
%cumsum (累積和) of explained
TP_PCA_expl = b3TCW138R2204058fullD1contextTPPCAexplained;%import data
TP_PCA_expl_csum = cumsum(TP_PCA_expl);

TNrand_PCA_expl = b3TCW138R2204058fullD1contextTNrand10000PCAexplain2022052620566;%import data
TNrand_PCA_expl_ave = mean(TNrand_PCA_expl,2);
TNrand_PCA_expl_ave_csum = cumsum(TNrand_PCA_expl_ave);
TNrand_PCA_expl_sd = std(TNrand_PCA_expl,0,2);%https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=std&s_tid=srchtitle_std_1
TNrand_PCA_expl_se = TNrand_PCA_expl_sd./sqrt(size(TNrand_PCA_expl,2));

figure
plot(x,TP_PCA_expl_csum,'-o','MarkerSize',5,'MarkerEdgeColor',[0 0 0], 'MarkerFaTPShuf2ceColor',[1 .6 .6]);
hold on
%plot(x,TNrand_PCA_expl_ave_csum, '-o','MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[.2 .2 .2]);
errorbar(x,TNrand_PCA_expl_ave_csum,TNrand_PCA_expl_se);


% contributed cells (from coeff)
SD_PCA_coeff = 2;

TP_PCA_coeff = b3TCW138R2204058fullD1contextTPPCAcoeff;
TP_PCA_coeff_mean = mean(TP_PCA_coeff);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
TP_PCA_coeff_sd= std(TP_PCA_coeff);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
sd_TP_PCA_coeff = (TP_PCA_coeff - TP_PCA_coeff_mean)./TP_PCA_coeff_sd;%各データの細胞ごとのSD
SDbinary_TP_PCA_coeff = sd_TP_PCA_coeff>SD_PCA_coeff;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
CellNum_TP_PCA_coeff = sum(SDbinary_TP_PCA_coeff);
%PCAでは独立になる様に各PCを決めているから独立になりやすい？各PC1-3(SD = 2), 0-1(SD=3)
%raw data(Takekawa system)をそのままPCAにかけてもensembleを抽出できない? or tdTomato+
%cellはensembleを形成していない?

TNrand_PCA_coeff = b3TCW138R2204058fullD1contextTNrand1pcaCoeff;
TNrand_PCA_coeff_mean = mean(TNrand_PCA_coeff);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
TNrand_PCA_coeff_sd= std(TNrand_PCA_coeff);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
sd_TNrand_PCA_coeff = (TNrand_PCA_coeff - TNrand_PCA_coeff_mean)./TNrand_PCA_coeff_sd;%各データの細胞ごとのSD
SDbinary_TNrand_PCA_coeff = sd_TNrand_PCA_coeff>SD_PCA_coeff;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
CellNum_TNrand_PCA_coeff = sum(SDbinary_TNrand_PCA_coeff);

TNrand_PCA_coeff = b3TCW138R2204058fullD1contextTNrand2pcaCoeff;
TNrand_PCA_coeff_mean = mean(TNrand_PCA_coeff);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
TNrand_PCA_coeff_sd= std(TNrand_PCA_coeff);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
sd_TNrand_PCA_coeff = (TNrand_PCA_coeff - TNrand_PCA_coeff_mean)./TNrand_PCA_coeff_sd;%各データの細胞ごとのSD
SDbinary_TNrand_PCA_coeff2 = sd_TNrand_PCA_coeff>SD_PCA_coeff;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
CellNum_TNrand_PCA_coeff2 = sum(SDbinary_TNrand_PCA_coeff2);
