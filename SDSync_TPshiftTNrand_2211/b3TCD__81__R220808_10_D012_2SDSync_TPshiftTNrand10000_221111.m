%% 0. clear
clear;

%このcodeはstep1/2のデータ名に関する情報とデータに関する情報をマニュアルで入力したら、step3以降は、"最後まで実行"でOK
%221024のprogressのデータを作ったコード
%% 1. Create folder for storing analysis data (Manual)
% 1-1.フォルダ名の作成
%https://jp.mathworks.com/matlabcentral/answers/479642-
cd 'C:\Users\Inokuchi_Lab\Documents\MATLAB\Asai_';%folder作成場所の指定, https://jp.mathworks.com/help/matlab/ref/cd.html

mouseID = char('b3TCD__81__');%フォルダ名, データ保存名, https://jp.mathworks.com/help/matlab/ref/char.html
dinfo = char('R220808_10_D012_');%フォルダ名
anltype = char('Rem_cell_2SDSync__');
ftime = char(datetime('now','Format','yyMMdd'));%フォルダ名, (データ名)

mkdir([mouseID,dinfo,anltype,ftime]) % フォルダの作成
cd ([mouseID,dinfo,anltype,ftime]);%今作成したフォルダに移動
pwd;%https://jp.mathworks.com/help/matlab/ref/pwd.html?searchHighlight=pwd&s_tid=srchtitle_pwd_1

%% 2. Data import > divide session by session > TPTN > TPTN matrix (Manual)

% 2-1.csv(Excel) fileをMatlab右の空白(Workspace)にドラッグ > Numeric
% matrixで、time(1列目)とcell ID/Tom+- (1or2行目)を含んだ状態でimport (:, manual)
Data_imp = b3TCD81007vtimeTomsort;%rename (manual, copy and paste from workspace)
Data_raw = Data_imp(3:end,2:end);%remove time and Tom+/-

% 2-2. the number of frames in each sessoin (Manual)
Day0_frame = 12000;
Day1pre_frame = 12000;
Day1context_frame = 12000;
Day1after_frame = 6000;
Day1post_frame = 12000;
Day2pre_frame = 12000;
Day2context_frame = 11999;
Day2after_frame = 6000;
Day2post_frame = 6000;

% first and last frames of each session
Day1pre_first = Day0_frame + 1;Day1pre_last = Day0_frame + Day1pre_frame;
Day1context_first = Day1pre_last + 1;Day1context_last = Day1pre_last + Day1context_frame;
Day1after_first = Day1context_last + 1;Day1after_last = Day1context_last + Day1after_frame;
Day1post_first = Day1after_last + 1;Day1post_last = Day1after_last + Day1post_frame;

Day2pre_first = Day1post_last + 1;Day2pre_last = Day1post_last + Day2pre_frame;
Day2context_first = Day2pre_last + 1;Day2context_last = Day2pre_last + Day2context_frame;
Day2after_first = Day2context_last + 1;Day2after_last = Day2context_last + Day2after_frame;
Day2post_first = Day2after_last + 1;Day2post_last = Day2after_last + Day2post_frame;

%Data_D0 = Data_raw(1:Day0_frame,:);
dFData_D1pre = Data_raw(Day1pre_first:Day1pre_last,:);
dFData_D1context = Data_raw(Day1context_first:Day1context_last,:);
dFData_D1ctx = dFData_D1context;
dFData_D1after = Data_raw(Day1after_first:Day1after_last,:);
dFData_D1post = Data_raw(Day1post_first:Day1post_last,:);

dFData_D2pre = Data_raw(Day2pre_first:Day2pre_last,:);
dFData_D2context = Data_raw(Day2context_first:Day2context_last,:);
dFData_D2after = Data_raw(Day2after_first:Day2after_last,:);
dFData_D2post = Data_raw(Day2post_first:Day2post_last,:);


% 2-3. the number of Tom+ and Tom- (Manual)
%Tom+= (1:TP_Cell,:) TP+? = (TP+1:TP+TP?, :) Tom-=(TN_Cell:endm,:)
TP_Cell = 4;TP_last=TP_Cell;
TPfP_Cell = 3; TPfP_first = 1+TP_last; TPfP_last = TPfP_Cell + TP_last;
TN_first = 1+TPfP_last;TN_Cell = size(Data_raw,2) - TPfP_last;
All_Cell = size(Data_raw,2); 

% dFData_D1pre_TP = dFData_D1pre(:,1:TP_last);
% dFData_D1pre_TN = dFData_D1pre(:,TN_first:end);
% dFData_D1pre_fPTN = dFData_D1pre(:,TPfP_first:end);

%% 3. cellular activity (SD, each sesion, All=TP+TPfP+TN), Day1

 % set parameter
  % set parameter
    SD_thres = 2;%binaryの時のSDのthreshold
    
    % data setup (step 2)
    
    mkdir('3_cellular2SD') % フォルダの作成
%     cd ('3_cellular2SD');%今作成したフォルダに移動, save直前に移動、function.codeがなくなり計算ができなくなるのを避ける

    % basic infomation (mean/SD/z-score) of this data > binarization based on SD
    %D1ctx. all
    cd 'C:\Users\Inokuchi_Lab\Documents\MATLAB\Asai_';%folder作成場所の指定, https://jp.mathworks.com/help/matlab/ref/cd.html
    SDbinaryData_D1ctx = func_SDbinaryData_221105(dFData_D1context,SD_thres);
        %     %%%%% func_SDbinaryData_221105の中身 %%%%%
        %     df_D1ctx_mean = mean(dFData_D1context);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
        %     df_D1ctx_sd= std(dFData_D1context);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
        %     sdData_D1ctx = (dFData_D1context - df_D1ctx_mean)./df_D1ctx_sd;%各データの細胞ごとのz-score
        %     SDbinaryData_D1ctx = sdData_D1ctx>SD_thres;%binarize, 2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
    
    % save the binarized data for correlation matrix using SD-binarized data
    cd ([mouseID,dinfo,anltype,ftime]);%今作成したフォルダに移動
    cd ('3_cellular2SD');%今作成したフォルダに移動

        dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-
        d_content = char('2SDbinar_');
        celltype = char('All_');

        session = char('D1context_');      
        writematrix(SDbinaryData_D1ctx,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%

    %%%%%%%%%%%%%%%%% save raster plot figure, Day1 context %%%%%%%%%%%%%%%%%%%%
        sesion_title = char('Day1 context');
        d_content = char('2SDbinar_raster_');
       
        %Day1context, all
        figure;imagesc(SDbinaryData_D1ctx');
        cmap_raster_all=[1 1 1;0 0 0];%white(0), black(1)
        colormap(cmap_raster_all);
%         title('Day1 context')
        title(sesion_title);
        xlabel('Time (min)');ylabel('Cell IDs');
        xticks(0:1200:12000);yticklabels(0:1:10);%
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%         saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'jpeg');  
        
        celltype = char('TP_');
        %Day1context, Tom+
        figure;imagesc(SDbinaryData_D1ctx(:,1:TP_last)');
        cmap_raster_TP=[1 1 1;1 0 1];%white(0), magenta(1)
        colormap(cmap_raster_TP);
        title(sesion_title);
        xlabel('Time (min)');ylabel('Cell IDs');
        xticks(0:1200:12000);xticklabels(0:1:10);
        yticks(1:TP_Cell);yticklabels(1:1:TP_Cell);%
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%         saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'jpeg');  
        %これにsynchroの情報を加える

        celltype = char('TN_');
        %Day1context, Tom-
        figure;imagesc(SDbinaryData_D1ctx(:,TN_first:end)');
        cmap_raster_TN=[1 1 1;0 1 1];%white(0), magenta(1)
        colormap(cmap_raster_TN);
        title(sesion_title);
        xlabel('Time (min)');ylabel('Cell IDs');
        xticks(0:1200:12000);xticklabels(0:1:10);
        saveas(gcf,[mouseID,session,celltype,d_content,dtime],'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%         saveas(gcf,[mouseID,session,celltype,d_content,dtime],'tif');
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'jpeg');  
        
        celltype = char('fPTN_');
        %Day1context, fP+Tom-
        figure;imagesc(SDbinaryData_D1ctx(:,TPfP_first:end)');
        cmap_raster_TN=[1 1 1;.3 .3 .3];%white(0), magenta(1)
        colormap(cmap_raster_TN);
        title(sesion_title);
        xlabel('Time (min)');ylabel('Cell IDs');
        xticks(0:1200:12000);yticklabels(0:1:10);%
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%         saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'jpeg');  
        
    % the number of active frame/cell
    %         All_Cell_Totalactive_Mean = mean(All_Cell_Totalactive);%全細胞の平均活動量
    All_D1ctx_Cell_Totalactive = sum(SDbinaryData_D1ctx, 1);

    % save cellular activity for combining these data from all mice
    dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-
    
    celltype = char('All_');%All, TP, TN, TNrand, 
    d_content = char('ActiveFrame_');
       
    session = char('D1context_');      
    writematrix(All_D1ctx_Cell_Totalactive,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
 
    %figure (the number of total active frames)
       
       d_content = char('ActiveFrame_');
       celltype = char('All_');
      
        session = char('D1context_');
    figure;bar(All_D1ctx_Cell_Totalactive)
    title('Active frame Day1 context')
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'tif');
    
% figure;histogram(All_D1ctx_Cell_Totalactive);

    cd ..;%上のフォルダに移動

%% 4. synchronous activity (SD, All), Day1
%この解析は、synchronous activityの簡易的なもの

 mkdir('4_2SDSynch') % フォルダの作成
    cd ('4_2SDSynch');%今作成したフォルダに移動

    % Synchronous activity (results from binarized df matrix)
    
    %Day1 context, All
    All_D1ctx_binarySD_SyncCells= sum(SDbinaryData_D1ctx, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
        NumCell_dF_All_D1ctx =size(dFData_D1ctx, 2);%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1ctx = size(dFData_D1context,1);%Time(frame)数
    D1ctx_SyncCells_Ratio = All_D1ctx_binarySD_SyncCells./NumCell_dF_All_D1ctx;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    %%%%%%%%%%%%%% SD synchronize (TP/TN) %%%%%%%%%%%%%%%%%%%%
       
    %Day1 context, TP/TN
    session = char('D1ctx_');
    d_content = char('2SDSync_');
    celltype = char('TP_');
    dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-
    
    TP_D1ctx_binarySD_SyncCells= sum(SDbinaryData_D1ctx(:,1:TP_last), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP_D1ctx =size(dFData_D1context(:,1:TP_last), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1ctx = size(dFData_D1context,1);%Time(frame)数
    TP_D1ctx_SyncCells_Ratio = TP_D1ctx_binarySD_SyncCells./NumCell_dF_TP_D1ctx;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
    writematrix(TP_D1ctx_binarySD_SyncCells,[mouseID,dinfo,celltype,session,d_content,dtime,'.txt'],'delimiter','\t');
    
%     sum(TP_D1ctx_binarySD_SyncCells)

            celltype = char('TN_');
            TN_D1ctx_binarySD_SyncCells= sum(SDbinaryData_D1ctx(:,TN_first:end), 2);
            NumCell_dF_TN_D1ctx =size(dFData_D1context(:,TN_first:end), 2);
                %         Timeframe_dF_D1ctx = size(dFData_D1ctx,1);
            TN_D1ctx_SyncCells_Ratio = TN_D1ctx_binarySD_SyncCells./NumCell_dF_TN_D1ctx;
            writematrix(TN_D1ctx_binarySD_SyncCells,[mouseID,dinfo,celltype,session,d_content,dtime,'.txt'],'delimiter','\t');
            

% figure (synchronous active cells and time) and save (All/TP/TN)
    dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-
    
    d_content = char('SynccellFrame_');

    %Day1 context
    session = char('D1context_');
% celltype = char('AllTPTN_');
% figure;plot(All_D1ctx_binarySD_SyncCells,"Color",[.5 .5 .5]);
% ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day1 context','FontSize',12);
celltype = char('TPTN_');
figure;
plot(TP_D1ctx_binarySD_SyncCells,"Color",'magenta');
hold on
plot(TN_D1ctx_binarySD_SyncCells,"Color",'cyan');
ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day1 context','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

% figure (histogram, synchronous active cells) (All) 
%     d_content = char('HistSynccellFrame_');
% 
%     %Day1 context
%       session = char('D1context_');
%         % celltype = char('All_');    
%         % figure;histogram(All_D1ctx_binarySD_SyncCells,"FaceColor",[.5 .5 .5]);
%         % xlabel('# Synchronous active cells');ylabel('# Frames');title('Day1 context','FontSize',12);
%         % saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%         % saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% celltype = char('TPTN_');
% figure;histogram(TP_D1ctx_binarySD_SyncCells,"FaceColor",'magenta');
% hold on
% histogram(TN_D1ctx_binarySD_SyncCells,"FaceColor",'cyan');
% xlabel('# Synchronous active cells');ylabel('# Frames');title('Day1 context','FontSize',12);
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

cd ..;%上のフォルダに移動

%% 5-2. synchronous SD, Day1 context TP_random_shift

%%% 4-6 (4-2 circshift).　Create Shuffled data (Shuffle type 5, raster(binary) across time, shuffle in each cells)
%細胞ごとにshiftしているので、mean/SDは変わらない
% Preparation/setup >> Create random IDs for Shuffle > random data & save > cal > summarize >>
% save
%この計算で比較する前に、全細胞の活動頻度がそれほど高くないことを確認した方がいい (あまり活動頻度が高いと意味がない)

% Folders for saving the shifted data
ftime1 = char(datetime('now','Format','HHmmss'));%フォルダ名, (データ名)
mkdir(['5_SDsync_TPShift_',ftime1]) % フォルダの作成
% cd (['5_SDsync_TPShift_',ftime1]);

%set parameters
iter_TPShift = 10000;%5000-10000;%1000:p=0.01を検出するためには最低でもこれくらい必要なのでは？ Num_Time_frames_TP/10;Num_Cell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)
% SD_thres = 3;
% NumCell_dF_TP =size(dFData_TP, 2);

%%%%%%%%% data setup for random shift %%%%%%%%%%%%%
dFData_D1ctx_TP = dFData_D1context(:,1:TP_last);%221024
        % dFData_D1ctx_All = dFData_D1context;%221024
        % dFData_D1ctx_TN = dFData_D1context(:,TN_first:end);%221024

%%%%%%%%% prepare matrix for shifted data and IDs %%%%%%%%%%%
NumCell_dF_TP = TP_Cell;%221024
% dFData_D1ctx_TP = dFData_D1context(:,1:TP_last);%221024
TPShift_IDs = zeros(iter_TPShift,NumCell_dF_TP);%shiftの履歴として一応残しておく
dfData_TPShift_D1ctx = dFData_D1ctx_TP;%zeros(Timeframe_dF_TP,NumCell_dF_TP);

    NumCell_TPShift = size(dfData_TPShift_D1ctx,2);%=size(dFData_TP,2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
    Timeframe_dF_TPShift_D1ctx = size(dfData_TPShift_D1ctx,1);%=size(dFData_TP,1);%Time(frame)数%Time(frame)数
    Timeframe_dF_TP_D1ctx = Timeframe_dF_TPShift_D1ctx;
    
% 4-2-. Matrix for summarizing repeated data
List_df_TPShift_D1ctx_mean = zeros(iter_TPShift,NumCell_TPShift);%shift dataの細胞の特性、各細胞のdf meanは変わらないはず
List_TPShift_D1ctx_Cell_Totalactive = zeros(iter_TPShift,NumCell_TPShift);%shift dataの細胞の特性、各細胞のdf meanは変わらないはず
List_TPShift_D1ctx_binarySD_Max = zeros(iter_TPShift,NumCell_TPShift);%shift dataの細胞の特性、活動していない細胞の有無確認
List_TPShift_D1ctx_binarySD_SyncCells = zeros(Timeframe_dF_TPShift_D1ctx,iter_TPShift);%for..endの開始前に配置
List_TPShift_D1ctx_SyncCells_Ratio = zeros(Timeframe_dF_TPShift_D1ctx,iter_TPShift);%for..endの開始前に配置、他のマウスと比較する時に使うかも
          
%%%%%%%%%%% create shifted data > cal. > list
tic;
for n = 1:iter_TPShift%shift回数
dfData_TPShift_D1ctx = dFData_D1ctx_TP;%re-set;

    % Create random IDs for Shuffle in each cells (repeat for total cells)
    for i = 1:NumCell_TPShift%shift the time of all cells one by one
        random_TPshift_D1ctx = randperm(Timeframe_dF_TPShift_D1ctx,1);%randperm(%TPの時間,1);
        dfData_TPShift_D1ctx(:,i)= circshift(dfData_TPShift_D1ctx(:,i),random_TPshift_D1ctx);%create shifted data
        TPShift_IDs(n,i) = random_TPshift_D1ctx;%shiftの履歴として一応残しておく
    end
   
    %%%%% cal. (the same flow as raw data, 3-1/2/3) %%%%%%
    % basic infomation (mean/SD) of this data > binarization based on SD
    % threshpold
    % Day1 context
%     SDbinaryData_TPShift_D1ctx = func_SDbinaryData_221105(dfData_TPShift_D1ctx);
    df_TPShift_D1ctx_mean = mean(dfData_TPShift_D1ctx);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_TPShift_D1ctx_sd= std(dfData_TPShift_D1ctx);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_TPShift_D1ctx = (dfData_TPShift_D1ctx - df_TPShift_D1ctx_mean)./df_TPShift_D1ctx_sd;%各データの細胞ごとのSD
    SDbinaryData_TPShift_D1ctx = sdData_TPShift_D1ctx>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1

    TPShift_D1ctx_binarySD_Max = max(SDbinaryData_TPShift_D1ctx);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
    TPShift_D1ctx_Cell_Totalactive = sum(SDbinaryData_TPShift_D1ctx, 1);%各細胞の活動量 (SD以上のframe数)
    TPShift_D1ctx_Cell_Totalactive_Mean = mean(TPShift_D1ctx_Cell_Totalactive);%全細胞の平均活動量


    % results from binarized df matrix
    TPShift_D1ctx_binarySD_SyncCells= sum(SDbinaryData_TPShift_D1ctx, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    TPShift_D1ctx_SyncCells_Ratio = TPShift_D1ctx_binarySD_SyncCells./NumCell_TPShift;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
    
   
    %%%%%% Summarize of synchronous active cells in each bin
    List_df_TPShift_D1ctx_mean(n,:) = df_TPShift_D1ctx_mean;%%TPShift6ではTPと変わらないようにshuffleしている
    List_TPShift_D1ctx_binarySD_Max(n,:)= TPShift_D1ctx_binarySD_Max;%10000回の各細胞の活動有無のまとめ
    List_TPShift_D1ctx_Cell_Totalactive(n,:) = TPShift_D1ctx_Cell_Totalactive;%10000回の各細胞のactive frame数 %TPShiftではTPと変わらないようにshuffleしている
        
    List_TPShift_D1ctx_binarySD_SyncCells(:,n) = TPShift_D1ctx_binarySD_SyncCells;%簡易同期のまとめ
    List_TPShift_D1ctx_SyncCells_Ratio(:,n) = TPShift_D1ctx_SyncCells_Ratio;%
    
end
toc
%%%%%% save the summary (individual data)
cd (['5_SDsync_TPShift_',ftime1]);
celltype = char('TPshift');
writematrix(List_df_TPShift_D1ctx_mean,[mouseID,dinfo,celltype,num2str(n),'_summary_dF_mean.txt'], 'delimiter','\t');%
writematrix(List_TPShift_D1ctx_binarySD_Max, [mouseID,dinfo,celltype,num2str(n),'_summary_binarySD_Max.txt'], 'delimiter','\t');%
writematrix(List_TPShift_D1ctx_Cell_Totalactive, [mouseID,dinfo,celltype,num2str(n),'_summary_Cell_Totalactive.txt'], 'delimiter','\t');%

writematrix(List_TPShift_D1ctx_binarySD_SyncCells,[mouseID,dinfo,celltype,num2str(n),'_summary_binarySD_SyncCells.txt'], 'delimiter','\t');%
writematrix(List_TPShift_D1ctx_SyncCells_Ratio, [mouseID,dinfo,celltype,num2str(n),'_summary_SyncCell_Ratio.txt'], 'delimiter','\t');%
% cd ..
%%%%%% cal mean/SD/SEM of iteration for crating Graph
    %         % cal. Average(mean)/SD/SEM for graph(10000 ramdom shiftの各回の各細胞のdF averageのmean/SD/SEM, rawと同じはず)
    %         Ave_List_df_TPShift_mean = mean(List_df_TPShift_D1ctx_mean);
    %         SD_Lisr_df_TPShift_mean = std(List_df_TPShift_D1ctx_mean);
    %         SEM_List_df_TPShift_mean = SD_Lisr_df_TPShift_mean./sqrt(iter_TPShift);
% 
            %10000回の各細胞の活動有無のmean/SD/SEM (不要)
            % Ave_List_TPShift_binarySD_Max = mean(List_TPShift_D1ctx_binarySD_Max);
            % SD_List_TPShift_binarySD_Max = std(List_TPShift_D1ctx_binarySD_Max);
            % SEM_List_TPShift_binarySD_Max = SD_List_TPShift_binarySD_Max./sqrt(iter_TPShift);
%     
            %10000回の各細胞のtotal active frameのmean/SD/SEM (ほぼ不要)
            % Ave_List_TPShift_Cell_Totalactive = mean(List_TPShift_D1ctx_Cell_Totalactive);
            % SD_List_TPShift_Cell_Totalactive = std(List_TPShift_D1ctx_Cell_Totalactive);
            % SEM_List_TPShift_Cell_Totalactive = SD_List_TPShift_Cell_Totalactive./sqrt(iter_TPShift);
% 
        %10000回の簡易同期のmean/SD/SEM
        %raw
%         Ave_List_TPShift_binarySD_SyncCells = mean(List_TPShift_D1ctx_binarySD_SyncCells,2);%各frameにおける10000回の簡易同期の平均
%         SD_List_TPShift_binarySD_SyncCells = std(List_TPShift_D1ctx_binarySD_SyncCells,2);%全体平均は、イベント回数をそろえているので意味がない(差がない)
%         SEM_List_TPShift_binarySD_SyncCells = SD_List_TPShift_binarySD_SyncCells./sqrt(iter_TPShift);%全体平均は、イベント回数をそろえているので意味がない(差がない)
% 
%         Sum_List_TPShift_binarySD_SyncCells =sum(List_TPShift_D1ctx_binarySD_SyncCells);%10000回の各回の合計event数、当然rawと同じになる
% 
        %10000回の簡易同期率(TP細胞中の活動している細胞の割合)のmean/SD/SEM
        % Ave_List_TPShift_SyncCells_Ratio = mean(List_TPShift_D1ctx_SyncCells_Ratio);
        % SD_List_TPShift_SyncCells_Ratio = std(List_TPShift_D1ctx_SyncCells_Ratio);
        % SEM_List_TPShift_SyncCells_Ratio = SD_List_TPShift_SyncCells_Ratio./sqrt(iter_TPShift);
%%%%%% cal mean/SD/SEM of iteration for crating Graph

%%%%%%%%%%% cal/summary of the number of frame every number of active cells
%TP raw
cd 'C:\Users\Inokuchi_Lab\Documents\MATLAB\Asai_';%folder作成場所の指定, https://jp.mathworks.com/help/matlab/ref/cd.html
Res_TP_D1ctx_binarySD_SyncCount = func_countSDsync_221105_2(TP_D1ctx_binarySD_SyncCells);
session = char('D1ctx_');
celltype = char('TP_');
d_content = char('SDsyncCount_');
cd ([mouseID,dinfo,anltype,ftime]);%今作成したフォルダに移動
cd (['5_SDsync_TPShift_',ftime1]);
writematrix(Res_TP_D1ctx_binarySD_SyncCount, [mouseID,dinfo,celltype,session,d_content,'pr.txt'], 'delimiter','\t');%

%TP_random shift
celltype = char('TPshift');
cd 'C:\Users\Inokuchi_Lab\Documents\MATLAB\Asai_';%folder作成場所の指定, https://jp.mathworks.com/help/matlab/ref/cd.html
Res_TPshift_D1ctx_binarySD_SyncCount = func_countSDsync_221105_2(List_TPShift_D1ctx_binarySD_SyncCells);
cd ([mouseID,dinfo,anltype,ftime]);%今作成したフォルダに移動
cd (['5_SDsync_TPShift_',ftime1]);
writematrix(Res_TPshift_D1ctx_binarySD_SyncCount, [mouseID,dinfo,celltype,num2str(n),session,d_content,'pr.txt'], 'delimiter','\t');%
%  figure (raster TP raw + TP SDsync)
% celltype = char('TP_');
%         %Day1context, Tom+
%         figure;imagesc(SDbinaryData_D1ctx(:,1:TP_last)');
%         cmap_raster_TP=[1 1 1;1 0 1];%white(0), magenta(1)
%         colormap(cmap_raster_TP);
%         hold on
%         imagesc(TP_D1ctx_binarySD_SyncFrame')
%         title(sesion_title);
%         xlabel('Time (min)');ylabel('Cell IDs');
%         xticks(0:1200:12000);xticklabels(0:1:10);
%         yticks(1:TP_Cell);yticklabels(1:1:TP_Cell);%
%         saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %         saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
%         saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'jpeg');  
%         %これにsynchroの情報を加える TP_D1ctx_binarySD_SyncFrame = TP_D1ctx_binarySD_SyncCells>=2;
cd ..
%% 5-2. synchronous SD, Day1 context TN_random_sample

%create rondom data%%%%%%%%%%%%%%%%%%%%%%%%%
% 5-1-0. save the shuffled data
ftime2 = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
celltype = char('TNrand_');
mkdir(['6_SDsync_TNrand_',ftime2]); % フォルダの作成
% cd(['6_SDsync_TNrand_',ftime2]);

% 5-1-1. Setting parameters
% for creating randomized TN data
dFData_TN = Data_raw(:, TN_first:end);%extract Tom- cells/remove Tom+ cells
NumCell_TN = size(dFData_TN,2);%the number of Tom- cells
Timeframe_dF_TN = size(dFData_TN,1);%Time(frame)数

    NumCell_dF_TP = TP_Cell;% the number of cells selected from Tom- cells
dFData_TNrand = zeros(Timeframe_dF_TN,NumCell_dF_TP);%TNrandのdF dataのformat作成
    NumCell_TNrand = size(dFData_TNrand,2);%the number of cells,(TNrandはTPに数をそろえているので、TPの数で代用してもいいが、念のため)
    Timeframe_dF_TNrand = size(dFData_TNrand,1);%Time(frame)数, 全session

% % %  dFData_D1ctx_TN = dFData_D1context(:,TN_first:end);%221024

% the number of shuffle/random sampling
iter_TNrand = 10000;%5000;%1000:p=0.01を検出するためには最低でもこれくらい必要なのでは？ Num_Time_frames_TP/10;Num_Cell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)

% paremeters for cal.
% SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく

% 5-1-2. Matrix for summarizing repeated data

%%%%% Day1 context %%%%%
List_df_TNrand_D1ctx_mean = zeros(iter_TNrand,NumCell_TNrand);%template for all sessions
List_TNrand_D1ctx_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
List_TNrand_D1ctx_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);

Timeframe_dF_TNrand = Day1context_frame;
List_TNrand_D1ctx_binarySD_SyncCells = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
List_TNrand_D1ctx_SyncCells_Ratio = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置

% 5-1-3. Create data sets of rondomized Tom Nega > cal. mean etc. >
% Summarize the results
tic;
for j = 1:iter_TNrand
    % Create random IDs for Shuffle
    TNrand_IDs = randperm(NumCell_TN,NumCell_TN);
    
    % Create Shuffled TN data
    TNrandIDs_DataTN = [TNrand_IDs; dFData_TN];%concatenate shuffled ID with raw data of Tom Nega
    TNrandIDs_DataTN_sort = sortrows(TNrandIDs_DataTN')';%sort, sortrawsで列ごとにsortするために、転置(')してshuffle IDを(1列目から)1行目に持ってきて、sort(rows)したのちに、再転置('), https://jp.mathworks.com/help/matlab/ref/double.sortrows.html
    Data_TNShuf = TNrandIDs_DataTN_sort(2:end,:);%remove shuffled IDs
    dFData_TNrand = Data_TNShuf(:,1:NumCell_dF_TP);
    
% devide session
dFData_TNrand_D1ctx = dFData_TNrand(Day1context_first:Day1context_last,:);%remove time and Tom+/-

%%%%%cal.%%%%%%%%%%%
%%%%% Day1 context %%%%%
   % basic infomation (mean/SD) of this data > binarization based on SD
   % threshold
    df_TNrand_D1ctx_mean = mean(dFData_TNrand_D1ctx);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_TNrand_D1ctx_sd= std(dFData_TNrand_D1ctx);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_TNrand_D1ctx = (dFData_TNrand_D1ctx - df_TNrand_D1ctx_mean)./df_TNrand_D1ctx_sd;%各データの細胞ごとのSD
    SDbinaryData_TNrand_D1ctx = sdData_TNrand_D1ctx>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1


    TNrand_D1ctx_binarySD_Max = max(SDbinaryData_TNrand_D1ctx);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
    TNrand_D1ctx_Cell_Totalactive = sum(SDbinaryData_TNrand_D1ctx, 1);%各細胞の活動量 (2SD以上のframe数)
    TNrand_D1ctx_Cell_Totalactive_Mean = mean(TNrand_D1ctx_Cell_Totalactive);%全細胞の平均活動量
    % results from binarized df matrix
    TNrand_D1ctx_binarySD_SyncCells= sum(SDbinaryData_TNrand_D1ctx, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    TNrand_D1ctx_SyncCells_Ratio = TNrand_D1ctx_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
  

    % Summarize results from randomized data
    List_df_TNrand_D1ctx_mean(j,:) = df_TNrand_D1ctx_mean;%10000回のランダムサンプリングの各細胞のdf mean
    List_TNrand_D1ctx_binarySD_Max(j,:)= TNrand_D1ctx_binarySD_Max;%10000回のランダムサンプリングの各細胞の活動有無のまとめ
    List_TNrand_D1ctx_Cell_Totalactive(j,:) = TNrand_D1ctx_Cell_Totalactive;%10000回のランダムサンプリングの各細胞のactive frame数
    %%%%%% Summarize of synchronous active cells in each bin
    List_TNrand_D1ctx_binarySD_SyncCells(:,j) = TNrand_D1ctx_binarySD_SyncCells;%簡易同期
    List_TNrand_D1ctx_SyncCells_Ratio(:,j) = TNrand_D1ctx_SyncCells_Ratio;%

end         
toc
%%%%%% save the summary (individual data)
cd(['6_SDsync_TNrand_',ftime2]); % フォルダの作成
celltype = char('TNrand');
writematrix(List_df_TNrand_D1ctx_mean,[mouseID,dinfo,celltype,num2str(j),'_summary_dF_mean.txt'], 'delimiter','\t');%
writematrix(List_TNrand_D1ctx_binarySD_Max, [mouseID,dinfo,celltype,num2str(j),'_summary_binarySD_Max.txt'], 'delimiter','\t');%
writematrix(List_TNrand_D1ctx_Cell_Totalactive, [mouseID,dinfo,celltype,num2str(j),'_summary_Cell_Totalactive.txt'], 'delimiter','\t');%

writematrix(List_TNrand_D1ctx_binarySD_SyncCells,[mouseID,dinfo,celltype,num2str(j),'_summary_binarySD_SyncCells.txt'], 'delimiter','\t');%
writematrix(List_TNrand_D1ctx_SyncCells_Ratio, [mouseID,dinfo,celltype,num2str(j),'_summary_SyncCell_Ratio.txt'], 'delimiter','\t');%
cd ..
%%%%%%%%%%% cal/summary of the number of frame every number of active cells
%TP raw
% Res_TP_D1ctx_binarySD_SyncCount = func_countSDsync_221105(TP_D1ctx_binarySD_SyncCells);
% celltype = char('TP_');
% d_content = char('SDsyncCount_');
% writematrix(Res_TP_D1ctx_binarySD_SyncCount, [mouseID,dinfo,celltype,d_content,'pr.txt'], 'delimiter','\t');%
%TP_random shift
celltype = char('TNrand');
d_content = char('SDsyncCount_');
cd 'C:\Users\Inokuchi_Lab\Documents\MATLAB\Asai_';%folder作成場所の指定, https://jp.mathworks.com/help/matlab/ref/cd.html
Res_TNrand_D1ctx_binarySD_SyncCount = func_countSDsync_221105_2(List_TNrand_D1ctx_binarySD_SyncCells);
cd ([mouseID,dinfo,anltype,ftime]);%今作成したフォルダに移動
cd(['6_SDsync_TNrand_',ftime2]); % フォルダの作成
writematrix(Res_TNrand_D1ctx_binarySD_SyncCount, [mouseID,dinfo,celltype,num2str(j),session,d_content,'pr.txt'], 'delimiter','\t');%
