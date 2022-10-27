%% 0. clear
clear;

%このcodeはstep1/2のデータ名に関する情報とデータに関する情報をマニュアルで入力したら、step3以降は、"最後まで実行"でOK
%221024のprogressのデータを作ったコード
%% 1. Create folder for storing analysis data (Manual)
% 1-1.フォルダ名の作成
%https://jp.mathworks.com/matlabcentral/answers/479642-
cd 'C:\Users\Inokuchi_Lab\Documents\MATLAB\Asai_';%folder作成場所の指定, https://jp.mathworks.com/help/matlab/ref/cd.html

mouseID = char('b3TCW__126__');%フォルダ名, データ保存名, https://jp.mathworks.com/help/matlab/ref/char.html
dinfo = char('R220214_15_D12_');%フォルダ名
anltype = char('Rem_cell_D1cxCorMat_SD_try_');
ftime = char(datetime('now','Format','yyMMdd'));%フォルダ名, (データ名)

mkdir([mouseID,dinfo,anltype,ftime]) % フォルダの作成
cd ([mouseID,dinfo,anltype,ftime]);%今作成したフォルダに移動
pwd;%https://jp.mathworks.com/help/matlab/ref/pwd.html?searchHighlight=pwd&s_tid=srchtitle_pwd_1

%% 2. Data import > divide session by session > TPTN > TPTN matrix (Manual)

% 2-1.csv(Excel) fileをMatlab右の空白(Workspace)にドラッグ > Numeric
% matrixで、time(1列目)とcell ID/Tom+- (1or2行目)を含んだ状態でimport (:LE, manual)
Data_imp = b3TCW126007vtiemTomOvsortmcor;%rename (manual, copy and paste from workspace)
Data_raw = Data_imp(3:end,2:end);%remove time and Tom+/-

% 2-2. the number of frames in each sessoin (Manual)
Day0_frame = 0;
Day1pre_frame = 11997;
Day1context_frame = 12000;
Day1after_frame = 0;
Day1post_frame = 11999;
Day2pre_frame = 11997;
Day2context_frame = 11997;
Day2after_frame = 0;
Day2post_frame = 11999;

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
TP_Cell = 3;TP_last=TP_Cell;
TPfP_Cell = 0; TPfP_first = 1+TP_last; TPfP_last = TPfP_Cell + TP_last;
TN_first = 1+TPfP_last;TN_Cell = size(Data_raw,2) - TPfP_last;

% dFData_D1pre_TP = dFData_D1pre(:,1:TP_last);
% dFData_D1pre_TN = dFData_D1pre(:,TN_first:end);
% dFData_D1pre_fPTN = dFData_D1pre(:,TPfP_first:end);



%% 3. cellular activity (SD, each sesion, All=TP+TPfP+TN), Day1

 % set parameter
    SD_thres = 2;%binaryの時のSDのthreshold
    
    % data setup (step 2)
    
    mkdir('3_cellular2SD') % フォルダの作成
    cd ('3_cellular2SD');%今作成したフォルダに移動

    % basic infomation (mean/SD/z-score) of this data > binarization based on SD
    %D1pre. all
    df_D1pre_mean = mean(dFData_D1pre);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_D1pre_sd= std(dFData_D1pre);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_D1pre = (dFData_D1pre - df_D1pre_mean)./df_D1pre_sd;%各データの細胞ごとのz-score
    SDbinaryData_D1pre = sdData_D1pre>SD_thres;%binarize, 2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
    %D1ctx. all
    df_D1ctx_mean = mean(dFData_D1context);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_D1ctx_sd= std(dFData_D1context);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_D1ctx = (dFData_D1context - df_D1ctx_mean)./df_D1ctx_sd;%各データの細胞ごとのz-score
    SDbinaryData_D1ctx = sdData_D1ctx>SD_thres;%binarize, 2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
    %D1after. all
    df_D1after_mean = mean(dFData_D1after);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_D1after_sd= std(dFData_D1after);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_D1after = (dFData_D1after - df_D1after_mean)./df_D1after_sd;%各データの細胞ごとのz-score
    SDbinaryData_D1after = sdData_D1after>SD_thres;%binarize, 2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
    %D1post. all    
    df_D1post_mean = mean(dFData_D1post);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_D1post_sd= std(dFData_D1post);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_D1post = (dFData_D1post - df_D1post_mean)./df_D1post_sd;%各データの細胞ごとのz-score
    SDbinaryData_D1post = sdData_D1post>SD_thres;%binarize, 2SD以上を1, それ未満を0に。s://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
    
        % save the binarized data for correlation matrix using SD-binarized data
        dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-
        d_content = char('2SDbinar_');
        celltype = char('All_');
%         SDbinaryData_Time = [Data_imp(503:end,1),SDbinaryData_D1ctx];%add time
        %dlmwrite([mouseID,dinfo,'TP_2SDbinary_Time.txt'], SDbinaryData_TP_Time, 'delimiter','\t');%
        session = char('D1pre_');      
        writematrix(SDbinaryData_D1pre,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
        session = char('D1context_');      
        writematrix(SDbinaryData_D1ctx,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
        session = char('D1after_');      
        writematrix(SDbinaryData_D1after,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
        session = char('D1post_');      
        writematrix(SDbinaryData_D1post,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%

        %save raster plot figure
        session = char('D1context_');
        sesion_title = char('Day1 context');
        d_content = char('2SDbinar_raster_');
       
        celltype = char('All_');

        %Day1context, all
        figure;imagesc(SDbinaryData_D1ctx');
        cmap_raster_all=[1 1 1;0 0 0];%white(0), black(1)
        colormap(cmap_raster_all);
%         title('Day1 context')
        title(sesion_title);
        xlabel('Time (min)');ylabel('Cell IDs');
        xticks(0:1200:12000);xticklabels(0:1:10);
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
        yticks(1:TP_Cell);xticklabels(1:1:TP_Cell);
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%         saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'jpeg');  
        
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
        xticks(0:1200:12000);xticklabels(0:1:10);
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%         saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'jpeg');  
        
        % the number of active frame
    % All_D1pre_binarySD__Max = max(SDbinaryData_D1pre);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
    All_D1pre_Cell_Totalactive = sum(SDbinaryData_D1pre, 1);%各細胞の活動量 (2SD以上のframe数)
    %         All_Cell_Totalactive_Mean = mean(All_Cell_Totalactive);%全細胞の平均活動量
    All_D1ctx_Cell_Totalactive = sum(SDbinaryData_D1ctx, 1);
    All_D1after_Cell_Totalactive = sum(SDbinaryData_D1after, 1);
    All_D1post_Cell_Totalactive = sum(SDbinaryData_D1post, 1);

    % save cellular activity for combining these data from all mice
    dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
    
    celltype = char('All_');%All, TP, TN, TNrand, 
    d_content = char('ActiveFrame_');
       
    session = char('D1pre_');
    writematrix(All_D1pre_Cell_Totalactive,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
    session = char('D1context_');      
    writematrix(All_D1ctx_Cell_Totalactive,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
    session = char('D1after_');      
    writematrix(All_D1after_Cell_Totalactive,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
    session = char('D1post_');      
    writematrix(All_D1post_Cell_Totalactive,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
    
    %figure (the number of total active frames)
       
       d_content = char('ActiveFrame_');
       celltype = char('All_');
    
       session = char('D1pre_');

    figure;bar(All_D1pre_Cell_Totalactive)
    title('Active frame Day1 pre')
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'tif');
    
        session = char('D1context_');
    figure;bar(All_D1ctx_Cell_Totalactive)
    title('Active frame Day1 context')
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'tif');
    
        session = char('D1after_');
    figure;bar(All_D1after_Cell_Totalactive);
    title('Active frame Day1 after');
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'tif');
    
        session = char('D1post_');
    figure;bar(All_D1post_Cell_Totalactive)
    title('Active frame Day1 post');
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'tif');
    
% figure;histogram(All_D1pre_Cell_Totalactive);
% figure;histogram(All_D1ctx_Cell_Totalactive);
% figure;histogram(All_D1after_Cell_Totalactive);
% figure;histogram(All_D1post_Cell_Totalactive);


%     cd ..;%上のフォルダに移動

%% 3. cellular activity (SD, each sesion, All=TP+TPfP+TN), Day2 (copy Day1 and paste new sheet, and replace D1 (Day1) to D2 (Day2)

%  % set parameter
%     SD_thres = 3;%binaryの時のSDのthreshold
%     
%     % data setup (step 2)
%     
%     mkdir('3_cellular activity_SD') % フォルダの作成
%     cd ('3_cellular activity_SD');%今作成したフォルダに移動

    % basic infomation (mean/SD/z-score) of this data > binarization based on SD
    %D2pre. all
    df_D2pre_mean = mean(dFData_D2pre);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_D2pre_sd= std(dFData_D2pre);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_D2pre = (dFData_D2pre - df_D2pre_mean)./df_D2pre_sd;%各データの細胞ごとのz-score
    SDbinaryData_D2pre = sdData_D2pre>SD_thres;%binarize, 2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
    %D2ctx. all
    df_D2ctx_mean = mean(dFData_D2context);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_D2ctx_sd= std(dFData_D2context);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_D2ctx = (dFData_D2context - df_D2ctx_mean)./df_D2ctx_sd;%各データの細胞ごとのz-score
    SDbinaryData_D2ctx = sdData_D2ctx>SD_thres;%binarize, 2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
    %D2after. all
    df_D2after_mean = mean(dFData_D2after);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_D2after_sd= std(dFData_D2after);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_D2after = (dFData_D2after - df_D2after_mean)./df_D2after_sd;%各データの細胞ごとのz-score
    SDbinaryData_D2after = sdData_D2after>SD_thres;%binarize, 2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
    %D2post. all    
    df_D2post_mean = mean(dFData_D2post);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_D2post_sd= std(dFData_D2post);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
    sdData_D2post = (dFData_D2post - df_D2post_mean)./df_D2post_sd;%各データの細胞ごとのz-score
    SDbinaryData_D2post = sdData_D2post>SD_thres;%binarize, 2SD以上を1, それ未満を0に。s://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
    
        % save the binarized data for correlation matrix using SD-binarized data
        d_content = char('2SDbinar_');
        celltype = char('All_');
%         SDbinaryData_Time = [Data_imp(503:end,1),SDbinaryData_D2ctx];%add time
        %dlmwrite([mouseID,dinfo,'TP_2SDbinary_Time.txt'], SDbinaryData_TP_Time, 'delimiter','\t');%
        session = char('D2pre_');      
        writematrix(SDbinaryData_D2pre,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
        session = char('D2context_');      
        writematrix(SDbinaryData_D2ctx,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
        session = char('D2after_');      
        writematrix(SDbinaryData_D2after,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
        session = char('D2post_');      
        writematrix(SDbinaryData_D2post,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%

        %save raster plot figure
        session = char('D2context_');
        sesion_title = char('Day2 context');
        d_content = char('2SDbinar_raster_');
       
        celltype = char('All_');

        %Day2context, all
        figure;imagesc(SDbinaryData_D2ctx');
        cmap_raster_all=[1 1 1;0 0 0];%white(0), black(1)
        colormap(cmap_raster_all);
%         title('Day2 context')
        title(sesion_title);
        xlabel('Time (min)');ylabel('Cell IDs');
        xticks(0:1200:12000);xticklabels(0:1:10);
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

        celltype = char('TP_');

        %Day2context, Tom+
        figure;imagesc(SDbinaryData_D2ctx(:,1:TP_last)');
        cmap_raster_TP=[1 1 1;1 0 1];%white(0), magenta(1)
        colormap(cmap_raster_TP);
        title(sesion_title);
        xlabel('Time (min)');ylabel('Cell IDs');
        xticks(0:1200:12000);xticklabels(0:1:10);
        yticks(1:TP_Cell);xticklabels(1:1:TP_Cell);
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

        celltype = char('TN_');

        %Day2context, Tom-
        figure;imagesc(SDbinaryData_D2ctx(:,TN_first:end)');
        cmap_raster_TN=[1 1 1;0 1 1];%white(0), magenta(1)
        colormap(cmap_raster_TN);
        title(sesion_title);
        xlabel('Time (min)');ylabel('Cell IDs');
        xticks(0:1200:12000);xticklabels(0:1:10);
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

        celltype = char('fPTN_');
        %Day2context, fP+Tom-
        figure;imagesc(SDbinaryData_D2ctx(:,TPfP_first:end)');
        cmap_raster_TN=[1 1 1;.3 .3 .3];%white(0), magenta(1)
        colormap(cmap_raster_TN);
        title(sesion_title);
        xlabel('Time (min)');ylabel('Cell IDs');
        xticks(0:1200:12000);xticklabels(0:1:10);
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
        saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

        % the number of active frame
    % All_D2pre_binarySD__Max = max(SDbinaryData_D2pre);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
    All_D2pre_Cell_Totalactive = sum(SDbinaryData_D2pre, 1);%各細胞の活動量 (2SD以上のframe数)
    %         All_Cell_Totalactive_Mean = mean(All_Cell_Totalactive);%全細胞の平均活動量
    All_D2ctx_Cell_Totalactive = sum(SDbinaryData_D2ctx, 1);
    All_D2after_Cell_Totalactive = sum(SDbinaryData_D2after, 1);
    All_D2post_Cell_Totalactive = sum(SDbinaryData_D2post, 1);

    % save cellular activity for combining these data from all mice
%     dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
    
    celltype = char('All_');%All, TP, TN, TNrand, 
    d_content = char('ActiveFrame_');
       
    session = char('D2pre_');
    writematrix(All_D2pre_Cell_Totalactive,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
    session = char('D2context_');      
    writematrix(All_D2ctx_Cell_Totalactive,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
    session = char('D2after_');      
    writematrix(All_D2after_Cell_Totalactive,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
    session = char('D2post_');      
    writematrix(All_D2post_Cell_Totalactive,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
    
    %figure (the number of total active frames)
       
       d_content = char('ActiveFrame_');
       celltype = char('All_');
    
       session = char('D2pre_');

    figure;bar(All_D2pre_Cell_Totalactive)
    title('Active frame Day2 pre')
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'tif');
    
        session = char('D2context_');
    figure;bar(All_D2ctx_Cell_Totalactive)
    title('Active frame Day2 context')
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'tif');
    
        session = char('D2after_');
    figure;bar(All_D2after_Cell_Totalactive);
    title('Active frame Day2 after');
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'tif');
    
        session = char('D2post_');
    figure;bar(All_D2post_Cell_Totalactive)
    title('Active frame Day2 post');
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,dinfo,session,d_content,celltype,dtime], 'tif');
    
% figure;histogram(All_D2pre_Cell_Totalactive);
% figure;histogram(All_D2ctx_Cell_Totalactive);
% figure;histogram(All_D2after_Cell_Totalactive);
% figure;histogram(All_D2post_Cell_Totalactive);


%     cd ..;%上のフォルダに移動
%% 3-2 cellular activity event (SD, All/TPvsTN), Day1

% preparation for counting events
D1pre_1 = Day1pre_frame - 1;
D1ctx_1 = Day1context_frame - 1;
D1after_1 = Day1after_frame -1;
D1post_1 = Day1post_frame -1;
D2pre_1 = Day2pre_frame -1;
D2ctx_1 = Day2context_frame-1;
D2after_1 = Day2after_frame -1;
D2post_1 = Day2post_frame -1;


% active event (Shift > count start > save
    dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-
    celltype = char('All_');%All, TP, TN, TNrand, 
    d_content = char('TotalEvent_');
    
    session = char('D1pre_');
    
    %Day1 pre
    SDbinaryData_D1pre_ev = SDbinaryData_D1pre(2:end,:) - SDbinaryData_D1pre(1:D1pre_1,:);
    SDbinaryData_D1pre_evst = SDbinaryData_D1pre_ev > 0;
    All_D1pre_Cell_Totalev = sum(SDbinaryData_D1pre_evst, 1);%各細胞の活動回数
% save cellular activity for combining these data from all mice
writematrix(All_D1pre_Cell_Totalev,[mouseID,dinfo,session,d_content,celltype,dtime,'.txt'],'delimiter','\t');%
    
    session = char('D1context_');
    %Day1 context
    SDbinaryData_D1ctx_ev = SDbinaryData_D1ctx(2:end,:) - SDbinaryData_D1ctx(1:D1ctx_1,:);
    SDbinaryData_D1ctx_evst = SDbinaryData_D1ctx_ev > 0;
    All_D1ctx_Cell_Totalev = sum(SDbinaryData_D1ctx_evst, 1);%各細胞の活動回数
writematrix(All_D1ctx_Cell_Totalev,[mouseID,dinfo,session,d_content,celltype,dtime,'.txt'],'delimiter','\t');
    
    session = char('D1after_');
    %Day1 after
    SDbinaryData_D1after_ev = SDbinaryData_D1after(2:end,:) - SDbinaryData_D1after(1:D1after_1,:);
    SDbinaryData_D1after_evst = SDbinaryData_D1after_ev > 0;
    All_D1after_Cell_Totalev = sum(SDbinaryData_D1after_evst, 1);%各細胞の活動回数
writematrix(All_D1after_Cell_Totalev,[mouseID,dinfo,session,d_content,celltype,dtime,'.txt'],'delimiter','\t');%
    
    session = char('D1post_');
    %Day1 post
    SDbinaryData_D1post_ev = SDbinaryData_D1post(2:end,:) - SDbinaryData_D1post(1:D1post_1,:);
    SDbinaryData_D1post_evst = SDbinaryData_D1post_ev > 0;
    All_D1post_Cell_Totalev = sum(SDbinaryData_D1post_evst, 1);%各細胞の活動回数
writematrix(All_D1post_Cell_Totalev,[mouseID,dinfo,session,d_content,celltype,dtime,'.txt'],'delimiter','\t');%

    %figure (event)
    celltype = char('All_');%All, TP, TN, TNrand, 
    d_content = char('TotalEvent_');
    %Day1 pre
    session = char('D1pre_');
    figure;bar(All_D1pre_Cell_Totalev);
    title('Event Day1 pre')
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
   
    %Day1 context
    session = char('D1context_');
    figure;bar(All_D1ctx_Cell_Totalev);
    title('Event Day1 context')
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
    
    %Day1 after
    session = char('D1after_');
    figure;bar(All_D1after_Cell_Totalev);
    title('Event Day1 after')
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
    
    %Day1 post
    session = char('D1post_');
    figure;bar(All_D1post_Cell_Totalev);
    title('Event Day1 post')
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');


%figure (2次元plot active frame-event (TP-TN))
    d_content = char('ScatActiveFrEv_');
    %Day1 pre
    session = char('D1pre_');
    celltype = char('TPvfPTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D1pre_Cell_Totalactive(:,TPfP_first:end),All_D1pre_Cell_Totalev(:,TPfP_first:end),25,'MarkerEdgeColor',[.5 .5 .5]);
hold on 
scatter(All_D1pre_Cell_Totalactive(:,1:TP_last),All_D1pre_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day1 pre');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    celltype = char('TPvTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D1pre_Cell_Totalactive(:,TN_first:end),All_D1pre_Cell_Totalev(:,TN_first:end),25,'MarkerEdgeColor',"cyan");
hold on 
scatter(All_D1pre_Cell_Totalactive(:,1:TP_last),All_D1pre_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day1 pre');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

% figure;scatter(All_D1pre_Cell_Totalactive(:,TN_first:end),All_D1pre_Cell_Totalev(:,TN_first:end),"cyan");
% hold on 
% scatter(All_D1pre_Cell_Totalactive(:,TPfP_first:TPfP_last),All_D1pre_Cell_Totalev(:,TPfP_first:TPfP_last),"yellow");
% scatter(All_D1pre_Cell_Totalactive(:,1:TP_last),All_D1pre_Cell_Totalev(:,1:TP_last),"magenta");
% hold off
% title('Day1 pre');xlabel('Frame');ylabel('Event');
% saveas(gcf, [mouseID,'_D1pre_TPfPTN_ScatTotalactiveFrameEvent',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_D1pre_TPfPTN_ScatTotalactiveFrameEvent',dtime], 'tif');

    %Day1 context
    session = char('D1context_');
    celltype = char('TPvfPTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D1ctx_Cell_Totalactive(:,TPfP_first:end),All_D1ctx_Cell_Totalev(:,TPfP_first:end),25,'MarkerEdgeColor',[.5 .5 .5]);
hold on 
scatter(All_D1ctx_Cell_Totalactive(:,1:TP_last),All_D1ctx_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day1 context');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

celltype = char('TPvTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D1ctx_Cell_Totalactive(:,TN_first:end),All_D1ctx_Cell_Totalev(:,TN_first:end),25,"cyan");
hold on 
scatter(All_D1ctx_Cell_Totalactive(:,1:TP_last),All_D1ctx_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day1 context');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

%Day1 after
    session = char('D1after_');
    celltype = char('TPvfPTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D1after_Cell_Totalactive(:,TPfP_first:end),All_D1after_Cell_Totalev(:,TPfP_first:end),25,'MarkerEdgeColor',[.5 .5 .5]);
hold on 
scatter(All_D1after_Cell_Totalactive(:,1:TP_last),All_D1after_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day1 after');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

celltype = char('TPvTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D1after_Cell_Totalactive(:,TN_first:end),All_D1after_Cell_Totalev(:,TN_first:end),25,'MarkerEdgeColor',"cyan");
hold on 
scatter(All_D1after_Cell_Totalactive(:,1:TP_last),All_D1after_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day1 after');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

%Day1 post
    session = char('D1post_');
    celltype = char('TPvfPTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
%figure;scatter(All_D1pre_Cell_Totalactive,All_D1pre_Cell_Totalev);
figure;scatter(All_D1post_Cell_Totalactive(:,TPfP_first:end),All_D1post_Cell_Totalev(:,TPfP_first:end),25,'MarkerEdgeColor',[.5 .5 .5]);
hold on 
scatter(All_D1post_Cell_Totalactive(:,1:TP_last),All_D1post_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day1 post');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

celltype = char('TPvTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D1post_Cell_Totalactive(:,TN_first:end),All_D1post_Cell_Totalev(:,TN_first:end),25,'MarkerEdgeColor',"cyan");
hold on 
scatter(All_D1post_Cell_Totalactive(:,1:TP_last),All_D1post_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day1 post');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

%  cd ..;%上のフォルダに移動

%% cellular activity event (SD, All/TPvsTN), Day2 (copy Day1 and paste new sheet, and replace D1 (Day1) to D2 (Day2)

% preparation for counting events

% active event (Shift > count start > save
%     dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%     celltype = char('All_');%All, TP, TN, TNrand, 
    d_content = char('TotalEvent_');
    
    session = char('D2pre_');
    
    %Day2 pre
    SDbinaryData_D2pre_ev = SDbinaryData_D2pre(2:end,:) - SDbinaryData_D2pre(1:D2pre_1,:);
    SDbinaryData_D2pre_evst = SDbinaryData_D2pre_ev > 0;
    All_D2pre_Cell_Totalev = sum(SDbinaryData_D2pre_evst, 1);%各細胞の活動回数
% save cellular activity for combining these data from all mice
writematrix(All_D2pre_Cell_Totalev,[mouseID,dinfo,session,d_content,celltype,dtime,'.txt'],'delimiter','\t');%
    
    session = char('D2context_');
    %Day2 context
    SDbinaryData_D2ctx_ev = SDbinaryData_D2ctx(2:end,:) - SDbinaryData_D2ctx(1:D2ctx_1,:);
    SDbinaryData_D2ctx_evst = SDbinaryData_D2ctx_ev > 0;
    All_D2ctx_Cell_Totalev = sum(SDbinaryData_D2ctx_evst, 1);%各細胞の活動回数
writematrix(All_D2ctx_Cell_Totalev,[mouseID,dinfo,session,d_content,celltype,dtime,'.txt'],'delimiter','\t');
    
    session = char('D2after_');
    %Day2 after
    SDbinaryData_D2after_ev = SDbinaryData_D2after(2:end,:) - SDbinaryData_D2after(1:D2after_1,:);
    SDbinaryData_D2after_evst = SDbinaryData_D2after_ev > 0;
    All_D2after_Cell_Totalev = sum(SDbinaryData_D2after_evst, 1);%各細胞の活動回数
writematrix(All_D2after_Cell_Totalev,[mouseID,dinfo,session,d_content,celltype,dtime,'.txt'],'delimiter','\t');%
    
    session = char('D2post_');
    %Day2 post
    SDbinaryData_D2post_ev = SDbinaryData_D2post(2:end,:) - SDbinaryData_D2post(1:D2post_1,:);
    SDbinaryData_D2post_evst = SDbinaryData_D2post_ev > 0;
    All_D2post_Cell_Totalev = sum(SDbinaryData_D2post_evst, 1);%各細胞の活動回数
writematrix(All_D2post_Cell_Totalev,[mouseID,dinfo,session,d_content,celltype,dtime,'.txt'],'delimiter','\t');%

    %figure (event)
    celltype = char('All_');%All, TP, TN, TNrand, 
    d_content = char('TotalEvent_');
    %Day2 pre
    session = char('D2pre_');
    figure;bar(All_D2pre_Cell_Totalev);
    title('Event Day2 pre')
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
   
    %Day2 context
    session = char('D2context_');
    figure;bar(All_D2ctx_Cell_Totalev);
    title('Event Day2 context')
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
    
    %Day2 after
    session = char('D2after_');
    figure;bar(All_D2after_Cell_Totalev);
    title('Event Day2 after')
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
    
    %Day2 post
    session = char('D2post_');
    figure;bar(All_D2post_Cell_Totalev);
    title('Event Day2 post')
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
    saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');


%figure (2次元plot active frame-event (TP-TN))
    d_content = char('ScatActiveFrEv_');
    %Day2 pre
    session = char('D2pre_');
    celltype = char('TPvfPTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D2pre_Cell_Totalactive(:,TPfP_first:end),All_D2pre_Cell_Totalev(:,TPfP_first:end),25,'MarkerEdgeColor',[.5 .5 .5]);
hold on 
scatter(All_D2pre_Cell_Totalactive(:,1:TP_last),All_D2pre_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day2 pre');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    celltype = char('TPvTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D2pre_Cell_Totalactive(:,TN_first:end),All_D2pre_Cell_Totalev(:,TN_first:end),25,'MarkerEdgeColor',"cyan");
hold on 
scatter(All_D2pre_Cell_Totalactive(:,1:TP_last),All_D2pre_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day2 pre');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

% figure;scatter(All_D2pre_Cell_Totalactive(:,TN_first:end),All_D2pre_Cell_Totalev(:,TN_first:end),"cyan");
% hold on 
% scatter(All_D2pre_Cell_Totalactive(:,TPfP_first:TPfP_last),All_D2pre_Cell_Totalev(:,TPfP_first:TPfP_last),"yellow");
% scatter(All_D2pre_Cell_Totalactive(:,1:TP_last),All_D2pre_Cell_Totalev(:,1:TP_last),"magenta");
% hold off
% title('Day2 pre');xlabel('Frame');ylabel('Event');
% saveas(gcf, [mouseID,'_D2pre_TPfPTN_ScatTotalactiveFrameEvent',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_D2pre_TPfPTN_ScatTotalactiveFrameEvent',dtime], 'tif');

    %Day2 context
    session = char('D2context_');
    celltype = char('TPvfPTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D2ctx_Cell_Totalactive(:,TPfP_first:end),All_D2ctx_Cell_Totalev(:,TPfP_first:end),25,'MarkerEdgeColor',[.5 .5 .5]);
hold on 
scatter(All_D2ctx_Cell_Totalactive(:,1:TP_last),All_D2ctx_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day2 context');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

celltype = char('TPvTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D2ctx_Cell_Totalactive(:,TN_first:end),All_D2ctx_Cell_Totalev(:,TN_first:end),25,"cyan");
hold on 
scatter(All_D2ctx_Cell_Totalactive(:,1:TP_last),All_D2ctx_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day2 context');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

%Day2 after
    session = char('D2after_');
    celltype = char('TPvfPTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D2after_Cell_Totalactive(:,TPfP_first:end),All_D2after_Cell_Totalev(:,TPfP_first:end),25,'MarkerEdgeColor',[.5 .5 .5]);
hold on 
scatter(All_D2after_Cell_Totalactive(:,1:TP_last),All_D2after_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day2 after');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

celltype = char('TPvTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D2after_Cell_Totalactive(:,TN_first:end),All_D2after_Cell_Totalev(:,TN_first:end),25,'MarkerEdgeColor',"cyan");
hold on 
scatter(All_D2after_Cell_Totalactive(:,1:TP_last),All_D2after_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day2 after');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

%Day2 post
    session = char('D2post_');
    celltype = char('TPvfPTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
%figure;scatter(All_D2pre_Cell_Totalactive,All_D2pre_Cell_Totalev);
figure;scatter(All_D2post_Cell_Totalactive(:,TPfP_first:end),All_D2post_Cell_Totalev(:,TPfP_first:end),25,'MarkerEdgeColor',[.5 .5 .5]);
hold on 
scatter(All_D2post_Cell_Totalactive(:,1:TP_last),All_D2post_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day2 post');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

celltype = char('TPvTN_');%All, TP(magenta[1 0 1]), TN(cyan[0 1 1]), TNrand(cyan[0 1 1]), fPTN(grey[.5 .5 .5]) 
figure;scatter(All_D2post_Cell_Totalactive(:,TN_first:end),All_D2post_Cell_Totalev(:,TN_first:end),25,'MarkerEdgeColor',"cyan");
hold on 
scatter(All_D2post_Cell_Totalactive(:,1:TP_last),All_D2post_Cell_Totalev(:,1:TP_last),25,"magenta");
hold off
title('Day2 post');xlabel('Time (Frame)');ylabel('Events');
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

 cd ..;%上のフォルダに移動

%% 4. synchronous activity (SD, All), Day1
%この解析は、synchronous activityの補助的なもの
    mkdir('4_synchronous2SD') % フォルダの作成
    cd ('4_synchronous2SD');%今作成したフォルダに移動

    % Synchronous activity (results from binarized df matrix)
    %Day1 pre, All
    All_D1pre_binarySD_SyncCells = sum(SDbinaryData_D1pre, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
        NumCell_dF_All_D1pre =size(dFData_D1pre, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1pre = size(dFData_D1pre,1);%Time(frame)数
    D1pre_SyncCells_Ratio = All_D1pre_binarySD_SyncCells./NumCell_dF_All_D1pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
    
    %Day1 context, All
    All_D1ctx_binarySD_SyncCells= sum(SDbinaryData_D1ctx, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%         NumCell_dF_All_D1pre =size(dFData_D1pre, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1ctx = size(dFData_D1context,1);%Time(frame)数
    D1ctx_SyncCells_Ratio = All_D1pre_binarySD_SyncCells./NumCell_dF_All_D1pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    %Day1 after, All
    All_D1after_binarySD_SyncCells= sum(SDbinaryData_D1after, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%         NumCell_dF_All_D1pre =size(dFData_D1pre, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1after = size(dFData_D1after,1);%Time(frame)数
    D1after_SyncCells_Ratio = All_D1after_binarySD_SyncCells./NumCell_dF_All_D1pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    %Day1 post, All
    All_D1post_binarySD_SyncCells= sum(SDbinaryData_D1post, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%         NumCell_dF_All_D1pre =size(dFData_D1pre, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1post = size(dFData_D1post,1);%Time(frame)数
    D1post_SyncCells_Ratio = All_D1post_binarySD_SyncCells./NumCell_dF_All_D1pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    % SD synchronize (TP vs TN)
    %Day1 pre, TP/TN
    TP_D1pre_binarySD_SyncCells= sum(SDbinaryData_D1pre(:,1:TP_last), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP_D1pre =size(dFData_D1pre(:,1:TP_last), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1pre = size(dFData_D1pre,1);%Time(frame)数
    TP_D1pre_SyncCells_Ratio = TP_D1pre_binarySD_SyncCells./NumCell_dF_TP_D1pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
        
    TN_D1pre_binarySD_SyncCells= sum(SDbinaryData_D1pre(:,TN_first:end), 2);
    NumCell_dF_TN_D1pre =size(dFData_D1pre(:,TN_first:end), 2);
        %         Timeframe_dF_D1pre = size(dFData_D1pre,1);%Time(frame)数
    TN_D1pre_SyncCells_Ratio = TN_D1pre_binarySD_SyncCells./NumCell_dF_TN_D1pre;
        
    %Day1 context, TP/TN
    TP_D1ctx_binarySD_SyncCells= sum(SDbinaryData_D1ctx(:,1:TP_last), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP_D1ctx =size(dFData_D1context(:,1:TP_last), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1ctx = size(dFData_D1context,1);%Time(frame)数
    TP_D1ctx_SyncCells_Ratio = TP_D1ctx_binarySD_SyncCells./NumCell_dF_TP_D1ctx;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
        
    TN_D1ctx_binarySD_SyncCells= sum(SDbinaryData_D1ctx(:,TN_first:end), 2);
    NumCell_dF_TN_D1ctx =size(dFData_D1context(:,TN_first:end), 2);
        %         Timeframe_dF_D1ctx = size(dFData_D1ctx,1);
    TN_D1ctx_SyncCells_Ratio = TN_D1ctx_binarySD_SyncCells./NumCell_dF_TN_D1ctx;

    %Day1 after, TP/TN
    TP_D1after_binarySD_SyncCells= sum(SDbinaryData_D1after(:,1:TP_last), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP_D1after =size(dFData_D1after(:,1:TP_last), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1after = size(dFData_D1after,1);%Time(frame)数
    TP_D1after_SyncCells_Ratio = TP_D1after_binarySD_SyncCells./NumCell_dF_TP_D1after;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
        
    TN_D1after_binarySD_SyncCells= sum(SDbinaryData_D1after(:,TN_first:end), 2);
    NumCell_dF_TN_D1after =size(dFData_D1after(:,TN_first:end), 2);
        %         Timeframe_dF_D1after = size(dFData_D1after,1);
    TN_D1after_SyncCells_Ratio = TN_D1after_binarySD_SyncCells./NumCell_dF_TN_D1after;

%Day1 post, TP/TN
    TP_D1post_binarySD_SyncCells= sum(SDbinaryData_D1post(:,1:TP_last), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP_D1post =size(dFData_D1post(:,1:TP_last), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D1post = size(dFData_D1post,1);%Time(frame)数
    TP_D1post_SyncCells_Ratio = TP_D1post_binarySD_SyncCells./NumCell_dF_TP_D1post;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
        
    TN_D1post_binarySD_SyncCells= sum(SDbinaryData_D1post(:,TN_first:end), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TN_D1post =size(dFData_D1post(:,TN_first:end), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        %         Timeframe_dF_D1post = size(dFData_D1post,1);%Time(frame)数
    TN_D1post_SyncCells_Ratio = TN_D1post_binarySD_SyncCells./NumCell_dF_TN_D1post;%同期活動している細胞の割合=同期活動している細胞数/全細胞数


% figure (synchronous active cells and time) and save
    dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-
    
    celltype = char('AllTPTN_');
    d_content = char('SynccellFrame_');

    %Day1 pre
    session = char('D1pre_');
figure;plot(All_D1pre_binarySD_SyncCells,"Color",[.5 .5 .5]);
hold on;
figure;plot(TP_D1pre_binarySD_SyncCells,"Color",'magenta');
plot(TN_D1pre_binarySD_SyncCells,"Color",'cyan');
hold off;
ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day1 pre','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

% figure;plot(TP_D1pre_SyncCells_Ratio,"Color",'magenta');
% hold on;
% plot(TN_D1pre_SyncCells_Ratio,"Color",'cyan');
% plot(D1pre_SyncCells_Ratio,"Color",[.5 .5 .5]);
% hold off;

    %Day1 context
    session = char('D1context_');
figure;plot(All_D1ctx_binarySD_SyncCells,"Color",[.5 .5 .5]);
hold on;
figure;plot(TP_D1ctx_binarySD_SyncCells,"Color",'magenta');
plot(TN_D1ctx_binarySD_SyncCells,"Color",'cyan');
ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day1 context','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day1 after
    session = char('D1after_');
figure;plot(All_D1after_binarySD_SyncCells,"Color",[.5 .5 .5]);
hold on;
figure;plot(TP_D1after_binarySD_SyncCells,"Color",'magenta');
plot(TN_D1after_binarySD_SyncCells,"Color",'cyan');
ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day1 after','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day1 post
    session = char('D1post_');
figure;plot(All_D1post_binarySD_SyncCells,"Color",[.5 .5 .5]);
hold on;
figure;plot(TP_D1post_binarySD_SyncCells,"Color",'magenta');
plot(TN_D1post_binarySD_SyncCells,"Color",'cyan');
ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day1 post','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');


% figure (histogram, synchronous active cells)
    celltype = char('All_');
    d_content = char('HistSynccellFrame_');

    %Day1 pre
    session = char('D1pre_');
figure;histogram(All_D1pre_binarySD_SyncCells,"FaceColor",[.5 .5 .5]);
xlabel('# Synchronous active cells');ylabel('# Frames');title('Day1 pre','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day1 context
    session = char('D1context_');
figure;histogram(All_D1ctx_binarySD_SyncCells,"FaceColor",[.5 .5 .5]);
xlabel('# Synchronous active cells');ylabel('# Frames');title('Day1 context','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day1 after
    session = char('D1after_');
figure;histogram(All_D1after_binarySD_SyncCells,"FaceColor",[.5 .5 .5]);
xlabel('# Synchronous active cells');ylabel('# Frames');title('Day1 after','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day1 post
    session = char('D1post_');
figure;histogram(All_D1post_binarySD_SyncCells,"FaceColor",[.5 .5 .5]);
xlabel('# Synchronous active cells');ylabel('# Frames');title('Day1 post','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

% cd ..;%上のフォルダに移動

%% 4. synchronous activity (SD, All), Day2 (copy Day1 and paste new sheet, and replace D1 (Day1) to D2 (Day2)
%この解析は、synchronous activityの補助的なもの
%     mkdir('4_synchronous activity_SD') % フォルダの作成
%     cd ('4_synchronous activity_SD');%今作成したフォルダに移動

    % Synchronous activity (results from binarized df matrix)
    %Day2 pre, All
    All_D2pre_binarySD_SyncCells = sum(SDbinaryData_D2pre, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
        NumCell_dF_All_D2pre =size(dFData_D2pre, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D2pre = size(dFData_D2pre,1);%Time(frame)数
    D2pre_SyncCells_Ratio = All_D2pre_binarySD_SyncCells./NumCell_dF_All_D2pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
    
    %Day2 context, All
    All_D2ctx_binarySD_SyncCells= sum(SDbinaryData_D2ctx, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%         NumCell_dF_All_D2pre =size(dFData_D2pre, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D2ctx = size(dFData_D2context,1);%Time(frame)数
    D2ctx_SyncCells_Ratio = All_D2pre_binarySD_SyncCells./NumCell_dF_All_D2pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    %Day2 after, All
    All_D2after_binarySD_SyncCells= sum(SDbinaryData_D2after, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%         NumCell_dF_All_D2pre =size(dFData_D2pre, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D2after = size(dFData_D2after,1);%Time(frame)数
    D2after_SyncCells_Ratio = All_D2after_binarySD_SyncCells./NumCell_dF_All_D2pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    %Day2 post, All
    All_D2post_binarySD_SyncCells= sum(SDbinaryData_D2post, 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%         NumCell_dF_All_D2pre =size(dFData_D2pre, 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D2post = size(dFData_D2post,1);%Time(frame)数
    D2post_SyncCells_Ratio = All_D2post_binarySD_SyncCells./NumCell_dF_All_D2pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    % SD synchronize (TP vs TN)
    %Day2 pre, TP/TN
    TP_D2pre_binarySD_SyncCells= sum(SDbinaryData_D2pre(:,1:TP_last), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP_D2pre =size(dFData_D2pre(:,1:TP_last), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D2pre = size(dFData_D2pre,1);%Time(frame)数
    TP_D2pre_SyncCells_Ratio = TP_D2pre_binarySD_SyncCells./NumCell_dF_TP_D2pre;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
        
    TN_D2pre_binarySD_SyncCells= sum(SDbinaryData_D2pre(:,TN_first:end), 2);
    NumCell_dF_TN_D2pre =size(dFData_D2pre(:,TN_first:end), 2);
        %         Timeframe_dF_D2pre = size(dFData_D2pre,1);%Time(frame)数
    TN_D2pre_SyncCells_Ratio = TN_D2pre_binarySD_SyncCells./NumCell_dF_TN_D2pre;
        
    %Day2 context, TP/TN
    TP_D2ctx_binarySD_SyncCells= sum(SDbinaryData_D2ctx(:,1:TP_last), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP_D2ctx =size(dFData_D2context(:,1:TP_last), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D2ctx = size(dFData_D2context,1);%Time(frame)数
    TP_D2ctx_SyncCells_Ratio = TP_D2ctx_binarySD_SyncCells./NumCell_dF_TP_D2ctx;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
        
    TN_D2ctx_binarySD_SyncCells= sum(SDbinaryData_D2ctx(:,TN_first:end), 2);
    NumCell_dF_TN_D2ctx =size(dFData_D2context(:,TN_first:end), 2);
        %         Timeframe_dF_D2ctx = size(dFData_D2ctx,1);
    TN_D2ctx_SyncCells_Ratio = TN_D2ctx_binarySD_SyncCells./NumCell_dF_TN_D2ctx;

    %Day2 after, TP/TN
    TP_D2after_binarySD_SyncCells= sum(SDbinaryData_D2after(:,1:TP_last), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP_D2after =size(dFData_D2after(:,1:TP_last), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D2after = size(dFData_D2after,1);%Time(frame)数
    TP_D2after_SyncCells_Ratio = TP_D2after_binarySD_SyncCells./NumCell_dF_TP_D2after;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
        
    TN_D2after_binarySD_SyncCells= sum(SDbinaryData_D2after(:,TN_first:end), 2);
    NumCell_dF_TN_D2after =size(dFData_D2after(:,TN_first:end), 2);
        %         Timeframe_dF_D2after = size(dFData_D2after,1);
    TN_D2after_SyncCells_Ratio = TN_D2after_binarySD_SyncCells./NumCell_dF_TN_D2after;

%Day2 post, TP/TN
    TP_D2post_binarySD_SyncCells= sum(SDbinaryData_D2post(:,1:TP_last), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TP_D2post =size(dFData_D2post(:,1:TP_last), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        Timeframe_dF_D2post = size(dFData_D2post,1);%Time(frame)数
    TP_D2post_SyncCells_Ratio = TP_D2post_binarySD_SyncCells./NumCell_dF_TP_D2post;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
        
    TN_D2post_binarySD_SyncCells= sum(SDbinaryData_D2post(:,TN_first:end), 2);%dFのSDで２値化したmatrixで、簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    NumCell_dF_TN_D2post =size(dFData_D2post(:,TN_first:end), 2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
        %         Timeframe_dF_D2post = size(dFData_D2post,1);%Time(frame)数
    TN_D2post_SyncCells_Ratio = TN_D2post_binarySD_SyncCells./NumCell_dF_TN_D2post;%同期活動している細胞の割合=同期活動している細胞数/全細胞数


% figure (synchronous active cells and time)
%     dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-

    celltype = char('AllTPTN_');
    d_content = char('SynccellFrame_');

    %Day2 pre
    session = char('D2pre_');
figure;plot(All_D2pre_binarySD_SyncCells,"Color",[.5 .5 .5]);
hold on;
figure;plot(TP_D2pre_binarySD_SyncCells,"Color",'magenta');
plot(TN_D2pre_binarySD_SyncCells,"Color",'cyan');
hold off;
ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day2 pre','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

% figure;plot(TP_D2pre_SyncCells_Ratio,"Color",'magenta');
% hold on;
% plot(TN_D2pre_SyncCells_Ratio,"Color",'cyan');
% plot(D2pre_SyncCells_Ratio,"Color",[.5 .5 .5]);
% hold off;

    %Day2 context
    session = char('D2context_');
figure;plot(All_D2ctx_binarySD_SyncCells,"Color",[.5 .5 .5]);
hold on;
figure;plot(TP_D2ctx_binarySD_SyncCells,"Color",'magenta');
plot(TN_D2ctx_binarySD_SyncCells,"Color",'cyan');
ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day2 context','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day2 after
    session = char('D2after_');
figure;plot(All_D2after_binarySD_SyncCells,"Color",[.5 .5 .5]);
hold on;
figure;plot(TP_D2after_binarySD_SyncCells,"Color",'magenta');
plot(TN_D2after_binarySD_SyncCells,"Color",'cyan');
ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day2 after','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day2 post
    session = char('D2post_');
figure;plot(All_D2post_binarySD_SyncCells,"Color",[.5 .5 .5]);
hold on;
figure;plot(TP_D2post_binarySD_SyncCells,"Color",'magenta');
plot(TN_D2post_binarySD_SyncCells,"Color",'cyan');
ylabel('# Synchronous active cells');xlabel('Time (Frame)');title('Day2 post','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');


% figure (histogram, synchronous active cells)
    celltype = char('All_');
    d_content = char('HistSynccellFrame_');

    %Day2 pre
    session = char('D2pre_');
figure;histogram(All_D2pre_binarySD_SyncCells,"FaceColor",[.5 .5 .5]);
xlabel('# Synchronous active cells');ylabel('# Frames');title('Day2 pre','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day2 context
    session = char('D2context_');
figure;histogram(All_D2ctx_binarySD_SyncCells,"FaceColor",[.5 .5 .5]);
xlabel('# Synchronous active cells');ylabel('# Frames');title('Day2 context','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day2 after
    session = char('D2after_');
figure;histogram(All_D2after_binarySD_SyncCells,"FaceColor",[.5 .5 .5]);
xlabel('# Synchronous active cells');ylabel('# Frames');title('Day2 after','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

    %Day2 post
    session = char('D2post_');
figure;histogram(All_D2post_binarySD_SyncCells,"FaceColor",[.5 .5 .5]);
xlabel('# Synchronous active cells');ylabel('# Frames');title('Day2 post','FontSize',12);
saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');

cd ..;%上のフォルダに移動

    %% 5. Correlation matrix (All, for TP_all vs TN all), Day1

    mkdir('5_dFcorrelation') % フォルダの作成
%     cd ('5_dFcorrelation');%今作成したフォルダに移動
    
    %CorrMat (Parameter set > CorrMat > Figure)
% Correlation matrix from Alan's(original name = out_matrix.m)
% How to use: CorrMat_Alan(binned_datamatrix(time, neurons),shift_window(4))
shift_cormat = 5;%parameter setupに移動
window_cormat = 99;%4;
% TP_CorrMat_Alan = CorrMat_Alan(dFData_TP,shift_cormat);
% writematrix(TP_CorrMat_Alan, [mouseID,dinfo,'_TPshuf2_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
% TP_CorrMat_Alan_sum = sum(TP_CorrMat_Alan);
% TP_CorrMat_sum_sum = sum(TP_CorrMat_Alan_sum);
% Ave_List_TP_CorrMat_sum_sum = mean(TP_CorrMat_sum_sum);%Shuffleの名残
% SD_List_TP_CorrMat_sum_sum = std(TP_CorrMat_sum_sum);
% SEM_List_TP_CorrMat_sum_sum = SD_List_TP_CorrMat_sum_sum./sqrt(iter_TPshuf2);%
% % writematrix(TP_CorrMat_sum_sum, [mouseID,dinfo,'_TP_',num2str(n),'_CorrMat_SumSum.txt'], 'delimiter','\t');%

% figure;
% imagesc(TP_CorrMat_Alan);
% dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
% %一応もう一度時間をとる
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');

dFData_D1ctx_TP = dFData_D1context(:,1:TP_last);%221024
dFData_D1ctx_All = dFData_D1context;%221024
dFData_D1ctx_TN = dFData_D1context(:,TN_first:end);%221024

%TP
tic;
TP_D1ctx_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_D1ctx_TP,window_cormat,shift_cormat);
toc
celltype = char('TP_');
d_content = char('CorMatSkip_');
writematrix(TP_D1ctx_CorrMat_skip, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat),'.txt'], 'delimiter','\t');%

TP_D1ctx_CorrMat_skip_sum = sum(TP_D1ctx_CorrMat_skip);
TP_D1ctx_CorrMat_skip_sumsum = sum(sum(TP_D1ctx_CorrMat_skip));

writematrix(TP_D1ctx_CorrMat_skip_sum, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat),'.txt'], 'delimiter','\t');%

dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-

figure;
imagesc(TP_D1ctx_CorrMat_skip);
colorbar;caxis([0 0.01]);
saveas(gcf, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat)], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat)], 'tif');

%All
tic;
All_D1ctx_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_D1ctx_All,window_cormat,shift_cormat);
toc
celltype = char('All_');
d_content = char('CorMatSkip_');
writematrix(All_D1ctx_CorrMat_skip, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat),'.txt'], 'delimiter','\t');%

All_D1ctx_CorrMat_skip_sum = sum(All_D1ctx_CorrMat_skip,'omitnan');
All_D1ctx_CorrMat_skip_sumsum = sum(sum(All_D1ctx_CorrMat_skip,'omitnan'));

writematrix(All_D1ctx_CorrMat_skip_sum, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat),'.txt'], 'delimiter','\t');%

figure;
imagesc(All_D1ctx_CorrMat_skip);
colorbar;caxis([0 0.01]);
saveas(gcf, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat)], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat)], 'tif');

%TN
tic;
TN_D1ctx_CorrMat_skip = CorrMat_Alan_As220506_skip(dFData_D1ctx_TN,window_cormat,shift_cormat);
toc
celltype = char('TN_');
d_content = char('CorMatSkip_');
writematrix(TN_D1ctx_CorrMat_skip, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat),'.txt'], 'delimiter','\t');%

TN_D1ctx_CorrMat_skip_sum = sum(TN_D1ctx_CorrMat_skip,'omitnan');
TN_D1ctx_CorrMat_skip_sumsum = sum(sum(TN_D1ctx_CorrMat_skip,'omitnan'));

writematrix(TN_D1ctx_CorrMat_skip_sum, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat),'.txt'], 'delimiter','\t');%

figure;
imagesc(TN_D1ctx_CorrMat_skip);
colorbar;caxis([0 0.01]);
saveas(gcf, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat)], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
saveas(gcf, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat)], 'tif');

% summary of results
Cell = [{'Tom+'};{'Tom-'};{'All'}];
List_D1ctx_CorrMat_skip_sumsum = [TP_D1ctx_CorrMat_skip_sumsum;TN_D1ctx_CorrMat_skip_sumsum;All_D1ctx_CorrMat_skip_sumsum];
Res_D1ctx_CorrMat_skip_sumsum = table(Cell,List_D1ctx_CorrMat_skip_sumsum);



% % %並列処理を試したが、計算時間は短縮されなかった。またはエラーのため、使わない
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

% saveas(gcf, [mouseID,'_CorrMat_As',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% saveas(gcf, [mouseID,'_CorrMat_As',dtime], 'tif');


%% 5. Correlation matrix (TP shift, for TP_all vs TN all), Day1

%%% 4-6 (4-2 circshift).　Create Shuffled data (Shuffle type 5, raster(binary) across time, shuffle in each cells)
%細胞ごとにshiftしているので、mean/SDは変わらない
% Preparation/setup >> Create random IDs for Shuffle > random data & save > cal > summarize >>
% save
%この計算で比較する前に、全細胞の活動頻度がそれほど高くないことを確認した方がいい (あまり活動頻度が高いと意味がない)

% Folders for saving the shifted data
ftime = char(datetime('now','Format','HHmmss'));%フォルダ名, (データ名)
mkdir([mouseID,dinfo,'_TPShift_',ftime]) % フォルダの作成
% % cd ([mouseID,dinfo,'_shuf2_',ftime]);%今作成したフォルダに移動
% mkdir([mouseID,dinfo,'_TPShift_time',ftime]) % フォルダの作成

%set parameters
iter_TPShift = 100;%5000-10000;%1000:p=0.01を検出するためには最低でもこれくらい必要なのでは？ Num_Time_frames_TP/10;Num_Cell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)
% NumCell_dF_TP =size(dFData_TP, 2);
% for cal (for...end内)
% SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
%for correlation matrix
% shift_cormat = 19;

% SD_thres = 2;%221024

% % % data setup
dFData_D1ctx_TP = dFData_D1context(:,1:TP_last);%221024
dFData_D1ctx_All = dFData_D1context;%221024
dFData_D1ctx_TN = dFData_D1context(:,TN_first:end);%221024

% prepare matrix for shifted data and IDs
NumCell_dF_TP = TP_Cell;%221024
dFData_D1ctx_TP = dFData_D1context(:,1:TP_last);%221024
TPShift_IDs = zeros(iter_TPShift,NumCell_dF_TP);%shiftの履歴として一応残しておく
dfData_TPShift_D1ctx = dFData_D1ctx_TP;%zeros(Timeframe_dF_TP,NumCell_dF_TP);

    NumCell_TPShift = size(dfData_TPShift_D1ctx,2);%=size(dFData_TP,2);%= Num_TomPosi;%細胞数, https://jp.mathworks.com/help/matlab/ref/size.html?searchHighlight=size&s_tid=srchtitle_size_1
    Timeframe_dF_TPShift_D1ctx = size(dfData_TPShift_D1ctx,1);%=size(dFData_TP,1);%Time(frame)数%Time(frame)数
    Timeframe_dF_TP_D1ctx = Timeframe_dF_TPShift_D1ctx;
% 4-2-. Matrix for summarizing repeated data
List_df_TPShift_D1ctx_mean = zeros(iter_TPShift,NumCell_TPShift);
List_TPShift_D1ctx_binarySD_SyncCells = zeros(iter_TPShift,Timeframe_dF_TPShift_D1ctx);%for..endの開始前に配置
List_TPShift_D1ctx_SyncCells_Ratio = zeros(iter_TPShift,Timeframe_dF_TPShift_D1ctx);%for..endの開始前に配置
List_TPShift_D1ctx_binarySD_Max = zeros(iter_TPShift,NumCell_TPShift);
List_TPShift_D1ctx_Cell_Totalactive = zeros(iter_TPShift,NumCell_TPShift);

Size_TP_D1ctx_CorMatSkip_sum = size(TP_D1ctx_CorrMat_skip_sum,2);
List_TPShift_D1ctx_CorMatSkip_sum = zeros(iter_TPShift,Size_TP_D1ctx_CorMatSkip_sum);

% create shifted data > cal. > list
for n = 1:iter_TPShift%shift回数
dfData_TPShift_D1ctx = dFData_D1ctx_TP;%re-set;

    % Create random IDs for Shuffle in each cells (repeat for total cells)
    for i = 1:NumCell_TPShift%shift the time of all cells one by one
        random_TPshift_D1ctx = randperm(Timeframe_dF_TPShift_D1ctx,1);%randperm(%TPの時間,1);
        dfData_TPShift_D1ctx(:,i)= circshift(dfData_TPShift_D1ctx(:,i),random_TPshift_D1ctx);%create shifted data
        TPShift_IDs(n,i) = random_TPshift_D1ctx;%shiftの履歴として一応残しておく
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
    df_TPShift_D1ctx_mean = mean(dfData_TPShift_D1ctx);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
    df_TPShift_D1ctx_sd= std(dfData_TPShift_D1ctx);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1

    sdData_TPShift_D1ctx = (dfData_TPShift_D1ctx - df_TPShift_D1ctx_mean)./df_TPShift_D1ctx_sd;%各データの細胞ごとのSD
    SDbinaryData_TPShift_D1ctx = sdData_TPShift_D1ctx>SD_thres;%3SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1

% %         % save the binarized data for correlation matrix using SD-binarized data
% %         SDbinaryData_TPShift6_Time = [Data_imp(503:end,1),SDbinaryData_TPShift6];%add time
% %         % dlmwrite([mouseID,dinfo,'TPShift6_3SDbinary_Time.txt'], SDbinaryData_TPShift6_Time, 'delimiter','\t');%
% %         writematrix([mouseID,dinfo,'TPShift6_SDbinary_Time.txt'], SDbinaryData_TPShift6_Time, 'delimiter','\t');%

    % results from binarized df matrix
    TPShift_D1ctx_binarySD_SyncCells= sum(SDbinaryData_TPShift_D1ctx, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
    %細胞数、frame数は既に取得済み
    TPShift_D1ctx_SyncCells_Ratio = TPShift_D1ctx_binarySD_SyncCells./NumCell_TPShift;%同期活動している細胞の割合=同期活動している細胞数/全細胞数

    TPShift_D1ctx_binarySD_Max = max(SDbinaryData_TPShift_D1ctx);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
    TPShift_D1ctx_Cell_Totalactive = sum(SDbinaryData_TPShift_D1ctx, 1);%各細胞の活動量 (3SD以上のframe数)
    TPShift_D1ctx_Cell_Totalactive_Mean = mean(TPShift_D1ctx_Cell_Totalactive);%全細胞の平均活動量

    % Summarize 
    List_df_TPShift_D1ctx_mean(n,:) = df_TPShift_D1ctx_mean;%%TPShift6ではTPと変わらないようにshuffleしている
    List_TPShift_D1ctx_binarySD_SyncCells(n,:) = TPShift_D1ctx_binarySD_SyncCells;%
    List_TPShift_D1ctx_SyncCells_Ratio(n,:) = TPShift_D1ctx_SyncCells_Ratio;%
    List_TPShift_D1ctx_binarySD_Max(n,:)= TPShift_D1ctx_binarySD_Max;
    List_TPShift_D1ctx_Cell_Totalactive(n,:) = TPShift_D1ctx_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
      
     %Correlation matrix from Alan's(original name = out_matrix.m)
    %How to use: CorrMat_Alan(binned_datamatrix(time, neurons),shift_window(4))
%     shift_cormat = 4;%parameter setupに移動

%     CorrMat_Alan = CorrMat_Alan(dfData_TPShift6,shift_cormat);
%     writematrix(CorrMat_Alan, [mouseID,dinfo,'_TPshift6_',num2str(n),'_CorrMat.txt'], 'delimiter','\t');%
    celltype = char('TPshift_');
    d_content = char('CorMatSkip_');

    tic;
    TPshift_D1ctx_CorrMat_skip = CorrMat_Alan_As220506_skip(dfData_TPShift_D1ctx,window_cormat,shift_cormat);
    toc
    writematrix(TPshift_D1ctx_CorrMat_skip, [mouseID,dinfo,celltype,d_content,num2str(window_cormat),'_',num2str(shift_cormat),'.txt'], 'delimiter','\t');%
    
%     figure;%全shuffleをfigureにしているとデータ量が大変なので、最後のデータだけ一例として保存してみる(for...endの後にある)。
%     imagesc(CorrMat);
% 
%     dtime = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));%https://jp.mathworks.com/matlabcentral/answers/479642-
%     %一応もう一度時間をとる
%     saveas(gcf, [mouseID,'_CorrMat',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_CorrMat',dtime], 'tif');
    TPshift_D1ctx_CorrMat_sum = sum(TPshift_D1ctx_CorrMat_skip,'omitnan');%sum(CorrMat_Alan);
    List_TPShift_D1ctx_CorMatSkip_sum(n,:) = TPshift_D1ctx_CorrMat_sum;

end

% cal. Average(mean)/SD/SEM for graph
Ave_List_df_TPShift6_mean = mean(List_df_TPShift_D1ctx_mean);
SD_Lisr_df_TPShift6_mean = std(List_df_TPShift_D1ctx_mean);
SEM_List_df_TPShift6_mean = SD_Lisr_df_TPShift6_mean./sqrt(iter_TPShift);

Ave_List_TPShift6_binarySD_SyncCells = mean(List_TPShift_D1ctx_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
SD_List_TPShift6_binarySD_SyncCells = std(List_TPShift_D1ctx_binarySD_SyncCells);%全体平均は、イベント回数をそろえているので意味がない(差がない)
SEM_List_TPShift6_binarySD_SyncCells = SD_List_TPShift6_binarySD_SyncCells./sqrt(iter_TPShift);%全体平均は、イベント回数をそろえているので意味がない(差がない)

Ave_List_TPShift6_SyncCells_Ratio = mean(List_TPShift_D1ctx_SyncCells_Ratio);
SD_List_TPShift6_SyncCells_Ratio = std(List_TPShift_D1ctx_SyncCells_Ratio);
SEM_List_TPShift6_SyncCells_Ratio = SD_List_TPShift6_SyncCells_Ratio./sqrt(iter_TPShift);

Ave_List_TPShift6_binarySD_Max = mean(List_TPShift_D1ctx_binarySD_Max);
SD_List_TPShift6_binarySD_Max = std(List_TPShift_D1ctx_binarySD_Max);
SEM_List_TPShift6_binarySD_Max = SD_List_TPShift6_binarySD_Max./sqrt(iter_TPShift);

Ave_List_TPShift6_Cell_Totalactive = mean(List_TPShift_D1ctx_Cell_Totalactive);
SD_List_TPShift6_Cell_Totalactive = std(List_TPShift_D1ctx_Cell_Totalactive);
SEM_List_TPShift6_Cell_Totalactive = SD_List_TPShift6_Cell_Totalactive./sqrt(iter_TPShift);

% save the summary and average/SD/SEM
% check first session
writematrix(List_df_TPShift_D1ctx_mean,[mouseID,dinfo,celltype,num2str(n),'_summary_dF_mean.txt'], 'delimiter','\t');%
writematrix(List_TPShift_D1ctx_binarySD_SyncCells,[mouseID,dinfo,celltype,num2str(n),'_summary_binarySD_SyncCells.txt'], 'delimiter','\t');%
writematrix(List_TPShift_D1ctx_SyncCells_Ratio, [mouseID,dinfo,celltype,num2str(n),'_summary_SyncCell_Ratio.txt'], 'delimiter','\t');%
writematrix(List_TPShift_D1ctx_binarySD_Max, [mouseID,dinfo,celltype,num2str(n),'_summary_binarySD_Max.txt'], 'delimiter','\t');%
writematrix(List_TPShift_D1ctx_Cell_Totalactive, [mouseID,dinfo,celltype,num2str(n),'_summary_Cell_Totalactive.txt'], 'delimiter','\t');%
% % % for comparison between original TP vs shuffled TP
% % Summary_df_TPShift6_TP_mean = [df_TP_mean;List_df_TPShift6_mean];%

% Create figures
% parameters for figure
xc = linspace(1,NumCell_dF_TP,NumCell_dF_TP);
xt = linspace(1,Timeframe_dF_TP_D1ctx,Timeframe_dF_TP_D1ctx);

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

% figure;
% bar_val = [df_TP_mean;Ave_List_df_TPShift6_mean];
% bar(xc,bar_val);
% legend('raw','shuffle');
% title('dF average');

figure;
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
plot(xt,TP_D1ctx_binarySD_SyncCells)
hold on
plot(xt,Ave_List_TPShift6_binarySD_SyncCells);
legend
title('Synchronous active cells');

    for frame = [1, 1000, 2000, 3000, 4000, 5000, 6000] 
        figure;
        histogram(List_TPShift_D1ctx_binarySD_SyncCells(:,frame));
        hold on
        %histogram(TP_binarySD_SyncCells(1,:));
        bar(TP_D1ctx_binarySD_SyncCells(frame,:),2500,0.01);%
        legend('shuffle','raw')
        title('Synchronous active cells');
    end

figure;
% xt = linspace(1,Timeframe_dF_TP,Timeframe_dF_TP);
plot(xt,TP_D1ctx_SyncCells_Ratio,'LineWidth',1);
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
        plot(xt,TP_D1ctx_binarySD_SyncCells)
        hold on
        plot(xt,List_TPShift_D1ctx_binarySD_SyncCells(1:3,:),'LineWidth',0.75);
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


% %%
% 
%     
%     
%     [R_D1pre,P_D1pre] = corrcoef(dFData_D1pre);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
%     [R_D1ctx,P_D1ctx] = corrcoef(dFData_D1context);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
%     [R_D1after,P_D1after] = corrcoef(dFData_D1after);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
%     [R_D1post,P_D1post] = corrcoef(dFData_D1post);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
% 
%     % subtact self correlation (1)
% R_D1pre_sub = R_D1pre - eye(size(R_D1pre));
% R_D1ctx_sub = R_D1ctx - eye(size(R_D1ctx));
% R_D1after_sub = R_D1after - eye(size(R_D1after));
% R_D1post_sub = R_D1post - eye(size(R_D1post));
% 
% %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% R_D1pre_sub_0 = R_D1pre_sub;
% R_D1pre_sub_0(R_D1pre_sub_0~=0) = 1;
% 
% R_D1pre_sub_0_TP_sumsum = sum(sum(R_D1pre_sub_0(1:TP_last,1:TP_last)));
% R_D1pre_sub_0_TN_sumsum = sum(sum(R_D1pre_sub_0(TN_first:end,TN_first:end)));
% R_D1pre_sub_0_All_sumsum = sum(sum(R_D1pre_sub_0));
% 
% 
% R_D1ctx_sub_0 = R_D1ctx_sub;
% R_D1ctx_sub_0(R_D1ctx_sub_0~=0) = 1;
% 
% R_D1ctx_sub_0_TP_sumsum = sum(sum(R_D1ctx_sub_0(1:TP_last,1:TP_last)));
% R_D1ctx_sub_0_TN_sumsum = sum(sum(R_D1ctx_sub_0(TN_first:end,TN_first:end)));
% R_D1ctx_sub_0_All_sumsum = sum(sum(R_D1ctx_sub_0));
% 
% 
% R_D1after_sub_0 = R_D1after_sub;
% R_D1after_sub_0(R_D1after_sub_0~=0) = 1;
% 
% R_D1after_sub_0_TP_sumsum = sum(sum(R_D1after_sub_0(1:TP_last,1:TP_last)));
% R_D1after_sub_0_TN_sumsum = sum(sum(R_D1after_sub_0(TN_first:end,TN_first:end)));
% R_D1after_sub_0_All_sumsum = sum(sum(R_D1after_sub_0));
% 
% 
% R_D1post_sub_0 = R_D1post_sub;
% R_D1post_sub_0(R_D1post_sub_0~=0) = 1;
% 
% R_D1post_sub_0_TP_sumsum = sum(sum(R_D1post_sub_0(1:TP_last,1:TP_last)));
% R_D1post_sub_0_TN_sumsum = sum(sum(R_D1post_sub_0(TN_first:end,TN_first:end)));
% R_D1post_sub_0_All_sumsum = sum(sum(R_D1post_sub_0));
% 
% %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% R_D1pre_sub_TP_sumsum = sum(sum(R_D1pre_sub(1:TP_last,1:TP_last)));
% R_D1pre_sub_TN_sumsum = sum(sum(R_D1pre_sub(TN_first:end,TN_first:end)));
% R_D1pre_sub_All_sumsum = sum(sum(R_D1pre_sub));
% 
% R_D1ctx_sub_TP_sumsum = sum(sum(R_D1ctx_sub(1:TP_last,1:TP_last)));
% R_D1ctx_sub_TN_sumsum = sum(sum(R_D1ctx_sub(TN_first:end,TN_first:end)));
% R_D1ctx_sub_All_sumsum = sum(sum(R_D1ctx_sub));
% 
% R_D1after_sub_TP_sumsum = sum(sum(R_D1after_sub(1:TP_last,1:TP_last)));
% R_D1after_sub_TN_sumsum = sum(sum(R_D1after_sub(TN_first:end,TN_first:end)));
% R_D1after_sub_All_sumsum = sum(sum(R_D1after_sub));
% 
% R_D1post_sub_TP_sumsum = sum(sum(R_D1post_sub(1:TP_last,1:TP_last)));
% R_D1post_sub_TN_sumsum = sum(sum(R_D1post_sub(TN_first:end,TN_first:end)));
% R_D1post_sub_All_sumsum = sum(sum(R_D1post_sub));
% 
% %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% R_D1pre_sub_TP_average = R_D1pre_sub_TP_sumsum/R_D1pre_sub_0_TP_sumsum;
% R_D1pre_sub_TN_average = R_D1pre_sub_TN_sumsum/R_D1pre_sub_0_TN_sumsum;
% R_D1pre_sub_All_average = R_D1pre_sub_All_sumsum/R_D1pre_sub_0_All_sumsum;
% 
% R_D1ctx_sub_TP_average = R_D1ctx_sub_TP_sumsum/R_D1ctx_sub_0_TP_sumsum;
% R_D1ctx_sub_TN_average = R_D1ctx_sub_TN_sumsum/R_D1ctx_sub_0_TN_sumsum;
% R_D1ctx_sub_All_average = R_D1ctx_sub_All_sumsum/R_D1ctx_sub_0_All_sumsum;
% 
% R_D1after_sub_TP_average = R_D1after_sub_TP_sumsum/R_D1after_sub_0_TP_sumsum;
% R_D1after_sub_TN_average = R_D1after_sub_TN_sumsum/R_D1after_sub_0_TN_sumsum;
% R_D1after_sub_All_average = R_D1after_sub_All_sumsum/R_D1after_sub_0_All_sumsum;
% 
% R_D1post_sub_TP_average = R_D1post_sub_TP_sumsum/R_D1post_sub_0_TP_sumsum;
% R_D1post_sub_TN_average = R_D1post_sub_TN_sumsum/R_D1post_sub_0_TN_sumsum;
% R_D1post_sub_All_average = R_D1post_sub_All_sumsum/R_D1post_sub_0_All_sumsum;
% 
% %対角以外のますのMax/MinをセッションごとにTP/TN/Allで求める(add 221017)
% %Day1 pre
% R_D1pre_sub_TP_Max = max(max(R_D1pre_sub(1:TP_last,1:TP_last)));
% R_D1pre_sub_TN_Max = max(max(R_D1pre_sub(TN_first:end,TN_first:end)));
% R_D1pre_sub_All_Max = max(max(R_D1pre_sub));
% 
% R_D1pre_sub_TP_Min = min(min(R_D1pre_sub(1:TP_last,1:TP_last)));
% R_D1pre_sub_TN_Min = min(min(R_D1pre_sub(TN_first:end,TN_first:end)));
% R_D1pre_sub_All_Min = min(min(R_D1pre_sub));
% %Day1 context
% R_D1ctx_sub_TP_Max = max(max(R_D1ctx_sub(1:TP_last,1:TP_last)));
% R_D1ctx_sub_TN_Max = max(max(R_D1ctx_sub(TN_first:end,TN_first:end)));
% R_D1ctx_sub_All_Max = max(max(R_D1ctx_sub));
% 
% R_D1ctx_sub_TP_Min = min(min(R_D1ctx_sub(1:TP_last,1:TP_last)));
% R_D1ctx_sub_TN_Min = min(min(R_D1ctx_sub(TN_first:end,TN_first:end)));
% R_D1ctx_sub_All_Min = min(min(R_D1ctx_sub));
% %Day1 after
% R_D1after_sub_TP_Max = max(max(R_D1after_sub(1:TP_last,1:TP_last)));
% R_D1after_sub_TN_Max = max(max(R_D1after_sub(TN_first:end,TN_first:end)));
% R_D1after_sub_All_Max = max(max(R_D1after_sub));
% 
% R_D1after_sub_TP_Min = min(min(R_D1after_sub(1:TP_last,1:TP_last)));
% R_D1after_sub_TN_Min = min(min(R_D1after_sub(TN_first:end,TN_first:end)));
% R_D1after_sub_All_Min = min(min(R_D1after_sub));
% %Day1 post
% R_D1post_sub_TP_Max = max(max(R_D1post_sub(1:TP_last,1:TP_last)));
% R_D1post_sub_TN_Max = max(max(R_D1post_sub(TN_first:end,TN_first:end)));
% R_D1post_sub_All_Max = max(max(R_D1post_sub));
% 
% R_D1post_sub_TP_Min = min(min(R_D1post_sub(1:TP_last,1:TP_last)));
% R_D1post_sub_TN_Min = min(min(R_D1post_sub(TN_first:end,TN_first:end)));
% R_D1post_sub_All_Min = min(min(R_D1post_sub));
% 
% 
% %save the matrix
%     dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-
% 
%     celltype = char('All_');
%     d_content = char('dFcorr_');
% 
%     %Day1 pre
%     session = char('D1pre_');
% writematrix(R_D1pre_sub,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
%     %Day1 context
%     session = char('D1context_');
% writematrix(R_D1pre_sub,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
%     %Day1 after
%     session = char('D1after_');
% writematrix(R_D1pre_sub,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
%     %Day1 post
%     session = char('D1post_');
% writematrix(R_D1pre_sub,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
% 
% 
% % summary of results
% Cell = [{'Tom+'};{'Tom-'};{'All'}];
% R_D1pre_TPTNTotal_average = [R_D1pre_sub_TP_average;R_D1pre_sub_TN_average;R_D1pre_sub_All_average];
% R_D1ctx_TPTNTotal_average = [R_D1ctx_sub_TP_average;R_D1ctx_sub_TN_average;R_D1ctx_sub_All_average];
% R_D1after_TPTNTotal_average = [R_D1after_sub_TP_average;R_D1after_sub_TN_average;R_D1after_sub_All_average];
% R_D1post_TPTNTotal_average = [R_D1post_sub_TP_average;R_D1post_sub_TN_average;R_D1post_sub_All_average];
% 
% R_D1pre_TPTNTotal_Max = [R_D1pre_sub_TP_Max;R_D1pre_sub_TN_Max;R_D1pre_sub_All_Max];
% R_D1ctx_TPTNTotal_Max = [R_D1ctx_sub_TP_Max;R_D1ctx_sub_TN_Max;R_D1ctx_sub_All_Max];
% R_D1after_TPTNTotal_Max = [R_D1after_sub_TP_Max;R_D1after_sub_TN_Max;R_D1after_sub_All_Max];
% R_D1post_TPTNTotal_Max = [R_D1post_sub_TP_Max;R_D1post_sub_TN_Max;R_D1post_sub_All_Max];
% 
% R_D1pre_TPTNTotal_Min = [R_D1pre_sub_TP_Min;R_D1pre_sub_TN_Min;R_D1pre_sub_All_Min];
% R_D1ctx_TPTNTotal_Min = [R_D1ctx_sub_TP_Min;R_D1ctx_sub_TN_Min;R_D1ctx_sub_All_Min];
% R_D1after_TPTNTotal_Min = [R_D1after_sub_TP_Min;R_D1after_sub_TN_Min;R_D1after_sub_All_Min];
% R_D1post_TPTNTotal_Min = [R_D1post_sub_TP_Min;R_D1post_sub_TN_Min;R_D1post_sub_All_Min];
% 
% 
% % Day1 pre/context/after/post ver.
%     celltype = char('TpTnAll_');
%     d_content = char('dFcorrAve_');
% 
%     session = char('Day1_');
% Result_R_average = table(Cell, R_D1pre_TPTNTotal_average,R_D1ctx_TPTNTotal_average,R_D1after_TPTNTotal_average,R_D1post_TPTNTotal_average);
% Result_R_Max = table(Cell, R_D1pre_TPTNTotal_Max,R_D1ctx_TPTNTotal_Max,R_D1after_TPTNTotal_Max,R_D1post_TPTNTotal_Max);
% Result_R_Min = table(Cell, R_D1pre_TPTNTotal_Min,R_D1ctx_TPTNTotal_Min,R_D1after_TPTNTotal_Min,R_D1post_TPTNTotal_Min);
% %     % afterなしver. day1 pre/context/post ver
% % Result_R_average = table(Cell, R_D1pre_TPTNTotal_average,R_D1ctx_TPTNTotal_average,R_D1post_TPTNTotal_average);
% writetable(Result_R_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMax_');
% writetable(Result_R_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMin_');
% writetable(Result_R_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');%
% 
% % figure(average of correlation)
%     x_corr_All =  categorical({'Tom+','Tom-','All'});
%     x_corr_All = reordercats(x_corr_All,{'Tom+','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
% 
%     celltype = char('TpTnAll_');
%     d_content = char('dFcorrAve_');
%     session = char('D1pre_');
% 
% figure;bar(x_corr_All,R_D1pre_TPTNTotal_average);
% title("Day1 pre")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     session = char('D1context_');
% figure;bar(x_corr_All,R_D1ctx_TPTNTotal_average);
% title("Day1 context")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     session = char('D1after_');
% figure;bar(x_corr_All,R_D1after_TPTNTotal_average);
% title("Day1 after")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     session = char('D1post_');
% figure;bar(x_corr_All,R_D1post_TPTNTotal_average);
% title("Day1 post")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     % figure (correlation) > save
%     cmax = 0.10;cmin = -0.10;
%    
%     celltype = char('All_');
%     d_content = char('dFcorr_');
%     session = char('D1pre_');
% 
%     figure;imagesc(R_D1pre)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 pre");
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorr_sub_');
%     figure;imagesc(R_D1pre_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 pre")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
% %     cmax = 0.10;cmin = -0.10;
%     d_content = char('dFcorr_'); 
%     session = char('D1context_');
% 
%     figure;imagesc(R_D1ctx)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 context")
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorr_sub_'); 
%     figure;imagesc(R_D1ctx_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 context")
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
% %     cmax = 0.20;cmin = -0.20;
%     d_content = char('dFcorr_'); 
%     session = char('D1after_');
% 
%     figure;imagesc(R_D1after)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 after")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorr_sub_'); 
%     figure;imagesc(R_D1after_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 after");
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorr_'); 
%     session = char('D1post_');
%     figure;imagesc(R_D1post)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 post")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
%     
%     d_content = char('dFcorr_sub_'); 
%     figure;imagesc(R_D1post_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 post")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
% %     cd ..;%上のフォルダに移動
% 
%     %% 5. correlation (All, for TP_all vs TN all), Day2 (copy Day1 and paste new sheet, and replace D1 (Day1) to D2 (Day2)
% 
% %     mkdir('5_dFcorrelation') % フォルダの作成
% %     cd ('5_dFcorrelation');%今作成したフォルダに移動
%     
%     [R_D2pre,P_D2pre] = corrcoef(dFData_D2pre);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
%     [R_D2ctx,P_D2ctx] = corrcoef(dFData_D2context);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
%     [R_D2after,P_D2after] = corrcoef(dFData_D2after);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
%     [R_D2post,P_D2post] = corrcoef(dFData_D2post);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
% 
%     % subtact self correlation (1)
% R_D2pre_sub = R_D2pre - eye(size(R_D2pre));
% R_D2ctx_sub = R_D2ctx - eye(size(R_D2ctx));
% R_D2after_sub = R_D2after - eye(size(R_D2after));
% R_D2post_sub = R_D2post - eye(size(R_D2post));
% 
% %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% R_D2pre_sub_0 = R_D2pre_sub;
% R_D2pre_sub_0(R_D2pre_sub_0~=0) = 1;
% 
% R_D2pre_sub_0_TP_sumsum = sum(sum(R_D2pre_sub_0(1:TP_last,1:TP_last)));
% R_D2pre_sub_0_TN_sumsum = sum(sum(R_D2pre_sub_0(TN_first:end,TN_first:end)));
% R_D2pre_sub_0_All_sumsum = sum(sum(R_D2pre_sub_0));
% 
% 
% R_D2ctx_sub_0 = R_D2ctx_sub;
% R_D2ctx_sub_0(R_D2ctx_sub_0~=0) = 1;
% 
% R_D2ctx_sub_0_TP_sumsum = sum(sum(R_D2ctx_sub_0(1:TP_last,1:TP_last)));
% R_D2ctx_sub_0_TN_sumsum = sum(sum(R_D2ctx_sub_0(TN_first:end,TN_first:end)));
% R_D2ctx_sub_0_All_sumsum = sum(sum(R_D2ctx_sub_0));
% 
% 
% R_D2after_sub_0 = R_D2after_sub;
% R_D2after_sub_0(R_D2after_sub_0~=0) = 1;
% 
% R_D2after_sub_0_TP_sumsum = sum(sum(R_D2after_sub_0(1:TP_last,1:TP_last)));
% R_D2after_sub_0_TN_sumsum = sum(sum(R_D2after_sub_0(TN_first:end,TN_first:end)));
% R_D2after_sub_0_All_sumsum = sum(sum(R_D2after_sub_0));
% 
% 
% R_D2post_sub_0 = R_D2post_sub;
% R_D2post_sub_0(R_D2post_sub_0~=0) = 1;
% 
% R_D2post_sub_0_TP_sumsum = sum(sum(R_D2post_sub_0(1:TP_last,1:TP_last)));
% R_D2post_sub_0_TN_sumsum = sum(sum(R_D2post_sub_0(TN_first:end,TN_first:end)));
% R_D2post_sub_0_All_sumsum = sum(sum(R_D2post_sub_0));
% 
% %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% R_D2pre_sub_TP_sumsum = sum(sum(R_D2pre_sub(1:TP_last,1:TP_last)));
% R_D2pre_sub_TN_sumsum = sum(sum(R_D2pre_sub(TN_first:end,TN_first:end)));
% R_D2pre_sub_All_sumsum = sum(sum(R_D2pre_sub));
% 
% R_D2ctx_sub_TP_sumsum = sum(sum(R_D2ctx_sub(1:TP_last,1:TP_last)));
% R_D2ctx_sub_TN_sumsum = sum(sum(R_D2ctx_sub(TN_first:end,TN_first:end)));
% R_D2ctx_sub_All_sumsum = sum(sum(R_D2ctx_sub));
% 
% R_D2after_sub_TP_sumsum = sum(sum(R_D2after_sub(1:TP_last,1:TP_last)));
% R_D2after_sub_TN_sumsum = sum(sum(R_D2after_sub(TN_first:end,TN_first:end)));
% R_D2after_sub_All_sumsum = sum(sum(R_D2after_sub));
% 
% R_D2post_sub_TP_sumsum = sum(sum(R_D2post_sub(1:TP_last,1:TP_last)));
% R_D2post_sub_TN_sumsum = sum(sum(R_D2post_sub(TN_first:end,TN_first:end)));
% R_D2post_sub_All_sumsum = sum(sum(R_D2post_sub));
% 
% %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% R_D2pre_sub_TP_average = R_D2pre_sub_TP_sumsum/R_D2pre_sub_0_TP_sumsum;
% R_D2pre_sub_TN_average = R_D2pre_sub_TN_sumsum/R_D2pre_sub_0_TN_sumsum;
% R_D2pre_sub_All_average = R_D2pre_sub_All_sumsum/R_D2pre_sub_0_All_sumsum;
% 
% R_D2ctx_sub_TP_average = R_D2ctx_sub_TP_sumsum/R_D2ctx_sub_0_TP_sumsum;
% R_D2ctx_sub_TN_average = R_D2ctx_sub_TN_sumsum/R_D2ctx_sub_0_TN_sumsum;
% R_D2ctx_sub_All_average = R_D2ctx_sub_All_sumsum/R_D2ctx_sub_0_All_sumsum;
% 
% R_D2after_sub_TP_average = R_D2after_sub_TP_sumsum/R_D2after_sub_0_TP_sumsum;
% R_D2after_sub_TN_average = R_D2after_sub_TN_sumsum/R_D2after_sub_0_TN_sumsum;
% R_D2after_sub_All_average = R_D2after_sub_All_sumsum/R_D2after_sub_0_All_sumsum;
% 
% R_D2post_sub_TP_average = R_D2post_sub_TP_sumsum/R_D2post_sub_0_TP_sumsum;
% R_D2post_sub_TN_average = R_D2post_sub_TN_sumsum/R_D2post_sub_0_TN_sumsum;
% R_D2post_sub_All_average = R_D2post_sub_All_sumsum/R_D2post_sub_0_All_sumsum;
% 
% %対角以外のますのMax/MinをセッションごとにTP/TN/Allで求める(add 221017)
% %Day2 pre
% R_D2pre_sub_TP_Max = max(max(R_D2pre_sub(1:TP_last,1:TP_last)));
% R_D2pre_sub_TN_Max = max(max(R_D2pre_sub(TN_first:end,TN_first:end)));
% R_D2pre_sub_All_Max = max(max(R_D2pre_sub));
% 
% R_D2pre_sub_TP_Min = min(min(R_D2pre_sub(1:TP_last,1:TP_last)));
% R_D2pre_sub_TN_Min = min(min(R_D2pre_sub(TN_first:end,TN_first:end)));
% R_D2pre_sub_All_Min = min(min(R_D2pre_sub));
% %Day2 context
% R_D2ctx_sub_TP_Max = max(max(R_D2ctx_sub(1:TP_last,1:TP_last)));
% R_D2ctx_sub_TN_Max = max(max(R_D2ctx_sub(TN_first:end,TN_first:end)));
% R_D2ctx_sub_All_Max = max(max(R_D2ctx_sub));
% 
% R_D2ctx_sub_TP_Min = min(min(R_D2ctx_sub(1:TP_last,1:TP_last)));
% R_D2ctx_sub_TN_Min = min(min(R_D2ctx_sub(TN_first:end,TN_first:end)));
% R_D2ctx_sub_All_Min = min(min(R_D2ctx_sub));
% %Day2 after
% R_D2after_sub_TP_Max = max(max(R_D2after_sub(1:TP_last,1:TP_last)));
% R_D2after_sub_TN_Max = max(max(R_D2after_sub(TN_first:end,TN_first:end)));
% R_D2after_sub_All_Max = max(max(R_D2after_sub));
% 
% R_D2after_sub_TP_Min = min(min(R_D2after_sub(1:TP_last,1:TP_last)));
% R_D2after_sub_TN_Min = min(min(R_D2after_sub(TN_first:end,TN_first:end)));
% R_D2after_sub_All_Min = min(min(R_D2after_sub));
% %Day2 post
% R_D2post_sub_TP_Max = max(max(R_D2post_sub(1:TP_last,1:TP_last)));
% R_D2post_sub_TN_Max = max(max(R_D2post_sub(TN_first:end,TN_first:end)));
% R_D2post_sub_All_Max = max(max(R_D2post_sub));
% 
% R_D2post_sub_TP_Min = min(min(R_D2post_sub(1:TP_last,1:TP_last)));
% R_D2post_sub_TN_Min = min(min(R_D2post_sub(TN_first:end,TN_first:end)));
% R_D2post_sub_All_Min = min(min(R_D2post_sub));
% 
% %save the matrix
%     dtime = char(datetime('now', 'Format', 'yyMMddHHmm'));%https://jp.mathworks.com/matlabcentral/answers/479642-
% 
%     celltype = char('All_');
%     d_content = char('dFcorr_');
% 
%     %Day2 pre
%     session = char('D2pre_');
% writematrix(R_D2pre_sub,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
%     %Day2 context
%     session = char('D2context_');
% writematrix(R_D2pre_sub,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
%     %Day2 after
%     session = char('D2after_');
% writematrix(R_D2pre_sub,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
%     %Day2 post
%     session = char('D2post_');
% writematrix(R_D2pre_sub,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'],'delimiter','\t');%
% 
% % summary of results
% Cell = [{'Tom+'};{'Tom-'};{'All'}];
% R_D2pre_TPTNTotal_average = [R_D2pre_sub_TP_average;R_D2pre_sub_TN_average;R_D2pre_sub_All_average];
% R_D2ctx_TPTNTotal_average = [R_D2ctx_sub_TP_average;R_D2ctx_sub_TN_average;R_D2ctx_sub_All_average];
% R_D2after_TPTNTotal_average = [R_D2after_sub_TP_average;R_D2after_sub_TN_average;R_D2after_sub_All_average];
% R_D2post_TPTNTotal_average = [R_D2post_sub_TP_average;R_D2post_sub_TN_average;R_D2post_sub_All_average];
% 
% R_D2pre_TPTNTotal_Max = [R_D2pre_sub_TP_Max;R_D2pre_sub_TN_Max;R_D2pre_sub_All_Max];
% R_D2ctx_TPTNTotal_Max = [R_D2ctx_sub_TP_Max;R_D2ctx_sub_TN_Max;R_D2ctx_sub_All_Max];
% R_D2after_TPTNTotal_Max = [R_D2after_sub_TP_Max;R_D2after_sub_TN_Max;R_D2after_sub_All_Max];
% R_D2post_TPTNTotal_Max = [R_D2post_sub_TP_Max;R_D2post_sub_TN_Max;R_D2post_sub_All_Max];
% 
% R_D2pre_TPTNTotal_Min = [R_D2pre_sub_TP_Min;R_D2pre_sub_TN_Min;R_D2pre_sub_All_Min];
% R_D2ctx_TPTNTotal_Min = [R_D2ctx_sub_TP_Min;R_D2ctx_sub_TN_Min;R_D2ctx_sub_All_Min];
% R_D2after_TPTNTotal_Min = [R_D2after_sub_TP_Min;R_D2after_sub_TN_Min;R_D2after_sub_All_Min];
% R_D2post_TPTNTotal_Min = [R_D2post_sub_TP_Min;R_D2post_sub_TN_Min;R_D2post_sub_All_Min];
% 
% 
% % Day2 pre/context/after/post ver.
%     celltype = char('TpTnAll_');
%     d_content = char('dFcorrAve_');
% 
%     session = char('Day2_');
% Result_R_average = table(Cell, R_D2pre_TPTNTotal_average,R_D2ctx_TPTNTotal_average,R_D2after_TPTNTotal_average,R_D2post_TPTNTotal_average);
% Result_R_Max = table(Cell, R_D2pre_TPTNTotal_Max,R_D2ctx_TPTNTotal_Max,R_D2after_TPTNTotal_Max,R_D2post_TPTNTotal_Max);
% Result_R_Min = table(Cell, R_D2pre_TPTNTotal_Min,R_D2ctx_TPTNTotal_Min,R_D2after_TPTNTotal_Min,R_D2post_TPTNTotal_Min);
% %     % afterなしver. Day2 pre/context/post ver
% % Result_R_average = table(Cell, R_D2pre_TPTNTotal_average,R_D2ctx_TPTNTotal_average,R_D2post_TPTNTotal_average);
% writetable(Result_R_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMax_');
% writetable(Result_R_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMin_');
% writetable(Result_R_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');%
% 
% % figure(average of correlation)
%     x_corr_All =  categorical({'Tom+','Tom-','All'});
%     x_corr_All = reordercats(x_corr_All,{'Tom+','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
% 
%     celltype = char('TpTnAll_');
%     d_content = char('dFcorrAve_');
%     session = char('D2pre_');
% 
% figure;bar(x_corr_All,R_D2pre_TPTNTotal_average);
% title("Day2 pre")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     session = char('D2context_');
% figure;bar(x_corr_All,R_D2ctx_TPTNTotal_average);
% title("Day2 context")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     session = char('D2after_');
% figure;bar(x_corr_All,R_D2after_TPTNTotal_average);
% title("Day2 after")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     session = char('D2post_');
% figure;bar(x_corr_All,R_D2post_TPTNTotal_average);
% title("Day2 post")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     % figure (correlation) > save
%     cmax = 0.10;cmin = -0.10;
%    
%     celltype = char('All_');
%     d_content = char('dFcorr_');
%     session = char('D2pre_');
% 
%     figure;imagesc(R_D2pre)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 pre");
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorr_sub_');
%     figure;imagesc(R_D2pre_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 pre")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
% %     cmax = 0.10;cmin = -0.10;
%     d_content = char('dFcorr_'); 
%     session = char('D2context_');
% 
%     figure;imagesc(R_D2ctx)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 context")
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorr_sub_'); 
%     figure;imagesc(R_D2ctx_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 context")
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
% %     cmax = 0.20;cmin = -0.20;
%     d_content = char('dFcorr_'); 
%     session = char('D2after_');
% 
%     figure;imagesc(R_D2after)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 after")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorr_sub_'); 
%     figure;imagesc(R_D2after_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 after");
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorr_'); 
%     session = char('D2post_');
%     figure;imagesc(R_D2post)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 post")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
%     
%     d_content = char('dFcorr_sub_'); 
%     figure;imagesc(R_D2post_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 post")  
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,dinfo,session,celltype,d_content,dtime], 'tif');
% 
% %     cd ..;%上のフォルダに移動
% 
% 
% %% Correlation (TNrand, for TP vs TNrand, TN, All)
% %     mkdir('5_dFcorrelation') % フォルダの作成
% %     cd ('5_dFcorrelation');%今作成したフォルダに移動
%     
% %create rondom data%%%%%%%%%%%%%%%%%%%%%%%%%
% % 5-1-0. save the shuffled data
% ftime2 = char(datetime('now','Format','yyyyMMdd'));%フォルダ名, (データ名)
% celltype = char('TNrand_');
% mkdir([mouseID,dinfo,'_TNrand_',ftime2]); % フォルダの作成
% cd ([mouseID,dinfo,'_TNrand_',ftime2]);
% 
% % 5-1-1. Setting parameters
% % for creating randomized TN data
% dFData_TN = Data_raw(:, TN_first:end);%extract Tom- cells/remove Tom+ cells
% NumCell_TN = size(dFData_TN,2);%the number of Tom- cells
% Timeframe_dF_TN = size(dFData_TN,1);%Time(frame)数
% 
%     NumCell_dF_TP = TP_Cell;% the number of cells selected from Tom- cells
% dFData_TNrand = zeros(Timeframe_dF_TN,NumCell_dF_TP);
%     NumCell_TNrand = size(dFData_TNrand,2);%the number of cells,(TNrandはTPに数をそろえているので、TPの数で代用してもいいが、念のため)
%     Timeframe_dF_TNrand = size(dFData_TNrand,1);%Time(frame)数
% 
% % the number of shuffle/random sampling
% iter_TNrand = 10000;%5000;%1000:p=0.01を検出するためには最低でもこれくらい必要なのでは？ Num_Time_frames_TP/10;Num_Cell_TP;%30;%Shuffleの回数(shuffleデータの数、この平均をshuffleデータとする予定。)
% 
% % paremeters for cal.
% SD_thres = 3;%binaryの時のSDのthreshold, このsessionだけで完結させるように一応入れておく
% 
% % 5-1-2. Matrix for summarizing repeated data
% 
% %%%%% Day1 pre %%%%%
% List_df_TNrand_D1pre_mean = zeros(iter_TNrand,NumCell_TNrand);%template for all sessions
% 
% Timeframe_dF_TNrand = Day1pre_frame;
% List_TNrand_D1pre_binarySD_SyncCells = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D1pre_SyncCells_Ratio = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D1pre_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
% List_TNrand_D1pre_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);
% 
% List_R_D1pre_TNrand_sub_average = zeros(iter_TNrand, 1);
% List_R_D1pre_TNrand_sub_Max = zeros(iter_TNrand, 1);
% List_R_D1pre_TNrand_sub_Min = zeros(iter_TNrand, 1);
% 
% %%%%% Day1 context %%%%%
% List_df_TNrand_D1ctx_mean = zeros(iter_TNrand,NumCell_TNrand);%template for all sessions
% 
% Timeframe_dF_TNrand = Day1context_frame;
% List_TNrand_D1ctx_binarySD_SyncCells = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D1ctx_SyncCells_Ratio = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D1ctx_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
% List_TNrand_D1ctx_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);
% 
% List_R_D1ctx_TNrand_sub_average = zeros(iter_TNrand, 1);
% List_R_D1ctx_TNrand_sub_Max = zeros(iter_TNrand, 1);
% List_R_D1ctx_TNrand_sub_Min = zeros(iter_TNrand, 1);
% 
% %%%%% Day1 after %%%%%
% List_df_TNrand_D1after_mean = zeros(iter_TNrand,NumCell_TNrand);%template for all sessions
% 
% Timeframe_dF_TNrand = Day1after_frame;
% List_TNrand_D1after_binarySD_SyncCells = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D1after_SyncCells_Ratio = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D1after_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
% List_TNrand_D1after_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);
% 
% List_R_D1after_TNrand_sub_average = zeros(iter_TNrand, 1);
% List_R_D1after_TNrand_sub_Max = zeros(iter_TNrand, 1);
% List_R_D1after_TNrand_sub_Min = zeros(iter_TNrand, 1);
% 
% %%%%% Day1 post %%%%%
% List_df_TNrand_D1post_mean = zeros(iter_TNrand,NumCell_TNrand);%template for all sessions
% 
% Timeframe_dF_TNrand = Day1post_frame;
% List_TNrand_D1post_binarySD_SyncCells = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D1post_SyncCells_Ratio = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D1post_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
% List_TNrand_D1post_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);
% 
% List_R_D1post_TNrand_sub_average = zeros(iter_TNrand, 1);
% List_R_D1post_TNrand_sub_Max = zeros(iter_TNrand, 1);
% List_R_D1post_TNrand_sub_Min = zeros(iter_TNrand, 1);
% 
% %%%%% Day2 pre %%%%%
% List_df_TNrand_D2pre_mean = zeros(iter_TNrand,NumCell_TNrand);%template for all sessions
% 
% Timeframe_dF_TNrand = Day2pre_frame;
% List_TNrand_D2pre_binarySD_SyncCells = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D2pre_SyncCells_Ratio = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D2pre_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
% List_TNrand_D2pre_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);
% 
% List_R_D2pre_TNrand_sub_average = zeros(iter_TNrand, 1);
% List_R_D2pre_TNrand_sub_Max = zeros(iter_TNrand, 1);
% List_R_D2pre_TNrand_sub_Min = zeros(iter_TNrand, 1);
% 
% %%%%% Day2 context %%%%%
% List_df_TNrand_D2ctx_mean = zeros(iter_TNrand,NumCell_TNrand);%template for all sessions
% 
% Timeframe_dF_TNrand = Day2context_frame;
% List_TNrand_D2ctx_binarySD_SyncCells = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D2ctx_SyncCells_Ratio = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D2ctx_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
% List_TNrand_D2ctx_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);
% 
% List_R_D2ctx_TNrand_sub_average = zeros(iter_TNrand, 1);
% List_R_D2ctx_TNrand_sub_Max = zeros(iter_TNrand, 1);
% List_R_D2ctx_TNrand_sub_Min = zeros(iter_TNrand, 1);
% 
% %%%%% Day2 after %%%%%
% List_df_TNrand_D2after_mean = zeros(iter_TNrand,NumCell_TNrand);%template for all sessions
% 
% Timeframe_dF_TNrand = Day2after_frame;
% List_TNrand_D2after_binarySD_SyncCells = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D2after_SyncCells_Ratio = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D2after_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
% List_TNrand_D2after_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);
% 
% List_R_D2after_TNrand_sub_average = zeros(iter_TNrand, 1);
% List_R_D2after_TNrand_sub_Max = zeros(iter_TNrand, 1);
% List_R_D2after_TNrand_sub_Min = zeros(iter_TNrand, 1);
% 
% %%%%% Day2 post %%%%%
% List_df_TNrand_D2post_mean = zeros(iter_TNrand,NumCell_TNrand);%template for all sessions
% 
% Timeframe_dF_TNrand = Day2post_frame;
% List_TNrand_D2post_binarySD_SyncCells = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D2post_SyncCells_Ratio = zeros(Timeframe_dF_TNrand,iter_TNrand);%for..endの開始前に配置
% List_TNrand_D2post_binarySD_Max = zeros(iter_TNrand,NumCell_TNrand);
% List_TNrand_D2post_Cell_Totalactive = zeros(iter_TNrand,NumCell_TNrand);
% 
% List_R_D2post_TNrand_sub_average = zeros(iter_TNrand, 1);
% List_R_D2post_TNrand_sub_Max = zeros(iter_TNrand, 1);
% List_R_D2post_TNrand_sub_Min = zeros(iter_TNrand, 1);
% 
% % 5-1-3. Create data sets of rondomized Tom Nega > cal. mean etc. >
% % Summarize the results
% tic;
% for i = 1:iter_TNrand
%     % Create random IDs for Shuffle
%     TNrand_IDs = randperm(NumCell_TN,NumCell_TN);
%     
%     % Create Shuffled TN data
%     TNrandIDs_DataTN = [TNrand_IDs; dFData_TN];%concatenate shuffled ID with raw data of Tom Nega
%     TNrandIDs_DataTN_sort = sortrows(TNrandIDs_DataTN')';%sort, sortrawsで列ごとにsortするために、転置(')してshuffle IDを(1列目から)1行目に持ってきて、sort(rows)したのちに、再転置('), https://jp.mathworks.com/help/matlab/ref/double.sortrows.html
%     Data_TNShuf = TNrandIDs_DataTN_sort(2:end,:);%remove shuffled IDs
%     dFData_TNrand = Data_TNShuf(:,1:NumCell_dF_TP);
%     
% % devide session
% dFData_TNrand_D1pre = dFData_TNrand(Day1pre_first:Day1pre_last,:);%remove time and Tom+/-
% % dFData_TNrand_D1context = dFData_TNrand(Day1context_first:Day1context_last,:);%remove time and Tom+/-
% dFData_TNrand_D1ctx = dFData_TNrand(Day1context_first:Day1context_last,:);%remove time and Tom+/-
% dFData_TNrand_D1after = dFData_TNrand(Day1after_first:Day1after_last,:);%remove time and Tom+/-
% dFData_TNrand_D1post = dFData_TNrand(Day1post_first:Day1post_last,:);%remove time and Tom+/-
% 
% dFData_TNrand_D2pre = dFData_TNrand(Day2pre_first:Day2pre_last,:);%remove time and Tom+/-
% dFData_TNrand_D2ctx = dFData_TNrand(Day2context_first:Day2context_last,:);%remove time and Tom+/-
% dFData_TNrand_D2after = dFData_TNrand(Day2after_first:Day2after_last,:);%remove time and Tom+/-
% dFData_TNrand_D2post = dFData_TNrand(Day2post_first:Day2post_last,:);%remove time and Tom+/-
% 
% %%%%%cal.%%%%%%%%%%%
% %%%%% Day1 pre %%%%%
% % threshold
%     df_TNrand_D1pre_mean = mean(dFData_TNrand_D1pre);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
%     df_TNrand_D1pre_sd= std(dFData_TNrand_D1pre);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% 
%     sdData_TNrand_D1pre = (dFData_TNrand_D1pre - df_TNrand_D1pre_mean)./df_TNrand_D1pre_sd;%各データの細胞ごとのSD
%     SDbinaryData_TNrand_D1pre = sdData_TNrand_D1pre>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% 
%     % results from binarized df matrix
%     TNrand_D1pre_binarySD_SyncCells= sum(SDbinaryData_TNrand_D1pre, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%     %細胞数、frame数は既に取得済み
%     TNrand_D1pre_SyncCells_Ratio = TNrand_D1pre_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
% 
%     TNrand_D1pre_binarySD_Max = max(SDbinaryData_TNrand_D1pre);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
%     TNrand_D1pre_Cell_Totalactive = sum(SDbinaryData_TNrand_D1pre, 1);%各細胞の活動量 (2SD以上のframe数)
%     TNrand_D1pre_Cell_Totalactive_Mean = mean(TNrand_D1pre_Cell_Totalactive);%全細胞の平均活動量
% 
%     % Summarize results from randomized data
%     List_df_TNrand_D1pre_mean(i,:) = df_TNrand_D1pre_mean;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
%     List_TNrand_D1pre_binarySD_SyncCells(:,i) = TNrand_D1pre_binarySD_SyncCells;%
%     List_TNrand_D1pre_SyncCells_Ratio(:,i) = TNrand_D1pre_SyncCells_Ratio;%
%     List_TNrand_D1pre_binarySD_Max(i,:)= TNrand_D1pre_binarySD_Max;
%     List_TNrand_D1pre_Cell_Totalactive(i,:) = TNrand_D1pre_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
% %%%%% Day1 context %%%%%
% % threshold
%     df_TNrand_D1ctx_mean = mean(dFData_TNrand_D1ctx);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
%     df_TNrand_D1ctx_sd= std(dFData_TNrand_D1ctx);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% 
%     sdData_TNrand_D1ctx = (dFData_TNrand_D1ctx - df_TNrand_D1ctx_mean)./df_TNrand_D1ctx_sd;%各データの細胞ごとのSD
%     SDbinaryData_TNrand_D1ctx = sdData_TNrand_D1ctx>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% 
%     % results from binarized df matrix
%     TNrand_D1ctx_binarySD_SyncCells= sum(SDbinaryData_TNrand_D1ctx, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%     %細胞数、frame数は既に取得済み
%     TNrand_D1ctx_SyncCells_Ratio = TNrand_D1ctx_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
% 
%     TNrand_D1ctx_binarySD_Max = max(SDbinaryData_TNrand_D1ctx);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
%     TNrand_D1ctx_Cell_Totalactive = sum(SDbinaryData_TNrand_D1ctx, 1);%各細胞の活動量 (2SD以上のframe数)
%     TNrand_D1ctx_Cell_Totalactive_Mean = mean(TNrand_D1ctx_Cell_Totalactive);%全細胞の平均活動量
% 
%     % Summarize results from randomized data
%     List_df_TNrand_D1ctx_mean(i,:) = df_TNrand_D1ctx_mean;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
%     List_TNrand_D1ctx_binarySD_SyncCells(:,i) = TNrand_D1ctx_binarySD_SyncCells;%
%     List_TNrand_D1ctx_SyncCells_Ratio(:,i) = TNrand_D1ctx_SyncCells_Ratio;%
%     List_TNrand_D1ctx_binarySD_Max(i,:)= TNrand_D1ctx_binarySD_Max;
%     List_TNrand_D1ctx_Cell_Totalactive(i,:) = TNrand_D1ctx_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
% % %%%%% Day1 after %%%%%
% % % threshold
% %     df_TNrand_D1after_mean = mean(dFData_TNrand_D1after);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
% %     df_TNrand_D1after_sd= std(dFData_TNrand_D1after);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% % 
% %     sdData_TNrand_D1after = (dFData_TNrand_D1after - df_TNrand_D1after_mean)./df_TNrand_D1after_sd;%各データの細胞ごとのSD
% %     SDbinaryData_TNrand_D1after = sdData_TNrand_D1after>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% % 
% %     % results from binarized df matrix
% %     TNrand_D1after_binarySD_SyncCells= sum(SDbinaryData_TNrand_D1after, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
% %     %細胞数、frame数は既に取得済み
% %     TNrand_D1after_SyncCells_Ratio = TNrand_D1after_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
% % 
% %     TNrand_D1after_binarySD_Max = max(SDbinaryData_TNrand_D1after);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
% %     TNrand_D1after_Cell_Totalactive = sum(SDbinaryData_TNrand_D1after, 1);%各細胞の活動量 (2SD以上のframe数)
% %     TNrand_D1after_Cell_Totalactive_Mean = mean(TNrand_D1after_Cell_Totalactive);%全細胞の平均活動量
% % 
% %     % Summarize results from randomized data
% %     List_df_TNrand_D1after_mean(i,:) = df_TNrand_D1after_mean;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
% %     List_TNrand_D1after_binarySD_SyncCells(:,i) = TNrand_D1after_binarySD_SyncCells;%
% %     List_TNrand_D1after_SyncCells_Ratio(:,i) = TNrand_D1after_SyncCells_Ratio;%
% %     List_TNrand_D1after_binarySD_Max(i,:)= TNrand_D1after_binarySD_Max;
% %     List_TNrand_D1after_Cell_Totalactive(i,:) = TNrand_D1after_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
% %%%%% Day1 post %%%%%
% % threshold
%     df_TNrand_D1post_mean = mean(dFData_TNrand_D1post);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
%     df_TNrand_D1post_sd= std(dFData_TNrand_D1post);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% 
%     sdData_TNrand_D1post = (dFData_TNrand_D1post - df_TNrand_D1post_mean)./df_TNrand_D1post_sd;%各データの細胞ごとのSD
%     SDbinaryData_TNrand_D1post = sdData_TNrand_D1post>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% 
%     % results from binarized df matrix
%     TNrand_D1post_binarySD_SyncCells= sum(SDbinaryData_TNrand_D1post, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%     %細胞数、frame数は既に取得済み
%     TNrand_D1post_SyncCells_Ratio = TNrand_D1post_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
% 
%     TNrand_D1post_binarySD_Max = max(SDbinaryData_TNrand_D1post);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
%     TNrand_D1post_Cell_Totalactive = sum(SDbinaryData_TNrand_D1post, 1);%各細胞の活動量 (2SD以上のframe数)
%     TNrand_D1post_Cell_Totalactive_Mean = mean(TNrand_D1post_Cell_Totalactive);%全細胞の平均活動量
% 
%     % Summarize results from randomized data
%     List_df_TNrand_D1post_mean(i,:) = df_TNrand_D1post_mean;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
%     List_TNrand_D1post_binarySD_SyncCells(:,i) = TNrand_D1post_binarySD_SyncCells;%
%     List_TNrand_D1post_SyncCells_Ratio(:,i) = TNrand_D1post_SyncCells_Ratio;%
%     List_TNrand_D1post_binarySD_Max(i,:)= TNrand_D1post_binarySD_Max;
%     List_TNrand_D1post_Cell_Totalactive(i,:) = TNrand_D1post_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
% %%%%% Day2 pre %%%%%
% % threshold
%     df_TNrand_D2pre_mean = mean(dFData_TNrand_D2pre);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
%     df_TNrand_D2pre_sd= std(dFData_TNrand_D2pre);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% 
%     sdData_TNrand_D2pre = (dFData_TNrand_D2pre - df_TNrand_D2pre_mean)./df_TNrand_D2pre_sd;%各データの細胞ごとのSD
%     SDbinaryData_TNrand_D2pre = sdData_TNrand_D2pre>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% 
%     % results from binarized df matrix
%     TNrand_D2pre_binarySD_SyncCells= sum(SDbinaryData_TNrand_D2pre, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%     %細胞数、frame数は既に取得済み
%     TNrand_D2pre_SyncCells_Ratio = TNrand_D2pre_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
% 
%     TNrand_D2pre_binarySD_Max = max(SDbinaryData_TNrand_D2pre);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
%     TNrand_D2pre_Cell_Totalactive = sum(SDbinaryData_TNrand_D2pre, 1);%各細胞の活動量 (2SD以上のframe数)
%     TNrand_D2pre_Cell_Totalactive_Mean = mean(TNrand_D2pre_Cell_Totalactive);%全細胞の平均活動量
% 
%     % Summarize results from randomized data
%     List_df_TNrand_D2pre_mean(i,:) = df_TNrand_D2pre_mean;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
%     List_TNrand_D2pre_binarySD_SyncCells(:,i) = TNrand_D2pre_binarySD_SyncCells;%
%     List_TNrand_D2pre_SyncCells_Ratio(:,i) = TNrand_D2pre_SyncCells_Ratio;%
%     List_TNrand_D2pre_binarySD_Max(i,:)= TNrand_D2pre_binarySD_Max;
%     List_TNrand_D2pre_Cell_Totalactive(i,:) = TNrand_D2pre_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
% %%%%% Day2 context %%%%%
% % threshold
%     df_TNrand_D2ctx_mean = mean(dFData_TNrand_D2ctx);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
%     df_TNrand_D2ctx_sd= std(dFData_TNrand_D2ctx);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% 
%     sdData_TNrand_D2ctx = (dFData_TNrand_D2ctx - df_TNrand_D2ctx_mean)./df_TNrand_D2ctx_sd;%各データの細胞ごとのSD
%     SDbinaryData_TNrand_D2ctx = sdData_TNrand_D2ctx>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% 
%     % results from binarized df matrix
%     TNrand_D2ctx_binarySD_SyncCells= sum(SDbinaryData_TNrand_D2ctx, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%     %細胞数、frame数は既に取得済み
%     TNrand_D2ctx_SyncCells_Ratio = TNrand_D2ctx_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
% 
%     TNrand_D2ctx_binarySD_Max = max(SDbinaryData_TNrand_D2ctx);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
%     TNrand_D2ctx_Cell_Totalactive = sum(SDbinaryData_TNrand_D2ctx, 1);%各細胞の活動量 (2SD以上のframe数)
%     TNrand_D2ctx_Cell_Totalactive_Mean = mean(TNrand_D2ctx_Cell_Totalactive);%全細胞の平均活動量
% 
%     % Summarize results from randomized data
%     List_df_TNrand_D2ctx_mean(i,:) = df_TNrand_D2ctx_mean;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
%     List_TNrand_D2ctx_binarySD_SyncCells(:,i) = TNrand_D2ctx_binarySD_SyncCells;%
%     List_TNrand_D2ctx_SyncCells_Ratio(:,i) = TNrand_D2ctx_SyncCells_Ratio;%
%     List_TNrand_D2ctx_binarySD_Max(i,:)= TNrand_D2ctx_binarySD_Max;
%     List_TNrand_D2ctx_Cell_Totalactive(i,:) = TNrand_D2ctx_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
% % %%%%% Day2 after %%%%%
% % % threshold
% %     df_TNrand_D2after_mean = mean(dFData_TNrand_D2after);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
% %     df_TNrand_D2after_sd= std(dFData_TNrand_D2after);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% % 
% %     sdData_TNrand_D2after = (dFData_TNrand_D2after - df_TNrand_D2after_mean)./df_TNrand_D2after_sd;%各データの細胞ごとのSD
% %     SDbinaryData_TNrand_D2after = sdData_TNrand_D2after>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% % 
% %     % results from binarized df matrix
% %     TNrand_D2after_binarySD_SyncCells= sum(SDbinaryData_TNrand_D2after, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
% %     %細胞数、frame数は既に取得済み
% %     TNrand_D2after_SyncCells_Ratio = TNrand_D2after_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
% % 
% %     TNrand_D2after_binarySD_Max = max(SDbinaryData_TNrand_D2after);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
% %     TNrand_D2after_Cell_Totalactive = sum(SDbinaryData_TNrand_D2after, 1);%各細胞の活動量 (2SD以上のframe数)
% %     TNrand_D2after_Cell_Totalactive_Mean = mean(TNrand_D2after_Cell_Totalactive);%全細胞の平均活動量
% % 
% %     % Summarize results from randomized data
% %     List_df_TNrand_D2after_mean(i,:) = df_TNrand_D2after_mean;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
% %     List_TNrand_D2after_binarySD_SyncCells(:,i) = TNrand_D2after_binarySD_SyncCells;%
% %     List_TNrand_D2after_SyncCells_Ratio(:,i) = TNrand_D2after_SyncCells_Ratio;%
% %     List_TNrand_D2after_binarySD_Max(i,:)= TNrand_D2after_binarySD_Max;
% %     List_TNrand_D2after_Cell_Totalactive(i,:) = TNrand_D2after_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
% %%%%% Day2 post %%%%%
% % threshold
%     df_TNrand_D2post_mean = mean(dFData_TNrand_D2post);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
%     df_TNrand_D2post_sd= std(dFData_TNrand_D2post);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
% 
%     sdData_TNrand_D2post = (dFData_TNrand_D2post - df_TNrand_D2post_mean)./df_TNrand_D2post_sd;%各データの細胞ごとのSD
%     SDbinaryData_TNrand_D2post = sdData_TNrand_D2post>SD_thres;%2SD以上を1, それ未満を0に。https://jp.mathworks.com/help/matlab/matlab_prog/find-array-elements-that-meet-a-condition.html?searchHighlight=%E6%9D%A1%E4%BB%B6&s_tid=srchtitle_%25E6%259D%25A1%25E4%25BB%25B6_1
% 
%     % results from binarized df matrix
%     TNrand_D2post_binarySD_SyncCells= sum(SDbinaryData_TNrand_D2post, 2);%簡易同期チェック, https://jp.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum&s_tid=srchtitle_sum_1
%     %細胞数、frame数は既に取得済み
%     TNrand_D2post_SyncCells_Ratio = TNrand_D2post_binarySD_SyncCells./NumCell_TNrand;%同期活動している細胞の割合=同期活動している細胞数/全細胞数
% 
%     TNrand_D2post_binarySD_Max = max(SDbinaryData_TNrand_D2post);%活動していない細胞の有無確認 (shuffle dataで活動していない細胞ができていないか確認のため, just in case), https://jp.mathworks.com/help/matlab/ref/max.html?searchHighlight=%E6%9C%80%E5%A4%A7%E5%80%A4%20max&s_tid=srchtitle_%25E6%259C%2580%25E5%25A4%25A7%25E5%2580%25A4%20max_1
%     TNrand_D2post_Cell_Totalactive = sum(SDbinaryData_TNrand_D2post, 1);%各細胞の活動量 (2SD以上のframe数)
%     TNrand_D2post_Cell_Totalactive_Mean = mean(TNrand_D2post_Cell_Totalactive);%全細胞の平均活動量
% 
%     % Summarize results from randomized data
%     List_df_TNrand_D2post_mean(i,:) = df_TNrand_D2post_mean;%,,,,%TPShift6ではTPと変わらないようにshuffleしている
%     List_TNrand_D2post_binarySD_SyncCells(:,i) = TNrand_D2post_binarySD_SyncCells;%
%     List_TNrand_D2post_SyncCells_Ratio(:,i) = TNrand_D2post_SyncCells_Ratio;%
%     List_TNrand_D2post_binarySD_Max(i,:)= TNrand_D2post_binarySD_Max;
%     List_TNrand_D2post_Cell_Totalactive(i,:) = TNrand_D2post_Cell_Totalactive;%%TPShift6ではTPと変わらないようにshuffleしている
% 
% 
% %%%%% Correlation %%%%%%%%%%%%%
% %%%%% Day1 pre %%%%%
% [R_D1pre_TNrand,P_D1pre_TNrand] = corrcoef(dFData_TNrand_D1pre);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
% % subtact self correlation (1)
% R_D1pre_TNrand_sub = R_D1pre_TNrand - eye(size(R_D1pre_TNrand));
% %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% R_D1pre_TNrand_sub_0 = R_D1pre_TNrand_sub;
% R_D1pre_TNrand_sub_0(R_D1pre_TNrand_sub_0~=0) = 1;
% % R_D1pre_TNrand_sub_0_TN_sumsum = sum(sum(R_D1pre_TNrand_sub_0(TN_first:end,TN_first:end)));
% R_D1pre_TNrand_sub_0_All_sumsum = sum(sum(R_D1pre_TNrand_sub_0));
% %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% % R_D1pre_TNrand_sub_TN_sumsum = sum(sum(R_D1pre_TNrand_sub(TN_first:end,TN_first:end)));
% R_D1pre_TNrand_sub_All_sumsum = sum(sum(R_D1pre_TNrand_sub));
% %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% % R_D1pre_TNrand_sub_TN_average = R_D1pre_TNrand_sub_TN_sumsum/R_D1pre_TNrand_sub_0_TN_sumsum;
% R_D1pre_TNrand_sub_All_average = R_D1pre_TNrand_sub_All_sumsum/R_D1pre_TNrand_sub_0_All_sumsum;
% %対角以外のますの数のMax/MinをセッションごとにTP/TN/Allで求める
% R_D1pre_TNrand_sub_All_Max = max(max(R_D1pre_TNrand_sub));
% R_D1pre_TNrand_sub_All_Min = min(min(R_D1pre_TNrand_sub));
% %save the matrix
% writematrix(R_D1pre_TNrand_sub,[mouseID,dinfo,'_D1pre_dFCorr',num2str(i),'.txt'],'delimiter','\t');
% % List_R_D1pre_TNrand_sub_TNrand_average(i,:) = R_D1pre_TNrand_sub_TN_average;
% List_R_D1pre_TNrand_sub_average(i,:) = R_D1pre_TNrand_sub_All_average;
% List_R_D1pre_TNrand_sub_Max(i,:) = R_D1pre_TNrand_sub_All_Max;
% List_R_D1pre_TNrand_sub_Min(i,:) = R_D1pre_TNrand_sub_All_Min;
% %%%%% Day1 context %%%%%
% [R_D1ctx_TNrand,P_D1ctx_TNrand] = corrcoef(dFData_TNrand_D1ctx);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
%      
% % subtact self correlation (1)
% R_D1ctx_TNrand_sub = R_D1ctx_TNrand - eye(size(R_D1ctx_TNrand));
% 
% %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% R_D1ctx_TNrand_sub_0 = R_D1ctx_TNrand_sub;
% R_D1ctx_TNrand_sub_0(R_D1ctx_TNrand_sub_0~=0) = 1;
% 
% % R_D1ctx_TNrand_sub_0_TN_sumsum = sum(sum(R_D1ctx_TNrand_sub_0(TN_first:end,TN_first:end)));
% R_D1ctx_TNrand_sub_0_All_sumsum = sum(sum(R_D1ctx_TNrand_sub_0));
% 
% %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% % R_D1ctx_TNrand_sub_TN_sumsum = sum(sum(R_D1ctx_TNrand_sub(TN_first:end,TN_first:end)));
% R_D1ctx_TNrand_sub_All_sumsum = sum(sum(R_D1ctx_TNrand_sub));
% 
% %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% % R_D1ctx_TNrand_sub_TN_average = R_D1ctx_TNrand_sub_TN_sumsum/R_D1ctx_TNrand_sub_0_TN_sumsum;
% R_D1ctx_TNrand_sub_All_average = R_D1ctx_TNrand_sub_All_sumsum/R_D1ctx_TNrand_sub_0_All_sumsum;
% 
% %対角以外のますの数のMax/MinをセッションごとにTP/TN/Allで求める
% R_D1ctx_TNrand_sub_All_Max = max(max(R_D1ctx_TNrand_sub));
% R_D1ctx_TNrand_sub_All_Min = min(min(R_D1ctx_TNrand_sub));
% 
% %save the matrix
% writematrix(R_D1ctx_TNrand_sub,[mouseID,dinfo,'_D1ctx_dFCorr',num2str(i),'.txt'],'delimiter','\t');
% 
% 
% % List_R_D1ctx_TNrand_sub_TNrand_average(i,:) = R_D1ctx_TNrand_sub_TN_average;
% List_R_D1ctx_TNrand_sub_average(i,:) = R_D1ctx_TNrand_sub_All_average;
% List_R_D1ctx_TNrand_sub_Max(i,:) = R_D1ctx_TNrand_sub_All_Max;
% List_R_D1ctx_TNrand_sub_Min(i,:) = R_D1ctx_TNrand_sub_All_Min;
% % %%%%% Day1 after %%%%%
% % [R_D1after_TNrand,P_D1after_TNrand] = corrcoef(dFData_TNrand_D1after);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
% % % subtact self correlation (1)
% % R_D1after_TNrand_sub = R_D1after_TNrand - eye(size(R_D1after_TNrand));
% % %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% % R_D1after_TNrand_sub_0 = R_D1after_TNrand_sub;
% % R_D1after_TNrand_sub_0(R_D1after_TNrand_sub_0~=0) = 1;
% % % R_D1after_TNrand_sub_0_TN_sumsum = sum(sum(R_D1after_TNrand_sub_0(TN_first:end,TN_first:end)));
% % R_D1after_TNrand_sub_0_All_sumsum = sum(sum(R_D1after_TNrand_sub_0));
% % %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% % % R_D1after_TNrand_sub_TN_sumsum = sum(sum(R_D1after_TNrand_sub(TN_first:end,TN_first:end)));
% % R_D1after_TNrand_sub_All_sumsum = sum(sum(R_D1after_TNrand_sub));
% % %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% % % R_D1after_TNrand_sub_TN_average = R_D1after_TNrand_sub_TN_sumsum/R_D1after_TNrand_sub_0_TN_sumsum;
% % R_D1after_TNrand_sub_All_average = R_D1after_TNrand_sub_All_sumsum/R_D1after_TNrand_sub_0_All_sumsum;
% % %対角以外のますの数のMax/MinをセッションごとにTP/TN/Allで求める
% % R_D1after_TNrand_sub_All_Max = max(max(R_D1after_TNrand_sub));
% % R_D1after_TNrand_sub_All_Min = min(min(R_D1after_TNrand_sub));
% % %save the matrix
% % writematrix(R_D1after_TNrand_sub,[mouseID,dinfo,'_D1after_dFCorr',num2str(i),'.txt'],'delimiter','\t');
% % % List_R_D1after_TNrand_sub_TNrand_average(i,:) = R_D1after_TNrand_sub_TN_average;
% % List_R_D1after_TNrand_sub_average(i,:) = R_D1after_TNrand_sub_All_average;
% % List_R_D1after_TNrand_sub_Max(i,:) = R_D1after_TNrand_sub_All_Max;
% % List_R_D1after_TNrand_sub_Min(i,:) = R_D1after_TNrand_sub_All_Min;
% %%%%% Day1 post %%%%%
% [R_D1post_TNrand,P_D1post_TNrand] = corrcoef(dFData_TNrand_D1post);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
% % subtact self correlation (1)
% R_D1post_TNrand_sub = R_D1post_TNrand - eye(size(R_D1post_TNrand));
% %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% R_D1post_TNrand_sub_0 = R_D1post_TNrand_sub;
% R_D1post_TNrand_sub_0(R_D1post_TNrand_sub_0~=0) = 1;
% % R_D1post_TNrand_sub_0_TN_sumsum = sum(sum(R_D1post_TNrand_sub_0(TN_first:end,TN_first:end)));
% R_D1post_TNrand_sub_0_All_sumsum = sum(sum(R_D1post_TNrand_sub_0));
% %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% % R_D1post_TNrand_sub_TN_sumsum = sum(sum(R_D1post_TNrand_sub(TN_first:end,TN_first:end)));
% R_D1post_TNrand_sub_All_sumsum = sum(sum(R_D1post_TNrand_sub));
% %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% % R_D1post_TNrand_sub_TN_average = R_D1post_TNrand_sub_TN_sumsum/R_D1post_TNrand_sub_0_TN_sumsum;
% R_D1post_TNrand_sub_All_average = R_D1post_TNrand_sub_All_sumsum/R_D1post_TNrand_sub_0_All_sumsum;
% %対角以外のますの数のMax/MinをセッションごとにTP/TN/Allで求める
% R_D1post_TNrand_sub_All_Max = max(max(R_D1post_TNrand_sub));
% R_D1post_TNrand_sub_All_Min = min(min(R_D1post_TNrand_sub));
% %save the matrix
% writematrix(R_D1post_TNrand_sub,[mouseID,dinfo,'_D1post_dFCorr',num2str(i),'.txt'],'delimiter','\t');
% % List_R_D1post_TNrand_sub_TNrand_average(i,:) = R_D1post_TNrand_sub_TN_average;
% List_R_D1post_TNrand_sub_average(i,:) = R_D1post_TNrand_sub_All_average;
% List_R_D1post_TNrand_sub_Max(i,:) = R_D1post_TNrand_sub_All_Max;
% List_R_D1post_TNrand_sub_Min(i,:) = R_D1post_TNrand_sub_All_Min;
% 
% %%%%% Day2 pre %%%%%
% [R_D2pre_TNrand,P_D2pre_TNrand] = corrcoef(dFData_TNrand_D2pre);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
% % subtact self correlation (1)
% R_D2pre_TNrand_sub = R_D2pre_TNrand - eye(size(R_D2pre_TNrand));
% %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% R_D2pre_TNrand_sub_0 = R_D2pre_TNrand_sub;
% R_D2pre_TNrand_sub_0(R_D2pre_TNrand_sub_0~=0) = 1;
% % R_D2pre_TNrand_sub_0_TN_sumsum = sum(sum(R_D2pre_TNrand_sub_0(TN_first:end,TN_first:end)));
% R_D2pre_TNrand_sub_0_All_sumsum = sum(sum(R_D2pre_TNrand_sub_0));
% %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% % R_D2pre_TNrand_sub_TN_sumsum = sum(sum(R_D2pre_TNrand_sub(TN_first:end,TN_first:end)));
% R_D2pre_TNrand_sub_All_sumsum = sum(sum(R_D2pre_TNrand_sub));
% %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% % R_D2pre_TNrand_sub_TN_average = R_D2pre_TNrand_sub_TN_sumsum/R_D2pre_TNrand_sub_0_TN_sumsum;
% R_D2pre_TNrand_sub_All_average = R_D2pre_TNrand_sub_All_sumsum/R_D2pre_TNrand_sub_0_All_sumsum;
% %対角以外のますの数のMax/MinをセッションごとにTP/TN/Allで求める
% R_D2pre_TNrand_sub_All_Max = max(max(R_D2pre_TNrand_sub));
% R_D2pre_TNrand_sub_All_Min = min(min(R_D2pre_TNrand_sub));
% %save the matrix
% writematrix(R_D2pre_TNrand_sub,[mouseID,dinfo,'_D2pre_dFCorr',num2str(i),'.txt'],'delimiter','\t');
% % List_R_D2pre_TNrand_sub_TNrand_average(i,:) = R_D2pre_TNrand_sub_TN_average;
% List_R_D2pre_TNrand_sub_average(i,:) = R_D2pre_TNrand_sub_All_average;
% List_R_D2pre_TNrand_sub_Max(i,:) = R_D2pre_TNrand_sub_All_Max;
% List_R_D2pre_TNrand_sub_Min(i,:) = R_D2pre_TNrand_sub_All_Min;
% %%%%% Day2 context %%%%%
% [R_D2ctx_TNrand,P_D2ctx_TNrand] = corrcoef(dFData_TNrand_D2ctx);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
%      
% % subtact self correlation (1)
% R_D2ctx_TNrand_sub = R_D2ctx_TNrand - eye(size(R_D2ctx_TNrand));
% 
% %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% R_D2ctx_TNrand_sub_0 = R_D2ctx_TNrand_sub;
% R_D2ctx_TNrand_sub_0(R_D2ctx_TNrand_sub_0~=0) = 1;
% 
% % R_D2ctx_TNrand_sub_0_TN_sumsum = sum(sum(R_D2ctx_TNrand_sub_0(TN_first:end,TN_first:end)));
% R_D2ctx_TNrand_sub_0_All_sumsum = sum(sum(R_D2ctx_TNrand_sub_0));
% 
% %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% % R_D2ctx_TNrand_sub_TN_sumsum = sum(sum(R_D2ctx_TNrand_sub(TN_first:end,TN_first:end)));
% R_D2ctx_TNrand_sub_All_sumsum = sum(sum(R_D2ctx_TNrand_sub));
% 
% %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% % R_D2ctx_TNrand_sub_TN_average = R_D2ctx_TNrand_sub_TN_sumsum/R_D2ctx_TNrand_sub_0_TN_sumsum;
% R_D2ctx_TNrand_sub_All_average = R_D2ctx_TNrand_sub_All_sumsum/R_D2ctx_TNrand_sub_0_All_sumsum;
% 
% %対角以外のますの数のMax/MinをセッションごとにTP/TN/Allで求める
% R_D2ctx_TNrand_sub_All_Max = max(max(R_D2ctx_TNrand_sub));
% R_D2ctx_TNrand_sub_All_Min = min(min(R_D2ctx_TNrand_sub));
% 
% %save the matrix
% writematrix(R_D2ctx_TNrand_sub,[mouseID,dinfo,'_D2ctx_dFCorr',num2str(i),'.txt'],'delimiter','\t');
% 
% 
% % List_R_D2ctx_TNrand_sub_TNrand_average(i,:) = R_D2ctx_TNrand_sub_TN_average;
% List_R_D2ctx_TNrand_sub_average(i,:) = R_D2ctx_TNrand_sub_All_average;
% List_R_D2ctx_TNrand_sub_Max(i,:) = R_D2ctx_TNrand_sub_All_Max;
% List_R_D2ctx_TNrand_sub_Min(i,:) = R_D2ctx_TNrand_sub_All_Min;
% % %%%%% Day2 after %%%%%
% % [R_D2after_TNrand,P_D2after_TNrand] = corrcoef(dFData_TNrand_D2after);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
% % % subtact self correlation (1)
% % R_D2after_TNrand_sub = R_D2after_TNrand - eye(size(R_D2after_TNrand));
% % %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% % R_D2after_TNrand_sub_0 = R_D2after_TNrand_sub;
% % R_D2after_TNrand_sub_0(R_D2after_TNrand_sub_0~=0) = 1;
% % % R_D2after_TNrand_sub_0_TN_sumsum = sum(sum(R_D2after_TNrand_sub_0(TN_first:end,TN_first:end)));
% % R_D2after_TNrand_sub_0_All_sumsum = sum(sum(R_D2after_TNrand_sub_0));
% % %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% % % R_D2after_TNrand_sub_TN_sumsum = sum(sum(R_D2after_TNrand_sub(TN_first:end,TN_first:end)));
% % R_D2after_TNrand_sub_All_sumsum = sum(sum(R_D2after_TNrand_sub));
% % %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% % % R_D2after_TNrand_sub_TN_average = R_D2after_TNrand_sub_TN_sumsum/R_D2after_TNrand_sub_0_TN_sumsum;
% % R_D2after_TNrand_sub_All_average = R_D2after_TNrand_sub_All_sumsum/R_D2after_TNrand_sub_0_All_sumsum;
% % %対角以外のますの数のMax/MinをセッションごとにTP/TN/Allで求める
% % R_D2after_TNrand_sub_All_Max = max(max(R_D2after_TNrand_sub));
% % R_D2after_TNrand_sub_All_Min = min(min(R_D2after_TNrand_sub));
% % %save the matrix
% % writematrix(R_D2after_TNrand_sub,[mouseID,dinfo,'_D2after_dFCorr',num2str(i),'.txt'],'delimiter','\t');
% % % List_R_D2after_TNrand_sub_TNrand_average(i,:) = R_D2after_TNrand_sub_TN_average;
% % List_R_D2after_TNrand_sub_average(i,:) = R_D2after_TNrand_sub_All_average;
% % List_R_D2after_TNrand_sub_Max(i,:) = R_D2after_TNrand_sub_All_Max;
% % List_R_D2after_TNrand_sub_Min(i,:) = R_D2after_TNrand_sub_All_Min;
% %%%%% Day2 post %%%%%
% [R_D2post_TNrand,P_D2post_TNrand] = corrcoef(dFData_TNrand_D2post);%https://jp.mathworks.com/help/matlab/ref/corrcoef.html, Rは相関係数の行列Pは統計値 (0.05以下なら有意だが、、、)
% % subtact self correlation (1)
% R_D2post_TNrand_sub = R_D2post_TNrand - eye(size(R_D2post_TNrand));
% %対角の数を求める > 対角以外を1に置換 > 対角以外のますの数をセッションごとにTP/TN/Allで求める
% R_D2post_TNrand_sub_0 = R_D2post_TNrand_sub;
% R_D2post_TNrand_sub_0(R_D2post_TNrand_sub_0~=0) = 1;
% % R_D2post_TNrand_sub_0_TN_sumsum = sum(sum(R_D2post_TNrand_sub_0(TN_first:end,TN_first:end)));
% R_D2post_TNrand_sub_0_All_sumsum = sum(sum(R_D2post_TNrand_sub_0));
% %対角以外のますの数の和をセッションごとにTP/TN/Allで求める
% % R_D2post_TNrand_sub_TN_sumsum = sum(sum(R_D2post_TNrand_sub(TN_first:end,TN_first:end)));
% R_D2post_TNrand_sub_All_sumsum = sum(sum(R_D2post_TNrand_sub));
% %対角以外のますの数のaverage(和/数)をセッションごとにTP/TN/Allで求める
% % R_D2post_TNrand_sub_TN_average = R_D2post_TNrand_sub_TN_sumsum/R_D2post_TNrand_sub_0_TN_sumsum;
% R_D2post_TNrand_sub_All_average = R_D2post_TNrand_sub_All_sumsum/R_D2post_TNrand_sub_0_All_sumsum;
% %対角以外のますの数のMax/MinをセッションごとにTP/TN/Allで求める
% R_D2post_TNrand_sub_All_Max = max(max(R_D2post_TNrand_sub));
% R_D2post_TNrand_sub_All_Min = min(min(R_D2post_TNrand_sub));
% %save the matrix
% writematrix(R_D2post_TNrand_sub,[mouseID,dinfo,'_D2post_dFCorr',num2str(i),'.txt'],'delimiter','\t');
% % List_R_D2post_TNrand_sub_TNrand_average(i,:) = R_D2post_TNrand_sub_TN_average;
% List_R_D2post_TNrand_sub_average(i,:) = R_D2post_TNrand_sub_All_average;
% List_R_D2post_TNrand_sub_Max(i,:) = R_D2post_TNrand_sub_All_Max;
% List_R_D2post_TNrand_sub_Min(i,:) = R_D2post_TNrand_sub_All_Min;
% 
% end
% toc
% % summary of results
% 
% % % day1 pre/context/after/post ver.
% % Result_R_average = table(Cell, R_D1pre_TPTNTotal_average,R_D1ctx_TNrand_TPTNTotal_average,R_D1after_TPTNTotal_average,R_D1post_TPTNTotal_average);
% % %     % afterなしver. day1 pre/context/post ver
% % %     Result_R_average = table(Cell, R_D1pre_TPTNTotal_average,R_D1ctx_TNrand_TPTNTotal_average,R_D1post_TPTNTotal_average);
% % writetable(Result_R_average,[mouseID,dinfo,'_Day1_CellCorr_Tnrand_',dtime,'.txt'], 'delimiter','\t');%
% 
% %%%%% Day1 pre %%%%%
% session = char('D1pre_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(List_R_D1pre_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMax_');
% writematrix(List_R_D1pre_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMin_');
% writematrix(List_R_D1pre_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% 
% Mean_List_R_D1pre_TNrand_sub_average = mean(List_R_D1pre_TNrand_sub_average);
% Mean_List_R_D1pre_TNrand_sub_Max = mean(List_R_D1pre_TNrand_sub_Max);
% Mean_List_R_D1pre_TNrand_sub_Min = mean(List_R_D1pre_TNrand_sub_Min);
% 
% Cell = [{'Tom+'};{'Tom- random'};{'Tom-'};{'All'}];
% R_D1pre_TPTNrandTNTotal_average = [R_D1pre_sub_TP_average;Mean_List_R_D1pre_TNrand_sub_average;R_D1pre_sub_TN_average;R_D1pre_sub_All_average];
% R_D1pre_TPTNrandTNTotal_Max = [R_D1pre_sub_TP_Max;Mean_List_R_D1pre_TNrand_sub_Max;R_D1pre_sub_TN_Max;R_D1pre_sub_All_Max];
% R_D1pre_TPTNrandTNTotal_Min = [R_D1pre_sub_TP_Min;Mean_List_R_D1pre_TNrand_sub_Min;R_D1pre_sub_TN_Min;R_D1pre_sub_All_Min];
% 
% session = char('D1pre_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(R_D1pre_TPTNrandTNTotal_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMax_');
% writematrix(R_D1pre_TPTNrandTNTotal_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMin_');
% writematrix(R_D1pre_TPTNrandTNTotal_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% 
%     x =  categorical({'Tom+','Tom- random','Tom- ','All'});
%     x = reordercats(x,{'Tom+','Tom- random','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
%     d_content = char('dFcorrAv_');
%     figure;bar(x,R_D1pre_TPTNrandTNTotal_average);
%     title("Day1 pre");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
%     
%     d_content = char('dFcorrMax_');
%     figure;bar(x,R_D1pre_TPTNrandTNTotal_Max);
%     title("Day1 pre");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorrMin_');
%     figure;bar(x,R_D1pre_TPTNrandTNTotal_Min);
%     title("Day1 pre");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     % figure (correlation), example
%     cmax = 0.10;cmin = -0.10;
% 
%     session = char('D1pre_');
%     celltype = char('TNrand_');
%     d_content = char('dFcorrAv_');
%     
%     figure;imagesc(R_D1pre_TNrand)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 pre")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
%     d_content = char('dFcorrAv_sub_');
%     figure;imagesc(R_D1pre_TNrand_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 pre")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
% 
%     %%%%%%%%%%%%%%%%　別のセクションから　%%%%%%%%%%%%%%%%%%%%%
% 
% figure;imagesc(R_D1pre(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day1 pre")  
%     saveas(gcf, [mouseID,'_D1pre_TP_CorrCell',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D1pre_TP_CorrCell',dtime], 'tif');
% 
%     figure;imagesc(R_D1pre_sub(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day1 pre")  
%     saveas(gcf, [mouseID,'_D1pre_TP_CorrCell_sub',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D1pre_TP_CorrCell_sub',dtime], 'tif');
% 
% %%%%% Day1 context %%%%%
% session = char('D1context_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(List_R_D1ctx_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMax_');
% writematrix(List_R_D1ctx_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMin_');
% writematrix(List_R_D1ctx_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% 
% Mean_List_R_D1ctx_TNrand_sub_average = mean(List_R_D1ctx_TNrand_sub_average);
% Mean_List_R_D1ctx_TNrand_sub_Max = mean(List_R_D1ctx_TNrand_sub_Max);
% Mean_List_R_D1ctx_TNrand_sub_Min = mean(List_R_D1ctx_TNrand_sub_Min);
% 
% Cell = [{'Tom+'};{'Tom- random'};{'Tom-'};{'All'}];
% R_D1ctx_TPTNrandTNTotal_average = [R_D1ctx_sub_TP_average;Mean_List_R_D1ctx_TNrand_sub_average;R_D1ctx_sub_TN_average;R_D1ctx_sub_All_average];
% R_D1ctx_TPTNrandTNTotal_Max = [R_D1ctx_sub_TP_Max;Mean_List_R_D1ctx_TNrand_sub_Max;R_D1ctx_sub_TN_Max;R_D1ctx_sub_All_Max];
% R_D1ctx_TPTNrandTNTotal_Min = [R_D1ctx_sub_TP_Min;Mean_List_R_D1ctx_TNrand_sub_Min;R_D1ctx_sub_TN_Min;R_D1ctx_sub_All_Min];
% 
% session = char('D1context_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(R_D1ctx_TPTNrandTNTotal_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMax_');
% writematrix(R_D1ctx_TPTNrandTNTotal_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMin_');
% writematrix(R_D1ctx_TPTNrandTNTotal_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% 
%     x =  categorical({'Tom+','Tom- random','Tom- ','All'});
%     x = reordercats(x,{'Tom+','Tom- random','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
%     d_content = char('dFcorrAv_');
%     figure;bar(x,R_D1ctx_TPTNrandTNTotal_average);
%     title("Day1 context");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
%     
%     d_content = char('dFcorrMax_');
%     figure;bar(x,R_D1ctx_TPTNrandTNTotal_Max);
%     title("Day1 context");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorrMin_');
%     figure;bar(x,R_D1ctx_TPTNrandTNTotal_Min);
%     title("Day1 context");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     % figure (correlation), example
%     cmax = 0.10;cmin = -0.10;
% 
%     session = char('D1context_');
%     celltype = char('TNrand_');
%     d_content = char('dFcorrAv_');
%     
%     figure;imagesc(R_D1ctx_TNrand)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 context")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
%     d_content = char('dFcorrAv_sub_');
%     figure;imagesc(R_D1ctx_TNrand_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 context")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
% 
%     %%%%%%%%%%%%%%%%　別のセクションから　%%%%%%%%%%%%%%%%%%%%%
% 
% figure;imagesc(R_D1ctx(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day1 context")  
%     saveas(gcf, [mouseID,'_D1context_TP_CorrCell',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D1context_TP_CorrCell',dtime], 'tif');
% 
%     figure;imagesc(R_D1ctx_sub(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day1 context")  
%     saveas(gcf, [mouseID,'_D1context_TP_CorrCell_sub',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D1context_TP_CorrCell_sub',dtime], 'tif');
% % %%%%% Day1 after %%%%%
% % session = char('D1after_');
% % celltype = char('TNrand_');
% % d_content = char('dFcorrAv_');
% % writematrix(List_R_D1after_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% % d_content = char('dFcorrMax_');
% % writematrix(List_R_D1after_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% % d_content = char('dFcorrMin_');
% % writematrix(List_R_D1after_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% % 
% % Mean_List_R_D1after_TNrand_sub_average = mean(List_R_D1after_TNrand_sub_average);
% % Mean_List_R_D1after_TNrand_sub_Max = mean(List_R_D1after_TNrand_sub_Max);
% % Mean_List_R_D1after_TNrand_sub_Min = mean(List_R_D1after_TNrand_sub_Min);
% % 
% % Cell = [{'Tom+'};{'Tom- random'};{'Tom-'};{'All'}];
% % R_D1after_TPTNrandTNTotal_average = [R_D1after_sub_TP_average;Mean_List_R_D1after_TNrand_sub_average;R_D1after_sub_TN_average;R_D1after_sub_All_average];
% % R_D1after_TPTNrandTNTotal_Max = [R_D1after_sub_TP_Max;Mean_List_R_D1after_TNrand_sub_Max;R_D1after_sub_TN_Max;R_D1after_sub_All_Max];
% % R_D1after_TPTNrandTNTotal_Min = [R_D1after_sub_TP_Min;Mean_List_R_D1after_TNrand_sub_Min;R_D1after_sub_TN_Min;R_D1after_sub_All_Min];
% % 
% % session = char('D1after_');
% % celltype = char('TNrand_');
% % d_content = char('dFcorrAv_');
% % writematrix(R_D1after_TPTNrandTNTotal_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% % d_content = char('dFcorrMax_');
% % writematrix(R_D1after_TPTNrandTNTotal_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% % d_content = char('dFcorrMin_');
% % writematrix(R_D1after_TPTNrandTNTotal_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% % 
% %     x =  categorical({'Tom+','Tom- random','Tom- ','All'});
% %     x = reordercats(x,{'Tom+','Tom- random','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
% %     d_content = char('dFcorrAv_');
% %     figure;bar(x,R_D1after_TPTNrandTNTotal_average);
% %     title("Day1 after");
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% %     
% %     d_content = char('dFcorrMax_');
% %     figure;bar(x,R_D1after_TPTNrandTNTotal_Max);
% %     title("Day1 after");
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% % 
% %     d_content = char('dFcorrMin_');
% %     figure;bar(x,R_D1after_TPTNrandTNTotal_Min);
% %     title("Day1 after");
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% % 
% %     % figure (correlation), example
% %     cmax = 0.10;cmin = -0.10;
% % 
% %     session = char('D1after_');
% %     celltype = char('TNrand_');
% %     d_content = char('dFcorrAv_');
% %     
% %     figure;imagesc(R_D1after_TNrand)
% %     colorbar;caxis([cmin cmax]);
% %     title("Day1 after")  
% %     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% % 
% %     d_content = char('dFcorrAv_sub_');
% %     figure;imagesc(R_D1after_TNrand_sub)
% %     colorbar;caxis([cmin cmax]);
% %     title("Day1 after")  
% %     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% % 
% % 
% %     %%%%%%%%%%%%%%%%　別のセクションから　%%%%%%%%%%%%%%%%%%%%%
% % 
% % figure;imagesc(R_D1after(1:TP_last,1:TP_last))
% %     colorbar;caxis([cmin cmax]);
% %     title("Day1 after")  
% %     saveas(gcf, [mouseID,'_D1after_TP_CorrCell',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,'_D1after_TP_CorrCell',dtime], 'tif');
% % 
% %     figure;imagesc(R_D1after_sub(1:TP_last,1:TP_last))
% %     colorbar;caxis([cmin cmax]);
% %     title("Day1 after")  
% %     saveas(gcf, [mouseID,'_D1after_TP_CorrCell_sub',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,'_D1after_TP_CorrCell_sub',dtime], 'tif');
% 
% %%%%% Day1 post %%%%%
% session = char('D1post_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(List_R_D1post_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMax_');
% writematrix(List_R_D1post_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMin_');
% writematrix(List_R_D1post_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% 
% Mean_List_R_D1post_TNrand_sub_average = mean(List_R_D1post_TNrand_sub_average);
% Mean_List_R_D1post_TNrand_sub_Max = mean(List_R_D1post_TNrand_sub_Max);
% Mean_List_R_D1post_TNrand_sub_Min = mean(List_R_D1post_TNrand_sub_Min);
% 
% Cell = [{'Tom+'};{'Tom- random'};{'Tom-'};{'All'}];
% R_D1post_TPTNrandTNTotal_average = [R_D1post_sub_TP_average;Mean_List_R_D1post_TNrand_sub_average;R_D1post_sub_TN_average;R_D1post_sub_All_average];
% R_D1post_TPTNrandTNTotal_Max = [R_D1post_sub_TP_Max;Mean_List_R_D1post_TNrand_sub_Max;R_D1post_sub_TN_Max;R_D1post_sub_All_Max];
% R_D1post_TPTNrandTNTotal_Min = [R_D1post_sub_TP_Min;Mean_List_R_D1post_TNrand_sub_Min;R_D1post_sub_TN_Min;R_D1post_sub_All_Min];
% 
% session = char('D1post_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(R_D1post_TPTNrandTNTotal_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMax_');
% writematrix(R_D1post_TPTNrandTNTotal_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMin_');
% writematrix(R_D1post_TPTNrandTNTotal_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% 
%     x =  categorical({'Tom+','Tom- random','Tom- ','All'});
%     x = reordercats(x,{'Tom+','Tom- random','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
%     d_content = char('dFcorrAv_');
%     figure;bar(x,R_D1post_TPTNrandTNTotal_average);
%     title("Day1 post");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
%     
%     d_content = char('dFcorrMax_');
%     figure;bar(x,R_D1post_TPTNrandTNTotal_Max);
%     title("Day1 post");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorrMin_');
%     figure;bar(x,R_D1post_TPTNrandTNTotal_Min);
%     title("Day1 post");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     % figure (correlation), example
%     cmax = 0.10;cmin = -0.10;
% 
%     session = char('D1post_');
%     celltype = char('TNrand_');
%     d_content = char('dFcorrAv_');
%     
%     figure;imagesc(R_D1post_TNrand)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 post")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
%     d_content = char('dFcorrAv_sub_');
%     figure;imagesc(R_D1post_TNrand_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day1 post")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
% 
%     %%%%%%%%%%%%%%%%　別のセクションから　%%%%%%%%%%%%%%%%%%%%%
% 
% figure;imagesc(R_D1post(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day1 post")  
%     saveas(gcf, [mouseID,'_D1post_TP_CorrCell',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D1post_TP_CorrCell',dtime], 'tif');
% 
%     figure;imagesc(R_D1post_sub(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day1 post")  
%     saveas(gcf, [mouseID,'_D1post_TP_CorrCell_sub',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D1post_TP_CorrCell_sub',dtime], 'tif');
% %%%%% Day2 pre %%%%%
% session = char('D2pre_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(List_R_D2pre_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMax_');
% writematrix(List_R_D2pre_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMin_');
% writematrix(List_R_D2pre_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% 
% Mean_List_R_D2pre_TNrand_sub_average = mean(List_R_D2pre_TNrand_sub_average);
% Mean_List_R_D2pre_TNrand_sub_Max = mean(List_R_D2pre_TNrand_sub_Max);
% Mean_List_R_D2pre_TNrand_sub_Min = mean(List_R_D2pre_TNrand_sub_Min);
% 
% Cell = [{'Tom+'};{'Tom- random'};{'Tom-'};{'All'}];
% R_D2pre_TPTNrandTNTotal_average = [R_D2pre_sub_TP_average;Mean_List_R_D2pre_TNrand_sub_average;R_D2pre_sub_TN_average;R_D2pre_sub_All_average];
% R_D2pre_TPTNrandTNTotal_Max = [R_D2pre_sub_TP_Max;Mean_List_R_D2pre_TNrand_sub_Max;R_D2pre_sub_TN_Max;R_D2pre_sub_All_Max];
% R_D2pre_TPTNrandTNTotal_Min = [R_D2pre_sub_TP_Min;Mean_List_R_D2pre_TNrand_sub_Min;R_D2pre_sub_TN_Min;R_D2pre_sub_All_Min];
% 
% session = char('D2pre_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(R_D2pre_TPTNrandTNTotal_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMax_');
% writematrix(R_D2pre_TPTNrandTNTotal_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMin_');
% writematrix(R_D2pre_TPTNrandTNTotal_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% 
%     x =  categorical({'Tom+','Tom- random','Tom- ','All'});
%     x = reordercats(x,{'Tom+','Tom- random','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
%     d_content = char('dFcorrAv_');
%     figure;bar(x,R_D2pre_TPTNrandTNTotal_average);
%     title("Day2 pre");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
%     
%     d_content = char('dFcorrMax_');
%     figure;bar(x,R_D2pre_TPTNrandTNTotal_Max);
%     title("Day2 pre");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorrMin_');
%     figure;bar(x,R_D2pre_TPTNrandTNTotal_Min);
%     title("Day2 pre");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     % figure (correlation), example
%     cmax = 0.10;cmin = -0.10;
% 
%     session = char('D2pre_');
%     celltype = char('TNrand_');
%     d_content = char('dFcorrAv_');
%     
%     figure;imagesc(R_D2pre_TNrand)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 pre")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
%     d_content = char('dFcorrAv_sub_');
%     figure;imagesc(R_D2pre_TNrand_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 pre")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
% 
%     %%%%%%%%%%%%%%%%　別のセクションから　%%%%%%%%%%%%%%%%%%%%%
% 
% figure;imagesc(R_D2pre(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day2 pre")  
%     saveas(gcf, [mouseID,'_D2pre_TP_CorrCell',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D2pre_TP_CorrCell',dtime], 'tif');
% 
%     figure;imagesc(R_D2pre_sub(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day2 pre")  
%     saveas(gcf, [mouseID,'_D2pre_TP_CorrCell_sub',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D2pre_TP_CorrCell_sub',dtime], 'tif');
% 
% %%%%% Day2 context %%%%%
% session = char('D2context_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(List_R_D2ctx_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMax_');
% writematrix(List_R_D2ctx_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMin_');
% writematrix(List_R_D2ctx_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% 
% Mean_List_R_D2ctx_TNrand_sub_average = mean(List_R_D2ctx_TNrand_sub_average);
% Mean_List_R_D2ctx_TNrand_sub_Max = mean(List_R_D2ctx_TNrand_sub_Max);
% Mean_List_R_D2ctx_TNrand_sub_Min = mean(List_R_D2ctx_TNrand_sub_Min);
% 
% Cell = [{'Tom+'};{'Tom- random'};{'Tom-'};{'All'}];
% R_D2ctx_TPTNrandTNTotal_average = [R_D2ctx_sub_TP_average;Mean_List_R_D2ctx_TNrand_sub_average;R_D2ctx_sub_TN_average;R_D2ctx_sub_All_average];
% R_D2ctx_TPTNrandTNTotal_Max = [R_D2ctx_sub_TP_Max;Mean_List_R_D2ctx_TNrand_sub_Max;R_D2ctx_sub_TN_Max;R_D2ctx_sub_All_Max];
% R_D2ctx_TPTNrandTNTotal_Min = [R_D2ctx_sub_TP_Min;Mean_List_R_D2ctx_TNrand_sub_Min;R_D2ctx_sub_TN_Min;R_D2ctx_sub_All_Min];
% 
% session = char('D2context_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(R_D2ctx_TPTNrandTNTotal_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMax_');
% writematrix(R_D2ctx_TPTNrandTNTotal_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMin_');
% writematrix(R_D2ctx_TPTNrandTNTotal_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% 
%     x =  categorical({'Tom+','Tom- random','Tom- ','All'});
%     x = reordercats(x,{'Tom+','Tom- random','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
%     d_content = char('dFcorrAv_');
%     figure;bar(x,R_D2ctx_TPTNrandTNTotal_average);
%     title("Day2 context");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
%     
%     d_content = char('dFcorrMax_');
%     figure;bar(x,R_D2ctx_TPTNrandTNTotal_Max);
%     title("Day2 context");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorrMin_');
%     figure;bar(x,R_D2ctx_TPTNrandTNTotal_Min);
%     title("Day2 context");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     % figure (correlation), example
%     cmax = 0.10;cmin = -0.10;
% 
%     session = char('D2context_');
%     celltype = char('TNrand_');
%     d_content = char('dFcorrAv_');
%     
%     figure;imagesc(R_D2ctx_TNrand)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 context")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
%     d_content = char('dFcorrAv_sub_');
%     figure;imagesc(R_D2ctx_TNrand_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 context")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
% 
%     %%%%%%%%%%%%%%%%　別のセクションから　%%%%%%%%%%%%%%%%%%%%%
% 
% figure;imagesc(R_D2ctx(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day2 context")  
%     saveas(gcf, [mouseID,'_D2context_TP_CorrCell',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D2context_TP_CorrCell',dtime], 'tif');
% 
%     figure;imagesc(R_D2ctx_sub(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day2 context")  
%     saveas(gcf, [mouseID,'_D2context_TP_CorrCell_sub',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D2context_TP_CorrCell_sub',dtime], 'tif');
% % %%%%% Day2 after %%%%%
% % session = char('D2after_');
% % celltype = char('TNrand_');
% % d_content = char('dFcorrAv_');
% % writematrix(List_R_D2after_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% % d_content = char('dFcorrMax_');
% % writematrix(List_R_D2after_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% % d_content = char('dFcorrMin_');
% % writematrix(List_R_D2after_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% % 
% % Mean_List_R_D2after_TNrand_sub_average = mean(List_R_D2after_TNrand_sub_average);
% % Mean_List_R_D2after_TNrand_sub_Max = mean(List_R_D2after_TNrand_sub_Max);
% % Mean_List_R_D2after_TNrand_sub_Min = mean(List_R_D2after_TNrand_sub_Min);
% % 
% % Cell = [{'Tom+'};{'Tom- random'};{'Tom-'};{'All'}];
% % R_D2after_TPTNrandTNTotal_average = [R_D2after_sub_TP_average;Mean_List_R_D2after_TNrand_sub_average;R_D2after_sub_TN_average;R_D2after_sub_All_average];
% % R_D2after_TPTNrandTNTotal_Max = [R_D2after_sub_TP_Max;Mean_List_R_D2after_TNrand_sub_Max;R_D2after_sub_TN_Max;R_D2after_sub_All_Max];
% % R_D2after_TPTNrandTNTotal_Min = [R_D2after_sub_TP_Min;Mean_List_R_D2after_TNrand_sub_Min;R_D2after_sub_TN_Min;R_D2after_sub_All_Min];
% % 
% % session = char('D2after_');
% % celltype = char('TNrand_');
% % d_content = char('dFcorrAv_');
% % writematrix(R_D2after_TPTNrandTNTotal_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% % d_content = char('dFcorrMax_');
% % writematrix(R_D2after_TPTNrandTNTotal_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% % d_content = char('dFcorrMin_');
% % writematrix(R_D2after_TPTNrandTNTotal_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% % 
% %     x =  categorical({'Tom+','Tom- random','Tom- ','All'});
% %     x = reordercats(x,{'Tom+','Tom- random','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
% %     d_content = char('dFcorrAv_');
% %     figure;bar(x,R_D2after_TPTNrandTNTotal_average);
% %     title("Day2 after");
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% %     
% %     d_content = char('dFcorrMax_');
% %     figure;bar(x,R_D2after_TPTNrandTNTotal_Max);
% %     title("Day2 after");
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% % 
% %     d_content = char('dFcorrMin_');
% %     figure;bar(x,R_D2after_TPTNrandTNTotal_Min);
% %     title("Day2 after");
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% % 
% %     % figure (correlation), example
% %     cmax = 0.10;cmin = -0.10;
% % 
% %     session = char('D2after_');
% %     celltype = char('TNrand_');
% %     d_content = char('dFcorrAv_');
% %     
% %     figure;imagesc(R_D2after_TNrand)
% %     colorbar;caxis([cmin cmax]);
% %     title("Day2 after")  
% %     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% % 
% %     d_content = char('dFcorrAv_sub_');
% %     figure;imagesc(R_D2after_TNrand_sub)
% %     colorbar;caxis([cmin cmax]);
% %     title("Day2 after")  
% %     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% % 
% % 
% %     %%%%%%%%%%%%%%%%　別のセクションから　%%%%%%%%%%%%%%%%%%%%%
% % 
% % figure;imagesc(R_D2after(1:TP_last,1:TP_last))
% %     colorbar;caxis([cmin cmax]);
% %     title("Day2 after")  
% %     saveas(gcf, [mouseID,'_D2after_TP_CorrCell',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,'_D2after_TP_CorrCell',dtime], 'tif');
% % 
% %     figure;imagesc(R_D2after_sub(1:TP_last,1:TP_last))
% %     colorbar;caxis([cmin cmax]);
% %     title("Day2 after")  
% %     saveas(gcf, [mouseID,'_D2after_TP_CorrCell_sub',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
% %     saveas(gcf, [mouseID,'_D2after_TP_CorrCell_sub',dtime], 'tif');
% 
% %%%%% Day2 post %%%%%
% session = char('D2post_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(List_R_D2post_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMax_');
% writematrix(List_R_D2post_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% d_content = char('dFcorrMin_');
% writematrix(List_R_D2post_TNrand_sub_average,[mouseID,dinfo,session,celltype,d_content,'summary',num2str(i),dtime,'.txt'], 'delimiter','\t');%
% 
% Mean_List_R_D2post_TNrand_sub_average = mean(List_R_D2post_TNrand_sub_average);
% Mean_List_R_D2post_TNrand_sub_Max = mean(List_R_D2post_TNrand_sub_Max);
% Mean_List_R_D2post_TNrand_sub_Min = mean(List_R_D2post_TNrand_sub_Min);
% 
% Cell = [{'Tom+'};{'Tom- random'};{'Tom-'};{'All'}];
% R_D2post_TPTNrandTNTotal_average = [R_D2post_sub_TP_average;Mean_List_R_D2post_TNrand_sub_average;R_D2post_sub_TN_average;R_D2post_sub_All_average];
% R_D2post_TPTNrandTNTotal_Max = [R_D2post_sub_TP_Max;Mean_List_R_D2post_TNrand_sub_Max;R_D2post_sub_TN_Max;R_D2post_sub_All_Max];
% R_D2post_TPTNrandTNTotal_Min = [R_D2post_sub_TP_Min;Mean_List_R_D2post_TNrand_sub_Min;R_D2post_sub_TN_Min;R_D2post_sub_All_Min];
% 
% session = char('D2post_');
% celltype = char('TNrand_');
% d_content = char('dFcorrAv_');
% writematrix(R_D2post_TPTNrandTNTotal_average,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMax_');
% writematrix(R_D2post_TPTNrandTNTotal_Max,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% d_content = char('dFcorrMin_');
% writematrix(R_D2post_TPTNrandTNTotal_Min,[mouseID,dinfo,session,celltype,d_content,dtime,'.txt'], 'delimiter','\t');
% 
%     x =  categorical({'Tom+','Tom- random','Tom- ','All'});
%     x = reordercats(x,{'Tom+','Tom- random','Tom- ','All'});%これを入れないと、アルファベット順になってしまう
%     d_content = char('dFcorrAv_');
%     figure;bar(x,R_D2post_TPTNrandTNTotal_average);
%     title("Day2 post");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
%     
%     d_content = char('dFcorrMax_');
%     figure;bar(x,R_D2post_TPTNrandTNTotal_Max);
%     title("Day2 post");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     d_content = char('dFcorrMin_');
%     figure;bar(x,R_D2post_TPTNrandTNTotal_Min);
%     title("Day2 post");
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,dtime], 'tif');
% 
%     % figure (correlation), example
%     cmax = 0.10;cmin = -0.10;
% 
%     session = char('D2post_');
%     celltype = char('TNrand_');
%     d_content = char('dFcorrAv_');
%     
%     figure;imagesc(R_D2post_TNrand)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 post")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
%     d_content = char('dFcorrAv_sub_');
%     figure;imagesc(R_D2post_TNrand_sub)
%     colorbar;caxis([cmin cmax]);
%     title("Day2 post")  
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,session,celltype,d_content,num2str(i),dtime], 'tif');
% 
% 
%     %%%%%%%%%%%%%%%%　別のセクションから　%%%%%%%%%%%%%%%%%%%%%
% 
% figure;imagesc(R_D2post(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day2 post")  
%     saveas(gcf, [mouseID,'_D2post_TP_CorrCell',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D2post_TP_CorrCell',dtime], 'tif');
% 
%     figure;imagesc(R_D2post_sub(1:TP_last,1:TP_last))
%     colorbar;caxis([cmin cmax]);
%     title("Day2 post")  
%     saveas(gcf, [mouseID,'_D2post_TP_CorrCell_sub',dtime], 'fig');%https://jp.mathworks.com/help/matlab/matlab_prog/convert-between-datetime-arrays-numbers-and-strings.html
%     saveas(gcf, [mouseID,'_D2post_TP_CorrCell_sub',dtime], 'tif');
% 
% 
% cd ..;%上のフォルダに移動
%    
