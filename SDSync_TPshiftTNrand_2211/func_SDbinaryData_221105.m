function SDbinaryData = func_SDbinaryData_221105(X,SD_thres);
%raw df data matrixから、各細胞ごとにmean, sdを算出し、z-sscoreに変換、閾値でbinarize
dfX_mean = mean(X);%各細胞のmean dF, https://jp.mathworks.com/help/matlab/ref/mean.html?searchHighlight=mean&s_tid=srchtitle_mean_1
dfX_sd= std(X);%各細胞のSD dF, https://jp.mathworks.com/help/matlab/ref/std.html?searchHighlight=%E6%A8%99%E6%BA%96%E5%81%8F%E5%B7%AE&s_tid=srchtitle_%25E6%25A8%2599%25E6%25BA%2596%25E5%2581%258F%25E5%25B7%25AE_1
sdDataX = (X - dfX_mean)./dfX_sd;%各データの細胞ごとのz-score
SDbinaryData = sdDataX>SD_thres;%thresholdで２値化
