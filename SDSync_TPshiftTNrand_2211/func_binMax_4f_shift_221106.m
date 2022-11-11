function binning_slide = func_binMax_4f_shift_221106(X)

[n_X, m_X] = size(X);
%n_X: total time frames
%m_X: total number of neurons

bin_width = 4;%frames
slide = bin_width/2;%2frame shift
X_bin = zeros(floor(n_X/slide)-1,m_X);
for i = 1:floor(n_X/slide)-1
    X_bin(i,:) = max(X(slide*i-slide+1:slide*i+slide,:));
end

binning_slide = X_bin;
end
