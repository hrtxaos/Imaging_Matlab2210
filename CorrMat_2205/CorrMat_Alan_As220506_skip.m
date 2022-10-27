%% Do not delelte this script, for PVD calculation
%correlationとそのp-valueを得られるように改良しようとしたが、correlation
%coefficincyの計算がcorr/corrcoefではなかったので諦めて、この値でshuffleやcontrolと比較することにした
%correlatoinの値の計算を少しだけ改変
function CorrMat_Alan_As220506_skip = Alan_cor_matching_new(X, shift_win,shift)
  %X:binned data matrix(time/neurons)
  %shift_win: time window (default = 5 binned frames > shift = 4でないと5frameにならないだろう)
  %X: raw dF matrix, shift = 4(frame step)
  [n_X, m_X] = size(X);
  % n_X:total time frames
  % m_X: total number of neurons
  combination = factorial(m_X)/factorial(2)/factorial(m_X-2);%Asai_220502
  if n_X <= shift_win
    error('The shift is too large.')
  end
  % Create matrix for substituting the correlation coefficency
%   CorrMat_Alan_As220506 = zeros(n_X-shift_win, n_X-shift_win);
  CorrMat_Alan_As220506_skip = zeros(floor((n_X-shift_win)/shift), floor((n_X-shift_win)/shift));
%   (切り捨て)floor;%https://jp.mathworks.com/help/matlab/ref/floor.html

    %repeat from 1st frame to last
  for i = 1:floor(n_X/shift)
    
    if (4*i+shift_win) > n_X
      break
    end
  
    for j = 1:floor(n_X/shift)
      
      if (4*j+shift_win) > n_X
        break
      end
      
      tmp_1 = corr(X(4*i-3:(4*i-3+shift_win), :)) - eye(m_X);%time "i"での神経活動の同期性を総当たりでcorrelation matrixに (eyeで対角成分の1を0に。これ尾wしないと、最後にsumを取って、correlationとするときに、negative correlationがかき消されるのだろう)
      tmp_2 = corr(X(4*j-3:(4*j-3+shift_win), :)) - eye(m_X);%time "j"での神経活動の同期性を総当たりでcorrelation matrixに
      %between time"i" and "j"
%       CorrMat_Alan_220502(i,j) = sum(sum(tmp_1 .* tmp_2)) / (m_X^2);%Alan's original correlationこれは、同じ細胞由来の成分を2回足しているのに対し、細胞数の2乗で割っているので、細胞数が多いほど厳しくなりすぎている気がする。
%       CorrMat_Alan_220502(i,j) = sum(sum(tmp_1 .* tmp_2)) / (m_X*2);%これだと細胞数が多くなると、組み合わせが多くなるので、correlationが1を超えてくるはず。
      CorrMat_Alan_As220506_skip(i,j) = sum(sum(tmp_1 .* tmp_2)) / (combination*2);%Asai_220502, 組み合わせの数で割る、つまり平均を取る感じ。
%       j = j + shift;%As_220506, xx 
    end
%       i = i + shift;%As_220506, xx
  end
end
