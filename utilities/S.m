% l1-prox function = soft thresholding
function y = S(x, t)
mx = max(abs(x) - t, 0);
y = sign(x).*mx;
