function LD = LoDFun(V)

% Data from MSL paper, 0.24 roughly 
v = [ 6000
    4992.512
3996.452
2519.6536
1890.0615
1352.2175
1203.1013
517.1026
0];

LoD = [0.2375
    0.23756167
0.24493533
0.26297614
0.27370104
0.29312414
0.27746123
0.33886394
0.33886];

scale = 0.29/0.24;
scale = 1; % for MSL-class 

LD = interp1(v, LoD*scale, V, 'linear');