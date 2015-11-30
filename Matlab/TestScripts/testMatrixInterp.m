clear; clc;
size = [20,20];
n = 100;
for i = 1:n
    M{i} = (i-1)*2/n*ones(size) + .01*rand(size);
end
m = cat(3,M{:});
t = linspace(0,1,n);
ti = linspace(0,1,10000);
% ti = rand();

tic
Mi = MatrixInterp(t,M,ti);
toc

% for i = 1:size(1)
%     for j = 1:size(2)
%         plot(t,squeeze(m(i,j,:)))
%         hold all
%     end
% end