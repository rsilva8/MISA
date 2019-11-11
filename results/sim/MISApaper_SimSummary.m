% pwd
% ans =
% \\loki\export\mialab\users\rsilva\projects\MultivariateICA\MISA\fresh\output\S_case12_A_fake\r002\tuln_noinit_
cd ..\..\r002\con_bound_\
inlist = dir('*.mat');
out = [];
for rr = 1:length(inlist)
    load(inlist(rr).name,'mISI')
    out = [out mISI];
end
median(out)
cd ..\tuln_\
inlist = dir('*.mat');
out = [];
for rr = 1:length(inlist)
    load(inlist(rr).name,'mISI')
    out = [out mISI];
end
median(out)
cd ..\tuln_noinit_\
inlist = dir('*.mat');
out = [];
for rr = 1:length(inlist)
    load(inlist(rr).name,'mISI')
    out = [out mISI];
end
median(out)


% pwd
% ans =
% \\loki\export\mialab\users\rsilva\projects\MultivariateICA\MISA\fresh\output\S_case14_A_fake\r002\qle_
cd ../../r002/con_bound_/
inlist = dir('*.mat');
out = [];
for rr = 1:length(inlist)
load(inlist(rr).name,'mISI')
out = [out mISI];
end
median(out)
cd ../../r002/dlahat_/
inlist = dir('*.mat');
out = [];
for rr = 1:length(inlist)
load(inlist(rr).name,'mISI')
out = [out mISI];
end
median(out)
cd ../../r002/qle_/
inlist = dir('*.mat');
out = [];
for rr = 1:length(inlist)
load(inlist(rr).name,'mISI')
out = [out mISI];
end
median(out)