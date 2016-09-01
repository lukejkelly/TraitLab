% Check reconstruction of Xi
% RJR 21/09/2009
% Assumes file is correctly formatted

xifile='synthdataMissingInBlocks3OUTXI.txt';

%leafnumbers=[2   	1   	3   	14   	7   	8   	9   	15   	4   	18   	19   	20   	13   	10   	16   	17   	12   	11   	6   	5 ];
%leafnumbers=[20   	1   	11   	12   	2   	4   	17   	18   	3   	5   	13   	14   	19   	15   	16   	8   	10   	9   	7   	6];
leafnumbers=[4   	6   	5   	16   	15   	1   	13   	14   	3   	18   	2   	17   	12   	11   	10   	8   	9   	7   	19   	20];

fid=fopen(xifile);
allstats = textscan(fid,'%f','headerlines',6);
fclose(fid);

fid=fopen(xifile);
xx=textscan(fid,'%s%[^\n]','delimiter','\n','headerlines',6);
fclose(fid);

cols = length(str2num(xx{1}{1}));

allstats=allstats{1};

Nsamp=length(allstats)/cols;
allstats = reshape(allstats,cols,Nsamp);

%allxi=allstats;

for i=1:(cols-1)
    allxi(leafnumbers(i),:)=allstats(i+1,:);
end


%%
% Draw
%truexi=[0.4232    0.9890    0.4935   0.9983   0.7352    0.9968   0.8554   0.8974   0.2671    0.6470   0.8406    0.7045   0.9861     0.6965    0.1101    0.9994  0.9667    0.9996    0.9030 0.5681];
%truexi=[    0.9827     0.9795    0.9751    0.9998    0.8031    0.9931    0.9978    0.6650    0.9732    0.2403    0.9928    0.3722    0.9995    0.9664    0.9990   0.6168    0.2796    0.4442    0.8121    0.4368];
truexi=[0.4670    0.7869    0.1287    0.1111    0.2743    0.9631    0.9110    0.5064    0.2306    0.2333    0.1880    0.4904    0.9463    0.8119    0.1040   0.6639    0.9754    0.9947    0.9999    0.9999];

good=zeros(1,20);

[sortedxi, index]=sort(truexi);

figure
hold on

for i=1:(cols-1)
    [low high]=hpd(allxi(index(i),100:end),5);
    plot([i+.3 i+.3],[low high],'LineWidth',2)
    plot(i,sortedxi(i),'xr','LineWidth',2)
    if sortedxi(i)>=low && sortedxi(i)<=high, good(i)=1; end
end
set(gca,'xtick',[])

XLIM([0 21]);
