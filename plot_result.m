close all;
%{
err_sample = [0.017250722	0.019648144	0.018000791;
0.00385007	0.005094558	0.006140878;
0.002225255	0.002776267	0.003445298;
0.001537398	0.001600408	0.002493266;
0.001208245	0.001446583	0.001940691];
%}
err_sample = [0.408820504 0.494650474 0.425018655;
0.16393496 0.201438674 0.257363612;
0.124376857 0.154421389 0.195215316;
0.09956037 0.114968409 0.162751583;
0.089689923 0.092422929 0.152381577];
stdev_sample = [0.133738946 0.149578183 0.115595246;
0.048538731 0.076259572 0.073876311;
0.023935037 0.094124045 0.06032811;
0.019430989 0.05857937 0.03871613;
0.016362442 0.017823443 0.047129658];














err_sample_x = [100
300
500
700
900];


err_source = [0.00120767	0.001308857	0.001956101;
0.002284646	0.002812135	0.003585946;
0.002690336	0.003630127	0.004036694;
0.003481484	0.004180403	0.004709381;
0.004568755	0.005506949	0.005200658];

err_source_x = [2
6
10
14
18];

err_g_noise = [0.002387645	0.002779743	0.00373534;
0.003268225	0.003905077	0.004265354;
0.008942933	0.009433814	0.009139924
];

err_g_noise_x = [1.00E-03
1.00E-02
1.00E-01];

err_u_noise = [0.002160402	0.002676054	0.003551076;
0.002562269	0.00304364	0.003921631;
0.004933391	0.005405321	0.005660125;
];

err_u_noise_x = [2.00E-03
2.00E-02
2.00E-01];

err_g_sample = [0.002623379	0.003230932	0.003928479
0.002997336	0.003848104	0.004416748
0.003644267	0.005387648	0.004924151
0.006394441	0.007781368	0.007163664
0.006627499	0.008414133	0.006997259];

err_g_sample_x = [0.1
    0.2
    0.3
    0.4
    0.5];

global_rank = [0.141347656	0.083007813;
0.050933228	0.025604248;
0.023641129	0.011960983;
0.016957039	0.009182513];

global_rank_x = [10 15 20 25];


figure('Color',[1 1 1]);
plot(err_sample_x,err_sample(:,1),'k-o',...
err_sample_x,err_sample(:,2),'k--s',...
err_sample_x,err_sample(:,3),'k-.x','LineWidth',1.2);

errorbar(err_sample_x,err_sample(:,1),stdev_sample(:,1),'rx');
hold on;
errorbar(err_sample_x,err_sample(:,2),stdev_sample(:,2),'bo');
hold on;
errorbar(err_sample_x,err_sample(:,3),stdev_sample(:,3),'g*');
hold on;

legend('Lp-ICA-F','Lp-ICA-G','FastICA');
set(gca, 'XTick',  [100 300 500 700 900]);
axis([0 1000 0 0.025]);
xlabel('Number of samples');
ylabel('MSE');

figure('Color',[1 1 1]);
plot(err_source_x,err_source(:,1),'k-o',...
err_source_x,err_source(:,2),'k--s',...
err_source_x,err_source(:,3),'k-.x','LineWidth',1.2);
legend('Lp-ICA-F','Lp-ICA-G','FastICA');
set(gca, 'XTick',  [2 6 10 14 18]);
set(gca, 'YTick',  [0.002 0.004 0.006 0.008 0.01]);
set(gca, 'YTickLabel',  {'0.002', '0.004', '0.006', '0.008', '0.01'});
axis([0 20 0 0.006]);
xlabel('Number of sources');
ylabel('MSE');

figure('Color',[1 1 1]);
plot(err_g_noise_x,err_g_noise(:,1),'k-o',...
err_g_noise_x,err_g_noise(:,2),'k--s',...
err_g_noise_x,err_g_noise(:,3),'k-.x','LineWidth',1.2);
legend('Lp-ICA-F','Lp-ICA-G','FastICA');
set(gca, 'Xscale',  'log');
set(gca, 'XTick',  [1e-3 1e-2 1e-1]);
set(gca, 'YTick',  [0.002 0.004 0.006 0.008 0.01]);
axis([5e-4 0.3 0 0.01]);
xlabel('Standard deviation of Gaussian noise');
ylabel('MSE');

figure('Color',[1 1 1]);
plot(err_u_noise_x,err_u_noise(:,1),'k-o',...
err_u_noise_x,err_u_noise(:,2),'k--s',...
err_u_noise_x,err_u_noise(:,3),'k-.x','LineWidth',1.2);
legend('Lp-ICA-F','Lp-ICA-G','FastICA');
set(gca, 'Xscale',  'log');
set(gca, 'XTick',  [2e-3 2e-2 2e-1]);
set(gca, 'YTick',  [0.002 0.004 0.006 0.008 0.01]);
axis([1e-3 0.3 0 0.01]);
xlabel('Interval of uniform noise');
ylabel('MSE');

figure('Color',[1 1 1]);
plot(err_g_sample_x,err_g_sample(:,1),'k-o',...
err_g_sample_x,err_g_sample(:,2),'k--s',...
err_g_sample_x,err_g_sample(:,3),'k-.x','LineWidth',1.2);
legend('Lp-ICA-F','Lp-ICA-G','FastICA');
set(gca, 'YTick',  [0.002 0.004 0.006 0.008 0.01]);
axis([0 0.6 0 0.01]);
xlabel('Ratio of samples from Gaussian distribution');
ylabel('MSE');

figure('Color',[1 1 1]);
plot(global_rank_x,global_rank(:,1),'k-o',...
global_rank_x,global_rank(:,2),'k--s','LineWidth',1.2);
legend('mean rank','median rank');
set(gca, 'XTick',  [10 15 20 25]);
axis([5 30 0 0.16]);
xlabel('Number of samples');
ylabel('Rank ratio');