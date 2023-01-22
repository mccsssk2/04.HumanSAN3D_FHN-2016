% Make the 2 cell plots. x axis = GJC (col 2), y axis = q_value (col 3)
% std0 = col 5, std1 = col10
close all;
close all;
clear all;
clear all;
data     = load('total_data.data');
gjc      = data(:,2) ;
q_value  = data(:,3) ;
mean0    = data(:,4) ;
std0     = data(:,5) ;
num_cls0 = data(:,6) ;
mean1    = data(:,9) ;
std1     = data(:,10);
num_cls1 = data(:,11);
%
%
% remove parts that are oscillation period 1, so that the Arnold's tongue shows clearly.
% CL0
j = 1;
for i=1:1:length(gjc)
	if num_cls0(i) > 1 % condition that num_cls0 does not have period 1, which is the data you want.
		gjc_ap0(j)      = gjc(i)     ;
		q_value_ap0(j)  = q_value(i) ;
		num_cls0_ap(j) = num_cls0(i);
		j = j + 1;
	end;
end;
figure;
h1 = scatter3(gjc_ap0, q_value_ap0,num_cls0_ap, 5, num_cls0_ap, 'filled');
view(0,90);
grid off;
xlabel('GJC','FontSize',20);
ylabel('q','FontSize',20);
title('distinct cycle lengths, Cell #1','FontSize',12);
colorbar;
axis([0 0.015 0 2000]);
saveas(h1,'cell1gjcq_cl_panel1.png','png');
%
% now draw the CL1.
j = 1;
for i=1:1:length(gjc)
	if num_cls1(i) > 1 % condition that num_cls0 does not have period 1, which is the data you want.
		gjc_ap1(j)      = gjc(i)     ;
		q_value_ap1(j)  = q_value(i) ;
		num_cls1_ap(j) = num_cls1(i);
		j = j + 1;
	end;
end;
figure;
h2 = scatter3(gjc_ap1, q_value_ap1, num_cls1_ap, 5, num_cls1_ap, 'filled');
view(0,90);
grid off;
xlabel('GJC','FontSize',20);
ylabel('q','FontSize',20);
title('distinct cycle lengths, Cell #2','FontSize',12);
colorbar;
axis([0 0.015 0 2000]);
saveas(h2,'cell2gjcq_cl_panel2.png','png');
%
% measure the synchronisation by taking the difference.

clear gjc_ap0, q_value_ap0;
j = 1;
for i=1:1:length(gjc)
	if abs(num_cls0(i) - num_cls1(i)) > 0 % condition that num_cls0 does not have period 1, which is the data you want.
		gjc_ap0(j)      = gjc(i)     ;
		q_value_ap0(j)  = q_value(i) ;
		synced(j)       = abs(num_cls0(i) - num_cls1(i));
		j = j + 1;
	end;
end;
figure;
h3 = scatter3(gjc_ap0, q_value_ap0, synced, 5, synced, 'filled');
view(0,90);
grid off;
xlabel('GJC','FontSize',20);
ylabel('q','FontSize',20);
title('CL Cell#1 - CL Cell#2','FontSize',12);
colorbar;
axis([0 0.015 0 2000]);
saveas(h3,'synced_or_not_panel3.png','png');
