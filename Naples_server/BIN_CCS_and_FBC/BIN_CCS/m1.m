icount=1;

% setfile1=['lcsf_temp/output_ref.dat'];
% setfile2=['lcsf_temp/lcfs_plot.dat']; 
 setfile1=['output_ref.dat'];
 setfile2=['lcfs_plot.dat']; 


A=load(setfile1);
C=load(setfile2);

[m,n]=size(C);

 
param1=['Time=',num2str(A(icount,1)+20),'sec'];

hCCS = plot(C(2*icount-1,2:A(icount,4)+1),C(2*icount,2:A(icount,4)+1), C(2*icount-1,A(icount,4)+2:A(icount,4)+A(icount,5)),C(2*icount,A(icount,4)+2:A(icount,4)+A(icount,5)))
xlim([1.5 5.0])
ylim([-3 3])
hold on
set(hCCS, 'LineWidth', 0.5)
set(hCCS, 'Color', 'r')
set(hCCS, 'LineStyle', '--')

%  xlabel({'R[m]';' ';param1})
xlabel('R[m]')
ylabel('Z[m]')



