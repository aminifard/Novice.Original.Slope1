
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% driver_cs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear;
close all;
format long
global n n_S m_S ub1 A b Mn options 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'book2.xlsx';
sheet = 1;
xlRange = 'A4:C363';
St = xlsread(filename,sheet,xlRange);
[n,m]=size(St);
s=zeros(72,95);
for i=1:n
    s(St(i,1),St(i,2))=St(i,3);
end
S=s;
[m_S,n_S] = size(S);
A = [S zeros(m_S,n_S) zeros(m_S,n_S);eye(n_S) eye(n_S) zeros(n_S);-eye(n_S) zeros(n_S) eye(n_S)];
[m,n] = size(A);
Lt=xlsread(filename,'A410:B450');
[nl,ml]=size(Lt);
l1=zeros(n_S,1);
for i=1:nl
    l1(Lt(i,1))=Lt(i,2);
end
lb1=l1;
Ut=xlsread(filename,'A367:B406');
[nu,mu]=size(Ut);
u1=zeros(n_S,1);
for i=1:nu
    u1(Ut(i,1))=Ut(i,2);
end
ub1=u1;
b = [zeros(m_S,1);ub1;-lb1];
Mn = [lb1;ub1-lb1;zeros(n_S,1)];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mychat = xlsread('optimizer.data.xlsx');
% S = mychat(102:111,1:100);
% [m_S,n_S] = size(S);
% A = [S zeros(m_S,n_S) zeros(m_S,n_S);eye(n_S) eye(n_S) zeros(n_S);-eye(n_S) zeros(n_S) eye(n_S)];
% [m,n] = size(A);
% lb1 = mychat(1:100,1);
% ub1 = mychat(113:212,1);
% b = [zeros(m_S,1);ub1;-lb1];
% %b = zeros(m_S,1);
% Mn = [lb1;ub1-lb1;zeros(n_S,1)];
% length(b)
Slist    = {'fpcbbtest'}; % contain the list of solvers

%%%%%%%%%%%%%%%%%%  
% select options %
%%%%%%%%%%%%%%%%%%

mu             = 2^(-8); 
options.mu     = mu;
options.init   = 3; 
options.xtol   = 1E-10;
options.ftol   = 1E-10;
options.gtol   = 0.2;
options.mxitr  = 10000;
options.eta    = 4;
options.fullMu = false;
options.gamma  = 0.85; 
options.c      = 1E-3; 
options.beta   = 0.5;
options.scale  = true;
options.stopCriterion = 1;

%%%%%%%%%%%%%%%%%%%%
% generate problem %
%%%%%%%%%%%%%%%%%%%%

alpha    = 0.5; 
seed     = 0;
full     = true;
                                        
   for j = 1:length(Slist) 

       fprintf(['Running ',Slist{j},'...\n'])
                
        tic
        eval(['Output=',Slist{j},''])      
        toc
               
        eval(['itern.',Slist{j},'=Output.iterfinal;'])
        eval(['Time.',Slist{j},'=Output.cpu;'])
        eval(['funccount.',Slist{j},'=Output.nf;'])
        eval(['Nonzeros.',Slist{j},'=Output.nz;'])
        eval(['x.',Slist{j},'=Output.x;'])
       
              
   end
  d=0;
  %size(Output.x(i));
   for i=1:n_S
       if Output.x(i)<lb1(i) || Output.x(i)>ub1(i)
           i
          Output.x(i)
          lb1(i)
          ub1(i)
       end
   end
vout=x.fpcbbtest(1:95,1)
sout=sum(vout==0) 
filename = 'testdata.xlsx';
A = [vout];
xlswrite(filename,A)