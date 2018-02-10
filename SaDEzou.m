%%%%%%%%%44444444444444444444444444444444       SADE    44444444444444444444444444444444444
function [v,i]=SaDEzou(Nt,Ni,popsize,Nd,F,CR,ranges,fit)
NI=3000;     %设置进化代数
NI4=NI;


F4=0.6;CR4=0.3;
N=38;
PD=6000;

xmin=[220   220   200   200   200   200   200   200   114   114   114   114   110    90    82   120    65    65    65   120   120   110    80    10    60    55    35    20    20    20    20    20    25    18     8    25    20    20];
xmax=[550   550   500   500   500   500   500   500   500   500   500   500   500   365   365   325   315   315   315   272   272   260   190   150   125   110    75    70    70    70    70    60    60    60    60    60    38    38];
c=[64782       64782       64670       64670       64670       64670       64670       64670      172832      172832      176003      173028       91340       63440       65468       77282      190928      285372      271676       39197       45576       28770       36902      105510       22233       30953       17044       81079      124767      121915      120780      104441       83224      111281       64142      103519       13547       13518];
b=[796.9 796.9 795.5 795.5 795.5 795.5 795.5 795.5 915.7 915.7 884.2 884.2 1250.1 1298.6 1298.6 1290.8 238.1 1149.5 1269.1 696.1 660.2 803.2 818.2 33.5 805.4 707.1 833.6 2188.7 1024.4 837.1 1305.2 716.6 1633.9 969.6 2625.8 1633.9 694.7 655.9];
a=[0.3133 0.3133 0.3127 0.3127 0.3127 0.3127 0.3127 0.3127 0.7075 0.7075 0.7515 0.7083 0.4211 0.5145 0.5691 0.5691 2.5881 3.8734 3.6842 0.4921 0.5728 0.3572 0.9415 52.123 1.1421 2.0275 3.0744 16.765 26.355 30.575 25.098 33.722 23.915 32.562 18.362 23.915 8.482 9.693];

% F=0.6;CR=0.3;
lamb=10^5;
NC=2;


best4=zeros(NC,N+3);fbest4=zeros(NC,NI);t4=zeros(1,NC);Tau1=0.1;Tau2=0.1;Fl=0.1;Fu=0.9;


tic
for nr=1:NC    
pop=zeros(popsize,N+3);pop1=zeros(popsize,N+3);
for i=1:popsize
    for j=1:N
        pop(i,j)=xmin(j)+rand()*(xmax(j)-xmin(j));
    end
    Psum=sum(pop(i,1:N));
        w=abs(Psum-PD);
        if(w~=0)
            if(Psum<PD)
                n=1:N;nc=size(n,2);
                j0=ceil(rand()*nc);jc=n(j0);
                while(pop(i,jc)==xmax(jc))
                    n(n==jc)=[];nc=size(n,2);
                    j0=ceil(rand()*nc);jc=n(j0);
                end
                pop(i,jc)=min(pop(i,jc)+w,xmax(jc));
            else
                n=1:N;nc=size(n,2);
                j0=ceil(rand()*nc);jc=n(j0);
                while(pop(i,jc)==xmin(jc))
                    n(n==jc)=[];nc=size(n,2);
                    j0=ceil(rand()*nc);jc=n(j0);
                end
                pop(i,jc)=max(pop(i,jc)-w,xmin(jc));
            end
        end
        pop(i,N+2)=0;
        for j=1:N
            pop(i,N+2)=pop(i,N+2)+a(j)*pop(i,j)^2+b(j)*pop(i,j)+c(j);
        end
        Psum=sum(pop(i,1:N));
        pop(i,N+3)=abs(Psum-PD);
        pop(i,N+1)=pop(i,N+2)+lamb*pop(i,N+3);
end
for ni=1:NI4
    %???????????????????
    if rand()<Tau1
        F4=Fl+rand()*Fu;
    end
    if rand()<Tau2
        CR4=rand();
    end
    pop1=pop;
    for i=1:popsize
        i1=ceil(rand()*popsize);
            while(i1==i)
                i1=ceil(rand()*popsize);
            end
            i2=ceil(rand()*popsize);
            while(i2==i||i2==i1)
                i2=ceil(rand()*popsize);
            end
            i3=ceil(rand()*popsize);
            while(i3==i||i3==i1||i3==i2)
                i3=ceil(rand()*popsize);
            end
        pop1(i,1:N)=pop(i1,1:N)+F4*(pop(i2,1:N)-pop(i3,1:N));%F后面乘以rand()效果更好
        j0=ceil(rand()*N);
        for j=1:N
            if(rand()<CR4||j==j0)
                pop1(i,j)=pop1(i,j);
            else
                pop1(i,j)=pop(i,j);
            end
            if pop1(i,j)>xmax(j)
                pop1(i,j)=xmax(j);
            elseif pop1(i,j)<xmin(j)
                pop1(i,j)=xmin(j);
            end
        end
        Psum=sum(pop1(i,1:N));
        w=abs(Psum-PD);
        if(w~=0)
            if(Psum<PD)
                n=1:N;nc=size(n,2);
                j0=ceil(rand()*nc);jc=n(j0);
                while(pop1(i,jc)==xmax(jc))
                    n(n==jc)=[];nc=size(n,2);
                    j0=ceil(rand()*nc);jc=n(j0);
                end
                pop1(i,jc)=min(pop1(i,jc)+w,xmax(jc));
            else
                n=1:N;nc=size(n,2);
                j0=ceil(rand()*nc);jc=n(j0);
                while(pop1(i,jc)==xmin(jc))
                    n(n==jc)=[];nc=size(n,2);
                    j0=ceil(rand()*nc);jc=n(j0);
                end
                pop1(i,jc)=max(pop1(i,jc)-w,xmin(jc));
            end
        end
        pop1(i,N+2)=0;
%         for j=1:N
%             pop1(i,N+2)=pop1(i,N+2)+a(j)*pop1(i,j)^2+b(j)*pop1(i,j)+c(j);
%         end
%         Psum=sum(pop1(i,1:N));
%         pop1(i,N+3)=abs(Psum-PD);
%         pop1(i,N+1)=pop1(i,N+2)+lamb*pop1(i,N+3);
          pop1(i,N+1)= sum(pop1(i,1:N).^2);
%         pop(i,N+1)=pop(i,N+2)+lamb*pop(i,N+3);
        if pop1(i,N+1)>pop(i,N+1)
            pop1(i,:)=pop(i,:);
        end
    end
    pop=pop1;
    fbest4(nr,ni)=min(pop(:,N+1));
end
ind1=find(pop(:,N+1)==min(pop(:,N+1)));ind=ind1(1);
best4(nr,:)=pop(ind,:);
t4(nr)=toc;
end
Averagefbest4=sum(fbest4)/NC;

fmin4=min(best4(:,N+1));
fmax4=max(best4(:,N+1));
fmean4=mean(best4(:,N+1));
fmedian4=median(best4(:,N+1));
fstd4=std(best4(:,N+1));
fdata4=[fmin4 fmax4 fmean4 fmedian4 fstd4];
% fdata4(1:4)=roundn(fdata4(1:4),-6)

x=1:NI;
x=x*popsize;
figure,semilogy(x,Averagefbest4,'b');%axis([1 NI*popsize 15274.930392 15276.930392]);
legend('SADE');
xlabel('Number of function evaluations');ylabel('Average function value');%title('Sphere function');
end