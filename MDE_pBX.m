% NI the number of Iteration
% popsize the number of particles
% D the number of dimension
% F
% CR
% fit the function of fitness
function [fitBests,nItersions]=MDE_pBX(Nt,Ni,Np,Nd,F,CR,range,fit)
%     flag=0;
%     valueBestOld=zeros(1,1);
tic;

disp(fit);
xmin= -range;
xmax = range;
R3s=zeros(Np,3);
fitBests=zeros(Nt,Ni);
nItersions=zeros(Nt,1);
Fmin=0.0;
Fmax=1;
CRmin=0.0;
CRmax=1;

q=0.15;
fsn=1.5;


for k=1:Nt
    XG = xmin+(xmax-xmin)*rand(Np,Nd);  %产生Np个D维向量
    XG_V=zeros(Np,Nd);
    Fm=0.5*ones(Np,1);
    CRm=0.6*ones(Np,1);
    Fs=cauchypdf(Fm,0.1);
    CRs=normrnd(CRm,0.1);
    Fs(Fs>Fmax)=Fmin+(Fmax-Fmin)*rand;
    Fs(Fs<Fmin)=Fmin+(Fmax-Fmin)*rand;
    CRs(CRs>CRmax)=CRmin+(CRmax-CRmin)*rand;
    CRs(CRs<CRmin)=CRmin+(CRmax-CRmin)*rand;
    for G=1:Ni
        temp_fit=fit(XG);
        [~,temp_fit_index]=sort(temp_fit,1);% 排序
        xp=ceil(0.5*Np*(1-(G-1)/Ni));
        topIndex=temp_fit_index(1:xp); 
        baseVectors=topIndex(randi(xp,Np,1));
        
%         QXG=XG(randperm(Np,q*Np),:);
%         [~,i]=min(fit(QXG));
        for n=1:Np

%         baseVector=topIndex(randperm(xp,1));
        
        QXG=XG(randperm(Np,q*Np),:);
        [~,i]=min(fit(QXG));
        R3s=randperm(Np,3);
        %% %%%%%%%%%%%%%%%%%%%%%%----变异操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
        son = XG(baseVectors(n),:) + Fs(n).*(QXG(i,:)-XG(baseVectors(n),:)+XG(R3s(:,2),:) - XG(R3s(:,1),:));
        XG_V(n,:)=son;
        end
        XG_V(XG_V<xmin)=(xmax - xmin)*rand(1) + xmin;
        XG_V(XG_V>xmax)=(xmax - xmin)*rand(1) + xmin;
        
        %% %%%%%%%%%%%%%%%%%%%%%---交叉操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             randR=randi([1,Nd],Np,1);
        %             dV=repmat(1:Nd,Np,1);
        %
        %             tempX=(rand(Np,Nd)<=CRs) | (dV==randR);
        %             XG_C=XG_V.*tempX-XG.*(tempX-1);
        
        randR=randi([1,Nd],Np,1);
        dV=repmat(1:Nd,Np,1);
        
        tt=(rand(Np,Nd)<CRs);
        ttt= (dV==randR);
        tempX= tt|ttt;
        XG_C=XG_V.*tempX-XG.*(tempX-1);
        
        %% %%%%%%%%%%%%%%%%----选择操作---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tempC=fit(XG_C)<fit(XG);
        tempSuccesslogic=tempC;
        tempC=repmat(tempC,1,Nd);
        XG_next=XG_C.*tempC-XG.*(tempC-1);
        
        
        %% 测试进化后的粒子
        XG = XG_next;
        %从中指找出最小的fitness 和 particle的indox
        [valueBest,indexBest] = min(fit(XG));
        %保存本次迭代最小fitness
        fitBests(k,G) = valueBest;
        nItersions(k)= G;
        
        %
        if valueBest==0
            break;
        end
        usable=(sum(tempSuccesslogic)/Np);
        
        if usable==0
            Fsuccess= Fsuccess_old;
            CRsuccess=CRsuccess_old;
        else
            Fsuccess=Fs.*tempSuccesslogic;
            CRsuccess=CRs.*tempSuccesslogic;
        end
        
        Fsuccess_old=Fsuccess;
        wf=0.8+0.2*rand(Np,1);
        Fsuccess(Fsuccess==0)=[];
        M=Fsuccess.^(fsn);
        N=length(Fsuccess);
        sumM=sum(M);
        MiN=(sumM/N)^(1/fsn);
        meanPowFsuccess=MiN;
        
        Fm=wf.*Fm+(1-wf)*(meanPowFsuccess);
        F2m=repmat(Fm,2,1);
        Fs =myCauchy(F2m,0.1);
        Fs=Fs(Fs>=0&Fs<=1);
        [r,~]=size(Fs);
        while(r<Np)
        Fs =myCauchy(F2m,0.1);
        Fs=Fs(Fs>=0&Fs<=1);
        [r,~]=size(Fs);
        end
        Fs=Fs(1:Np);
        
        
        
        CRsuccess_old=CRsuccess;
        wCR=0.9+0.1*rand(Np,1);
        
        CRsuccess(CRsuccess==0)=[];
        M=CRsuccess.^(fsn);
        N=length(CRsuccess);
        sumM=sum(M);
        MiN=(sumM/N)^(1/fsn);
        meanPowCRsuccess=MiN;
        
        CRm=wCR.*CRm+(1-wCR)*(meanPowCRsuccess);
        CR2m=repmat(Fm,2,1);
       
        CRs =normrnd(CR2m,0.1,2*Np,1);
        CRs=CRs(CRs>=0&CRs<=1);
        [r,~]=size(CRs);
        while(r<Np)
           CRs =normrnd(CR2m,0.1,2*Np,1);
        CRs=CRs(CRs>=0&CRs<=1);
        [r,~]=size(CRs);
        end
        CRs=CRs(1:Np);
          
    end
    disp(['MDE_pBX第',num2str(k),'次测试 在第',num2str(G),'次迭代找到最优值：',num2str(fitBests(k,G))]);
end

toc

end

function y=setWithInAre(x,a)
c=(x>=a(1))|(x<=a(2));
y=x.*c+(a(1)+(a(2)-a(1))*rand(size(x))).*(1-c);
end
