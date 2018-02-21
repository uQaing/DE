% version 1.1.1
% NI the number of Iteration
% popsize the number of particles
% D the number of dimension
% F
% CR
% fit the function of fitness
function [fitBests,nItersions]=SSDE19(Nt,Ni,Npp,Nd,Kk,Ww,upLow,fit)


disp(fit);
Fmin=0.1;
Fmax=0.8;
CRmin=0.3;
CRmax=0.9;
PBest=0.9;
BBest=0.5;
Np=Npp-1;
xmin= -upLow;
xmax = upLow;
q=0.15;
fsn=1.5;
fitBests=zeros(Nt,Ni);
nItersions=zeros(Nt,1);

for t=1:Nt
    Fm=0.5*ones(Np,1);
    CRm=0.6*ones(Np,1);
    plan1=0;
    plan2=0;    
    XG = xmin+(xmax-xmin)*rand(Np,Nd);  %产生Np个D维向量
    XG_V=zeros(Np,Nd);
    % Fs and CRs in first iterasion.
    Fs=Fmin+(Fmax-Fmin)*rand(Np,1);
    CRs=CRmin+(CRmax-CRmin)*rand(Np,Nd);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ps=PBest+0.1*rand(Np,1);
    Bs=BBest+0.1*rand(Np,1);
   xd_Min=min(XG,[],1);
   xd_Max=max(XG,[],1);
    Fs_success=Fs;

     for G=1:Ni

        XG_n=xd_Min+(xd_Max-xd_Min).*rand(Np,Nd);
        tempN=fit(XG_n)<fit(XG);
        tempN=repmat(tempN,1,Nd);
        XG=XG_n.*tempN-XG.*(tempN-1);
        
        temp_fit=fit(XG);
        [~,temp_fit_index]=sort(temp_fit,1);% 排序
        xp=ceil(0.5*Np*((1-(G-1)/Ni)).^1.5);
        topIndex=temp_fit_index(1:xp); 
        baseVectors=topIndex(randi(xp,Np,1));
        
        QXG=XG(randperm(Np,q*Npp),:);
        [~,i]=min(fit(QXG));
        %% %%%%%%%%%%%%%%%%%%%%%%----变异操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        for n=1:Np
        R3=randperm(Np,5);
        
            if rand <Ps(n)
                if rand <Bs(n)
%                     son=XG(indexBest,:)+ Fs(n).*(XG(R3(2),:) - XG(R3(1),:));
                      son = XG(R3(3),:) + Fs(n)*(QXG(i,:)-XG(R3(3),:) +XG(R3(2),:) - XG(R3(1),:));
                else
                    son = XG(baseVectors(n),:) + Fs(n)*(XG(R3(2),:) - XG(R3(1),:));
%                       son = XG(baseVector,:) + Fs(n)*(XG(R3(2),:) - XG(R3(1),:));
                end
            else
                    son = XG(R3(3),:) + (0.1+0.7*rand)*(XG(R3(2),:) - XG(R3(1),:));
                    tempCRs=normrnd(0.5,1);
                    while((tempCRs<0)||(tempCRs>1))
                        tempCRs=normrnd(0.5,0.1);
                    end
                    CRs(n,:)= tempCRs;
            end           
            XG_V(n,:)=son;
        end
          XG_V=setWithInAre(XG_V,[xmin,xmax]);
        %% %%%%%%%%%%%%%%%%%%%%%---交叉操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        randR=randi([1,Nd],Np,1);
        dV=repmat(1:Nd,Np,1);
        
        tt=(rand(Np,Nd)<CRs);
        ttt= (dV==randR);
        tempX= tt|ttt;
        XG_C=XG_V.*tempX-XG.*(tempX-1);
        %% %%%%%%%%%%%%%%%%----选择操作---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tempC_i=fit(XG_C)<fit(XG);
        tempC=repmat(tempC_i,1,Nd);
        XG_next=XG_C.*tempC-XG.*(tempC-1);
        
        
        %%
        XG = XG_next;
        fits=fit(XG);
        %从中指找出最小的fitness 和 particle的indox
        [valueBest,indexBest] = min(fits);
        %保存本次迭代最小fitness
        fitBests(t,G) = valueBest;
        nItersions(t)=G;
        if valueBest==0
            break;
        end
        
        
        %% ready next itermsion

        %%create next P
        usable=(sum(tempC_i)/Np);
        
     K=0.05;
        if  usable>K
            plan1=plan1+1;     
            Fs_success=Fs.*tempC_i;
            
            
            temp=tempC_i.*Ps;
            u=sum(temp);
            d=sum(temp.^2);
            Ps=tempC_i.*Ps+(1-tempC_i).*(d/u);
            
            temp=tempC_i.*Bs;
            u=sum(temp);
            d=sum(temp.^2);
            Bs=tempC_i.*Bs+(1-tempC_i).*(d/u);
               
            temp=tempC_i.*Fs;
            Fsuccess=temp;
%         wf=0.8+0.2*rand(Np,1);
        Fsuccess(temp==0)=[];
        M=Fsuccess.^(fsn);
        N=length(Fsuccess);
        sumM=sum(M);
        MiN=(sumM/N)^(1/fsn);
        meanPowFsuccess=MiN;
        
        Fm=tempC_i.*Fm+(1-tempC_i)*(meanPowFsuccess);
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
        
        temp=tempC_i.*CRs;
        CRsuccess=temp;
%         wCR=0.9+0.1*rand(Np,1);
        CRsuccess(temp==0)=[];
        M=CRsuccess.^(fsn);
        N=length(CRsuccess);
        sumM=sum(M);
        MiN=(sumM/N)^(1/fsn);
        meanPowCRsuccess=MiN;
        
        CRm=tempC_i.*CRm+(1-tempC_i)*(meanPowCRsuccess);
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
                       
        else
            plan2=plan2+1;
            
            Ps=PBest+0.1*rand(Np,1);
            Bs=BBest+0.1*rand(Np,1);
            Fs=normrnd(mean(Fs_success),0.1,Np,1);
            
            
            avgFit=mean(fits);
            minFit=min(fits);
            maxFit=max(fits);
            %
            avgFits=repmat(avgFit,Np,1);
            minFits=repmat(minFit,Np,1);
            maxFits=repmat(maxFit,Np,1);
            
            W=1;
            J=100;
            t3=((fits-avgFits)./(maxFits-minFits));
            t3=t3.^W;
            t3=log(J*t3+1)/log(J+1);
            tempCR=(fits<=avgFits);
            CRs=CRs.*tempCR-(CRmin+(CRmax-CRmin)*t3).*(1-tempCR);
            
            Fs=setWithInAre(Fs,[Fmin,Fmax]);
            CRs=setWithInAre(CRs,[CRmin,CRmax]);
        end
            XGG=XG;
            XG_end_fromMean=mean(XG);
            if fit(XG_end_fromMean)<fit(XG(Np,:))
                XGG(Npp,:)=XG_end_fromMean;
            end
            xd_Min=min(XGG,[],1);
            xd_Max=max(XGG,[],1); 
        
     end
    disp(['SSDE第',num2str(t),'次测试 在第',num2str(G),'次迭代找到最优值：',num2str(fitBests(t,G)),' plan1= ',num2str(plan1),'plan2= ',num2str(plan2)]);
end

end
function y=setWithInAre(x,a)
c=(x>=a(1))&(x<=a(2));
y=x.*c+(a(1)+(a(2)-a(1))*rand(size(x))).*(1-c);
end