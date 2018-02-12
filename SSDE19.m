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

fitBests=zeros(Nt,Ni);
nItersions=zeros(Nt,1);
R3s=zeros(Np,3);
for t=1:Nt



    plan1=0;
    plan2=0;
    
    XG = xmin+(xmax-xmin)*rand(Np,Nd);  %产生Np个D维向量
    XG_V=zeros(Np,Nd);
%     [~,indexBest] = min(fit(XG));
    
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
        xp=ceil(0.5*Np*(1-(G-1)/Ni));
        topIndex=temp_fit_index(1:xp); 
%         baseVector=topIndex(randperm(xp,1));
        baseVectors=topIndex(randi(xp,Np,1));
        
        QXG=XG(randperm(Np,q*Npp),:);
        [~,i]=min(fit(QXG));
        %% %%%%%%%%%%%%%%%%%%%%%%----变异操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        for n=1:Np
%         baseVector=topIndex(randperm(xp,1));
        
%         QXG=XG(randperm(Np,q*Npp),:);
%         [~,i]=min(fit(QXG));
        
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
%                 son=xmin+(xmax-xmin)*rand(1,Nd);
                
                    son = XG(R3(3),:) + (0.1+0.7*rand)*(XG(R3(2),:) - XG(R3(1),:));
                    %                       CRs(n,:)=CRmin+(CRmax-CRmin)*rand(1,Nd);
                    tempCRs=normrnd(0.5,1);
                    while((tempCRs<0)||(tempCRs>1))
                        tempCRs=normrnd(0.5,0.1);
                    end
                    CRs(n,:)= tempCRs;


%                       son = XG(baseVector,:) + Fs(n)*(XG(R3(2),:) - XG(R3(1),:));
%              son=  XG(n,:)+ (0.1+0.7*rand)*(XG(R3(1),:) - XG(n,:))+Fs(n)*(XG(R3(2),:) - XG(R3(3),:));
%              son = XG(R3(1),:) + Fs(n)*(XG(R3(2),:) - XG(R3(3),:))+Fs(n)*(XG(R3(4),:) - XG(R3(5),:));
            end           
            XG_V(n,:)=son;
        end



%        for n=1:Np
%        R3s(n,:)=randperm(Np,3);
%        end
%        sonPTB = XG(R3s(3),:) + Fs.*(QXG(i,:)-XG(R3s(3),:)+XG(R3s(:,2),:) - XG(R3s(:,1),:));
%        sonPFB = XG(baseVector,:) + Fs.*(XG(R3s(2),:) - XG(R3s(1),:));
%        sonFP = XG(R3s(3),:) + Fs.*(XG(R3s(2),:) - XG(R3s(1),:));
%        
%        randP=rand(Np,1);
%        randB=rand(Np,1);
%        logicP=randP<Ps;
%        logicB=randB<Bs;
%        
%        logicPTB=logicP&logicB;
%        logicPFB=logicP&(~logicB);
%        logicFP=(~logicP);
%        
%        templogic=logicPTB+logicPFB+logicFP;
%        
%        sonNextPTB=logicPTB.*sonPTB;
%        sonNextPFB=logicPFB.*sonPFB;
%        sonNextFP=logicFP.*sonFP;
%           
%         XG_V=sonNextPTB+sonNextPFB+sonNextFP;
       
%         XG_V(XG_V<xmin)=(xmax - xmin)*rand(1) + xmin;
%         XG_V(XG_V>xmax)=(xmax - xmin)*rand(1) + xmin;
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
            Ps_success=Ps.*tempC_i;
            Bs_success=Bs.*tempC_i;
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
            u=sum(temp);
            d=sum(temp.^2);
            Fs=tempC_i.*Fs+(1-tempC_i).*(d/u);
           
            temp=tempC_i.*CRs;
            u=sum(temp,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            d=sum(temp.^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            CRs=tempC_i.*CRs+(1-tempC_i).*(d./u);
                       
        else
            plan2=plan2+1;
            
            Ps=PBest+0.1*rand(Np,1);
            Bs=BBest+0.1*rand(Np,1);
%             Fs=gamrnd(mean(Fs_success),0.1,Np,1);
            Fs=normrnd(mean(Fs_success),0.1,Np,1);
            
            
            % Fs and CRs in first iterasion.
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
            
        end
             Ps(Ps>PBest+0.1)=PBest+0.1*rand;
             Ps(Ps<PBest)=PBest+0.1*rand;
             
             Bs(Ps>BBest+0.1)=BBest+0.1*rand;
             Bs(Ps<BBest)=BBest+0.1*rand;
             
             Fs(Fs>Fmax)=Fmin+(Fmax-Fmin)*rand;
             Fs(Fs<Fmin)=Fmin+(Fmax-Fmin)*rand;
             
            %
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