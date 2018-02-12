% NI the number of Iteration 
% popsize the number of particles
% D the number of dimension
% F
% CR
% fit the function of fitness

function [fitBests,nItersions]=MDE(Nt,Ni,Npp,Nd,F,CR,upLow,fit)

disp(fit);
Np=Npp-1;
Fmin=0.5;
Fmax=1;
CRmin=0.8; 
CRmax=1;
Mf=0.75;
Qf=0.1;
Mcr=0.9;
Qcr=0.1;

xmin= -upLow;
xmax = upLow;
    
    fitBests=zeros(Nt,Ni);
    nItersions=zeros(Nt,1);
    for t=1:Nt
        G=0;
        nn=0;
         test0=0;
         test1=0;
         
         testBest=0;
         testRand=0;
         testRe=0;
 
        XG = xmin+(xmax-xmin)*rand(Np,Nd);  %����Np��Dά����       
        XG_V=zeros(Np,Nd);
        [valueBest,indexBest] = min(fit(XG));
        
        Fs=normrnd(Mf,Qf,Np,1);
            
        CRs=CRmin+(CRmax-CRmin)*rand(Np,Nd);
        xd_Min=min(XG,[],1);
        xd_Max=max(XG,[],1);  
       for G=1:Ni
            nn=nn+1;
            
            XG_n=xd_Min+(xd_Max-xd_Min).*rand(Np,Nd);          
            tempN=fit(XG_n)<fit(XG);
            tempN=repmat(tempN,1,Nd);
            XG=XG_n.*tempN-XG.*(tempN-1);
            
            r=(G/Ni)^(1/4);
            Fs=normrnd(mean(Fs),Qf,Np,1);
            meanCR=mean(CRs,2);
            CRs=(meanCR-Qcr)+(2*Qcr)*rand;
            Fs=setWithInAre(Fs,[Fmin,Fmax]);
%             Fs(Fs>Fmax)=Fmin+(Fmax-Fmin)*rand;
%             Fs(Fs<Fmin)=Fmin+(Fmax-Fmin)*rand;
            CRs=setWithInAre(CRs,[Fmin,Fmax]);
%             CRs(CRs>CRmax)=CRmin+(CRmax-CRmin)*rand;
%             CRs(CRs<CRmin)=CRmin+(CRmax-CRmin)*rand;
            %% %%%%%%%%%%%%%%%%%%%%%%----�������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
            for n=1:Np
                R3=randperm(Np,3);
                    if rand <r
                        son=XG(indexBest,:)+ Fs(n)*(XG(R3(2),:) - XG(R3(1),:));
                        testBest=testBest+1;
                    else
                        son = XG(R3(3),:) + Fs(n)*(XG(R3(2),:) - XG(R3(1),:));
                        testRand=testRand+1;
                    end
               
                    testRe=testRe+1;
                son=setWithInAre(son,[xmin,xmax]);
%                 son(son<xmin)=(xmax - xmin)*rand(1) + xmin;
%                 son(son>xmax)=(xmax - xmin)*rand(1) + xmin;
                XG_V(n,:)=son;
            end

            %% %%%%%%%%%%%%%%%%%%%%%---�������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            randR=randi([1,Nd],Np,1);
            dV=repmat(1:Nd,Np,1);

            tempX=(rand(Np,Nd)<=CRs) | (dV==randR);
            XG_C=XG_V.*tempX-XG.*(tempX-1);
            %% %%%%%%%%%%%%%%%%----ѡ�����---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tempC_i=fit(XG_C)<fit(XG);
            tempC=repmat(tempC_i,1,Nd);
            XG_next=XG_C.*tempC-XG.*(tempC-1);

            %% 
            XG = XG_next;
            fits=fit(XG);
            %����ָ�ҳ���С��fitness �� particle��indox
            [valueBest,indexBest] = min(fits);
            %���汾�ε�����Сfitness
            fitBests(t,G) = valueBest;
            nItersions(t)=G;
            if valueBest==0
                break;
            end
            
            
            %% ready next itermsion
            XGG=XG;
            XG_end_fromMean=mean(XG);
            if fit(XG_end_fromMean)<fit(XG(Np,:))
                XGG(Npp,:)=XG_end_fromMean;
            end
            xd_Min=min(XGG,[],1);
            xd_Max=max(XGG,[],1); 
           
       end
        disp(['MDE��',num2str(t),'�β��� �ڵ�',num2str(G),'�ε����ҵ�����ֵ��',num2str(fitBests(t,G))]);
    end

end

function y=setWithInAre(x,a)
c=(x>=a(1))|(x<=a(2));
y=x.*c+(a(1)+(a(2)-a(1))*rand(size(x))).*(1-c);

end

