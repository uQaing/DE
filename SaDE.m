function [v,i]=SaDE(Nt,Ni,Np,Nd,F,CR,range,fit)
tic;
disp(fit);
    xmin= -range;
    xmax = range;
    Fmin=-0.4;
    Fmax=1.4;
    
    for k=1:Nt
     R5s=zeros(Np,5);
    fitBests=zeros(Nt,Ni);
    nItersions=zeros(Nt,1);
    LP=50;
    lp=1;
    SuccessMemory=zeros(LP,4);
    FailureMemory=zeros(LP,4);
    CRMemory=zeros(Np,4);
    CRm=[0.5,0.5,0.5,0.5];
    CRs=normrnd(0.5,0.1,Np,4);
    P=[1/4,1/4,1/4,1/4];
        
        XG = xmin+(xmax-xmin)*rand(Np,Nd);  %����Np��Dά����
        [~,indexBest] = min(fit(XG));
        for G=1:Ni
            alphabet = [1 2 3 4]; 
            stratetyFlag=randsrc(1,1,[alphabet; P]);
            %% %%%%%%%%%%%%%%%%%%%%%%----�������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Fs=normrnd(0.5,0.3,Np,1);
            Fs(Fs<Fmin|Fs>Fmax)=(Fmax -Fmin)*rand(1) + Fmin;
            for n=1:Np
                R5s(n,:)=randperm(Np,5);
            end
            if stratetyFlag==4
            son=  XG+ rand(Np,1).*(XG(R5s(:,1),:) - XG)+Fs.*(XG(R5s(:,2),:) - XG(R5s(:,3),:));
            son(son<xmin)=(xmax - xmin)*rand(1) + xmin;
            son(son>xmax)=(xmax - xmin)*rand(1) + xmin;
            XG_C=son;
            else
                switch stratetyFlag
                    case 1
                        son = XG(R5s(:,3),:) + Fs.*(XG(R5s(:,2),:) - XG(R5s(:,1),:));
                    case 2
                        son = XG + Fs.*(XG(indexBest,:) - XG)+Fs.*(XG(R5s(:,1),:) - XG(R5s(:,2),:))+Fs.*(XG(R5s(:,3),:) - XG(R5s(:,4),:));
                    case 3
                        son = XG(R5s(:,1),:) + Fs.*(XG(R5s(:,2),:) - XG(R5s(:,3),:))+Fs.*(XG(R5s(:,4),:) - XG(R5s(:,5),:));
                end
            
            son(son<xmin)=(xmax - xmin)*rand(1) + xmin;
            son(son>xmax)=(xmax - xmin)*rand(1) + xmin;
            XG_V=son;
            %% %%%%%%%%%%%%%%%%%%%%%---�������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            randR=randi([1,Nd],Np,1);
            dV=repmat(1:Nd,Np,1);

            tempX=(rand(Np,Nd)<=CRs(:,stratetyFlag)) | (dV==randR);
            XG_C=XG_V.*tempX-XG.*(tempX-1);
 
            end
            %% %%%%%%%%%%%%%%%%----ѡ�����---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tempC=fit(XG_C)<fit(XG);  
            
            successNum=sum(tempC);
            failureNum=Np-successNum;
            if lp==LP
                lp=1;
            end
            SuccessMemory(lp,stratetyFlag)=successNum;
            FailureMemory(lp,stratetyFlag)=failureNum;
            CRMemory(:,stratetyFlag)=CRs(:,stratetyFlag);
            if G<LP
            CRs=normrnd(0.5,0.1,Np,4);
            else
                tt=repmat(median(CRMemory),Np,1);
                
                CRs=normrnd(tt,0.1,Np,4);
            end
            lp=lp+1;
           SuccessSum=sum(SuccessMemory,1);
           FailureSum=sum(FailureMemory,1);
           S=SuccessSum./(SuccessSum+FailureSum);
           S(isnan(S))=0;
           S=S+0.01;
           P=S./sum(S,2);
%            sp= sum(P,2);
%            if sp~=1
%                 error('sum P is not 1!!! ');  
%            end
            tempC=repmat(tempC,1,Nd);
            XG_next=XG_C.*tempC-XG.*(tempC-1);
            

            %% ���Խ����������
            XG = XG_next;
            %����ָ�ҳ���С��fitness �� particle��indox
            [valueBest,indexBest] = min(fit(XG));
            %���汾�ε�����Сfitness
            fitBests(k,G) = valueBest;    
            nItersions(k)= G;
         
            if valueBest==0
                break;
            end

                
        end
        disp(['SaDE��',num2str(k),'�β��� �ڵ�',num2str(G),'�ε����ҵ�����ֵ��',num2str(fitBests(k,G))]);
    end

toc
end