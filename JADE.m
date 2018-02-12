% NI the number of Iteration
% popsize the number of particles
% D the number of dimension
% F
% CR
% fit the function of fitness
function [fitBests,nItersions]=JADE(Nt,Ni,Np,Nd,F,CR,range,fit)
tic;
disp(fit);

xmin= -range;
xmax = range;


fitBests=zeros(Nt,Ni);
nItersions=zeros(Nt,1);

c=0.1;    
for k=1:Nt
    uCR=0.5;
    uF=0.5;
    successindex=1;
    inferiorIndex=1;
    XG = xmin+(xmax-xmin)*rand(Np,Nd);  %����Np��Dά���� 
    XGNext=zeros(Np,Nd); 
    XA=[];
    for g=1:Ni
        %% %%%%%%%%%%%%%%%%%%%%%%----�������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        temp_fit=fit(XG);
        [B,I]=sort(temp_fit);% ����
        
        CRis =normrnd(uCR,0.1,2*Np,1);
        CRis=CRis(CRis>=0&CRis<=1);
        [r,~]=size(CRis);
        while(r<Np)
           CRis =normrnd(uCR,0.1,2*Np,1);
        CRis=CRis(CRis>=0&CRis<=1);
        [r,~]=size(CRis);
        end
        CRis=CRis(1:Np);
        
        Fis=cauchycdf(ones(Np,1)*uF,0.1);
        Fis((Fis<=0)|(Fis>=1))=1;
        for i=1:Np
%             CRi =normrnd(uCR,0.1);
%             while((CRi<0)||(CRi>1))
%                 CRi =normrnd(uCR,0.1);
%             end
%             Fi=cauchycdf(uF,0.1);
%             Fi((Fi<=0)||(Fi>=1))=1;
            
            XGp_bestg=XG(I(randperm(5,1)),:);
            XGr1g=XG(randperm(Np,1),:);
            PA=[XG;XA];
            XGr2g=PA(randperm(Np,1),:);
            son = XG(i,:) +Fis(i)*(XGp_bestg-XG(i,:))+Fis(i)*(XGr1g-XGr2g);
            son((son<xmin)|(son>xmax))=(xmax - xmin)*rand(1) + xmin;
            XG_V=son;            
        %% %%%%%%%%%%%%%%%%%%%%%---�������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        tempX=(rand(1,Nd)<=CRis(i)) | (randperm(Nd,1)==[1:Nd]);
        XG_C=XG_V.*tempX-XG(i,:).*(tempX-1);
       %% %%%%%%%%%%%%%%%%----ѡ�����---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Nfit=fit(XG_C);
       Sfit=fit(XG(i,:));
        if(Nfit<Sfit)
            fitBests(k,g)=Nfit;
            XGNext(i,:)=XG_C;
            SF(successindex)=Fis(i);
            SCR(successindex)=CRis(i);
            successindex=successindex+1;
            XA(inferiorIndex,:)=XG(i,:);
        else
            fitBests(k,g)=Sfit;
            XGNext(i,:)=XG(i,:);
        end
           if(inferiorIndex<=Np)
            inferiorIndex=inferiorIndex+1;
            else
                inferiorIndex=randperm(Np,1);
            end 
        end %end Np
            uF=(1-c)*uF+c*(sum(SF.^2)/sum(SF));                    
            uCR=(1-c)*uCR+c*(mean(SCR));
            XG=XGNext;
    end %end Ni
        
        disp(['JADE��',num2str(k),'�β��� �ڵ�',num2str(g),'�ε����ҵ�����ֵ��',num2str(fitBests(k,g))]);
end %end Nt
    toc
end



