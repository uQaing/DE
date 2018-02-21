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
    XG = xmin+(xmax-xmin)*rand(Np,Nd);  %产生Np个D维向量 
    XGNext=zeros(Np,Nd); 
    XA=[];
    for g=1:Ni
%         SF=[];
%         SCR=[];
        %% %%%%%%%%%%%%%%%%%%%%%%----变异操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        temp_fit=fit(XG);
        [B,I]=sort(temp_fit);% 升序
        
        CRis =normrnd(uCR,0.1,2*Np,1);
        CRis=CRis(CRis>=0&CRis<=1);
        [r,~]=size(CRis);
        while(r<Np)
           CRis =normrnd(uCR,0.1,2*Np,1);
        CRis=CRis(CRis>=0&CRis<=1);
        [r,~]=size(CRis);
        end
        CRis=CRis(1:Np);

        
        Fis=myCauchy(ones(Np,1)*uF,0.1);
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
            %             son((son<xmin)|(son>xmax))=(xmax - xmin)*rand(1) + xmin;
            son=setWithInAre(son,[xmin,xmax]);
            XG_V=son;
        %% %%%%%%%%%%%%%%%%%%%%%---交叉操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        tempX=(rand(1,Nd)<=CRis(i)) | (randperm(Nd,1)==[1:Nd]);
        XG_C=XG_V.*tempX-XG(i,:).*(tempX-1);
       %% %%%%%%%%%%%%%%%%----选择操作---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        disp(['JADE第',num2str(k),'次测试 在第',num2str(g),'次迭代找到最优值：',num2str(fitBests(k,g))]);
end %end Nt
    toc
end
function y=setWithInAre(x,a)
c=(x>=a(1))&(x<=a(2));
y=x.*c+(a(1)+(a(2)-a(1))*rand(size(x))).*(1-c);
end


