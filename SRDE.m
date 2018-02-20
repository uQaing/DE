% NI the number of Iteration 
% popsize the number of particles
% D the number of dimension
% F
% CR
% fit the function of fitness
function [fitBests,nItersions]=SRDE(Nt,Ni,Np,Nd,F,CR,range,fit)
disp(fit);
    
Fmin=0.1;
Fmax=0.8;
CRmin=0.1; 
CRmax=0.9;
P=0.9;
B=0.5;
xmin= -range;
xmax = range;
    
    fitBests=zeros(Nt,Ni);
    nItersions=zeros(Nt,1);
XG_V=zeros(Np,Nd);
    for k=1:Nt
        
        XG = xmin+(xmax-xmin)*rand(Np,Nd);  %产生Np个D维向量
        
        [valueBest,indexBest] = min(fit(XG));
        
        % Fs and CRs in first iterasion.
        Fs=Fmin+(Fmax-Fmin)*rand(Np,1);
        CRs=CRmin+(CRmax-CRmin)*rand(Np,1);
        
        for G=1:Ni
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%----变异操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
            for n=1:Np
                R3=randperm(Np,3);
                if rand <P
                    if rand <B
                        son=XG(indexBest,:)+ Fs(n)*(XG(R3(2),:) - XG(R3(1),:));
                    else
                        son = XG(R3(3),:) + Fs(n)*(XG(R3(2),:) - XG(R3(1),:));
                    end
                else
                    son=xmin+(xmax-xmin)*rand(1,Nd);
                end
                son=setWithInAre(son,[xmin,xmax]);
%                 son(son<xmin)=(xmax - xmin)*rand(1) + xmin;
%                 son(son>xmax)=(xmax - xmin)*rand(1) + xmin;
                XG_V(n,:)=son;
            end

               
            %% %%%%%%%%%%%%%%%%%%%%%---交叉操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            randR=randi([1,Nd],Np,1);
            dV=repmat(1:Nd,Np,1);

            tempX=(rand(Np,Nd)<=CRs) | (dV==randR);
            XG_C=XG_V.*tempX-XG.*(tempX-1);
            %% %%%%%%%%%%%%%%%%----选择操作---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tempC=fit(XG_C)<fit(XG);
            tempC=repmat(tempC,1,Nd);
            XG_next=XG_C.*tempC-XG.*(tempC-1);

            %% 
            XG = XG_next;
            fits=fit(XG);
            %从中指找出最小的fitness 和 particle的indox
            [valueBest,indexBest] = min(fits);
            %保存本次迭代最小fitness
            fitBests(k,G) = valueBest;
            nItersions(k)=G;
            if valueBest==0
                break;
            end
            
            %% create next Fs
            Fs=Fmin+(Fmax-Fmin)*rand(Np,1);
            
            %% create next CRs
            avgFit=mean(fits);
            minFit=min(fits);
            maxFit=max(fits);
            
            avgFits=repmat(avgFit,Np,1);
            minFits=repmat(minFit,Np,1);
            maxFits=repmat(maxFit,Np,1);
            tt=((fits-minFits)./(maxFits-minFits));
            tempCR=(fits<=avgFits);
            CRs=CRs.*tempCR-(CRmin+(CRmax-CRmin)*tt).*(1-tempCR);
            
        end
        disp(['SRDE第',num2str(k),'次测试 在第',num2str(G),'次迭代找到最优值：',num2str(fitBests(k,G))]);
    end


end

function y=setWithInAre(x,a)
c=(x>=a(1))&(x<=a(2));
y=x.*c+(a(1)+(a(2)-a(1))*rand(size(x))).*(1-c);
end
