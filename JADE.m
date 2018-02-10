% NI the number of Iteration
% popsize the number of particles
% D the number of dimension
% F
% CR
% fit the function of fitness
function [fitBests,nItersions]=JADE(Nt,Ni,Np,Nd,F,CR,range,fit)
disp(fit);

xmin= -range;
xmax = range;
R3s=zeros(Np,3);
index=zeros(Np,1);
fitBests=zeros(Nt,Ni);
nItersions=zeros(Nt,1);
uCRss=zeros(Np,Nd);
uFs=zeros(Np,1);
uCRss0=0.5;
uF0=0.5;

    kmmk=100;
    if Np<100
        kmmk=Np;
    end
    
for k=1:Nt
    xmin= -range;
    xmax = range;
    R3s=zeros(Np,3);
    index=zeros(Np,1);
    uCRss=zeros(Np,Nd);
    uFs=zeros(Np,1);
    uCRss0=0.5;
    uF0=0.5;
    
    
    XG = xmin+(xmax-xmin)*rand(Np,Nd);  %产生Np个D维向量
    CRss=normrnd(uCRss0,0.1,Np,Nd);
    Fs=normrnd(uF0,0.1,Np,1);
    [valueBest,indexBest] = min(fit(XG));

        
    for G=1:Ni
        %% %%%%%%%%%%%%%%%%%%%%%%----变异操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp_fit=fit(XG);
        temp_fit_sort=sort(temp_fit,1);% 排序
        xp=ceil(rand*kmmk);
        xp(xp==0)=1;
        temp_fit_min100=temp_fit_sort(1:xp);
        [l,ll]=size(temp_fit_min100);
        
        for n=1:Np
            
            
            R3s(n,:)=randperm(Np,3);
            
            ttmm=temp_fit_min100(randperm(l,1));
            tmtm=find(temp_fit==ttmm);
            index(n)= tmtm(1);
        end
        MM=(XG(index,:)-XG);
        
        
        %              while index(1)==R3s(1)||index(1)==R3s(2)
        %              R3s=randperm(Np,2);
        %              end
        
        son = XG(R3s(:,3),:) +Fs.*MM+Fs.*(XG(R3s(:,2),:) - XG(R3s(:,1),:));
        %                 son=XG(indexBest,:)+Fs.*(XG(R3s(:,2),:) - XG(R3s(:,1),:));
        son(son<xmin)=(xmax - xmin)*rand(1) + xmin;
        son(son>xmax)=(xmax - xmin)*rand(1) + xmin;
        XG_V=son;
        
        
        %% %%%%%%%%%%%%%%%%%%%%%---交叉操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        randR=randi([1,Nd],Np,1);
        dV=repmat(1:Nd,Np,1);
        
        tempX=(rand(Np,Nd)<=CRss) | (dV==randR);
        XG_C=XG_V.*tempX-XG.*(tempX-1);
        
        %% %%%%%%%%%%%%%%%%----选择操作---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tempC_i=fit(XG_C)<fit(XG);
        tempC=repmat(tempC_i,1,Nd);
        XG_next=XG_C.*tempC-XG.*(tempC-1);
        
        XA=XG.*tempC;
        
        %% 测试进化后的粒子
        XG = XG_next;
        %从中指找出最小的fitness 和 particle的indox
        [valueBest,indexBest] = min(fit(XG));
        %保存本次迭代最小fitness
        fitBests(k,G) = valueBest;
        nItersions(k)= G;
        
        
        usable=(sum(tempC_i)/Np);
        if usable > 0
            temp=tempC_i.*Fs;
            u=sum(temp);
            d=sum(temp.^2);
            c=rand(Np,1);
            uFs=(1-c).*uFs+c.*(d/u);
            
            
            d=sum(tempC_i.*CRss);
            u=sum(tempC_i);
            c=rand(Np,Nd);
            uCRss=(1-c).*uCRss+c.*(d/u);
            
            %                 CRss=normcdf(uCRss,0.1);
            %                 Fs=normcdf(uFs,0.1);
            CRss=uCRss;
            Fs=uFs;
            
        else
            CRss=normrnd(uCRss0,0.1,Np,Nd);
            Fs=normrnd(uF0,0.1,Np,1);
        end
        
        
    end
    disp(['JADE第',num2str(k),'次测试 在第',num2str(G),'次迭代找到最优值：',num2str(fitBests(k,G))]);
end


end

