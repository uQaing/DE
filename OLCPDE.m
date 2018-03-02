
function [fitBests,nItersions]= OLCPDE(Nt,Ni,Np,Nd,F,CR,upLow,fit)
disp(fit);
Fmin=0.1;
Fmax=0.8;
CRmin=0.3;
CRmax=0.9;
PBest=0.9;
BBest=0.5;
u=4;
xmin= -upLow;
xmax = upLow;

fitBests=zeros(Nt,Ni);
nItersions=zeros(Nt,1);



for k=1:Nt
    Y=zeros(Np,Nd);
Y(1,:)=rand(1,Nd);
for np=2:Np
    Y(np,:)=u*Y(np-1,:).*(1-Y(np-1,:));
end
X=xmin+Y.*(xmax-xmin);

temp_fit=fit(X);
[~,temp_fit_index]=sort(temp_fit,1);% 排序

elite_Index=temp_fit_index(1:(Np/2));
normal_Index=temp_fit_index((Np/2)+1:end);

XG_elite=X(elite_Index,:);
XG_normal=X(normal_Index,:);
XG=[XG_elite;XG_normal;XG_elite];
    for G=1:Ni
        %% %%%%%%%%%%%%%%%%%%%%%%----变异操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for n=1:Np/2
            R2s(n,:)=randperm((Np/2),2);
            R3s(n,:)=randperm((Np/2),3);
            R5s(n,:)=randperm((Np/2),5);
        end
        [~,best_v_index]=min(fit(XG_elite));
        son_elite = XG_elite(best_v_index,:) + F*(XG_elite(R3s(:,2),:) - XG_elite(R3s(:,1),:));
        son_normal=XG_normal(R5s(:,1),:)+F*(XG_normal(R5s(:,2),:)+XG_normal(R5s(:,3),:)-XG_normal(R5s(:,4),:)-XG_normal(R5s(:,5),:));
        son_parallel=XG_elite(best_v_index,:) + F*(XG_normal(R2s(:,2),:) - XG_normal(R2s(:,1),:));
        
        
        XG_V_elite=setWithInAre(son_elite,[xmin,xmax]);
        XG_V_normal=setWithInAre(son_normal,[xmin,xmax]);
        XG_V_parallel=setWithInAre(son_parallel,[xmin,xmax]);
        XG_V=[XG_V_elite;XG_V_normal;XG_V_parallel];
        
        randR=randi([1,Nd],Np*1.5,1);
        dV=repmat(1:Nd,Np*1.5,1);
        
        tempX=(rand(Np*1.5,Nd)<=CR) | (dV==randR);
        XG_C=XG_V.*tempX-XG.*(tempX-1);
        
        % inv study 选择
        XG_C_normal=XG_C((Np/2)+1:Np,:);
        XG_normal_inv=max(XG_normal)+min(XG_normal)-XG_normal;
        XG_C_normal_inv=max(XG_C_normal)+min(XG_C_normal)-XG_C_normal;
        
        XG_inv=[XG_normal;XG_normal_inv;XG_C_normal;XG_C_normal_inv];
        temp_fit=fit(XG_inv);
        [~,temp_fit_index]=sort(temp_fit,1);% 排序
        
        inv_Index=temp_fit_index(1:(Np/2));
        XG_next_inv_elite=XG_inv(inv_Index,:);
        
        
        %% %%%%%%%%%%%%%%%%----选择操作 只对精英和并行的选择---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tempC=fit(XG_C)<fit(XG);
        tempC=repmat(tempC,1,Nd);
        XG_next=XG_C.*tempC-XG.*(tempC-1);
        
        XG_next((Np/2)+1:Np,:)= XG_next_inv_elite ;          %% 测试进化后的粒子
        XG = XG_next;
        %从中指找出最小的fitness 和 particle的indox
        [valueBest,indexBest] = min(fit(XG));
        %保存本次迭代最小fitness
        fitBests(k,G) = valueBest;
        nItersions(k)= G;
        
        
        if valueBest==0
            break;
        end
        temp_fit=fit(XG);
        [~,temp_fit_index]=sort(temp_fit,1);% 排序
        
        elite_Index=temp_fit_index(1:(Np/2));
        normal_Index=temp_fit_index(Np+1:end);
        
        XG_elite=XG(elite_Index,:);
        XG_normal=XG(normal_Index,:);
        XG=[XG_elite;XG_normal;XG_elite];
    end
    disp(['OLCPDE第',num2str(k),'次测试 在第',num2str(G),'次迭代找到最优值：',num2str(fitBests(k,G))]);
    
end
o=0;
end


