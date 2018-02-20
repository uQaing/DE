% NI the number of Iteration 
% popsize the number of particles
% D the number of dimension
% F
% CR
% fit the function of fitness
function [fitBests,nItersions]=STDE(Nt,Ni,Np,Nd,F,CR,range,fit)
%     flag=0;
%     valueBestOld=zeros(1,1);
tic;
disp(fit);
    xmin= -range;
    xmax = range;
    R3s=zeros(Np,3);
    fitBests=zeros(Nt,Ni);
    nItersions=zeros(Nt,1);
    for k=1:Nt
%         if Nd==2
%             x1= -5+15*rand(Np,1);
%             x2= 0+15*rand(Np,1);
%             XG=[x1,x2];
%         else
%         XG = xmin+(xmax-xmin)*rand(Np,Nd);  %产生Np个D维向量
% end
       XG = xmin+(xmax-xmin)*rand(Np,Nd);  %产生Np个D维向量
        
        
        for G=1:Ni
            %% %%%%%%%%%%%%%%%%%%%%%%----变异操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for n=1:Np
                R3s(n,:)=randperm(Np,3);
            end
                son = XG(R3s(:,3),:) + F*(XG(R3s(:,2),:) - XG(R3s(:,1),:));
                son=setWithInAre(son,[xmin,xmax]);
%                 son(son<xmin)=(xmax - xmin)*rand(1) + xmin;
%                 son(son>xmax)=(xmax - xmin)*rand(1) + xmin;
                XG_V=son;
            

            %% %%%%%%%%%%%%%%%%%%%%%---交叉操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            randR=randi([1,Nd],Np,1);
            dV=repmat(1:Nd,Np,1);

            tempX=(rand(Np,Nd)<=CR) | (dV==randR);
            XG_C=XG_V.*tempX-XG.*(tempX-1);
 
            %% %%%%%%%%%%%%%%%%----选择操作---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tempC=fit(XG_C)<fit(XG);
            tempC=repmat(tempC,1,Nd);
            XG_next=XG_C.*tempC-XG.*(tempC-1);
            

            %% 测试进化后的粒子
            XG = XG_next;
            %从中指找出最小的fitness 和 particle的indox
            [valueBest,indexBest] = min(fit(XG));
            %保存本次迭代最小fitness
            fitBests(k,G) = valueBest;    
            nItersions(k)= G;
%             valueBestOld(flag==0)=valueBest;
%             if valueBest<=endBestValue
%                 break;
%             elseif valueBestOld==valueBest
%                 flag=flag+1;
%             elseif valueBestOld ~= valueBest
%                 flag=0;
%                  valueBestOld = valueBest;
%             end
%             
            if valueBest==0
                break;
            end

                
        end
        disp(['STDE第',num2str(k),'次测试 在第',num2str(G),'次迭代找到最优值：',num2str(fitBests(k,G))]);
    end

toc
end
function y=setWithInAre(x,a)
c=(x>=a(1))&(x<=a(2));
y=x.*c+(a(1)+(a(2)-a(1))*rand(size(x))).*(1-c);
end
