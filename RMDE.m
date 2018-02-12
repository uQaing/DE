% NI the number of Iteration 
% popsize the number of particles
% D the number of dimension
% F
% CR
% fit the function of fitness
function [fitBests,object]=RMDE(Nt,Ni,Np,Nd,F,CR,range,fit)
disp(fit);

    xmin= -range;
    xmax = range;
    
    P=0.9;
    B=0.5;

    fitBests=zeros(Nt,Ni);
    XG_V=zeros(Np,Nd);
    object=zeros(Nt,1);
    for t=1:Nt
        XG = xmin+(xmax-xmin)*rand(Np,Nd);  %����Np��Dά����
		[valueBest,indexBest] = min(fit(XG));
        for i=1:Ni
            %% %%%%%%%%%%%%%%%%%%%%%%----�������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for n=1:Np
                R3=randperm(Np,3);
                if rand <P
                    if rand <B
                        son=XG(indexBest,:)+ F*(XG(R3(2),:) - XG(R3(1),:));
                    else
                        son = XG(R3(3),:) + F*(XG(R3(2),:) - XG(R3(1),:));
                    end
                else
                    son=xmin+(xmax-xmin)*rand(1,Nd);
                end
                son=setWithInAre(son,[xmin,xmax]);
%                 son(son<xmin)=(xmax - xmin)*rand(1) + xmin;
%                 son(son>xmax)=(xmax - xmin)*rand(1) + xmin;
                XG_V(n,:)=son;
            end
           
            %% %%%%%%%%%%%%%%%%%%%%%---�������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            randR=randi([1,Nd],Np,1);
            dV=repmat(1:Nd,Np,1);

            tempX=(rand(Np,Nd)<=CR) | (dV==randR);
            XG_C=XG_V.*tempX-XG.*(tempX-1);

            %% %%%%%%%%%%%%%%%%----ѡ�����---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tempC=fit(XG_C)<fit(XG);
            tempC=repmat(tempC,1,Nd);
            XG_next=XG_C.*tempC-XG.*(tempC-1);

            %% ���Խ����������
            XG = XG_next;
            %����ָ�ҳ���С��fitness �� particle��indox
            [valueBest,indexBest] = min(fit(XG));
            %���汾�ε�����Сfitness
            fitBests(t,i) = valueBest;
            object(t)=i;
            if valueBest==0
                break;
            end
            
        end
        disp(['RMDE��',num2str(t),'�β��� �ڵ�',num2str(i),'�ε����ҵ�����ֵ��',num2str(fitBests(t,i))]);
    end
end

function y=setWithInAre(x,a)
c=(x>=a(1))|(x<=a(2));
y=x.*c+(a(1)+(a(2)-a(1))*rand(size(x))).*(1-c);
end
