% NI the number of Iteration
% popsize the number of particles
% D the number of dimension
% F
% CR
% fit the function of fitness
function [fitBests,nItersions]=MRDE(Nt,Ni,Np,Nd,Kk,Ww,upLow,fit)


disp(fit);

Fmin=0.1;
Fmax=0.8;
CRmin=0.3;
CRmax=0.9;
PBest=0.9;
BBest=0.5;

xmin= -upLow;
xmax = upLow;

fitBests=zeros(Nt,Ni);
nItersions=zeros(Nt,1);
for t=1:Nt
    G=0;
    plan1=0;
plan2=0;
    nn=0;
    test0=0;
    test1=0;
    
    testBest=0;
    testRand=0;
    testRe=0;
    
    XG = xmin+(xmax-xmin)*rand(Np,Nd);  %产生Np个D维向量
    XG_V=zeros(Np,Nd);
    [valueBest,indexBest] = min(fit(XG));
    
    % Fs and CRs in first iterasion.
    Fs=Fmin+(Fmax-Fmin)*rand(Np,1);
    CRs=CRmin+(CRmax-CRmin)*rand(Np,Nd);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ps=PBest+0.1*rand(Np,1);
    Bs=BBest+0.1*rand(Np,1);
    
    Ps_success=Ps;
    Bs_success=Bs;
    Fs_success=Fs;
    CRs_success=CRs;
    %         for G=1:Ni
    while G<Ni
        G=G+1;
        nn=nn+1;
        xd_Min=min(XG,[],1);
        xd_Max=max(XG,[],1);
        XG_n=xd_Min+(xd_Max-xd_Min).*rand(Np,Nd);
        tempN=fit(XG_n)<fit(XG);
        tempN=repmat(tempN,1,Nd);
        XG=XG_n.*tempN-XG.*(tempN-1);
        
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%----变异操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for n=1:Np
            R3=randperm(Np,3);
            if rand <Ps(n)
                if rand <Bs(n)
                    son=XG(indexBest,:)+ Fs(n)*(XG(R3(2),:) - XG(R3(1),:));
                    testBest=testBest+1;
                else
                    son = XG(R3(3),:) + Fs(n)*(XG(R3(2),:) - XG(R3(1),:));
                    testRand=testRand+1;
                end
            else
                son=xmin+(xmax-xmin)*rand(1,Nd);
                testRe=testRe+1;
            end
            son(son<xmin)=(xmax - xmin)*rand(1) + xmin;
            son(son>xmax)=(xmax - xmin)*rand(1) + xmin;
            XG_V(n,:)=son;
        end
        
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
        %             K=0;
        K=0.9+0.05*rand;
        %                 K=0.05;
        if usable > K
            plan1=plan1+1;
            
            test0=test0+1;
            
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
            u=sum(temp,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            d=sum(temp.^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            CRs=tempC_i.*CRs+(1-tempC_i).*(d./u);
            Ps_success=Ps;
            Bs_success=Bs;
            Fs_success=Fs;
            CRs_success=CRs;
            
        else
            %                 G = 1;
            plan2=plan2+1;
            test1=test1+1;
            
            Ps=normrnd(mean(Ps_success),0.1,Np,1);
            Bs=normrnd(mean(Bs_success),0.1,Np,1);
            Fs=normrnd(mean(Fs_success),0.1,Np,1);
            
            
            % Fs and CRs in first iterasion.
            %         Fs=Fmin+(Fmax-Fmin)*rand(Np,1);
            %         CRs=CRmin+(CRmax-CRmin)*rand(Np,1);
            %         Ps=0.9+0.1*rand(Np,1);
            %         Bs=0.5+0.5*rand(Np,1);
            %
            %
            avgFit=mean(fits);
            minFit=min(fits);
            maxFit=max(fits);
            %
            avgFits=repmat(avgFit,Np,1);
            minFits=repmat(minFit,Np,1);
            maxFits=repmat(maxFit,Np,1);
            
            W=1;
            J=100;
            %         q=((maxFits-minFits).^W)./(((maxFits-avgFits).^W)-((minFits-avgFits).^W));
            %         t1=q.*((((fits-avgFits)./(maxFits-minFits)).^W)-(((minFits-avgFits)./(maxFits-minFits)).^W));
            %         t1=real(t1);
            %
            %         tempCR=(fits<=avgFits);
            %         CRs=CRs.*tempCR-(CRmin+(CRmax-CRmin)*t1).*(1-tempCR);
            
            %         CRs=(CRmin+(CRmax-CRmin)*t1);
            
            
%                     t2=((fits-avgFits)./(maxFits-minFits));
%                     t2=t2.^W;
%                     tempCR=(fits<=avgFits);
%                     CRs=CRs.*tempCR-(CRmin+(CRmax-CRmin)*t2).*(1-tempCR);
            
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
        
    end
    %         test0
    %         test1
    %         testBest
    %         testRand
    %         testRe
    %         nn
    disp(['MRDE第',num2str(t),'次测试 在第',num2str(G),'次迭代找到最优值：',num2str(fitBests(t,G)),' plan1= ',num2str(plan1),'plan2= ',num2str(plan2)]);
    %         disp(sm);
end

end
