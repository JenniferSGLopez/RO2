format compact;
close all;
fclose all;

nq = 200;
gq = nq;
hq = nq;

MCtrials=1000; %1000

data = importdata('Level3_RO2SufficientSubSampleDetermination_9Nov2018.txt',' ',1);

D=data.data;
w=D(:,6)<.8;
fprintf(' %d cases with insufficient power.\n',sum(w));
D(:,1)=0;
D(:,6)=0;
D=unique(D(w,:),'rows');
fprintf(' %d cases to run.\n',size(D,1));

% Run first part of sample set
D=D(1:floor(size(D,1)/4),:);
fprintf(' %d cases to run.\n',size(D,1));

% % Run second part of sample set
% D=D((floor(size(D,1)/4)+1):((2*floor(size(D,1)/4)+1)),:);
% fprintf(' %d cases to run.\n',size(D,1));
%
% % Run third part of sample set
% D=D((2*floor(size(D,1)/4)+2):((3*floor(size(D,1)/4)+2)),:);
% fprintf(' %d cases to run.\n',size(D,1));
%
% % Run fourth part of sample set
% D=D=D((3*floor(size(D,1)/4)+3):((4*floor(size(D,1)/4))),:);
% fprintf(' %d cases to run.\n',size(D,1));

gamma_v=D(:,1);
Case_v=D(:,2);
scen_v=D(:,3);
ES_v=D(:,4);
e_var_v=D(:,5);
Power_v=D(:,6);
n_v=D(:,7);
g_v=D(:,8);
h_v=D(:,9);

sigma2_a_v=1*(Case_v==1) + 0.5*(Case_v==2) + 1*(Case_v==3);
sigma2_x_v=1*(Case_v==1) + 1*(Case_v==2) + 0.5*(Case_v==3);
sigma2_w_v=1*(Case_v==1) + 1*(Case_v==2) + 0.5*(Case_v==3);

%sigma2_a_v=1*(Case_v==1) + 0.5*(Case_v==2) + 1*(Case_v==3) + 2*(Case_v==4);
%sigma2_x_v=1*(Case_v==1) + 1*(Case_v==2) + 0.5*(Case_v==3) + 2*(Case_v==4);

sigma2_e_v=4*(e_var_v==1) + 1*(e_var_v==2) + 2*(e_var_v==3);
sigma2_r_v= sigma2_e_v.*(2*(scen_v==1)+ 1*(scen_v==2)+ 1*(scen_v==3));
sigma2_u_v= sigma2_e_v.*(1*(scen_v==1)+ 1*(scen_v==2)+ 1*(scen_v==3));

%sigma2_e_v=0.5*(e_var_v==1) + 2*(e_var_v==2);
%sigma2_e_v=0.5*(e_var_v==1) + 1*(e_var_v==2) + 2*(e_var_v==3);
%sigma2_r_v= sigma2_e_v.*(2*(scen_v==1)+.5*(scen_v==2)+(scen_v==3));

%% Set up Monte Carlo variables
MC_history=nan(MCtrials*size(D,1),47);

rng('shuffle');

for i=1:size(D,1)
    Case=Case_v(i);
    scen=scen_v(i);
    ES=ES_v(i);
    e_var=e_var_v(i);
    n=n_v(i);
    g=g_v(i);
    h=h_v(i);
    sigma2_a=sigma2_a_v(i);
    sigma2_x=sigma2_x_v(i);
    sigma2_w=sigma2_w_v(i);
    sigma2_e=sigma2_e_v(i);
    sigma2_r=sigma2_r_v(i);
    sigma2_u=sigma2_u_v(i);
    
    % Run Monte Carlo trials
    for MCk=1:MCtrials
        fprintf(1,'Now running Monte Carlo iteration %d of %d', MCk,MCtrials);
        tic
        % Generate Data based on Small Samples
        a = randn(n*g*h,1)*sqrt(sigma2_a);
        x= repmat(randn(1,g*h),n,1)*sqrt(sigma2_x);
        x=x(:);
        w= repmat(randn(1,h),n*g,1)*sqrt(sigma2_w);
        w=w(:);
        
        e = randn(n*g*h,1)*sqrt(sigma2_e);
        %r = normrnd(mu,sigma)
        r= repmat(randn(1,g*h),n,1)*sqrt(sigma2_r);
        r=r(:);
        u= repmat(randn(1,h),n*g,1)*sqrt(sigma2_u);
        u=u(:);
        
        %         Y = g_000 + g_001.*w + g_010.*x + (g_011.*x).*w + g_100.*a ...
        %                + (g_101.*a).*w + (g_110.*a).*x + ((g_111.*a).*x).*w ...
        %                + u_00k + u_01k.*x + u_10k.*a + (u_11k.*a).*x + r_0jk ...
        %                + r_1jk.*a + e_ijk;
        y = ES + ES.*w + ES.*x + ES.*x.*w + ES.*a + ES.*a.*w + ES.*a.*x + ...
            ES.*a.*x.*w + u + u.*x + u.*a + u.*a.*x + r + r.*a + e;
        
        % Sample Statistics
        a_mean = mean(a);
        a_var = var(a);
        x_mean = mean(x);
        x=reshape(x,[n g*h]);
        x_var=mean(var(x,0,2));
        %x_var = var(x);
        w_mean = mean(w);
        w=reshape(w,[n*g h]);
        w_var=mean(var(w,0,2));
        %w_var = var(w);
        %y_mean = mean(y);
        %y_var = var(y);
        
        a = randn(nq*gq*hq,1)*sqrt(a_var);
        x= repmat(randn(1,gq*hq),nq,1)*sqrt(x_var);
        x=x(:);
        w= repmat(randn(1,hq),nq*gq,1)*sqrt(w_var);
        w=w(:);
        
        
        %         % Multivariate Nomral Distribution Sample Statisitcs
        %         N = n*g*h -1;
        %         covAY = (1/N)*(a-a_mean)'*(y-y_mean);
        %         covXY = (1/N)*(x-x_mean)'*(y-y_mean);
        %         covWY = (1/N)*(w-w_mean)'*(y-y_mean);
        %
        %         mu = [a_mean,x_mean,w_mean,y_mean];
        %
        %         sigma = [a_var 0     0      covAY;
        %                  0     x_var 0      covXY;
        %                  0     0     w_var  covWY;
        %                  covAY covXY covWY  y_var];
        %         m=min(eig(sigma));
        %         if m<0
        %             m;
        %             sigma=sigma-m*eye(4);
        %         end
        
        ngh_reml=nq^3;
        %ngh_reml=[30,50,100]
        %             M = mvnrnd(mu,sigma,ngh_reml);
        %             a = M(:,1);
        %             x = M(:,2);
        %             w = M(:,3);
        %             y = M(:,4);
        e = randn(nq*gq*hq,1)*sqrt(sigma2_e);
        r= repmat(randn(1,gq*hq),nq,1)*sqrt(sigma2_r);
        r=r(:);
        u= repmat(randn(1,hq),nq*gq,1)*sqrt(sigma2_u);
        u=u(:);
        
        %% REML
        % Response
        %            Y = g_000 + g_001.*w + g_010.*x + (g_011.*x).*w + g_100.*a ...
        %                + (g_101.*a).*w + (g_110.*a).*x + ((g_111.*a).*x).*w ...
        %                + u_00k + u_01k.*x + u_10k.*a + (u_11k.*a).*x + r_0jk ...
        %                + r_1jk.*a + e_ijk;
        
        g_000 = ES;
        g_001 = ES;
        g_010 = ES;
        g_011 = ES;
        g_100 = ES;
        g_101 = ES;
        g_110 = ES;
        g_111 = ES;
        
        Y = g_000 + g_001.*w + g_010.*x + (g_011.*x).*w + g_100.*a ...
            + (g_101.*a).*w + (g_110.*a).*x + ((g_111.*a).*x).*w ...
            + u + u.*x + u.*a + (u.*a).*x + r + r.*a + e;
        
        % Design Matricies
        
        N = nq*gq*hq;
        X = [ones(N,1) w x x.*w a a.*w a.*x a.*x.*w];
        %            Z = {[ones(N,1) x a a.*x] [ones(N,1) a]};
        Z = {ones(N,1) x a a.*x ones(N,1) a};
        
        Level3 = repmat(1:hq,nq*gq,1);
        Level3 = Level3(:);
        
        Level2 = repmat(1:hq*gq,nq,1);
        Level2 = Level2(:);
        
        %             Level1 = repmat(1:(gq),nq,1);
        %             Level1 = Level1(:);
        %             %Level1 = nominal(Level1);
        %             %GG = Level1;
        %             GG = {Level3 Level2 Level1};
        
        %            GG = {Level3,Level2};
        
        GG = {Level3,Level3,Level3,Level3,Level2,Level2};
        
        %% Linear Mixed Model Fit Method
        %            optoptions=optimoptions('fminunc','Display','iter','Algorithm','quasi-newton');
        optoptions=optimoptions('fminunc','Algorithm','quasi-newton');
        
        %% INSERT PART 7 of Level-3 RO2 Page 4
        lme = fitlmematrix(X,Y,Z,GG,'FixedEffectPredictors',...
            {'g_000', 'g_001','g_010','g_011','g_100','g_101','g_110','g_111'},'RandomEffectPredictors',...
            {{'u_00k'},{'u_01k'},{'u_10k'},{'u_11k'},{'r_0j'},{'r_1j'}},'RandomEffectGroups', {'Unit','Unit','Unit','Unit','Team','Team'},...
            'FitMethod','REML','Optimizer','fminunc','OptimizerOptions',optoptions);
        
        MC_history((i-1)*MCtrials +MCk,:)=[n g h nthroot(ngh_reml,3) Case scen ES e_var ...
            lme.Coefficients.Estimate' lme.Coefficients.SE' ...
            lme.Coefficients.tStat' lme.Coefficients.pValue' sqrt(cell2mat(lme.covarianceParameters))' sqrt(lme.MSE)];
        
        %            MC_history=[MC_history; n g h ngh_reml Case scen ES e_var ...
        %                lme.Coefficients.Estimate' lme.Coefficients.SE' ...
        %                lme.Coefficients.tStat' lme.Coefficients.pValue' lme.mse];
        toc
    end
end


csvwrite('Level3_RO2_REML200part1Block_11Dec2018.csv', MC_history);

