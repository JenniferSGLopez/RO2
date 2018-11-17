format compact;
close all;
fclose all;

nq = 50;
gq = nq;
hq = nq;
fq = nq;

MCtrials=1000; %1000

data = importdata('Level4_RO2SufficientSubSampleDetermination_10Oct2018.txt',' ',1);

D=data.data;
w=D(:,6)<.8;
fprintf(' %d cases with insufficient power.\n',sum(w));
D(:,1)=0;
D(:,6)=0;
D=unique(D(w,:),'rows');
fprintf(' %d cases to run.\n',size(D,1));

% Part 1
D=D((1:10),:); %Use small sample for testing
fprintf(' %d cases to run.\n',size(D,1));

% Part 2
%D=D((11:20),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 3
%D=D((21:30),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 4
%D=D((31:40),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 5
%D=D((41:50),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 6
%D=D((51:60),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 7
%D=D((61:70),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 8
%D=D((71:80),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 9
%D=D((81:90),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 10
%D=D((91:100),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 11
%D=D((101:110),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 12
%D=D((111:120),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 13
%D=D((121:130),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 14
%D=D((131:140),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 15
%D=D((141:150),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 16
%D=D((151:160),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 17
%D=D((161:170),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 18
%D=D((171:180),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 19
%D=D((181:190),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

% Part 20
%D=D((191:198),:); %Use small sample for testing
%fprintf(' %d cases to run.\n',size(D,1));

delta_v=D(:,1);
Case_v=D(:,2);
scen_v=D(:,3);
ES_v=D(:,4);
e_var_v=D(:,5);
Power_v=D(:,6);
n_v=D(:,7);
g_v=D(:,8);
h_v=D(:,9);
f_v=D(:,10);

sigma2_a_v=1*(Case_v==1) + 0.5*(Case_v==2) + 1*(Case_v==3);
sigma2_x_v=1*(Case_v==1) + 1*(Case_v==2) + 0.5*(Case_v==3);
sigma2_w_v=1*(Case_v==1) + 1*(Case_v==2) + 0.5*(Case_v==3);
sigma2_z_v=1*(Case_v==1) + 1*(Case_v==2) + 0.5*(Case_v==3);

%sigma2_a_v=1*(Case_v==1) + 0.5*(Case_v==2) + 1*(Case_v==3) + 2*(Case_v==4);
%sigma2_x_v=1*(Case_v==1) + 1*(Case_v==2) + 0.5*(Case_v==3) + 2*(Case_v==4);

sigma2_e_v=8*(e_var_v==1) + 1*(e_var_v==2) + 2*(e_var_v==3);
sigma2_r_v= sigma2_e_v.*(4*(scen_v==1)+ 1*(scen_v==2)+ 1*(scen_v==3));
sigma2_u_v= sigma2_e_v.*(2*(scen_v==1)+ 1*(scen_v==2)+ 1*(scen_v==3));
sigma2_v_v= sigma2_e_v.*(1*(scen_v==1)+ 1*(scen_v==2)+ 1*(scen_v==3));

%sigma2_e_v=0.5*(e_var_v==1) + 1*(e_var_v==2) + 2*(e_var_v==3);
%sigma2_r_v= sigma2_e_v.*(2*(scen_v==1)+.5*(scen_v==2)+(scen_v==3));


%% Set up Monte Carlo variables
MC_history=nan(MCtrials*size(D,1),88); % NEED TO CHANGE 74

rng('shuffle')

for i=1:size(D,1)
    Case=Case_v(i);
    scen=scen_v(i);
    ES=ES_v(i);
    e_var=e_var_v(i);
    n=n_v(i);
    g=g_v(i);
    h=h_v(i);
    f=f_v(i);
    sigma2_a=sigma2_a_v(i);
    sigma2_x=sigma2_x_v(i);
    sigma2_w=sigma2_w_v(i);
    sigma2_z=sigma2_z_v(i);
    sigma2_e=sigma2_e_v(i);
    sigma2_r=sigma2_r_v(i);
    sigma2_u=sigma2_u_v(i);
    sigma2_v=sigma2_v_v(i);
    
    
    % Run Monte Carlo trials
    for MCk=1:MCtrials
        fprintf(1,'Now running Monte Carlo iteration %d of %d\n', MCk,MCtrials);
        tic
        % Generate Data based on Small Samples
        a = randn(n*g*h*f,1)*sqrt(sigma2_a);
        x= repmat(randn(1,g*h*f),n,1)*sqrt(sigma2_x);
        x=x(:);
        w= repmat(randn(1,h*f),n*g,1)*sqrt(sigma2_w);
        w=w(:);
        z= repmat(randn(1,f),n*g*h,1)*sqrt(sigma2_z);
        z=z(:);
        
        e = randn(n*g*h*f,1)*sqrt(sigma2_e);
        %r = normrnd(mu,sigma)
        r= repmat(randn(1,g*h*f),n,1)*(sigma2_r);
        r=r(:);
        u= repmat(randn(1,h*f),n*g,1)*sqrt(sigma2_u);
        u=u(:);
        v= repmat(randn(1,f),n*g*h,1)*sqrt(sigma2_v);
        v=v(:);
        
        % MODEL
        %y = d0000 + d1000a + d0100x + d0010w + d0001z + d1100ax +
        %    d1010aw + d1001az + d0110xw + d0101xz + d0011wz + d1110axw +
        %    d1101axz + d1011awz + d0111xwz + d1111axwz +
        %    + v000l + v100la + v010lx + v001lw + v110lax + v101law
        %    + v011lxw + v111laxw
        %    + u_00kl + u_01kl.*x + u_10kl.*a + (u_11kl.*a).*x + r_0jkl ...
        %    + r_1jkl.*a + e_ijkl;
        
        y = ES + ES*a + ES*x + ES*w + ES*z + ES*a.*x + ES*a.*w + ES*a.*z ...
            + ES*x.*w + ES*x.*z + ES*w.*z + ES*a.*x.*w + ES*a.*x.*z ...
            + ES*a.*w.*z + ES*x.*w.*z + ES*a.*x.*w.*z ...
            + v + v.*a + v.*x + v.*w + v.*a.*x + v.*a.*w + v.*x.*w + v.*a.*x.*w ...
            + u + u.*x + u.*a + u.*a.*x ...
            + r + r.*a + e;
        
        % Sample Statistics
        a_mean = mean(a);
        a_var = var(a);
        x_mean = mean(x);
        x=reshape(x,[n g*h*f]);
        x_var=mean(var(x,0,2));
        %x_var = var(x);
        w_mean = mean(w);
        w=reshape(w,[n*g h*f]);
        w_var=mean(var(w,0,2));
        %w_var = var(w);
        z_mean = mean(z);
        z=reshape(z,[n*g*h f]);
        z_var=mean(var(z,0,2));
        %z_var = var(z);
        %y_mean = mean(y);
        %y_var = var(y);
        
        a = randn(nq*gq*hq*fq,1)*sqrt(a_var);
        x= repmat(randn(1,gq*hq*fq),nq,1)*sqrt(x_var);
        x=x(:);
        w= repmat(randn(1,hq*fq),nq*gq,1)*sqrt(w_var);
        w=w(:);
        z= repmat(randn(1,fq),nq*gq*hq,1)*sqrt(z_var);
        z=z(:);
        
        
        % Multivariate Nomral Distribution Sample Statisitcs
        %         N = n*g*h*f -1;
        %         covAY = (1/N)*(a-a_mean)'*(y-y_mean);
        %         covXY = (1/N)*(x-x_mean)'*(y-y_mean);
        %         covWY = (1/N)*(w-w_mean)'*(y-y_mean);
        %         covZY = (1/N)*(z-z_mean)'*(y-y_mean);
        
        %         mu = [a_mean,x_mean,w_mean,z_mean,y_mean];
        %
        %         sigma = [a_var  0     0      0     covAY;
        %             0      x_var 0      0     covXY;
        %             0      0     w_var  0     covWY;
        %             0      0     0      z_var covZY;
        %             covAY  covXY covWY  covZY y_var];
        %         m=min(eig(sigma));
        %         if m<0
        %             m;
        %             sigma=sigma-m*eye(5);
        %         end
        
        %         nq = 30;
        %         gq = nq;
        %         hq = nq;
        %         fq = nq;
        
        nghf_reml=nq^4;
        %ngh_reml=[30,50,100]
        %         M = mvnrnd(mu,sigma,nghf_reml);
        %         a = M(:,1);
        %         x = M(:,2);
        %         w = M(:,3);
        %         z = M(:,4);
        %         y = M(:,5);
        %
        e = randn(nq*gq*hq*fq,1)*sqrt(sigma2_e);
        r= repmat(randn(1,gq*hq*fq),nq,1)*sqrt(sigma2_r);
        r=r(:);
        u= repmat(randn(1,hq*fq),nq*gq,1)*sqrt(sigma2_u);
        u=u(:);
        v= repmat(randn(1,fq),nq*gq*hq,1)*sqrt(sigma2_v);
        v=v(:);
        
        %% REML
        % Response
        %y = d0000 + d1000a + d0100x + d0010w + d0001z + d1100ax +
        %    d1010aw + d1001az + d0110xw + d0101xz + d0011wz + d1110axw +
        %    d1101axz + d1011awz + d0111xwz + d1111axwz +
        %    + v000l + v100la + v010lx + v001lw + v110lax + v101law
        %    + v011lxw + v111laxw
        %    + u_00kl + u_01kl.*x + u_10kl.*a + (u_11kl.*a).*x + r_0jkl ...
        %    + r_1jkl.*a + e_ijkl;
        
        d_0000 = ES;
        d_1000 = ES;
        d_0100 = ES;
        d_0010 = ES;
        d_0001 = ES;
        d_1100 = ES;
        d_1010 = ES;
        d_1001 = ES;
        d_0110 = ES;
        d_0101 = ES;
        d_0011 = ES;
        d_1110 = ES;
        d_1101 = ES;
        d_1011 = ES;
        d_0111 = ES;
        d_1111 = ES;
        
        Y = d_0000 + d_1000*a + d_0100*x + d_0010*w + d_0001*z + d_1100*a.*x ...
            + d_1010*a.*w + d_1001*a.*z + d_0110*x.*w + d_0101*x.*z + d_0011*w.*z ...
            + d_1110*a.*x.*w + d_1101*a.*x.*z + d_1011*a.*w.*z + d_0111*x.*w.*z ...
            + d_1111*a.*x.*w.*z ...
            + v + v.*a + v.*x + v.*w + v.*a.*x ...
            + v.*a.*w + v.*x.*w + v.*a.*x.*w ...
            + u + u.*x + u.*a + u.*a.*x ...
            + r + r.*a + e;
        
        % Design Matricies
        N = nq*gq*hq*fq;
        X = [ones(N,1) a x w z a.*x a.*w a.*z x.*w x.*z w.*z a.*x.*w a.*x.*z a.*w.*z x.*w.*z a.*x.*w.*z];
        Z = {ones(N,1) a x w a.*x a.*w x.*w a.*x.*w ones(N,1) x a a.*x ones(N,1) a};
        %Z = {[ones(N,1) a x w a.*x a.*w x.*w a.*x.*w] [ones(N,1) x a a.*x] [ones(N,1) a]};
        
        Level4 = repmat(1:fq,nq*gq*hq,1);
        Level4 = Level4(:);
        
        Level3 = repmat(1:hq*fq,nq*gq,1);
        Level3 = Level3(:);
        
        Level2 = repmat(1:gq*hq*fq,nq,1);
        Level2 = Level2(:);
        
        %Level1 = repmat(1:(gq),nq,1);
        %Level1 = Level1(:);
        %Level1 = nominal(Level1);
        %GG = Level1;
        GG = {Level4,Level4,Level4,Level4,Level4,Level4,Level4,Level4,Level3,Level3,Level3,Level3,Level2,Level2};
        
        %% Linear Mixed Model Fit Method
        %optoptions=optimoptions('fminunc','Display','iter','Algorithm','quasi-newton');
        optoptions=optimoptions('fminunc','Algorithm','quasi-newton');
        
        lme = fitlmematrix(X,Y,Z,GG,'FixedEffectPredictors',...
            {'d_0000' 'd_1000' 'd_0100' 'd_0010' 'd_0001' 'd_1100' 'd_1010' 'd_1001' 'd_0110' 'd_0101' 'd_0011' 'd_1110' 'd_1101' 'd_1011' 'd_0111' 'd_1111'},...
            'RandomEffectPredictors',...
            {{'v_000l'},{'v_100l'},{'v_010l'},{'v_001l'},{'v_110l'},{'v_101l'},{'v_011l'},{'v_111l'},{'u_00kl'},{'u_01kl'},{'u_10kl'},{'u_11kl'},{'r_0j'},{'r_1j'}},'RandomEffectGroups', {'Organization','Organization','Organization','Organization','Organization','Organization','Organization','Organization','Unit','Unit','Unit','Unit','Team','Team'},...
            'FitMethod','REML','Optimizer','fminunc','OptimizerOptions',optoptions);
        
        %             lme = fitlmematrix(X,Y,Z,GG,'FixedEffectPredictors',...
        %                 {'d_0000' 'd_1000' 'd_0100' 'd_0010' 'd_0001' 'd_1100' 'd_1010' 'd_1001' 'd_0110' 'd_0101' 'd_0011' 'd_1110' 'd_1101' 'd_1011' 'd_0111' 'd_1111'},...
        %                 'RandomEffectPredictors',...
        %                 {{'v_000l','v_100l','v_010l','v_001l','v_110l','v_101l','v_011l','v_111l'},{'u_00kl','u_01kl','u_10kl','u_11kl'},{'r_0j','r_1j'}},'RandomEffectGroups', {'Organization' 'Unit' 'Team'},...
        %                 'FitMethod','REML','Optimizer','fminunc','OptimizerOptions',optoptions);
        
        
        MC_history((i-1)*MCtrials +MCk,:)=[n g h f nthroot(nghf_reml,4) Case scen ES e_var ...
            lme.Coefficients.Estimate' lme.Coefficients.SE' ...
            lme.Coefficients.tStat' lme.Coefficients.pValue' sqrt(cell2mat(lme.covarianceParameters))' sqrt(lme.MSE)];
        
        %            MC_history=[MC_history; n g h ngh_reml Case scen ES e_var ...
        %                lme.Coefficients.Estimate' lme.Coefficients.SE' ...
        %                lme.Coefficients.tStat' lme.Coefficients.pValue' lme.mse];
        toc
    end
end


csvwrite('Level4_RO2_REML30BlockPart2_14Nov2018.csv', MC_history);


