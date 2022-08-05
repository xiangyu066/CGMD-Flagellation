%**************************************************************************
%       FlagellaGrowthSimulation_MD_V13.m
%      ===================================
% CREATED: December 28, 2016
% MODIFIED: April 7, 2017
% EDITOR: 嚇死你的眉毛(X.Y. ZHUANG)
%
% DESCRITION:
%   <...>
%
% NOTES:
%   - 
%   - MD計算只計算前後一個, 因此使用矩陣向量化使matlab運算加速
%   - 將動畫向量化
%
% COMPILE:
%   - WIN7 SP1(64-bit)
%   - MATLAB R2013a(VER.8.1)
%**************************************************************************
close all, clear all
tic

%% set paramters of system
KT = 4.114; % [pN*nm]
site_step = 0.47; % [nm]
sigma = 56; % diameter, [nm]
epsilon = 1*KT; 
dt = 1e-4;
DIFF_Const = 1e+3; % [nm^2/s]
F_loading = 3; % pN
drag_coef = KT/DIFF_Const;
T_loading = 0.14; % generation time for one flagellin
N = 200; % 預設通道內最後會有的flagellin數量

%% initialization
output = ['V13_D',num2str(DIFF_Const),'_F',num2str(F_loading),'_E',num2str(epsilon/KT),'_T',num2str(T_loading),'_dt',num2str(dt),'.txt'];
site = 0; % boundary of flagella
monitor_timestep = 0.5;

x = zeros(1,N);
Vx = zeros(1,N);
flag_flagellin = zeros(1,N); % 是否有裝載的flagellin
Fx_loading = zeros(1,N); % pump施力

flag_pump = 0; % 監控pump上是否有flagellin
flag_pump2 = 0 ;% 監控固定loading time抓Flagellin

theta = linspace(0,2*pi,100); % 用來畫圓的

%%
tcount = 0; % 計算目前是第幾個timestep
ncount = 0; % 計算目前有幾個flagellin在通道內

fid = fopen(output,'w');
fprintf(fid,'%.1f\t%.1f\t%d\n',0.0,0.0,0);

while (site < 8100)
    
    tcount = tcount+1;

    %%%%%%% adding new partical into tube %%%%%%%
    if (rem(tcount,round(T_loading/dt))==0 && flag_pump2 ==0)
        flag_pump2 = 1;
    end
    
    if (flag_pump ==0 && tcount==1)
        ncount = ncount+1;
        
        x(ncount) = -sigma/2;
        Vx(ncount) = 0;
        flag_flagellin(ncount) = 1;
        Fx_loading(ncount) = F_loading;
        
        flag_pump = 1;
    elseif (flag_pump == 0 && flag_pump2==1)
        ncount = ncount+1;
        
        x(ncount) = -sigma/2;
        Vx(ncount) = 0;
        flag_flagellin(ncount) = 1;
        Fx_loading(ncount) = F_loading;
        
        flag_pump = 1;
        flag_pump2 = 0;
    end
    
    if (ncount>0) % 若通道內有flagellin才計算模型
        
        %%%%%%% calculate force interaction (MD) %%%%%%%
        dx1 = x(flag_flagellin>0)-circshift(x(flag_flagellin>0),[0,1]);
        dx2 = x(flag_flagellin>0)-circshift(x(flag_flagellin>0),[0,-1]);
        
        if (length(dx1)>1)
            rij1 = sqrt(dx1.^2);
            rij2 = sqrt(dx2.^2);
            
            n1 = (dx1>0)-(dx1<0);
            n2 = (dx2>0)-(dx2<0);
            
            Fx = n1.*-4*epsilon.*((6*sigma^6)./rij1.^7 - (12*sigma^12)./rij1.^13)...
                +n2.*-4*epsilon.*((6*sigma^6)./rij2.^7 - (12*sigma^12)./rij2.^13)...
                +Fx_loading(flag_flagellin>0);
        else
            Fx = Fx_loading(flag_flagellin>0);
        end
        
        % update next position and velocity
        Vx(flag_flagellin>0) = Fx/drag_coef;
        if (Vx(ncount)<0 && flag_pump ==1) % 避免在pump上的flagellin倒退
            Vx(ncount) = 0;
        end
        
        delta_x = sqrt(2*DIFF_Const*dt)*randn(1,N); % generate random step
        if (flag_pump ==1) % 在pump上的flagellin不受熱擾動
            delta_x(ncount) = 0;
        end
         
        x(flag_flagellin>0) = x(flag_flagellin>0)+delta_x(flag_flagellin>0)+Vx(flag_flagellin>0).*dt;
        
        % 若不是在pump上的flagellin要倒退回通道, 則強制停在sigma/2的地方
        bw = (x(flag_flagellin>0)>sigma/2);
        if (flag_pump == 1) 
            bw(ncount) = true;
        end
        
        x(flag_flagellin>0) = x(flag_flagellin>0).*bw + (sigma/2)*~bw;
        
        % 檢查是否有可以變成軌道的粒子並重置矩陣
        Num_flagellin2site = sum(x(flag_flagellin>0)>site-sigma/2);

        if (Num_flagellin2site>0)
            site = site+site_step*Num_flagellin2site;
            x = [x(Num_flagellin2site+1:end),zeros(1,Num_flagellin2site)];
            Vx = [Vx(Num_flagellin2site+1:end),zeros(1,Num_flagellin2site)];
            Fx_loading = [Fx_loading(Num_flagellin2site+1:end),zeros(1,Num_flagellin2site)];
            flag_flagellin = [flag_flagellin(Num_flagellin2site+1:end),zeros(1,Num_flagellin2site)];
            ncount =  ncount-Num_flagellin2site;
        end
    end
    
    % 判斷脫離pump上是否還有flagellin
    if (ncount==0)
        flag_pump = 0;
    elseif(x(ncount)>sigma/2 && flag_pump ==1)
        flag_pump = 0;
        Fx_loading(ncount) = 0;
    end
    
    %%%%%%% output file and display simulation %%%%%%%
    if (rem(tcount,round(monitor_timestep/dt))==0)
        
%         A = repmat(sigma/2*cos(theta'),[1,sum(flag_flagellin>0)]);
%         A2 = repmat(sigma/2*sin(theta'),[1,sum(flag_flagellin>0)]);
%         B = repmat(x(flag_flagellin>0),[length(theta),1]);
%         xx = A+B;
%         yy = A2+sigma/2;
%         
%         if (flag_pump==1)
%             fill(xx(:,end),yy(:,end),'r',xx(:,1:end-1),yy(:,1:end-1),'b');
%         else
%             fill(xx,yy,'b');
%         end
%         
%         title(num2str((tcount)*dt))
%         set(gca,'DataAspectRatio',[1 1 1]); % set aspect ratio x:y to 1:1
%         axis([-0.01,site,0,sigma]);
%         drawnow, clf(gcf);
        
        fprintf(fid,'%.1f\t%.1f\t%d\n',tcount*dt,site,ncount);
        disp(['T = ',num2str(tcount*dt), '[s] / ','flagella = ',num2str(site),' [nm] / ncount = ', num2str(ncount), ' / flag = ', num2str(flag_pump),' / flag2 = ',num2str(flag_pump2)])
    end
    
end
fclose(fid);

%%
toc

