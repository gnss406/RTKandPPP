%*****************************************************************
% Kalman滤波
function [xe, v,  K,Pv]=akf_m(d,x,A,H,G,P,Q ,R)
% d     观测数据
% x     初值
% A H G 系统矩阵
% P     初值的方差
% Q R   噪声协方差矩阵
% flag  1 平滑  0 只进行Kalman滤波   缺省值：为1

m=length(x);

if( nargin < 9)  %进行平滑
    flag=1;
end
if( nargin < 10)  %进行平滑
    B=zeros(size(x,1),size(x,1));
    u=zeros(size(x,1),size(d,2));
end

nn=size(d,1);


x0=x;
%normal kalman filtering
for k=1:size(d,2)

    %Time update
    %predict state
    x=A*x+B*u(:,k);
    %predict covariance
    P=A*P*A'+G*Q*G';
    xminus(:,k)=x;
    Pminus(:,:,k)=P;   
    %compute Kalman gain
    Pv=H*P*H'+R;
    K=P*H'/(H*P*H'+R);
    %  P3=P;
    %update state 
    v(:,k)=d(:,k)-H*x;   %新息%
    vv=v(:,k);
    md(k)=vv'*Pv*vv;
    
    K=P*H'/(H*P*H'+R);
    x=x+K*(d(:,k)-H*x);
      
    %update covariance
    %P=(eye(m)-K*H)*P;
    P=(eye(m)-K*H)*P*(eye(m)-K*H)'+K*R*K';

    Pplus(:,:,k)=P;
    xplus(:,k)=x;

    xe(:,k)=x;
end

x=x0;
% compute 极限误差
xalpha=(median(abs(v))/0.6745*3)^2;

%akf
for k=1:size(d,2)

    %Time update
    %predict state
    x=A*x+B*u(:,k);
    
    lamda=1;
    pv=lamda*H*A*P*A'*H'+H*G*Q*G'*H'+R;
    v(:,k)=d(:,k)-H*x;   %新息%
    vv=v(:,k);
    md(k)=vv'*inv(Pv)*vv;%
    %迭代次数
    nt=0;
    P0=P;
    P=A*P*A'+G*Q*G';
    if k==2937
        k
    end
    % 参考An adaptive fading Kalman filter based on Mahalanobis distance
    while md(k)>xalpha&&nt<200
        gamma=md(k);
        
        p_minus=lamda*P;
%         lamda1=lamda+(gamma-xalpha)/(vv'*inv(pv)*(H*p_minus*H')*inv(pv)*vv);
        lamda1=lamda+(gamma-xalpha)/(vv'*inv(pv)*(H*P*H')*inv(pv)*vv);
        lamda=lamda1;
        pv=H*p_minus*H'+R;
        md(k)=vv'*inv(pv)*vv;
        nt=nt+1;
    end
    P=P0;

    %predict covariance
    P=lamda*A*P*A'+G*Q*G';

    %compute Kalman gain
    K=P*H'/(H*P*H'+R);

    x=x+K*(d(:,k)-H*x);
      
    %update covariance
    P=(eye(m)-K*H)*P*(eye(m)-K*H)'+K*R*K';

    xe(:,k)=x;
end






end
% %********************************************
% function [M,Pyk]=MSHjuli(y,y1)
% 
% %马氏距离作为判别量
% %y为观测量；y1为新息；
% %Pyk为关于新息的概率密度函数；
% %Pyk1为测量协方差矩阵；
% Pyk=(exp(-1/2*(y-y1)'*inv(Pyk1)*(y-y1)))/sqrt((2*pi)^2*abs(Pyk1));
% 
% M=(y-y1)'*inv(Pyk1)*(y-y1);%求出每个y到y1的马氏距离
% end
% 
% %********************************************
% function [  ]=KaFjianyan(y)
% %卡方检验的作为检验量
% %以下部分进总体均值未知时,方差的假设检验中的右边检验
% n=length(y(:)); %样本容量n
% alf=0.05;        %自由度
% c=32.9;         %从分布表查的左端值
% d=8.91;         %从分布表查的右端值
% 
% y1=means(y);
% 
% for i=1:n-10
% x(i,:)=y(i:i+9);
% 
% x1(i)=means(x(i,:));  
% end
% for i=1:n-1
%     b(i)=(y(i)-y1)^2+(y(i+1)-y1)^2;
% end
% for t=1:10
%     e(i)=(x(i)-x1(i))^2+(x(i+1)-x1(i))^2;
% end
% 
%     rmsv=sqrt(b(end)/n-1);   %子样中误差
%     
%     KF=sqrt(b(end))/rmsv;  %卡方分布统计量
%    if (c<KF&&KF<d)
%        [xe v s P3 K PP]=kalmansmoother(d,x,A,H,G,P,Q,R,flag,B,u)
%    else 
%    end
% end
%   function [chi2 ]=KaFjianyan1(y)
% Alf=0.05;       %置信区间95%
% n=length(y(:)); %样本容量n
% df=n-1;         %自由度df
% y1=means(y);
% for i=1:n-1
%     b(i)=(y(i)-y1)^2+(y(i+1)-y1)^2;
% end
%  rmsv=b(end)/n-1;   %子样中误差
%  
% Py(i)=pdf(norm,y(i),yi,rmsv); %在y的地i个数的概率密度函数；
% chi2=chi2inv(df,Py);  %利用卡方分布函数求出接受域的右边界
% end
% %*********************************************
% %t检验进行检验
% function [h,sig,ci3]=Tjianyan(Y)
% Alf=0.05;  %显著性水平
% u=mean(Y); %求均值
% [h, sig, ci3]=ttest(Y, u, Alf, tail);%对μ进行T检验
% %h=1,显著性水平下拒绝原假设；h=0,显著性水平下接受原假设；
% % sig 表示在 X 的均值等于 u 的原假设下较大或者统计意义下较大的概率值 
% %ci3表示返回一个置信度为 100(1-Alf)％的均值的置信区间
% %tail=0, 表示备择假设:μ≠m (默认 , 双边检验 );
% %tail=1, 表示备择假设: μ>m(右边检验 );
% %tail=-1, 表示备择假设 :μ<m(左边检验 ).
% end
% %*********************************************
% %牛顿迭代法
% 
% %*******************************************
