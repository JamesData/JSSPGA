%% [MinVal,P] = JSP(T,Jm,NIND,MAXGEN,GGAP,XOVR,MUTR)
% 参数说明：
% T          - 各工件各工序使用的时间
% Jm         - 各工件各工序使用的机器
% NIND       - 个体数目(Number of individuals)
% MAXGEN     - 最大遗传代数(Maximum number of generations)
% GGAP       - 代沟(Generation gap)
% XOVR       - 交叉率
% MUTR       - 变异率
%%

function [MinVal,P] = JSP(T,Jm,NIND,MAXGEN,GGAP,XOVR,MUTR)
% 比如 MinVal = JSP(T,Jm,40,50,0.9,0.8,0.6)
% NIND=60;        个体数目(Number of individuals)
% MAXGEN=500;     最大遗传代数(Maximum number of generations)
% GGAP=0.9;       代沟(Generation gap)
% XOVR=0.8;       交叉率
% MUTR=0.6;       变异率



[PNumber,MNumber]=size(T);              % PNumber工件数量   % MNumber单个工件工序的数量
gen=0;                                  % 代计数器
trace=zeros(2, MAXGEN);                 % 寻优结果的初始值

% 初始化群
% 基于调度优先级的编码方法：每个基因对应一道工序，代表该工序在进行调度操作时的处理优先级。
% 例如一个3!3问题，则染色体用1-3的整数表示为： [1 3 2 1 1 3 3 2 2]
% 基因与工序的对应关系如下：
%  P11  P31  P21  P12  P13  P32  P33  P22  P23
%  [1    3    2    1    1    3    3    2    2 ]
% 其中Pij表示第i个零件的第j道工序
WNumber=PNumber*MNumber;                % 工序总数量
Chrom=zeros(NIND,WNumber);              % 个体集合矩阵
Number=zeros(1,PNumber);
for i=1:PNumber
    Number(i)=MNumber;
end

% 随机生成工序序列
for j=1:NIND
    WPNumberTemp=Number;
    for i=1:WNumber
        val=unidrnd(PNumber);          % unidrnd: Discrete uniform random numbers，one to N
        
        while WPNumberTemp(val)==0     % 判断该工件是否编码完
            val=unidrnd(PNumber);
        end
        
        Chrom(j,i)= val;
        WPNumberTemp(val)=WPNumberTemp(val)-1;
    end
end

% 计算目标函数值
[PVal, ObjV, P]=cal(Chrom,T,Jm);
% PVal      - 为最佳调度工序时间
% ObjV      - 为群中每个个体的调度工序时间
% P         - 为最佳输出的调度工序

while gen<MAXGEN
    FitnV=ranking(ObjV);                                 % 分配适应度值(Assign fitness values)
    SelCh=select('sus', Chrom, FitnV, GGAP);             % 选择
    SelCh=across(SelCh,XOVR,PNumber);                    % 交叉
    SelCh=aberrance(SelCh,MUTR);                         % 变异
    [PVal, ObjVSel, P]=cal(SelCh,T,Jm);                  % 计算目标函数值
    [Chrom, ObjV] =reins(Chrom, SelCh,1, 1, ObjV, ObjVSel);            % 重插入子代的新种群
    gen=gen+1;                                           % 代计数器增加
    
    % 输出最优解及其序号，并在目标函数图像中标出，Y为最优解,I为种群的序号
    % 遗传算法性能跟踪
    trace(1, gen)=min(ObjV);
    trace(2, gen)=mean(ObjV);
    
    % 初始化
    if gen==1
        Val1=PVal;
        Val2=P;
        MinVal=min(ObjV);                  % 初始化最小调度工序时间
    end
    % 记录最佳调度工序时间、最佳输出的调度工序、最小的调度工序时间
    if MinVal> trace(1,gen)
        Val1=PVal;
        Val2=P;
        MinVal=trace(1,gen);
    end
end

PVal=Val1;                  % 工序时间
P=Val2;                     % 工序

% 计算解的变化
figure(1);
hold on;
plot([0,0],[0,0]);
BestV=plot(trace(1,:),'r');
hold on;
MeanV=plot(trace(2,:),'b');
grid;
legend([BestV,MeanV],'解的变化','种群均值的变化');
% 显示结果
figure(2);
for i=1:WNumber
    val= P(1,i);
    a=(mod(val,10))+1;                             % Modulus after division
    b=((val-a+1)/10);
    mText=Jm(b,a);
    PlotRec(PVal(1,i),PVal(2,i),mText);
    hold on;
    mPoint1=PVal(1,i);
    mPoint2=PVal(2,i);
    x1=mPoint1;
    y1=mText-0.5;
    x2=mPoint2;
    y2=mText-0.5;
    x3=mPoint2;
    y3=mText;
    x4=mPoint1;
    y4=mText;
    fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,0.5,1]);
    word=num2str(P(1,i));
    %text(0.5*mPoint1+0.5*mPoint2,mText-0.5,word);
    text(mPoint1,mText-0.7,word);
end
