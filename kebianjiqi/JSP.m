
function [MinVal,P] = JSP(T,Jm,JmNumber,NIND,MAXGEN,GGAP,XOVR,MUTR)

%  NIND=40;        %个体数目(Number of individuals)
%  MAXGEN=50;      %最大遗传代数(Maximum number of generations)
%  GGAP=0.9;       %代沟(Generation gap)
%  XOVR=0.8;       %交叉率
%  MUTR=0.6;       %变异率
 tic
 gen=0;  %代计数器
 
 [PNumber MNumber]=size(Jm);  %  PNumber 工件个数   工序个数 
 trace=zeros(2, MAXGEN);      %寻优结果的初始值
 WNumber=PNumber*MNumber;     %工序总个数

 %初始化
 Number=zeros(1,PNumber);
 for i=1:PNumber
   Number(i)=MNumber;
 end

% 初始化群 生成 2层代码 第一层代码表示工序  第2层 代码 表示机器
Chrom=zeros(NIND,2*WNumber);
for j=1:NIND
    
   WPNumberTemp=Number;    
   for i=1:WNumber    
       
     %随机产成工序
     val=unidrnd(PNumber); 
     while WPNumberTemp(val)==0     
         val=unidrnd(PNumber);
     end 
       
       %第一层代码表示工序
       Chrom(j,i)= val;          
       WPNumberTemp(val)=WPNumberTemp(val)-1;
       
       %第2层 代码 表示机器
       Temp=Jm{val,MNumber-WPNumberTemp(val)};
       SizeTemp=length(Temp);
       %随机产成工序机器
       Chrom(j,i+WNumber)= unidrnd(SizeTemp);
   end       
end

%     Chrom=Chrom;
 
 %计算目标函数值
   [PVal ObjV P S]=cal(Chrom,JmNumber,T,Jm);  

while gen<MAXGEN
    FitnV=ranking(ObjV);                                 %分配适应度值(Assign fitness values)         
     SelCh=select('rws', Chrom, FitnV, GGAP);               %选择
%   
    SelCh=across(SelCh,XOVR,Jm,T);            %交叉工序，机器不变化 
  %   SelCh=aberrance(SelCh,MUTR,Jm,T);        %变异工序，机器不变化 
     SelCh=aberranceJm(SelCh,MUTR,Jm,T);                %变异 机器    
    
    [PVal ObjVSel P S]=cal(SelCh,JmNumber,T,Jm);   %计算目标函数值
    [Chrom ObjV] =reins(Chrom, SelCh,1, 1, ObjV, ObjVSel);                %重插入子代的新种群
    
    gen=gen+1;                                             %代计数器增加
    
    %输出最优解及其序号，并在目标函数图像中标出，Y为最优解,I为种群的序号
    trace(1, gen)=min(ObjV);                               %遗传算法性能跟踪
    trace(2, gen)=mean(ObjV);  
    
    % 初始化
    if gen==1
        Val1=PVal;
        Val2=P;
        MinVal=min(ObjV);%最小时间
        STemp=S;
    end    
    %记录 最小的工序
    if MinVal> trace(1,gen)
        Val1=PVal;
        Val2=P;
        MinVal=trace(1,gen);
        STemp=S;
    end   
end

PVal=Val1; %工序时间
P=Val2;  %工序 
S=STemp; %调度基因含机器基因

%计算解的变化
hold on;
plot([0,0],[0,0]);
plot(trace(1,:));
hold on;
plot(trace(2,:),'-.');grid;
legend('解的变化','种群均值的变化');  
%显示结果
figure(2);

MP=S(1,PNumber*MNumber+1:PNumber*MNumber*2); 

for i=1:WNumber
    
    val= P(1,i);
    a=(mod(val,100)); %工序
    b=((val-a)/100); %工件
    Temp=Jm{b,a};
    mText=Temp(MP(1,i));
    
    x1=PVal(1,i);
    x2=PVal(2,i);
    
    y1=mText-1;
    y2=mText;
    PlotRec(x1,x2,mText);
    
    PlotRec(PVal(1,i),PVal(2,i),mText);
    hold on;
    
     fill([x1,x2,x2,x1],[y1,y1,y2,y2],[1-1/b,1/b,b/PNumber]);
    text((x1+x2)/2,mText-0.25,num2str(P(i)));  
     
end 
toc


