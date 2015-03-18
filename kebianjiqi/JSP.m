
function [MinVal,P] = JSP(T,Jm,JmNumber,NIND,MAXGEN,GGAP,XOVR,MUTR)

%  NIND=40;        %������Ŀ(Number of individuals)
%  MAXGEN=50;      %����Ŵ�����(Maximum number of generations)
%  GGAP=0.9;       %����(Generation gap)
%  XOVR=0.8;       %������
%  MUTR=0.6;       %������
 tic
 gen=0;  %��������
 
 [PNumber MNumber]=size(Jm);  %  PNumber ��������   ������� 
 trace=zeros(2, MAXGEN);      %Ѱ�Ž���ĳ�ʼֵ
 WNumber=PNumber*MNumber;     %�����ܸ���

 %��ʼ��
 Number=zeros(1,PNumber);
 for i=1:PNumber
   Number(i)=MNumber;
 end

% ��ʼ��Ⱥ ���� 2����� ��һ������ʾ����  ��2�� ���� ��ʾ����
Chrom=zeros(NIND,2*WNumber);
for j=1:NIND
    
   WPNumberTemp=Number;    
   for i=1:WNumber    
       
     %������ɹ���
     val=unidrnd(PNumber); 
     while WPNumberTemp(val)==0     
         val=unidrnd(PNumber);
     end 
       
       %��һ������ʾ����
       Chrom(j,i)= val;          
       WPNumberTemp(val)=WPNumberTemp(val)-1;
       
       %��2�� ���� ��ʾ����
       Temp=Jm{val,MNumber-WPNumberTemp(val)};
       SizeTemp=length(Temp);
       %������ɹ������
       Chrom(j,i+WNumber)= unidrnd(SizeTemp);
   end       
end

%     Chrom=Chrom;
 
 %����Ŀ�꺯��ֵ
   [PVal ObjV P S]=cal(Chrom,JmNumber,T,Jm);  

while gen<MAXGEN
    FitnV=ranking(ObjV);                                 %������Ӧ��ֵ(Assign fitness values)         
     SelCh=select('rws', Chrom, FitnV, GGAP);               %ѡ��
%   
    SelCh=across(SelCh,XOVR,Jm,T);            %���湤�򣬻������仯 
  %   SelCh=aberrance(SelCh,MUTR,Jm,T);        %���칤�򣬻������仯 
     SelCh=aberranceJm(SelCh,MUTR,Jm,T);                %���� ����    
    
    [PVal ObjVSel P S]=cal(SelCh,JmNumber,T,Jm);   %����Ŀ�꺯��ֵ
    [Chrom ObjV] =reins(Chrom, SelCh,1, 1, ObjV, ObjVSel);                %�ز����Ӵ�������Ⱥ
    
    gen=gen+1;                                             %������������
    
    %������Ž⼰����ţ�����Ŀ�꺯��ͼ���б����YΪ���Ž�,IΪ��Ⱥ�����
    trace(1, gen)=min(ObjV);                               %�Ŵ��㷨���ܸ���
    trace(2, gen)=mean(ObjV);  
    
    % ��ʼ��
    if gen==1
        Val1=PVal;
        Val2=P;
        MinVal=min(ObjV);%��Сʱ��
        STemp=S;
    end    
    %��¼ ��С�Ĺ���
    if MinVal> trace(1,gen)
        Val1=PVal;
        Val2=P;
        MinVal=trace(1,gen);
        STemp=S;
    end   
end

PVal=Val1; %����ʱ��
P=Val2;  %���� 
S=STemp; %���Ȼ��򺬻�������

%�����ı仯
hold on;
plot([0,0],[0,0]);
plot(trace(1,:));
hold on;
plot(trace(2,:),'-.');grid;
legend('��ı仯','��Ⱥ��ֵ�ı仯');  
%��ʾ���
figure(2);

MP=S(1,PNumber*MNumber+1:PNumber*MNumber*2); 

for i=1:WNumber
    
    val= P(1,i);
    a=(mod(val,100)); %����
    b=((val-a)/100); %����
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


