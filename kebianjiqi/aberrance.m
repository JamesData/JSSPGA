function ChromNew=aberrance(Chrom,MUTR,Jm,T)


%  Chrom=[3,2,2,3,1,1,2,2];
%  MUTR=0.4;

%��ʼ��
[NIND,WNumber]=size(Chrom);
WNumber=WNumber/2;

ChromNew=Chrom;

[PNumber MNumber]=size(Jm);
Number=zeros(1,PNumber);
for i=1:PNumber
  Number(i)=1;
end

for i=1:NIND    
    %�Ƿ����
    if MUTR>rand;
        
        %����λ��
        Pos1=unidrnd(WNumber);
        Pos2=unidrnd(WNumber);
        
        %����λ�ò���ͬ
        while Pos1==Pos2      
            Pos2=unidrnd(WNumber);
        end 
        
        %ȡһ������
        S=Chrom(i,:); 
        
        %����
        temp=S(Pos1);
        S(Pos1)=S(Pos2);
        S(Pos2)=temp;
           
%         % ���潻��ǰ�Ļ��� ����
%         STemp=Chrom(i,:);       
%         for k=1:WNumber            
%             Pos1=Find(S(k),STemp);           
%             S(WNumber+k)=STemp(WNumber+Pos1);
%             STemp(Pos1)=0;            
%         end
          
        temp=S(Pos1+WNumber);
        S(Pos1+WNumber)=S(Pos2+WNumber);
        S(Pos2+WNumber)=temp;

        WPNumberTemp=Number; 
        for j=1:WNumber
           
          JMTemp=Jm{S(j), WPNumberTemp(S(j))};
          SizeTemp=length(JMTemp);
                                  
         %ѡ������� �ӹ�ʱ���ٵ�ѡ���ʴ�
         if SizeTemp<S(j+WNumber)      
                 S(j+WNumber)=selectJm(S(j++WNumber),T{S(j),WPNumberTemp(S(j))});
         end
           
            WPNumberTemp(S(j))=WPNumberTemp(S(j))+1;
        end   
        

        %���ݷ�����Ⱥ
        ChromNew(i,:)=S;
    end
end




 