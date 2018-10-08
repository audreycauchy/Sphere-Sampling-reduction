#include <NTL/matrix.h>
#include <NTL/vec_vec_RR.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <fstream>
#include <iostream>
using namespace std;
using namespace NTL;
void GSO(mat_RR A);
mat_RR AO;
mat_RR U;
vec_RR DIS;
mat_ZZ A11;
int insert_sign;

int main()
{mat_RR B;
 int isre;

 vec_RR list;//build a list of normal distribute
 list.SetLength(60001);
 ifstream min("normallist1.txt");
 min>>list;


 int left;
 int right;
 int mid;
 int mm;
 vec_RR memory;
 vec_RR memory1;
 vec_RR proj_mem;
 int count;
 int new1;
 RR short_length;
 mat_RR mat0;
 RR pow_a;//随机数生成底数
 RR pow_b;//随机数生成指数
 RR pow_c;
 RR randmax;
 RR mul;
 RR detL;
 RR ghlbound;
 RR ghlbound1;
 mul=0.48;
 int cm;

	

 vec_RR ini1;
 int num0;
 int nn;
 int cnum;
 int rounds;
 int num2;
 int num3;
 int num4;
 int num5;
 RR sum_S;
 int insertnum;
 int insertnum1;
 RR ran_length;
 vec_RR sum_GSO;
 vec_RR pro_coe;
 vec_RR cond;
 RR  pownum;
 RR ran_pass;
 ZZ detA;
 mat_ZZ A1;
 mat_ZZ A;
 vec_RR con_length;
 vec_RR con;
 vec_RR short_vec;
 vec_RR GSOd_pow;//正交基长度的平方
 vec_RR proj_vec;
 vec_RR ini;
 RR g_length1;
 int roundnormal;

 int out_sign;
 out_sign=0;
 pow_a=2;
 pow_b=53;
 count=1;
 pow_c=-1;
 detL=1;

 ifstream fin("svpchallengedim100seed0.txt");//文件读入基
 fin>>A;
 int on;
 on=0;
 int nc=A.NumCols();
 int nr=A.NumRows();
 AO.SetDims(nr,nc);
 U.SetDims(nr,nc);
 DIS.SetLength(nc);
 ini.SetLength(nc);
 ini1.SetLength(nc+1);
 memory.SetLength(nc);
 //proj_sum.SetLength(nc);
 con_length.SetLength(nr);
 //pre_vec.SetLength(nc);
 A1.SetDims(nr+1,nc);
 B.SetDims(nc,nr);
 con.SetLength(nc);
 pro_coe.SetLength(nc);
 sum_GSO.SetLength(nc);
 GSOd_pow.SetLength(nc);
 proj_vec.SetLength(nc+1);
 A11.SetDims(nc+1,nc);
 memory1.SetLength(nc);
 //random_b=conv<long>(pow_b);
 randmax=pow(pow_a,pow_b)-1;
 for(new1=1;new1<=nr;new1++)//初始化空向量用于更新
	{ini(new1)=0;
         ini1(new1)=0;

	}
	ini1(nr+1)=0;
 short_vec=ini;

 
 double delta=0.99;//BKZ参数
  long BlockSize=20;//BKZ参数
  long prune=10;//BKZ参数
  LLLCheckFct check = 0; //BKZ参数
  long verbose = 0;//BKZ参数
LLL(detA,A,0);
for(BlockSize=3;BlockSize<=40;BlockSize++)
{
 BKZ_FP(A,delta,BlockSize,prune,check,verbose);//BKZ执行 约化基
}

for(num0=1;num0<=nc;num0++)
	{B(num0)=conv<vec_RR>(A(num0));
	}

	GSO(B);
pownum=1.0/nr;

for(num5=1;num5<=nc;num5++)
{
detL=detL*DIS(num5);
}
ghlbound=1.05*pow(conv<RR>(tgamma(conv<double>(nc*0.5+1))),pownum)/sqrt(ComputePi_RR())*pow(detL,pownum);
ghlbound1=ghlbound*ghlbound;
cout<<ghlbound1<<endl;



//insert_sign=nc;
for(nn=1;nn<=nc;nn++)//求正交基长度的平方
{GSOd_pow(nn)=DIS(nn)*DIS(nn);}
proj_mem=GSOd_pow;
short_length=GSOd_pow(1)-1;
g_length1=DIS(1);//可调参数 控制半径

for(cnum=1;cnum<=nc;cnum++ )//避免重复计算 从单位圆采样还原到以第一个正交基为半径的圆上
	{pro_coe(cnum)=g_length1/DIS(cnum);}

int a=1;


while((out_sign>60)||(out_sign==0))
{



//while(insert_sign==nc)
//for(rounds=1;rounds<=10000;rounds++)//采出100000000个

ran_length=pow(conv<RR>(RandomBits_ZZ(53))/randmax,pownum);
for(roundnormal=1;roundnormal<=nc;roundnormal++)
	{	
		ran_pass=(conv<RR>(RandomBits_ZZ(53)))/randmax;
	    right=60001;
	   	left=1;
         mid=floor((right+left)*0.5);
       

	       		while(mm!=1)
	       {
	       	if(list(mid)>=ran_pass)
	           {
	           	if(mid-left==1)
	               {mm=1;
	               	
	               	con_length(roundnormal)=((left+mid)*0.5-1)*0.0001-3;
	               	
	               }
	             else
	                 {right=mid;
	                 	
	                  mid=floor((left+right)*0.5);}}
	       else
	       	{
	       		if(right-mid==1)
	       	    {mm=1;
	       	 
	       	    	con_length(roundnormal)=((mid+right)*0.5-1)*0.0001-3;
	       	    	
	       	    }
	       	 else
	                 {left=mid;
	                 	
	                 mid=floor((left+right)*0.5);}}
	   }
	     
	    mm=0;}
	    sum_S=con_length*con_length;

//cout<<con_length<<"\n";
con_length=con_length*(1/sqrt(sum_S))*ran_length;//单位圆上采样乘上半径上的随机分布
		
//


for(num2=1;num2<=nc;num2++)//映射到原球面上采
	 	 	{	 	 con(num2)=con_length(num2)*pro_coe(num2)*mul;

	 	 con(num2)=round(con(num2));

	 	 
	 	                     }
	 	                     //cout<<con<<endl;

//cout<<con<<"\n";
if(con!=ini)//判断插入位置和还原在原基下的坐标
{
out_sign=0;
cond=con;
sum_GSO=ini;
proj_vec=ini1;


	for(num3=nr;num3>=1;num3--)
		{  
			for(num4=nr;num4>num3;num4--)
      	{
  sum_GSO(num3)=sum_GSO(num3)+cond(num4)*U(num4)(num3);
       }
       cond(num3)=cond(num3)-sum_GSO(num3);
       cond(num3)=ceil(cond(num3)-0.5);//原基下坐标
       sum_GSO(num3)=cond(num3)+sum_GSO(num3);
       proj_vec(num3)=sum_GSO(num3)*sum_GSO(num3)*GSOd_pow(num3)+proj_vec(num3+1);//投影长度



       	if((proj_vec(num3)<GSOd_pow(num3))&&(proj_vec(num3)!=0))//比该投影小 更新插入编号

       		{
       			
       	 out_sign=num3;

       	
       	
       	



       	 
       	 }
       }
       if((out_sign!=0)&&(cond*cond!=1))
       {
       	
       		
       		
       	 memory=cond*B;



     
      
       	
//insert_sign=out_sign

       	}
       	else
       		{out_sign=0;}
       	        
           }}
       	



      //cout<<out_sign<<endl;
     // cout<<memory<<endl;
//out_sign=0;	









           while(isre!=1)
           	{
           	  if((out_sign==1)&&(memory*memory<ghlbound1))
           		//if(insert_sign==1)
               {cout<<"short_vec="<<memory<<endl;
           		break;}
           cm=0;		

           count++;

           

           

for(insertnum=1;insertnum<out_sign;insertnum++)
	{
		A1(insertnum)=conv<vec_ZZ>(B(insertnum));
	}

A1(out_sign)=conv<vec_ZZ>(memory);
for(insertnum1=out_sign+1;insertnum1<=nc+1;insertnum1++)
    {A1(insertnum1)=conv<vec_ZZ>(B(insertnum1-1));}
//cout<<A1-A11<<endl;
//A11=A1;
//cout<<A1(insert_sign)<<endl;

//cout<<"B"<<endl;




LLL(detA,A1,0);
BlockSize=40;
 BKZ_FP(A1,delta,BlockSize,prune,check,verbose);//BKZ执行 约化基

for(num0=2;num0<=nc+1;num0++)
	{B(num0-1)=conv<vec_RR>(A1(num0));

		
	}

	GSO(B);
	if(count%50==0)
{cout<<"round"<<count<<endl;
	ofstream fout("100revbkz40.txt", ios::app|ios::ate);
           fout<<"round"<<count<<endl;
           fout<<"insert_sign="<<out_sign<<endl;
       		//fout<<"bound="<<proj_mem(insert_sign)<<endl;
       		//fout<<"pro_vec"<<proj_vec(insert_sign)<<endl;
           fout<<DIS<<endl;
           fout<<memory<<endl;
       		fout<<"mul"<<mul<<endl;}

g_length1=DIS(1);
for(cnum=1;cnum<=nc;cnum++ )//避免重复计算 从单位圆采样还原到以第一个正交基为半径的圆上
	{pro_coe(cnum)=g_length1/DIS(cnum);}

	for(nn=1;nn<=nc;nn++)//求正交基长度的平方
{GSOd_pow(nn)=DIS(nn)*DIS(nn);}

proj_mem=GSOd_pow;
         
         
         //fout<<"proj_mem="<<proj_mem<<"\n";	
insert_sign=nc;
out_sign=0;
memory=ini;
while((out_sign==0)||(out_sign>50))
//for(rounds=1;rounds<=10000;rounds++)//采出100000000个
{   cm=cm+1;
	if(cm>20000)
		{if(mul<0.58)
			{
			mul=mul+0.005;
			
			}
		 else
		 	{mul=0.48;}
			cm=0;
		if(insert_sign<70)
			{memory=memory1;
				out_sign=insert_sign;
				break;}
		}

ran_length=pow(conv<RR>(RandomBits_ZZ(53))/randmax,pownum);
for(roundnormal=1;roundnormal<=nc;roundnormal++)
	{	
		ran_pass=(conv<RR>(RandomBits_ZZ(53)))/randmax;
	    right=60001;
	   	left=1;
         mid=floor((right+left)*0.5);
       

	       		while(mm!=1)
	       {
	       	if(list(mid)>=ran_pass)
	           {
	           	if(mid-left==1)
	               {mm=1;
	               	
	               	con_length(roundnormal)=((left+mid)*0.5-1)*0.0001-3;
	               	
	               }
	             else
	                 {right=mid;
	                 	
	                  mid=floor((left+right)*0.5);}}
	       else
	       	{
	       		if(right-mid==1)
	       	    {mm=1;
	       	 
	       	    	con_length(roundnormal)=((mid+right)*0.5-1)*0.0001-3;
	       	    	
	       	    }
	       	 else
	                 {left=mid;
	                 	
	                 mid=floor((left+right)*0.5);}}
	   }
	     
	    mm=0;}
	    sum_S=con_length*con_length;
//cout<<con_length<<"\n";

con_length=con_length*(1/sqrt(sum_S))*ran_length;//单位圆上采样乘上半径上的随机分布
		
//cout<<con_length<<"\n";


for(num2=1;num2<=nc;num2++)//映射到原球面上采
	 	 	{	 	 con(num2)=con_length(num2)*pro_coe(num2)*mul;

	 	 con(num2)=round(con(num2));

	 	 
	 	                     }
	 	                     //cout<<con<<endl;

//cout<<con<<"\n";

if(con!=ini)//判断插入位置和还原在原基下的坐标
{
out_sign=0;
cond=con;
sum_GSO=ini;
proj_vec=ini1;


	for(num3=nr;num3>=1;num3--)
		{  
			for(num4=nr;num4>num3;num4--)
      	{
  sum_GSO(num3)=sum_GSO(num3)+cond(num4)*U(num4)(num3);
       }
       cond(num3)=cond(num3)-sum_GSO(num3);
       cond(num3)=ceil(cond(num3)-0.5);//原基下坐标
       sum_GSO(num3)=cond(num3)+sum_GSO(num3);
       proj_vec(num3)=sum_GSO(num3)*sum_GSO(num3)*GSOd_pow(num3)+proj_vec(num3+1);//投影长度



       	if((proj_vec(num3)<GSOd_pow(num3))&&(proj_vec(num3)!=0))//比该投影小 更新插入编号

       		{
       			
       	 out_sign=num3;

       	
       	
       	



       	 
       	 }
       }
       if((out_sign!=0)&&(cond*cond!=1))
       {
       	 
       		
       		
       	 memory=cond*B;


if(insert_sign>out_sign)

     {insert_sign=out_sign;
     	memory1=memory;
     }
      
       	
//insert_sign=out_sign

       	}
       	else
       		{out_sign=0;}
       	        
           }


       }
      //cout<<"insert_sign="<<insert_sign<<endl;
      //fout<<"cond="<<cond<<endl;
      //cout<<out_sign<<endl;

      //cout<<proj_vec(out_sign)<<"A"<<endl;
      //cout<<GSOd_pow(out_sign)<<"B"<<endl;
      //fout<<memory<<endl;
      //cout<<memory<<endl;

           	}
     







//}

}

void GSO(mat_RR A)
{//mat_RR U,AO;
 //vec_RR DIS;

 int i,ii,j,ar,ac,m;
 RR x,y,sqrt;
 sqrt=0.5;

 
  ar=A.NumRows();
 ac=A.NumCols();
AO=A;
 for(j=2;j<=ar;j++)
	 for(ii=1;ii<=(j-1);ii++)
		{   //cout<<ii<<"\n";
			InnerProduct(x,A(j),AO(ii));
 			InnerProduct(y,AO(ii),AO(ii));
					U(j)(ii)=x/y;
                                       AO(j)=AO(j)-U(j)(ii)*AO(ii);

}
for(m=1;m<=ar;m++)
{U(m)(m)=1;
}



for(i=1;i<=ac;i++)//更新施密特正交化长度
{InnerProduct(DIS(i),AO(i),AO(i));
pow(DIS(i),DIS(i),sqrt);}


//return GSO1;




}
