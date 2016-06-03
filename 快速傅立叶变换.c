/*这个程序最多512个点，如果需要处理更多的数据，可以改动R的值*/
#include   <stdio.h>
#include   <math.h>
#define   R   512
#define   PI   3.141592653
typedef struct//定义一个复数结构
     {
      double   real;//由于三角运算结果是double型，所以不能用float定义
      double   img;
      }complex;
complex   x[R], W[R];//FFT输入和输出的存储空间
void   fft(); //快速傅里叶变换算法实现
void   change();

int   N=0;//定义并初始化序列的长度，必须是2的自然数次幂，不能超过N
int   L=0;
void   add(complex   a,complex   b,complex   *c)   //复数加法
           {   
             c->real=a.real+b.real;   
             c->img=a.img+b.img;   
            }   
    
 void   mul(complex   a,complex   b,complex   *c)   //复数乘法
             {   
               c->real=a.real*b.real   -   a.img*b.img;   
               c->img=a.real*b.img   +   a.img*b.real;   
              }   
 void   sub(complex   a,complex   b,complex   *c)   //复数减法
              {   
                c->real=a.real-b.real;   
                c->img=a.img-b.img;   
               }   
int   main()
        {
         int   i;
         printf("序列长度为2的L次方，请输入L的值（1~9）:\n");//输入指数值（1~9）确定序列的长度
         scanf("%d",&L);
		 N=(int) pow(2,L);
			 if(L<1)
				 printf("L值输入错误");
			 else if(N>R)
				 printf("L值过大");
			 else												//只有在L值合法情况才执行主程序
				{
				printf("请依次输入X(n)的复数序列值(实部 虚部):\n");                        
				for(i=0;i<N;i++)
					{
						scanf("%lf %lf",&x[i].real,&x[i].img);
              
					}
				fft();  
				printf("结果如下\n");   
            for(i=0;i<N;i++)									//输出结果
               {   
                 printf("%.4f",x[i].real);   
                 if(x[i].img>=0.0001)   
                 printf("+%.4fj\n",x[i].img);   
                 else   if(fabs(x[i].img)<0.0001)   
                 printf("\n");   
                 else     
                 printf("%.4fj\n",x[i].img);   
               }   
		        return   0;
				} 
		}
            void   fft()   
                {   
                  int   a,i=0,j=0,k=0,l=0;   
                  complex   HI,LO,OU;   //定义奇偶两部分和结果的复数结构
                  change();   //倒序列 
					for(a=0;a<N;a++)   
					{   
						W[a].real=cos(2*PI/N*a);   //指数从0到N-1的所有旋转因子
						W[a].img=-1*sin(2*PI/N*a);   
					}   
				  for(i=0;i<L;i++)
                       {       
                         l=1<<i;   //最外层循环，每循环一次左移一位，即*2，控制旋转因子的数量
                         for(j=0;j<N;j+= 2*l   )//在第L级中，同一个旋转因子对应间隔为2的L次方的蝶形
                           {                           
                             for(k=0;k<l;k++)
                                {                 
                                  mul(x[j+k+l],W[N*k/2/l],&OU);  //一次蝶形运算的过程 
                                  add(x[j+k],OU,&HI);   
                                  sub(x[j+k],OU,&LO); 
                                  x[j+k]=HI;   x[j+k+l]=LO;   
                                   }   
                              }   
                         }   
                    } 

void   change()   //序列的倒序
           {   
				complex   temp;  
				int i=0,j=0,k=0;
				for(i=0;i<N-1;i++)		//最后一个序列不需要交换，只需要循环N-1次
                {                   
					k=N/2;
					if(j>i)				//不相等时变换序列
						{   
							temp=x[i];   
							x[i]=x[j];   
							x[j]=temp;   
						}  
					while(k<=j)			//k是N/2，比K大则j最高位是1
                       {   
							j=j-k;		//j-k=j-N/2就是把j最高位清0
							k=k/2;		//比较次高位					
                        } 
				  j=j+k;				//0时置1
                   
                   }   
               }
        


