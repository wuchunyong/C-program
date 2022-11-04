#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654

float **creat_grid(int nx,int nz);//创建nx行nz列的网格
void free_grid(float **u,int nx);//释放nx行的网格
void update_grid(float **u1,float **u2,float **u3,float **v,float **hy,
                 int nx,int nz,int fn,float dt,float dx,float f,int t);//更新网格
void swap_grid(float ***x,float ***y);//交换两个网格
void write_data(float **u2,int nx,int nz);//写入nx行，nz列网格数据到csv文件
float focus1(float f,float dt,int t);//不含空间衰减函数的震源函数，雷克子波
float focus2(float **hy,float f,float dt,int t,int x,int z,int fn);//含空间衰减函数的震源函数，雷克子波，x，z为当前点坐标

int main()
{
    int nx=800,nz=800,nt=3000,fn=2;//nx，nz为网格行列数，nt为循环次数，fn为震源的个数
    float **u1,**u2,**u3,**v,**hy;//u1，u2，u3分别为前一时间，当前时间，后一时间的的位移网格数据，v为速度网格数据，hy为fn行两列的震源位置
    float dt=0.001,dx=10.0,f=12.5;//均匀网格，dx，dz统一用dx代替
    int i,j,t;
    u1=creat_grid(nx,nz);
    u2=creat_grid(nx,nz);
    u3=creat_grid(nx,nz);
    v=creat_grid(nx,nz);
    hy=creat_grid(fn,2);

    //确定各震源的位置
    for (i=0;i<fn;i++)
    {
        if(i<6)
        {
            hy[i][0]=400;
            hy[i][1]=300+50*i;
        }
        else
        {
            hy[i][0]=500;
            hy[i][1]=300+20*(i-10);
        }
        
    }

    //初始化位移网格，速度网格
    for(i=0;i<nx;i++)
    {
        for(j=0;j<nz;j++)
        {
            u1[i][j]=0;
            u2[i][j]=0;
            u3[i][j]=0;
            if(j<200)
            v[i][j]=1500.0;
            else if (j<600)
            v[i][j]=2000.0;
            else
            v[i][j]=2500.0;
            
        }
    }
    remove("1.csv");//删除当前文件夹下的1.csv文件
    for(t=0;t<nt;t++)//更新nt次网格
    {
        update_grid(u1,u2,u3,v,hy,nx,nz,fn,dt,dx,f,t);
        swap_grid(&u1,&u2);
        swap_grid(&u2,&u3);
        if(t%100==0)
        {
            write_data(u2,nx,nz);//每一百次写入一次数据
        }
    }
    free_grid(u1,nx);//释放网格
    free_grid(u2,nx);
    free_grid(u3,nx);
    free_grid(v,nx);
    free_grid(hy,fn);
    printf("OK");
    return 0;
}

float **creat_grid(int nx,int nz)
{
    int i;
    float **u;
    u=(float**)malloc(sizeof(float*)*nx);
    for(i=0;i<nx;i++)
        u[i]=(float*)malloc(sizeof(float)*nz);
    return u;
}

void free_grid(float **u,int nx)
{
    int i;
    for(i=0;i<nx;i++)
    free(u[i]);
    free(u);
    return;
}

void update_grid(float **u1,float **u2,float **u3,float **v,float **hy,
                int nx,int nz,int fn,float dt,float dx,float f,int t)
{
    int i,j;
    float p,mixed,second;
    for(i=1;i<nx-1;i++)
    {
        for(j=1;j<nz-1;j++)
        {
            p=v[i][j]*dt/dx;
            u3[i][j]=2*u2[i][j]-u1[i][j]+p*p*(u2[i+1][j]+u2[i-1][j]+u2[i][j+1]+u2[i][j-1]-4*u2[i][j]);
            u3[i][j]=u3[i][j]+focus2(hy,f,dt,t,i,j,fn);
            // if ((i==100&&j==100)||(i==600&&j==700))
            /*if(i==400&&j==400)
            {
                u3[i][j]=u3[i][j]+focus1(f,dt,t);
            }*/
        }
    }
    //上边界
    j=0;
    for(i=1;i<nx-1;i++)
    {
        p=v[i][j]*dt/dx;
        mixed=(u2[i][j+1]-u2[i][j])-(u1[i][j+1]-u1[i][j]);
        second=u2[i-1][j]-2*u2[i][j]+u2[i+1][j];
        u3[i][j]=2*u2[i][j]-u1[i][j]+p*mixed+0.5*p*p*second;
    }
    //下边界
    j=nz-1;
    for(i=1;i<nx-1;i++)
    {
        p=v[i][j]*dt/dx;
        mixed=-(u2[i][j]-u2[i][j-1])+(u1[i][j]-u1[i][j-1]);
        second=u2[i-1][j]-2*u2[i][j]+u2[i+1][j];
        u3[i][j]=2*u2[i][j]-u1[i][j]+p*mixed+0.5*p*p*second;
    }
    //左边界
    i=0;
    for(j=1;j<nz-1;j++)
    {
        p=v[i][j]*dt/dx;
        mixed=(u2[i+1][j]-u2[i][j])-(u1[i+1][j]-u1[i][j]);
        second=u2[i][j-1]-2*u2[i][j]+u2[i][j+1];
        u3[i][j]=2*u2[i][j]-u1[i][j]+p*mixed+0.5*p*p*second;
    }
    //右边界
    i=nx-1;
    for(j=1;j<nz-1;j++)
    {
        p=v[i][j]*dt/dx;
        mixed=-(u2[i][j]-u2[i-1][j])+(u1[i][j]-u1[i-1][j]);
        second=u2[i][j-1]-2*u2[i][j]+u2[i][j+1];
        u3[i][j]=2*u2[i][j]-u1[i][j]+p*mixed+0.5*p*p*second;
    }
    //四角
    u3[0][0]=0.5*(u2[0][1]+u2[1][0]);
    u3[0][nz-1]=0.5*(u2[0][nz-2]+u2[1][nz-1]);
    u3[nx-1][0]=0.5*(u2[nx-2][0]+u2[nx-1][1]);
    u3[nx-1][nz-1]=0.5*(u2[nx-2][nz-1]+u2[nx-1][nz-2]);
}

void swap_grid(float ***x,float ***y)
{
    float **temp;
    temp=*x;
    *x=*y;
    *y=temp;
    return;
}

void write_data(float **u2,int nx,int nz)
{
    FILE *file;
    int i,j;
    file=fopen("1.csv","a");
    //fseek(file, 0, SEEK_END);
    for(i=0;i<nx;i++)
    {
        for(j=0;j<nz;j++)
        {
            fprintf(file,"%.3f,",u2[i][j]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
    return;
}

float focus1(float f,float dt,int t)
{
    float s,t0=50;
    s=(1-2*pow(pi*f*(t-t0)*dt,2))*exp(-pow(pi*f*(t-t0)*dt,2));
    return s;
}
float focus2(float **hy,float f,float dt,int t,int x,int z,int fn)
{
    float s,h,sum=0,t0=50,a=0.2;
    int i,x0,z0;
    for(i=0;i<fn;i++)
    {
        x0=hy[i][0];
        z0=hy[i][1];
        s=(1-2*pow(pi*f*(t-t0)*dt,2))*exp(-pow(pi*f*(t-t0)*dt,2));
        h=exp(-a*a*(pow((x-x0),2)+pow((z-z0),2)));
        sum=sum+s*h;
    }
    return sum;
}