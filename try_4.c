#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
float **creat_grid(int nx,int nz);//创建网格
void free_grid(float **u,int nx);
void update_grid(float **u1,float **u2,float **u3,float **v,int nx,int nz,float dt,float dx,float f,int t);
void swap_grid(float ***x,float ***y);
void write_data();
float focus(float f,float dt,int t);


int main()
{
    int nx=800,nz=800,nt=3000;
    float **u1,**u2,**u3,**v;
    float dt=0.001,dx=10.0,f=12.5;
    int i,j,t;
    u1=creat_grid(nx,nz);
    u2=creat_grid(nx,nz);
    u3=creat_grid(nx,nz);
    v=creat_grid(nx,nz);
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
            v[i][j]=2500.0;
            else
            v[i][j]=3500.0;
            
        }
    }
    remove("1.csv");
    for(t=0;t<nt;t++)
    {
        update_grid(u1,u2,u3,v,nx,nz,dt,dx,f,t);
        swap_grid(&u1,&u2);
        swap_grid(&u2,&u3);
        if(t%100==0)
        {
            write_data(u2,nx,nz);
        }
    }
    free_grid(u1,nx);
    free_grid(u2,nx);
    free_grid(u3,nx);
    free_grid(v,nx);
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

void update_grid(float **u1,float **u2,float **u3,float **v,int nx,int nz,float dt,float dx,float f,int t)
{
    int i,j;
    float p,mixed,second;
    for(i=1;i<nx-1;i++)
    {
        for(j=1;j<nz-1;j++)
        {
            p=v[i][j]*dt/dx;
            u3[i][j]=2*u2[i][j]-u1[i][j]+p*p*(u2[i][j+1]+u2[i][j-1]+u2[i+1][j]+u2[i-1][j]-4*u2[i][j]);
            //if ((i==100&&j==100)||(i==600&&j==700))
            if(i==400&&j==400)
            {
                u3[i][j]=u3[i][j]+focus(f,dt,t);
            }
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

float focus(float f,float dt,int t)
{
    float x,t0=0;
    x=(1-2*pow(pi*f*(t-t0)*dt,2))*exp(-pow(pi*f*(t-t0)*dt,2));
    return x;
}