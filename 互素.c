#include <stdio.h>
int main()
{
	int i,x,y,t,flag=0;
	scanf("%d %d",&x,&y);
	if(x<y)
	{
		t=x;
		x=y;
		y=t;
	}
	printf("%d %d\n",x,y);
	for(i=2;i<=y;i++)
	{
		if(x%i==0&&y%i==0)
		{
			printf("%d\n",i);
			flag=1;
		}
	}
	printf("%d",flag);
	return 0;
}
