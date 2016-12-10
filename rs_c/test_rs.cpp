#include "rs.h"

#include <stdlib.h> 
#include <stdio.h> 
#include <time.h> 

void main()
{
	int i,j;
	unsigned char data[kk], data1[kk], recd[nn], bb[nn-kk], dataout[kk];

	for(j = 0;j<1000;j++)
	{
		/* generate data sent */
		srand( (unsigned)time( NULL ) );
		for (i=0; i<kk; i++)
		{
			data[i] = rand()%256 ;
			data1[i] = data[i] ;
		}

			

		i =0;
		
		encode_rs_8(recd, data, bb);
		
		recd[rand()%15] = rand()%256;
		recd[rand()%10] = rand()%256;
		
		decode_rs_8(recd, dataout);
		
		for (i=0; i<kk; i++)
		{
			//printf("%5d",data[i]);
			if(data[i]!=dataout[i])
			{	
				printf("error\n");
			}
		}
	}
}
