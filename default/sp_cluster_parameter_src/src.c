#include<stdio.h>

#define test_string "This is a test string/abc"
#define tar_file "sp_cluster_knot.txt"

int main(int argc,char*argv[])
{
    char fn_buffer_1[1024];
    char fn_buffer_2[1024];
    char fn_buffer_3[1024];

    FILE* fp = NULL;

    double var;
    

    sprintf(&fn_buffer_1[0],"sp_cluster_knot.txt");

    fp = fopen(&fn_buffer_1[0],"r");
    fscanf(fp,"%*d");

    for(int i=0;i<178;i++)
    {   fscanf(fp,"%lf\n",&var);
        printf("%lf\n",var);
    }

    fclose(fp);


//    printf("%s\n",test_string);

  //  sprintf(&fn_buffer_2[0],test_string);

    //printf("%s\n",fn_buffer_2);
    
  //  fp = fopen(tar_file,"r");


   // for(int i=0;i<178;i++)
   // {   fscanf(fp,"%lf\n",&var);
     //   printf("%lf\n",var);
    //}


    return 0;
}
