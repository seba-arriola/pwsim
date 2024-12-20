#include "pwsim.h"

int isEmpty(const char *str)
{
    char ch;

    do
    {
        ch = *(str++);

        if(ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r' && ch != '\0')
            return 0;

    } while (ch != '\0');

    return 1;
}


int removeEmptyLines(char *path){
	
	FILE *srcFile;
    FILE *tempFile;
    
    srcFile  = fopen(path, "r");
    
    if(srcFile == NULL)
    {
		printf("File %s cannot be opened\n", path);
		return -1;
	}
    
    char temp_file[500];
    
    strcpy(temp_file,path);
    strcat(temp_file,"_tmp");
    
    tempFile = fopen(temp_file, "w");
	
    char buffer[BUFFER_SIZE];

    while ((fgets(buffer, BUFFER_SIZE, srcFile)) != NULL)
    {
        if(!isEmpty(buffer))
            fputs(buffer, tempFile);
    }
	
	fclose(srcFile);
    fclose(tempFile);

    remove(path);
    rename(temp_file, path);
	
	return 0;
	
}

/* function that counts lines in a text file
 * Receives the name of the text file
 * */
int countLines(char *filename)
{
	FILE *fp; 
    int count = 0;  // Line counter 
	char line[400];
	if(removeEmptyLines(filename)<0)
		return 0;
    fp = fopen(filename, "r");
    
    if(fp == NULL)
    {
		return 0;
	}
	
    while (!feof (fp))
		if (fgets(line, sizeof (line), fp))
			count++;

    fclose(fp); 
    return count;
}

