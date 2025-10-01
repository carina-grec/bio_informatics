#include <stdio.h>
#include <string.h>

#define MAX_LENGTH 1000

void find_alphabet(char *s){
    int counts[256] = {0};
    int total = 0;
    for(int i = 0; s[i] != '\0'; i++){
        unsigned char c = s[i];
        counts[c]++;
        total++;
    }
    printf("Component percentages:\n");
    for(int i = 0; i < 256; i++){
        if(counts[i] > 0){
            double percent = (counts[i] * 100.0) / total;
            printf("%c: %.2f%%\n", i, percent);
        }
    }
}

int main(){
    char filename[256];
    char sequence[MAX_LENGTH * 10] = ""; 
    printf("Enter FASTA filename: ");
    fgets(filename, sizeof(filename), stdin);
    filename[strcspn(filename, "\n")] = 0;
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("Error opening file.\n");
        return 1;
    }
    char line[MAX_LENGTH];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') continue; 
        line[strcspn(line, "\n")] = 0;
        strcat(sequence, line);
    }
    fclose(fp);
    find_alphabet(sequence);
    return 0;
}