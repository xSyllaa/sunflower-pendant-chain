/* Independent naive triple-loop checker (V3R1T4S run 2026-10-19-1011). Input: lines of space-separated point
   indices (0-based), one member per line; '#' comments. Reports duplicates, Sperner violations, 3-sunflowers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
int main(int argc,char**argv){
  FILE*f=argc>1?fopen(argv[1],"r"):stdin; if(!f){perror("open");return 2;}
  size_t cap=1024,n=0; uint64_t*s=malloc(cap*8); char line[4096];
  while(fgets(line,sizeof line,f)){ if(line[0]=='#')continue; char*p=line; uint64_t m=0; int any=0;
    while(*p){ while(*p==' '||*p=='\t'||*p==',')p++; if(*p<'0'||*p>'9')break; long v=strtol(p,&p,10); if(v<0||v>63){fprintf(stderr,"bad point\n");return 2;} m|=1ULL<<v; any=1;}
    if(!any && strspn(line," \t\r\n")==strlen(line)) { /* blank line = empty set only if marked */ continue; }
    if(n==cap){cap*=2;s=realloc(s,cap*8);} s[n++]=m; }
  long dup=0,sp=0,sf=0;
  for(size_t i=0;i<n;i++)for(size_t j=0;j<n;j++){ if(i==j)continue; if(s[i]==s[j]){if(i<j)dup++;} else if((s[i]&s[j])==s[i])sp++; }
  for(size_t i=0;i<n;i++)for(size_t j=i+1;j<n;j++){ uint64_t I=s[i]&s[j];
    for(size_t k=j+1;k<n;k++){ if((s[i]&s[k])==I && (s[j]&s[k])==I){ if(sf<3)fprintf(stderr,"sunflower %zu %zu %zu\n",i,j,k); sf++; } } }
  printf("members %zu duplicates %ld sperner_violations %ld sunflowers %ld\n",n,dup,sp,sf);
  return (dup||sp||sf)?1:0; }
