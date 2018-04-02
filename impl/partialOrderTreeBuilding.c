#include <stdio.h>
#include <regex.h>
#include "sonLib.h"

struct partition {
  stList *threadSets;
};

struct threadSet {
  stList *threads;
}



stHash *parsePO(FILE *poFile) {
  stHash *partitions = stHash_construct();
  regex_t regex;
  regcomp(&regex, "[ASL]\d+");

  char *line;
  while ((line = getLine(poFile) != NULL)) {
    char base;
    char *nodeInfo;
    int num = sscanf(line, "%c:%s\n", &base, nodeInfo);
    if (num != 2) continue;
    int match = regexec(&regex, nodeInfo, 0, NULL, 0);
    assert(match != REG_NOMATCH);

    fprintf(stderr, "Parsed node with base %c\n", base);

  }


int main(int argc, char **argv) {

}
