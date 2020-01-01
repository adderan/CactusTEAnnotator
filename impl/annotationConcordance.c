#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "sonLib.h"

typedef struct {
    char *chrom;
    char *family;
    int64_t start;
    int64_t end;
    int64_t id;
} Annotation;


static int gffStartCompare(const void *a, const void *b) {
    Annotation *A = (Annotation*) a;
    Annotation *B = (Annotation*) b;
    if (A->start < B->start) return -1;
    if (A->start > B->start) return 1;
    if (A->end < B-> end) return -1;
    if (A->end > B-> end) return 1;
    return 0;
}

stSortedSet *readRepeatMaskerGff(char *gffFilename) {
    stSortedSet *annotations = stSortedSet_construct3(gffStartCompare, NULL);
    FILE *gff = fopen(gffFilename, "r");
    char a[10], b[10], c[10], d[10];
    char _chrom[100];
    int64_t _start;
    int64_t _end;
    char score[100];
    char strand;
    char _family[100];
    char e[10], f[10], g[10], h[10];
    int64_t _id;
    char *line = NULL;
    size_t nBytes;
    while(getline(&line, &nBytes, gff) != EOF) {
        int64_t ret = sscanf(line, "%s %s %s %s %s %ld %ld %s %c %s %s %s %s %s %ld\n",
            a, b, c, d, _chrom, &_start, &_end, score, &strand, _family, e, f, g, h, &_id);

        if (ret != 15) continue;
        Annotation *annotation = malloc(sizeof(Annotation));
        annotation->chrom = calloc(strlen(_chrom) + 1, sizeof(char));
        strcpy(annotation->chrom, _chrom);
        annotation->family = calloc(strlen(_family) + 1, sizeof(char));
        strcpy(annotation->family, _family);
        annotation->start = _start;
        annotation->end = _end;
        annotation->id = _id;

        stSortedSet_insert(annotations, annotation);

        fprintf(stderr, "%s %ld %ld %ld %s\n", annotation->chrom, annotation->start, annotation->end,
            annotation->id, annotation->family);
    }
    return annotations;
}

double overlap(Annotation *query, Annotation *target) {
    assert(target->start < target->end);
    assert(query->start < query->end);
    if ((target->start > query->end) || (query->start > target->end))
        return 0.0;
    int64_t overlapStart = (target->start > query->start) ? target->start : query->start;
    int64_t overlapEnd = (target->end < query->end) ? target->end : query->end;

    int64_t overlappingBases = overlapEnd - overlapStart + 1;
    double overlapFraction = (double) overlappingBases/ (double) (target->end - target->start);

    return overlapFraction;
}

int main(int argc, char **argv) {
    char *queryFilename = NULL;
	char *targetFilename = NULL;
    while (1) {
        static struct option long_options[] = {
            { "query", required_argument, 0, 'a' }, 
			{ "target", required_argument, 0, 'b'},
            { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                queryFilename = strdup(optarg);
                break;
			case 'b':
				targetFilename = strdup(optarg);
				break;
            default:
                return 1;
        }
    }
    
    stSortedSet *targetAnnotations = readRepeatMaskerGff(targetFilename);
    stSortedSet *queryAnnotations = readRepeatMaskerGff(queryFilename);

    
    stSortedSetIterator *targetIterator = stSortedSet_getIterator(targetAnnotations);
    stSortedSetIterator *queryIterator = stSortedSet_getIterator(queryAnnotations);

    fprintf(stderr, "Found %ld query annotations\n", stSortedSet_size(queryAnnotations));
    fprintf(stderr, "Found %ld target annotations\n", stSortedSet_size(targetAnnotations));
    stList *matchingAnnotations_query = stList_construct();
    stList *matchingAnnotations_target = stList_construct();

    //Find annotations that overlap by more than specified amount
    Annotation *target = stSortedSet_getNext(targetIterator);
    Annotation *query = stSortedSet_getNext(queryIterator);
    while(target && query) {
        if(overlap(target, query) > 0.5) {
            stList_append(matchingAnnotations_query, query);
            stList_append(matchingAnnotations_target, target);
            target = stSortedSet_getNext(targetIterator);
            query = stSortedSet_getNext(queryIterator);
        }
        else if (target->start < query->start) {
            target = stSortedSet_getNext(targetIterator);
        }
        else {
            query = stSortedSet_getNext(queryIterator);
        }
    }
    
    fprintf(stderr, "Found %ld overlapping annotations\n", stList_length(matchingAnnotations_query));

    int64_t nConcordant = 0;
    int64_t total = 0;
    for (int64_t i = 0; i < stList_length(matchingAnnotations_query); i++) {
        if (!(i % 50 == 0)) continue;
        Annotation *query_i = stList_get(matchingAnnotations_query, i);
        Annotation *target_i = stList_get(matchingAnnotations_target, i);
        for (int64_t j = 0; j < i; j++) { 
            Annotation *query_j = stList_get(matchingAnnotations_query, j);
            Annotation *target_j = stList_get(matchingAnnotations_target, j);
            bool a = (strcmp(query_i->family, query_j->family) == 0);
            bool b = (strcmp(target_i->family, target_j->family) == 0);
            if ((a && b) || (!a && !b)) {
                nConcordant++;
            }
            total++;
        }
    }
    fprintf(stderr, "Sampled %ld pairs\n", total);
    double concordance = (double) nConcordant / (double) total;
    fprintf(stderr, "Concordance = %lf\n", concordance);
    stSortedSet_destruct(targetAnnotations);
    stSortedSet_destruct(queryAnnotations);
}