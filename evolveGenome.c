/*

evolveGenome.c
Written by: Craig Lowe (lowec@stanford.edu)

*/

#include "common.h"
#include "linefile.h"
#include "options.h"
#include "memalloc.h"
#include "bed.h"
#include "fa.h"

/*---------------------------------------------------------------------------*/

static struct optionSpec optionSpecs[] =
/* command line option specifications */
{
	{"randomize", OPTION_BOOLEAN},
	{"maskDels", OPTION_BOOLEAN},
	{"dups", OPTION_INT},
	{"minDup", OPTION_INT},
	{"maxDup", OPTION_INT},
	{"timesToDup", OPTION_INT},
	{"dels", OPTION_INT},
	{"minDel", OPTION_INT},
	{"maxDel", OPTION_INT},
	{"mutRate", OPTION_DOUBLE},
	{"delsFile", OPTION_STRING},
	{NULL, 0}
};

boolean optRandomize = FALSE;
boolean optMaskDels = FALSE;
int optDups = 0;
int optMinDup = 1000;
int optMaxDup = 500000;
int optTimesToDup = 1;
int optDels = 0;
int optMinDel = 100;
int optMaxDel = 10000;
double optMutRate = 0.0;
char *optDelsFile = NULL;

/*---------------------------------------------------------------------------*/

void usage()
/* Explain usage and exit. */
{
errAbort(
	 "evolveGenome - randomly create amplifications and deletions in a genome\n"
	 "usage:\n"
	 "   evolveGenome in.fa noGap.bed out.fa ampsAndDels.bed\n"
	 "options:\n"
	 "    -randomize     (FALSE)     Seed the random number generator with the current system time\n"
	 "    -dups          (0)         Number of duplications\n"
	 "    -minDup        (1000)      Minimum size of a duplication\n"
	 "    -maxDup        (1000000)   Maximum size of a duplication\n"
	 "    -timesToDup    (1)         Number of times to duplicate the region\n"
	 "    -dels          (0)         Number of deletions\n"
	 "    -minDel        (100)       Minimum size of a deletion\n"
	 "    -maxDel        (10000)     Maximum size of a deletion\n"
	 "    -delsFile      (NULL)      Bed file of regions to be deleted (instead of -dels, -minDel, -maxDel)\n"
	 "    -mutRate       (0.0)       Rate of bp mutation\n"
	 "    -maskDels      (FALSE)     Instead of deleting the bases, turn them into Ns\n"
	 "notes:\n"
	 );
}

/*---------------------------------------------------------------------------*/

double possiblePositions(struct bed *bedList, unsigned int windowLength)
{
	double sum = 0;
	struct bed *futon = NULL;
	for(futon=bedList; futon!=NULL; futon=futon->next)
        {
		if(windowLength < futon->chromEnd - futon->chromStart)
		{
			sum += (double)(futon->chromEnd - futon->chromStart - windowLength);
		}
	}
	return(sum);
}

struct bed *randomlyPlaceNewBed(unsigned int placingLength, char *placingName, double randStart, struct bed *regionBedList)
{
	double length = 0;
	struct bed *answer = NULL, *futon = NULL;
	
	for(futon=regionBedList; futon!=NULL; futon=futon->next)
	{
		verbose(3, "start %e %s\t%u\t%u\n", randStart, futon->chrom, futon->chromStart, futon->chromEnd);
		verbose(3, "subtract: %u %u\n", futon->chromEnd - futon->chromStart, placingLength);
		if(placingLength > futon->chromEnd - futon->chromStart){length = 0;}
		else{length = (futon->chromEnd - futon->chromStart) - placingLength;}
		if(randStart < length)
		{
			verbose(3, "found: %f %f\n", randStart, length);
			AllocVar(answer);
			answer->chrom = cloneString(futon->chrom);
			answer->chromStart = futon->chromStart + randStart;
			answer->chromEnd = answer->chromStart  + placingLength;
			answer->name = cloneString(placingName);
			return(answer);
		}
		else
		{
			verbose(3, "a %f %f\n", randStart, length);
			randStart -= length;
			verbose(3, "b %f %f\n", randStart, length);
		}
	}
	errAbort("Error: at end of region list, but did not find a location to place the bed\n");
	return(NULL);
}


unsigned int pickInClosedInterval(unsigned int min, unsigned int max)
{
	return((double)(max-min+1) * rand() / (RAND_MAX + 1.0) + min);
}


double pickDoubleZeroToOne()
{
	return((double)rand() / (double)(RAND_MAX + 1.0));
}


boolean bedOverlap(struct bed *a, struct bed *b)
{
	assert(a != NULL);
	assert(b != NULL);
	if(sameString(a->chrom,b->chrom))
	{
		if(min(a->chromEnd,b->chromEnd) > max(a->chromStart,b->chromStart))
		{
	        return(TRUE);
	        }
	}
	return(FALSE);
}


boolean bedOverlapList(struct bed *a, struct bed *head)
{
	struct bed *curr = NULL;
	for(curr=head; curr != NULL; curr=curr->next)
	{
		if(bedOverlap(a, curr)){return(TRUE);}
	}
	return(FALSE);
}


struct bed *makeRandomlyPlacedBeds(int numberOfBeds, unsigned int minSize, unsigned int maxSize, char *name, struct bed *noGapBedList, struct bed *avoidMe)
{
	int i = 0;
	struct bed *answer = NULL, *randBed = NULL;
	unsigned int size = 0;
	double randStart = 0;

	i=0;
	while(i<numberOfBeds)
	{
		size = pickInClosedInterval(minSize, maxSize);
		randStart = possiblePositions(noGapBedList, size) * rand() / ( RAND_MAX + 1.0 );
		verbose(3, "size start positions: %u %e %e\n", size, randStart, possiblePositions(noGapBedList, size));
		randBed = randomlyPlaceNewBed(size, name, randStart, noGapBedList);
		if(!bedOverlapList(randBed, answer) && !bedOverlapList(randBed, avoidMe))
		{
			randBed->next = answer;
			answer = randBed;
			i++;
		}
	}
	slReverse(&answer);
	return(answer);
}


void bedOutAll(struct bed *printMe, char *outfile)
{
	struct bed *curr = NULL;
	FILE *f = mustOpen(outfile, "w");
	for(curr=printMe; curr != NULL; curr=curr->next)
	{
		bedTabOutN(curr, 4, f);
	}
	carefulClose(&f);
}


void printArray(char *buf, unsigned int length)
{
	unsigned int x = 0;
	for(x=0;x<length;x++)
	{
		uglyf("*%c", buf[x]);
	}
	uglyf("*\n");
}

void appendString(char *dest, char *src, unsigned int length)
{
	unsigned int x = 0;
	for(x=0;x<length;x++)
	{
		dest[x] = src[x];
	}
}


unsigned int basesCovered(struct bed *bedList)
{
	unsigned int total = 0;
	struct bed *futon = NULL;
	for(futon=bedList; futon!=NULL; futon=futon->next)
	{
		total += futon->chromEnd - futon->chromStart;
	}
	return(total);
}


struct dnaSeq *duplicateRegions(struct dnaSeq *seqList, struct bed *dupMeList, unsigned int copiesToMake)
{
	struct dnaSeq *currSeq = NULL, *answer = NULL;
	struct bed *currRegion = NULL;
	unsigned int seqSize=0, i=0, dupedBases = 0, oldSize=0;
	DNA *src = NULL, *dest = NULL, *newChrom = NULL;

	dupedBases = basesCovered(dupMeList) * copiesToMake;
	if(dupedBases == 0){return(NULL);}
	AllocArray(newChrom, dupedBases);

	for(currRegion = dupMeList; currRegion != NULL; currRegion = currRegion->next)
	{
		for(currSeq = seqList; currSeq != NULL; currSeq = currSeq->next)
		{
			if(sameString(currRegion->chrom, currSeq->name))
			{
				seqSize = currRegion->chromEnd - currRegion->chromStart;
				for(i=0; i<copiesToMake; i++)
				{
					if(answer == NULL)
					{
						answer = newDnaSeq(newChrom, dupedBases, "duplication");
					}
					src = currSeq->dna + currRegion->chromStart;
					dest = answer->dna + oldSize;
					appendString(dest, src, seqSize);
					oldSize += seqSize;
				}
			}
		}
	}
	return(answer);
}


void maskRegion(struct dnaSeq *seq, struct bed *region)
{
	unsigned int i = 0;

	for(i=region->chromStart; i < region->chromEnd; i++)
	{
		seq->dna[i] = 'N';
	}
}


void deleteRegion(struct dnaSeq *seq, struct bed *region)
{
	unsigned int i = 0, regionSize = 0;

	regionSize = region->chromEnd - region->chromStart;
	seq->size = seq->size - regionSize;
	for(i=region->chromStart; i < seq->size; i++)
	{
		seq->dna[i] = seq->dna[i+regionSize];
	}
}


void deleteRegions(struct dnaSeq *seqList, struct bed *delMeListSafe)
{
	struct bed *currRegion = NULL, *delMeList = NULL;
	struct dnaSeq *currSeq = NULL;

	delMeList = cloneBedList(delMeListSafe);
	slSort(&delMeList, bedCmp);
	slReverse(&delMeList);
	for(currRegion = delMeList, currSeq = seqList; currRegion != NULL; currRegion = currRegion->next)
	{
		verbose(3, "Deleting: %s %u %u\n", currRegion->chrom, currRegion->chromStart, currRegion->chromEnd);
		if(differentString(currRegion->chrom, currSeq->name))
		{
			for(currSeq = seqList; currSeq != NULL && differentString(currRegion->chrom, currSeq->name); currSeq = currSeq->next)
			{
			}
			if(currSeq == NULL){errAbort("%s is in the list of deletions, but not the fa file\n", currRegion->chrom);}
		}
		if(optMaskDels){maskRegion(currSeq, currRegion);}
		else{deleteRegion(currSeq, currRegion);}
	}
}


void mutateBase(struct dnaSeq *dna, unsigned int pos)
{
	unsigned int draw = 0;

	draw = pickInClosedInterval(0, 2);
	if(draw == 0)
	{
		if(dna->dna[pos] == 'A' || dna->dna[pos] == 'a'){dna->dna[pos] = 'T';}
		else{dna->dna[pos] = 'A';}
	}
	else if(draw == 1)
	{
		if(dna->dna[pos] == 'C' || dna->dna[pos] == 'c'){dna->dna[pos] = 'T';}
		else{dna->dna[pos] = 'C';}
	}
	else if(draw == 2)
	{
		if(dna->dna[pos] == 'G' || dna->dna[pos] == 'g'){dna->dna[pos] = 'T';}
		else{dna->dna[pos] = 'G';}
	}
}


void mutateGenome(struct dnaSeq *seqList, double mutationRate)
{
	struct dnaSeq *curr = NULL;
	unsigned int x = 0;

	for(curr=seqList; curr != NULL; curr=curr->next)
	{
		for(x = 0; x < curr->size; x++)
		{
			if(pickDoubleZeroToOne() <= mutationRate)
			{
				mutateBase(curr, x);
			}
		}
	}
}


/*---------------------------------------------------------------------------*/

void evolveGenome(char *inFaFile, char *noGapBedFile, char *outFaFile, char *outOperationsFile)
{
struct bed *noGapBedList = NULL, *deletions = NULL, *duplications = NULL, *operations = NULL;
struct dnaSeq *seqList = NULL, *newChroms = NULL;

if(optRandomize){srand(time(NULL));}

verbose(2, "reading in fa file and no gap file\n");
seqList = faReadAllDna(inFaFile);
noGapBedList = bedLoadNAll(noGapBedFile, 3);

verbose(2, "calcing where dups and dels will be\n");
duplications = makeRandomlyPlacedBeds(optDups, optMinDup, optMaxDup, "duplications", noGapBedList, NULL);
if(optDelsFile == NULL)
{
	deletions = makeRandomlyPlacedBeds(optDels, optMinDel, optMaxDel, "deletions", noGapBedList, duplications);
}
else
{
	deletions = bedLoadNAll(optDelsFile, 4);
}
verbose(2, "calculated %d dups and %d dels\n", slCount(duplications), slCount(deletions));

verbose(2, "making duplications\n");
newChroms = duplicateRegions(seqList, duplications, optTimesToDup);
verbose(2, "making deletions\n");
deleteRegions(seqList, deletions);
verbose(2, "calculated %d dups and %d dels\n", slCount(duplications), slCount(deletions));

verbose(2, "writing fa file\n");
seqList = slCat(seqList, newChroms);

if(optMutRate > 0)
{
	mutateGenome(seqList, optMutRate);
}

faWriteAll(outFaFile, seqList);

verbose(2, "writing bed file\n");
verbose(2, "calculated %d dups and %d dels\n", slCount(duplications), slCount(deletions));
operations = slCat(duplications, deletions);
verbose(2, "calculated %d operations\n", slCount(operations));

bedOutAll(operations, outOperationsFile);

}

/*---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, optionSpecs);
if (argc != 5)
    usage();

optRandomize = optionExists("randomize");
optMaskDels = optionExists("maskDels");
optDels = optionInt("dels", optDels);
optMinDel = optionInt("minDel", optMinDel);
optMaxDel = optionInt("maxDel", optMaxDel);
optDups = optionInt("dups", optDups);
optMinDup = optionInt("minDup", optMinDup);
optMaxDup = optionInt("maxDup", optMaxDup);
optTimesToDup = optionInt("timesToDup", optTimesToDup);
optMutRate = optionDouble("mutRate", optMutRate);
optDelsFile = optionVal("delsFile", optDelsFile);

if(optDelsFile != NULL && optDels != 0){errAbort("Error: can not do both -dels > 0 and -delsFile\n");}

evolveGenome(argv[1],argv[2],argv[3],argv[4]);
return(0);
}

