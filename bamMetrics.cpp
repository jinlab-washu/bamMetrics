#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <libgen.h>
#include <vector>
#include <string>

using namespace std;

#define DEFAULT_REF_FILE "/home/bioinfo/software/knightlab/genomes/hg38/grch38/grch38_primary_plus_decoy.fasta"
#define DEFAULT_BED_FILE "/home/bioinfo/software/knightlab/genomes/hg38/bed_files/hg38_coding_padded15_Jan2017.bed"

int numdepths = 11;
int depths[] = { 1, 2, 4, 8, 10, 15, 20, 30, 40, 50, 100 };

typedef struct {
	string chr, anno;
	int chrnum, start, end;
} BED_REGION;

typedef struct {
	string chr;
	off_t offset;
	int seqlen;
} CHR_INFO;

void readBedFile(char *filename, vector<BED_REGION>& regions);
void setAlign(char **fields, int& start, int& end, int &alen, char *qalign, char *ralign, char *qscore);
void setChromInfo(char *reffile, char *covfile, char *depthfile);
void resetChrSeq(char *chr, int curchrnum, int chrnum);
void settabs(char *s, char **fields);
int chr2num(char *chr);

vector<CHR_INFO> chromlist;

int reflen = 0;
char *refseq = NULL;
int *refcnts = NULL;
int *refcntsmapq0 = NULL;
char *refflags = NULL;
FILE *reffp = NULL;
vector<BED_REGION> regions;

FILE *covfp = NULL;
FILE *depthfp = NULL;

bool genomeFlag = false;
bool oneFlag = false;

uint64_t targetBases = 0;
uint64_t targetDepths[1000];
uint64_t targetDepthsMapq0[1000];

int genomeTargetReads = 0;
int genomeXReads = 0;
int genomeYReads = 0;
int genomeOtherReads = 0;

int readLenBins[1000];

int main(int argc, char *argv[])
{
	char *bamfile = NULL;
	char *bedfile = NULL;
	char *reffile = NULL;
	char *outfile = NULL;
	char *covfile = NULL;
	char *depthfile = NULL;

	bool countDupsFlag = false;
	bool longReadFlag = false;
	int minqscore = 0;

	int numThreads = 8;

	int arg = 1;
	while (arg < argc) {
		if (arg + 1 < argc && strcmp(argv[arg], "-o") == 0) {
			outfile = argv[arg+1];
			arg += 2;
		} else if (arg + 1 < argc && strcmp(argv[arg], "-b") == 0) {
			bedfile = argv[arg+1];
			arg += 2;
		} else if (arg + 1 < argc && strcmp(argv[arg], "-r") == 0) {
			reffile = argv[arg+1];
			arg += 2;
		} else if (arg + 1 < argc && strcmp(argv[arg], "-c") == 0) {
			covfile = argv[arg+1];
			arg += 2;
		} else if (arg + 1 < argc && strcmp(argv[arg], "-d") == 0) {
			depthfile = argv[arg+1];
			arg += 2;
		} else if (arg + 1 < argc && strcmp(argv[arg], "-q") == 0) {
			if (sscanf(argv[arg+1], "%d", &minqscore) != 1) {
				fprintf(stderr, "Error:  Invalid -q option value:  %s\n", argv[arg+1]);
				exit(-1);
			}
			minqscore += 33;
			arg += 2;
		} else if (arg + 1 < argc && strcmp(argv[arg], "-t") == 0) {
			if (sscanf(argv[arg+1], "%d", &numThreads) != 1) {
				fprintf(stderr, "Error:  Invalid -t option value:  %s\n", argv[arg+1]);
				exit(-1);
			}
			arg += 2;
		} else if (strcmp(argv[arg], "-g") == 0) {
			genomeFlag = true;
			arg += 1;
		} else if (strcmp(argv[arg], "-1") == 0) {
			oneFlag = true;
			arg += 1;
		} else if (strcmp(argv[arg], "--countdups") == 0) {
			countDupsFlag = true;
			arg += 1;
		} else if (strcmp(argv[arg], "--long") == 0) {
			longReadFlag = true;
			arg += 1;
		} else {
			break;
		}
	}

	if (genomeFlag && (bedfile || covfile)) {
		fprintf(stderr, "Error:  The whole genome mode cannot be used with the -b or -c options.\n");
		exit(-1);
	}

	if (arg + 1 != argc) {
		fprintf(stderr, "Usage:  bamMetrics [-g] [-1] [-b bedFile] [-r refFile] [-o outputFile] [-c coverageOutputFile] [-d depthOutputFile] [--countdups] [--long] [-q #] bamFile\n");
		exit(-1);
	}

	bamfile = argv[arg];
	if (bedfile == NULL) {
		bedfile = strdup(DEFAULT_BED_FILE);
	}
	if (reffile == NULL) {
		reffile = strdup(DEFAULT_REF_FILE);
	}

	FILE *outfp = NULL;
	if (outfile == NULL) {
		outfp = stdout;
	} else {
		outfp = fopen(outfile, "w");
		if (outfp == NULL) {
			fprintf(stderr, "Error:  Cannot write stats file:  %s\n", outfile);
			exit(-1);
		}
	}

	int curChrNum = 0;
	refseq = (char *) malloc(300000000 * sizeof(char));
	refcnts = (int *) malloc(300000000 * sizeof(int));
	refcntsmapq0 = (int *) malloc(300000000 * sizeof(int));
	refflags = (char *) malloc(300000000 * sizeof(char));

	setChromInfo(reffile, covfile, depthfile);

	char qalign[1000000], ralign[1000000], qscore[1000000];

	if (!genomeFlag) {
		readBedFile(bedfile, regions);
	}

	int numReads = 0;
	uint64_t numBases = 0;
	uint64_t numUniqueBases = 0;

	int numNondupReads = 0;
	int numNondupBases = 0;

	int numDuplicate = 0;
	int numUnmapped = 0;
	int numRepeat = 0;
	int numUnique = 0;

	uint64_t basecnt = 0;
	uint64_t errcnt = 0;

	int numTargetReads = 0;
	uint64_t numTargetBases = 0;

	for (int i=0; i < 1000; ++i) {
		targetDepths[i] = 0;
		targetDepthsMapq0[i] = 0;
	}

	for (int i=0; i < 1000; ++i) {
		readLenBins[i] = 0;
	}

	char cmd[1000];
	sprintf(cmd, "/opt/samtools/bin/samtools view -@ %d %s", numThreads, bamfile);
	FILE *fp = popen(cmd, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error:  Unable to read bamfile:  %s\n", bamfile);
		exit(-1);
	}

	int progressInterval = (longReadFlag ? 50000 : 1000000);

	int curchrnum = 0;

	char *fields[1000];
	char line[1000000];
	int lcnt = 0;
	while (fgets(line, 1000000, fp)) {
		if (line[0] == '@') {
			continue;
		}

		char origline[1000000];
		strcpy(origline, line);

		settabs(line, fields);
		
		int readlen = strlen(fields[9]);
		int flag = atoi(fields[1]);

		if ((flag & 0x100) || (flag & 0x900)) {
			continue;
		}

		++numReads;
		numBases += (uint64_t) readlen;

		if (readlen < 500000) {
			int bin = (int) (readlen / 1000.0);
			readLenBins[bin] += 1;
		} else {
			int bin = (int) (500000 / 1000);
			readLenBins[bin] += 1;
		}

		if (numReads % progressInterval == 0) {
			fprintf(stderr, "   -> %d reads, at %s:%s\n", numReads, fields[2], fields[3]);
			fflush(stderr);
		}

		// Set read counts (num, PCR dup, repeat)

		if (flag & 4) {
			++numUnmapped;
			continue;
		}

		if (flag & 0x400) {
			++numDuplicate;
			if (!countDupsFlag) {
				continue;
			}
		}
		
		int mapq = atoi(fields[4]);
		if (mapq == 0) {
			++numRepeat;
		} else {
			++numUnique;
		}

		char *chr = fields[2];
		int chrnum = chr2num(chr);
		if (chrnum == -1 || (genomeFlag && chrnum > 22)) {
			if (mapq >= 20) {
				if (chrnum == 24) {
					genomeXReads++;
				} else if (chrnum == 25) {
					genomeYReads++;
				} else {
					genomeOtherReads++;
				}
			}
			continue;
		}
		if (mapq >= 20) {
			genomeTargetReads++;
		}

		if (chrnum != curchrnum) {
			resetChrSeq(chr, curchrnum, chrnum);
			curchrnum = chrnum;
		}

		// Get alignment
		int start = 0;
		int end = 0;
		int alen = 0;
		setAlign(fields, start, end, alen, qalign, ralign, qscore);

		// Set coverage counts, error counts, on target counts
		int tflag = 0;
		for (int i=start; i <= end; ++i) {
			// ++refcntsmapq0[i];
			if (mapq == 0) {
				continue;
			}

			// ++refcnts[i];
			
			if (genomeFlag || refflags[i]) {
				tflag = 1;
				++numTargetBases;
			}
		}
		if (tflag) {
			++numTargetReads;
		}

		if (mapq == 0) {
			int rpos = start;
			for (int i=0; i < alen; ++i) {
				if (ralign[i] != '-') {
					if ((int) qscore[i] >= minqscore) {
						++refcntsmapq0[rpos];
					}
					++rpos;
				}
			}
			continue;
		}

		int rpos = start;
		for (int i=0; i < alen; ++i) {
			if (qalign[i] != '-') {
				++numUniqueBases;
			}
			++basecnt;
			if (qalign[i] != ralign[i]) {
				++errcnt;
			}
			if (ralign[i] != '-') {
				if ((int) qscore[i] >= minqscore) {
					++refcnts[rpos];
					++refcntsmapq0[rpos];
				}
				++rpos;
			}
		}

		/*
		printf("%s\n", origline);
		printf("Read %s, bc=%d  ec=%d, align %s:%d..%d\n   %s\n   %s\n",
			   fields[0], basecnt, errcnt, chr, start, end, qalign, ralign);
		printf("Continue? ");
		fflush(stdout);
		char resp[1000];
		fgets(resp, 1000, stdin);
		if (resp[0] == 'n') {
			pclose(fp);
			exit(0);
		}
		*/
	}
	pclose(fp);

	resetChrSeq(NULL, curchrnum, -1);

	fprintf(outfp, "Read Length:\t%d\n", (int) (numBases * 1.0 / numReads + 0.5));
	fprintf(outfp, "Num reads (M):\t%.1f\n", numReads / 1000000.0);
	fprintf(outfp, "Num bases (G):\t%.1f\n", numBases / 1000000000.0);
	
	uint64_t depthTotal = 0;
	uint64_t cumDepth = 0;
	uint64_t cumDepths[1000];
	uint64_t cumDepthMapq0 = 0;
	uint64_t cumDepthsMapq0[1000];
	for (int i=999; i >= 0; --i) {
		cumDepth += targetDepths[i];
		cumDepths[i] = cumDepth;
		depthTotal += targetDepths[i] * ((uint64_t) i);
		cumDepthMapq0 += targetDepthsMapq0[i];
		cumDepthsMapq0[i] = cumDepthMapq0;
	}
	uint64_t halfDepth = cumDepth / 2;
	int medianDepth = 0;
	for (int i=0; i < 1000; ++i) {
		if (cumDepths[i] < halfDepth) {
			break;
		}
		medianDepth = i;
	}
	
	fprintf(outfp, "Mean coverage:  \t%.1f\n", depthTotal * 1.0 / cumDepth);
	fprintf(outfp, "Median coverage:\t%d\n", medianDepth);

	fprintf(outfp, "PCR duplicates: \t%.2f%%\n", (numDuplicate * 100.0) / numReads);
	fprintf(outfp, "Multiply mapped:\t%.2f%%\n", (numRepeat * 100.0) / numReads);
	fprintf(outfp, "Unmapped:       \t%.2f%%\n", (numUnmapped * 100.0) / numReads);
	fprintf(outfp, "Reads on-target:\t%.2f%%\n", (numTargetReads * 100.0) / numReads);
	if (!genomeFlag) {
		fprintf(outfp, "Bases on-target:\t%.2f%%\n", (numTargetBases * 100.0) / numBases);
	} else {
		fprintf(outfp, "Genome Target Reads:\t%.2f%%\n", (genomeTargetReads * 100.0) / numReads);
		fprintf(outfp, "Genome ChrX Reads:\t%.2f%%\n", (genomeXReads * 100.0) / numReads);
		fprintf(outfp, "Genome ChrY Reads:\t%.2f%%\n", (genomeYReads * 100.0) / numReads);
		fprintf(outfp, "Genome MT/Decoy Reads:\t%.2f%%\n", (genomeOtherReads * 100.0) / numReads);
	}
	fprintf(outfp, "Mean error rate:\t%.2f%%\n", (errcnt * 100.0) / basecnt);

	for (int i=0; i < numdepths; ++i) {
		fprintf(outfp, "%dx target base coverage:\t%.1f%%\n",
				depths[i], (cumDepths[depths[i]] * 100.0) / cumDepths[0]);
	}
	
	for (int i=0; i < numdepths; ++i) {
		fprintf(outfp, "%dx target cov, incl MQ 0:\t%.1f%%\n",
				depths[i], (cumDepthsMapq0[depths[i]] * 100.0) / cumDepthsMapq0[0]);
	}

	if (longReadFlag) {
		fprintf(outfp, "Read Length Distribution:\n");
		for (int i=0; i < 50; ++i) {
			fprintf(outfp, "   %dk-%dk:\t%d\n", i, i+1, readLenBins[i]);
		}
		for (int i=50; i < 100; i += 10) {
			int cnt = 0;
			for (int j=i; j < i + 10; ++j)  {
				cnt += readLenBins[j];
			}
			fprintf(outfp, "   %dk-%dk:\t%d\n", i, i+10, cnt);
		}
		for (int i=100; i < 500; i += 50) {
			int cnt = 0;
			for (int j=i; j < i + 50; ++j)  {
				cnt += readLenBins[j];
			}
			fprintf(outfp, "   %dk-%dk:\t%d\n", i, i+50, cnt);
		}
		fprintf(outfp, "\t500k+:\t%d\n", readLenBins[500]);
	}

	
	fclose(outfp);

	if (depthfp != NULL) {
		fprintf(depthfp, "Depth\tCovered Base %\tDepth Distribution\n");
		for (int i=0; i <= 200; ++i) {
			fprintf(depthfp, "%d\t%.4f\t%.4f\n", i,
					(cumDepths[i] * 1.0) / cumDepths[0], 
					(i + 1 == 200 ? 0.0 : ( (cumDepths[i] * 1.0 / cumDepths[0]) - (cumDepths[i+1] * 1.0 / cumDepths[0]) )));
		}
		fclose(depthfp);
	}

	return 0;
}

void settabs(char *s, char **fields)
{
	int i = 0;

	fields[i] = s;
	for ( ; *s && *s != '\n'; ++s) {
		if (*s == '\t') {
			++i;
			fields[i] = s+1;
			*s = '\0';
		}
	}
	++i;
	fields[i] = NULL;
}

void setAlign(char **fields, int& start, int& end, int &alen, char *qalign, char *ralign, char *qscore)
{
	start = atoi(fields[3]);
	end = start - 1;
	alen = 0;

	char *readp = fields[9];
	char *scorep = fields[10];

	char *cigar = fields[5];

	int refpos = start - 1;
	for (char *s=cigar; *s; ++s) {
		int num = 0;
		for ( ; *s && isdigit(*s); ++s) {
			num = num * 10 + (*s - '0');
		}
		if (*s == 'M') {
			for (int i=0; i < num; ++i) {
				qalign[alen] = *readp++;
				qscore[alen] = *scorep++;
				ralign[alen] = refseq[refpos++];
				++end;
				++alen;
			}
		} else if (*s == 'D') {
			for (int i=0; i < num; ++i) {
				qalign[alen] = '-';
				qscore[alen] = '\0';
				ralign[alen] = refseq[refpos++];
				++end;
				++alen;
			}
		} else if (*s == 'I') { 
			for (int i=0; i < num; ++i) {
				qalign[alen] = *readp++;
				qscore[alen] = *scorep++;
				ralign[alen] = '-';
				++alen;
			}
		} else if (*s == 'S') {
			for (int i=0; i < num; ++i) {
				readp++;
				scorep++;
			}
		} else if (*s == 'H') {
		} else {
			fprintf(stderr, "Error:  Invalid cigar string:  %s\n", cigar);
			exit(-1);
		}
	}
	qalign[alen] = '\0';
	ralign[alen] = '\0';
	qscore[alen] = '\0';
	
	int last = -1;
	int lastscore = 0;
	for (int i=0; i < alen; ++i) {
		if (qscore[i] == '\0') {
			continue;
		}
		if (last + 1 < i) {
			int ls = (int) lastscore - 33;
			int ns = (int) qscore[i] - 33;
			if (lastscore == 0) {
				ls = ns;
			}

			int newscore = (ls + ns) / 2;
			char score = (char) (newscore + 33);
			while (last + 1 < i) {
				qscore[last+1] = score;
				++last;
			}
		}
		last = i;
	}
	while (last + 1 < alen) {
		qscore[last+1] = qscore[last];
		++last;
	}
}

void setChromInfo(char *reffile, char *covfile, char *depthfile)
{
	reffp = fopen64(reffile, "r");
	if (reffp == NULL) {
		fprintf(stderr, "Error:  Cannot open reference file:  %s\n", reffile);
		exit(-1);
	}

	char faifile[10000];
	sprintf(faifile, "%s.fai", reffile);
	FILE *fp = fopen(faifile, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error:  Cannot open .fai file:  %s\n", faifile);
		exit(-1);
	}
	char line[10000];
	char *fields[1000];
	while (fgets(line, 10000, fp)) {
		settabs(line, fields);
		char *s = fields[0];
		for ( ; *s && !isspace(*s); ++s) ;
		*s = '\0';

		CHR_INFO info;
		info.chr = fields[0];
		info.seqlen = atoi(fields[1]);
		info.offset = (uint64_t) 0;
		for (char *s=fields[2]; *s; ++s) {
			info.offset = (info.offset * (off_t) 10) + (off_t) (*s - '0');
		}
		chromlist.push_back(info);
	}
	fclose(fp);

	if (covfile != NULL) {
		covfp = fopen(covfile, "w");
		if (covfp == NULL) {
			fprintf(stderr, "Error:  Cannot write coverage file:  %s\n", covfile);
			exit(-1);
		}
		fprintf(covfp, "Chrom\tStart\tEnd\tAvg. Depth\tAvg. Depth (incl. MapQ 0)");
		for (int i=0; i < numdepths; ++i) {
			fprintf(covfp, "\t%dx", depths[i]);
		}
		fprintf(covfp, "\n");
	}

	if (depthfile != NULL) {
		depthfp = fopen(depthfile, "w");
		if (depthfp == NULL) {
			fprintf(stderr, "Error:  Cannot write coverage file:  %s\n", depthfile);
			exit(-1);
		}
	}
}

void resetChrSeq(char *chr, int curchrnum, int chrnum)
{
	// Compute previous chr coverage stats.
	if (curchrnum != 0) {
		if (genomeFlag) {
			if (curchrnum <= 22) {
				for (int pos=1; pos <= reflen; ++pos) {
					if (refseq[pos-1] == 'N') {
						continue;
					}
					int depth = refcnts[pos];
					++targetBases;
					targetDepths[min(depth,999)] += 1;
					int depthmapq0 = refcntsmapq0[pos];
					targetDepthsMapq0[min(depthmapq0,999)] += 1;
				}
			}
		} else {
			for (int i=0; i < (int) regions.size(); ++i) {
				if (regions[i].chrnum != curchrnum) {
					continue;
				}
				int regionBases = 0;
				int regionTotalDepth = 0;
				int regionTotalDepthMapq0 = 0;
				int regionDepths[1000];
				for (int j=0; j < 1000; ++j) {
					regionDepths[j] = 0;
				}
				for (int pos=regions[i].start; pos < regions[i].end; ++pos) {
					int depth = refcnts[pos];
					int depthmapq0 = refcntsmapq0[pos];

					++regionBases;
					regionDepths[min(depth,999)] += 1;

					regionTotalDepth += depth;
					regionTotalDepthMapq0 += depthmapq0;

					if (curchrnum != -1) {
						++targetBases;
						targetDepths[min(depth,999)] += 1;
						targetDepthsMapq0[min(depthmapq0,999)] += 1;
					}
				}
				float meanDepth = regionTotalDepth * 1.0 / regionBases;
				float meanDepthMapq0 = regionTotalDepthMapq0 * 1.0 / regionBases;

				if (covfp != NULL) {
					int cumDepth = 0;
					int cumDepths[1000];
					for (int j=999; j >= 0; --j) {
						cumDepth += regionDepths[j];
						cumDepths[j] = cumDepth;
					}

					fprintf(covfp, "%s\t%d\t%d\t%.1f\t%.1f", regions[i].chr.c_str(), regions[i].start, regions[i].end,
							meanDepth, meanDepthMapq0);
					for (int j=0; j < numdepths; ++j) {
						fprintf(covfp, "\t%.1f%%", (cumDepths[depths[j]] * 100.0) / cumDepths[0]);
					}
					fprintf(covfp, "\n");
					fflush(covfp);
				}

				/*
				if (depthfp != NULL) {
					fprintf(depthfp, ">%s\t%d\t%d\t%.1f", regions[i].chr.c_str(), regions[i].start, regions[i].end,
							meanDepth);
					for (int pos=regions[i].start; pos < regions[i].end; ++pos) {
						fprintf(depthfp, "%d\n", refcnts[pos]);
					}
				}
				*/
			}

			if (oneFlag) {
				fclose(covfp);
				exit(0);
			}
		}
	}

	if (chrnum == -1) {
		return;
	}

	string c = chr;

	int idx = 0;
	for ( ; idx < chromlist.size(); ++idx) {
		if (c == chromlist[idx].chr) {
			break;
		}
	}
	if (idx == chromlist.size()) {
		fprintf(stderr, "Error:  Chromosome not found in reference file:  %s\n", chr);
		exit(-1);
	}

	fseeko(reffp, chromlist[idx].offset, SEEK_SET);

	reflen = 0;
	char line[10000];
	while (reflen < chromlist[idx].seqlen) {
		if (!fgets(line, 10000, reffp)) {
			fprintf(stderr, "Error:  Reach EOF of reference file when reading chromosome:  %s\n", chr);
			exit(-1);
		}
		for (char *s=line; *s != '\n'; ++s) {
			refseq[reflen++] = toupper(*s);
		}
	}
	if (reflen != chromlist[idx].seqlen) {
		fprintf(stderr, "Error:  Could not ref chromosome from reference file:  %s\n", chr);
		exit(-1);
	}
	for (int i=0; i <= reflen; ++i) {
		refcnts[i] = 0;
		refcntsmapq0[i] = 0;
		refflags[i] = 0;
	}
	for (int i=0; i < (int) regions.size(); ++i) {
		if (regions[i].chrnum != chrnum) {
			continue;
		}
		for (int pos=regions[i].start; pos < regions[i].end; ++pos) {
			refflags[pos] = 1;
		}
	}
}
	
void readBedFile(char *filename, vector<BED_REGION>& regions)
{
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error:  Cannot read bed file:  %s\n", filename);
		exit(-1);
	}
	char line[10000];
	char *fields[1000];
	while (fgets(line, 10000, fp)) {
		settabs(line, fields);
		BED_REGION r;
		r.chr = fields[0];
		r.chrnum = chr2num(fields[0]);
		r.start = atoi(fields[1]);
		r.end = atoi(fields[2]);
		regions.push_back(r);
	}
	fclose(fp);
}

int chr2num(char *chr)
{
	int idx = (chr[0] == 'c' ? 3 : 0);
	if (isdigit(chr[idx]) && (chr[idx+1] == '\0' || (isdigit(chr[idx+1]) && chr[idx+2] == '\0'))) {
		return atoi(chr+idx);
	} else if (chr[idx] == 'M' && (chr[idx+1] == '\0' || (chr[idx+1] == 'T' && chr[idx+2] == '\0'))) {
		return 26;
	} else if (chr[idx] == 'X' && chr[idx+1] == '\0') {
		return 24;
	} else if (chr[idx] == 'Y' && chr[idx+1] == '\0') {
		return 25;
	} else {
		return -1;
	}
}
