# bamMetrics
Modified copy of Jim Knight's bamMetrics tool


bamMetrics came from /home/bioinfo/software/knightlab/soft/gatkutil_Apr2019/bamMetrics.cpp at Yale

Compiled using the command on compute0: ```g++ -O2 -o bamMetrics_docker bamMetrics_mgibio_docker.cpp```

***NOTE: The base bamMetrics.cpp has been modified on line 199 to point to the location of the samtools program in the mgibio/samtools-cwl:1.0.0 docker image.

The compiled program is labeled ```bamMetrics_docker``` and can be found here: ```/gscmnt/gc2698/jin810/programs```. See below for changes.

Modified:
```
spencer.king@virtual-workstation2:/gscmnt/gc2698/jin810/programs$ diff bamMetrics.cpp bamMetrics_mgibio_docker.cpp 
198c198
< 	sprintf(cmd, "/home/jk2269/soft/samtools-1.8/bin/samtools view -@ %d %s", numThreads, bamfile);
---
> 	sprintf(cmd, "/opt/samtools/bin/samtools view -@ %d %s", numThreads, bamfile);
```
This is because the mgibio/samtools-cwl:1.0.0 Docker image puts samtools at /opt/samtools/bin/samtools.

*If you would like to run bamMetrics_docker outside of a cwl pipeline, you will need to run it with the mgibio/samtools-cwl:1.0.0 docker image.

*If you would like to run bamMetrics using a different location of samtools (or without docker), you will need to change ```/opt/samtools/bin/samtools``` in line 199 to point to your samtools location.
