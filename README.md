# bamMetrics
## Docker location: https://hub.docker.com/repository/docker/sam16711/bam_metrics
### Docker files: [bamMetrics/docker_image_files](./docker_image_files)
> If you need to rebuild image, see files above.  
  
**NOTE: Docker image is used in exome pipeline on Washu Compute0 and Compute1. If the docker image isn't working, see below for standalone version and buliding.**  


Modified copy of Jim Knight's bamMetrics tool

bamMetrics came from /home/bioinfo/software/knightlab/soft/gatkutil_Apr2019/bamMetrics.cpp at Yale

Compiled using the command on compute0: ```g++ -O2 -o bamMetrics_docker bamMetrics.cpp```

The compiled program is labeled ```bamMetrics_docker``` and can be found here: ```/gscmnt/gc2698/jin810/programs```. See below for changes.

The Makefile specifies the path to samtools. This is defaulted to ```/opt/samtools/bin/samtools``` since this is where the mgibio/samtools-cwl:1.0.0 Docker image puts samtools.

*If you would like to run bamMetrics_docker outside of a cwl pipeline, you will need to run it with the mgibio/samtools-cwl:1.0.0 docker image.

## Building
Clone the repository and run ```make```.

If you've installed samtools somewhere other than ```/opt/samtools/bin/samtools``` you will need to change the ```sampath``` variable in the Makefile.

If you've made changes after building run ```make clean``` before making a new build.
