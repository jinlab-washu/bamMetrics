#!/bin/bash

#Check for pre-existing compute0 or compute1 reference cache for samtools cram support

if [ -d /gscmnt ] #check for compute0 root
then
	export REF_CACHE=/gscmnt/gc2698/jin810/references/samtools_ref/cache/%2s/%2s/%s
	export REF_PATH=/gscmnt/gc2698/jin810/references/samtools_ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
else
	if [ -d /storage1 ] #check for storage1 root (compute1)
	then
		export REF_CACHE=/storage1/fs1/jin810/Active/references/samtools_ref/cache/%2s/%2s/%s
		export REF_PATH=/storage1/fs1/jin810/Active/references/samtools_ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
	else
		echo "
		WARNING!!!	
		No REF_CACHE or REF_PATH set. You must create and set a ref_cache for using Samtools with CRAMS.
		Export paths using the following commands:
		export REF_PATH=/some_dir/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
		export REF_CACHE=/some_dir/cache/%2s/%2s/%s
		
		Use '/opt/bin/seq_cache_populate.pl -root /some_dir/cache yeast.fasta' to create MD5s cache for cram use.
		See here for more details: http://www.htslib.org/workflow/
		"
	fi
fi

exec "$@"
