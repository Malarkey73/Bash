package SomaticCallerConfig;

require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(%Config);

=pod

=head1 NAME

 SomaticCallerConfig


=head1 DESCRIPTION

=head1 SYNOPSIS

 use SomaticCallerConfig qw(%Config);

 my $caveman_jar = $Config{CAVEMAN_JAR};

=head1 COPYRIGHT

 Copyright (c) 2011 Genome Research Ltd.

 Author: Alistair G. Rust <first_initial last_initial 12@sanger.ac.uk>
         Mamun Rashid     <last_initial first_initial  8@sanger.ac.uk>

 This file is part of Cake.

 Cake is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published
 by the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this program.  If not, see
 <http://www.gnu.org/licenses/>.
=cut

sub SCRIPTS_DIR() {'/shared/ucl/apps/cake/1.0/trunk/scripts/'}
#sub SOFTWARE_DIR() {'/Users/rm8/Software/Bioinformatics/'}

%Config = (

	   #################################
	   #  -- Sequencing data type  --  #
       #################################
	   #SEQUENCING_DATA_TYPE => 'exome',
	   SEQUENCING_DATA_TYPE => 'genome',  

	   
	   ########################################################################
	   #  							-- System Type   --					 	 ##  
       # 																	 ##
       # --------------------------------------------------------------------##
       ## Specifiy the Architecture type . Cake can be run both both 		 ##
       ## standalone and distributed system  								 ##
       #																	 ##
       #  System type Can be CLUSTER | STANDALONE							 ##
       #  Cluster type can only be LSF 
       ########################################################################
       
	   SYSTEM_TYPE => 'STANDALONE',
	   #CLUSTER_TYPE => 'LSF',


	   ## Specify maximum number of sample to run in one go ##
	   #MAX_NUM_SAMPLES_RUNNING => 10,
	   
	   ## Certain LSF system assigns submitted jobs to a particular
	   ## group with cerain resources. If you have such system, specify 
	   ## your group id, otherwise comment it out.
	   
	   #GROUP => 'team113',
	   

	   ########################################################################
	   #  		-- Queue Selection for Callers / Variant Filtering  --		 ##  
       # 																	 ##
       # --------------------------------------------------------------------##
       ## If you are using the distributed system [ CLUSTER ; LSF ] specify  ##
       ## the lsf queue system for different callers. Default is set 
       ## according to different running time 								 ## 
       ########################################################################

	   BAMBINO_QUEUE => 'normal',

	   MPILEUP_QUEUE => 'long',

	   VARSCAN_QUEUE => 'normal',

	   CAVEMAN_QUEUE => 'basement',
	   
	   SOMATICSNIPER_QUEUE => 'long',

	   ############################################
	   # -- Queue Selection for Post-Processing - #
	   ############################################
    
	   POST_PRO_QUEUE => 'normal',
	   
	   
	   #########################################################################
	   # -- 						Temporary Directory 			    	-- #
	   #########################################################################
		# compute nodes will use $TMPDIR for this
	  	TMP_DIR => '$ENV{'TMPDIR'}',
                

	   #########################################################################
	   # -- 					Post process filter flag			    	-- #
	   # -------------------------------------------------------------------- ##
	   # These Filters switch on / off different Variant Filtering Modules    ##
	   # Each filter is independent and produces a filtered VCF file. 		  ##
	   #########################################################################
	   
	   EXONIC_FILTER => 'TRUE',
	   SNP_FILTER => 'TRUE',
	   INDEL_FILTER => 'TRUE',
	   GERMLINE_FILTER => 'TRUE',
	   BASE_POS_FILTER => 'TRUE',
	   SOMATIC_COVERAGE_FILTER => 'TRUE',
	   ANNOTATE_CONSEQUENCE => 'FALSE',
	   FILTER_CONSEQUENCE => 'TRUE',
	   ALLELE_RATIO_FILTER => 'TRUE',
	   COSMIC_ANNOTATION_FLAG => 'TRUE',
	   
	   #########################################################################
	   # -- 					Relevant Software / Tools			    	-- #
	   # -------------------------------------------------------------------- ##
	   #  This segment contains executable paths [ bin / jar ] to all the 	  ## 
	   #  softwares required to successfully run Cake Variant Calling and 	  ##
	   #  Variant filtering pipeline										  ##	
	   #  																	  ##
	   #  *** Automatic Installation Script ***								  ##
	   #  Cake is distributed with a software / module installation script.	  ##
	   #  The role of the script is to download and install the necessary	  ##	
	   #  software and modules. The script also updates the links / Paths     ##
	   #  in this file. 													  ##
	   # 
	   #  *** Manual Installation  ***										  ##
	   #  If your are following the manual installation process, please do    ##   
	   #  remember to update the PATHs yourself.							  ##
	   #########################################################################
	  
       PERL_PATH => '/usr/bin/perl',
           
	   VARSCAN_JAR =>'/shared/ucl/apps/cake/1.0/support/VarScan.v2.3.6.jar',
	
	   #CAVEMAN_JAR => 'Set_the_Path_to_caveman_jar',
	   
	   SOMATIC_SNIPER =>'/shared/ucl/apps/cake/1.0/support/somatic-sniper-1.0.4/bin/bam-somaticsniper',

	   BAMBINO_JAR   =>'/shared/ucl/apps/cake/1.0/support/bambino_bundle_1.06.jar',

	   SAM_TOOLS_BIN =>'/shared/ucl/apps/samtools/0.1.19/samtools',
	   
	   # Required for somatic sniper  ##
	   SAM_TOOLS_1_6_BIN =>'/shared/ucl/apps/cake/1.0/support/somatic-sniper-1.0.4/vendor/samtools/samtools',

	   BCF_TOOLS_BIN =>'/shared/ucl/apps/samtools/0.1.19/bcftools/bcftools',

	   VCF_TOOLS_PATH =>'/shared/ucl/apps/cake/1.0/support/vcftools_0.1.12a/bin',

	   BGZIP =>'/shared/ucl/apps/cake/1.0/support/tabix-0.2.6/bgzip',

	   TABIX =>'/shared/ucl/apps/cake/1.0/support/tabix-0.2.6/tabix',
	   
	   #########################################################################
	   # -- 					Somatic calling scripts			    		-- #
	   # -------------------------------------------------------------------- ##
	   # Path to various Cake scrips. Please change the SCRIPTS_DIR variable  ##
	   #########################################################################

	   BAMBINO_CALLER => &SCRIPTS_DIR . 'bambino_caller.pl',

	   VARSCAN_CALLER => &SCRIPTS_DIR . 'varscan_caller.pl',

	   CAVEMAN_CALLER => &SCRIPTS_DIR . 'caveman_caller.pl',

	   MPILEUP_CALLER => &SCRIPTS_DIR . 'mpileup_caller.pl',

	   SOMATICSNIPER_CALLER => &SCRIPTS_DIR . 'somaticsniper_caller.pl',

	   SOMATIC_POST_PROCESSING_CALLER => &SCRIPTS_DIR . 'somatic_post_processing_caller.pl',

	   ENSEMBL_VARIANT_EFFECT_PREDICTOR => '/shared/ucl/apps/cake/1.0/support/ensembl-tools/scripts/variant_effect_predictor/variant_effect_predictor.pl',
	   
	   VCF_SUBSET_POSITION	=> &SCRIPTS_DIR . 'vcf-subset-positions.pl',
			
	   
	   #########################################################################
	   # -- 					Reference Genome Files			    		-- #
	   # -------------------------------------------------------------------- ##
	   # Specify location of the reference fasta file [ indexed ] for the	  ## 
	   # species you are going to analyse.   								  ##
	   #
	   # Make sure the the Fasta file used for BAM aligned (build, assembly)  ##
	   # is same as the one you are using for analysis. 					  ##
	   #																	  ##
	   # If you do not have indexed fasta [ .fai ] the following command      ## 
	   # samtools faidx FASTA_FILE 								
	   # 
	   #
	   # *** A sample fasta file containing first 5000000 bases of chr 19
	   # *** can be downloaded from 
	   # http://sourceforge.net/projects/cakesomatic/files/Relevant_Files/Other_Required_Files/GRCm38_um_just_sliced_chr19.fa
	   # http://sourceforge.net/projects/cakesomatic/files/Relevant_Files/Other_Required_Files/GRCm38_um_just_sliced_chr19.fa.fai
	   #
	   #########################################################################
	   REFERENCE_FASTA =>
	   {
	    # hg19
	    HUMAN => '/home/rmgzshd/Scratch/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa',
	    
	    # NCBIM37 / 38
	    #MOUSE => '/Users/rm8/Sanger/NGS_Analysis/Reference_Genome_Files/mouse/reference_genome/NCBIM37_um_mod1.fa',
	    # NCBIM38 / 38 - Chr 19 
	    MOUSE => '/shared/ucl/apps/cake/1.0/required/GRCm38_um_just_sliced_chr19.fa',
	   },



	   #########################################################################
	   # -- 					   Known SNP Files			    		    -- #
	   # -------------------------------------------------------------------- ##
	   # 
 	   # Cake uses user provided known SNP positions to filter out known 	  ##
 	   # single nucleotide polymorphisms.   								  ##
	   # 																	  ##
	   #																	  ##	
	   # The files must be in [ Chr \t Position ] format. Please refer to the ##
	   # Cake user guide for description and downlaod link for these knwon SNP##
	   # files. 															  ##
	   # 																	  ##
	   # Example : 															  ##
	   # Human Known SNPs :
	   # 	1. 1000 genome : http://www.1000genomes.org/data/		  		  ##	
	   # 	2. Washington Exomes : http://evs.gs.washington.edu/			  ##
       #	       evs_bulk_data/ESP6500SI-V2.coverage.all_sites.txt.tar.gz   ##
	   #																      ##
	   # 																      ##
	   # Mouse Known SNPs : 
	   #   1. Mouse Genome Project : 
	   #    ftp://ftp-mouse.sanger.ac.uk/REL-1211-SNPs_Indels/
	   #         mgp.v2.snps.annot.reformat.vcf.gz 
	   #
	   # If downloaded files are in VCF format Convert the vcf files in to    ##
	   # position file using the                  							  ##
	   #   scripts/vcf2bed.pl script in the script directory                  ##
	   #
	   #
	   # *** A test sample file can be downloaded from 
	   # *** http://sourceforge.net/projects/cakesomatic/files/Relevant_Files/Other_Required_Files/Test_SNP_File.pos
	   #
	   #
	   #########################################################################
	 
	   GENOMES =>
	   {
	    	MOUSE =>
	    	{
 			'KNOWN_MOUSE_SNP' => '/shared/ucl/apps/cake/1.0/required/Test_SNP_File.pos',
				
	    	},
	    	
	    	HUMAN =>
	    	{
	    		'KNOWN_HUMAN_SNP' =>'/home/rmgzshd/Scratch/Homo_sapiens/1000genomes/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.pos',
			},
	   },
	

	   #########################################################################
	   # -- 			Exonic co-ordinate / target interval file			-- #
	   # -------------------------------------------------------------------- ##
	   # 
 	   # Cake uses user provided Exonic co-ordinate / Gene co-ordinate files  ##
 	   # to filter out unwanted variants. 									  ##
	   # 																	  ##
	   # File format :														  ##	
	   # [ Chromosome \t Start \t End ]										  ##
	   #
	   # Exonic / Gene co-ordinate files can be downloaded from Ensembl/UCSC ##
	   # User may prefer to use custom bait design files for focusing only on ##
	   # region of interest
	   # 
	   # A brief discussion can be found here to obtain Ensembl/UCSC exonic   ##
	   # co-ordinates
	   # 1. http://www.biostars.org/p/6391/
	   # 2. Biomart is another good resource to obtain these co-ordinates     ##
	   # http://www.ensembl.org/biomart/martview/2efa45a9db099f5ce774533c3f5b7e59
	   #																	  ##
	   # If downloaded files are in VCF format Convert the vcf files in to    ##
	   # position file using the                  							  ##
	   #   scripts/vcf2bed.pl script in the script directory                  ##
	   
	   #
	   # *** A sample test file can be downloaded from
	   # http://sourceforge.net/projects/cakesomatic/files/Relevant_Files/Other_Required_Files/Test_exonic_cordinate_file.interval
	   #
	   #
	   #########################################################################
	   
	   
	   EXON_COORDINATES_FILE =>
	   {
	  		# Human Exon hg19
	   		HUMAN => '/home/rmgzshd/Scratch/Homo_sapiens/UCSC/hg19/Annotation/Exons/UCSCexons.txt',

			# Mouse Exonic Co-Ord
	   		MOUSE => '/shared/ucl/apps/cake/1.0/required/Test_exonic_cordinate_file.interval',
	   },


	   #########################################################################
	   # -- 						Known indel Regions						-- #
	   # -------------------------------------------------------------------- ##
	   # 
 	   # Variants falling into known indel regions are often unreliable and   ##
 	   # difficult to interpret. Cake uses user provided Known indel region   ##
 	   # file to discard any somatic variant occuring in these regions 		  ##
	   # 																	  ##
	   # File format :														  ##	
	   # [ Chromosome \t Start \t End ]										  ##
	   #
	   # Known indel regions can be obtained from 
	   # 1. Human : ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/
	   #     20121024_phase1_exome_indel_regenotyped/ALL.phase1_exome_indel_consensus_regenotyped_withAA.20110521.indels.exome.genotypes.vcf.gz
	   # 2. Mouse : ftp://ftp-mouse.sanger.ac.uk/REL-1211-SNPs_Indels/mgp.v2.indels.annot.reformat.vcf.gz
	   #  
	   # 
	   # If downloaded files are in VCF format Convert the vcf files in to    ##
	   # position file using the                  							  ##
	   #   scripts/vcf2bed.pl script in the script directory                  ##
	   #
	   #
	   # *** A sample test file can be downloaded from 
	   # *** http://sourceforge.net/projects/cakesomatic/files/Relevant_Files/Other_Required_Files/Test_Indel_interval_file.interval
	   #
	   #########################################################################
	   
	   INDEL_FILE =>
	   {
	   	# NCBIM37
	    MOUSE => '/shared/ucl/apps/cake/1.0/required/Test_Indel_interval_file.interval',
	
	    #hg19 see the REAME.txt in the folder with this
	    HUMAN => '/home/rmgzshd/Scratch/Homo_sapiens/UCSC/hg19/Annotation/Indels/indels.txt',
	    
	   },


	   #########################################################################
	   # -- 						Cosmic Annotation Data					-- #
	   # -------------------------------------------------------------------- ##
	   # 
 	   # Cake annotates human somatic variants with somatic mutation identified 
 	   # in COSMIC database					
	   # 																	  ##
	   # File format :														  ##	
	   # [ Chromosome \t Start \t End ]										  ##
	   #
	   # Cosmic annotation file :  
	   #  
	   # *** A Sample test file can be downloaded from 
	   # *** http://sourceforge.net/projects/cakesomatic/files/Relevant_Files/Other_Required_Files/Test_Cosmic_Region_Annotation_File.interval
	   #########################################################################
	   
	   #COSMIC_MUTATIONS_FILE => 'Set_the_Path_to_Human_Known_Indel_Location',



	   #########################################################################
	   # -- 						Caveman specific Data					-- #
	   # -------------------------------------------------------------------- ##
	   # 
 	   # Cake annotates human somatic variants with somatic mutation identified 
 	   # in COSMIC database					
	   # 																	  ##
	   # File format :														  ##	
	   # [ Chromosome \t Start \t End ]										  ##
	   #																	  ##
	   # Cosmic annotation file :  											  ##
	   #  																	  ##
	   # 																	  
	   #########################################################################
	   
	   #IGNORE_REGIONS_FILE =>
	   #{
	   # 	HUMAN => '/lustre/scratch103/sanger/rm8/GenomeFiles/requiredFiles/dummyIgnoreRegions.txt',
	   # 	MOUSE => '/lustre/scratch103/sanger/rm8/GenomeFiles/requiredFiles/dummyIgnoreRegions.txt',
	   #},
	   #
	   #
	   #COPY_NUMBER_FILES =>
	   #{
	   # 	HUMAN =>
	   # 	{
	   #  		TUMOUR => '/lustre/scratch103/sanger/rm8/GenomeFiles/requiredFiles/Human_Exome_Tumour_Default_CN.txt',
	   #  		NORMAL => '/lustre/scratch103/sanger/rm8/GenomeFiles/requiredFiles/Human_Exome_Normal_Default_CN.txt',
	   # 	},
	   # 	MOUSE =>
	   # 	{
	   #  		TUMOUR => '/lustre/scratch103/sanger/rm8/GenomeFiles/requiredFiles/Mouse_Exome_Tumour_Default_CN.txt',
 	   #		NORMAL => '/lustre/scratch103/sanger/rm8/GenomeFiles/requiredFiles/Mouse_Exome_Normal_Default_CN.txt',
	   # 	},
	   #},

	   
	   #########################################################################
	   # -- 						Variant Calling Paramters				-- #
	   # -------------------------------------------------------------------- ##
	   # 
 	   # Cake combines a number of somatic / variant calling algorithm and    ##
 	   # each of these callers has a number of parameters. Most of these      ##
 	   # paramters are common sequence quality paramters, however each variant##
 	   # caller has some paramter exclusive to themselves.                    ##
 
	   # *** Please refer to Table 2 in Cake user guide for more detail       ##
	   #     description about these paramter and their default values        ##
	   																	      ##
	   # 																	  ##
	   #########################################################################
	   
	   
	   MIN_MAPPING_QUALITY => 15,
	   MIN_BASE_POSITION_FREQUENCY => 0.6,
	   MIN_COVERAGE           => 4,
	   MIN_NUCLEOTIDE_QUALITY => 10,
	   
	   ## -- Bambino Specific parameter --
	   MIN_FLANKING_QUALITY   => 15,
	   MIN_ALT_ALLELE_COUNT    => 2,
	   MIN_MINOR_FREQUENCY    => 0,
	   MMF_MAX_HQ_MISMATCHES  => 5,
	   MMF_MIN_HQ_THRESHOLD   => 15,
	   MMF_MAX_LQ_MISMATCHES  => 6,
	   UNIQUE_FILTER_COVERAGE => 2,

	   ## -- Mpileup specific parameter --
	   BWA_DOWNGRADE_COFF     => 50,
	   NO_OF_READS_TO_CONSIDER_REALIGNMENT => 3,
	   FREQ_OF_READS          => 0.0002,

	   MPILEUP_QUALITY_THRESHOLD => 10,

	   ## -- Caveman specific parameter
	   NO_OF_BASES_TO_INCREMENT => 250000,


	   ## Somatic Sniper specific parameter ##
	   SOMATIC_QUALITY => 15,
	   

	   
	   

       #########################################################################
	   # -- 						Variant Filtering Paramters				-- #
	   # -------------------------------------------------------------------- ##
	   # 
 	   # Cake combines a number of variant filtering modules, which follows   ##
 	   # the variant calling pipeline. [ SNP filtering, Base Position filter, ##
 	   # Consequence Annotation and filtering, Variant Allele Ratio filter ]  ## 
 	   # 
	   # *** Please refer to Table 2 in Cake user guide for more detail       ##
	   #     description about these paramter and their default values        ##
	   																	      ##
	   # 																	  ##
	   #########################################################################
       
       #########################################################################
       ##	
       ##    NO_VCF_TO_INTERSECT 											  ##
	   ## -- No of Variant Caller output to intersect to indentify a variant  ##
	   ## -- as somatic mutation. By Default, for a variant to be somatic it  ##
	   ## -- needs to be identified as somatic by at least two of the variant ##
	   ## -- callers. 
	   
	   ##    NO_GENOTYPE_TO_INTERSECT 										  ##
	   ## -- After genomic position base intersection Cake intersects variants## 
	   ## -- based on genotype.  
	   
	   #########################################################################
	   
	   NO_VCF_TO_INTERSECT    => 2,
	   NO_GENOTYPE_TO_INTERSECT => 2,
	   
	   
	   
	   ###########################################################################
	   ## 				     -- Genotype Threshold --		 			   -- 	##
	   ##																	  	##
	   ## 	Many Variant caller does not produce any genotype for variants.   	##
	   ##   They merely spits out the number of reads showing both alleles	  	##
	   ##   Moreover Cancer genomes are often unstable, hence require
	   ##   more careful threhold selection is required to capture the 
	   ##   allelic variation. 
	   ## 																	  	##
	   ##	So we have set two different sets of thresholds 				    ##
	   ##			    Homozygous Ref | Heterozygous | Homozygous Alternative  ##
	   ##							   |		      |
	   ##   Tumor Minor                |              |
	   ## Allele Percentage  0 - 10%   |   11% - 50%  | 	51% - 100% 			##
	   ##		-----------------------|--------------|------------------------ ##
	   ##   Normal Minor               |              |
	   ## Allele Percentage  0 - 25%   |   26% - 75%  | 	76% - 100% 			##
	   ##																		##
	   ###########################################################################
	   
	   TUMOR_MIN_MINOR_ALLELE => 0.10, 
	   TUMOR_MAX_MINOR_ALLELE => 0.50,
	   
	   NORMAL_MIN_MINOR_ALLELE => 0.25,
	   NORMAL_MAX_MINOR_ALLELE => 0.75,
	   
	   
	   
	   #########################################################################
	   # -- 						Somatic Coverage Filter 				-- #
	   # -------------------------------------------------------------------- ##
	   # 
 	   # Cake detects somatic variant comparing Tumour and Normal reads. 	  ## 
 	   # So it is essential for both tumor and normal sample to have	      ## 
 	   # sufficient coverage at the genomic location under question           ##
	   # 																      ##
	   # So we apply a somatic coverage filter to make sure both tumor and 	  ##
	   # normal sample has user specified minimum covergae [ coverage maay ]  ##
	   # significantly vary between studies
	   #
	   # 
	   #########################################################################
       
	   
	   TUMOR_MIN_DEPTH => 10,
	   NORMAL_MIN_DEPTH => 10,
	   
	   
	   #########################################################################
	   ## 				     -- Filter Variant Allele Ratio 			   -- ##
	   ##																	  ##
	   ## 	1. VAR :  Ratio of Normal to Tumor variant allele ratio values    ##
	   ##   2. NVAR : Variant allle ratio in Normal samples 				  ##
	   ## 
	   ## 	Possible filter field and respective values :-> 				  ##
	   ##	1. VAR  -> 0.10												      ##
	   ##	2. NVAR -> 0.10													  ##
	   ##																      ##
	   #########################################################################
	   
	   FILTER_FOR_RATIO => 'VAR',
	   RATIO_FILTER_THRESHOLD => 0.10,
	   
	   
	   #########################################################################
	   ## 				    -- Consequnce Annotation and Filter			   -- ##
	   ##																	  ##
	   ## 	Cake Annotates variants according to their effect in Amino Acid 
	   ##   alteration using ensembl variant prediction tools. Variants with 
	   ##   non interesting [ user defined ] outcomes are discarded.          ##
	   ##    
	   ##
	   ##   Default Ensembl data base version used is 71. If you are using a  ##
	   ##	older version of the database and API, please comment the new     ##  
	   ##   ensembl CONSEQUENCE_TO_FILTER field and open the old one  
	   ##
	   ##	Also feel free to add or remove consequnces from the list accoring##
	   ##	to your need.													  ##
	   ##																	  ##	
	   ##
	   ##  ENSEMBL_DATABASE_VERSION 
	   ##  Specify which ensembl data base version you want to use
	   ##  Please make sure your ensembl database, API and VEP script version
	   ##  is in sync
	   ##
	   ##  ENSEMBL_CONNECTION_TYPE
	   ##  By default this is 'database' i.e. it will try to connect to ensembl
	   ##  over the internet and retrieve data. If you have a local cache
	   ##  configured, switch to cache mode by changing it to 'cache'
	   ##  http://www.ensembl.org/info/docs/variation/vep/vep_script.html#cache
	   ##
	   ## 
	   #########################################################################
	   
	   ## -- Consequence Filter -- ##
	   ENSEMBL_DATABASE_VERSION => 75,
	   ENSEMBL_CONNECTION_TYPE => 'database',
	   
	   ## Consequence for new ensembl version ##
	   CONSEQUENCE_TO_FILTER => 'synonymous_variant,intron_variant,3_prime_UTR_variant,5_prime_UTR_variant,upstream_gene_variant,synonymous_variant_NMD_transcript_variant,non_coding_exon_variant_nc_transcript_variant,intron_variant_nc_transcript_variant,intergenic_variant,downstream_gene_variant,',
	   
	   ## Consequence for old ensembl version ##
	   #CONSEQUENCE_TO_FILTER => 'SYNONYMOUS_CODING,INTRONIC,3PRIME_UTR,5PRIME_UTR,UPSTREAM,NMD_TRANSCRIPT,WITHIN_NON_CODING_GENE,INTERGENIC,DOWNSTREAM,',       
          
      	   
);


