package Parser::parse_options;

use strict;
use Exporter 'import';
use Getopt::Long::Descriptive;

our @EXPORT_OK = qw(parse_options);

# Function that will parse the user options 
sub parse_options {
    my user_args = @_;

    my ($opts, $usage)= describe_options(
        "Usage: compadre_kickoff.pl %o",                                  
            
            "For usage and documentation:",  
            [ 'help|h|?', "Brief help message", { shortcircuit => 1 } ],      
            [ 'man', "Print manual page", { shortcircuit => 1 } ],    

            "Required options (one of the following):",
            [ 'plink_ibd|p=s', "Specify path to a .genome IBD estimates file produced by PLINK" ],
            [ 'input|i=s{1,9}', "Specify path to an IBD estimates file and additional column information (Ex. FILE=test.txt IBD0=5 IBD1=6 )", {
                test => sub
                    { 
                        my $val = shift;
                        my $opt_name = shift;
                        my @possible_keys = qw(FILE FID1 IID1 FID2 IID2 IBD0 IBD1 IBD2 PI_HAT RELATEDNESS);
                        my ($key,$value) = split(/=/, $val);
                        $key = uc($key);
                        if(grep ($_ eq $key, @possible_keys ))
                        {
                            $ibd_estimates{$key} = $value;
                        }
                        else
                        {
                            die "\n\n[COMPADRE] Error: Invalid key before \"=\" for --$opt_name: $key\n" . "Must be one of these: " . join(', ', @possible_keys) .  "\n\n";
                        }
                    },
                } 
            ],

            "    (or both):"
            [ 'file=s', "Path to PLINK formatted data without the file extensions; behaves the same as in PLINK meaning this expects the file path up to the file prefix (requires --genome)" ], 
            [ 'genome', "Read in --file and calculate IBD estimates using PLINK" ],  

            "COMPADRE options (new):"
            [ 'segment_data=s', "Shared segments data in format readable by ERSA (see full documentation)" ], 
            [ 'port_number=i', "Port number for additional Python computation ", { default => 6000 , default_needed => 1 } ],  
            [ 'run_padre', "Run PADRE after standard pedigree reconstruction is complete", { default => 0 } ],

            "General options:",
            [ 'rel_threshold|t=f', "Set the minimum level of relatedness for two people to be considered related", {  default => 0.09375, default_needed => 1 } ],
            [ 'degree_rel_cutoff=i', "Set the maximum degree of relatedness for two people to be considered related", { default => 3, default_needed => 1} ],  
            [ 'output_dir|o=s', "Specify path to the output directory for all results" ],   
            [ 'verbose|v=i', "Set verbosity level (0=none; 1=default; 2=more; 3=max)", { default => 1, default_needed => 1 } ], 

            "prePRIMUS IBD estimation options:",
            [ 'plink_ex=s', "Path to the plink executable file (searches environment variables by default)" ],
            [ 'ref_pops=s', "Comma separated list of 1000 Genomes populations used for reference allele freqs (overrides default method)" {
                # Make sure that we have valid 1KG populations
                test => sub 
                    {
                        my $val = shift;
                        @ref_pops = split(/,/,$val);
                        @ref_pops = uc(@ref_pops);
                        foreach my $pop (@ref_pops)
                        {
                            next if $pop =~ /none/i;
                            
                            if(!grep ($_ eq $pop, @onekg_pops ))
                            {
                                die "\n\n[COMPADRE] Error: Invalid 1KG population: $pop\n" . 
                                "Must be a comma seperated list these: @onekg_pops\n" . 
                                "For example: --ref_pops CEU,TSI,YRI\n\n";
                                
                            }
                        }
                        return 1; # Return true to pass validation
                    },
                } 
            ], 
            [ 'no_automatic_IBD', "Turn off automatic selection of the HapMap3 populations for reference allele freqs (On by default)", { default => 1 } ],    
            [ 'remove_AIMs', "Automatically remove ancestry informative markers", { default => 0 } ],  
            [ 'keep_AIMs', "Do not remove ancestry informative markers", { default => 0 } ],   
            [ 'internal_ref', "Use the dataset provided in --file to get reference allele frequencies", { default = 0 } ],
            [ 'alt_ref_stem=s', "Path to PLINK formatted data (no file extensions) used for allele frequencies" ],   
            [ 'keep_inter_files', "Keep prePRIMUS intermediate files" ],  
            [ 'no_mito=i', "Disable mitochondrial data (1 to disable, 0 to enable)", { default => 1 } ],  
            [ 'no_y=i', "Disable y data (1 to disable, 0 to enable)", { default => 1 } ],  
            [ 'min_pihat_threshold=f', "Set a minimum pi-hat threshold that will be used in the plink", { default => 0, default_needed => 1 } ], 
            [ 'gzip_genome', "Request gzipped output for the .genome generated by the PLINK ibd detection (off by default)" ], 
            [ 'samples=s', "A file containing a list of all samples (e.g., a .fam file) to ensure individuals with no relatives are included" ],   
            [ 'max_memory=i', "Specify amount of memory to be used in PLINK prePRIMUS commands (in MB)" ],
            [ 'no_PCA_plot', "Forgo running pca and getting PCA plots in order to get a faster runtime", { default => 0, default_needed => 1} ],   

            "Identification of maximum unrelated set (IMUS) options:",
            [ 'no_IMUS', "Don't identify a maximum unrelated set", { default => 1} ], 
            [ 'missing_val=f', "Set value that denotes missing data in IBD file", { default => 1, default_needed => 1} ],  
            [ 'size', "Specify to weight on set size (Default unless a binary trait is specified first)" ], 
            [ 'high_btrait=s', "File with FID, IID, and binary trait to weight for the higher value" ],  
            [ 'low_btrait=s', "File with FID, IID, and binary trait to weight for the lower value" ],                           
            [ 'high_qtrait=s', "File with FID, IID, and quantitative trait to weight for higher values" ],                   
            [ 'low_qtrait=s', "File with FID, IID, and quantitative trait to weight for lower values" ],                     
            [ 'mean_qtrait=s', "File with FID, IID, and quantitative trait to weight towards the mean value" ],                   
            [ 'tails_qtrait=s', "File with FID, IID, and quantitative trait to weight against the middle values" ],   
            [ 'generate_likelihoods_only', "only generate the likelihood vectors", { default => 0 }],
            ["int_likelihood_cutoff=f", "initial cutoff for likelihood ratios", { default => 0.3, default_needed => 1 }],
            ["missing_val=f", "value that represents missing data within the provided trait file"],

            "Pedigree Reconstruction (PR) options:",
            [ 'no_PR', "Don't reconstruct pedigrees" ],   
            [ 'max_gens=i', "Max number of generations sampled in reconstructed pedigree" ],
            [ 'max_gen_gap=i', "Max number of generations between two people that have a child", { default => 0 } ],  
            [ 'age_file=s', "Specify path to the file containing the age of each sample" ],   
            [ 'ages=s{1,4}', "Like --age_file but requires FILE=[file], optional specification of file columns", {
                test => sub
                    { 
                        my $val = shift;
                        my $opt_name = shift;
                        my @possible_keys = qw(FILE FID IID AGE);
                        my ($key,$value) = split(/=/,$val);
                        $key = uc($key);
                        if(grep ($_ eq $key, @possible_keys ))
                        {
                            $ages{$key} = $value;
                        }
                        else
                        {
                            die "\n\n[COMPADRE] Error: Invalid key '$key' for --$opt_name option.\n" . "Must be one of these: " . join(', ', @possible_keys) . "\n\n";

                        }
                        return 1;
                    },
                } 
            ], 
            [ 'sex_file=s', "Specify path to the file containing the sex of each sample" ],  
            [ 'sexes=s{1,6}', "Like --sex_file but requires FILE=[file], optional specification of file columns", {
                test => sub
                    { 
                        my $val = shift;
                        my $opt_name = shift;
                        my @possible_keys = qw(FILE FID IID SEX MALE FEMALE);
                        my ($key,$value) = split(/=/,$val);
                        $key = uc($key);
                        if(grep ($_ eq $key, @possible_keys ))
                        {
                            $sexes{$key} = $value;
                        }
                        else
                        {
                            die "\n\n[COMPADRE] Error: Invalid key '$key' for --$opt_name option.\n" ."Must be one of these: " . join(', ', @possible_keys) . "\n\n";
                        }
                        return 1;
                    },
                } 
            ],    
            [ 'mito_data|mito_matches|mito=s{1,7}', "Path to mito matching status for each pair of samples (requires FILE=[file])", {
                test => sub 
                    {
                        my $val = shift;
                        my $opt_name = shift;
                        my @possible_keys = qw(FILE FID1 IID1 FID2 IID2 MATCH MATCH_VAL);
                        my ($key,$value) = split(/=/,$val);
                        $key = uc($key);
                        if(grep ($_ eq $key, @possible_keys ))
                        {
                            $mito_matches{$key} = $value;
                        }
                        else
                        {
                            die "\n\n[COMPADRE] Error: Invalid key before \"=\" for --$opt_name: $val\n" . "Must be one of these: " . join(', ', @possible_keys) . "\n\n";
                        }
                        return 1;
                    },
                } 
            ],  
            [ 'y_data|y_matches|y=s{1,7}', "Path to y matching status for each pair of samples (requires FILE=[file])", {
                test => sub
                    { 
                        my $val = shift;
                        my $opt_name = shift;
                        my @possible_keys = qw(FILE FID1 IID1 FID2 IID2 MATCH MATCH_VAL);
                        my ($key,$value) = split(/=/, $val);
                        $key = uc($key);
                        if(grep ($_ eq $key, @possible_keys ))
                        {
                            $y_matches{$key} = $value;
                        }
                        else
                        {
                            die "\n\n[COMPADRE] Error: Invalid key before \"=\" for --$opt_name: $val\n" . "Must be one of these: " . join(', ', @possible_keys) . "\n\n";
                        }
                    },
                } 
            ],
            [ 'MT_error_rate=f', "Proportion of the MT sequence that must not match to be called a non-match", { default => 0.01, default_needed => 1 } ],
            [ 'Y_error_rate=f', "Proportion of the Y sequence that must not match to be called a non-match", { default => 0.01, default_needed => 1 } ],   
            [ 'use_mito_match=i', "Ensures that the final reconstructed pedigree complies with the MT matching rules", { default => 0 }],
            [ 'use_y_match=i', "Ensures that the final reconstructed pedigree complies with the y chromosome matching rules", { default => 0 }],  
            [ 'affection_file=s', "Specify path to the file containing the affection status of each sample" ], 
            [ 'affections=s{1,5}', "Like --affection_file; need FILE=[file], optional specification of file columns", {
                test => sub
                    { 
                        my $val = shift;
                        my $opt_name = shift;
                        my @possible_keys = qw(FILE FID IID AFFECTION AFFECTION_VALUE);
                        my ($key,$value) = split(/=/,@_[1]);
                        uc($key);
                        if(grep ($_ eq $key, @possible_keys ))
                        {
                            $affections{$key} = $value;
                        }
                        else
                        {
                            die "\n\n[COMPADRE] Error: Invalid key before \"=\" for --affections option: --$opt_name\n" .
                            "Must be one of these: " . join(', ', @possible_keys) . "\n\n";
                        }
                        return 1;
                    },	
                } 
            ], 
            [ 'int_likelihood_cutoff=f', "Initial minimum likelihood for a relationship to reconstruction", { default => 0.1, default_needed => 1 } ],  

            "ERSA options:",
            [ 'min_cm=f', "minimum segment size to consider", { default => 2.5, default_needed => 1 } ], 
            [ 'max_cm=f', "maximum segment size to consider for estimating the exponential distribution of segment sizes in the population", { default => 10.0, default_needed => 1 } ],
            [ 'max_meioses=i', "maximum number of meioses to consider", { default => 40, default_needed => 1 } ], 
            [ 'rec_per_meioses=f', "expected number of recombination events per meioses", { default => 35.2548101, default => 1 } ],   
            [ 'ascertained_chromosome=s', "Ascertained chromosome", { default => 'no_ascertainment', default_needed => 1 } ],  
            [ 'ascertained_position=i', "chromosomal position of ascertained disease locus" ],
            [ 'control_files=s', "GERMLINE or Beagle fibd output file(s) for population control" ],  
            [ 'control_sample_size=i', "Sample size of control population. Used with --control_files flag" ], 
            [ 'exp_mean=f', "Mean of exp distribution of shared segment size in population", { default => 3.197036753, default_needed => 1} ],
            [ 'pois_mean=f', "Mean of the Poisson distribution of the number of segments shared between a pair of individuals in the population", { default => 13.73, default_needed => 1 } ],
            [ 'pair_file=s', "Restrict pairwise comparisons to the ID pairs specified i this file" ],
            [ 'single_pair', "Restrict pairwise comparisons to the pairs specified in this flag (id1:id2)" ],
            [ 'number_of_ancestors=i', "Restrict relationships to [1] one parent (half-sibs/cousins), [2] two parents (full-sibs/cousins), or [0] (parent-offspring/grandparent-granchild).", { default => 'None' } ], 
            [ 'number_of_chromosomes=i', "Number of chromosomes", { default => 22, default_needed = 1 } ], 
            [ 'parent_offspring_option', "Option to evaluate potential parent-offspring and sibling relationships based on total proportion of the genome that is shared ibd1", { default => 1, default_needed => 1 } ],                                                      
            [ 'parent_offspring_zscore=f', "Z-score for rejecting a sibling relationship in favor of a parent-offspring relationship (alpha=0.01)", { default => 2.33, default_needed => 1 } ], 
            [ 'adjust_pop_dist', "Option to adjust the population distribution of shared segments downward for segments that could not be detected due to recent ancestry", { default => 0, default_needed => 1 } ],  
            [ 'confidence_level=f', "Confidence level for confidence interval around the estimated degree of relationship.", { default => 0.95, default_needed => 1 } ], 
            [ 'mask_common_shared_regions', "excludes chromosomal regions that are commonly shared from evaluation. Used only when the control_files or mask_region_file parameter is specified", { default => 0, default_needed => 1 } ],  
            [ 'mask_region_cross_length=i', "length in base pairs that a shared segment must extend past a masked segment in order to avoid truncation. Used only when mask_common_shared_regions parameter is specified ", { default => 1000000, default_needed => 1 } ],   
            [ 'mask_region_file=s', "file containing chromosomal regions to exclude from evaluation." ], 
            [ 'mask_region_threshold=f', "Threshold for the ratio of observed vs. expected segment sharing in controls before a region will be masked.", { default => 4.0, default_needed => 1 } ], 
            [ 'mask_region_sim_count=i', "number of simulations performed of the null distribution of shared segment locations in controls; results written to output_file.sim", { default => 0, default_needed => 1 } ], 
            [ 'recombination_files=s', "file containing genetic distances for all chromosomes. This parameter must be specified with Beagle fibd input files" ], 
            [ 'beagle_markers_files=s', "Beagle marker files (one file required for each chromosome, wildcards required, ex: chr*beagle.marker). Each filename must begin with the chromosome name followed by a period. This parameter must be specified with Beagle fibd input files" ],     
            [ 'ersa_model_output=s', "ERSA will generate this file when the --model_output_file option specified on the commandline when running ERSA. For each pair of individuals in the dataset, it will contain the likelihood that they are related as 1st through 40th degree relatives sharing 0, 1, or two parents in common at each degree. PRIMUS uses these likelihhods and the likelihood that they are unrelated (obtained from the file specified with the --ersa_results option) to find the most likely way each family network is connected." ],  
            [ 'ersa_results=s', "ERSA's main output file (specified with the --output_file option) needs to be provided. PRIMUS + ERSA uses the likelihood that a pair of individuals are unrelated in its algorithm. Provide the path to the file here." ],
            [ 'ersa_results=s', "ERSA's main output file (specified with the --output_file option) needs to be provided. PRIMUS + ERSA uses the likelihood that a pair of individuals are unrelated in its algorithm. Provide the path to the file here." ],                             
            [ 'project_summary_file=s', "The project level summary file is in the main output directory which we call the Project level directory (default is *_PRIMUS). In the main ouput directory, there are two summary files. One is a summary containing the pairwise relationship table (*_pairwuse_table.txt), and this is NOT the file you want to provide with this option. You want to provide the other summary file." ],             
            [ 'PADRE_multiple_test_correct', "PADRE multiple testing correction", { default => 0, default_needed => 1} ],   
            { argv => \@args }                                                             
        );
    # Handle cases where the user passed the help flag or they passed the man flag
    if ($opt->$help) {
        print $usage->text;
        exit 0;
    }
    # We are going to run the man documentation for the compadre_kickoff.pl 
    if ($opt->man) {
        require Pod::Usage;
        Pod::Usage::pod2usage({ -input => $0, -verbose => 2 }); 
    }

    # Lets valid that the right arguments were parse before constructing the configuration object


    # Initialize your default configuration object
    my $config = {
        global => {
            verbose         =>  $opt->verbose,
            study_name      => "",
            output_dir      => ,
            ersa_data       => $opt->segment_data,
            port_number     => $opt->port_number,
            run_padre       => $opt->run_padre,
            reference_pop   => "",
            log_file        => undef,
            debug           => 0,
            test            => 0,
            cluster         => 0,
            public_html_dir => "",
            usage           => $usage, # Store the usage object right here
        },
        imus => {
            sexes                            => { "FID" => 1, "IID" => 2, "SEX" => 3, "MALE" => 1, "FEMALE" => 2, "FILE" => $opt->sex_file },
            affections                       => { "FID" => 1, "IID" => 2, "AFFECTION" => 3, "AFFECTION_VALUE" => 2, "FILE" => $opt->affection_file },
            ages                             => { "FID" => 1, "IID" => 2, "AGE" => 3, "FILE" => $opt->age_file},
            traits                           => {},
            trait_order                      => [],
            relatedness_threshold            => $opt->rel_threshold,
            degree_rel_cutoff                => $opt->degree_rel_cutoff,
            max_memory                       => $opt->max_memory,
            initial_likelihood_cutoff        => $opt->int_likelihood_cutoff,
            max_generations                  => "none",
            max_generation_gap               => $opt->max_gen_gap,
            missing_data_value               => $opt->missing_val,
            get_max_unrelated_set            => $opt->no_IMUS,
            reconstruct_pedigrees            => $opt->no_PR,
            generate_likelihood_vectors_only => $opt->generate_likelihoods_only,
            ibd_estimates                    => { "FID1" => 1, "IID1" => 2, "FID2" => 3, "IID2" => 4, "IBD0" => 7, "IBD1" => 8, "IBD2" => 9, "PI_HAT" => 10, "FILE" => $opt->plink_ibd},
            mito_matches                     => { "FID1" => 1, "IID1" => 2, "FID2" => 3, "IID2" => 4, "MATCH" => 5, "MATCH_VAL" => 1 },
            y_matches                        => { "FID1" => 1, "IID1" => 2, "FID2" => 3, "IID2" => 4, "MATCH" => 5, "MATCH_VAL" => 1 },
        },
        preprimus => {
            run_prePRIMUS                     => 0,
            data_stem                         => "",
            plink_path                        => "plink",
            no_automatic_IBD                  => $opt->no_automatic_IBD,
            remove_AIMs                       => $opt->remove_AIMs,
            keep_AIMs                         => %opt->keep_AIMs,
            no_mito                           => $opt->no_mito,
            no_y                              => $opt->no_y,
            use_mito_match                    => $opt->use_mito_match,
            use_y_match                       => $opt->use_y_match,
            MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH => $opt->MT_error_rate,
            Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH  => $opt->Y_error_rate,
            internal_ref                      => $opt->internal_ref,
            alt_ref_stem                      => "",
            keep_prePRIMUS_intermediate_files => $opt->keep_inter_files,
            no_PCA_plot                       => $opt->no_PCA_plot,
            ref_pops                          => $opt->ref_pops,
            rerun                             => 0,
            min_pihat_threshold               => $opt->min_pihat_threshold,
            samples_file                      => $opt->samples,
            gzip_genome                       => $opt->gzip_genome,
        },
        padre => {
            run_PRIMUS_plus_ERSA        => 0,
            ersa_model_output           => undef,
            ersa_results                => undef,
            project_summary_file        => undef,
            PADRE_multiple_test_correct => 0,
        },
        ersa => {
            min_cm                      => $opt->min_cm,
            max_cm                      => $opt->max_cm,
            max_meioses                 => $opt->max_meioses,
            rec_per_meioses             => $opt->rec_per_meioses,
            ascertained_chromosome      => $opt->ascertained_chromosome,
            ascertained_position        => $opt->ascertained_position,
            control_files               => $opt->control_files,
            control_sample_size         => $opt->control_sample_size,
            exp_mean                    => $opt->exp_mean,
            pois_mean                   => $opt->pois_mean,
            number_of_ancestors         => $opt->number_of_ancestors,
            number_of_chromosomes       => $opt->number_of_chromosomes,
            parent_offspring_option     => $opt->parent_offspring_option,
            parent_offspring_zscore     => $opt->parent_offspring_zscore,
            adjust_pop_dist             => $opt->adjust_pop_dist,
            confidence_level            => $opt->confidence_level,
            mask_common_shared_regions  => $opt->mask_common_shared_regions,
            mask_region_cross_length    => $opt->mask_region_cross_length,
            mask_region_file            => $opt->mask_region_file,
            mask_region_threshold       => $opt->mask_region_threshold,
            mask_region_sim_count       => $opt->mask_region_sim_count,
            recombination_files         => $opt->recombination_files,
            beagle_markers_files        => $opt->beagle_markers_files,
        }
    };                                                     
}