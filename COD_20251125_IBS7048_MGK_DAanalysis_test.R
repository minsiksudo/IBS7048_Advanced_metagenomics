
# Before running DA analysis, you need to double check the distribution of your data

set.seed(seed)


# I did CLR transform with `microbiome` package ---------------------------------

my_phyloseq_da <- my_phyloseq %>% # set your file here
        subset_taxa(., Kingdom %in% "Bacteria") %>%
        prune_taxa(taxa_sums(.) > 0, .) %>%
        microbiome::transform(., 'clr')

#dummy <- compositions::clr(otu_table(phyloseq_metaphlan$phy_microbe_rel)) %>% data.frame()
#dummy2 <- otu_table(physeq_mic_clr_school) %>% data.frame()

# LME4 for CLR normalized data --------------------------------------------

lmer_function <- function(phyloseq_lmer) {
        
        lmer_taxa_rel_clr <- data.frame()
        
        all <- phyloseq_lmer %>% taxa_sums() %>% length()
        
        for(i in 1:all) {
                #Creating a data frame that includes CLR transformed data of i-th bug.
                #making differnt otu tables for each sample type
                otu_table <- phyloseq_lmer %>% otu_table
                
                #tax table for different sample type
                taxa_data <- otu_table[,i] %>% data.frame()
                
                #Making a merged dataframe having sample data and CLR transformed output
                lme_data <- merge(taxa_data, sample_data(phyloseq_lmer), by = 0) %>% column_to_rownames("Row.names")
                
                #generating a character of formula.
                #Here, as the taxa are already CLR transformed, I did not make model at log-scale.
                lme4_formula <- paste(names(taxa_data), "~", "your model") # <-- modify your model here
                
                
                class_result <- tryCatch({
                        model <- lmerTest::lmer(formula = lme4_formula, data = lme_data)
                        summary(model) %>%
                                .$coefficients %>%
                                data.frame(check.names = FALSE) %>%
                                mutate(feature = colnames(otu_table)[i]) %>%
                                rownames_to_column("value") %>%
                                subset(value %in% c("your variable 1", "your variable 2"))  # <-- modify based on your model here
                } %>% suppressMessages(), 
                error = function(e) {
                        message(paste("Error in model for taxon", i, ":", e$message))
                        return(NULL)
                }) 
                
                #row binding all the associations of i-th taxa to one data frame
                lmer_taxa_rel_clr <- rbind(lmer_taxa_rel_clr,
                                                 class_result) %>%
                        remove_rownames()
                
                
        } 
        lmer_taxa_rel_clr
}


da_result <- lmer_function(my_physeq_da)# set your file here


da_result <- da_result %>%
        mutate(qval = qvalue::qvalue(`Pr(>|t|)`, lambda = 0)$qvalue,
               p_bh = p.adjust(p = `Pr(>|t|)`, method = "BH")) %>%
        suppressMessages()
