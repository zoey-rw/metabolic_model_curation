#' Functions to aid in balancing equations within genome-scale metabolic models
#' To be used in conjunction with the Python package memote "find_charge_unbalanced_reactions()"
#' 
#' TO DO: output from "get_bigg_charges" could be presented in a cleaner/more useful way
#' 
#' possibly useful: https://github.com/SBRG/bigg_models/blob/cf16b53ab77ea699a1e132ce10a3eca1690e0aee/bin/load_metanetx



#' DOWNLOAD REFERENCE DATA
#' 
#' Get reference data on names and charges for chemical reactions, metabolites, and cellular compartments from MetanetX and BiGG
#'  Files can be downloaded from Metanetx/BiGG (e.g. using wget) or called directly from site via URL
#' @param reac_xref_path 
#' @param chem_xref_path
#' @param comp_xref_path
#' @param reac_prop_path
#' @param chem_prop_path
#' @param bigg_met_path 
#' @param bigg_rxn_path 
#'
#' @return returns a list with 7 named dataframes for reference
#' @export
#'
#' @examples
#' ref_data <- get_reference_data()
get_reference_data <- function(threads =10, 
															 reac_xref_path = "/projectnb/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/reference_data/reac_xref.tsv", 
														 chem_xref_path = "/projectnb/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/reference_data/chem_xref.tsv",
														 comp_xref_path = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/comp_xref.tsv",
														 reac_prop_path = "/projectnb/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/reference_data/reac_prop.tsv", 
														 chem_prop_path = "/projectnb/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/reference_data/chem_prop.tsv", 
														 bigg_met_path = "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt",
														 bigg_rxn_path = "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt") {
	
	# Cross-referencing tables
	reac_xref <- data.table::fread(reac_xref_path, fill=T, skip=351, nThread = threads) %>% 
		separate(col = "#source", into = c("source","source_id"), sep = ":")
	chem_xref <- data.table::fread(chem_xref_path, fill=T, skip=351, nThread = threads) %>% 
		separate(col = "#source", into = c("source","source_id"), sep = ":")
	comp_xref <- data.table::fread(comp_xref_path, fill=T, skip=351, nThread = threads) %>% 
		separate(col = "#source", into = c("source","source_id"), sep = ":")
	
	# Property tables
	reac_prop <- data.table::fread(reac_prop_path, fill=T, skip=351, nThread = threads)
	chem_prop <- data.table::fread(chem_prop_path, fill=T, skip=351, nThread = threads)
	
	# BiGG data from url
	bigg_met <- data.table::fread(bigg_met_path, fill=T)
	bigg_rxn <- data.table::fread(bigg_rxn_path, fill=T) %>% 
		separate(reaction_string, sep = " <-> ", into = c("left","right"), remove = F,fill = "right")
	
	out <- list("chem_xref" = chem_xref,
							"comp_xref" = comp_xref,
							"reac_xref" = reac_xref,
							"reac_prop" = reac_prop,
							"chem_prop" = chem_prop,
							"bigg_met" = bigg_met,
							"bigg_rxn" = bigg_rxn) 
	return(out)
}

# TO DO: could be vectorized to be faster
get_bigg_met <- function(metanetxID = "MNXM3453",
														 reference_data = NULL) {
	chem_xref = reference_data$chem_xref
	bigg_met = reference_data$bigg_met
	
	chem_met_info <- chem_xref %>% filter(ID == !!metanetxID & source == "biggM" & !grepl("M_", source_id)) %>% select(source_id)
	bigg_met_ID <- unique(chem_met_info$source_id)
	return(bigg_met_ID)
	
}

get_bigg_met(metanetxID = "MNXM3453", reference_data = ref_data)

#' INVESTIGATE SPECIFIC REACTIONS
#' 
#' Get data on names and charges for specific chemical reaction from MetanetX and BiGG
#' @param metanetxID Optional - ID of problematic/unbalanced chemical reaction
#' @param biggID Optional - ID of problematic/unbalanced chemical reaction
#' @param reference_data
#'
#' @return returns series of messages describing the metabolite(s) missing charges in a reaction, and the typical charges for that metabolite in BiGG models
#' @export
#'
#' @examples
#' get_bigg_charges(metanetxID = "MNXR94775",
#' 								 reference_data = ref_data) 
#' 
#' get_bigg_charges(metanetxID = NULL, 
#' 								 biggID = "2DDARAA",
#' 								 reference_data = ref_data)
get_bigg_charges <- function(metanetxID = NULL, 
														 biggID = NULL,
														 reference_data = NULL) {
	
	pacman::p_load(tidyverse, minval, curl) 
	
	if (is.null(reference_data)) {
		message("First run get_reference_data() and pass the output as the argument to reference_data")
	} else {
		chem_xref = reference_data$chem_xref
	comp_xref = reference_data$comp_xref
	reac_xref = reference_data$reac_xref
	reac_prop = reference_data$reac_prop
	chem_prop = reference_data$chem_prop
	bigg_met = reference_data$bigg_met
	bigg_rxn = reference_data$bigg_rxn
	}
	
	# Get metanetID from BiGG ID, if necessary
	if (!is.null(biggID)) {
		rxn_info <- reac_xref %>% filter(source_id == !!biggID) 
		metanetxID <- unique(rxn_info$ID)
	}
	
	rxn <- reac_prop %>% filter(`#ID` == !!metanetxID)
	
	rxn_split <- strsplit(rxn$mnx_equation, " = ") %>% unlist()
	
	#side <- "left"
	for (side in c("left","right")) {
		rxn_split <- strsplit(rxn$mnx_equation, " = ") %>% unlist()
		rxn_side <- switch(side,
											 "left" = rxn_split[[1]],
											 "right" = rxn_split[[2]])
		
		met <- minval::metabolites(rxn_side) %>% 
			lapply(., function(x) gsub("@MNX[CD][0-9]","", x)) %>% 
			unlist()
		
		compartments <- minval::metabolites(rxn_side) %>% 
			sapply(., function(x) {strsplit(x, "@") %>% 
					sapply(tail, 1 )})
		names(compartments) <- met
		
		init_info <- list()
		
		for (met in met) {
			db_entry <-  chem_prop %>% filter(`#ID` == met)
			init_info[[met]] <- list("charge" = db_entry$charge,
															 "mass" = db_entry$mass)
		}
		init_charge <- lapply(init_info, "[[", 1)
		init_mass <- lapply(init_info, "[[",2)
		
		missing_charges <- init_charge[is.na(init_charge)]
		bigg_charges <- init_charge
		names(bigg_charges) <- lapply(names(bigg_charges), function(x) get_bigg_met(metanetxID = x, reference_data = reference_data)) %>% unlist()
		
		if (length(missing_charges)==0) {
			message("MetanetX has the following charges on ", side, " of equation: ")
			print(unlist(bigg_charges))
		}
		
		#met <- "MNXM3453" # for testing
		for (met in names(missing_charges)){
			bigg_met_id <- get_bigg_met(metanetxID = met, reference_data = reference_data)
			
			message("MetanetX is missing charge on ", side, " side of equation for: \nMetanetX metabolite: ", met, "\nBiGG metabolite: ", bigg_met_id)
			

			# Get list of BiGG reactions that involve this metabolite on the same side
			bigg_ids <- chem_xref %>% filter(ID == met) %>% select(source_id) %>% unlist()
			bigg_universal_id <- bigg_met %>% filter(universal_bigg_id %in% bigg_ids) %>% select(universal_bigg_id) %>% unique() %>% unlist()
			if (side == "left") {
			rxn_list <- bigg_rxn %>% filter(grepl(bigg_universal_id, left)) %>% select(bigg_id) %>% unlist()
			} else {
				rxn_list <- bigg_rxn %>% filter(grepl(bigg_universal_id, right)) %>% select(bigg_id) %>% unlist()
			}
			
			#rxn_of_interest <- "2AGPGAT161"
			rxn_of_interest <- "2AGPG161tipp"
			for (rxn_of_interest in rxn_list) {
				
				message("Metabolite used in reaction: ", rxn_of_interest)
				rxn_info <- bigg_rxn %>% filter(bigg_id == !!rxn_of_interest) 
				
				model_list <- rxn_info$model_list %>% strsplit("; ") %>% unlist()
				api_url <- paste0("http://bigg.ucsd.edu/api/v2/models/", model_list[[1]], "/reactions/", rxn_of_interest)
				
				#if (RCurl::url.exists(api_url)) {
				json_out <- jsonlite::fromJSON(api_url, flatten = T)
				charges <- json_out$metabolites %>% filter(bigg_id == !!bigg_universal_id) 
				
				# Get metabolite compartment for decision-making
				comp <- compartments[met]
				comp_desc <- comp_xref[comp_xref$ID==comp,]$description

				
				# Report output
				for (k in 1:nrow(charges)){
					message("Observed charge for: ", met, " is ", charges[k,]$stoichiometry, " in ", "compartment: ", charges[k,]$compartment_bigg_id)
				}
			}
			# After the final loop, report the actual compartment to decide which of the charges is most appropriate
			message("Metabolite ", met, " is in compartment: ", comp, ", which has the following description from MetanetX: ", comp_desc)
		}
	}
}



ref_data <- get_reference_data()

get_bigg_charges(metanetxID = "MNXR94775",
								 reference_data = ref_data) 

get_bigg_charges(metanetxID = NULL, 
								 biggID = "2AGPGAT160",
								 reference_data = ref_data)




get_bigg_charges(metanetxID = "MNXR94774", 
								 reference_data = ref_data)
