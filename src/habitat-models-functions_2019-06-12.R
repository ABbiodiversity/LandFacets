#
# Title: Functions to execute species habitat models
# Created: June 13th, 2019
# Last Updated: August 16th, 2019
# Author: Brandon Allen
# Objectives: Soil and land facet habitat model functions
# Keywords: Southern models, Grid predictions, Site predictions, AUC calculation
#

###################
# Southern models # 
###################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

southern_models <- function (data.analysis, results.store, landscape.models, prediction.matrix, species.ID) {
        
        # Convert data to probabilities based on the # of quadrants sampled
        data.analysis["pcount"] <- data.analysis[, species.ID] / data.analysis$nQuadrant
        
        # Creating lists for storage of all models
        space.climate.store <- landscape.store <- list(NULL) 
        
        for (model in 1:length(landscape.models)) {
                
                landscape.store[[model]] <- try(glm(landscape.models[[model]], family = "binomial", data = data.analysis, weights = visit))
                
        }
        
        nModels <- length(landscape.store)
        
        # AIC calculation  (I'm using AIC here, because this is primarily for prediction, rather than finding a minimial best model)
        aic.ta <- rep(999999999, (nModels))
        
        for (i in 1:(nModels)) {
                
                if (!is.null(landscape.store[[i]]) & class(landscape.store[[i]])[1] != "try-error") {  # last part is to not used non-converged models, unless none converged
                        aic.ta[i] <- AICc(landscape.store[[i]])
                }
                
        }
        
        aic.delta <- aic.ta - min(aic.ta)
        aic.exp <- exp(-1 / 2 * aic.delta)
        aic.wt.ta <- aic.exp / sum(aic.exp)
        
        # Store coefficients
        # Prediction Matrix for 100% of each type of HF/Veg
        
        p1 <- p1.se <- array(0, c(nModels, nrow(pm))) # Predictions for each model for each soil and HF type type.  These are the main coefficients
        p.site1 <- array(0, c(nModels, nrow(data.analysis)))  # Predictions for each model and each site.  These are used as the offsets in the next stage
        p1.aspen <- p1.se.aspen <- rep(0, nModels)  # Predictions for each model for aspen coefficient
        
        for (i in 1:nModels) {
                
                if (class(landscape.store[[i]])[1] != "try-error") {   # Prediction is 0 if model failed, but this is not used because AIC wt would equal 0
                        
                        p <- predict(landscape.store[[i]], newdata = data.frame(prediction.matrix, paspen = 0), se.fit = TRUE)  # For each type.  Predictions made at 0% Aspen.  All predictions made with new protocol.  Aspen effect added later, and plotted as separate points
                        p1[i,] <- p$fit
                        p1.se[i,] <- p$se.fit
                        p.site1[i,] <- predict(landscape.store[[i]])  # For each site (using the original d.sp data frame)
                        if ("paspen" %in% names(coef(landscape.store[[i]]))) {
                                p1.aspen[i] <- coefficients(summary(landscape.store[[i]]))["paspen","Estimate"]
                                p1.se.aspen[i] <- coefficients(summary(landscape.store[[i]]))["paspen","Std. Error"]
                        }
                        
                        
                }
                
        }
        
        tTypeMean <- tTypeVar <- NULL  # Logit-scaled model averages for each veg type
        
        aic.wt.ta.adj <- ifelse(aic.wt.ta < 0.01, 0, aic.wt.ta)  # Use adjusted weight to avoid poor fitting models with extreme values for certain types (usually numerically equivalent to 0's or 1's)
        aic.wt.ta.adj <- aic.wt.ta.adj / sum(aic.wt.ta.adj)
        
        tpAspenMean <- sum(aic.wt.ta * p1.aspen)
        tpAspenVar <- sum(aic.wt.ta * sqrt(p1.se.aspen^2 + (rep(tpAspenMean, nModels) - p1.aspen)^2))^2  # AIC-weighted variance of mean...
        
        for (i in 1:nrow(prediction.matrix)) {
                
                tTypeMean[i] <- sum(aic.wt.ta.adj * p1[, i])  # Mean of nModels for each veg type, on transformed scale
                tTypeVar[i] <- sum(aic.wt.ta.adj * sqrt(p1.se[, i]^2+(rep(tTypeMean[i], nModels) - p1[, i])^2))^2  # AIC-weighted variance of mean...
                
        }
        
        names(tTypeMean) <- names(tTypeVar)<- colnames(prediction.matrix) # Rename the final intercept coefficient to productive
        
        data.analysis$prediction <- colSums(p.site1 * aic.wt.ta.adj)  # Add logit-scaled model average prediction to data frame
        
        # Put coefficients in Coef matrix (and SE's)
        results.store$landscape.coef[species.ID, na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))] <- plogis(tTypeMean[colnames(results.store$landscape.coef)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.coef)))]])  # On ordinal scale
        results.store$landscape.se[species.ID, na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))] <- sqrt(tTypeVar[colnames(results.store$landscape.se)[na.omit(match(names(tTypeMean), colnames(results.store$landscape.se)))]])  # On logit scale
        results.store$paspen.coef[species.ID, "paspen"] <- tpAspenMean
        results.store$paspen.se[species.ID, "paspen"] <- sqrt(tpAspenVar)
        
        # Do Not Need To Adjust! HardLin is a classification and is combined with Urban Industrial as Alienating footprint
        # 3. Do the adjustments to coefficients
        # # 3.2 Variance-weighted average of hard linear (often poorly estimated) with urban/industrial (no early seral to adjust soft linear in south)
        results.store$landscape.coef[species.ID, "HardLin"] <- plogis( (qlogis(results.store$landscape.coef[species.ID, "HardLin"]) / results.store$landscape.se[species.ID, "HardLin"]^2 + qlogis(results.store$landscape.coef[species.ID, "UrbInd"])/results.store$landscape.se[species.ID, "UrbInd"]^2) / (1 / results.store$landscape.se[species.ID, "HardLin"]^2 + 1 / results.store$landscape.se[species.ID, "UrbInd"]^2) )  # Inverse-variance weighted, done on logit scale then converted back
        results.store$landscape.se[species.ID, "HardLin"] <- sqrt(1 / (1 / results.store$landscape.se[species.ID, "HardLin"]^2 + 1 / results.store$landscape.se[species.ID, "UrbInd"]^2))
        
        # 4. Residual variation due to surrounding HF 
        # This does climate and space covariates first, so that best model from those can be used in main models of residual surrounding-HF effects (and so that surrounding HF models can be model-averaged)
        # 4.1 First, find best climate and/or spatial covariates.  
        # Climate and spatial variables sets are tried, with surrounding THD and THD^2 as initial HF covariates
        p1 <- colSums(results.store$landscape.coef[species.ID, ] * t(data.analysis[ ,colnames(results.store$landscape.coef)]))  # Multiply coefficients by soil type proportion for each site to get offset
        
        p1 <- plogis(qlogis(p1) + results.store$paspen.coef[species.ID, ] * data.analysis$paspen)  # And add in effect of pAspen on logit scale
        
        wt1 <- data.analysis$nQuadrant * data.analysis$visit  # Need to do this to use model.matrix below for plotting
        space.climate.store[[1]] <- try(glm(pcount ~ offset(qlogis(0.998 * p1 + 0.001)), data = data.analysis, family = "binomial", weights = wt1))  # Prediction can be 0 (water), so need to scale a bit to avoid infinities in logit offset. Protocol not included here, because already accounted for in offset values 
        space.climate.store[[2]] <- try(update(space.climate.store[[1]], .~.+ PET))
        space.climate.store[[3]] <- try(update(space.climate.store[[1]], .~.+ AHM))
        space.climate.store[[4]] <- try(update(space.climate.store[[1]], .~.+ MAT))
        space.climate.store[[5]] <- try(update(space.climate.store[[1]], .~.+ FFP))
        space.climate.store[[6]] <- try(update(space.climate.store[[1]], .~.+ MAP + FFP))
        space.climate.store[[7]] <- try(update(space.climate.store[[1]], .~.+ MAP +FFP + MAP:FFP))
        space.climate.store[[8]] <- try(update(space.climate.store[[1]], .~.+ MAT + MAP + PET + AHM))
        space.climate.store[[9]] <- try(update(space.climate.store[[1]], .~.+ MAT + MAP + PET + AHM + MAP:PET + MAT:AHM))
        space.climate.store[[10]] <- try(update(space.climate.store[[1]] ,.~.+ MAT + MAP))
        space.climate.store[[11]] <- try(update(space.climate.store[[1]] ,.~.+ MWMT + MCMT))
        space.climate.store[[12]] <- try(update(space.climate.store[[1]] ,.~.+ AHM + PET))
        space.climate.store[[13]] <- try(update(space.climate.store[[1]] ,.~.+ MAT + I(MAT^2) + MWMT + I(MWMT^2)))
        space.climate.store[[14]] <- try(update(space.climate.store[[1]] ,.~.+ MWMT + MCMT + FFP + MAT))
        for (i in 1:14) space.climate.store[[i + 14]] <- try(update(space.climate.store[[i]], .~.+ Lat + Long + Lat:Long)) 
        for (i in 1:14) space.climate.store[[i + 28]] <- try(update(space.climate.store[[i]], .~.+ Lat + Long + Lat:Long + I(Lat^2) + I(Long^2)))
        
        nModels.sc<-length(space.climate.store)
        # BIC calculation to select best covariate set Uses BIC for more conservative variable set
        bic.sc <- rep(999999999, nModels.sc)
        for (i in 1:(nModels.sc)) {
                
                if (!is.null(space.climate.store[[i]]) & class(space.climate.store[[i]])[1] != "try-error") {
                        bic.sc[i] <- BIC(space.climate.store[[i]])
                }
                
        }
        
        bic.delta.sc <- bic.sc - min(bic.sc)
        bic.exp.sc <- exp(-1 / 2 * bic.delta.sc)
        bic.wt.sc <- bic.exp.sc / sum(bic.exp.sc)
        best.model.sc <- space.climate.store[[which.max(bic.wt.sc)]]
        
        # And the subset of climate and/or spatial variables
        vnames <- names(coef(space.climate.store[[which.max(bic.wt.sc)]]))
        vnames1 <- ifelse(vnames == "(Intercept)", "Intercept", vnames)  # This is for the names used in km2.res (where "(", "^", etc. can't be used)
        vnames1 <- ifelse(vnames1 == "I(Lat^2)", "Lat2", vnames1)
        vnames1 <- ifelse(vnames1 == "I(Long^2)", "Long2", vnames1)
        vnames1 <- ifelse(vnames1 == "I(MAT^2)", "MAT2", vnames1)
        vnames1 <- ifelse(vnames1 == "I(MWMT^2)", "MWMT2", vnames1)
        vnames1 <- gsub(":", "", vnames1)
        
        # Storing of climate coefficients fails
        coef.match <- colnames(results.store$climate.coef) %in% vnames1
        coef.match <- match(colnames(results.store$climate.coef), vnames1, nomatch = 0)
        results.store$climate.coef[species.ID, coef.match != 0] <- coef(best.model.sc)[coef.match]
        
        return(results.store)
}


####################
# Site predictions # 
####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

site_predictions <- function(model.coef, landscape.data, climate.data, paspen.data, species.data, species.ID) {

        veg.prediction <- colSums(model.coef$landscape.coef[species.ID, ]*t(landscape.data[ ,colnames(model.coef$landscape.coef)])) # Prediction based on veg types only.
        
        # Calculate cimate and spatial prediction for each species
        
        climate.prediction <- colSums(model.coef$climate.coef[species.ID, c(FALSE, as.logical(!is.na(model.coef$climate.coef[species.ID, -1])))] * t(climate.data[ ,colnames(model.coef$climate.coef)[c(FALSE, as.logical(!is.na(model.coef$climate.coef[species.ID, -1])))]])) # Prediction of residual (climate and spatial) effect
        
        # Take the difference between unqiue species intercept and the predicted climate residuals
        
        climate.prediction <- climate.prediction + model.coef$climate.coef[species.ID, 1]
        
        # Calculate the pAspen prediction for each species, same for both current and reference condition.
        
        aspen.prediction <- colSums(model.coef$paspen.coef[species.ID, "paspen"] * t(paspen.data))
        
        # Combine the soil and climate/space based predictions
        
        site.prediction <- plogis(qlogis(0.998*veg.prediction + 0.001) + climate.prediction + aspen.prediction) # Apply same transformation as used in fitting residual model
        
        auc.predict <- rep(site.prediction, species.data$nQuadrant)
        
        auc.obs <- NULL
        
        for (i in 1:nrow(species.data)) {
                
                if (species.data[i, species.ID] == 0) {
                        
                        occ.temp <- rep(0, species.data[i, "nQuadrant"])
                        
                }else {
                        
                        occ.temp <- c(rep(1, species.data[i, species.ID]), rep(0, (species.data[i, "nQuadrant"] - species.data[i, species.ID])))
                }
                
                auc.obs <- c(auc.obs, occ.temp)
                
                rm(occ.temp)
        }

        # Storage of matrix required for AUC
        predict.out <- matrix(c(auc.obs, auc.predict), nrow = length(auc.predict), 
                              ncol = 2, byrow = FALSE, 
                              dimnames = list(NULL, c(paste(species.ID, "obs", sep = "_"), paste(species.ID, "pred", sep = "_"))))
        
        return(predict.out)
        
        rm(auc.obs, auc.predict, auc.out, veg.prediction, climate.prediction, aspen.prediction, site.prediction)
        
}

#########################
# Landscape predictions # 
#########################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

landscape_only_predictions <- function(model.coef, landscape.data, species.data, species.ID) {
        
        veg.prediction <- colSums(model.coef$landscape.coef[species.ID, match(colnames(model.coef$landscape.coef), colnames(landscape.data))[!is.na(match(colnames(model.coef$landscape.coef), colnames(landscape.data)))]]*t(landscape.data[ , match(colnames(landscape.data), colnames(model.coef$landscape.coef))])) # Prediction based on veg types only.
        
        auc.predict <- rep(veg.prediction, species.data$nQuadrant)
        
        auc.obs <- NULL
        
        for (i in 1:nrow(species.data)) {
                
                if (species.data[i, species.ID] == 0) {
                        
                        occ.temp <- rep(0, species.data[i, "nQuadrant"])
                        
                }else {
                        
                        occ.temp <- c(rep(1, species.data[i, species.ID]), rep(0, (species.data[i, "nQuadrant"] - species.data[i, species.ID])))
                }
                
                auc.obs <- c(auc.obs, occ.temp)
                
                rm(occ.temp)
        }
        
        # Storage of matrix required for AUC
        predict.out <- matrix(c(auc.obs, auc.predict), nrow = length(auc.predict), 
                              ncol = 2, byrow = FALSE, 
                              dimnames = list(NULL, c(paste(species.ID, "obs", sep = "_"), paste(species.ID, "pred", sep = "_"))))
        
        return(predict.out)
        
        rm(auc.obs, auc.predict, auc.out, veg.prediction)
        
}


####################
# Grid predictions # 
####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grid_predictions <- function(model.coef, landscape.data, climate.data, paspen.data, species.ID) {
        
        veg.prediction <- colSums(model.coef$landscape.coef[species.ID, match(colnames(model.coef$landscape.coef), colnames(landscape.data))[!is.na(match(colnames(model.coef$landscape.coef), colnames(landscape.data)))]]*t(landscape.data[ , match(colnames(landscape.data), colnames(model.coef$landscape.coef))])) # Prediction based on veg types only.
        
        # Calculate cimate and spatial prediction for each species
        
        climate.prediction <- colSums(model.coef$climate.coef[species.ID, c(FALSE, as.logical(!is.na(model.coef$climate.coef[species.ID, -1])))] * t(climate.data[ ,colnames(model.coef$climate.coef)[c(FALSE, as.logical(!is.na(model.coef$climate.coef[species.ID, -1])))]])) # Prediction of residual (climate and spatial) effect
        
        # Take the difference between unqiue species intercept and the predicted climate residuals
        
        climate.prediction <- climate.prediction + model.coef$climate.coef[species.ID, 1]
        
        # Calculate the pAspen prediction for each species, same for both current and reference condition.
        
        aspen.prediction <- colSums(model.coef$paspen.coef[species.ID, "paspen"] * t(paspen.data))
        
        # Combine the soil and climate/space based predictions
        
        site.prediction <- plogis(qlogis(0.998*veg.prediction + 0.001) + climate.prediction + aspen.prediction) # Apply same transformation as used in fitting residual model
        
        return(site.prediction)
        
        rm(veg.prediction, climate.prediction, aspen.prediction, site.prediction)
        
}

###################
# AUC Calculation # 
###################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model_auc <- function(soil.matrix, facet.matrix, obs, pred, species.ID) {
        
        # Calcualte AUC of both the soil and facet matrix
        auc.1 <- roc(facet.matrix[, obs], facet.matrix[, pred])$auc
        auc.2 <- roc(soil.matrix[, obs], soil.matrix[, pred])$auc
        return(matrix(data = c(auc.1,
                               auc.2,
                               roc.test(auc.1, auc.2)$p.value),
                      nrow = 1, 
                      ncol = 3, 
                      dimnames = list(species.ID, c("facetAUC", "soilAUC", "pvalue"))))
        
}