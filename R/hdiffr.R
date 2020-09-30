#' R functions for heterogeneous diffusion model
#'
#' hdiffr performs an event history diffusion model that is modified in order to
#' account for heterogeneity in the diffusion process. In addition to modeling
#' chracteristics that influence the direct propensity of a focal actor to adopt,
#' hdiffr can also model characteristics of the focal actor that makes it more
#' susceptible to the prior adoptions of other actors, as well as characteristics
#' of prior-adopting actors that make them more or less influential on others, and
#' characteristics of the social structure among actors that influence the transmission
#' of influence between actors.
#' @param data data with SS (startstate), ES (endstate), ET (endtime), ST (starttime) columns (required)
#' @param xvars propensity varlist
#' @param vvars susceptibility varlist
#' @param wvars infectiousness varlist
#' @param zdvars distance-based proximity varlist
#' @param zgvars group-based proximity varlist
#' @param zmvar data for matrix-based proximity
#' @param zmvar2 second data for matrix-based proximity
#' @param idvar focal ID
#' @param multispell option for multiple time spells per observation
#' @param repeatevents to use with multispell to allow IDs to be at-risk of adoption repeatedly
#' @param repeatvar to use with multispell and reveatvars to add a coefficient testing influence of repeat events in infectiousness
#' @param vintercept to test simple contagion
#' @param robust cluster errors based on ID
#' @param hrno to display output using coefficients rather than hazard ratios
#' @param dist to specify distribution for survreg (default is exponential)
#'
#' @return survival model results
#' @export
#' @examples
#' library(hdiffr)
#' data(exampleData)
#' result <- hdiffr(data = data, xvars = c('lnsale', 'roa', 'vote'),
#'       vvars = 'activistgrp', wvars = 'lnsale', zgvars = 'sic',
#'       idvar = 'gvkey', multispell = 1, vintercept = 1, hrno = 1)
hdiffr <- function(data, xvars, vvars = NULL, wvars = NULL, zdvars = NULL, zgvars = NULL, zmvar = NULL, zmvar2 = NULL, idvar, multispell = 0,
    repeatevents = 0, repeatvar = 0, vintercept = 1, robust = 0, hrno = 1, dist = "exponential") {

    data_wf <- data[, c("ES", "ET", "SS", "ST", idvar, xvars, vvars, wvars, zdvars, zgvars)]

    data_wf$all_case <- 1

    data_wf$dups <- duplicated(data_wf[[idvar]])
    spellmax <- max(aggregate(data_wf$dups, by = list(idvar = data_wf[[idvar]]), FUN = sum)$x)
    if (spellmax > 0) {
        warning("Warning: More than one observation per ID detected.")
    }
    data_wf <- within(data_wf, rm(dups))

    dups2 <- duplicated(data_wf[[idvar]][which(data_wf$ES == 1)])
    data_wf$dups2 <- NA
    data_wf$dups2[which(data_wf$ES == 1)] <- dups2
    eventmax <- max(aggregate(data_wf$dups2, by = list(idvar = data_wf[[idvar]]), FUN = sum, na.rm = T)$x)
    if (eventmax > 0) {
        warning("Warning: More than one observation per ID detected.")
    }
    data_wf <- within(data_wf, rm(dups2))

    # save(data_wf, file = 'working_file.RData')


    data_gi <- data_wf[which(data_wf$ES == 1), ]
    colnames(data_gi)[which(colnames(data_gi) == "ET")] <- "ET_j"
    colnames(data_gi)[which(colnames(data_gi) == idvar)] <- paste(idvar, "j", sep = "_")

    if (exists("wvars") & !is.null(wvars)) {
        for (i in 1:length(wvars)) {
            data_gi[[paste(wvars[i], "j", sep = "_")]] <- data_gi[[wvars[i]]]
        }
    }
    if (exists("zdvars") & !is.null(zdvars)) {
        for (i in 1:length(zdvars)) {
            data_gi[[paste(zdvars[i], "j", sep = "_")]] <- data_gi[[zdvars[i]]]
        }
    }
    if (exists("zgvars") & !is.null(zgvars)) {
        for (i in 1:length(zgvars)) {
            data_gi[[paste(zgvars[i], "j", sep = "_")]] <- data_gi[[zgvars[i]]]
        }
    }
    gi_include <- grep("[(_j$)]", colnames(data_gi))
    data_gi_2 <- cbind(all_case_2 = data_gi[, c("all_case")], data_gi[, c(gi_include)])

    data_gi_2 <- data_gi_2[order(data_gi_2$ET_j), ]

    data_gi_2$et_max <- max(data_gi_2$ET_j)

    event_dupl <- aggregate(all_case ~ ET_j, data_gi_2, FUN = "length")
    colnames(event_dupl)[2] <- "event_dupl"

    data_gi_2 <- merge(x = data_gi_2, y = event_dupl, by = "ET_j", all.x = TRUE)

    # save(data_gi_2, file = 'general_influence.RData')

    # load(file = 'working_file.RData')

    data_wf <- data_wf[order(data_wf[[idvar]]), ]

    if (exists("multispell") & multispell == 1) {
        data_wf <- data_wf[order(data_wf[[idvar]], data_wf$ST), ]
        data_wf$i_spell <- with(data_wf, ave(rep(1, nrow(data_wf)), eval(as.symbol(idvar)), FUN = seq_along))
    }

    data_wf_gi <- merge(x = data_wf, y = data_gi_2, by = "all_case", all = TRUE)

    # save(data_wf_gi, file = 'working_file.RData')

    data_wf_gi <- data_wf_gi[order(data_wf_gi[[idvar]], data_wf_gi$ST, data_wf_gi$ET_j, data_wf_gi[[paste(idvar, "j", sep = "_")]]), ]

    data_wf_gi$seqnum_dupl <- c(1:dim(data_wf_gi)[1])

    # save(data_wf_gi, file = 'tempfile_1.RData')

    if (exists("multispell") & multispell == 1) {
        data_wf_gi$sp_lose <- 1
        data_wf_gi$sp_lose[which(data_wf_gi$ET_j <= data_wf_gi$ST)] <- 0
        data_wf_gi$sp_lose[which(data_wf_gi$ET_j > data_wf_gi$ET)] <- 0
        sp_notlost <- aggregate(data_wf_gi$sp_lose, by = list(idvar = data_wf_gi[[idvar]], i_spell = data_wf_gi$i_spell), FUN = sum)
        colnames(sp_notlost)[which(colnames(sp_notlost) == "idvar")] <- eval(substitute(idvar))
        colnames(sp_notlost)[which(colnames(sp_notlost) == "x")] <- "sp_notlost"
        data_wf_gi <- merge(x = data_wf_gi, y = sp_notlost, by = c(eval(substitute(idvar)), "i_spell"), all.x = TRUE)
        data_wf_gi <- data_wf_gi[order(data_wf_gi[[idvar]], data_wf_gi$i_spell, data_wf_gi$ST, data_wf_gi$ET_j, data_wf_gi[[paste(idvar,
            "j", sep = "_")]]), ]
        # if i-spell has no j-events: drop all but first event within spell
        if (length(which(data_wf_gi$sp_notlost == 0 & c(FALSE, diff(data_wf_gi[[idvar]]) == 0) & c(FALSE, diff(data_wf_gi$i_spell) == 0))) >
            0) {
            data_wf_gi <- data_wf_gi[-which(data_wf_gi$sp_notlost == 0 & c(FALSE, diff(data_wf_gi[[idvar]]) == 0) & c(FALSE, diff(data_wf_gi$i_spell) ==
                0)), ]
        }
        # if i-spell has j-events in it: drop too-early & too-late events within spell
        data_wf_gi <- data_wf_gi[-which(data_wf_gi$sp_lose == 0 & data_wf_gi$sp_notlost > 0), ]
        # if i-spell has j-events in it (but not at exact end of i-spell): add spell from last j-event to i-spell endtime but prevent that added
        # spell from counting toward adoption
        if (length(which(c(diff(data_wf_gi[[idvar]]) == 0, FALSE) & c(diff(data_wf_gi$i_spell) != 0, TRUE) & data_wf_gi$sp_notlost > 0 &
            data_wf_gi$ET_j < data_wf_gi$ET)) > 0) {
            data_wf_gi <- rbind(data_wf_gi, data_wf_gi[which(c(diff(data_wf_gi[[idvar]]) == 0, FALSE) & c(diff(data_wf_gi$i_spell) != 0,
                TRUE) & data_wf_gi$sp_notlost > 0 & data_wf_gi$ET_j < data_wf_gi$ET), ])

        }
        data_wf_gi <- data_wf_gi[order(data_wf_gi[[idvar]], data_wf_gi$i_spell, data_wf_gi$ST, data_wf_gi$ET_j, data_wf_gi[[paste(idvar,
            "j", sep = "_")]]), ]

        data_wf_gi$all_case[which(c(diff(data_wf_gi[[idvar]]) == 0, FALSE) & c(diff(data_wf_gi$i_spell) != 0, TRUE) & c(FALSE, diff(data_wf_gi$seqnum_dupl) ==
            0))] <- 0
    }

    # first spell for ID starts with overall i's starttime and other spells start at end of last spell (last adoption)
    require(dplyr)
    spell_st <- data.frame(data_wf_gi %>% group_by(eval(as.symbol(idvar))) %>% mutate(V1 = dplyr::lag(ET_j, n = 1, default = NA)), stringsAsFactors = FALSE)
    data_wf_gi <- cbind(data_wf_gi, spell_st = spell_st$V1)
    data_wf_gi$spell_st[which(c(TRUE, diff(data_wf_gi[[idvar]]) != 0))] <- data_wf_gi$ST[which(c(TRUE, diff(data_wf_gi[[idvar]]) != 0))]
    # multispell data first j-event in i-spell takes i-spell's starttime
    if (exists("multispell") & multispell == 1) {
        data_wf_gi <- data_wf_gi[order(data_wf_gi[[idvar]], data_wf_gi$i_spell, data_wf_gi$ST, data_wf_gi$ET_j, data_wf_gi[[paste(idvar,
            "j", sep = "_")]]), ]
        data_wf_gi$spell_st[which(c(FALSE, diff(data_wf_gi[[idvar]]) == 0) & c(TRUE, diff(data_wf_gi$i_spell) != 0))] <- data_wf_gi$ST[which(c(FALSE,
            diff(data_wf_gi[[idvar]]) == 0) & c(TRUE, diff(data_wf_gi$i_spell) != 0))]
    }

    data_wf_gi$spell_et <- data_wf_gi$ET_j

    # multispell data spell added above needs proper endtime
    if (exists("multispell") & multispell == 1) {
        data_wf_gi$spell_et[which(c(diff(data_wf_gi[[idvar]]) == 0, FALSE) & c(diff(data_wf_gi$i_spell) != 0, TRUE) & c(FALSE, diff(data_wf_gi$seqnum_dupl) ==
            0))] <- data_wf_gi$ET[which(c(diff(data_wf_gi[[idvar]]) == 0, FALSE) & c(diff(data_wf_gi$i_spell) != 0, TRUE) & c(FALSE, diff(data_wf_gi$seqnum_dupl) ==
            0))]
        # if i-spells has no j-events: use i-spell's starttime and endtime and do not count toward adoption
        data_wf_gi$spell_st[which(data_wf_gi$sp_notlost == 0)] <- data_wf_gi$ST[which(data_wf_gi$sp_notlost == 0)]
        data_wf_gi$spell_et[which(data_wf_gi$sp_notlost == 0)] <- data_wf_gi$ET[which(data_wf_gi$sp_notlost == 0)]
        data_wf_gi$all_case[which(data_wf_gi$sp_notlost == 0)] <- 0
    }

    data_wf_gi$spell_es <- 0
    data_wf_gi$spell_es[which(data_wf_gi$spell_et == data_wf_gi$ET & data_wf_gi$ES == 1)] <- 1

    # data is now in spell form divided by event times including 0 length spells when duplicate events occur there's a final fill of event
    # times from et_max to right censor after create covariates save(data_wf_gi, file = 'tempfile_2.RData')

    if (length(which(data_wf_gi$spell_et > data_wf_gi$ET)) > 0) {
        data_wf_gi <- data_wf_gi[-which(data_wf_gi$spell_et > data_wf_gi$ET), ]
    }

    if (exists("multispell") & multispell == 1) {
        data_wf_gi <- rbind(data_wf_gi, data_wf_gi[which(c(diff(data_wf_gi[[idvar]]) != 0, TRUE) & data_wf_gi$spell_es == 0), ])
        data_wf_gi <- data_wf_gi[order(data_wf_gi[[idvar]], data_wf_gi$ST, data_wf_gi$ET_j, data_wf_gi[[paste(idvar, "j", sep = "_")]]),
            ]
        data_wf_gi$spell_st[which(c(diff(data_wf_gi[[idvar]]) != 0, TRUE) & data_wf_gi$spell_es == 0)] <- data_wf_gi$et_max[which(c(diff(data_wf_gi[[idvar]]) !=
            0, TRUE) & data_wf_gi$spell_es == 0)]
        data_wf_gi$spell_et[which(c(diff(data_wf_gi[[idvar]]) != 0, TRUE) & data_wf_gi$spell_es == 0)] <- data_wf_gi$ET[which(c(diff(data_wf_gi[[idvar]]) !=
            0, TRUE) & data_wf_gi$spell_es == 0)]
    }

    # create diffusion counts
    data_wf_gi <- data_wf_gi[order(data_wf_gi[[idvar]], data_wf_gi$spell_et, data_wf_gi$spell_st), ]

    # count j-events up to last (in i-spell) and other spells start at end of last spell (last adoption)
    data_wf_gi$nprior_last <- with(data_wf_gi, ave(rep(1, nrow(data_wf_gi)), eval(as.symbol(idvar)), FUN = seq_along))
    nprior <- data.frame(data_wf_gi %>% group_by(eval(as.symbol(idvar))) %>% mutate(V1 = dplyr::lag(nprior_last, n = 1, default = NA)), stringsAsFactors = FALSE)
    data_wf_gi <- cbind(data_wf_gi, nprior = nprior$V1)
    data_wf_gi$nprior[is.na(data_wf_gi$nprior)] <- 0

    # since each spell links cases i to a specific prior adopter j, can form diffusion counts through equivalence relations between i and j
    # below zg_eq_`var' = 1 when there's a match, and zg_cnt_`var' cumulate these over time
    if (exists("wvars") & !is.null(wvars)) {
        for (i in 1:length(wvars)) {
            wt_var <- data.frame(data_wf_gi %>% mutate(V1 = dplyr::lag(eval(as.symbol(paste(wvars[i], "j", sep = "_"))), n = 1, default = NA)),
                stringsAsFactors = FALSE)$V1
            data_wf_gi <- cbind(data_wf_gi, wt_var)
            colnames(data_wf_gi)[dim(data_wf_gi)[2]] <- paste("wt", wvars[i], sep = "_")
            data_wf_gi[[paste("wt", wvars[i], sep = "_")]][which(c(TRUE, diff(data_wf_gi[[idvar]]) != 0))] <- NA
            data_wf_gi[[paste("wt", wvars[i], sep = "_")]][which(c(FALSE, diff(data_wf_gi$nprior) == 0))] <- NA
            data_wf_gi <- data.frame(data_wf_gi %>% group_by(idvar = eval(as.symbol(idvar))) %>% mutate(nprior_diff = c(FALSE, diff(nprior))),
                stringsAsFactors = FALSE)
            data_wf_gi[[paste("wt", wvars[i], sep = "_")]][which(data_wf_gi$nprior_diff == 0)] <- NA
            data_wf_gi[[paste("w", wvars[i], sep = "_")]] <- ave(data_wf_gi[[paste("wt", wvars[i], sep = "_")]], data_wf_gi[[idvar]], FUN = function(x) cumsum(tidyr::replace_na(x,
                0)))
            data_wf_gi <- within(data_wf_gi, rm(idvar, nprior_diff))
        }
    }

    if (exists("zgvars") & !is.null(zgvars)) {
        for (i in 1:length(zgvars)) {
            ddt4 <- data.table::data.table(data_wf_gi)
            zgvar_lag <- data.frame(data_wf_gi %>% mutate(V1 = dplyr::lag(eval(as.symbol(paste(zgvars[i], "j", sep = "_"))), n = 1, default = NA)),
                stringsAsFactors = FALSE)$V1
            data_wf_gi[[paste("zg_eq", zgvars[i], sep = "_")]] <- NA
            data_wf_gi[[paste("zg_eq", zgvars[i], sep = "_")]][which(data_wf_gi[[zgvars[i]]] == zgvar_lag)] <- 1
            idvar_lag <- data.frame(data_wf_gi %>% mutate(V1 = dplyr::lag(eval(as.symbol(idvar)), n = 1, default = NA)), stringsAsFactors = FALSE)$V1
            data_wf_gi[[paste("zg_eq", zgvars[i], sep = "_")]][which(data_wf_gi[[idvar]] != idvar_lag)] <- 0
            data_wf_gi <- data.frame(data_wf_gi %>% group_by(idvar = eval(as.symbol(idvar))) %>% mutate(nprior_diff = c(FALSE, diff(nprior))),
                stringsAsFactors = FALSE)
            data_wf_gi[[paste("zg_eq", zgvars[i], sep = "_")]][which(data_wf_gi$nprior_diff == 0)] <- NA
            data_wf_gi[[paste("zg", zgvars[i], sep = "_")]] <- ave(data_wf_gi[[paste("zg_eq", zgvars[i], sep = "_")]], data_wf_gi[[idvar]],
                FUN = function(x) cumsum(tidyr::replace_na(x, 0)))
            data_wf_gi <- within(data_wf_gi, rm(idvar, nprior_diff))
        }
    }
    # save(data_wf_gi, file = 'tempfile_3.RData')

    # delete spells for duplicate events: all covariates have now been generated from them for multispell data, keep i-spells without
    # j-events as well
    if (exists("multispell") & multispell == 1) {
        data_wf_gi <- data_wf_gi[which(data_wf_gi$spell_et > data_wf_gi$spell_st | (data_wf_gi$spell_et == data_wf_gi$spell_st & data_wf_gi$sp_notlost ==
            0)), ]
    } else if (!exists("multispell") | multispell == 0) {
        data_wf_gi <- data_wf_gi[which(data_wf_gi$spell_et > data_wf_gi$spell_st), ]
    }

    # fill last record for each id from max event time to right-censor time do not do this above because the max event time may be a
    # duplicate
    data_wf_gi$spell_et[which(data_wf_gi$spell_et >= data_wf_gi$et_max & data_wf_gi$ES == 0 & c(diff(data_wf_gi[[idvar]]) != 0, TRUE))] <- data_wf_gi$ET[which(data_wf_gi$spell_et >=
        data_wf_gi$et_max & data_wf_gi$ES == 0 & c(diff(data_wf_gi[[idvar]]) != 0, TRUE))] + 1

    # save(data_wf_gi, file = 'working_file.RData')

    # create susceptibility vector variables
    if (exists("vvars") & !is.null(vvars)) {
        for (i in 1:length(vvars)) {
            data_wf_gi[[paste("v", vvars[i], sep = "_")]] <- data_wf_gi[[vvars[i]]] * data_wf_gi$nprior
        }
    }

    if (exists("vintercept") & vintercept == 1) {
        data_wf_gi$v_intercept <- data_wf_gi$nprior
    }

    vvars <- c(vvars, "intercept")

    # run regression
    formula_1 <- "survival::Surv(spell_et, spell_es)"

    xx_xvars <- paste(xvars, collapse = " + ")

    if (exists("vvars") & !is.null(vvars)) {
        vv_vvars <- paste("v", vvars, sep = "_")
        vv_vvars <- paste(vv_vvars, collapse = " + ")
        formula_2 <- paste(xx_xvars, vv_vvars, sep = " + ")
    }
    if (exists("wvars") & !is.null(wvars)) {
        ww_wvars <- paste("w", wvars, sep = "_")
        ww_wvars <- paste(ww_wvars, collapse = " + ")
        formula_2 <- paste(formula_2, ww_wvars, sep = " + ")
    }
    if (exists("zdvars") & !is.null(zdvars)) {
        zz_zdvars <- paste("zd", zdvars, sep = "_")
        zz_zdvars <- paste(zz_zdvars, collapse = " + ")
        formula_2 <- paste(formula_2, zz_zdvars, sep = " + ")
    }
    if (exists("zgvars") & !is.null(zgvars)) {
        zz_zgvars <- paste("zg", zgvars, sep = "_")
        zz_zgvars <- paste(zz_zgvars, collapse = " + ")
        formula_2 <- paste(formula_2, zz_zgvars, sep = " + ")
    }

    formula <- paste(formula_1, formula_2, sep = " ~ ")

    fit_exp <- survival::survreg(as.formula(formula), data = data_wf_gi, dist = dist)

    return(fit_exp)
}
