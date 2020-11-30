###CURRENT ISSUES
####1) comparison.charac() needs overall mean for continuous variables

network.charac <- function(data.nma, outcome, N, type.outcome, time){
  
  # Bind variables to function
  event0 <- NULL
  one <- NULL
  s <- NULL
  n.arms <- NULL

  # Check if number of events and participants is integer
  tmp.check1 <- data.nma$arm.data %>% select(outcome)
  tmp.check2 <- data.nma$arm.data %>% select(N)
  
  if(type.outcome == "binomial" && all(tmp.check1%%1!=0)) {
     stop('The "outcome" variable (number of events) must be an integer')
    } 
  if(all(tmp.check2%%1!=0)) {
           stop('The "N" variable (number of participants) must be an integer')
    }
  
  outcome2 <- rlang::quo(!! as.name(outcome))#this may be redundant now
  N2 <- rlang::quo(!! as.name(N))
  
    n.interventions <- data.nma$arm.data %>% select(data.nma$varname.t) %>% unique() %>% nrow()
    
    n.studies <- data.nma$arm.data %>% select(data.nma$varname.s) %>% unique() %>% nrow()
    
    n.patients <- data.nma$arm.data %>% select(N) %>% sum() 
    
    tmp1 <- data.nma$arm.data %>% 
      select(data.nma$varname.s, data.nma$varname.t) %>% 
      nest(data=c(data.nma$varname.t))
    
    cnt <- data.nma$arm.data %>% 
      select(data.nma$varname.s, data.nma$varname.t) %>% 
      count(across(data.nma$varname.s))
    
    tmp1 <- suppressMessages(bind_cols(tmp1, cnt)) %>%
      filter(n>1)
    pairs <- tmp1[1,"data"] %>% unlist %>% sort %>% combn(2)
    for(i in 2:nrow(tmp1)){
      pairs <- tmp1[i,"data"] %>% unlist %>% sort %>% combn(2) %>% cbind(pairs)
    }
    
    n.pairwise.direct <- unique(pairs, MARGIN=2) %>% ncol()
    
    edgesANDnodes <- suppressMessages(network.structure(data.nma))
    edges <- edgesANDnodes[[1]]
    nodes <- edgesANDnodes[[2]]
    net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
    
    n.pairwise <- sum(sapply(clusters(net)$csize, function(n) n*(n-1)/2))
    
    tmp2 <- data.nma$arm.data %>% select(data.nma$varname.s) %>% table()
    n.2arm.studies <- sum(tmp2 == 2)
    
    n.multiarm.studies <- sum(tmp2 > 2)
    
    tble <- tibble(Characteristic=c("Number of Interventions",
                                    "Number of Studies",
                                    "Total Number of Patients in Network",
                                    "Total Possible Pairwise Comparisons",
                                    "Total Number of Pairwise Comparisons With Direct Data",
                                    "Is the network connected?",
                                    "Number of Two-arm Studies",
                                    "Number of Multi-Arms Studies"),
                   Value= c(n.interventions %>% as.character,
                             n.studies,
                             n.patients,
                             n.pairwise,
                             n.pairwise.direct,
                             is.connected(net),
                             n.2arm.studies,
                             n.multiarm.studies))
    
    if(type.outcome=="binomial"){

      n.events <- data.nma$arm.data %>% select(outcome) %>% sum 
      
      tmp3 <- suppressMessages(data.nma$arm.data %>% 
        select(data.nma$varname.s, outcome) %>% 
        mutate(event0=(as.integer((!! outcome2)==0)), one=1) %>% 
        group_by(across(data.nma$varname.s)) %>%
        summarise(s=sum(event0), n.arms=sum(one)) %>%
        mutate(no.0=(s==0), some.0=(s!=0), all.0=(s==n.arms)))
      
      n.studies.no.0 <- tmp3$no.0 %>% sum
      
      n.studies.some.0 <- tmp3$some.0 %>% sum
      
      n.studies.all.0 <- tmp3$all.0 %>% sum
      
      add.tble <- tibble(Characteristic=c("Total Number of Events in Network",
                                          "Number of Studies With  No Zero Events",
                                          "Number of Studies With At Least One Zero Event",
                                          "Number of Studies with All Zero Events"),
                         Value= c(n.events %>% as.character,
                                   n.studies.no.0,
                                   n.studies.some.0,
                                   n.studies.all.0))
      
      return(bind_rows(tble,add.tble))
      
    }
    if(type.outcome %in% c("rate", "rate2")){
      
      n.events <- data.nma$arm.data %>% select(outcome) %>% sum 
      
      tmp3 <- data.nma$arm.data %>% 
        select(data.nma$varname.s, outcome) %>% 
        mutate(event0=(as.integer((!! outcome2)==0)), one=1) %>% 
        group_by(across(data.nma$varname.s)) %>%
        summarise(s=sum(event0), n.arms=sum(one)) %>%
        mutate(no.0=(s==0), some.0=(s!=0), all.0=(s==n.arms))
      
      n.studies.no.0 <- tmp3$no.0 %>% sum
      
      n.studies.some.0 <- tmp3$some.0 %>% sum
      
      n.studies.all.0 <- tmp3$all.0 %>% sum
      
      mean.person.time.Fup <- data.nma$arm.data %>% select(time) %>% colMeans %>% round(digits=2)
      str.mean.person.time.Fup <- ifelse(type.outcome=="rate", "Mean person-time of studies", "Mean person follow up time")
      
      add.tble <- tibble(Characteristic=c("Total Number of Events in Network",
                                          "Number of Studies With  No Zero Events",
                                          "Number of Studies With At Least One Zero Event",
                                          "Number of Studies with All Zero Events",
                                          str.mean.person.time.Fup),
                         Value= c(n.events %>% as.character,
                                  n.studies.no.0,
                                  n.studies.some.0,
                                  n.studies.all.0,
                                  mean.person.time.Fup %>% as.character))
      
      return(bind_rows(tble,add.tble))
      
    }
    if (type.outcome=="continuous"){
      
      just.N <- data.nma$arm.data %>% select(N)
      just.outcome <- data.nma$arm.data %>% select(outcome)
      
      total <- sum(just.N*just.outcome)
      
      w.average <- total/n.patients
      
      add.tble <- tibble(Characteristic="Average Outcome",
                         Value=signif(w.average, digits=4) %>% as.character)
      
      return(bind_rows(tble,add.tble))
      }
}
 
intervention.charac <- function(data.nma, outcome, N, type.outcome, time=NULL) {
  
  # Bind variables to the function
  w.outcome <- NULL
  person.time.fup <- NULL
  
  outcome2 <- rlang::quo(!! as.name(outcome))
  N2 <- rlang::quo(!! as.name(N))
  
  if ("n" %in% colnames(data.nma$arm.data)) {
  n.studies <- data.nma$arm.data %>% select(-n) %>% count(across(data.nma$varname.t)) %>% rename(n.studies = n) 
  } else {
    n.studies <- data.nma$arm.data %>% count(across(data.nma$varname.t)) %>% rename(n.studies = n) 
  }
  
  n.patients <- data.nma$arm.data %>% 
    group_by(across(data.nma$varname.t)) %>% 
    summarise(n.patients = as.integer(sum(!! N2)))
  
  tmp.tbl <- left_join(n.studies, n.patients, by=data.nma$varname.t)

  if(type.outcome %in% c("bin","binom","binomial","binary")){  
  
  n.events <- data.nma$arm.data %>% 
    group_by(across(data.nma$varname.t)) %>% 
    summarise(n.events = as.integer(sum(!! outcome2))) 
  
   tmp.rate <- data.nma$arm.data %>% 
    mutate(tmp.rate = (!! outcome2) / (!! N2)) %>% 
    group_by(across(data.nma$varname.t)) %>% 
    summarise(min.outcome = min(tmp.rate), max.outcome = max(tmp.rate))
  
   tmp.tbl <- left_join(n.studies, n.events, by=data.nma$varname.t) %>% 
     left_join(., n.patients, by=data.nma$varname.t) %>%
     left_join(., tmp.rate, by=data.nma$varname.t) %>% 
     mutate(av.outcome=n.events/n.patients) 
  
  } else if(type.outcome %in% c("cont", "continuous")){
    
    tmp.outcome <- data.nma$arm.data %>% 
     mutate(tmp.outcome = (!! outcome2), w.outcome=(!! outcome2)*(!! N2)) %>% 
     group_by(across(data.nma$varname.t)) %>% 
     summarise(min.outcome = min(tmp.outcome), max.outcome = max(tmp.outcome), w.outcome=sum(w.outcome))
    
    tmp.tbl <- left_join(n.studies, n.patients, by=data.nma$varname.t) %>% 
     left_join(., tmp.outcome, by=data.nma$varname.t) %>%
     mutate(av.outcome=w.outcome/n.patients) %>%
     select(-w.outcome)
    
  } else if (type.outcome =="rate"){
    time2 <- rlang::quo(!! as.name(time))
    
    person.time <- data.nma$arm.data %>% 
      group_by(across(data.nma$varname.t)) %>% 
      summarise(person.time.fup = as.integer(sum(!! time2)))

    n.events <- data.nma$arm.data %>%
      group_by(across(data.nma$varname.t)) %>%
      summarise(n.events = as.integer(sum(!! outcome2)))
    
    tmp.proportion <- data.nma$arm.data %>% 
      mutate(tmp.proportion = (!! outcome2) / (!! N2)) %>% 
      group_by(across(data.nma$varname.t)) %>% 
      summarise(min.proportion = min(tmp.proportion), max.proportion = max(tmp.proportion))
    
    tmp.rate <- data.nma$arm.data %>% 
      mutate(tmp.rate = (!! outcome2) / (!! time2)) %>% 
      group_by(across(data.nma$varname.t)) %>% 
      summarise(min.event.rate = min(tmp.rate), max.event.rate = max(tmp.rate))
    
    tmp.tbl <- left_join(n.studies, n.events, by=data.nma$varname.t) %>% 
      left_join(., n.patients, by=data.nma$varname.t) %>%
      #left_join(., tmp.proportion, by=data.nma$varname.t) %>%
      left_join(., person.time, by=data.nma$varname.t) %>% 
      left_join(., tmp.rate, by=data.nma$varname.t) %>% 
      mutate(events.per.person=n.events/n.patients) %>%
      mutate(av.event.rate=n.events/person.time.fup)
  } else if (type.outcome =="rate2"){
      time2 <- quo(!! as.name(time))
      
      person.time <- data.nma$arm.data %>% 
        group_by(across(data.nma$varname.t)) %>% 
        summarise(person.time.fup = as.integer(sum((!! time2)*(!! N2))))
      
      n.events <- data.nma$arm.data %>%
        group_by(across(data.nma$varname.t)) %>%
        summarise(n.events = as.integer(sum(!! outcome2)))
      
      tmp.rate <- data.nma$arm.data %>% 
        mutate(tmp.rate = (!! outcome2) / ((!! time2)*(!! N2))) %>% 
        group_by(across(data.nma$varname.t)) %>% 
        summarise(min.event.rate = min(tmp.rate), max.event.rate = max(tmp.rate))
      
      tmp.tbl <- left_join(n.studies, n.events, by=data.nma$varname.t) %>% 
        left_join(., n.patients, by=data.nma$varname.t) %>%
        #left_join(., tmp.proportion, by=data.nma$varname.t) %>%
        left_join(., person.time, by=data.nma$varname.t) %>% 
        left_join(., tmp.rate, by=data.nma$varname.t) %>% 
        mutate(events.per.person=n.events/n.patients) %>%
        mutate(overall.event.rate=n.events/person.time.fup)
  }
  
  return(tmp.tbl)
}

comparison.charac <- function(data.nma, outcome, N, type.outcome, time=NULL) {
  
  # Binding Variables to function
  comparison <- NULL
  patient_time <- NULL

  tmp1 <- by.comparison(data.nma=data.nma, outcome=outcome, type.outcome, N=N, time=time)
  
  add.patients <- tmp1 %>% select(paste0(N,".e"), paste0(N,".c")) %>% rowSums()
  add.patients <- cbind.data.frame(tmp1$comparison, add.patients, stringsAsFactors = FALSE)
  names(add.patients) <- c("comparison", "n.patients")
  n.patients <- add.patients %>% group_by(comparison) %>% summarise(n.patients = as.integer(sum(n.patients)))
  
  add.outcomes <- tmp1 %>% select(paste0(outcome,".e"), paste0(outcome,".c")) %>% rowSums()
  add.outcomes <- cbind.data.frame(tmp1$comparison, add.outcomes, stringsAsFactors = FALSE)
  names(add.outcomes) <- c("comparison", "n.outcomes")
  n.outcomes <- add.outcomes %>% group_by(comparison) %>% summarise(n.outcomes = as.integer(sum(n.outcomes)))
  
  n.studies <- tmp1 %>% count(comparison) %>% rename(n.studies=n)
  
  tmp2 <- left_join (n.studies, n.patients, by="comparison") 
  
  if (type.outcome %in% c("bin","binom","binomial","binary")){  
    tmp2 %<>% left_join(., n.outcomes, by="comparison")%>% 
      mutate(proportion=n.outcomes/n.patients)
  } else if(type.outcome == "rate"){
    add.time <- tmp1 %>% select(paste0(time,".e"), paste0(time,".c")) %>% rowSums()
    add.time <- cbind.data.frame(tmp1$comparison, add.time, stringsAsFactors = FALSE)
    names(add.time) <- c("comparison", "n.time")
    
    n.time <- add.time %>% group_by(comparison) %>% summarise(patient_time = as.integer(sum(n.time)))
    
    tmp2 %<>% left_join(.,n.outcomes, by="comparison")%<>% 
      left_join(., n.time, by="comparison") %>%
      mutate(proportion=n.outcomes/n.patients)%>%
      mutate(event.rate=n.outcomes/patient_time)
    
  } else if(type.outcome == "rate2"){
    add.time <- tmp1 %>% select(paste0(time,".e"), paste0(time,".c")) %>% rowSums()
    add.time <- cbind.data.frame(tmp1$comparison, add.time, stringsAsFactors = FALSE)
    names(add.time) <- c("comparison", "n.time")
    
    n.time <- add.time %>% group_by(comparison) %>% summarise(patient_time = as.integer(sum(n.time)))
    
    tmp2 %<>% left_join(.,n.outcomes, by="comparison")%<>% 
      left_join(., n.time, by="comparison") %>%
      mutate(proportion=n.outcomes/n.patients)%>%
      mutate(patient_time=patient_time * n.patients)%>%
      mutate(event.rate=n.outcomes/patient_time)
    
  } else if (type.outcome =="continuous"){
    #total.outcome <- tmp1%>% 
    #mutate(total.outcome=(paste0(outcome,".e")*paste0(N,".e")+ paste0(outcome,".c")*paste0(N,".c"))/(paste0(N,".e")+paste0(N,".c") ))
    
    #tmp2 %<>% left_join(., n.outcomes, by="comparison")                    
    
  }
  return(tmp2)
}

#' Generate Network Characteristics
#' @description Generates tables of network characteristics
#' @param data A \code{BUGSnetData} object produced by \code{data.prep()}
#' @param outcome A string indicating the name of the outcome variable
#' @param N A string indicating the name of the variable containing the number of participants in each arm
#' @param type.outcome A string. Options are: "binomial", "continuous", "rate" (e.g # of events and # person-time reported), 
#' "rate2" (e.g # events and followup time reported)
#' @param time A string required when type.outcome = "rate" or "rate2". Name of variable 
#'   indicating person-time followup (e.g person years) or study followup time
#' @return \code{network} - Table of network characteristics
#' @return \code{intervention} - Summary statistics broken down by treatment
#' @return \code{comparison} - Summary statistics broken down by treatment comparison
#' @examples
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#' varname.t = "Treatment", 
#' varname.s = "Study")
#' 
#' network.char <- net.tab(data = diabetes.slr,
#' outcome = "diabetes",
#' N = "n",
#' type.outcome = "rate2",
#' time = "followup")
#' @seealso \code{\link{data.prep}}
#' @export



net.tab <- function(data, outcome, N, type.outcome, time=NULL){
  return(list(network = network.charac(data, outcome, N, type.outcome,time),
              intervention = intervention.charac(data, outcome, N, type.outcome, time),
              comparison = comparison.charac(data, outcome, N, type.outcome, time)))
  
  if(class(data) != "BUGSnetData")
    stop("\'data\' must be a valid BUGSnetData object created using the data.prep function.")
}


