context("Network plot")

data(diabetes)
prepped1 <- data.prep(diabetes, "Treatment", "Study")
struct1 <- network.structure(prepped1)
struct1$edges[,"from"] <- paste0(struct1$edges$from,"+bla")
struct1$edges[,"to"] <- paste0(struct1$edges$to,"+bla")
struct1$nodes[,"Treatment"] <- paste0(struct1$nodes$Treatment,"+bla")

diabetes$Treatment <- paste0(diabetes$Treatment,"+bla")
prepped2 <- data.prep(diabetes, "Treatment", "Study")
struct2 <- network.structure(prepped2)

test_that("Special characters produce correct network structure", 
          {expect_identical(struct1, struct2)})

data1 <- nma.model(prepped1, outcome = "diabetes", N="n", effects = "fixed", 
          family = "binomial", link = "logit", reference = "Placebo")$data

data2 <- nma.model(prepped2, outcome = "diabetes", N="n", effects = "fixed", 
                   family = "binomial", link = "logit", reference = "Placebo+bla")$data

test_that("Special characters produce correct nma.model data", 
          {expect_identical(data1, data2)})
