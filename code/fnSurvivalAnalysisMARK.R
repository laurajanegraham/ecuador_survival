# This file can be edited change the models that are in the analysis. 
SurvivalAnalysisMARK <- function(sp_proc, sp_ddl){
        require(RMark)
        Phi.dot <- list(formula=~1)
        Phi.habitat <- list(formula=~habitat)
        Phi.TSM <- list(formula=~TSM)
        Phi.TSMhabitat <- list(formula=~TSM+habitat)
        p.dot <- list(formula=~1)
        #p.habitat <- list(formula=~habitat)
        cml <- create.model.list("CJS")
        mark.wrapper(cml,data=sp_proc,ddl=sp_ddl,output=FALSE)
}