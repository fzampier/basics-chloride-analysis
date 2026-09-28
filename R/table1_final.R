# table1_final.R — rebuild Table 1 on the analytic cohort (baseline + >=1 follow-up). Aggregate only.
if (!exists("BASICS_PATHS", inherits = TRUE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grepl("^--file=", .args)][[1]])
  source(file.path(dirname(normalizePath(.file)), "project_paths.R"))
}
P <- get("BASICS_PATHS", inherits = TRUE)
source(file.path(P$root, "R", "exposure_construction.R"))
load(P$data_file)
d <- dados
d$arm <- ifelse(d$grupo=="Saline Solution","Saline","Plasma-Lyte")
d$age <- as.numeric(difftime(d$dtarrol,d$dtnasc,units="days"))/365.25
d$female <- as.integer(d$sexo=="Feminino")
d$admtype <- factor(ifelse(d$proc_sepse=="Planejada","Elective",
                    ifelse(d$proc_sepse=="Não planejada, sem sepse","Non-sepsis","Sepsis")),
                    levels=c("Elective","Non-sepsis","Sepsis"))
scr <- d$creabas
kap <- ifelse(d$female==1,0.7,0.9); alp <- ifelse(d$female==1,-0.241,-0.302); sxf <- ifelse(d$female==1,1.012,1)
d$egfr <- 142*pmin(scr/kap,1)^alp*pmax(scr/kap,1)^(-1.200)*0.9938^d$age*sxf
d$creat_umol <- d$creabas*88.4
d <- add_postbaseline_exposures(d)
d$studyfl <- d$study_vol3
d$totalfl <- d$total_vol3
d$load3 <- d$study_load3

fu <- !is.na(d$clobas)&(!is.na(d$clod1)|!is.na(d$clod2)|!is.na(d$clod3))   # analytic cohort
co <- d[fu,]
cat("ANALYTIC COHORT n =",nrow(co)," (",sprintf("%.0f%%",100*nrow(co)/nrow(d)),"of",nrow(d),"trial)\n")
cat("  by arm: Saline",sum(co$arm=="Saline")," Plasma-Lyte",sum(co$arm=="Plasma-Lyte"),"\n")
p2 <- !is.na(d$clobas)&!is.na(d$clod2)
cat("PAIRED DAY-2 n =",sum(p2)," (Saline",sum(p2&d$arm=="Saline"),", PL",sum(p2&d$arm=="Plasma-Lyte"),")\n")
cat("SELECTION: mean SOFA analytic",round(mean(co$sofa,na.rm=T),1),"vs rest",round(mean(d$sofa[!fu],na.rm=T),1),
    "| %Saline analytic",sprintf("%.1f",100*mean(co$arm=="Saline")),"vs rest",sprintf("%.1f",100*mean(d$arm[!fu]=="Saline")),"\n\n")

miq <- function(x){x<-x[!is.na(x)];sprintf("%.0f [%.0f-%.0f]",median(x),quantile(x,.25),quantile(x,.75))}
S <- co[co$arm=="Saline",]; P <- co[co$arm=="Plasma-Lyte",]
row <- function(lab,fS,fP) cat(sprintf("| %s | %s | %s |\n",lab,fS,fP))
cat(sprintf("--- MARKDOWN ROWS (Saline n=%d, PL n=%d) ---\n",nrow(S),nrow(P)))
row("Age, y", miq(S$age), miq(P$age))
row("Female, No. (%)", sprintf("%d (%.0f)",sum(S$female,na.rm=T),100*mean(S$female,na.rm=T)), sprintf("%d (%.0f)",sum(P$female,na.rm=T),100*mean(P$female,na.rm=T)))
cat("    [sex NA -> Saline",sum(is.na(S$female)),", PL",sum(is.na(P$female)),"]\n")
row("SOFA score", miq(S$sofa), miq(P$sofa))
row("APACHE II score", miq(S$apache), miq(P$apache))
row("Admission, elec/non-sep/sep", paste(table(S$admtype),collapse=" / "), paste(table(P$admtype),collapse=" / "))
cat("    [admtype NA -> Saline",sum(is.na(S$admtype)),", PL",sum(is.na(P$admtype)),"]\n")
row("eGFR CKD-EPI, mL/min/1.73m2", miq(S$egfr), miq(P$egfr))
row("Baseline chloride, mmol/L", miq(S$clobas), miq(P$clobas))
row("Baseline creatinine, umol/L", miq(S$creat_umol), miq(P$creat_umol))
row("Baseline sodium, mmol/L", miq(S$sodbas), miq(P$sodbas))
row("Body weight, kg", miq(S$peso), miq(P$peso))
row("Mean arterial pressure, mmHg", miq(S$pambas), miq(P$pambas))
vp <- function(s){ on <- as.character(s$norbas) %in% c("<= 0,1 mcg/kg/min","> 0,1 mcg/kg/min"); sprintf("%d (%.0f)",sum(on),100*sum(on)/nrow(s)) }
row("Vasopressor use, No. (%)", vp(S), vp(P))
cat("    [norbas levels:", paste(levels(factor(co$norbas)),collapse=" | "), "; NA Saline",sum(is.na(S$norbas)),"PL",sum(is.na(P$norbas)),"]\n")
row("Study fluid d1-3, mL", miq(S$studyfl), miq(P$studyfl))
row("Total fluid d1-3, mL", miq(S$totalfl), miq(P$totalfl))
row("Study-fluid chloride load d1-3, mmol", miq(S$load3), miq(P$load3))
