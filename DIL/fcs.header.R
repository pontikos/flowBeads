#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

#The function ‘read.FCSheader’ works with the output of the FACS
#machine software from a number of vendors (FCS 2.0, FCS 3.0 and
#List Mode Data LMD). The output of the function is the TEXT
#section of the FCS files. The user can specify some keywords to
#limit the output to the information of interest.

option_list <- list( 
make_option(c("-f","--fcs"), default=NULL, help = "fcsfile to parse"),
make_option(c("-l","--files"), default=NULL, help = "file containing list of fcs files"),
make_option(c("-k","--keywords"), default=NULL, help="list of parameters [default %default]")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

if (is.null(opt$fcs) && is.null(opt$files)) {
	stop("Not FCS file specified on command line!")
}

if (! is.null(opt$keywords)) {
	unlist(strsplit(opt$keywords, ",")) -> keywords
} else {
	keywords <- NULL
}

#FCSversion : 3 
#FIL : ptpn22_14th_pair_6beads_002.fcs 
#SYS : Windows XP 5.1 
#TOT : 5177                
#PAR : 16 
#MODE : L 
#CREATOR : BD FACSDiva Software Version 6.2 
#TUBE NAME : 6beads_002 
#SRC : ptpn22_14th_pair 
#EXPERIMENT NAME : tims_cd25_p6_plus_sort_cd31_plus_tim7_PLUSP8 
#GUID : 12c398ec-fdac-449e-bb2d-2dd9546768d7 
#DATE : 03-NOV-2011 
#BTIM : 19:41:27 
#ETIM : 19:41:50 
#CYT : LSRFortessa 
#CYTNUM : H7B200001 
#WINDOW EXTENSION : 4.00 
#EXPORT USER NAME : Administrator 
#EXPORT TIME : 06-NOV-2011-19:13:05 

headers <- c(
'FCSversion',
'FIL',
'SYS',
'TOT',
'PAR',
'MODE',
'CREATOR',
'TUBE NAME',
'SRC',
'EXPERIMENT NAME',
'GUID',
'DATE',
'BTIM',
'ETIM',
'CYT',
'CYTNUM',
'WINDOW EXTENSION',
'EXPORT USER NAME',
'EXPORT TIME'
)

#as.Date("07/14/2010 12:31:27 PM", format="%d/%M/%Y %H:%M:%S")
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

if (!is.null(opt$files)) {
	files <- as.character(read.csv(opt$files)[,1])
	suppressWarnings(read.FCSheader(files[[1]], keyword=keywords)) -> fcs.header
	fcs.header[[files[[1]]]] -> h
	names(h) <- sub('$', '', names(h),  fixed=T)
	parameter.names <- grep('^P\\d+N$',names(h),value=TRUE)
	cat(c(headers,parameter.names), sep="\t")
	cat('\n')
	for (fcs in files) {
		suppressWarnings(read.FCSheader(fcs, keyword=keywords)) -> fcs.header
		fcs.header[[fcs]] -> h
		names(h) <- sub('$', '', names(h),  fixed=T)
		cat(h[c(headers,parameter.names)], sep="\t")
		cat('\n')
	}
} else if (!is.null(opt$fcs)) {
	fcs <- opt$fcs
	suppressWarnings(read.FCSheader(fcs, keyword=keywords)) -> fcs.header
	fcs.header[[fcs]] -> h
	names(h) <- sub('$', '', names(h),  fixed=T)
	parameter.names <- grep('^P\\d+N$',names(h),value=TRUE)
	cat(c(headers,parameter.names), sep="\t")
	cat('\n')
	cat(h[c(headers,parameter.names)], sep="\t")
	cat('\n')
}


