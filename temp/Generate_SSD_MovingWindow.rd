 \name{Generate_SSD_MovingWindow}
 \alias{Generate_SSD_MovingWindow}
 \title{Generate SNP set data file (SSD) }
 \description{
 
 Generate a SNP set data file (SSD) from binary plink formated data files using moving window method.
 If you want to use plink formated data files, you must generate the SSD files first. 

 }
 \usage{
Generate_SSD_MovingWindow(File.Bed, File.Bim, File.Fam, File.SSD
, File.Info, WindowSize,Overlap)
 }
\arguments{
      \item{File.Bed}{ the name of the binary ped file (BED).}
      \item{File.Bim}{ the name of the binary map file (BIM).}
      \item{File.Fam}{ the name of the FAM file (FAM).}
      \item{File.SSD}{ the name of the SSD file generated. }
      \item{File.Info}{ the name of the SSD info file generated. }
      \item{WindowSize}{ Size of moving windows. }
      \item{Overlap}{ Overlap for each windows. }
}

\details{
 
 The SSD file is binary formated file with genotype informations. 
 The SSD infor file is a text file. The first 6 rows have general information of data and SNP sets. The information of each set can be found from the 8th row. 
                                                             
}


\author{Seunggeun Lee, Larisa Miropolsky}

