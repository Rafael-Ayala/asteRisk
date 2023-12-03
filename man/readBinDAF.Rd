\name{readBinDAF}
\alias{readBinDAF}
\title{Read a generic binary DAF file
}
\description{
DAF (Double Precision Array File) is a binary file architecture designed to store arrays
of double precision arrays used by SPICE, NAIF toolkit software library. The
architecture forms the basis of multiple file formats used to store different data
related to astrodynamics, such as SPK, PCK and CK files (https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html).

DAF files provide a generic architecture onto which more specific file formats
are implemented. They are organized in records of fixed length (1024 bytes) containing 
different information. The first record is called the file record, and contains 
global metadata for the file. This is followed by an optional block comprising 
any number of comment records. After this, the file consits of sets of summary records, 
name records and element records. These are structured as  blocks of 1 summary record 
(which contains multiple array summaries, providing metadata about each array), 
followed by 1 name record (comprising names for the corresponding arrays whose 
summaries were in the previous summary record) and finally by as many element 
records as required to store the arrays described in the corresponding summary 
records. For a detailed description of the DAF architecture, see NAIF
documentation (https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html).

This function allows to read any binary file with a DAF architecture in a generic way,
without applying any specific formatting to its contents. Note that this will result 
in just a list of the file global metadata, comments and different arrays (each
one with a summary, a name and a set of elements). The number of elements and
meaning of each included in each array, as well as the descriptor integers and doubles
contained in the summaries, vary for each specific file type. Therefore, reading
a DAF file generically is not likely to bring much meaningful information, unless
a precise understanding of the specific file format is taken into account later.
}
\usage{
readBinDAF(filename)
}
\arguments{
  \item{filename}{Path to the binary DAF file.
  }
}
\value{
A list with three elements. The first element, named \code{metadata}, corresponds
to the "file record", which is always the first record in a DAF file, and is a 
list with the following metadata elements:
  \item{fileType }{String of the format DAF/XXXX indicating the specific file format
  for the DAF file}
  \item{numDoublesSummary }{Number of double precision numbers in each array summary}
  \item{numIntsSummary }{Number of integers in each array summary}
  \item{summaryRecordSizeDoubles }{Size of each array summary in doubles}
  \item{numCharsName }{Number of characters in each array name}
  \item{description }{String with an internal name or description of the DAF file}
  \item{firstSummaryRecNum }{Integer indicating which is the 1st summary record
  This can be used to infer the number of comment files. For example, a value of 3
  indicates that the 3rd record is the 1st summary record. Since the 1st record
  is always the file record, this means there is 1 comment record}
  \item{lastSummaryRecNum }{Integer indicating which is the last summary record}
  \item{firstFreeAddress }{Integer indicating the first free address (in bytes)
  of the file, i.e., the address of the last byte of the last element record plus 1}
  \item{endianString }{String indicating the endianness of the file. This is
  automatically taken into account when reading the file. It can be either LTL-IEEE
  (little endianness) or BIG-IEEE (big endianness), and is determined by the architecture
  of the system where the file was written}
  \item{ftpString }{String used to verify integrity of the DAF file. It should be
  exactly equal to \code{"FTPSTR:\r:\n:\r\n:\r:\x81:\020\xce:ENDFTP"}}

The second element is named \code{comments}, and is a character vector where
each element is a line of comments.

The third element is named \code{arrays}, and is a nested list where each top-level
element represents one of the arrays stored in the DAF file and its associated 
metadata. Each of the top-level elements is itself a list with the following 3 
elements:
  \item{arrayName }{String with the name of the array}
  \item{arraySummary }{A list with the multiple doubles and integers that are stored
  in each array summary of summary records and which provide metadata describing
  each array. The elements are named as Double1, Double2, ..., DoubleN; Integer1,
  Integer2, ..., Integer(M-2) and finally initialArrayAddress and finalArrayAddress. 
  N and M are respectively the number of doubles and integers in each array summary, 
  and are given in elements \code{numDoublesSummary} and \code{numIntsSummary} of 
  the \code{metadata} element of the top-level list. Note that the number of doubles 
  and integers describing each array, as well as the meaning of each, varies between 
  different specific file formats, and therefore no exact meaning can be derived 
  when simply reading the file as a generic DAF file. The exception to this are 
  the last 2 integers, which always are respectively the initial and final addresses
  of the elements corresponding to the array within the DAF file, in double precision 
  numbers (and therefore, in order to obtain byte addresses it must be multiplied by 8 
  and subtract 7)}
  \item{arrayElements }{A numeric vector with all the elements of the array. Note
  that this includes potentially constants, actual data and additional array metadata.
  Furthermore, the number, order and meaning of the elements differs greatly between
  different specific types and subtypes of DAF files, and therefore it is hard to
  extract any meaningful information without knowledge of the internal organization 
  of each array}
}
\references{
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
}
\examples{
# The file vgr2_jup230.bsp distributed by NAIF includes information for the 
# Jupiter flyby of Voyager 2

testDAFFilePath <- tempfile()
download.file("https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/spk/vgr2_jup230.bsp", testDAFFilePath)
testDAF <- readBinDAF(testDAFFilePath)
testDAF$metadata
# The file seems to be of type SPK
testDAF$comments
length(testDAF$arrays)
# It contains a single array
}
