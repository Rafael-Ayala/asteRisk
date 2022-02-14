\name{readOEM}
\alias{readOEM}
\title{Read an Orbital Ephemeris Message file
}
\description{
OEM (Orbital Ephemeris Message) is one of the three standard file formats defined
by the CCSDS for transferring spacecraft orbit information. OEM files contain
the position and velocity of a given object at multiple times (epochs). They can
also contain optionally acceleration values, covariance matrixes that indicate
the uncertainty of the provided state vectors and other additional information.
This function reads OEM files, retrieving also the optional fields.
}
\usage{
readOEM(filename)
}
\arguments{
  \item{filename}{Path to the OEM file.
  }
}
\value{
A list with two elements. The first element, named \code{header}, is a list with
the following elements:
  \item{OMVersion }{Version of the OEM format used in the file}
  \item{creationDate }{Date of creation of the file}
  \item{creator }{Individual or organization that generated the file}

The second element is named \code{dataBlocks}, and it contains one element for
each ephemeris data block found in the OEM file. Each of these elements is a 
list with the following elements that provide information about the ephemerides
of object (note that some elements are not mandatory and therefore might not be
present in all OEM files; in these cases, their value is set to NULL):
  \item{objectName }{Name of the object}
  \item{objectID }{Object identifier for the object. Frequently, although not
  always, the identifier has the format YYYY-NNNP{PP}, where YYYY is the year
  of launch, NNN is the three-digit serial number specifying the launch number
  during year YYYY and P{PP} is a part specifier comprising 1 to 3 capital letters
  that indicate the part of the object put into space during the launch.}
  \item{referenceFrame }{Frame of reference in which ephemerides are provided.}
  \item{refFrameEpoch }{Epoch for the frame of reference, for cases where it is
  not intrinsic to the frame itself, such as TEME.}
  \item{centerName }{Name of the center of coordinates. For example, a celestial
  body, the barycenter of the entire Solar System or even other spacecraft.}
  \item{timeSystem }{Time system used for the ephemerides, covariance matrixes
  and all other time fields of the data block.}
  \item{startTime }{Start time of the time span covered by the ephemerides and
  covariance matrixes in this data block.}
  \item{endTime }{End time of the time span covered by the ephemerides and
  covariance matrixes in this data block.}
  \item{usableStartTime }{Start time of the usable time span covered by the 
  ephemerides and covariance matrixes in this data block.}
  \item{usableEndTime }{End time of the usable time span covered by the 
  ephemerides and covariance matrixes in this data block.}
  \item{interpolationMethod }{Recommended interpolation method to calculate
  ephemerides at epochs between those directly given in the data block.}
  \item{interpolationOrder }{Recommended interpolation degree to calculate
  ephemerides at epochs between those directly given in the data block.}
  \item{mass }{Mass in kg of the object.}
  \item{dragArea }{Effective area of the object subjected to drag, in square
  meters.}
  \item{dragCoefficient }{Drag coefficient of the object.}
  \item{solarRadArea }{Effective area of the object subjected to solar radiation
  pressure, in square meters.}
  \item{solarRadCoefficient }{Solar radiation pressure coefficient of the object.}
  \item{ephemerides }{Data frame with the 7 or 10 columns providing the ephemerides
  for the object. The 1st column provides the epochs for each ephemeris; columns
  2 to 4 provide the X, Y and Z components of position (in km), and columns 5 to
  7 provide the X, Y and Z components of velocity (in km/s). Columns 8 to 10 are
  optional, and if present provide the X, Y and Z components of acceleration (in
  km/s2)}
  \item{covarianceMatrixes }{List where each element is a 3-element list that
  provides a covariance matrix for this data block. Each of the 3-element lists
  corresponding to a covariance matrix contains the following elements:
  \itemize{
    \item{epoch }{Epoch of the navigation solution related to this covariance
    matrix.}
    \item{referenceFrame }{Reference frame for the covariance matrix. Frequently
    this is the same as for the ephemerides. In order to facilitate interpretation
    of the covariance matrix, conversion to a perifocal frame of reference might
    be advisable.}
    \item{covarianceMatrix }{Covariance matrix that provides information about
    the uncertainties of position and velocities. This is a symmetric 6x6 matrix
    where the values in the diagonal are the squared standard deviations of each 
    variable, and the other values are covariances between 2 variables (those 
    corresponding to the row and column of each value). The rows and columns
    correspond to the following variables, in the specified order: X position,
    Y position, Z position, X velocity, Y velocity and Z velocity.}
  }}
}
\references{
https://public.ccsds.org/Pubs/502x0b2c1e2.pdf
https://spotthestation.nasa.gov/trajectory_data.cfm
}
\examples{
# The file testOEM.txt provided with the package includes ephemerides data
# for the ISS publicly available

testOEM_ISS <- readOEM(paste0(path.package("asteRisk"), "/testOEM.txt"))
testOEM_ISS$header
testOEM_ISS$dataBlocks[[1]]$objectName
head(testOEM_ISS$dataBlocks[[1]]$ephemerides)
}
