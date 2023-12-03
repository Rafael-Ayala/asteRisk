\name{lambert}
\alias{lambert}
\title{Solve Lambert's problem to determine a transfer orbit
}
\description{
Given 2 position vectors and a time difference between them, calculates the
velocity that the object should have at the initial and final states to follow
a transfer orbit in the target time (Lambert's problem).
The transfer orbit can be elliptical, parabolic or hyperbolic, all of which are
dealt with. In the case of elliptical orbits, multi-revolution orbits are possible.
Additionally, low-path and high-path orbits are also possible. All possible solutions
are returned for the case of elliptical transfers.
The user must specify if a retrograde transfer is desired. By default, prograde
transfers are calculated.
Currently, the formulation by Dario Izzo is applied to solve Lambert's problem
(https://link.springer.com/article/10.1007/s10569-014-9587-y).
}
\usage{
lambert(initialPosition, finalPosition, initialTime, finalTime, retrogradeTransfer=FALSE, 
        centralBody="Earth", maxIterations=2000, atol=0.00001, rtol=0.000001)
}
\arguments{
  \item{initialPosition}{Vector with the 3 components of the initial position in 
  Cartesian coordinates, in meters.
  }
  \item{finalPosition}{Vector with the 3 components of the initial position in 
  Cartesian coordinates, in meters.
  }
  \item{initialTime}{Either date-time string in UTC indicating the time
  corresponding to the initial state vector of the satellite, or numeric value 
  indicating the starting time in seconds since an arbitrary reference instant. 
  If provided as a date-time string in UTC, \code{finalTime} can be provided 
  either as another date-time string in UTC (in which case the transfer time
  will be determined as the difference between the 2 date-time strings), or
  as a numeric value indicating the seconds ellapsed since \code{initialTime}.
  If provided as a numeric value, \code{finalTime} can only be provided as another
  numeric value indicating the number of seconds since the same arbitrary reference
  for objects in deep space, and also for objects near Earth if \code{targetTime}
  is provided as a date-time string.
  }
  \item{finalTime}{Either date-time string in UTC indicating the time
  corresponding to the final state vector of the satellite, or numeric value 
  indicating the final time in seconds. 
  If \code{initialTime} was provided as a date-time string in UTC, then \code{finalTime}
  can be provided as either of the two (date-time string or numeric value in seconds).
  In this case, providing \code{finalTime} as a numeric value will be interpreted
  as the number of seconds since the instant specified for \code{initialTime}.
  If \code{initialTime} was provided as a numeric value (indicating the time in
  seconds since an arbitrary reference point), then \code{finalTime} can only
  be provided as another numeric value, which will be interpreted as the number
  of seconds since the same reference instant as that used for \code{initialTime},
  and therefore the transfer time will be calculated as the difference between
  \code{initialTime} and \code{finalTime}.
  }
  \item{retrogradeTransfer}{Logical indicating if retrograde transfer orbits
  should be calculated, i.e., transfer orbits with an inclination higher than
  180º. By default, \code{retrogradeTransfer=FALSE}, and prograde transfer orbits
  are calculated.
  }
  \item{centralBody}{String indicating the central body around which the satellite
  is orbiting. Can be one of \code{c("Sun", "Mercury", "Venus", "Earth", "Moon", 
  "Mars", "Jupiter", "Saturn", "Uranus", "Neptune")} (case insensitive).
  }
  \item{maxIterations}{Maximum number of iterations to perform at the different
  root-finding steps when solving Lambert's problem.
  }
  \item{atol}{Absolute tolerance value used for the convergence criterion at the
  different root-finding steps when solving Lambert's problem.
  }
  \item{rtol}{Relative tolerance value used for the convergence criterion at the
  different root-finding steps when solving Lambert's problem.
  }
}
\value{
A list with a number of elements equal to the number of possible transfer orbits
found. For hyperbolic and parabolic transfer orbits, there will always be a single
possible transfer orbit. For elliptic transfer orbits, multi-revolution orbits
may exist if the transfer time is large enough. If they are possible, they will
be calculated together with the basic, single-revolution orbit.
Each element of the top-level list is in itself a list with the following elements:
  \item{numberRevs}{Number of complete revolutions that will be performed in this
  transfer orbit. Always equals 0 for parabolic and hyperbolic transfer orbits,
  as well as for the non-multirevolution orbit of elliptic transfers.
  }
  \item{path}{String indicating if the transfer orbit is of the high-path type
  (i.e., has its second focus located beyond the vector connecting the initial
  and final positions) or of the low-path type (second focus located between
  this vector and the first focus).
  }
  \item{orbitType}{String indicating the type of transfer orbit (elliptic,
  parabolic or hyperbolic).
  }
  \item{v1}{Velocity vector in m/s that the object should have at the start of
  the transfer orbit.
  }
  \item{v2}{Velocity vector in m/s that the object should have at the end of
  the transfer orbit.
  }
}
\references{
https://link.springer.com/article/10.1007/s10569-014-9587-y
}
\examples{
# Consider the following initial and final positions:
initialPosition <- c(15945.34, 0, 0) * 1000
finalPosition <- c(12214.83899, 10249.46731, 0) * 1000
# Given a time difference of 76 minutes between the 2 states, calculate the
# velocity that the spacecraft should have at the beginning and end of the
# transfer orbit
lambertSolution <- lambert(initialPosition, finalPosition, 0, 76*60)
length(lambertSolution)
lambertSolution[[1]]$orbitType
lambertSolution[[1]]$v1
lambertSolution[[1]]$v2

# A single transfer orbit is possible (elliptic, single-revolution orbit)
}
