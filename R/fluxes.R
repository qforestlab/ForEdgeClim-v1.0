#' Function to calculate ground heat flux G
#'
#' @param net_rad_ground Net radiation at ground surface
#' @return G Ground heat flux
#' @export
#'
# G is positive if E is entering the soil
# More dense voxels lose less heat to ground
calculate_G <- function(net_rad_ground){
  G = p_ground*(1-den)*net_rad_ground
  return(G)
}

#' Function to simulate air-to-air heat exchange in 3D.
#' Air temperatures are mapped to their correct spatial positions.
#'
#' @param temp_air Air temperature (K)
#' @param temp_soil Ground surface temperature (K)
#' @param x_coords X coordinates corresponding to input vectors
#' @param y_coords Y coordinates corresponding to input vectors
#' @param z_coords Z coordinates corresponding to input vectors
#' @return Net heat exchange per voxel (W)
#' @export
#'
calculate_D <- function(temp_air, temp_soil, x_coords, y_coords, z_coords){

  area <- voxel_length^2
  V <<- voxel_length^3

  nx <- x_dim
  ny <- y_dim
  nz <- z_dim

  # Map temperatures to their correct spatial positions
  idx <- as.integer(
    x_coords + nx * (y_coords - 1) + nx * ny * (z_coords - 1)
  )

  n_voxels <- nx * ny * nz

  if (length(temp_air) != n_voxels ||
      length(temp_soil) != n_voxels ||
      length(idx) != n_voxels ||
      anyNA(idx) ||
      any(idx < 1 | idx > n_voxels) ||
      anyDuplicated(idx)) {
    stop("Invalid voxel coordinates or incomplete grid.")
  }

  air_values <- numeric(n_voxels)
  soil_values <- numeric(n_voxels)

  air_values[idx] <- temp_air
  soil_values[idx] <- temp_soil

  T_air_grid <- array(air_values, dim = c(nx, ny, nz))
  T_soil_grid <- array(soil_values, dim = c(nx, ny, nz))

  # Initialise heat exchange arrays
  D_x <- array(0, dim = c(nx, ny, nz))
  D_y <- array(0, dim = c(nx, ny, nz))
  D_z <- array(0, dim = c(nx, ny, nz))

  # --------------------------------------------------
  # X direction: forest core and forest edge
  # --------------------------------------------------

  # Eastern forest edge (max X): exchange with macroclimate
  D_x[nx, , ] <- h * area *
    (T_air_grid[nx, , ] - macro_temp)

  # Western forest core (min X): zero-flux boundary
  D_x[1, , ] <- 0

  # Exchange between neighbouring voxels
  if (nx > 1) {
    D_x[2:nx, , ] <- D_x[2:nx, , ] +
      h * area *
      (T_air_grid[2:nx, , ] -
         T_air_grid[1:(nx-1), , ])

    D_x[1:(nx-1), , ] <- D_x[1:(nx-1), , ] +
      h * area *
      (T_air_grid[1:(nx-1), , ] -
         T_air_grid[2:nx, , ])
  }

  # --------------------------------------------------
  # Y direction: northern and southern forest core
  # --------------------------------------------------

  # Zero-flux boundary conditions at min Y and max Y
  D_y[, ny, ] <- 0
  D_y[, 1, ] <- 0

  if (ny > 1) {
    D_y[, 1:(ny-1), ] <- D_y[, 1:(ny-1), ] +
      h * area *
      (T_air_grid[, 1:(ny-1), ] -
         T_air_grid[, 2:ny, ])

    D_y[, 2:ny, ] <- D_y[, 2:ny, ] +
      h * area *
      (T_air_grid[, 2:ny, ] -
         T_air_grid[, 1:(ny-1), ])
  }

  # --------------------------------------------------
  # Z direction: upper canopy and ground
  # --------------------------------------------------

  # Upper canopy (max Z): exchange with macroclimate
  D_z[, , nz] <- h * area *
    (T_air_grid[, , nz] - macro_temp)

  # Ground (min Z): exchange with ground temperature
  D_z[, , 1] <- h * area *
    (T_air_grid[, , 1] - T_soil_grid[, , 1])

  if (nz > 1) {
    D_z[, , 2:nz] <- D_z[, , 2:nz] +
      h * area *
      (T_air_grid[, , 2:nz] -
         T_air_grid[, , 1:(nz-1)])

    D_z[, , 1:(nz-1)] <- D_z[, , 1:(nz-1)] +
      h * area *
      (T_air_grid[, , 1:(nz-1)] -
         T_air_grid[, , 2:nz])
  }

  # Total net heat exchange per voxel
  # Positive values indicate heat leaving the voxel
  total_diffusion <- D_x + D_y + D_z

  # Restore original voxel ordering
  return(as.vector(total_diffusion)[idx])
}

#' Function to calculate sensible heat transfer between forest surface and air
#'
#' @param T_surf Forest surface temperature
#' @param T_air Air temperature
#' @return H Sensible heat flux
#' @export
#'
calculate_H <- function(T_surf,T_air){
  # H is simulated via 'convective conductance' where g functions as an effective conductance that accounts for
  # both conduction within the boundary layer and convection around it.
  # H is positive when energy is lost from the surface to the air.
  H = den * g_forest * (T_surf - T_air)

  return(H)

}

#' Function for saturated vapor pressure; this is the empirical equation of Tetens
#'
#' @param temp Forest surface temperature
#' @return Saturated vapor pressure
#' @export
#'
saturated_vapor_pressure <- function(temp) {
  0.6108 * exp(17.27 * temp / (temp + 237.3)) # temp in °C, result in kPa
}

#' Function to calculate latent heat flux LE (~ evapotranspiration, ET: LE = L_v . ET with L_v = latent heat of vaporization of water (J/kg))
#'
#' @param temperature Forest surface temperature
#' @param net_rad Net radiation
#' @return LE Latent heat flux
#' @export
#'
calculate_LE <- function(temperature, net_rad, G) {
  # slope of the saturation pressure curve; temp in °C; slope in kPa/K = kPa/°C
  slope = 4098 * saturated_vapor_pressure(temperature) / (temperature + 237.3)^2
  # Latent heat flux by the empirical method of Priestly-Taylor
  LE = den * alpha_PT * (net_rad - G) * slope / (slope + gamma_psy)
  # LE cannot be negative
  LE[LE<0] = 0

  return(LE) # LE in W/m²
}
