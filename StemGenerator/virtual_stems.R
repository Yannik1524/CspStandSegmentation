library(mgcv)

random_circle <- function(x = 0, y = 0, r = 0.2, theta = 0.05, n = 1000){
  # calculate 8 points with theta noise
  shape <- conicfit::calculateCircle(x,y,r,steps = 9)[-1,] + rnorm(16, 0, theta)
  shape <- data.frame(x = shape[,1], y = shape[,2])
  angles <- atan2(shape$y - y, shape$x - x)
  dists <- sqrt((shape$x - x)^2 + (shape$y - y)^2)
  # sort by angle
  shape <- shape[order(angles), ]
  spline <- mgcv::gam(dists ~ s(angles, bs = "cc", k = 8), knots = list(angles = c(-pi, pi)))
  # predict full circle
  fc_angles <- sort(runif(n, -pi, pi))
  fc_dists <- predict(spline, newdata = data.frame(angles = fc_angles))
  # convert back to cartesian
  fc_x <- x + fc_dists * cos(fc_angles)
  fc_y <- y + fc_dists * sin(fc_angles)
  return(data.frame(x = fc_x, y = fc_y))
}

circle2cylinder <- function(circle, zmin, zmax){
  # convert circle to cylinder
  z <- runif(nrow(circle), zmin, zmax)
  cylinder <- data.frame(x = circle$x, y = circle$y, z = z)
  return(cylinder)
}

add_outside_noise <- function(circle, noise = 0.01){
  # add noise to the outside of the circle using an uniform offset from the center
  center <- colMeans(circle[, c("x", "y")])
  angles <- atan2(circle$y - center[2], circle$x - center[1])
  dists <- sqrt((circle$x - center[1])^2 + (circle$y - center[2])^2)
  noise_dists <- dists + abs(rnorm(nrow(circle), 0, noise))
  circle$x <- center[1] + noise_dists * cos(angles)
  circle$y <- center[2] + noise_dists * sin(angles)

  return(circle)
}

add_scene_noise <- function(cylinder, n, sd, max_distance = 2, prob = 0.5,x = 0, y = 0, r = 0.2){
  # add noise to the scene
  if(runif(1) > prob) return(cylinder)
  center <- c(x, y)
  noise_pos <- center + runif(-max_distance, max_distance, n = 2)
  noise <- data.frame(x = rnorm(n, noise_pos[1], sd = sd),
                                    y = rnorm(n, noise_pos[2], sd = sd))
  dists <- sqrt((noise$x - center[1])^2 + (noise$y - center[2])^2)
  noise <- noise[dists > r, ]
  noise$z <- runif(nrow(noise), min(cylinder$z), max(cylinder$z))
  noise <- rbind(cylinder, noise)
  return(noise)
}

rotate_cylinder <- function(cylinder, max_angle = 0){
  ax <- runif(1, -max_angle, max_angle)
  ay <- runif(1, -max_angle, max_angle)

  center <- colMeans(cylinder[, c("x", "y", "z")])
  xyz <- as.matrix(cylinder[, c("x", "y", "z")]) - center

  Rx <- matrix(c(
    1, 0, 0,
    0, cos(ax), -sin(ax),
    0, sin(ax),  cos(ax)
  ), nrow = 3, byrow = TRUE)

  Ry <- matrix(c(
    cos(ay), 0, sin(ay),
    0, 1, 0,
    -sin(ay), 0, cos(ay)
  ), nrow = 3, byrow = TRUE)

  rotated <- xyz %*% t(Ry %*% Rx) + center

  data.frame(x = rotated[, 1], y = rotated[, 2], z = rotated[, 3])
}

# function to create a branch as a random walk of points in 3d space starting from a random 
# cylinder point and going in a limited random direction away from the cylinder, with a random 
# length and radius, and a random number of points, and then adding noise to the points
add_branch <- function(cylinder, max_length = 2, max_radius = 0.05, max_points = 300, noise_sd = 0.01, direction_sd = 0.01, direction_bias = 0.8, n = sample(0:4, 1)){
  # select a random point on the cylinder
  if(n == 0) return(cylinder)
  for(i in 1:n){

    idx <- sample(1:nrow(cylinder), 1)
    start_point <- as.numeric(cylinder[idx, ])
    
    # calculate directional vector from cylinder center at the points height to the point on the cylinder
    center <- colMeans(cylinder[, c("x", "y")])
    center_point <- c(center[1], center[2], start_point[3])

    # generate a random direction away from the cylinder
    direction <- start_point - center_point
    direction <- direction / sqrt(sum(direction^2)) # normalize
    
    # generate a random length and radius for the branch
    length <- runif(1, 0, max_length)
    radius <- runif(1, 0, max_radius)
    
    # generate a random number of points for the branch
    n_points <- sample(10:max_points, 1)
    
    # create the branch points as a random walk in the direction away from the cylinder
    branch_points <- matrix(NA, nrow = n_points, ncol = 3)
    branch_points[1, ] <- start_point
    
    for(i in 2:n_points){
      step_length <- runif(1, 0, length / n_points)
      # blend outward direction with random noise, then renormalize
      step_direction <- direction * direction_bias + rnorm(3) * (1 - direction_bias)
      step_direction <- step_direction / sqrt(sum(step_direction^2))
      branch_points[i, ] <- branch_points[i - 1, ] + step_direction * step_length
      # add noise to the point
      branch_points[i, ] <- branch_points[i, ] + rnorm(3, sd = noise_sd)
    }
    
    branch <- data.frame(x = branch_points[,1], y = branch_points[,2], z = branch_points[,3])
    cylinder <- rbind(cylinder, branch)
  }
  return(cylinder)
}

simulate_occlusion <- function(cylinder, occlusion_prob = 0.5, occlusion_angular_radius = 180, x = 0, y = 0){
  if(runif(1) > occlusion_prob) return(cylinder)
  
  # convert occlusion_angular_radius to radians
  occlusion_angular_radius <- occlusion_angular_radius * pi / 180
  
  # select a random point on the cylinder
  idx <- sample(1:nrow(cylinder), 1)
  occlusion_center <- as.numeric(cylinder[idx, ])
  
  # calculate angles from the center of the cylinder to the occlusion center
  center <- c(x, y)
  occlusion_angle <- atan2(occlusion_center[2] - center[2], occlusion_center[1] - center[1])

  #calculate angles for all points from center
  angles <- atan2(cylinder$y - center[2], cylinder$x - center[1])
  
  # calculate angular distance from occlusion center to all points
  angular_dists <- abs(angles - occlusion_angle)
  angular_dists <- pmin(angular_dists, 2 * pi - angular_dists) # wrap around

  # invert angular distance to get occlusion probability
  occlusion_probs <- (1 - (angular_dists / occlusion_angular_radius))^2
  # random remove points by occlusion probability
  cylinder <- cylinder[sample(1:nrow(cylinder), round(nrow(cylinder) * (occlusion_angular_radius / (2 * pi))), prob = occlusion_probs), ]

  return(cylinder)
}

simulate_registration_error <- function(cylinder, x = 0, y = 0, max_offset = 0.1, offset_prob = 0.5){
  if(runif(1) > offset_prob) return(cylinder)
  
  offset_angular_radius <- runif(1, 0, pi/2)

  # random x and y offset 
  x_off <- runif(1,0,max_offset)
  y_off <- runif(1,0,max_offset)

  # select a random point on the cylinder
  idx <- sample(1:nrow(cylinder), 1)
  offset_center <- as.numeric(cylinder[idx, ])
  offset_angle <- atan2(offset_center[2] - center[2], offset_center[1] - center[1])
  # calculate angles from the center of the cylinder to the offset center
  center <- c(x, y)

  #calculate angles for all points from center
  angles <- atan2(cylinder$y - center[2], cylinder$x - center[1])
  
  # calculate angular distance from offset center to all points
  angular_dists <- abs(angles - offset_angle)
  angular_dists <- pmin(angular_dists, 2 * pi - angular_dists) # wrap around

  # invert angular distance to get offset probability
  offset_probs <- (1 - (angular_dists / offset_angular_radius))^2
  offset_idx <- sample(1:nrow(cylinder), round(nrow(cylinder) * (offset_angular_radius / (2 * pi))), prob = offset_probs)
  cylinder[offset_idx, "x"] <- cylinder[offset_idx, "x"] + x_off
  cylinder[offset_idx, "y"] <- cylinder[offset_idx, "y"] + y_off
  return(cylinder)
}

# set.seed(10)
# cylinder <- random_circle() |> add_outside_noise(0.02) |> circle2cylinder(0, 1) |> simulate_occlusion() |> add_branch()  |> rotate_cylinder(0.1) ; cylinder[,1:2] |> plot(asp = 1) #|> add_scene_noise(1000, 0.5, 1.5, prob = 0.5)

# set.seed(10)
# random_circle() |> add_outside_noise(0.02) |> add_scene_noise(1000, 0.5, 1.5, prob = 0.5) |> conicfit::CircleFitByPratt()
trunk_generator <- function(x_range = c(-0.5, 0.5), 
                            y_range = c(-0.5, 0.5), 
                            z_range = c(0, 1), 
                            r_range = c(0.02, 1.1), 
                            n_points_range = c(80, 1000), 
                            noise_sd = 0.01,
                            noise_prob = 0.5,
                            max_rotation_angle = 0.1, 
                            max_branch_length = 2, 
                            occlusion_prob = 0.5, 
                            occlusion_angular_radius_range = c(90,180), 
                            offset_prob = 0.5, 
                            max_offset = 0.1){
  # generate a random trunk
  x <- runif(1, x_range[1], x_range[2])
  y <- runif(1, y_range[1], y_range[2])
  r <- runif(1, r_range[1], r_range[2])
  n_points <- round(runif(1, n_points_range[1], n_points_range[2]))
  cylinder <- random_circle(x = x, y = y, r = r, n = n_points) |> 
    add_outside_noise(noise_sd) |> 
    circle2cylinder(z_range[1], z_range[2]) |> 
    simulate_occlusion(occlusion_prob = occlusion_prob, occlusion_angular_radius = round(runif(1, occlusion_angular_radius_range[1], occlusion_angular_radius_range[2]))) |> 
    add_branch(max_branch_length = max_branch_length)  |> 
    rotate_cylinder(max_angle = max_rotation_angle) |> 
    add_scene_noise(round(n_points*0.5), max_distance = 1.5, prob = noise_prob, x = x, y = y, r = r) # maximum 50% noise points
  list(cylinder = cylinder, x = x, y = y, r = r) |> return()
}