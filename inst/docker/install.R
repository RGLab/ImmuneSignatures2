# Don't use package cache. https://rstudio.github.io/renv/articles/renv.html#cache-1
# Unneccesary since you only install packages once, when building container.
# Default cache adds lots of extra files under home directory which makes it really slow when
# running in LabKey.
renv::settings$use.cache(FALSE)

# For simplicity, install packages into default site library, so no need to worry about
# juggling libPaths. Not recommended outside of containerized R installations!
renv::restore(library = "/usr/local/lib/R/site-library")
remotes::install_github("rglab/immunesignatures2", dependencies = FALSE)
